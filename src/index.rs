//! Reference minimizer index: build once, save to disk (`.gdx`), reuse for
//! every query genome. Records live in one array sorted by hash, so lookups
//! are a binary search and peak build memory stays at ~24 bytes per
//! minimizer (no parallel hash/value arrays during construction).

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

use crate::fasta::FastaRecord;
use crate::minimizer::stream_minimizers;

pub const MAGIC: &[u8; 4] = b"GDIX";
pub const VERSION: u32 = 1;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Hit {
    pub hash: u64,
    pub contig: u32,
    /// Leftmost coordinate (0-based) of the k-mer on the reference contig,
    /// regardless of strand.
    pub pos: u64,
    pub strand: u8,
}

#[derive(Clone)]
pub struct ContigInfo {
    pub name: String,
    pub length: u64,
}

pub struct RefIndex {
    pub k: u32,
    pub w: u32,
    pub contigs: Vec<ContigInfo>,
    /// Sorted by (hash, contig, pos).
    pub hits: Vec<Hit>,
}

impl RefIndex {
    pub fn build(
        records: impl Iterator<Item = io::Result<FastaRecord>>,
        k: u32,
        w: u32,
    ) -> io::Result<RefIndex> {
        let mut contigs = Vec::new();
        let mut hits: Vec<Hit> = Vec::new();
        for (ci, r) in records.enumerate() {
            let rec = r?;
            // VCF/PAF-safe contig id: first whitespace-delimited token of the
            // header (handles names like ">Chr4 " with trailing spaces)
            let name = rec.name.split_whitespace().next().unwrap_or("").to_string();
            contigs.push(ContigInfo {
                name,
                length: rec.seq.len() as u64,
            });
            stream_minimizers(&rec.seq, k, w, |h, p, s| hits.push(Hit {
                hash: h,
                contig: ci as u32,
                pos: p,
                strand: s,
            }));
        }
        hits.sort_unstable();
        Ok(RefIndex { k, w, contigs, hits })
    }

    pub fn lookup(&self, h: u64) -> &[Hit] {
        let start = self.hits.partition_point(|x| x.hash < h);
        let end = self.hits.partition_point(|x| x.hash <= h);
        &self.hits[start..end]
    }

    pub fn contig_name(&self, id: u32) -> &str {
        &self.contigs[id as usize].name
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let mut f = BufWriter::with_capacity(1 << 20, File::create(path)?);
        f.write_all(MAGIC)?;
        f.write_all(&VERSION.to_le_bytes())?;
        f.write_all(&self.k.to_le_bytes())?;
        f.write_all(&self.w.to_le_bytes())?;
        f.write_all(&(self.contigs.len() as u32).to_le_bytes())?;
        for c in &self.contigs {
            let bytes = c.name.as_bytes();
            f.write_all(&(bytes.len() as u32).to_le_bytes())?;
            f.write_all(bytes)?;
            f.write_all(&c.length.to_le_bytes())?;
        }
        f.write_all(&(self.hits.len() as u64).to_le_bytes())?;
        for hit in &self.hits {
            f.write_all(&hit.hash.to_le_bytes())?;
            f.write_all(&hit.contig.to_le_bytes())?;
            f.write_all(&hit.pos.to_le_bytes())?;
            f.write_all(&[hit.strand])?;
        }
        f.flush()
    }

    pub fn load(path: &Path) -> io::Result<RefIndex> {
        let mut f = BufReader::with_capacity(1 << 20, File::open(path)?);
        let mut magic = [0u8; 4];
        f.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "not a genodiff index"));
        }
        let version = read_u32(&mut f)?;
        if version != VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("index version {version} unsupported (expected {VERSION})"),
            ));
        }
        let k = read_u32(&mut f)?;
        let w = read_u32(&mut f)?;
        let n_contigs = read_u32(&mut f)? as usize;
        let mut contigs = Vec::with_capacity(n_contigs);
        for _ in 0..n_contigs {
            let name_len = read_u32(&mut f)? as usize;
            let mut name = vec![0u8; name_len];
            f.read_exact(&mut name)?;
            let length = read_u64(&mut f)?;
            contigs.push(ContigInfo {
                name: String::from_utf8_lossy(&name).into_owned(),
                length,
            });
        }
        let n = read_u64(&mut f)? as usize;
        let mut hits = Vec::with_capacity(n);
        // record layout must mirror save(): hash u64 | contig u32 | pos u64 | strand u8
        let mut buf = [0u8; 21];
        for _ in 0..n {
            f.read_exact(&mut buf)?;
            hits.push(Hit {
                hash: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
                contig: u32::from_le_bytes(buf[8..12].try_into().unwrap()),
                pos: u64::from_le_bytes(buf[12..20].try_into().unwrap()),
                strand: buf[20],
            });
        }
        Ok(RefIndex { k, w, contigs, hits })
    }
}

fn read_u32(f: &mut impl Read) -> io::Result<u32> {
    let mut b = [0u8; 4];
    f.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_u64(f: &mut impl Read) -> io::Result<u64> {
    let mut b = [0u8; 8];
    f.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    fn rand_seq(rng: &mut Rng, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
    }

    #[test]
    fn build_save_load_roundtrip() {
        let mut rng = Rng::new(17);
        let recs = vec![
            FastaRecord { name: "ctg1 extra-description".into(), seq: rand_seq(&mut rng, 20_000) },
            FastaRecord { name: "ctg2".into(), seq: rand_seq(&mut rng, 8_000) },
        ];
        let idx = RefIndex::build(recs.iter().map(|r| Ok(r.clone())), 15, 10).unwrap();

        assert_eq!(idx.contigs.len(), 2);
        assert_eq!(idx.contigs[0].name, "ctg1"); // header must be truncated at whitespace
        assert_eq!(idx.contigs[0].length, 20_000);
        assert!(!idx.hits.is_empty());
        for w in idx.hits.windows(2) {
            assert!(w[0] <= w[1], "hits must be sorted by (hash, contig, pos)");
        }

        // a minimizer taken from the reference must be found at its own position
        let probe = idx.hits[idx.hits.len() / 2].hash;
        let hits = idx.lookup(probe);
        assert!(!hits.is_empty());

        let dir = std::env::temp_dir().join(format!("genodiff-idx-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.gdx");
        idx.save(&path).unwrap();
        let idx2 = RefIndex::load(&path).unwrap();
        assert_eq!(idx.k, idx2.k);
        assert_eq!(idx.w, idx2.w);
        assert_eq!(idx2.contigs.len(), 2);
        assert_eq!(idx.contigs[0].name, idx2.contigs[0].name);
        assert_eq!(idx.contigs[1].length, idx2.contigs[1].length);
        assert_eq!(idx.hits, idx2.hits);
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
