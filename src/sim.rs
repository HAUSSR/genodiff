//! Variant simulator: take a reference genome, inject known SNP/InDel/SV
//! edits, output a mutated query genome plus a truth VCF. This is the
//! ground-truth generator used to validate every pipeline stage.

use std::io;

use crate::fasta::FastaRecord;
use crate::minimizer::{base_code, revcomp};
use crate::rng::Rng;
use crate::vcf::VcfWriter;

/// minimum spacing between planned variants
pub const GAP: u64 = 100;
/// keep variants away from contig ends (also guarantees an anchor base for
/// explicit InDel VCF records)
pub const MARGIN: u64 = 200;

#[derive(Clone)]
pub enum VarKind {
    Snp {
        pos: u64,
        alt: u8,
    },
    /// delete [pos, pos+len)
    Del {
        pos: u64,
        len: u64,
    },
    /// insert `seq` between positions pos-1 and pos
    Ins {
        pos: u64,
        seq: Vec<u8>,
    },
    Inv {
        pos: u64,
        len: u64,
    },
    /// tandem duplication of [pos, pos+len)
    Dup {
        pos: u64,
        len: u64,
    },
}

impl VarKind {
    pub fn pos(&self) -> u64 {
        match self {
            VarKind::Snp { pos, .. }
            | VarKind::Del { pos, .. }
            | VarKind::Ins { pos, .. }
            | VarKind::Inv { pos, .. }
            | VarKind::Dup { pos, .. } => *pos,
        }
    }

    pub fn len(&self) -> u64 {
        match self {
            VarKind::Snp { .. } => 1,
            VarKind::Ins { seq, .. } => seq.len() as u64,
            VarKind::Del { len, .. } | VarKind::Inv { len, .. } | VarKind::Dup { len, .. } => *len,
        }
    }

    pub fn kind_str(&self) -> &'static str {
        match self {
            VarKind::Snp { .. } => "SNP",
            VarKind::Del { .. } => "DEL",
            VarKind::Ins { .. } => "INS",
            VarKind::Inv { .. } => "INV",
            VarKind::Dup { .. } => "DUP",
        }
    }
}

pub struct SimParams {
    pub snp: u64,
    /// small indels, 1..=49 bp
    pub indel: u64,
    /// SV deletions, 50..=10000 bp
    pub del: u64,
    /// SV insertions of random sequence, 50..=5000 bp
    pub ins: u64,
    /// inversions, 100..=10000 bp
    pub inv: u64,
    /// tandem duplications, 100..=5000 bp
    pub dup: u64,
    pub seed: u64,
}

impl Default for SimParams {
    fn default() -> Self {
        SimParams {
            snp: 500,
            indel: 50,
            del: 5,
            ins: 5,
            inv: 2,
            dup: 2,
            seed: 42,
        }
    }
}

pub struct SimResults {
    pub queries: Vec<FastaRecord>,
    /// (reference contig index, variant), sorted by (contig, pos)
    pub variants: Vec<(usize, VarKind)>,
}

pub fn simulate(refs: &[FastaRecord], p: &SimParams) -> SimResults {
    let mut rng = Rng::new(p.seed);
    // weight contigs by usable length (interior minus margins) so a variant
    // can actually be placed on the contig it is drawn for
    let weights: Vec<u64> = refs
        .iter()
        .map(|r| (r.seq.len() as u64).saturating_sub(2 * MARGIN))
        .collect();
    let total: u64 = weights.iter().sum();

    // Build all variant requests first, then place longest-first: dense small
    // variants carved up early would leave no clear span for kilobase-scale SVs.
    #[derive(Clone)]
    enum Req {
        Snp,
        Del(u64),
        Ins(u64),
        Inv(u64),
        Dup(u64),
    }
    let mut reqs: Vec<Req> = Vec::new();
    reqs.resize(p.snp as usize, Req::Snp);
    for _ in 0..p.indel {
        // 1..=49 bp: strictly below the 50 bp symbolic-SV threshold used in the VCF
        let len = 1 + rng.below(49);
        if rng.below(2) == 0 {
            reqs.push(Req::Del(len));
        } else {
            reqs.push(Req::Ins(len));
        }
    }
    for _ in 0..p.del {
        reqs.push(Req::Del(50 + rng.below(9_951)));
    }
    for _ in 0..p.ins {
        reqs.push(Req::Ins(50 + rng.below(4_951)));
    }
    for _ in 0..p.inv {
        reqs.push(Req::Inv(100 + rng.below(9_901)));
    }
    for _ in 0..p.dup {
        reqs.push(Req::Dup(100 + rng.below(4_901)));
    }
    reqs.sort_by_key(|r| {
        let len = match r {
            Req::Snp => 1,
            Req::Del(l) | Req::Ins(l) | Req::Inv(l) | Req::Dup(l) => *l,
        };
        std::cmp::Reverse(len)
    });

    let mut placed: Vec<Vec<(u64, u64)>> = vec![Vec::new(); refs.len()];
    let mut variants: Vec<(usize, VarKind)> = Vec::new();
    if total > 0 {
        for req in reqs {
            let vlen = match req {
                Req::Snp => 1,
                Req::Del(l) | Req::Ins(l) | Req::Inv(l) | Req::Dup(l) => l,
            };
            if let Some((ci, pos)) = pick_pos(&mut rng, refs, &weights, total, vlen, &mut placed) {
                let var = match req {
                    Req::Snp => {
                        let refbase = refs[ci].seq[pos as usize];
                        let code = base_code(refbase).min(3);
                        let alt = b"ACGT"[(code as usize + 1 + rng.below(3) as usize) % 4];
                        VarKind::Snp { pos, alt }
                    }
                    Req::Del(_) => VarKind::Del { pos, len: vlen },
                    Req::Ins(_) => {
                        let seq: Vec<u8> = (0..vlen).map(|_| rng.base()).collect();
                        VarKind::Ins { pos, seq }
                    }
                    Req::Inv(_) => VarKind::Inv { pos, len: vlen },
                    Req::Dup(_) => VarKind::Dup { pos, len: vlen },
                };
                variants.push((ci, var));
            }
        }
    }

    let mut by_contig: Vec<Vec<VarKind>> = vec![Vec::new(); refs.len()];
    for (ci, v) in &variants {
        by_contig[*ci].push(v.clone());
    }
    let mut queries = Vec::with_capacity(refs.len());
    for (ci, rec) in refs.iter().enumerate() {
        let mut edits = std::mem::take(&mut by_contig[ci]);
        edits.sort_by_key(|e| e.pos());
        let seq = apply_edits(&rec.seq, &edits);
        queries.push(FastaRecord {
            name: rec.name.clone(),
            seq,
        });
    }
    variants.sort_by_key(|(ci, v)| (*ci, v.pos()));

    SimResults { queries, variants }
}

/// Pick a random contig (weighted by usable length) and a non-overlapping
/// position for a variant of length `vlen`. Gives up after 100 attempts.
fn pick_pos(
    rng: &mut Rng,
    refs: &[FastaRecord],
    weights: &[u64],
    total_weight: u64,
    vlen: u64,
    placed: &mut [Vec<(u64, u64)>],
) -> Option<(usize, u64)> {
    for _ in 0..100 {
        let mut r = rng.below(total_weight);
        let mut ci = 0usize;
        for (i, wt) in weights.iter().enumerate() {
            if r < *wt {
                ci = i;
                break;
            }
            r -= wt;
        }
        let len = refs[ci].seq.len() as u64;
        if len < 2 * MARGIN + vlen {
            continue;
        }
        let pos = MARGIN + rng.below(len - 2 * MARGIN - vlen + 1);
        let (a, b) = (pos.saturating_sub(GAP), pos + vlen + GAP);
        if placed[ci].iter().any(|&(x, y)| a < y && x < b) {
            continue;
        }
        placed[ci].push((pos, pos + vlen));
        return Some((ci, pos));
    }
    None
}

/// Apply non-overlapping edits (any order) to `seq`.
pub fn apply_edits(seq: &[u8], edits: &[VarKind]) -> Vec<u8> {
    let mut sorted: Vec<&VarKind> = edits.iter().collect();
    sorted.sort_by_key(|e| e.pos());
    let mut out = Vec::with_capacity(seq.len());
    let mut cursor = 0usize;
    for e in sorted {
        match e {
            VarKind::Snp { pos, alt } => {
                let pos = *pos as usize;
                out.extend_from_slice(&seq[cursor..pos]);
                out.push(*alt);
                cursor = pos + 1;
            }
            VarKind::Del { pos, len } => {
                let (pos, len) = (*pos as usize, *len as usize);
                out.extend_from_slice(&seq[cursor..pos]);
                cursor = pos + len;
            }
            VarKind::Ins { pos, seq: iseq } => {
                let pos = *pos as usize;
                out.extend_from_slice(&seq[cursor..pos]);
                out.extend_from_slice(iseq);
                cursor = pos;
            }
            VarKind::Inv { pos, len } => {
                let (pos, len) = (*pos as usize, *len as usize);
                out.extend_from_slice(&seq[cursor..pos]);
                out.extend_from_slice(&revcomp(&seq[pos..pos + len]));
                cursor = pos + len;
            }
            VarKind::Dup { pos, len } => {
                let (pos, len) = (*pos as usize, *len as usize);
                out.extend_from_slice(&seq[cursor..pos + len]);
                out.extend_from_slice(&seq[pos..pos + len]);
                cursor = pos + len;
            }
        }
    }
    out.extend_from_slice(&seq[cursor..]);
    out
}

/// Small variants (<50 bp) are written as explicit REF/ALT records, larger
/// ones as symbolic ALT with SVTYPE/END/SVLEN.
pub fn write_truth_vcf<W: io::Write>(
    w: &mut VcfWriter<W>,
    refs: &[FastaRecord],
    variants: &[(usize, VarKind)],
) -> io::Result<()> {
    for c in refs {
        w.add_contig(&c.name, c.seq.len() as u64)?;
    }
    w.start_body()?;
    for (ci, v) in variants {
        let rec = &refs[*ci];
        let seq = &rec.seq;
        let chrom = rec.name.as_str();
        match v {
            VarKind::Snp { pos, alt } => {
                let p = *pos as usize;
                w.record(chrom, pos + 1, &seq[p..p + 1], &[*alt], ".")?;
            }
            VarKind::Del { pos, len } => {
                if *len < 50 {
                    let (p, l) = (*pos as usize, *len as usize);
                    w.record(chrom, *pos, &seq[p - 1..p + l], &seq[p - 1..p], ".")?;
                } else {
                    let p = *pos as usize;
                    let info = format!("SVTYPE=DEL;END={};SVLEN={}", pos + len, len);
                    w.record(chrom, pos + 1, &seq[p..p + 1], b"<DEL>", &info)?;
                }
            }
            VarKind::Ins { pos, seq: iseq } => {
                let p = *pos as usize;
                if iseq.len() < 50 {
                    let mut alt = Vec::with_capacity(1 + iseq.len());
                    alt.push(seq[p - 1]);
                    alt.extend_from_slice(iseq);
                    w.record(chrom, *pos, &seq[p - 1..p], &alt, ".")?;
                } else {
                    let info = format!("SVTYPE=INS;SVLEN={}", iseq.len());
                    w.record(chrom, *pos, &seq[p - 1..p], b"<INS>", &info)?;
                }
            }
            VarKind::Inv { pos, len } => {
                let p = *pos as usize;
                let info = format!("SVTYPE=INV;END={};SVLEN={}", pos + len, len);
                w.record(chrom, pos + 1, &seq[p..p + 1], b"<INV>", &info)?;
            }
            VarKind::Dup { pos, len } => {
                let p = *pos as usize;
                let info = format!("SVTYPE=DUP;END={};SVLEN={}", pos + len, len);
                w.record(chrom, pos + 1, &seq[p..p + 1], b"<DUP>", &info)?;
            }
        }
    }
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    fn rand_seq(rng: &mut Rng, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
    }

    #[test]
    fn apply_edits_matches_each_variant_kind() {
        let mut rng = Rng::new(11);
        let seq = rand_seq(&mut rng, 10_000);

        assert_eq!(apply_edits(&seq, &[]), seq);

        let alt = if seq[500] == b'A' { b'C' } else { b'A' };
        let q = apply_edits(&seq, &[VarKind::Snp { pos: 500, alt }]);
        assert_eq!(q.len(), seq.len());
        assert_eq!(q[500], alt);
        assert_eq!(&q[501..], &seq[501..]);

        let q = apply_edits(&seq, &[VarKind::Del { pos: 1000, len: 37 }]);
        assert_eq!(q.len(), seq.len() - 37);
        assert_eq!(&q[997..1000], &seq[997..1000]);
        assert_eq!(&q[1000..], &seq[1037..]);

        let ins = b"GGCCGGCC".to_vec();
        let q = apply_edits(
            &seq,
            &[VarKind::Ins {
                pos: 2000,
                seq: ins.clone(),
            }],
        );
        assert_eq!(q.len(), seq.len() + 8);
        assert_eq!(&q[1993..2000], &seq[1993..2000]);
        assert_eq!(&q[2000..2008], &ins[..]);
        assert_eq!(&q[2008..], &seq[2000..]);

        let q = apply_edits(
            &seq,
            &[VarKind::Inv {
                pos: 3000,
                len: 100,
            }],
        );
        assert_eq!(q.len(), seq.len());
        assert_eq!(&q[3000..3100], &revcomp(&seq[3000..3100])[..]);
        assert_eq!(&q[2900..3000], &seq[2900..3000]);

        let q = apply_edits(&seq, &[VarKind::Dup { pos: 4000, len: 60 }]);
        assert_eq!(q.len(), seq.len() + 60);
        assert_eq!(&q[4000..4060], &seq[4000..4060]);
        assert_eq!(&q[4060..4120], &seq[4000..4060]);
    }

    #[test]
    fn simulate_places_exact_counts_without_overlap() {
        let mut rng = Rng::new(21);
        let refs = vec![
            FastaRecord {
                name: "c1".into(),
                seq: rand_seq(&mut rng, 200_000),
            },
            FastaRecord {
                name: "c2".into(),
                seq: rand_seq(&mut rng, 100_000),
            },
        ];
        let p = SimParams {
            snp: 300,
            indel: 30,
            del: 3,
            ins: 3,
            inv: 2,
            dup: 2,
            seed: 9,
        };
        let res = simulate(&refs, &p);
        let n_of = |f: fn(&VarKind) -> bool| res.variants.iter().filter(|(_, v)| f(v)).count();
        assert_eq!(n_of(|v| matches!(v, VarKind::Snp { .. })), 300);
        assert_eq!(
            n_of(|v| matches!(v, VarKind::Del { len, .. } if *len < 50))
                + n_of(|v| matches!(v, VarKind::Ins { seq, .. } if seq.len() < 50)),
            30
        );
        assert_eq!(
            n_of(|v| matches!(v, VarKind::Del { len, .. } if *len >= 50)),
            p.del as usize
        );
        assert_eq!(
            n_of(|v| matches!(v, VarKind::Ins { seq, .. } if seq.len() >= 50)),
            3
        );
        assert_eq!(n_of(|v| matches!(v, VarKind::Inv { .. })), 2);
        assert_eq!(n_of(|v| matches!(v, VarKind::Dup { .. })), 2);

        // variants within a contig must respect the GAP separation
        let mut by_c: Vec<Vec<&VarKind>> = vec![Vec::new(); refs.len()];
        for (ci, v) in &res.variants {
            by_c[*ci].push(v);
        }
        for list in by_c {
            for w in list.windows(2) {
                assert!(
                    w[1].pos() >= w[0].pos() + w[0].len() + GAP,
                    "variants overlap: {:?} then {:?}",
                    (w[0].pos(), w[0].len()),
                    (w[1].pos(), w[1].len())
                );
            }
        }

        // query contigs keep reference names and plausible lengths
        assert_eq!(res.queries.len(), refs.len());
        for (q, r) in res.queries.iter().zip(&refs) {
            assert_eq!(q.name, r.name);
        }
    }

    #[test]
    fn truth_vcf_has_one_record_per_variant() {
        let mut rng = Rng::new(31);
        let refs = vec![FastaRecord {
            name: "c1".into(),
            seq: rand_seq(&mut rng, 50_000),
        }];
        let p = SimParams {
            snp: 40,
            indel: 10,
            del: 1,
            ins: 1,
            inv: 1,
            dup: 1,
            seed: 4,
        };
        let res = simulate(&refs, &p);
        let mut buf: Vec<u8> = Vec::new();
        {
            let mut vw = crate::vcf::VcfWriter::create(&mut buf).unwrap();
            write_truth_vcf(&mut vw, &refs, &res.variants).unwrap();
        }
        let text = String::from_utf8(buf).unwrap();
        assert!(text.starts_with("##fileformat=VCFv4.2"));
        assert!(text.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
        let records = text.lines().filter(|l| !l.starts_with('#')).count();
        assert_eq!(records, res.variants.len());
        assert!(text.contains("SVTYPE=INV"));
        assert!(text.contains("SVTYPE=DUP"));
    }
}
