//! Seed anchors: query minimizers looked up in the reference index.
//! This is the input stage for chaining — every anchor is a candidate
//! point on a collinear chain.

use crate::index::RefIndex;
use crate::minimizer::stream_minimizers;

#[derive(Clone, Copy)]
pub struct Anchor {
    pub ref_contig: u32,
    pub ref_pos: u64,
    pub query_pos: u64,
    /// 0 = query matches reference forward strand, 1 = reverse.
    /// ref and query minimizers share a canonical hash; both record whether
    /// they had to flip, and a match is reverse-strand iff exactly one did.
    pub strand: u8,
}

#[derive(Default)]
pub struct SeedStats {
    pub query_minimizers: u64,
    /// minimizers with at least one reference hit (regardless of --max-occ)
    pub matched_minimizers: u64,
    /// matched minimizers dropped because they occur too often in the reference
    pub repeats_skipped: u64,
    pub anchors: u64,
}

impl SeedStats {
    pub fn add(&mut self, o: &SeedStats) {
        self.query_minimizers += o.query_minimizers;
        self.matched_minimizers += o.matched_minimizers;
        self.repeats_skipped += o.repeats_skipped;
        self.anchors += o.anchors;
    }

    pub fn matched_fraction(&self) -> f64 {
        if self.query_minimizers == 0 {
            0.0
        } else {
            self.matched_minimizers as f64 / self.query_minimizers as f64
        }
    }
}

/// Collect anchors for one chunk of query sequence.
///
/// `skip` ignores k-mer starts before this local position: when a contig
/// is split into overlapping chunks, minimizers in the first `w-1` positions
/// of a chunk were already decided by the previous chunk, so they are
/// neither counted nor looked up (keeps stats and output exact).
///
/// `max_occ` = 0 means unlimited; otherwise hashes occurring more than
/// max_occ times in the reference are skipped (repeat control).
pub fn collect_anchors(
    idx: &RefIndex,
    qseq: &[u8],
    skip: u64,
    max_occ: usize,
    stats: &mut SeedStats,
) -> Vec<Anchor> {
    let mut anchors = Vec::new();
    stream_minimizers(qseq, idx.k, idx.w, |h, qpos, qstrand| {
        if qpos < skip {
            return;
        }
        stats.query_minimizers += 1;
        let hits = idx.lookup(h);
        if !hits.is_empty() {
            stats.matched_minimizers += 1;
        }
        if max_occ == 0 || hits.len() <= max_occ {
            stats.anchors += hits.len() as u64;
            for hit in hits {
                anchors.push(Anchor {
                    ref_contig: hit.contig,
                    ref_pos: hit.pos,
                    query_pos: qpos,
                    strand: hit.strand ^ qstrand,
                });
            }
        }
    });
    anchors.sort_by_key(|a| (a.ref_contig, a.ref_pos, a.query_pos));
    stats.anchors += anchors.len() as u64;
    anchors
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::FastaRecord;
    use crate::minimizer::revcomp;
    use crate::rng::Rng;
    use crate::sim::{simulate, SimParams, VarKind};

    fn rand_seq(rng: &mut Rng, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
    }

    fn build_test_index(refs: &[FastaRecord]) -> RefIndex {
        RefIndex::build(refs.iter().map(|r| Ok(r.clone())), 15, 10).unwrap()
    }

    #[test]
    fn snp_only_query_keeps_most_minimizers_matched() {
        let mut rng = Rng::new(2);
        let refs = vec![FastaRecord {
            name: "chr1".into(),
            seq: rand_seq(&mut rng, 100_000),
        }];
        let p = SimParams {
            snp: 100,
            indel: 0,
            del: 0,
            ins: 0,
            inv: 0,
            dup: 0,
            seed: 5,
        };
        let sim = simulate(&refs, &p);
        let idx = build_test_index(&refs);
        let mut stats = SeedStats::default();
        let anchors = collect_anchors(&idx, &sim.queries[0].seq, 0, 1000, &mut stats);
        assert!(!anchors.is_empty());
        assert!(
            stats.matched_fraction() >= 0.95,
            "matched fraction too low: {:.3}",
            stats.matched_fraction()
        );
    }

    #[test]
    fn reverse_complemented_query_hits_opposite_strand() {
        let mut rng = Rng::new(6);
        let refs = vec![FastaRecord {
            name: "chr1".into(),
            seq: rand_seq(&mut rng, 100_000),
        }];
        let p = SimParams {
            snp: 50,
            indel: 0,
            del: 0,
            ins: 0,
            inv: 0,
            dup: 0,
            seed: 5,
        };
        let sim = simulate(&refs, &p);
        let rcq = revcomp(&sim.queries[0].seq);
        let idx = build_test_index(&refs);
        let mut stats = SeedStats::default();
        let anchors = collect_anchors(&idx, &rcq, 0, 1000, &mut stats);
        assert!(stats.matched_fraction() >= 0.95);
        assert!(
            !anchors.is_empty(),
            "revcomp query must still hit the index"
        );
        // Nearly all anchors are reverse-strand. A handful of forward seeds are
        // expected: a random sequence contains some k-mer pairs that are each
        // other's reverse complement, which legitimately match forward.
        let rev = anchors.iter().filter(|a| a.strand == 1).count();
        assert!(
            rev * 100 >= anchors.len() * 99,
            "reverse-strand anchors: {}/{}",
            rev,
            anchors.len()
        );
    }

    #[test]
    fn inverted_region_flips_strand_of_local_anchors() {
        let mut rng = Rng::new(8);
        let refs = vec![FastaRecord {
            name: "chr1".into(),
            seq: rand_seq(&mut rng, 200_000),
        }];
        let p = SimParams {
            snp: 0,
            indel: 0,
            del: 0,
            ins: 0,
            inv: 1,
            dup: 0,
            seed: 12,
        };
        let sim = simulate(&refs, &p);
        let (_, inv) = sim
            .variants
            .iter()
            .find(|(_, v)| matches!(v, VarKind::Inv { .. }))
            .expect("one inversion");
        let (ipos, ilen) = (inv.pos(), inv.len());

        let idx = build_test_index(&refs);
        let mut stats = SeedStats::default();
        let anchors = collect_anchors(&idx, &sim.queries[0].seq, 0, 1000, &mut stats);

        // strictly inside the inversion: anchors must map to the reverse strand
        let inside: Vec<_> = anchors
            .iter()
            .filter(|a| a.query_pos > ipos + 100 && a.query_pos < ipos + ilen - 100)
            .collect();
        assert!(!inside.is_empty(), "inverted region must produce anchors");
        // a few seeds legitimately match the opposite strand where the
        // reference happens to contain a k-mer's reverse complement elsewhere,
        // so require an overwhelming majority rather than every anchor
        let ins_rev = inside.iter().filter(|a| a.strand == 1).count();
        assert!(
            ins_rev * 100 >= inside.len() * 99,
            "reverse-strand anchors inside inversion: {}/{}",
            ins_rev,
            inside.len()
        );
        // far outside: forward strand
        let outside: Vec<_> = anchors
            .iter()
            .filter(|a| a.query_pos + 15_000 < ipos || a.query_pos > ipos + ilen + 15_000)
            .collect();
        assert!(!outside.is_empty());
        let out_fwd = outside.iter().filter(|a| a.strand == 0).count();
        assert!(
            out_fwd * 100 >= outside.len() * 99,
            "forward-strand anchors outside inversion: {}/{}",
            out_fwd,
            outside.len()
        );
    }
}
