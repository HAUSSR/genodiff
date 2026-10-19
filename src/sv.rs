//! Structural variant calling from chain topology (M4).
//!
//! Chain junctions are the SV signal:
//!   - minus chain with plus-chain flanks on the same ref contig -> inversion
//!   - same contig, ref gap > query gap                        -> deletion
//!   - same contig, query gap > ref gap                        -> insertion
//!   - same-strand chains on overlapping ref spans             -> duplication
//!   - chains on different ref contigs (both substantial)      -> breakend
//!
//! Inversion flanks are searched globally rather than only among immediate
//! query neighbors, so intervening noise chains cannot hide the topology.
//! Breakpoint precision is bounded by chain ends (a few hundred bp).

use crate::call::{Var, VarRec};
use crate::chain::ChainRecord;

#[derive(Clone, Copy)]
pub struct SvParams {
    pub min_sv_len: u64,
    /// how far (query bp) a flank chain may sit from the minus chain
    pub max_near: u64,
    /// cross-contig junctions need both chains to carry at least this many
    /// anchors: filters paralog/scaffold noise out of BND calls (real
    /// translocations are backed by long chains on both sides)
    pub bnd_min_anchors: u32,
    /// cross-contig junctions also need the two chains to be close in query:
    /// a real translocation junction is tight; spurious scaffold/paralog
    /// junctions jump around
    pub bnd_max_qgap: u64,
}

impl Default for SvParams {
    fn default() -> Self {
        SvParams {
            min_sv_len: 50,
            max_near: 100_000,
            bnd_min_anchors: 100,
            bnd_max_qgap: 5_000,
        }
    }
}

/// Mark chains mapping into reference bins overloaded by mutually distant
/// query regions (repeat attractors: TE clusters, subtelomeric repeats).
/// Those chains are paralogous, not orthologous, and are excluded from
/// variant calling. A real collinear region maps once per query stretch;
/// five-plus chains from >=50 kb-apart query positions into one 16 kb ref
/// bin is the repeat signature.
pub fn flag_repeat_loci(
    chains: &[ChainRecord],
    bin: u64,
    min_chains: usize,
    spread: u64,
) -> Vec<bool> {
    use std::collections::HashMap;
    let mut groups: HashMap<(u32, u64), Vec<usize>> = HashMap::new();
    for (i, c) in chains.iter().enumerate() {
        let mid = (c.r_start + c.r_end) / 2;
        groups.entry((c.ref_contig, mid / bin)).or_default().push(i);
    }
    let mut flagged = vec![false; chains.len()];
    for idxs in groups.values() {
        if idxs.len() < min_chains {
            continue;
        }
        let qmin = idxs.iter().map(|&i| chains[i].q_start).min().unwrap();
        let qmax = idxs.iter().map(|&i| chains[i].q_start).max().unwrap();
        if qmax.saturating_sub(qmin) >= spread {
            for &i in idxs {
                flagged[i] = true;
            }
        }
    }
    flagged
}

/// True when a query-side and a reference-side gap between two chains are
/// consistent with "the same intervening sequence" (shared repeat or a
/// tight junction). Real inversions flanked by long repeats have symmetric
/// gaps; paralogs do not.
fn gaps_match(qgap: u64, rgap: u64, max_near: u64) -> bool {
    qgap <= max_near
        && rgap <= max_near
        && qgap.abs_diff(rgap) <= qgap.max(rgap) / 4 + 500
}

/// Call SVs from the chain layout of one query contig. `chains` must be in
/// query order (as produced by the chainer); `flagged` marks repeat-attractor
/// chains that are ignored.
pub fn call_svs(
    chains: &[ChainRecord],
    flagged: &[bool],
    p: &SvParams,
    ref_seqs: &[Vec<u8>],
    qseq: &[u8],
) -> Vec<VarRec> {
    use std::collections::HashMap;
    let mut out = Vec::new();

    // ---- inversions: a minus chain whose ref span is abutted by plus
    // chains on the same contig. Flanks are searched across ALL chains of
    // the contig (binary search over sorted endpoints), not just immediate
    // query neighbors, so intervening noise chains cannot hide the topology.
    let mut by_end: HashMap<u32, Vec<(u64, usize)>> = HashMap::new();
    let mut by_start: HashMap<u32, Vec<(u64, usize)>> = HashMap::new();
    for (i, c) in chains.iter().enumerate() {
        if flagged[i] || c.strand != 0 {
            continue;
        }
        by_end.entry(c.ref_contig).or_default().push((c.r_end, i));
        by_start.entry(c.ref_contig).or_default().push((c.r_start, i));
    }
    for v in by_end.values_mut() {
        v.sort_unstable();
    }
    for v in by_start.values_mut() {
        v.sort_unstable();
    }

    for (i, m) in chains.iter().enumerate() {
        if flagged[i] || m.strand != 1 {
            continue;
        }
        let len = m.r_end.saturating_sub(m.r_start);
        if len < p.min_sv_len {
            continue;
        }
        // left flank: the nearest plus chain ending at/before the minus span
        // start, in query order, with a symmetric query/ref gap (whatever
        // separates them — a tight junction or a shared repeat — appears on
        // both axes with the same length)
        let mut left_ok = false;
        if let Some(v) = by_end.get(&m.ref_contig) {
            let j = v.partition_point(|&(re, _)| re <= m.r_start);
            if let Some(ci) = j.checked_sub(1).map(|j| v[j].1) {
                let c = &chains[ci];
                left_ok = c.q_end <= m.q_start
                    && gaps_match(m.q_start - c.q_end, m.r_start - c.r_end, p.max_near);
            }
        }
        // right flank: the nearest plus chain starting at/after the minus
        // span end, in query order, with a symmetric gap
        let mut right_ok = false;
        if let Some(v) = by_start.get(&m.ref_contig) {
            let j = v.partition_point(|&(rs, _)| rs < m.r_end);
            if j < v.len() {
                let c = &chains[v[j].1];
                right_ok = c.q_start >= m.q_end
                    && gaps_match(c.q_start - m.q_end, c.r_start - m.r_end, p.max_near);
            }
        }
        if left_ok && right_ok {
            let s = m.r_start.max(1) as usize;
            let e = m.r_end as usize;
            let seq = ref_seqs[m.ref_contig as usize][s..e].to_vec();
            out.push(VarRec {
                contig: m.ref_contig,
                pos: m.r_start.max(1),
                var: Var::SvInv { seq },
                // the inverted block in the query: same bases, reverse order
                qspan: Some((m.q_start, m.q_end)),
            });
        }
    }

    // ---- pairwise junctions in query order ----
    for w in 0..chains.len().saturating_sub(1) {
        let (a, b) = (&chains[w], &chains[w + 1]);
        if flagged[w] || flagged[w + 1] {
            continue;
        }
        if b.q_start + 64 < a.q_end {
            continue; // query-overlapping pair: suppression leftovers
        }
        if a.ref_contig != b.ref_contig {
            let qgap = b.q_start.saturating_sub(a.q_end);
            if a.n_anchors >= p.bnd_min_anchors
                && b.n_anchors >= p.bnd_min_anchors
                && qgap <= p.bnd_max_qgap
            {
                out.push(VarRec {
                    contig: a.ref_contig,
                    pos: a.r_end.max(1),
                    var: Var::SvTra { chr2: b.ref_contig, pos2: b.r_start.max(1) },
                    qspan: None, // the far side is described by CHR2/POS2
                });
            }
            continue;
        }
        if a.strand != b.strand {
            continue; // orientation mix at a junction: leave to INV pattern
        }
        // duplication: same-strand chains on overlapping ref spans
        let ov = a.r_end.min(b.r_end).saturating_sub(a.r_start.max(b.r_start));
        if ov >= p.min_sv_len {
            let dup_start = a.r_start.max(b.r_start).max(1) as usize;
            let dup_end = a.r_end.min(b.r_end) as usize;
            let seq = ref_seqs[a.ref_contig as usize][dup_start..dup_end].to_vec();
            out.push(VarRec {
                contig: a.ref_contig,
                pos: dup_start as u64,
                var: Var::SvDup { seq },
                // both copies lie between the chains' query extents
                qspan: Some((a.q_start, b.q_end)),
            });
            continue;
        }
        let qgap = b.q_start.saturating_sub(a.q_end);
        // strand-aware ref advance from a to b: plus chains advance forward,
        // minus chains walk the reference backwards
        let rgap = if a.strand == 0 {
            b.r_start as i64 - a.r_end as i64
        } else {
            a.r_start as i64 - b.r_end as i64
        };
        let diff = rgap - qgap as i64;
        // a believable junction is tight on BOTH axes: a small query gap
        // paired with a megabase-scale ref gap (or the mirror image) is the
        // signature of a mis-mapped chain pair / large-scale rearrangement,
        // not a locus-level SV
        if diff >= p.min_sv_len as i64
            && qgap <= p.max_near
            && rgap <= p.max_near as i64
        {
            let start = if a.strand == 0 { a.r_end } else { b.r_end };
            let s = start.max(1) as usize;
            let del_seq = ref_seqs[a.ref_contig as usize][s..(s + diff as usize)].to_vec();
            out.push(VarRec {
                contig: a.ref_contig,
                pos: start.max(1),
                var: Var::SvDel { seq: del_seq },
                qspan: None, // deleted bases exist only in the reference
            });
        } else if -diff >= p.min_sv_len as i64
            && qgap <= p.max_near
            && rgap <= p.max_near as i64
        {
            let x = if a.strand == 0 { b.r_start } else { a.r_start };
            let x_idx = x.max(1) as usize;
            // inserted query bases are the gap between the chains, whatever
            // their strand
            let (qa, qb) = (a.q_end as usize, b.q_start as usize);
            let ins_seq = if qb > qa { qseq[qa..qb].to_vec() } else { Vec::new() };
            if !ins_seq.is_empty() {
                out.push(VarRec {
                    contig: a.ref_contig,
                    pos: x_idx as u64,
                    var: Var::SvIns { seq: ins_seq },
                    // exact query coordinates of the inserted bases
                    qspan: Some((qa as u64, qb as u64)),
                });
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chain(rc: u32, strand: u8, qs: u64, qe: u64, rs: u64, re: u64) -> ChainRecord {
        ChainRecord {
            ref_contig: rc,
            strand,
            q_start: qs,
            q_end: qe,
            r_start: rs,
            r_end: re,
            n_anchors: 100,
            anchors: Vec::new(),
        }
    }

    fn svs(chains: &[ChainRecord]) -> Vec<VarRec> {
        // dummy sequences long enough for all test coordinates
        let ref_seqs = vec![vec![b'A'; 600_000], vec![b'A'; 600_000], vec![b'A'; 600_000],
                            vec![b'A'; 600_000], vec![b'A'; 600_000], vec![b'A'; 600_000]];
        let qseq = vec![b'A'; 600_000];
        call_svs(chains, &vec![false; chains.len()], &SvParams::default(), &ref_seqs, &qseq)
    }

    #[test]
    fn junction_deletion_insertion_translocation() {
        // ref gap 5100 vs query gap 100 -> DEL 5000 at ref 10000
        let s = svs(&[
            chain(0, 0, 0, 10_000, 0, 10_000),
            chain(0, 0, 10_100, 20_000, 15_100, 25_000),
        ]);
        assert_eq!(s.len(), 1);
        assert_eq!(s[0].pos, 10_000);
        assert!(matches!(s[0].var, Var::SvDel { ref seq } if seq.len() == 5_000));
        // query gap 5100 vs ref gap 100 -> INS 5000 before ref 10100
        let s = svs(&[
            chain(0, 0, 0, 10_000, 0, 10_000),
            chain(0, 0, 15_100, 25_000, 10_100, 20_000),
        ]);
        assert_eq!(s.len(), 1);
        assert_eq!(s[0].pos, 10_100);
        assert!(matches!(s[0].var, Var::SvIns { ref seq } if seq.len() == 5_100));

        // different ref contigs, both chains substantial -> BND
        let s = svs(&[
            chain(0, 0, 0, 10_000, 0, 10_000),
            chain(1, 0, 10_100, 20_000, 100_000, 110_000),
        ]);
        assert_eq!(s.len(), 1);
        assert_eq!((s[0].contig, s[0].pos), (0, 10_000));
        assert!(matches!(s[0].var, Var::SvTra { chr2: 1, pos2: 100_000 }));
    }

    #[test]
    fn bnd_requires_substantial_chains_on_both_sides() {
        // a short spurious chain jumping to another contig must NOT become
        // a breakend (paralog noise); the same junction with a substantial
        // chain does
        let mut small = chain(1, 0, 10_100, 12_000, 100_000, 102_000);
        small.n_anchors = 20;
        let s = svs(&[chain(0, 0, 0, 10_000, 0, 10_000), small]);
        assert!(s.is_empty(), "{s:?}");

        let mut big = chain(1, 0, 10_100, 20_000, 100_000, 110_000);
        big.n_anchors = 100;
        let s = svs(&[chain(0, 0, 0, 10_000, 0, 10_000), big]);
        assert_eq!(s.len(), 1);
        assert!(matches!(s[0].var, Var::SvTra { .. }));
    }

    #[test]
    fn bnd_requires_tight_query_junction() {
        // two substantial chains on different contigs but 40kb apart in
        // query: scaffold/paralog noise, not a translocation junction
        let s = svs(&[
            chain(0, 0, 0, 10_000, 0, 10_000),
            chain(1, 0, 50_000, 60_000, 100_000, 110_000),
        ]);
        assert!(s.is_empty(), "{s:?}");
    }

    #[test]
    fn minus_strand_pair_uses_reversed_ref_gap() {
        // minus chains walk the reference backwards: a is ref-higher,
        // b is ref-lower; 5000 ref bases between them missing from query
        let s = svs(&[
            chain(0, 1, 0, 10_000, 20_000, 30_000),
            chain(0, 1, 10_100, 20_000, 5_000, 15_000),
        ]);
        // a.r_start=20000, b.r_end=15000 -> rgap = 5000; qgap = 100
        assert_eq!(s.len(), 1);
        assert_eq!(s[0].pos, 15_000);
        assert!(matches!(s[0].var, Var::SvDel { ref seq } if seq.len() == 4_900));
    }

    #[test]
    fn inversion_from_plus_minus_plus_flanks() {
        let s = svs(&[
            chain(0, 0, 0, 10_000, 0, 10_000),
            chain(0, 1, 10_100, 13_000, 10_100, 13_000),
            chain(0, 0, 13_100, 23_000, 13_100, 23_000),
        ]);
        let invs: Vec<_> = s.iter().filter(|v| matches!(v.var, Var::SvInv { .. })).collect();
        assert_eq!(invs.len(), 1);
        assert_eq!(invs[0].pos, 10_100);
        assert!(matches!(invs[0].var, Var::SvInv { ref seq } if seq.len() == 2_900));
        assert_eq!(invs[0].qspan, Some((10_100, 13_000)));
        assert_eq!(s.len(), 1, "flanking junctions must not add records: {s:?}");
    }

    #[test]
    fn inversion_survives_intervening_noise_chains() {
        // the + - + sandwich is interrupted in QUERY order by a small chain
        // from another contig (paralog noise); the flank search must still
        // see the reference-side topology
        let s = svs(&[
            chain(0, 0, 0, 10_000, 0, 10_000),
            chain(5, 1, 10_100, 10_500, 500_000, 500_400), // noise
            chain(0, 1, 10_600, 13_500, 10_100, 13_000),   // inversion
            chain(0, 0, 13_600, 23_000, 13_100, 22_500),   // right flank
        ]);
        let invs: Vec<_> = s.iter().filter(|v| matches!(v.var, Var::SvInv { .. })).collect();
        assert_eq!(invs.len(), 1, "inversion must be found despite the noise chain");
        assert_eq!(invs[0].pos, 10_100);
        assert!(matches!(invs[0].var, Var::SvInv { ref seq } if seq.len() == 2_900));
    }

    #[test]
    fn inversion_inside_a_larger_plus_chain_is_not_called() {
        // a spurious minus chain fully contained in one big plus chain's ref
        // span has no plus chain ENDING at its start: no INV
        let s = svs(&[
            chain(0, 0, 0, 100_000, 0, 100_000),
            chain(0, 1, 20_000, 22_000, 50_000, 52_000), // paralog, minus
            chain(0, 0, 100_100, 200_000, 100_100, 200_000),
        ]);
        assert!(
            !s.iter().any(|v| matches!(v.var, Var::SvInv { .. })),
            "contained minus chain is not an inversion: {s:?}"
        );
    }

    #[test]
    fn duplication_from_overlapping_ref_spans() {
        let s = svs(&[
            chain(0, 0, 0, 10_000, 1_000, 11_000),
            chain(0, 0, 10_500, 20_500, 1_200, 11_200),
        ]);
        assert_eq!(s.len(), 1);
        assert_eq!(s[0].pos, 1_200);
        assert!(matches!(s[0].var, Var::SvDup { ref seq } if seq.len() == 9_800));
        assert_eq!(s[0].qspan, Some((0, 20_500)));
    }

    #[test]
    fn junction_gaps_must_be_tight_on_both_axes() {
        // small query gap + megabase-scale ref gap: a mis-mapped chain pair
        // (rearrangement signature), not a deletion
        let s = svs(&[
            chain(0, 0, 0, 10_000, 0, 10_000),
            chain(0, 0, 10_100, 20_000, 500_000, 510_000),
        ]);
        assert!(s.is_empty(), "{s:?}");
        // mirror image: small ref gap + megabase-scale query gap is not an
        // insertion
        let s = svs(&[
            chain(0, 0, 0, 10_000, 0, 10_000),
            chain(0, 0, 500_100, 510_000, 10_100, 20_000),
        ]);
        assert!(s.is_empty(), "{s:?}");
    }

    #[test]
    fn minus_strand_insertion_takes_the_inter_chain_gap() {
        // minus chains: the inserted bases are the query gap BETWEEN the
        // chains, not the span from a.q_start across both chains to b.q_end
        let s = svs(&[
            chain(0, 1, 0, 10_000, 20_000, 30_000),
            chain(0, 1, 15_100, 25_000, 5_000, 15_000),
        ]);
        assert_eq!(s.len(), 1);
        assert_eq!(s[0].pos, 20_000);
        assert!(matches!(s[0].var, Var::SvIns { ref seq } if seq.len() == 5_100));
        assert_eq!(s[0].qspan, Some((10_000, 15_100)));
    }

    #[test]
    fn smooth_continuation_emits_nothing() {
        let s = svs(&[
            chain(0, 0, 0, 10_000, 0, 10_000),
            chain(0, 0, 10_050, 20_000, 10_050, 20_000),
            chain(0, 0, 20_050, 30_000, 20_050, 30_000),
        ]);
        assert!(s.is_empty(), "{s:?}");
    }

    #[test]
    fn repeat_flagged_chains_are_ignored() {
        let chains = vec![
            chain(0, 0, 0, 10_000, 0, 10_000),
            chain(0, 0, 10_100, 20_000, 15_100, 25_000),
        ];
        let flagged = vec![false, true];
        let ref_seqs = vec![vec![b'A'; 600_000]];
        let qseq = vec![b'A'; 600_000];
        assert!(call_svs(&chains, &flagged, &SvParams::default(), &ref_seqs, &qseq).is_empty());
    }

    #[test]
    fn repeat_loci_are_flagged_by_scattered_query_hits() {
        // one ref bin hit by 5 chains from query positions 1Mb apart
        let mut chains = Vec::new();
        for i in 0..5u64 {
            chains.push(chain(0, 1, i * 1_000_000, i * 1_000_000 + 3_000, 114_300_000, 114_303_000));
        }
        let flagged = flag_repeat_loci(&chains, 16_384, 5, 50_000);
        assert!(flagged.iter().all(|f| *f));

        // a unique locus hit twice (a real tandem dup) must NOT be flagged
        let chains = vec![
            chain(0, 0, 0, 10_000, 1_000, 11_000),
            chain(0, 0, 10_500, 20_500, 1_200, 11_200),
        ];
        let flagged = flag_repeat_loci(&chains, 16_384, 5, 50_000);
        assert!(flagged.iter().all(|f| !*f));
    }
}
