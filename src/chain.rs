//! Greedy collinear chaining of seed anchors into alignment blocks (M2).
//!
//! Streams query minimizers in query order and keeps at most one open chain,
//! so memory is O(1) regardless of genome size. Repeat copies that do not
//! continue the current chain are discarded on the fly — this is what keeps
//! the output small on repeat-rich genomes, where raw per-anchor output
//! explodes combinatorially.
//!
//! The greedy heuristic: a chain breaks when the next minimizer has no hit
//! within the drift tolerance, which happens at structural variants and at
//! differentially inserted repeats. Those breaks are exactly the signal
//! later stages (SV calling) consume. Shared repeats do NOT break chains:
//! their coordinates shift together in both genomes, so the drift stays ~0
//! and the chain bridges them.

use crate::index::{Hit, RefIndex};
use crate::minimizer::stream_minimizers;
use crate::seeds::SeedStats;

#[derive(Clone, Copy, Debug)]
pub struct ChainParams {
    /// base tolerance for |Δref − Δquery| between consecutive anchors
    pub max_drift: i64,
    /// hard cap on that tolerance (large query gaps must not chain blindly)
    pub tol_cap: i64,
    /// close a chain when the query jumps more than this between usable hits
    pub max_chain_gap: u64,
    /// chains with fewer anchors are dropped
    pub min_anchors: u32,
    /// chains whose query span is shorter than this are dropped
    pub min_span: u64,
    /// drop chains whose query span is mostly contained in a higher-scoring
    /// chain (repeat-copy noise lands on top of the primary alignment)
    pub suppress_secondary: bool,
    /// keep a chain open through this many query bp of non-extending
    /// minimizers before closing it. Mid-frequency repeats (a handful of
    /// copies in the reference) hit far-away loci and would otherwise kill
    /// an otherwise perfect chain; the true continuation is normally the
    /// next unique minimizer. Real discontinuities (>~max_drift indels)
    /// still break: every post-junction minimizer keeps failing, so the
    /// chain times out after this window.
    pub patience: u64,
}

impl Default for ChainParams {
    fn default() -> Self {
        ChainParams {
            max_drift: 150,
            tol_cap: 2000,
            max_chain_gap: 100_000,
            // 25 anchors at w=19 span ~475 bp: below this, a "chain" is
            // indistinguishable from repeat noise on plant-scale genomes
            min_anchors: 25,
            min_span: 500,
            suppress_secondary: true,
            patience: 1_000,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ChainRecord {
    pub ref_contig: u32,
    /// 0 = query matches reference forward strand, 1 = reverse
    pub strand: u8,
    pub q_start: u64,
    pub q_end: u64,
    /// ref span is always leftmost..rightmost on the forward reference,
    /// regardless of chain strand
    pub r_start: u64,
    pub r_end: u64,
    pub n_anchors: u32,
    /// per-anchor (query_pos, ref_pos), 0-based k-mer starts, in query order.
    /// Only populated when the chainer is created with retain_anchors.
    pub anchors: Vec<(u32, u32)>,
}

#[derive(Default, Clone, Copy, Debug)]
pub struct ChainStats {
    /// chains that passed min_anchors/min_span and survived suppression
    pub chains: u64,
    pub anchors_in_chains: u64,
    /// all closed chains, including ones too short to emit
    pub closed_chains: u64,
    /// closed chains rejected by min_anchors/min_span
    pub short_dropped: u64,
    /// emitted chains dropped by secondary suppression
    pub secondary_dropped: u64,
}

struct OpenChain {
    ref_contig: u32,
    strand: u8,
    last_q: u64,
    last_r: u64,
    q_start: u64,
    q_end: u64,
    r_min: u64,
    r_max: u64,
    n: u32,
    anchors: Vec<(u32, u32)>,
}

pub struct GreedyChainer {
    k: u64,
    params: ChainParams,
    retain_anchors: bool,
    cur: Option<OpenChain>,
    /// first minimizer that failed to extend the current chain (candidate
    /// start for the next one, so chain boundaries stay tight at junctions)
    pending: Option<(u64, u8, Hit)>,
    out: Vec<ChainRecord>,
    stats: ChainStats,
}

impl GreedyChainer {
    pub fn new(k: u32, params: ChainParams, retain_anchors: bool) -> Self {
        GreedyChainer {
            k: k as u64,
            params,
            retain_anchors,
            cur: None,
            pending: None,
            out: Vec::new(),
            stats: ChainStats::default(),
        }
    }

    /// Feed the reference hits of one query minimizer. `hits` must be
    /// non-empty and already occurrence-filtered; minimizers with no hits
    /// are simply not fed (the chain bridges over them until a later
    /// minimizer extends it or the gap cap closes it).
    pub fn feed(&mut self, qpos: u64, qstrand: u8, hits: &[Hit]) {
        let mut extend: Option<usize> = None;
        if let Some(cur) = self.cur.as_ref() {
            let dq = qpos as i64 - cur.last_q as i64;
            if dq > 0 && dq <= self.params.max_chain_gap as i64 {
                // forward chain: ref advances with query;
                // reverse chain: ref position decreases as query advances
                let predicted = if cur.strand == 0 {
                    cur.last_r as i64 + dq
                } else {
                    cur.last_r as i64 - dq
                };
                let tol = self.tol(dq);
                let mut best_d = i64::MAX;
                for (i, h) in hits.iter().enumerate() {
                    if h.contig != cur.ref_contig || (h.strand ^ qstrand) != cur.strand {
                        continue;
                    }
                    let d = (h.pos as i64 - predicted).abs();
                    if d <= tol && d < best_d {
                        best_d = d;
                        extend = Some(i);
                    }
                }
            }
        }
        match extend {
            Some(i) => {
                self.pending = None; // chain continued past the skipped region
                let h = &hits[i];
                let cur = self.cur.as_mut().unwrap();
                if h.pos < cur.r_min {
                    cur.r_min = h.pos;
                }
                if h.pos > cur.r_max {
                    cur.r_max = h.pos;
                }
                cur.last_q = qpos;
                cur.last_r = h.pos;
                cur.q_end = qpos + self.k;
                cur.n += 1;
                if self.retain_anchors {
                    cur.anchors.push((qpos as u32, h.pos as u32));
                }
            }
            None => {
                if let Some(cur) = self.cur.as_ref() {
                    if qpos.saturating_sub(cur.last_q) <= self.params.patience {
                        // within the patience window: skip this minimizer
                        // instead of killing the chain (mid-frequency repeat
                        // copies hit far-away loci; the next unique minimizer
                        // usually extends the chain again)
                        if self.pending.is_none() && qpos > cur.last_q {
                            self.pending = Some((qpos, qstrand, hits[0]));
                        }
                        return;
                    }
                }
                self.flush();
                // start the next chain from the first minimizer that failed
                // against the previous one (tight junction), else this one
                let (pq, ps, ph) = self.pending.take().unwrap_or((qpos, qstrand, hits[0]));
                let mut anchors = Vec::new();
                if self.retain_anchors {
                    anchors.push((pq as u32, ph.pos as u32));
                }
                self.cur = Some(OpenChain {
                    ref_contig: ph.contig,
                    strand: ph.strand ^ ps,
                    last_q: pq,
                    last_r: ph.pos,
                    q_start: pq,
                    q_end: pq + self.k,
                    r_min: ph.pos,
                    r_max: ph.pos,
                    n: 1,
                    anchors,
                });
            }
        }
    }

    fn tol(&self, dq: i64) -> i64 {
        self.params.max_drift.max(dq / 8).min(self.params.tol_cap)
    }

    fn flush(&mut self) {
        if let Some(c) = self.cur.take() {
            self.stats.closed_chains += 1;
            let span = c.q_end.saturating_sub(c.q_start);
            if c.n >= self.params.min_anchors && span >= self.params.min_span {
                self.stats.chains += 1;
                self.stats.anchors_in_chains += c.n as u64;
                self.out.push(ChainRecord {
                    ref_contig: c.ref_contig,
                    strand: c.strand,
                    q_start: c.q_start,
                    q_end: c.q_end,
                    r_start: c.r_min,
                    r_end: c.r_max + self.k,
                    n_anchors: c.n,
                    anchors: if self.retain_anchors {
                        c.anchors
                    } else {
                        Vec::new()
                    },
                });
            }
        }
    }

    pub fn finish(mut self) -> (Vec<ChainRecord>, ChainStats) {
        self.flush();
        (self.out, self.stats)
    }
}

/// Drop chains whose query span overlaps the union of higher-scoring chains'
/// query spans by at least half of their own span. Chains compete by
/// (n_anchors, span); survivors keep their original query order.
pub fn suppress_secondary(chains: Vec<ChainRecord>) -> (Vec<ChainRecord>, u64) {
    let mut order: Vec<usize> = (0..chains.len()).collect();
    order.sort_by_key(|&i| {
        let c = &chains[i];
        std::cmp::Reverse((c.n_anchors, c.q_end.saturating_sub(c.q_start)))
    });
    // merged kept query intervals, sorted and disjoint by construction
    let mut kept: Vec<(u64, u64)> = Vec::new();
    let mut drop = vec![false; chains.len()];
    for &i in &order {
        let c = &chains[i];
        let span = c.q_end.saturating_sub(c.q_start);
        let ov = overlap_with(&kept, c.q_start, c.q_end);
        if span > 0 && ov * 2 >= span {
            drop[i] = true;
        } else {
            insert_interval(&mut kept, c.q_start, c.q_end);
        }
    }
    let n_dropped = drop.iter().filter(|&&d| d).count() as u64;
    (
        chains
            .into_iter()
            .zip(drop)
            .filter(|(_, d)| !*d)
            .map(|(c, _)| c)
            .collect(),
        n_dropped,
    )
}

/// Total bp of [s,e) covered by the sorted, disjoint intervals in `kept`.
fn overlap_with(kept: &[(u64, u64)], s: u64, e: u64) -> u64 {
    let mut ov = 0u64;
    let mut i = kept.partition_point(|&(_, ke)| ke <= s);
    while i < kept.len() {
        let (ks, ke) = kept[i];
        if ks >= e {
            break;
        }
        ov += ke.min(e).saturating_sub(ks.max(s));
        i += 1;
    }
    ov
}

/// Insert [s,e) into the sorted interval list, merging overlaps.
fn insert_interval(kept: &mut Vec<(u64, u64)>, s: u64, e: u64) {
    let lo = kept.partition_point(|&(ks, _)| ks < s);
    let start = if lo > 0 && kept[lo - 1].1 >= s {
        lo - 1
    } else {
        lo
    };
    let mut end = start;
    while end < kept.len() && kept[end].0 < e {
        end += 1;
    }
    let mut ns = s;
    let mut ne = e;
    for &(ks, ke) in &kept[start..end] {
        ns = ns.min(ks);
        ne = ne.max(ke);
    }
    kept.splice(start..end, std::iter::once((ns, ne)));
}

/// Chain one whole query contig against the index. Returns emitted chains,
/// chain stats, and seed stats (minimizer/match counters). With
/// `retain_anchors`, chains carry per-anchor coordinates (needed by the
/// variant caller; map mode skips them to save memory).
pub fn chain_contig(
    idx: &RefIndex,
    seq: &[u8],
    max_occ: usize,
    params: ChainParams,
    retain_anchors: bool,
) -> (Vec<ChainRecord>, ChainStats, SeedStats) {
    let mut stats = SeedStats::default();
    let mut chainer = GreedyChainer::new(idx.k, params, retain_anchors);
    stream_minimizers(seq, idx.k, idx.w, |h, qpos, qstrand| {
        stats.query_minimizers += 1;
        let hits = idx.lookup(h);
        if hits.is_empty() {
            return;
        }
        stats.matched_minimizers += 1;
        if max_occ > 0 && hits.len() > max_occ {
            stats.repeats_skipped += 1;
            return;
        }
        stats.anchors += hits.len() as u64;
        chainer.feed(qpos, qstrand, hits);
    });
    let (chains, cstats) = chainer.finish();
    (chains, cstats, stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::FastaRecord;
    use crate::rng::Rng;
    use crate::sim::{apply_edits, simulate, SimParams, VarKind};

    fn rand_seq(rng: &mut Rng, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
    }

    fn build_idx(refs: &[FastaRecord]) -> RefIndex {
        RefIndex::build(refs.iter().map(|r| Ok(r.clone())), 19, 19).unwrap()
    }

    fn coverage(chains: &[ChainRecord]) -> u64 {
        // greedy chains are sequential in query coordinates: disjoint spans
        chains.iter().map(|c| c.q_end - c.q_start).sum()
    }

    #[test]
    fn snp_only_query_forms_one_long_chain() {
        let mut rng = Rng::new(1);
        let refs = vec![FastaRecord {
            name: "c1".into(),
            seq: rand_seq(&mut rng, 200_000),
        }];
        let sim = simulate(
            &refs,
            &SimParams {
                snp: 200,
                indel: 0,
                del: 0,
                ins: 0,
                inv: 0,
                dup: 0,
                seed: 2,
            },
        );
        let idx = build_idx(&refs);
        let (chains, cstats, sstats) = chain_contig(
            &idx,
            &sim.queries[0].seq,
            100,
            ChainParams::default(),
            false,
        );
        assert!(sstats.matched_fraction() > 0.9);
        assert_eq!(cstats.chains, chains.len() as u64);
        assert!(
            chains.len() <= 5,
            "unexpected chain fragmentation: {}",
            chains.len()
        );
        let cov = coverage(&chains);
        assert!(
            cov * 100 >= 90 * 200_000,
            "chain coverage too low: {cov}/200000 in {} chains",
            chains.len()
        );
        assert!(chains[0].n_anchors > 5_000);
    }

    #[test]
    fn large_deletion_breaks_the_chain() {
        let mut rng = Rng::new(4);
        let seq = rand_seq(&mut rng, 200_000);
        let refs = vec![FastaRecord {
            name: "c1".into(),
            seq: seq.clone(),
        }];
        let query = apply_edits(
            &seq,
            &[VarKind::Del {
                pos: 100_000,
                len: 5_000,
            }],
        );
        let idx = build_idx(&refs);
        let (chains, _, _) = chain_contig(&idx, &query, 100, ChainParams::default(), false);
        assert!(chains.len() >= 2, "a 5kb deletion must break the chain");
        assert!(chains.len() <= 4, "over-fragmented: {}", chains.len());
        let cov = coverage(&chains);
        assert!(cov * 100 >= 95 * 195_000, "coverage {cov}");
        // every emitted chain is forward here
        assert!(chains.iter().all(|c| c.strand == 0));
    }

    #[test]
    fn inversion_forms_reverse_strand_chain() {
        let mut rng = Rng::new(6);
        let seq = rand_seq(&mut rng, 200_000);
        let refs = vec![FastaRecord {
            name: "c1".into(),
            seq: seq.clone(),
        }];
        let (ipos, ilen) = (100_000u64, 3_000u64);
        let query = apply_edits(
            &seq,
            &[VarKind::Inv {
                pos: ipos,
                len: ilen,
            }],
        );
        let idx = build_idx(&refs);
        let (chains, _, _) = chain_contig(&idx, &query, 100, ChainParams::default(), false);

        let inv: Vec<_> = chains
            .iter()
            .filter(|c| c.strand == 1 && c.q_start < ipos + ilen - 100 && c.q_end > ipos + 100)
            .collect();
        assert!(
            !inv.is_empty(),
            "no reverse-strand chain covering the inversion"
        );
        // the inverted block itself chains on the minus strand
        assert!(inv[0].n_anchors > 100, "inversion chain too short");
        // flanking regions stay forward
        assert!(chains
            .iter()
            .any(|c| c.strand == 0 && c.q_end <= ipos + 100));
        assert!(chains
            .iter()
            .any(|c| c.strand == 0 && c.q_start >= ipos + ilen - 100));
    }

    #[test]
    fn tandem_repeat_array_chains_into_one_block() {
        let mut rng = Rng::new(5);
        let te = rand_seq(&mut rng, 5_000);
        let mut ctg2 = Vec::with_capacity(250_000);
        for _ in 0..50 {
            ctg2.extend_from_slice(&te);
        }
        let refs = vec![
            FastaRecord {
                name: "c1".into(),
                seq: rand_seq(&mut rng, 200_000),
            },
            FastaRecord {
                name: "c2".into(),
                seq: ctg2,
            },
        ];
        let idx = build_idx(&refs);
        // self-mapping: every TE minimizer hits 50 reference copies, yet the
        // tandem array must collapse into a single (or few) chain blocks
        let (chains, cstats, sstats) =
            chain_contig(&idx, &refs[1].seq, 100, ChainParams::default(), false);
        assert!(
            sstats.anchors > 1_000_000,
            "expected heavy repeat hit traffic"
        );
        assert!(
            chains.len() <= 3,
            "tandem array fragmented: {}",
            chains.len()
        );
        assert!(cstats.anchors_in_chains > 10_000);
        assert!(coverage(&chains) * 100 >= 99 * 250_000);
        assert!(chains.iter().all(|c| c.n_anchors >= 3));
    }

    #[test]
    fn unrelated_sequence_yields_no_chains() {
        let mut rng = Rng::new(9);
        let refs = vec![FastaRecord {
            name: "c1".into(),
            seq: rand_seq(&mut rng, 100_000),
        }];
        let query = rand_seq(&mut rng, 100_000);
        let idx = build_idx(&refs);
        let (chains, _, _) = chain_contig(&idx, &query, 100, ChainParams::default(), false);
        assert!(chains.is_empty());
    }

    #[test]
    fn secondary_chains_are_suppressed_by_query_overlap() {
        let mk = |qs, qe, n| ChainRecord {
            ref_contig: 0,
            strand: 0,
            q_start: qs,
            q_end: qe,
            r_start: 0,
            r_end: qe - qs,
            n_anchors: n,
            anchors: Vec::new(),
        };
        let chains = vec![
            mk(0, 1_000, 10),      // fully inside the primary -> dropped
            mk(0, 10_000, 100),    // primary: highest anchor count
            mk(5_000, 7_000, 20),  // 100% inside the primary -> dropped
            mk(9_500, 12_000, 30), // 500/2500 = 20% overlap -> kept
            mk(20_000, 21_000, 5), // disjoint -> kept (suppression ignores size)
        ];
        let (kept, nd) = suppress_secondary(chains);
        assert_eq!(nd, 2);
        assert_eq!(kept.len(), 3);
        // survivors stay in query order
        assert_eq!(kept[0].q_start, 0);
        assert_eq!(kept[0].n_anchors, 100);
        assert_eq!(kept[1].q_start, 9_500);
        assert_eq!(kept[2].q_start, 20_000);
    }

    #[test]
    fn inversion_survives_secondary_suppression() {
        // the minus-strand chain of an inversion overlaps the flanking plus
        // chains only at its boundaries; suppression must keep all three
        let mut rng = Rng::new(6);
        let seq = rand_seq(&mut rng, 200_000);
        let refs = vec![FastaRecord {
            name: "c1".into(),
            seq: seq.clone(),
        }];
        let query = apply_edits(
            &seq,
            &[VarKind::Inv {
                pos: 100_000,
                len: 3_000,
            }],
        );
        let idx = build_idx(&refs);
        let (chains, cstats, _) = chain_contig(&idx, &query, 100, ChainParams::default(), false);
        assert_eq!(
            cstats.secondary_dropped, 0,
            "no chain should overlap another by >=50%"
        );
        assert!(chains.iter().any(|c| c.strand == 1));
        assert!(chains.iter().any(|c| c.strand == 0));
    }

    #[test]
    fn chainer_skips_non_extending_minimizers_within_patience() {
        // collinear backbone q=r every 20bp; every 400bp a "repeat" minimizer
        // whose only hit is far away (must not kill the chain)
        let feed_stream = |ch: &mut GreedyChainer| {
            for q in (0..=2000u64).step_by(20) {
                if q % 400 == 200 {
                    ch.feed(
                        q,
                        0,
                        &[Hit {
                            hash: 999,
                            contig: 0,
                            pos: 150_000,
                            strand: 0,
                        }],
                    );
                } else {
                    ch.feed(
                        q,
                        0,
                        &[Hit {
                            hash: q,
                            contig: 0,
                            pos: q,
                            strand: 0,
                        }],
                    );
                }
            }
        };
        let params = ChainParams {
            patience: 1_000,
            min_anchors: 5,
            min_span: 1,
            ..Default::default()
        };
        let mut ch = GreedyChainer::new(19, params, false);
        feed_stream(&mut ch);
        let (chains, _) = ch.finish();
        assert_eq!(chains.len(), 1, "repeat decoys must not break the chain");
        // 101 backbone minimizers, 5 replaced by decoys: all backbone in
        assert_eq!(chains[0].n_anchors, 96);

        // patience=0: the old kill-on-first-failure behavior, must fragment
        let strict = ChainParams {
            patience: 0,
            min_anchors: 5,
            min_span: 1,
            ..Default::default()
        };
        let mut ch2 = GreedyChainer::new(19, strict, false);
        feed_stream(&mut ch2);
        let (chains2, _) = ch2.finish();
        assert!(
            chains2.len() >= 5,
            "without patience every decoy breaks the chain"
        );
    }

    #[test]
    fn large_insertions_break_chains_at_the_junction() {
        // 300bp novel insertions exceed max_drift: each becomes a chain break
        // (later an SV call), with tight junctions from the pending start
        let mut rng = Rng::new(13);
        let base = rand_seq(&mut rng, 200_000);
        let elem = rand_seq(&mut rng, 300);
        let refs = vec![FastaRecord {
            name: "c1".into(),
            seq: base.clone(),
        }];
        // insert a copy every 8 kb across [10k, 170k): 20 insertions
        let mut query: Vec<u8> = Vec::with_capacity(206_000);
        let mut next_ins = 10_000usize;
        for (i, b) in base.iter().enumerate() {
            if i == next_ins && next_ins < 170_000 {
                query.extend_from_slice(&elem);
                next_ins += 8_000;
            }
            query.push(*b);
        }
        let idx = build_idx(&refs);
        let (chains, _, _) = chain_contig(&idx, &query, 100, ChainParams::default(), false);
        assert!(chains.len() >= 20, "each insertion must break the chain");
        assert!(chains.len() <= 22, "over-fragmented: {}", chains.len());
        let cov = coverage(&chains);
        // the chain blocks tile everything except the inserted copies
        assert!(
            cov * 100 >= 95 * query.len() as u64,
            "coverage {cov}/{}",
            query.len()
        );
        assert!(chains.iter().all(|c| c.strand == 0));
    }
}
