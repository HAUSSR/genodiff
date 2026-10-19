//! Variant calling from chains (M3+M4): banded alignment of inter-anchor
//! gaps, extraction of SNP/InDel records, left normalization, and SV records
//! from chain topology (deletions/insertions/inversions/duplications/
//! breakends at chain junctions, via the `sv` module).

use crate::align::{align, suggest_band, AlignScores, Op};
use crate::chain::{chain_contig, ChainParams};
use crate::index::RefIndex;
use crate::minimizer::comp;
use crate::sv::{call_svs, flag_repeat_loci, SvParams};

#[derive(Clone)]
pub struct CallParams {
    pub chain: ChainParams,
    pub max_occ: usize,
    /// gaps longer than this on either side are not base-aligned
    pub max_align_len: usize,
    /// for symbolic SV emission from unalignable gaps: the small side must
    /// be at most this many bases
    pub sv_small_side: usize,
    pub min_sv_len: u64,
    /// VCF emission: reference-derived SVs (DEL/INV/DUP) are always
    /// position-described (symbolic allele + END/SVLEN) — their bases are
    /// recoverable from the reference. Only INS carries foreign sequence,
    /// and only up to this many bases; beyond it falls back to symbolic
    /// <INS> with SVLEN (safety valve, roughly the chain-gap bound)
    pub sv_seq_cap: usize,
}

impl Default for CallParams {
    fn default() -> Self {
        CallParams {
            chain: ChainParams::default(),
            max_occ: 100,
            max_align_len: 2000,
            sv_small_side: 100,
            min_sv_len: 50,
            sv_seq_cap: 100_000,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Var {
    Snp { ref_base: u8, alt_base: u8 },
    /// deletion of `len` reference bases starting at 0-based pos-1
    Del { len: u32 },
    /// insertion of `seq` (reference-forward strand) before 0-based pos
    Ins { seq: Vec<u8> },
    SvDel { seq: Vec<u8> },
    SvIns { seq: Vec<u8> },
    SvInv { seq: Vec<u8> },
    SvDup { seq: Vec<u8> },
    /// breakend: the query joins reference contig `chr2` around `pos2`
    SvTra { chr2: u32, pos2: u64 },
}

#[derive(Clone, Debug)]
pub struct VarRec {
    pub contig: u32,
    /// VCF POS (1-based; anchor-base convention for indels)
    pub pos: u64,
    pub var: Var,
    /// query-genome span of this variant (0-based half-open), so SV records
    /// can say exactly where the involved sequence lives in the QUERY
    /// assembly (VCF INFO QSEQ/QSTART/QEND); None = reference-only span
    pub qspan: Option<(u64, u64)>,
}

#[derive(Default, Clone, Copy, Debug)]
pub struct CallStats {
    pub chains_used: u64,
    /// chains skipped by the span-sanity filter (repeat noise)
    pub chains_junk: u64,
    pub gaps_aligned: u64,
    /// one-sided gaps: plain indel, emitted without the aligner
    pub gaps_pure: u64,
    /// seed overlaps inside the chain (consistent, nothing to do)
    pub gaps_overlapping: u64,
    /// indels recovered from unequal seed overlaps (length exact,
    /// position pinned by left normalization)
    pub gaps_slip: u64,
    /// inconsistent anchor geometry (skipped)
    pub gaps_weird: u64,
    pub gaps_sv: u64,
    /// unalignable and not cleanly symbolic (skipped)
    pub gaps_skipped: u64,
    /// mismatch columns with N on either side (skipped)
    pub columns_n: u64,
    pub snps: u64,
    pub dels: u64,
    pub ins: u64,
    pub sv_dels: u64,
    pub sv_ins: u64,
    pub sv_inv: u64,
    pub sv_dup: u64,
    pub sv_tra: u64,
    /// chains excluded as repeat-attractor (paralogous) mappings
    pub chains_repeat: u64,
}

impl CallStats {
    pub fn add(&mut self, o: &CallStats) {
        self.chains_used += o.chains_used;
        self.chains_junk += o.chains_junk;
        self.gaps_aligned += o.gaps_aligned;
        self.gaps_pure += o.gaps_pure;
        self.gaps_overlapping += o.gaps_overlapping;
        self.gaps_slip += o.gaps_slip;
        self.gaps_weird += o.gaps_weird;
        self.gaps_sv += o.gaps_sv;
        self.gaps_skipped += o.gaps_skipped;
        self.columns_n += o.columns_n;
        self.snps += o.snps;
        self.dels += o.dels;
        self.ins += o.ins;
        self.sv_dels += o.sv_dels;
        self.sv_ins += o.sv_ins;
        self.sv_inv += o.sv_inv;
        self.sv_dup += o.sv_dup;
        self.sv_tra += o.sv_tra;
        self.chains_repeat += o.chains_repeat;
    }
}

/// Call variants of one query contig against the indexed reference.
/// `ref_seqs` holds the reference contig sequences in index order.
pub fn call_contig(
    idx: &RefIndex,
    ref_seqs: &[Vec<u8>],
    qseq: &[u8],
    params: &CallParams,
) -> (Vec<VarRec>, CallStats) {
    let (chains, _, _) = chain_contig(idx, qseq, params.max_occ, params.chain, true);
    let mut vars = Vec::new();
    let mut st = CallStats::default();
    let sc = AlignScores::default();
    let k = idx.k as usize;

    // repeat attractors: ref bins hit by many scattered query regions are
    // paralogous mappings (TE clusters etc.) — exclude from calling
    let flagged = flag_repeat_loci(&chains, 16_384, 5, 50_000);
    st.chains_repeat = flagged.iter().filter(|&&f| f).count() as u64;

    // M4: SVs from chain-junction topology (with explicit sequences)
    for rec in call_svs(
        &chains,
        &flagged,
        &SvParams { min_sv_len: params.min_sv_len, ..SvParams::default() },
        ref_seqs,
        qseq,
    ) {
        push_qvar(&mut vars, &mut st, rec.contig, rec.pos, rec.var, rec.qspan);
    }

    for (ci, ch) in chains.iter().enumerate() {
        if flagged[ci] {
            continue;
        }
        // span sanity: a chain whose reference and query spans differ by more
        // than 25% is repeat noise (e.g. query TE copies piling onto one
        // reference locus), not a collinear block worth aligning
        let qspan = ch.q_end - ch.q_start;
        let rspan = ch.r_end - ch.r_start;
        if (rspan as i64 - qspan as i64).unsigned_abs() * 4 > rspan.max(qspan) {
            st.chains_junk += 1;
            continue;
        }
        st.chains_used += 1;
        let rseq = &ref_seqs[ch.ref_contig as usize];

        // oriented anchor coordinates: on the chain's own strand, both os and
        // r increase monotonically along the chain
        let olen = qspan as usize;
        let mut anch: Vec<(usize, usize)> = ch
            .anchors
            .iter()
            .map(|&(q, r)| {
                let os = if ch.strand == 0 {
                    (q as u64 - ch.q_start) as usize
                } else {
                    (ch.q_end - q as u64 - k as u64) as usize
                };
                (os, r as usize)
            })
            .collect();
        anch.sort_unstable();
        debug_assert!(anch.iter().all(|&(os, _)| os + k <= olen));

        // oriented query base at offset os (no revcomp materialization)
        let base = |os: usize| -> u8 {
            if ch.strand == 0 {
                qseq[(ch.q_start + os as u64) as usize]
            } else {
                comp(qseq[(ch.q_end - 1 - os as u64) as usize])
            }
        };
        // left-normalization floor: rolling may go anywhere with pos >= 1
        // (the anchor base at the final pos always exists in the reference)
        let floor = 1usize;

        for w in anch.windows(2) {
            let (os0, r0) = w[0];
            let (os1, r1) = w[1];
            // geometry check: a pair only makes sense on the chain's own
            // diagonal. Anchors salvaged from repeat noise near inversion
            // seams violate this and would emit bogus SNPs.
            let dgap = (r1 as i64 - r0 as i64) - (os1 as i64 - os0 as i64);
            if dgap.unsigned_abs() > params.chain.max_drift as u64 + k as u64 {
                st.gaps_weird += 1;
                continue;
            }
            let (gs, rs) = (os0 + k, r0 + k);
            let (ge, re) = (os1, r1);
            // signed gap lengths: negative = the next seed overlaps this one
            // (minimizers sit ~every (w+1)/2 bases, k-mers are k long, so
            // adjacent seeds usually overlap)
            let gq = ge as i64 - gs as i64;
            let gr = re as i64 - rs as i64;

            if gq >= 0 && gr >= 0 {
                let (gl, rl) = (gq as usize, gr as usize);
                if gl == 0 && rl == 0 {
                    continue; // adjacent anchors
                }
                if gl == 0 {
                    // pure deletion: ref[rs..re] absent from query
                    let (mut d, len) = (rs, rl);
                    normalize_del(rseq, &mut d, len, floor);
                    push_var(&mut vars, &mut st, ch.ref_contig, d as u64, Var::Del { len: len as u32 });
                    st.gaps_pure += 1;
                    continue;
                }
                if rl == 0 {
                    // pure insertion before rs; base() already yields
                    // reference-forward bases on either strand
                    let seq: Vec<u8> = (gs..ge).map(|i| base(i)).collect();
                    let (mut q, mut seq) = (rs, seq);
                    normalize_ins(rseq, &mut q, &mut seq, floor);
                    push_var(&mut vars, &mut st, ch.ref_contig, q as u64, Var::Ins { seq });
                    st.gaps_pure += 1;
                    continue;
                }
                // two-sided gap
                if gl.max(rl) > params.max_align_len {
                    let (mn, mx) = (gl.min(rl), gl.max(rl));
                    if mn <= params.sv_small_side && (mx - mn) as u64 >= params.min_sv_len {
                        if rl > gl {
                            let seq = rseq[rs..rs + (rl - gl)].to_vec();
                            push_var(&mut vars, &mut st, ch.ref_contig, rs as u64, Var::SvDel { seq });
                        } else {
                            let ins_start = gs + rl;
                            let seq: Vec<u8> = (ins_start..ge).map(|i| base(i)).collect();
                            // absolute query coordinates of the inserted bases
                            // (oriented offsets map directly on the plus strand
                            // and mirror on the minus strand)
                            let qspan = if ch.strand == 0 {
                                Some((ch.q_start + ins_start as u64, ch.q_start + ge as u64))
                            } else {
                                Some((ch.q_end - ge as u64, ch.q_end - ins_start as u64))
                            };
                            push_qvar(&mut vars, &mut st, ch.ref_contig, rs as u64,
                                      Var::SvIns { seq }, qspan);
                        }
                        st.gaps_sv += 1;
                    } else {
                        st.gaps_skipped += 1;
                    }
                    continue;
                }
                let oseq: Vec<u8> = (gs..ge).map(|i| base(i)).collect();
                if oseq == rseq[rs..re] {
                    st.gaps_aligned += 1; // identical gap: skip the aligner
                    continue;
                }
                let band = suggest_band(rl, gl);
                let ops = match align(&rseq[rs..re], &oseq, band, &sc) {
                    Some(o) => o,
                    None => match align(&rseq[rs..re], &oseq, band * 2, &sc) {
                        Some(o) => o,
                        None => {
                            st.gaps_skipped += 1;
                            continue;
                        }
                    },
                };
                st.gaps_aligned += 1;

                let mut ri = rs; // absolute reference cursor
                let mut oi = gs; // oriented query cursor
                let mut i = 0;
                while i < ops.len() {
                    match ops[i] {
                        Op::Match => {
                            ri += 1;
                            oi += 1;
                            i += 1;
                        }
                        Op::Mismatch => {
                            let (rb, ob) = (rseq[ri], base(oi));
                            if rb == b'N' || ob == b'N' {
                                st.columns_n += 1;
                            } else {
                                // base() is already reference-forward on both
                                // strands, so ob is the VCF ALT as-is
                                push_var(&mut vars, &mut st, ch.ref_contig, ri as u64 + 1, Var::Snp { ref_base: rb, alt_base: ob });
                            }
                            ri += 1;
                            oi += 1;
                            i += 1;
                        }
                        Op::Del => {
                            let d0 = ri;
                            while i < ops.len() && ops[i] == Op::Del {
                                ri += 1;
                                i += 1;
                            }
                            let (mut d, len) = (d0, ri - d0);
                            normalize_del(rseq, &mut d, len, floor);
                            push_var(&mut vars, &mut st, ch.ref_contig, d as u64, Var::Del { len: len as u32 });
                        }
                        Op::Ins => {
                            let ostart = oi;
                            while i < ops.len() && ops[i] == Op::Ins {
                                oi += 1;
                                i += 1;
                            }
                            let seq: Vec<u8> = (ostart..oi).map(|x| base(x)).collect();
                            let (mut q, mut seq) = (ri, seq);
                            normalize_ins(rseq, &mut q, &mut seq, floor);
                            push_var(&mut vars, &mut st, ch.ref_contig, q as u64, Var::Ins { seq });
                        }
                    }
                }
            } else {
                // overlapping seeds on at least one side
                if gq == gr {
                    st.gaps_overlapping += 1; // consistent overlap, nothing to do
                    continue;
                }
                // unequal overlap (or overlap on one side only): the two seed
                // claims are shifted against each other — an indel of exactly
                // |gq - gr| lives in the overlap. Length is exact; the exact
                // position is not, so left-normalization does the pinning.
                let d = gq - gr;
                let left = (r0 + k).min(r1); // overlap region start, ref side
                if d < 0 {
                    let (mut p, len) = (left, (-d) as usize);
                    normalize_del(rseq, &mut p, len, floor);
                    push_var(&mut vars, &mut st, ch.ref_contig, p as u64, Var::Del { len: len as u32 });
                } else {
                    // query carries d extra bases; take the d oriented bases
                    // just left of the next seed's start as the allele (any
                    // d bases of the overlap are equivalent up to repeats,
                    // and left-normalization pins the position)
                    let qe = ge;
                    let qs = qe.saturating_sub(d as usize);
                    if qe > qs && qe <= olen {
                        let seq: Vec<u8> = (qs..qe).map(|i| base(i)).collect();
                        let (mut p, mut seq) = (left, seq);
                        normalize_ins(rseq, &mut p, &mut seq, floor);
                        push_var(&mut vars, &mut st, ch.ref_contig, p as u64, Var::Ins { seq });
                    } else {
                        st.gaps_skipped += 1;
                    }
                }
                st.gaps_slip += 1;
            }
        }
    }
    (vars, st)
}

fn push_var(vars: &mut Vec<VarRec>, st: &mut CallStats, contig: u32, pos: u64, var: Var) {
    push_qvar(vars, st, contig, pos, var, None);
}

fn push_qvar(
    vars: &mut Vec<VarRec>,
    st: &mut CallStats,
    contig: u32,
    pos: u64,
    var: Var,
    qspan: Option<(u64, u64)>,
) {
    match &var {
        Var::Snp { .. } => st.snps += 1,
        Var::Del { .. } => st.dels += 1,
        Var::Ins { .. } => st.ins += 1,
        Var::SvDel { .. } => st.sv_dels += 1,
        Var::SvIns { .. } => st.sv_ins += 1,
        Var::SvInv { .. } => st.sv_inv += 1,
        Var::SvDup { .. } => st.sv_dup += 1,
        Var::SvTra { .. } => st.sv_tra += 1,
    }
    vars.push(VarRec { contig, pos, var, qspan });
}

/// Shift a deletion of `len` bases at `*d` as far left as the reference
/// allows (VCF left-alignment).
fn normalize_del(rseq: &[u8], d: &mut usize, len: usize, floor: usize) {
    while *d > floor && rseq[*d + len - 1] == rseq[*d - 1] {
        *d -= 1;
    }
}

/// Shift an insertion of `seq` before `*q` as far left as the reference
/// allows, rotating the sequence to keep the allele equivalent.
fn normalize_ins(rseq: &[u8], q: &mut usize, seq: &mut Vec<u8>, floor: usize) {
    while *q > floor && seq.last() == Some(&rseq[*q - 1]) {
        seq.rotate_right(1);
        *q -= 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::FastaRecord;
    use crate::minimizer::revcomp;
    use crate::rng::Rng;
    use crate::sim::{apply_edits, simulate, SimParams, VarKind};

    fn rand_seq(rng: &mut Rng, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
    }

    fn build_idx(refs: &[FastaRecord]) -> RefIndex {
        RefIndex::build(refs.iter().map(|r| Ok(r.clone())), 19, 19).unwrap()
    }

    fn ref_seqs(refs: &[FastaRecord]) -> Vec<Vec<u8>> {
        refs.iter().map(|r| r.seq.clone()).collect()
    }

    /// truth variants, left-normalized against the reference, as comparable tuples
    fn expected_set(refs: &[FastaRecord], sim: &crate::sim::SimResults) -> Vec<(u32, u64, Var)> {
        let mut out = Vec::new();
        for (ci, v) in &sim.variants {
            let rseq = &refs[*ci].seq;
            match v {
                VarKind::Snp { pos, alt } => out.push((
                    *ci as u32,
                    pos + 1,
                    Var::Snp { ref_base: rseq[*pos as usize], alt_base: *alt },
                )),
                VarKind::Del { pos, len } => {
                    let (mut d, len) = (*pos as usize, *len as usize);
                    normalize_del(rseq, &mut d, len, 1);
                    out.push((*ci as u32, d as u64, Var::Del { len: len as u32 }));
                }
                VarKind::Ins { pos, seq } => {
                    let (mut q, mut s) = (*pos as usize, seq.clone());
                    normalize_ins(rseq, &mut q, &mut s, 1);
                    out.push((*ci as u32, q as u64, Var::Ins { seq: s }));
                }
                _ => {} // SVs are handled by the sv_recall test
            }
        }
        out
    }

    fn called_set(vars: &[VarRec]) -> Vec<(u32, u64, Var)> {
        vars.iter().map(|v| (v.contig, v.pos, v.var.clone())).collect()
    }

    #[test]
    fn snp_and_small_indel_recall_precision() {
        let mut rng = Rng::new(20);
        let refs = vec![FastaRecord { name: "c1".into(), seq: rand_seq(&mut rng, 200_000) }];
        let sim = simulate(
            &refs,
            &SimParams { snp: 150, indel: 40, del: 0, ins: 0, inv: 0, dup: 0, seed: 5 },
        );
        let idx = build_idx(&refs);
        let (vars, st) = call_contig(&idx, &ref_seqs(&refs), &sim.queries[0].seq, &CallParams::default());

        let expected = expected_set(&refs, &sim);
        let called = called_set(&vars);
        let n_expected_snp = expected.iter().filter(|(_, _, v)| matches!(v, Var::Snp { .. })).count();
        let n_expected_indel = expected.len() - n_expected_snp;
        let hit = called.iter().filter(|c| expected.contains(c)).count();
        let n_called_snp = called.iter().filter(|(_, _, v)| matches!(v, Var::Snp { .. })).count();
        let n_called_indel = called.len() - n_called_snp;

        assert!(st.chains_junk == 0, "junk chains on a clean genome: {}", st.chains_junk);
        assert!(
            hit as u64 * 100 >= expected.len() as u64 * 95,
            "recall too low: {}/{} (stats: snps={} dels={} ins={} skipped={})",
            hit,
            expected.len(),
            st.snps,
            st.dels,
            st.ins,
            st.gaps_skipped
        );
        assert!(
            hit as u64 * 100 >= called.len() as u64 * 95,
            "precision too low: {}/{}",
            hit,
            called.len()
        );
        assert!(n_called_snp as u64 * 100 >= n_expected_snp as u64 * 95);
        assert!(n_called_indel as u64 * 100 >= n_expected_indel as u64 * 85);
    }

    #[test]
    fn snp_inside_inversion_called_on_correct_strand() {
        let mut rng = Rng::new(21);
        let seq = rand_seq(&mut rng, 200_000);
        let refs = vec![FastaRecord { name: "c1".into(), seq: seq.clone() }];
        // invert [100k, 103k): inside the inverted block, flip the query
        // base at query pos 101000; that position mirrors reference pos
        // ipos+ilen-1-(101000-ipos) = 102999
        let (ipos, ilen, qpos) = (100_000usize, 3_000usize, 101_000usize);
        let mut query = apply_edits(&seq, &[VarKind::Inv { pos: ipos as u64, len: ilen as u64 }]);
        let flipped = if query[qpos] == b'A' { b'C' } else { b'A' };
        query[qpos] = flipped;
        let mirror = ipos + ilen - 1 - (qpos - ipos);

        let idx = build_idx(&refs);
        let (vars, _) = call_contig(&idx, &ref_seqs(&refs), &query, &CallParams::default());
        // M4 emits the INV record itself; the SNP we injected must be the
        // only other variant inside the block
        let inside: Vec<_> = vars
            .iter()
            .filter(|v| {
                v.pos > (ipos + 100) as u64
                    && v.pos < (ipos + ilen - 100) as u64
                    && !matches!(v.var, Var::SvInv { .. })
            })
            .collect();
        assert_eq!(inside.len(), 1, "exactly one non-SV variant inside the inversion: {inside:?}");
        assert!(
            vars.iter().any(|v| matches!(v.var, Var::SvInv { .. })),
            "the inversion itself must be called as SvInv"
        );
        match &inside[0].var {
            Var::Snp { ref_base, alt_base } => {
                assert_eq!(inside[0].pos, mirror as u64 + 1, "POS must be the mirrored ref pos");
                assert_eq!(*ref_base, seq[mirror]);
                // VCF alleles are on the reference's forward strand: the
                // inverted query carries comp(flipped) at the mirrored pos
                assert_eq!(*alt_base, crate::minimizer::comp(flipped));
            }
            other => panic!("expected SNP, got {other:?}"),
        }
        // flanks: no variants (they are identical)
        let flank: Vec<_> = vars
            .iter()
            .filter(|v| {
                (v.pos + 10_000 < ipos as u64 || v.pos > (ipos + ilen + 10_000) as u64)
                    && !matches!(v.var, Var::SvInv { .. })
            })
            .collect();
        assert!(flank.is_empty(), "flanks must be variant-free: {flank:?}");
    }

    #[test]
    fn homopolymer_indels_are_left_normalized() {
        // build: unique flanks + a 30x A homopolymer at [5000,5030)
        let mut rng = Rng::new(22);
        let mut seq = rand_seq(&mut rng, 10_000);
        seq[4999] = b'C'; // non-A base right before the run
        for b in seq.iter_mut().take(5030).skip(5000) {
            *b = b'A';
        }
        let refs = vec![FastaRecord { name: "c1".into(), seq: seq.clone() }];

        // deletion of one A in the middle of the run -> POS must normalize
        // to 5000 (the anchor-base form of the leftmost representation)
        let mut q1 = seq.clone();
        q1.remove(5015);
        let idx = build_idx(&refs);
        let (vars, _) = call_contig(&idx, &ref_seqs(&refs), &q1, &CallParams::default());
        let dels: Vec<_> = vars.iter().filter(|v| matches!(v.var, Var::Del { .. })).collect();
        assert_eq!(dels.len(), 1, "dels: {dels:?}");
        assert_eq!(dels[0].pos, 5000);
        assert!(matches!(dels[0].var, Var::Del { len: 1 }));

        // insertion of one A inside the run -> same leftmost position
        let mut q2 = seq.clone();
        q2.insert(5022, b'A');
        let (vars, _) = call_contig(&idx, &ref_seqs(&refs), &q2, &CallParams::default());
        let inss: Vec<_> = vars.iter().filter(|v| matches!(v.var, Var::Ins { .. })).collect();
        assert_eq!(inss.len(), 1, "inss: {inss:?}");
        assert_eq!(inss[0].pos, 5000);
    }

    #[test]
    fn revcomp_genome_calls_complement_snps() {
        // whole-contig reverse complement: only the SNP we inject shows up;
        // VCF alleles are on the reference's forward strand, so the ALT is
        // the complement of the base the revcomp query carries
        let mut rng = Rng::new(23);
        let seq = rand_seq(&mut rng, 100_000);
        let refs = vec![FastaRecord { name: "c1".into(), seq: seq.clone() }];
        let mut query = revcomp(&seq);
        let p = 50_000usize;
        let flipped = if query[p] == b'G' { b'T' } else { b'G' };
        query[p] = flipped;
        let idx = build_idx(&refs);
        let (vars, st) = call_contig(&idx, &ref_seqs(&refs), &query, &CallParams::default());
        assert_eq!(vars.len(), 1, "only the injected SNP should be called");
        // query pos p mirrors ref pos n-1-p
        let ref_pos = seq.len() - 1 - p;
        match &vars[0].var {
            Var::Snp { ref_base, alt_base } => {
                assert_eq!(vars[0].pos, ref_pos as u64 + 1);
                assert_eq!(*ref_base, seq[ref_pos]);
                assert_eq!(*alt_base, crate::minimizer::comp(flipped));
            }
            other => panic!("expected SNP, got {other:?}"),
        }
        assert_eq!(st.chains_used, 1);
    }

    #[test]
    fn sv_recall_on_simulated_genome() {
        let mut rng = Rng::new(24);
        let refs = vec![
            FastaRecord { name: "c1".into(), seq: rand_seq(&mut rng, 1_200_000) },
            FastaRecord { name: "c2".into(), seq: rand_seq(&mut rng, 800_000) },
        ];
        let sim = simulate(
            &refs,
            &SimParams { snp: 300, indel: 30, del: 8, ins: 8, inv: 4, dup: 4, seed: 33 },
        );
        let idx = build_idx(&refs);
        let mut vars = Vec::new();
        for q in &sim.queries {
            let (v, _) = call_contig(&idx, &ref_seqs(&refs), &q.seq, &CallParams::default());
            vars.extend(v);
        }

        // truth SVs: the simulator's >=50bp classes
        let truth: Vec<(usize, &VarKind)> = sim
            .variants
            .iter()
            .filter(|(_, v)| match v {
                VarKind::Del { len, .. } => *len >= 50,
                VarKind::Ins { seq, .. } => seq.len() >= 50,
                VarKind::Inv { .. } | VarKind::Dup { .. } => true,
                _ => false,
            })
            .map(|(ci, v)| (*ci, v))
            .collect();
        assert_eq!(truth.len(), 24, "simulator must place all requested SVs");

        let kind_of = |v: &Var| -> u8 {
            match v {
                Var::SvDel { .. } => 0,
                Var::SvIns { .. } => 1,
                Var::SvInv { .. } => 2,
                Var::SvDup { .. } => 3,
                _ => 4,
            }
        };
        let truth_kind = |v: &VarKind| -> u8 {
            match v {
                VarKind::Del { .. } => 0,
                VarKind::Ins { .. } => 1,
                VarKind::Inv { .. } => 2,
                VarKind::Dup { .. } => 3,
                _ => 4,
            }
        };

        let mut tp = [0u32; 4];
        for (ci, t) in &truth {
            let (tpos, tlen) = (t.pos(), t.len());
            let hit = vars.iter().any(|v| {
                v.contig == *ci as u32
                    && kind_of(&v.var) == truth_kind(t)
                    && (v.pos as i64 - tpos as i64).abs() <= 2_000
                    && match &v.var {
                        // explicit sequences: check length consistency
                        Var::SvDel { seq } | Var::SvInv { seq } => {
                            tlen <= seq.len() as u64 + 2_000 && (seq.len() as u64) <= tlen + 2_000
                        }
                        Var::SvIns { seq } | Var::SvDup { seq } => {
                            tlen <= seq.len() as u64 + 2_000 && (seq.len() as u64) <= tlen + 2_000
                        }
                        _ => false,
                    }
            });
            tp[truth_kind(t) as usize] += hit as u32;
        }
        assert!(tp[0] >= 7, "DEL recall {}/8", tp[0]);
        assert!(tp[1] >= 7, "INS recall {}/8", tp[1]);
        assert!(tp[2] >= 3, "INV recall {}/4", tp[2]);
        assert!(tp[3] >= 3, "DUP recall {}/4", tp[3]);

        // false positives: called SVs far from any truth SV of any class
        let sv_calls: Vec<&VarRec> = vars
            .iter()
            .filter(|v| (0..4).contains(&kind_of(&v.var)))
            .collect();
        let fp = sv_calls
            .iter()
            .filter(|v| {
                !truth
                    .iter()
                    .any(|(ci, t)| *ci as u32 == v.contig && (v.pos as i64 - t.pos() as i64).abs() <= 5_000)
            })
            .count();
        assert!(fp <= 4, "{fp} SV false positives");
    }
}
