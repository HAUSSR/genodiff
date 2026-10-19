//! Banded global alignment with affine gaps (Gotoh) for near-identical
//! sequences — the gap aligner behind variant calling (M3).
//!
//! The band follows the straight line from (0,0) to (m,n), which keeps both
//! extreme gap distributions (all-indel-early / all-indel-late) inside the
//! band as long as band >= |m-n| + slack (see `suggest_band`).

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Op {
    Match,
    Mismatch,
    /// consume a reference base with no query base (deletion)
    Del,
    /// consume a query base with no reference base (insertion)
    Ins,
}

#[derive(Clone, Copy)]
pub struct AlignScores {
    pub match_: i32,
    pub mismatch: i32,
    /// charged once when a gap starts
    pub gap_open: i32,
    /// charged per additional gap base
    pub gap_extend: i32,
}

impl Default for AlignScores {
    fn default() -> Self {
        // 2 mismatches (-12) beat a spurious ins+del pair (-16); a real 1bp
        // indel (-8) beats miscalling it as flanking mismatches (-12)
        AlignScores {
            match_: 2,
            mismatch: -6,
            gap_open: -8,
            gap_extend: -2,
        }
    }
}

const NEG: i32 = i32::MIN / 4;
const UNREACHABLE: i32 = i32::MIN / 8;

/// Suggested half-band for aligning sequences of lengths m and n.
pub fn suggest_band(m: usize, n: usize) -> usize {
    m.abs_diff(n) + 24 + m.min(n) / 6
}

/// Global alignment of `a` (reference) and `b` (query) restricted to a
/// diagonal band. Returns the ops consuming all of a and b, or None when no
/// path fits the band. Empty sides are the caller's business.
pub fn align(a: &[u8], b: &[u8], band: usize, sc: &AlignScores) -> Option<Vec<Op>> {
    let (m, n) = (a.len(), b.len());
    if m == 0 || n == 0 {
        return None;
    }
    let w = 2 * band + 1;
    let mut los = Vec::with_capacity(n + 1);
    let mut his = Vec::with_capacity(n + 1);
    for j in 0..=n {
        let ic = j * m / n;
        los.push(ic.saturating_sub(band));
        his.push((ic + band).min(m));
    }
    // flat per-row storage, column = i - los[j]; states: 0=M, 1=X (del), 2=Y (ins)
    let mut s = vec![[NEG; 3]; w * (n + 1)];
    let mut d = vec![[0u8; 3]; w * (n + 1)];
    for j in 0..=n {
        let (lo, hi) = (los[j], his[j]);
        for i in lo..=hi {
            let c = j * w + (i - lo);
            if i == 0 && j == 0 {
                s[c][0] = 0;
            }
            if j == 0 && i >= 1 {
                s[c][1] = sc.gap_open + sc.gap_extend * (i as i32 - 1);
                d[c][1] = 1;
            }
            if i == 0 && j >= 1 {
                s[c][2] = sc.gap_open + sc.gap_extend * (j as i32 - 1);
                d[c][2] = 2;
            }
        }
    }
    for j in 1..=n {
        let (lo, hi) = (los[j], his[j]);
        let (plo, phi) = (los[j - 1], his[j - 1]);
        for i in lo..=hi {
            if i == 0 {
                continue; // insertion-only border cell already initialized
            }
            let c = j * w + (i - lo);
            // M: from (i-1, j-1)
            if i >= 1 && (i - 1) >= plo && (i - 1) <= phi {
                let pc = (j - 1) * w + ((i - 1) - plo);
                let mut best = NEG;
                let mut dm = 0u8;
                for (st, val) in [(0u8, s[pc][0]), (1, s[pc][1]), (2, s[pc][2])] {
                    if val > best {
                        best = val;
                        dm = st;
                    }
                }
                if best > UNREACHABLE {
                    let sub = if a[i - 1] == b[j - 1] {
                        sc.match_
                    } else {
                        sc.mismatch
                    };
                    s[c][0] = best + sub;
                    d[c][0] = dm;
                }
            }
            // X (deletion): from (i-1, j)
            if i >= lo + 1 {
                let pc = j * w + ((i - 1) - lo);
                let mut bx = UNREACHABLE;
                let mut dx = 0u8;
                if s[pc][0] + sc.gap_open > bx {
                    bx = s[pc][0] + sc.gap_open;
                    dx = 0;
                }
                if s[pc][1] + sc.gap_extend > bx {
                    bx = s[pc][1] + sc.gap_extend;
                    dx = 1;
                }
                if s[pc][2] + sc.gap_open > bx {
                    bx = s[pc][2] + sc.gap_open;
                    dx = 2;
                }
                s[c][1] = bx;
                d[c][1] = dx;
            }
            // Y (insertion): from (i, j-1)
            if i >= plo && i <= phi {
                let pc = (j - 1) * w + (i - plo);
                let mut by = UNREACHABLE;
                let mut dy = 0u8;
                if s[pc][0] + sc.gap_open > by {
                    by = s[pc][0] + sc.gap_open;
                    dy = 0;
                }
                if s[pc][1] + sc.gap_open > by {
                    by = s[pc][1] + sc.gap_open;
                    dy = 1;
                }
                if s[pc][2] + sc.gap_extend > by {
                    by = s[pc][2] + sc.gap_extend;
                    dy = 2;
                }
                s[c][2] = by;
                d[c][2] = dy;
            }
        }
    }
    let lo_n = los[n];
    debug_assert!(m >= lo_n && m <= his[n]);
    let c = n * w + (m - lo_n);
    let (best, mut state) = if s[c][0] >= s[c][1] && s[c][0] >= s[c][2] {
        (s[c][0], 0u8)
    } else if s[c][1] >= s[c][2] {
        (s[c][1], 1)
    } else {
        (s[c][2], 2)
    };
    if best <= UNREACHABLE {
        return None;
    }
    let mut ops = Vec::with_capacity(m + n);
    let (mut i, mut j) = (m, n);
    while i > 0 || j > 0 {
        let cc = j * w + (i - los[j]);
        match state {
            0 => {
                ops.push(if a[i - 1] == b[j - 1] {
                    Op::Match
                } else {
                    Op::Mismatch
                });
                state = d[cc][0];
                i -= 1;
                j -= 1;
            }
            1 => {
                ops.push(Op::Del);
                state = d[cc][1];
                i -= 1;
            }
            _ => {
                ops.push(Op::Ins);
                state = d[cc][2];
                j -= 1;
            }
        }
    }
    ops.reverse();
    Some(ops)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    fn ops_string(ops: &[Op]) -> String {
        ops.iter()
            .map(|o| match o {
                Op::Match => '=',
                Op::Mismatch => '*',
                Op::Del => 'D',
                Op::Ins => 'I',
            })
            .collect()
    }

    #[test]
    fn identical_sequences_all_match() {
        let a = b"ACGTACGTACGTAAAACCCCGGGGTTTT";
        let ops = align(
            a,
            a,
            suggest_band(a.len(), a.len()),
            &AlignScores::default(),
        )
        .unwrap();
        assert!(ops.iter().all(|o| *o == Op::Match));
    }

    #[test]
    fn single_mismatch_is_called() {
        let a = b"ACGTACGTACGTACGTACGT";
        let mut b = a.to_vec();
        b[10] = b'T';
        let ops = align(
            a,
            &b,
            suggest_band(a.len(), b.len()),
            &AlignScores::default(),
        )
        .unwrap();
        assert_eq!(ops_string(&ops), "==========*=========");
    }

    #[test]
    fn deletion_and_insertion_are_called() {
        let a = b"ACGTACGTACGTAAAACGTACGT";
        let mut b = a.to_vec();
        b.remove(12); // 1bp deletion
        let ops = align(
            a,
            &b,
            suggest_band(a.len(), b.len()),
            &AlignScores::default(),
        )
        .unwrap();
        assert_eq!(ops_string(&ops), "============D==========");
        let mut c = a.to_vec();
        c.insert(8, b'G'); // 1bp insertion
        let ops = align(
            a,
            &c,
            suggest_band(a.len(), c.len()),
            &AlignScores::default(),
        )
        .unwrap();
        assert_eq!(ops_string(&ops), "========I===============");
    }

    #[test]
    fn tight_band_returns_none_when_indel_present() {
        // 16bp deletion (ref 20 vs query 4): the band center jumps 5 cells
        // per query row, so band=1 leaves the rows disconnected (no path at
        // all -> None); a wide band aligns and reports the deletion
        let a = b"ACGTACGTACGTACGTACGT";
        let b = b"ACGT";
        assert!(align(a, b, 1, &AlignScores::default()).is_none());
        let ops = align(a, b, 20, &AlignScores::default()).unwrap();
        assert_eq!(ops.len(), 20);
        assert_eq!(ops.iter().filter(|o| **o == Op::Del).count(), 16);
    }

    /// Unbanded reference Gotoh for tiny inputs; with a huge band our aligner
    /// must return the exact same optimal score.
    fn brute_score(a: &[u8], b: &[u8], sc: &AlignScores) -> i32 {
        let (m, n) = (a.len(), b.len());
        let neg = i32::MIN / 4;
        let mut mm = vec![vec![neg; n + 1]; m + 1];
        let mut xx = vec![vec![neg; n + 1]; m + 1];
        let mut yy = vec![vec![neg; n + 1]; m + 1];
        mm[0][0] = 0;
        for i in 1..=m {
            xx[i][0] = sc.gap_open + sc.gap_extend * (i as i32 - 1);
        }
        for j in 1..=n {
            yy[0][j] = sc.gap_open + sc.gap_extend * (j as i32 - 1);
        }
        for i in 1..=m {
            for j in 1..=n {
                let best_prev = mm[i - 1][j - 1].max(xx[i - 1][j - 1]).max(yy[i - 1][j - 1]);
                if best_prev > neg / 2 {
                    mm[i][j] = best_prev
                        + if a[i - 1] == b[j - 1] {
                            sc.match_
                        } else {
                            sc.mismatch
                        };
                }
                xx[i][j] = (mm[i - 1][j] + sc.gap_open)
                    .max(xx[i - 1][j] + sc.gap_extend)
                    .max(yy[i - 1][j] + sc.gap_open);
                yy[i][j] = (mm[i][j - 1] + sc.gap_open)
                    .max(xx[i][j - 1] + sc.gap_open)
                    .max(yy[i][j - 1] + sc.gap_extend);
            }
        }
        mm[m][n].max(xx[m][n]).max(yy[m][n])
    }

    fn score_of(ops: &[Op], a: &[u8], b: &[u8], sc: &AlignScores) -> i32 {
        let (mut i, mut j, mut score) = (0usize, 0usize, 0);
        let mut in_del = false;
        let mut in_ins = false;
        for op in ops {
            match op {
                Op::Match => {
                    score += sc.match_;
                    in_del = false;
                    in_ins = false;
                    i += 1;
                    j += 1;
                }
                Op::Mismatch => {
                    score += sc.mismatch;
                    in_del = false;
                    in_ins = false;
                    i += 1;
                    j += 1;
                }
                Op::Del => {
                    score += if in_del { sc.gap_extend } else { sc.gap_open };
                    in_del = true;
                    in_ins = false;
                    i += 1;
                }
                Op::Ins => {
                    score += if in_ins { sc.gap_extend } else { sc.gap_open };
                    in_ins = true;
                    in_del = false;
                    j += 1;
                }
            }
        }
        assert_eq!(i, a.len(), "ops must consume all of a");
        assert_eq!(j, b.len(), "ops must consume all of b");
        score
    }

    #[test]
    fn matches_unbanded_optimum_on_random_pairs() {
        let mut rng = Rng::new(41);
        let sc = AlignScores::default();
        for _ in 0..200 {
            let m = 1 + rng.below(8) as usize;
            let n = 1 + rng.below(8) as usize;
            let a: Vec<u8> = (0..m).map(|_| b"ACGT"[rng.below(4) as usize]).collect();
            let b: Vec<u8> = (0..n).map(|_| b"ACGT"[rng.below(4) as usize]).collect();
            let huge = m.max(n) + 2; // band never binds
            let ops = align(&a, &b, huge, &sc).unwrap();
            let got = score_of(&ops, &a, &b, &sc);
            let want = brute_score(&a, &b, &sc);
            assert_eq!(got, want, "a={:?} b={:?}", a, b);
        }
    }

    #[test]
    fn suggested_band_suffices_for_diverged_random_pairs() {
        let mut rng = Rng::new(43);
        let sc = AlignScores::default();
        for trial in 0..50 {
            let n = 100 + rng.below(200) as usize;
            let a: Vec<u8> = (0..n).map(|_| b"ACGT"[rng.below(4) as usize]).collect();
            let mut b = a.clone();
            // ~3% divergence: mutations + a few short indels
            for _ in 0..(n / 30) {
                let p = rng.below(b.len() as u64 - 1) as usize;
                match rng.below(3) {
                    0 => b[p] = b"ACGT"[rng.below(4) as usize],
                    1 => {
                        b.remove(p);
                    }
                    _ => {
                        b.insert(p, b"ACGT"[rng.below(4) as usize]);
                    }
                }
            }
            let band = suggest_band(a.len(), b.len());
            let ops = align(&a, &b, band, &sc)
                .unwrap_or_else(|| panic!("trial {trial}: band {} too tight", band));
            // replay validity: every Match column must truly be equal
            let (mut i, mut j) = (0, 0);
            for op in &ops {
                match op {
                    Op::Match => {
                        assert_eq!(a[i], b[j]);
                        i += 1;
                        j += 1;
                    }
                    Op::Mismatch => {
                        i += 1;
                        j += 1;
                    }
                    Op::Del => i += 1,
                    Op::Ins => j += 1,
                }
            }
            assert_eq!((i, j), (a.len(), b.len()));
        }
    }
}
