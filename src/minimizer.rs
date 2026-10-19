//! Minimizer extraction with canonical (strand-neutral) k-mer hashing.
//!
//! k-mers are packed 2-bit; forward and reverse-complement values are
//! computed incrementally and the smaller one is canonical, so a sequence
//! and its reverse complement yield the same minimizer hashes.

pub const INVALID_CODE: u8 = 255;

/// 0..=3 for ACGT (either case), INVALID_CODE otherwise.
#[inline]
pub fn base_code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => INVALID_CODE,
    }
}

/// splitmix64 finalizer: spreads packed k-mer bits over the full u64 range.
#[inline]
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// complement of a single base (A<->T, C<->G); N passes through
#[inline]
pub fn comp(b: u8) -> u8 {
    match base_code(b) {
        0 => b'T',
        1 => b'G',
        2 => b'C',
        3 => b'A',
        _ => b'N',
    }
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    // complement in 2-bit code is 3 - code: A(0)<->T(3), C(1)<->G(2)
    seq.iter().rev().map(|&b| comp(b)).collect()
}

struct RollingKmer {
    k: u32,
    fwd: u64,
    rev: u64,
    mask: u64,
}

impl RollingKmer {
    fn new(k: u32) -> Self {
        let mask = if k >= 32 {
            !0u64
        } else {
            (1u64 << (2 * k)) - 1
        };
        RollingKmer { k, fwd: 0, rev: 0, mask }
    }

    /// Push the next base (code 0..=3). Callers feed a dummy 0 at N positions;
    /// validity is tracked separately so dummy hashes are never emitted.
    #[inline]
    fn push(&mut self, code: u8) {
        let c = code as u64;
        self.fwd = ((self.fwd << 2) | c) & self.mask;
        self.rev = (self.rev >> 2) | ((3 - c) << (2 * (self.k - 1)));
    }

    #[inline]
    fn hash(&self) -> u64 {
        mix64(self.fwd.min(self.rev))
    }

    /// 1 when the canonical orientation is the reverse complement.
    #[inline]
    fn strand(&self) -> u8 {
        (self.fwd > self.rev) as u8
    }
}

/// Stream (hash, kmer_start, strand) minimizers of `seq`.
/// Guarantees every window of w consecutive k-mers contains an emitted
/// minimizer; consecutive windows sharing the same minimizer emit it once.
pub fn stream_minimizers(seq: &[u8], k: u32, w: u32, mut emit: impl FnMut(u64, u64, u8)) {
    assert!((1..=32).contains(&k), "k must be in 1..=32");
    assert!(w >= 1, "w must be >= 1");
    let (k, w) = (k as usize, w as usize);
    let n = seq.len();
    if n < k + w - 1 {
        return;
    }
    let nk = n - k + 1;

    let mut rk = RollingKmer::new(k as u32);
    let mut ring_h = vec![u64::MAX; w];
    let mut ring_p = vec![0usize; w];
    let mut ring_s = vec![0u8; w];
    // k-mer start s is valid only once s >= valid_from; an N at position j
    // invalidates every start in (j-k, j], i.e. the next valid start is j+1.
    let mut valid_from = 0usize;
    let mut last_emitted: isize = -1;

    for i in 0..n {
        let code = base_code(seq[i]);
        if code == INVALID_CODE {
            rk.push(0);
            valid_from = i + 1;
        } else {
            rk.push(code);
        }
        if i + 1 < k {
            continue; // rolling state not yet a complete k-mer: hash/strand meaningless
        }
        let s = i + 1 - k; // k-mer start; rolling state now holds exactly this k-mer
        let slot = s % w;
        let h = if s >= valid_from { rk.hash() } else { u64::MAX };
        ring_h[slot] = h;
        ring_p[slot] = s;
        ring_s[slot] = rk.strand();

        if s + 1 >= w {
            // window covers k-mer starts [s+1-w, s]; the ring holds exactly them
            let mut bi = 0usize;
            let mut bh = u64::MAX;
            for j in 0..w {
                if ring_h[j] < bh || (ring_h[j] == bh && ring_p[j] < ring_p[bi]) {
                    bh = ring_h[j];
                    bi = j;
                }
            }
            if bh != u64::MAX {
                let p = ring_p[bi];
                if p as isize != last_emitted {
                    emit(bh, p as u64, ring_s[bi]);
                    last_emitted = p as isize;
                }
            }
        }
        if s + 1 == nk {
            break; // last k-mer processed; further bases would pollute the ring
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;
    use std::collections::HashSet;

    /// Brute-force canonical hash of the k-mer starting at `p`,
    /// u64::MAX if it contains a non-ACGT base. Must mirror RollingKmer:
    /// fwd packs the k-mer big-endian, rc packs its reverse complement
    /// (read back-to-front, complemented).
    fn brute_hash(seq: &[u8], k: usize, p: usize) -> u64 {
        let mut fwd = 0u64;
        let mut rc = 0u64;
        for j in 0..k {
            let c = base_code(seq[p + j]);
            if c == INVALID_CODE {
                return u64::MAX;
            }
            fwd = fwd * 4 + c as u64;
        }
        for j in (0..k).rev() {
            rc = rc * 4 + (3 - base_code(seq[p + j]) as u64);
        }
        mix64(fwd.min(rc))
    }

    fn rand_seq(rng: &mut Rng, n: usize, n_rate: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                if n_rate > 0 && rng.below(n_rate) == 0 {
                    b'N'
                } else {
                    b"ACGT"[rng.below(4) as usize]
                }
            })
            .collect()
    }

    #[test]
    fn window_property_holds_on_random_sequences() {
        let mut rng = Rng::new(7);
        for trial in 0..30 {
            let n = (200 + rng.below(400) as usize).max(60);
            let seq = rand_seq(&mut rng, n, if trial % 3 == 0 { 25 } else { 0 });
            let (k, w) = (15u32, 10u32);
            let mut emitted: HashSet<(u64, u64)> = HashSet::new();
            let mut positions: Vec<u64> = Vec::new();
            stream_minimizers(&seq, k, w, |h, p, _| {
                emitted.insert((h, p));
                positions.push(p);
            });
            // emitted positions are strictly increasing (dedup guarantee)
            for w2 in positions.windows(2) {
                assert!(w2[0] < w2[1], "positions not strictly increasing");
            }

            let nk = seq.len() - k as usize + 1;
            if nk < w as usize {
                assert!(positions.is_empty());
                continue;
            }
            for ws in 0..=(nk - w as usize) {
                // argmin (hash, position) over the window, brute force
                let mut best = (u64::MAX, usize::MAX);
                for j in ws..ws + w as usize {
                    let cand = (brute_hash(&seq, k as usize, j), j);
                    if cand.0 < best.0 || (cand.0 == best.0 && cand.1 < best.1) {
                        best = cand;
                    }
                }
                if best.0 != u64::MAX {
                    assert!(
                        emitted.contains(&(best.0, best.1 as u64)),
                        "trial {trial}: window {ws} minimizer at {} not emitted",
                        best.1
                    );
                }
            }
        }
    }

    #[test]
    fn revcomp_yields_identical_minimizer_hashes() {
        let mut rng = Rng::new(3);
        let seq = rand_seq(&mut rng, 2000, 0);
        let rc = revcomp(&seq);
        let mut a = Vec::new();
        let mut b = Vec::new();
        stream_minimizers(&seq, 15, 10, |h, _, _| a.push(h));
        stream_minimizers(&rc, 15, 10, |h, _, _| b.push(h));
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b, "minimizer hashes must be strand-symmetric");
        assert!(!a.is_empty());
    }

    #[test]
    fn revcomp_position_mapping_is_mirrored() {
        let mut rng = Rng::new(5);
        let n = 3000usize;
        let seq = rand_seq(&mut rng, n, 0);
        let rc = revcomp(&seq);
        let k = 15u32;
        let mut fwd = Vec::new();
        let mut rev = Vec::new();
        stream_minimizers(&seq, k, 10, |h, p, s| fwd.push((h, p, s)));
        stream_minimizers(&rc, k, 10, |h, p, s| rev.push((h, p, s)));
        let nk = (n - k as usize + 1) as u64;
        let set: HashSet<(u64, u64, u8)> = rev.into_iter().collect();
        for (h, p, s) in fwd {
            let mirrored = (h, nk - 1 - p, 1 - s);
            assert!(
                set.contains(&mirrored),
                "fwd minimizer ({h},{p},{s}) has no mirror {mirrored:?}"
            );
        }
    }

    #[test]
    fn too_short_sequence_emits_nothing() {
        let seq = b"ACGTACGTA".to_vec(); // 9 bases < k+w-1
        let mut n = 0;
        stream_minimizers(&seq, 15, 10, |_, _, _| n += 1);
        assert_eq!(n, 0);
    }
}
