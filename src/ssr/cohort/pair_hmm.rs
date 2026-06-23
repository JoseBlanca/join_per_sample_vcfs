//! The banded forward pair-HMM `align_subst` (arch `ssr_call_genotyping.md` §4,
//! spec §6, implementation plan §3.4, B3).
//!
//! `align_subst(obs | v)` is the probability `P(obs | v)` that an observed
//! repeat-tract sequence arose from one fixed placement variant `v` under a flat
//! per-base error `ε` — the *scoring* counterpart of the simulator's
//! `apply_substitutions` (the plan §B contract). Following HipSTR's repeat-block
//! rule it admits **in-tract substitutions only**; the only length-changing
//! operations (whole-unit slip + placement) are marginalized *outside* this term by
//! the likelihood's `Σ_Δ` / `Σ_v` (B2). Per-base gaps are confined to the **flanks**
//! (extraction boundary slop), never the tract interior — which is what keeps `ε`
//! (composition) and stutter (length) identifiable.
//!
//! Three regimes, fast to slow:
//! - `obs == v` byte-for-byte (the clean post-gate majority) → `(1−ε)^len`, no DP;
//! - equal length, differing bytes → the substitution closed form
//!   `(1−ε)^match · (ε/3)^mismatch`;
//! - differing length → a banded forward summing the few ways the length difference
//!   is absorbed as **flank** gaps (interior gaps are disallowed, so an in-tract
//!   indel is *not* competed — it scores far below a true flank indel).
//!
//! New code reusing the Stage-1 SSR pair-HMM's banded/scratch *pattern* — Stage-1 is
//! Viterbi (max-path); this is a forward **sum** (a probability).

/// Bases of slop tolerated at each end of the tract (extraction boundary). A gap
/// move is admitted only within this many positions of either end; the interior is
/// diagonal-only. Provisional — pinned in F2.
const FLANK_SLOP: usize = 2;

/// Reused forward-DP buffer, so scoring a locus's reads allocates once.
#[derive(Debug, Default)]
pub(crate) struct HmmScratch {
    cells: Vec<f64>,
}

impl HmmScratch {
    pub(crate) fn new() -> Self {
        Self::default()
    }
}

/// `P(obs | variant)` under a flat per-base error `ε`: substitutions in the tract,
/// gaps only in the flanks.
pub(crate) fn align_subst(obs: &[u8], variant: &[u8], eps: f64, scratch: &mut HmmScratch) -> f64 {
    if obs == variant {
        return (1.0 - eps).powi(obs.len() as i32);
    }
    if obs.len() == variant.len() {
        return substitution_product(obs, variant, eps);
    }
    banded_forward(obs, variant, eps, scratch)
}

/// The closed-form substitution score for two equal-length sequences.
fn substitution_product(obs: &[u8], variant: &[u8], eps: f64) -> f64 {
    let mismatch = eps / 3.0;
    obs.iter()
        .zip(variant)
        .map(|(a, b)| if a == b { 1.0 - eps } else { mismatch })
        .product()
}

/// Whether position `pos` (0-based) lies in either flank of a sequence of `len`.
fn in_flank(pos: usize, len: usize) -> bool {
    pos < FLANK_SLOP || pos + FLANK_SLOP >= len
}

/// Banded forward sum for unequal-length sequences, with gaps confined to the
/// flanks. `obs` has length `m`, `variant` length `n`; cells outside the band
/// `|i − j| ≤ |m − n| + FLANK_SLOP` stay zero.
fn banded_forward(obs: &[u8], variant: &[u8], eps: f64, scratch: &mut HmmScratch) -> f64 {
    let m = obs.len();
    let n = variant.len();
    let band = m.abs_diff(n) + FLANK_SLOP;
    let gap = eps; // a flank-slop base is error-like (flat emission); pinned in F2.
    let mismatch = eps / 3.0;
    let emit = |a: u8, b: u8| if a == b { 1.0 - eps } else { mismatch };

    let width = n + 1;
    scratch.cells.clear();
    scratch.cells.resize((m + 1) * width, 0.0);
    let at = |i: usize, j: usize| i * width + j;
    scratch.cells[at(0, 0)] = 1.0;

    for i in 0..=m {
        for j in 0..=n {
            if i == 0 && j == 0 {
                continue;
            }
            if i.abs_diff(j) > band {
                continue;
            }
            let mut acc = 0.0;
            if i > 0 && j > 0 {
                acc += scratch.cells[at(i - 1, j - 1)] * emit(obs[i - 1], variant[j - 1]);
            }
            // Insertion: consume obs[i-1] against a gap in the variant — flank only.
            if i > 0 && in_flank(i - 1, m) {
                acc += scratch.cells[at(i - 1, j)] * gap;
            }
            // Deletion: consume variant[j-1] against a gap in obs — flank only.
            if j > 0 && in_flank(j - 1, n) {
                acc += scratch.cells[at(i, j - 1)] * gap;
            }
            scratch.cells[at(i, j)] = acc;
        }
    }
    scratch.cells[at(m, n)]
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 0.01;

    #[test]
    fn exact_match_is_one_minus_eps_to_the_len() {
        let mut scratch = HmmScratch::new();
        let seq = b"CACACACA";
        let p = align_subst(seq, seq, EPS, &mut scratch);
        assert!((p - (1.0 - EPS).powi(8)).abs() < 1e-15);
    }

    #[test]
    fn one_substitution_scales_by_eps_over_three() {
        let mut scratch = HmmScratch::new();
        let v = b"CACACACA";
        let obs = b"CACAGACA"; // one substitution (C→G) at an interior position
        let got = align_subst(obs, v, EPS, &mut scratch);
        let expected = (1.0 - EPS).powi(7) * (EPS / 3.0);
        assert!(
            (got - expected).abs() < 1e-15,
            "got {got}, expected {expected}"
        );
        // i.e. one substitution costs exactly a factor (ε/3)/(1−ε) vs the exact match.
        let exact = align_subst(v, v, EPS, &mut scratch);
        assert!((got / exact - (EPS / 3.0) / (1.0 - EPS)).abs() < 1e-12);
    }

    #[test]
    fn a_single_flank_indel_is_scored() {
        let mut scratch = HmmScratch::new();
        let v = b"CACACACA";
        let obs = b"CACACACAC"; // one extra base at the trailing flank
        let got = align_subst(obs, v, EPS, &mut scratch);
        // Dominant path: diagonal over all of v, then one flank gap.
        let expected_leading = (1.0 - EPS).powi(8) * EPS;
        assert!(got > 0.0);
        // Within a small factor of the clean flank-gap path (other paths are tiny).
        assert!(
            (got / expected_leading - 1.0).abs() < 0.05,
            "flank-indel score {got} should track the flank-gap path {expected_leading}"
        );
    }

    #[test]
    fn an_interior_indel_is_not_competed_with_a_flank_indel() {
        let mut scratch = HmmScratch::new();
        let v = b"CACACACACACA"; // 12 bp
        let flank_indel = b"CACACACACACAC"; // extra base at the end (flank)
        let mut interior_indel = v.to_vec();
        interior_indel.insert(6, b'G'); // extra base in the middle of the tract
        // sanity: both are one base longer than v
        assert_eq!(flank_indel.len(), v.len() + 1);
        assert_eq!(interior_indel.len(), v.len() + 1);

        let flank = align_subst(flank_indel, v, EPS, &mut scratch);
        let interior = align_subst(&interior_indel, v, EPS, &mut scratch);
        assert!(
            interior < flank * 1e-3,
            "interior indel ({interior}) must score far below a flank indel ({flank})"
        );
    }

    #[test]
    fn longer_obs_returns_a_probability_not_exceeding_one() {
        let mut scratch = HmmScratch::new();
        let v = b"CACACA";
        let obs = b"CACACAC";
        let p = align_subst(obs, v, EPS, &mut scratch);
        assert!((0.0..=1.0).contains(&p), "p = {p}");
    }

    #[test]
    fn scratch_reuse_gives_identical_results() {
        let mut a = HmmScratch::new();
        let mut b = HmmScratch::new();
        let v = b"CACACACA";
        let obs = b"CACACACAC";
        // Prime `a` with an unrelated longer alignment, then reuse it.
        let _ = align_subst(b"CACACACACACA", b"CACA", EPS, &mut a);
        let reused = align_subst(obs, v, EPS, &mut a);
        let fresh = align_subst(obs, v, EPS, &mut b);
        assert_eq!(reused.to_bits(), fresh.to_bits());
    }
}
