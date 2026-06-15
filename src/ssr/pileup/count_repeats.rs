//! Fast path — exact motif counting for clean spanning reads (arch §4).
//!
//! The cheap bulk of Stage 1. A read whose tract (the bases between its two
//! clean flank anchors) is a **pure integer tiling** of the locus motif yields
//! one confident on-ladder repeat count plus a base-quality weight in
//! `O(tract length)` — no DP, no allocation. A pure tiling is on-ladder by
//! construction. A tract that is *not* a clean whole-number tiling returns
//! `None`: the read is not fast-path material and falls to the slow-path
//! pair-HMM instead (the gate's "pure tiling" clause, arch §2).
//!
//! This module owns only the tiling check and the weight. The triage-typed
//! wrapper (taking a classified read + its locus, anchoring the tract) lands
//! with `triage.rs`; keeping the core here lets it be tested on raw bytes.

use std::sync::LazyLock;

use crate::ssr::types::Motif;

/// Phred quality → probability the base is **correct**, `1 − 10^(−Q/10)`, as a
/// 256-entry lookup. Process-wide and built once (the per-Q math is identical
/// for every read); mirrors the BAQ engine's `Q2P` table pattern (arch §5.3) —
/// the pattern, not the code. Q0 maps to `0.0` (a base with no confidence), the
/// common Q30 to `0.999`.
static PHRED_CORRECT: LazyLock<[f32; 256]> = LazyLock::new(|| {
    let mut table = [0.0f32; 256];
    for (q, slot) in table.iter_mut().enumerate() {
        *slot = (1.0 - 10f64.powf(-(q as f64) / 10.0)) as f32;
    }
    table
});

#[inline]
fn phred_correct(q: u8) -> f32 {
    PHRED_CORRECT[q as usize]
}

/// The repeat-unit count if `tract` is a pure integer tiling of `motif`, else
/// `None` — a partial trailing unit or any interior base that breaks the tiling
/// (interrupted/impure). The tiling test alone, without the base-quality weight:
/// shared by the fast-path counter ([`count_pure_tiling`]) and the off-ladder
/// degenerate-case check (a tract that *is* a pure tiling is an on-ladder rung,
/// not off-ladder — arch §5.8).
pub(crate) fn pure_tiling_units(tract: &[u8], motif: &Motif) -> Option<u16> {
    let period = motif.period();
    // `Motif::new` guarantees `period >= 1`; the `period == 0` guard (which
    // short-circuits before `is_multiple_of`) keeps a future zero-period path a
    // clean `None` rather than a divide-by-zero.
    if period == 0 || !tract.len().is_multiple_of(period) {
        return None;
    }
    let units = tract.len() / period;
    // Repeat counts are stored as `u16` end-to-end (shared-types §2). A tract
    // long enough to overflow that is not a real spanning-read allele; reject
    // rather than truncate.
    if units > u16::MAX as usize {
        return None;
    }

    // Confirm the pure tiling: every base must match the motif at its phase.
    let unit = motif.as_bytes();
    if tract
        .iter()
        .enumerate()
        .any(|(i, &b)| b != unit[i % period])
    {
        return None;
    }
    Some(units as u16)
}

/// Count the repeat units in `tract` when it is a pure integer tiling of
/// `motif`, returning `(units, weight)`; `None` when the tract is not a clean
/// whole-number tiling (the read then takes the slow path, arch §2/§4).
///
/// `tract` is the read's tract bases and `quals` their matching Phred scores,
/// in read order, so `tract.len() == quals.len()` (a length mismatch is treated
/// as malformed input and returns `None`). The returned `units` is the repeat
/// count `L*`; `weight` is the **mean per-base probability the tract was read
/// correctly** (`1 − 10^(−Q/10)` averaged over the tract), a confidence
/// aggregate Stage 2 uses to down-weight a length whose confident support is all
/// low-quality (spec §4.3) — it is *not* a likelihood. An empty tract is a
/// vacuous clean tiling: `(0, 1.0)`.
pub(crate) fn count_pure_tiling(tract: &[u8], quals: &[u8], motif: &Motif) -> Option<(u16, f32)> {
    if tract.len() != quals.len() {
        return None;
    }
    let units = pure_tiling_units(tract, motif)?;

    let weight = if tract.is_empty() {
        1.0
    } else {
        let sum: f32 = quals.iter().map(|&q| phred_correct(q)).sum();
        sum / tract.len() as f32
    };
    Some((units, weight))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn motif(bytes: &[u8]) -> Motif {
        Motif::new(bytes).unwrap()
    }

    /// A `weight` is within float tolerance of the expected mean prob-correct.
    fn assert_weight_close(got: f32, want: f32) {
        assert!(
            (got - want).abs() < 1e-4,
            "weight {got} not close to expected {want}"
        );
    }

    #[test]
    fn counts_a_pure_dinucleotide_tiling() {
        // CACACA = (CA)×3, all Q30 → weight ≈ 0.999.
        let (units, weight) = count_pure_tiling(b"CACACA", &[30; 6], &motif(b"CA")).unwrap();
        assert_eq!(units, 3);
        assert_weight_close(weight, 0.999);
    }

    #[test]
    fn counts_a_pure_trinucleotide_tiling() {
        let (units, _) = count_pure_tiling(b"CAGCAGCAG", &[40; 9], &motif(b"CAG")).unwrap();
        assert_eq!(units, 3);
    }

    #[test]
    fn counts_a_homopolymer_tiling() {
        let (units, _) = count_pure_tiling(b"AAAAA", &[30; 5], &motif(b"A")).unwrap();
        assert_eq!(units, 5);
    }

    #[test]
    fn rejects_an_interrupted_tract() {
        // One interior base breaks the (CA) tiling → slow path.
        assert_eq!(count_pure_tiling(b"CACGCA", &[30; 6], &motif(b"CA")), None);
    }

    #[test]
    fn rejects_a_partial_trailing_unit() {
        // 5 bases is not a whole number of (CA) units.
        assert_eq!(count_pure_tiling(b"CACAC", &[30; 5], &motif(b"CA")), None);
    }

    #[test]
    fn empty_tract_is_a_vacuous_zero_unit_tiling() {
        assert_eq!(count_pure_tiling(b"", &[], &motif(b"CA")), Some((0, 1.0)));
    }

    #[test]
    fn weight_tracks_base_quality() {
        // Uniform Q10 → mean prob-correct 0.9.
        let (_, weight) = count_pure_tiling(b"CACA", &[10; 4], &motif(b"CA")).unwrap();
        assert_weight_close(weight, 0.9);
    }

    #[test]
    fn weight_is_the_mean_of_mixed_qualities() {
        // Two Q30 (0.999) + two Q10 (0.9) → mean 0.9495.
        let (_, weight) = count_pure_tiling(b"CACA", &[30, 30, 10, 10], &motif(b"CA")).unwrap();
        assert_weight_close(weight, (0.999 + 0.999 + 0.9 + 0.9) / 4.0);
    }

    #[test]
    fn mismatched_tract_and_quals_length_is_rejected() {
        assert_eq!(count_pure_tiling(b"CACA", &[30; 3], &motif(b"CA")), None);
    }
}
