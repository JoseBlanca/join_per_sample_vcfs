//! Shared read-to-allele attribution primitives.
//!
//! Two length-only stages still hard-attribute each spanning read to the nearest parent
//! allele by repeat-length distance and act on the signed slip `Δ = read_units − parent_units`:
//! the pre-pass confident-genotype stats ([`prepass::accumulate_locus`](super::prepass), a
//! confident *single* genotype so no same-length ambiguity arises). That length-only primitive
//! is [`nearest_parent`].
//!
//! The genotyping stages — the per-locus slip refit ([`em::attribute_locus`](super::em)) and
//! the allele-balance FP term (via the EM's final E-step, [`em`](super::em)) — instead key on
//! the candidate **sequence**, so two same-length called alleles (an interruption
//! polymorphism) are told apart by composition, not collapsed by length (spec §5.3, Mark-2 §6
//! amendment). [`nearest_called_by_sequence`] is the hard version (slip stats);
//! [`allele_responsibilities`] is the soft normalization the balance term consumes.

use smallvec::{SmallVec, smallvec};

use crate::ssr::cohort::pair_hmm::{HmmScratch, align_subst};

/// Hard-attribute a read of `read_units` repeat units to its nearest parent allele among
/// `parent_units`, returning `(index, delta)` where `delta = read_units − parent_units[index]`
/// (signed; `0` is a faithful read). A distance tie resolves to the first (lowest-index)
/// parent — the single shared tie-break so it cannot drift between call sites. `None` iff
/// `parent_units` is empty.
pub(crate) fn nearest_parent(read_units: i32, parent_units: &[u16]) -> Option<(usize, i32)> {
    parent_units
        .iter()
        .enumerate()
        .min_by_key(|&(_, &u)| (u as i32 - read_units).abs())
        .map(|(i, &u)| (i, read_units - u as i32))
}

/// Hard-attribute a read to its nearest **called** allele, breaking a length-distance tie by
/// composition. `called` is `(sequence, length in repeat units)` per called allele; returns
/// `(index into `called`, signed length slip Δ = read_units − called_units[index])`, or `None`
/// iff `called` is empty.
///
/// **Length distance first; composition only to break a *same-length* tie.** A slip *is* a
/// length change, so the parent is the length-nearest allele (and `align_subst` cannot even
/// score a multi-unit slip — the gaps exceed its flank band). The composition tie-break is used
/// **only among length-nearest alleles that share a length** — exactly the same-length case an
/// interruption polymorphism creates, where the `align_subst` substitution score (the metric
/// inside `Qᵣ`, so attribution and the likelihood agree, spec §5.3) picks the composition-
/// matching allele; the slip `Δ` is then identical whichever same-length allele wins, so the
/// `θ_locus` slip statistics are unchanged (effect-neutral). Two *different-length* alleles that
/// are equidistant from the read carry no composition signal comparable across lengths, so they
/// keep the deterministic **lowest-index** rule of [`nearest_parent`] (letting `align_subst`
/// decide there could flip the slip sign). A single nearest allele is taken directly — no
/// alignment is run, so a pure locus / length-separated het pays only integer comparisons.
pub(crate) fn nearest_called_by_sequence(
    obs: &[u8],
    read_units: i32,
    called: &[(&[u8], u16)],
    eps: f64,
    scratch: &mut HmmScratch,
) -> Option<(usize, i32)> {
    let min_dist = called
        .iter()
        .map(|&(_, units)| (i32::from(units) - read_units).abs())
        .min()?;
    // Indices of the length-nearest alleles, lowest first.
    let tied: SmallVec<[usize; 2]> = called
        .iter()
        .enumerate()
        .filter(|&(_, &(_, units))| (i32::from(units) - read_units).abs() == min_dist)
        .map(|(idx, _)| idx)
        .collect();
    let all_same_length = tied.iter().all(|&idx| called[idx].1 == called[tied[0]].1);
    let best_idx = if tied.len() == 1 || !all_same_length {
        // A single nearest, or a *different-length* equidistant tie: take the lowest index —
        // no `align_subst`. Matches `nearest_parent`.
        tied[0]
    } else {
        // ≥ 2 alleles tie at the *same* length → break on composition (`align_subst`); strict
        // `>` keeps the lowest-index maximum. Δ is identical across them, so this is
        // effect-neutral for the slip stats.
        tied.iter()
            .map(|&idx| (idx, align_subst(obs, called[idx].0, eps, scratch)))
            .reduce(|best, cur| if cur.1 > best.1 { cur } else { best })
            .map_or(tied[0], |(idx, _)| idx)
    };
    Some((best_idx, read_units - i32::from(called[best_idx].1)))
}

/// The soft per-allele responsibilities of one read: each called allele's `Qᵣ(obs | a)`
/// normalized over the called alleles (Mark-2 §6 amendment — the input to the allele-balance
/// deconvolution). The `Qᵣ` values are supplied by the caller (the EM, which owns the read
/// model; arch §3 / Q-I1), so this stays a pure, order-free, testable normalization.
///
/// On an all-zero input (a read that explains *none* of the called alleles — random junk the
/// `λ` outlier term handles upstream) it returns an even split so the responsibilities still
/// sum to 1; the balance caller skips such reads rather than crediting them, so this branch is
/// defensive.
pub(crate) fn allele_responsibilities(qr_per_allele: &[f64]) -> SmallVec<[f64; 2]> {
    let total: f64 = qr_per_allele.iter().sum();
    if total <= 0.0 {
        let n = qr_per_allele.len();
        return smallvec![if n == 0 { 0.0 } else { 1.0 / n as f64 }; n];
    }
    qr_per_allele.iter().map(|q| q / total).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn picks_the_nearest_parent_and_signs_the_delta() {
        // 9 units: nearest of {6, 10} is 10 (distance 1 < 3), delta = 9 - 10 = -1.
        assert_eq!(nearest_parent(9, &[6, 10]), Some((1, -1)));
        // 7 units: nearest of {6, 10} is 6, delta = +1.
        assert_eq!(nearest_parent(7, &[6, 10]), Some((0, 1)));
        // Faithful read on a parent → delta 0.
        assert_eq!(nearest_parent(6, &[6, 10]), Some((0, 0)));
    }

    #[test]
    fn breaks_a_distance_tie_toward_the_first_parent() {
        // 8 is equidistant from 6 and 10 → the first (index 0) wins.
        assert_eq!(nearest_parent(8, &[6, 10]), Some((0, 2)));
        // Order matters: the same tie with the parents swapped picks the other length,
        // still at index 0 — the shared "first parent" rule.
        assert_eq!(nearest_parent(8, &[10, 6]), Some((0, -2)));
    }

    #[test]
    fn returns_none_for_no_parents() {
        assert_eq!(nearest_parent(8, &[]), None);
    }

    // ── P1.4: sequence-aware attribution ──

    const EPS: f64 = 0.01;

    fn ca(units: u16) -> Vec<u8> {
        "CA".repeat(units as usize).into_bytes()
    }

    /// `CA×6` (12 bp) and its same-length interruption sibling (one interior base flipped).
    fn pure6() -> Vec<u8> {
        ca(6)
    }
    fn interrupted6() -> Vec<u8> {
        let mut s = pure6();
        s[5] = b'T';
        s
    }

    #[test]
    fn same_length_alleles_attributed_by_composition() {
        let (pure, interrupted) = (pure6(), interrupted6());
        let called: &[(&[u8], u16)] = &[(&pure, 6), (&interrupted, 6)];
        let mut scratch = HmmScratch::new();
        // A read matching the interrupted composition → the interrupted allele (index 1), Δ 0.
        assert_eq!(
            nearest_called_by_sequence(&interrupted, 6, called, EPS, &mut scratch),
            Some((1, 0))
        );
        // A read matching the pure composition → the pure allele (index 0), Δ 0.
        assert_eq!(
            nearest_called_by_sequence(&pure, 6, called, EPS, &mut scratch),
            Some((0, 0))
        );
    }

    #[test]
    fn slip_delta_is_length_based_and_effect_neutral_for_same_length_alleles() {
        let (pure, interrupted) = (pure6(), interrupted6());
        let called: &[(&[u8], u16)] = &[(&pure, 6), (&interrupted, 6)];
        let mut scratch = HmmScratch::new();
        // A +1 slip read (length 7): whichever same-length allele composition picks, Δ = +1.
        let read7 = ca(7);
        let (_, delta) = nearest_called_by_sequence(&read7, 7, called, EPS, &mut scratch).unwrap();
        assert_eq!(
            delta, 1,
            "a length-7 read is a +1 slip of a length-6 allele"
        );
        // A −3 slip read (length 3): align_subst cannot score a 3-unit slip (gaps exceed the
        // flank band), so both score ~0 → the tie falls to the lower index, but the length-
        // based Δ = −3 is still correct — the reason attribution stays length-first.
        let read3 = ca(3);
        assert_eq!(
            nearest_called_by_sequence(&read3, 3, called, EPS, &mut scratch),
            Some((0, -3))
        );
    }

    #[test]
    fn length_separated_alleles_match_the_length_only_primitive() {
        let (a6, b10) = (ca(6), ca(10));
        let called: &[(&[u8], u16)] = &[(&a6, 6), (&b10, 10)];
        let mut scratch = HmmScratch::new();
        // A length-7 read is nearest A (dist 1 < 3): index 0, Δ +1 — identical to nearest_parent.
        assert_eq!(
            nearest_called_by_sequence(&ca(7), 7, called, EPS, &mut scratch),
            Some((0, 1))
        );
        assert_eq!(nearest_parent(7, &[6, 10]), Some((0, 1)));
    }

    #[test]
    fn nearest_called_by_sequence_is_none_for_no_alleles() {
        let mut scratch = HmmScratch::new();
        assert_eq!(
            nearest_called_by_sequence(b"CACA", 2, &[], EPS, &mut scratch),
            None
        );
    }

    #[test]
    fn equidistant_different_length_tie_falls_to_lower_index_not_composition() {
        // A read equidistant (1 unit) from two DIFFERENT-length alleles whose compositions
        // differ: the length-6 allele mismatches the read, the length-8 one matches it.
        // Composition must NOT flip the attribution to the length-8 allele (that is only for
        // SAME-length ties); the deterministic lowest-index rule of `nearest_parent` holds, so
        // the slip Δ stays +1 (not −1). (Regression for the equidistant-slip-sign bug.)
        let a6: Vec<u8> = "CG".repeat(6).into_bytes(); // 12 bytes — mismatches the read
        let a8: Vec<u8> = "CA".repeat(8).into_bytes(); // 16 bytes — matches the read
        let read7: Vec<u8> = "CA".repeat(7).into_bytes(); // 14 bytes, 7 units
        let called: &[(&[u8], u16)] = &[(&a6, 6), (&a8, 8)];
        let mut scratch = HmmScratch::new();
        assert_eq!(
            nearest_called_by_sequence(&read7, 7, called, EPS, &mut scratch),
            Some((0, 1)),
            "equidistant different-length tie must not be flipped by composition"
        );
        assert_eq!(nearest_parent(7, &[6, 8]), Some((0, 1)));
    }

    #[test]
    fn responsibilities_normalize_over_the_called_alleles() {
        // Already-normalized input passes through.
        let r = allele_responsibilities(&[0.8, 0.2]);
        assert!((r[0] - 0.8).abs() < 1e-12 && (r[1] - 0.2).abs() < 1e-12);
        // Equal q_r → an even split (a genuine same-length het at balanced depth).
        let r = allele_responsibilities(&[0.4, 0.4]);
        assert!((r[0] - 0.5).abs() < 1e-12 && (r[1] - 0.5).abs() < 1e-12);
        // Unnormalized input is normalized by the total.
        let r = allele_responsibilities(&[3.0, 1.0]);
        assert!((r[0] - 0.75).abs() < 1e-12 && (r[1] - 0.25).abs() < 1e-12);
    }

    #[test]
    fn responsibilities_even_split_on_all_zero_evidence() {
        let r = allele_responsibilities(&[0.0, 0.0]);
        assert_eq!(r.len(), 2);
        assert!((r[0] - 0.5).abs() < 1e-12 && (r[1] - 0.5).abs() < 1e-12);
    }
}
