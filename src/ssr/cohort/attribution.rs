//! Shared read-to-allele attribution primitive.
//!
//! Several stages hard-attribute each spanning read to the nearest parent allele by
//! repeat-length distance and then act on the signed slip `Δ = read_units − parent_units`:
//! the per-locus slip refit ([`em::attribute_locus`](super::em)), the pre-pass confident-
//! genotype stats ([`prepass::accumulate_locus`](super::prepass)), and the apparent allele
//! balance ([`vcf_out::allele_balance`](super::vcf_out)). They must all break a distance tie
//! the same way, so the decision lives here once rather than being re-derived per site
//! (review Mi1). The soft per-read responsibility split is a deferred refinement that would
//! replace this hard assignment in every caller at once.

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
}
