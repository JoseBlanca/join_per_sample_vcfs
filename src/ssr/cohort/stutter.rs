//! Stutter-kernel types — the A1 noun for slip placement (arch
//! `ssr_call_genotyping.md` §3, implementation plan §3.4).
//!
//! The kernel `S_θ(Δ) = level × shape(Δ)`, the reachability enumeration
//! `reach_variants`, and the `θ_locus` M-step are written in B2; this file defines
//! only the placement-variant shape they produce.

/// One placement-distinct realization of `candidate ⊕ Δ` (arch §4, verify-fix #3).
///
/// A pure tract has a single placement of a `Δ`-unit slip; an impure tract (runs
/// separated by interruptions) has one variant per run the slip could land in. The
/// likelihood sums `align` over these with a uniform position prior.
///
/// (Shape, not final — B2's `reach_variants` may extend it; for now it carries the
/// realized tract sequence.)
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PlacementVariant {
    /// The realized tract sequence for this placement of the slip.
    pub(crate) seq: Box<[u8]>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn placement_variant_holds_its_sequence() {
        let v = PlacementVariant {
            seq: Box::from(b"ATATATAT".as_slice()),
        };
        assert_eq!(&*v.seq, b"ATATATAT");
    }
}
