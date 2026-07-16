//! Shared ng vocabulary — the domain newtypes cross-step code speaks. It starts as this
//! one file and splits into concept modules (`units`, `locus`, …) as clusters grow
//! (`doc/devel/ng/arch/module_layout.md` principle 3). Seeded here with only what the
//! `RefSeq` reference accessor needs.

/// Which reference sequence a coordinate refers to: an index into the reference contig
/// table ([`crate::fasta::ContigList`]), in `@SQ` / `.fai` order. Unconstrained — any
/// `u32` is a legal index at the type level, and an out-of-range id is caught at fetch
/// time — so the field is public and there is no checked constructor.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct ContigId(pub u32);

impl ContigId {
    #[inline]
    pub fn get(self) -> u32 {
        self.0
    }
}

// ---------------------------------------------------------------------
// Scalar newtypes — the domain quantities cross-step code speaks. Seeded
// here with only the scalars read filtering (ng step 1) touches; the rest
// of the `ng_step_interfaces.md` §1 vocabulary lands as later steps need
// it. See `doc/devel/ng/arch/read_filtering.md` §2.1.
// ---------------------------------------------------------------------

/// SAM mapping quality (MAPQ): the aligner's Phred-scaled confidence that the
/// read is placed at the right locus. `0` = "could be anywhere"; `60` = "as
/// sure as this aligner gets". MAPQ unavailable (SAM `0xFF`) is treated as `0`
/// by callers. Unconstrained — every `u8` is a legal value, so the field is
/// public and there is no checked constructor.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct MapQual(pub u8);

impl MapQual {
    #[inline]
    pub fn get(self) -> u8 {
        self.0
    }
}

/// A single base call's Phred quality (0–93). Unconstrained — any `u8` is a
/// legal value, so the field is public and there is no checked constructor.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct BaseQual(pub u8);

impl BaseQual {
    #[inline]
    pub fn get(self) -> u8 {
        self.0
    }
}

/// A length in base pairs — the generic length currency both the SNP/indel and
/// STR paths speak (only *repeat-unit* quantities carry the `Ssr` prefix). Here
/// it measures a read's decoded length. Unconstrained — any `u64` is a legal
/// value, so the field is public and there is no checked constructor.
///
/// `u64` since B2 (spec §4): ng speaks one width, so nothing narrows, nothing is
/// checked, and no off-by-width bug is possible. Ids stay `u32` — they index a
/// table, they are not positions.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Bp(pub u64);

impl Bp {
    #[inline]
    pub fn get(self) -> u64 {
        self.0
    }
}

/// A fraction of mismatched bases, constrained to `[0, 1]`. Unlike the
/// unconstrained newtypes above, an out-of-range value is *unrepresentable*:
/// the field is private and construction goes through the checked
/// [`Self::try_new`]. Read filtering uses it as the mismatch-fraction threshold
/// (filter #8), whose source is an untrusted CLI/config value — so the policy
/// is fail loudly, never silently coerce.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct MismatchFraction(f32);

impl MismatchFraction {
    /// The only constructor. A fraction outside `[0, 1]` is a user error —
    /// reject it rather than coerce.
    pub fn try_new(x: f32) -> Result<Self, DomainError> {
        (0.0..=1.0)
            .contains(&x)
            .then_some(Self(x))
            .ok_or(DomainError::MismatchFraction(x))
    }

    #[inline]
    pub fn get(self) -> f32 {
        self.0
    }
}

/// A domain-invariant violation — the ng-wide error raised when an untrusted
/// value falls outside a constrained newtype's range. Introduced here with its
/// first variant; later constrained types (`AlleleFreq`, `InbreedingF`,
/// `Theta`, …) add their own variants as they arrive. `#[non_exhaustive]` so
/// matchers accept those future variants without breaking.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum DomainError {
    /// A [`MismatchFraction`] was constructed from a value outside `[0, 1]`.
    #[error("mismatch fraction {0} is outside [0, 1]")]
    MismatchFraction(f32),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mismatch_fraction_accepts_boundary_values() {
        assert_eq!(MismatchFraction::try_new(0.0).unwrap().get(), 0.0);
        assert_eq!(MismatchFraction::try_new(1.0).unwrap().get(), 1.0);
        assert_eq!(MismatchFraction::try_new(0.10).unwrap().get(), 0.10);
    }

    #[test]
    fn mismatch_fraction_rejects_out_of_range() {
        assert_eq!(
            MismatchFraction::try_new(-0.01),
            Err(DomainError::MismatchFraction(-0.01))
        );
        assert!(matches!(
            MismatchFraction::try_new(1.01),
            Err(DomainError::MismatchFraction(_))
        ));
        // The `[0, 1]` range check rejects the infinities as well.
        assert!(MismatchFraction::try_new(f32::INFINITY).is_err());
        assert!(MismatchFraction::try_new(f32::NEG_INFINITY).is_err());
    }

    #[test]
    fn mismatch_fraction_rejects_nan() {
        // NaN is neither <= 1 nor >= 0, so the range check rejects it.
        assert!(MismatchFraction::try_new(f32::NAN).is_err());
    }

    #[test]
    fn unconstrained_newtypes_expose_their_value() {
        assert_eq!(MapQual(20).get(), 20);
        assert_eq!(BaseQual(93).get(), 93);
        assert_eq!(Bp(150).get(), 150);
    }
}
