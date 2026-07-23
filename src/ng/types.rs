//! Shared ng vocabulary — the domain newtypes cross-step code speaks. It starts as this
//! one file and splits into concept modules (`units`, `locus`, …) as clusters grow
//! (`doc/devel/ng/arch/module_layout.md` principle 3). Seeded here with only what the
//! `RefSeq` reference accessor needs.

use std::fmt;

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

/// A **1-based** reference position — the coordinate ng speaks everywhere
/// (`ng_step_interfaces.md` §1), matching VCF/SAM/IGV and the production engine.
/// Unconstrained: any `u64` is a legal value at the type level, so the field is
/// public and there is no checked constructor.
///
/// The base is the point. ng chose 1-based so that a coordinate printed in a log,
/// a VCF, or a bug report means the same thing everywhere — no mental `+ 1`, and
/// no class of off-by-one that only shows up as a wrong genotype. The exceptions
/// are named and local: `RepeatInterval` is 0-based because it indexes a byte
/// slice, and `regions::RegionSet` is production's (converted at
/// `GenomeRegions`).
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Position(pub u64);

impl Position {
    #[inline]
    pub fn get(self) -> u64 {
        self.0
    }
}

/// One base, genome-wide: which contig, and where along it. A [`Position`] on
/// its own does not identify a base — position 1000 exists on every contig —
/// which is why so many signatures thread `(contig, position)` as two
/// parameters. This is that pair given a name, and the point-shaped sibling of
/// [`GenomeRegion`].
///
/// **The derived [`Ord`] is genome order**: contig index first, then position
/// within it, because that is the field order. That is what lets the type serve
/// directly as a sort key wherever reads or loci are ordered along the
/// reference — the read-order guard and the k-way merge are its first uses
/// (`doc/devel/ng/arch/sample_reads.md` §1.1).
///
/// The ordering is only meaningful because every alignment file's `ref_id` was
/// proved equal to its [`ContigId`] when the file was opened
/// (`doc/devel/ng/spec/alignment_file.md` §3.1). Without that gate, contig
/// indices from different files would not be comparable at all.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct GenomePosition {
    pub contig: ContigId,
    pub position: Position,
}

/// A **physical** piece of DNA: a contig plus a 1-based **inclusive** range
/// `[start, end]`. No genetic claim — that is what distinguishes it from a
/// *locus* (`typed_regions.md` §1.1).
///
/// The consolidation `ng_step_interfaces.md` §1 reserved, and the typed-region
/// generator is its first real use: it replaces `regions::Region` (0-based,
/// `u32`, production's) and `bam::ContigInterval` for ng's purposes.
///
/// **Inclusive, not half-open**, and deliberately: a half-open `end` cannot name
/// the last base of a contig without arithmetic, and "the region ends at 21"
/// meaning base 20 is exactly the ambiguity the 1-based decision exists to
/// delete. The cost is one `+ 1` in [`Self::len`], which is where it belongs —
/// stated once, in the type.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct GenomeRegion {
    pub contig: ContigId,
    pub start: Position,
    pub end: Position,
}

impl GenomeRegion {
    /// The region's length in bases. **The one place the inclusive `+ 1` lives.**
    ///
    /// Saturating rather than panicking on an inverted region: this is a plain
    /// data type with public fields (no constructor to enforce `start <= end`),
    /// and a length of 0 is a truer answer for an empty span than a panic in a
    /// getter would be. Callers that require well-formedness say so themselves.
    #[inline]
    pub fn len(self) -> u64 {
        (self.end.get() + 1).saturating_sub(self.start.get())
    }

    /// Whether the region covers no bases (`end < start`).
    #[inline]
    pub fn is_empty(self) -> bool {
        self.end.get() < self.start.get()
    }

    /// Whether `pos` falls inside the region, bounds included.
    #[inline]
    pub fn contains(self, pos: Position) -> bool {
        self.start <= pos && pos <= self.end
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
    /// A flat emission model was built from a per-base error rate that is not a finite
    /// probability in `[0, 1]`.
    #[error("per-base error rate {0} is not a finite probability in [0, 1]")]
    ErrorRate(f64),
    /// A read's bases and its qualities were paired but differ in length.
    #[error("read has {bases} bases but {qualities} qualities")]
    ReadQualityLengthMismatch { bases: usize, qualities: usize },
}

// ---------------------------------------------------------------------
// The motif — STR domain vocabulary, shared across steps
// ---------------------------------------------------------------------

/// STR scope: a repeat unit (period) is between 1 and this many bases.
pub const MAX_MOTIF_LEN: usize = 6;

/// A motif's bytes were not a valid STR period.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum MotifError {
    /// Length is `0` or above [`MAX_MOTIF_LEN`] — outside the STR period range.
    #[error("motif length {len} is outside the STR period range 1..={MAX_MOTIF_LEN}")]
    BadLength { len: usize },
}

/// A tandem-repeat unit — the repeat's period, 1..=[`MAX_MOTIF_LEN`] bases.
///
/// The bytes are stored **verbatim**: the reference-strand, phase-faithful unit
/// exactly as it tiles the locus (e.g. `CAG`), *not* canonicalized. Rotating to
/// a canonical form (`CAG` → `AGC`) would break tiling, and reconstruction reads
/// phase-correct bytes off the reference anyway; the canonical *class* used for
/// stutter pooling is derived on demand downstream, never stored here.
///
/// Inline and `Copy`: a fixed 6-byte buffer plus a length, never heap-allocated.
/// Unused tail bytes are zero, so the derived `Eq`/`Hash` compare only the live
/// prefix (`0` is not a valid base, so it cannot collide with one).
///
/// INVARIANT (relied on by the derived `Eq`/`Hash`): every constructor MUST
/// zero-initialize the unused tail of `buf`. [`Motif::new`] is currently the
/// only one; any future constructor must uphold this or the derived impls will
/// treat equal motifs as distinct.
///
/// ## Why this is a port, when spec §4 said to reuse production's
///
/// Spec §4 kept `ssr::types::Motif` on the grounds that it *"carries no
/// coordinates and no width, and so has nothing to rebase — the Revision's
/// 'reuse where it costs production nothing' case exactly."* The first half is
/// true; **the conclusion was wrong**, and the compiler said so.
///
/// `ssr::types::Motif` is `pub(crate)`. ng's [`SsrSegment`] is `pub` (the ng-sibling
/// convention) and returns a motif, so reusing it trips rustc's
/// `private_interfaces` lint — a `pub` item leaking a `pub(crate)` type. The
/// three ways out: widen `Motif` in `src/ssr/types.rs` (**touching production —
/// forbidden**); demote ng's whole classification surface to `pub(crate)` (bends ng's
/// convention *and* buys `dead_code` warnings for every item until its Milestone
/// D consumer exists); or port the 40 lines. So reuse did **not** cost
/// production nothing — it cost a visibility compromise, which is precisely the
/// coupling "a fresh ng caller from scratch" (owner, 2026-07-16) exists to
/// avoid.
///
/// Ported, ng's `region_typing` names nothing from `src/ssr/` outside its
/// `#[cfg(test)]` differential — which is exactly where the dependency belongs.
/// The type is coordinate-free and trivially checkable, so the duplication is
/// cheap and the drift risk is near zero; the differential compares motifs by
/// bytes and would catch it anyway.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Motif {
    buf: [u8; MAX_MOTIF_LEN],
    len: u8,
}

impl Motif {
    /// Build a motif from its bytes, validating the STR period range.
    ///
    /// The bytes are taken verbatim (no canonicalization); the caller is
    /// responsible for supplying the phase-faithful, reference-strand unit.
    pub fn new(bytes: &[u8]) -> Result<Self, MotifError> {
        let len = bytes.len();
        if len == 0 || len > MAX_MOTIF_LEN {
            return Err(MotifError::BadLength { len });
        }
        let mut buf = [0u8; MAX_MOTIF_LEN];
        buf[..len].copy_from_slice(bytes);
        Ok(Self {
            buf,
            len: len as u8,
        })
    }

    /// The motif bytes.
    #[inline]
    pub fn as_bytes(&self) -> &[u8] {
        &self.buf[..self.len as usize]
    }

    /// The period (motif length, in bases).
    #[inline]
    pub fn period(&self) -> usize {
        self.len as usize
    }
}

impl fmt::Debug for Motif {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Bases are ASCII; render as text for readable test output / logs,
        // falling back to bytes if a motif ever held non-UTF-8.
        match std::str::from_utf8(self.as_bytes()) {
            Ok(s) => write!(f, "Motif({s:?})"),
            Err(_) => write!(f, "Motif({:?})", self.as_bytes()),
        }
    }
}

// Where `Motif` used to live, and why it moved. It was written in
// `region_typing::segment_criteria`, because that is where the STR classifier needed it. It
// then gained consumers in `locus_generation` and `alignment` — three modules across
// pipeline-stage boundaries — and `alignment` states in its own module doc that it is **not
// a pipeline step and knows no callers**, which an import from step 3 flatly contradicted,
// on ng's *public* surface. `module_layout.md` principle 3 already assigns shared STR domain
// vocabulary here. `segment_criteria` re-exports it, so step 3's own callers are untouched.

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
        assert_eq!(Position(1).get(), 1);
    }

    fn genome_position(contig: u32, position: u64) -> GenomePosition {
        GenomePosition {
            contig: ContigId(contig),
            position: Position(position),
        }
    }

    /// The whole reason `GenomePosition` is a struct with this field order and
    /// not a `(Position, ContigId)` pair: sorting must give **genome order** —
    /// contig-major, position-minor — with no comparator written at the call
    /// site. A transposed field order would still compile and still sort; it
    /// would just interleave contigs, which is the failure this test exists to
    /// catch.
    #[test]
    fn sorting_yields_contig_major_position_minor_order() {
        let mut shuffled = vec![
            genome_position(2, 5),
            genome_position(0, 900),
            genome_position(1, 1),
            genome_position(0, 7),
            genome_position(2, 1),
            genome_position(1, 1000),
            genome_position(0, 100),
        ];
        shuffled.sort();

        assert_eq!(
            shuffled,
            vec![
                genome_position(0, 7),
                genome_position(0, 100),
                genome_position(0, 900),
                genome_position(1, 1),
                genome_position(1, 1000),
                genome_position(2, 1),
                genome_position(2, 5),
            ],
            "contig index dominates: contig 2 position 1 sorts after contig 1 position 1000"
        );
    }

    /// Equal keys are legal and compare equal — the order guard rejects only a
    /// strict decrease (`spec/alignment_file.md` §3.2), so `Ord` must report
    /// `Equal` for a repeated position rather than imposing a tie-break of its
    /// own, and `Less` for any genuine advance along the genome.
    #[test]
    fn ordering_is_equal_for_a_repeated_position_and_strict_otherwise() {
        use std::cmp::Ordering;

        let read_start = genome_position(3, 42);
        assert_eq!(read_start.cmp(&genome_position(3, 42)), Ordering::Equal);
        assert!(
            read_start < genome_position(3, 43),
            "later position on same contig"
        );
        assert!(
            read_start < genome_position(4, 1),
            "later contig, earlier position"
        );
    }

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    /// **Inclusive is the whole point**, so `len` is where the `+ 1` lives — and
    /// the single-base case is the one that catches a half-open slip.
    #[test]
    fn region_len_counts_both_ends() {
        assert_eq!(region(1, 1).len(), 1, "a single base is length 1, not 0");
        assert_eq!(region(1, 10).len(), 10);
        assert_eq!(
            region(6, 21).len(),
            16,
            "the classification fixture's tract"
        );
        assert!(!region(1, 1).is_empty());
    }

    /// A region can be built inverted (public fields, no constructor), so the
    /// accessors must answer rather than panic.
    #[test]
    fn an_inverted_region_is_empty_not_a_panic() {
        assert!(region(10, 9).is_empty());
        assert_eq!(region(10, 9).len(), 0);
        // Saturating, not wrapping: a wildly inverted span is still 0, not ~u64::MAX.
        assert_eq!(region(u64::MAX, 1).len(), 0);
    }

    #[test]
    fn region_contains_includes_both_bounds() {
        let r = region(6, 21);
        assert!(r.contains(Position(6)), "start is inside");
        assert!(r.contains(Position(21)), "end is inside — inclusive");
        assert!(!r.contains(Position(5)));
        assert!(!r.contains(Position(22)));
    }

    /// Position 0 is representable but meaningless, since 1-based coordinates
    /// begin at one. Recorded rather than enforced: `Position` is unconstrained,
    /// as `MapQual` and `Bp` are, and every place where a 0 is a *bug* rejects it
    /// itself — `RefSeq::fetch_into` with `InvalidStart`, `Locus::new` with
    /// `BadCoordinates`, `classify` with an assert.
    #[test]
    fn position_zero_is_representable_and_rejected_where_it_matters() {
        assert_eq!(Position(0).get(), 0);
        assert!(region(0, 5).contains(Position(0)));
    }
}
