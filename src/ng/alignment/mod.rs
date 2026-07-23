//! ng alignment — lining a read up against a reference stretch.
//!
//! This is a **folder, not a file**, because it holds competing implementations that get
//! compared: the three aligner and normalizer traits will live here in `mod.rs` and each
//! algorithm sits in its own file beside its siblings — as does the shared [`emission`]
//! component (`doc/devel/ng/arch/alignment.md` §Module home). File names say
//! which trait an algorithm implements and whether it is repeat-aware, because both are
//! load-bearing — the best-path and marginal families implement *different* traits, and
//! every repeat-aware algorithm has an affine counterpart it must not be confused with.
//!
//! **It is not a pipeline step, and it knows none.** Step 2 (read preparation) and step 7
//! (read likelihood) both call it; it knows about neither. It takes sequences and returns
//! alignments or probabilities, which is why it sits outside the step-per-module rule
//! rather than breaking it (`doc/devel/ng/spec/alignment.md` §1).
//!
//! Landed so far: [`Alignment`] and [`RepeatSpan`], the two output shapes; [`RepeatGeometry`]
//! and [`RepeatContext`], which tell a repeat-aware algorithm where the repeat sits and how it
//! slips; the [`BestPathAligner`] trait they meet at; and the two shared components, the
//! [`emission`] scoring every aligner uses and the [`stutter`] model; and the first
//! algorithm, [`ssr_best_path_flat_gap`] — the repeat-aware delimiter, which is this
//! module's **only byte-parity oracle** against production (`delimit_parity`, test-only).
//! Still to come in the plan `doc/devel/ng/impl_plan/alignment_best_path.md`: banding
//! (Milestone C), the two-penalty aligner (D), and the affine aligner (E, gated).

#[cfg(test)]
mod delimit_parity;
pub mod emission;
pub mod ssr_best_path_flat_gap;
pub mod ssr_best_path_unit_slip;
pub mod stutter;

pub use emission::{BaseScores, Emission, FlatEmission, PerQualityEmission};
pub use stutter::{MAX_SLIP, StutterModel, StutterRates};

use crate::ng::types::{BaseQual, Bp, DomainError, Motif};
use crate::pileup::walker::CigarOp;
use std::ops::Range;

/// A read placed against a reference stretch: where the placement starts, and the
/// operations that spell it out from there.
///
/// This is the output shape of the **affine** best-path aligner — the general-purpose
/// one. The repeat-aware aligners return where a read's repeat sits instead, which is a
/// different shape and a different type; that is why the aligner trait's output is an
/// associated type rather than this struct (arch §2.1, §3).
///
/// # Reusing production's `CigarOp`
///
/// The operations are production's [`CigarOp`], **reused rather than re-minted**. A
/// parallel operation type would buy nothing and cost a conversion at every boundary
/// where an alignment meets code that already speaks CIGAR — and unlike ng's [`Motif`],
/// which had to be ported because production's is `pub(crate)` and would have leaked out
/// of a `pub` signature, `CigarOp` is already fully public, so reuse here costs
/// production nothing (arch §5).
///
/// Worth naming, because the import path misleads: `CigarOp` lives inside
/// `pileup::walker`, so this module appears to depend on the pileup walker, which it does
/// not — it knows no callers at all. `CigarOp` is really crate-wide CIGAR vocabulary
/// (`bam/`, `ssr/pileup/`, `pop_var_caller/cli/`, the benches and ng all consume it), and
/// its home inside a pipeline stage is a **pre-existing misplacement**. Lifting it to a
/// top-level peer is a production edit, and production is frozen (owner, 2026-07-16), so
/// this is recorded as known debt rather than fixed here; porting this module back is the
/// natural moment to do it.
///
/// [`Motif`]: crate::ng::types::Motif
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    /// Where the placement starts, measured in bases from the start of **the reference
    /// stretch the aligner was given** — not a genome position. The aligner is handed a
    /// slice and knows nothing about where that slice came from, so it cannot return a
    /// coordinate; turning this into a genome position is the caller's job, and it needs
    /// the stretch's own start to do it.
    ///
    /// `u64`, not `usize`: ng speaks one width for coordinates and lengths, so nothing
    /// narrows and no off-by-width bug is possible ([`Bp`], and
    /// [`RepeatInterval`] — the same shape, an offset into the slice its
    /// producer was handed, `u64` for the same reason).
    ///
    /// [`Bp`]: crate::ng::types::Bp
    /// [`RepeatInterval`]: crate::ng::tandem_repeat::RepeatInterval
    pub reference_offset: u64,
    /// The operations from `reference_offset` onward, in read order.
    ///
    /// Named `cigar` rather than `ops` to match how every other site in the crate spells
    /// this concept (`cigar` / `cigar_ops`); the architecture doc's sketch says `ops`, but
    /// it states that its signatures are illustrative and the contract is what ships.
    ///
    /// **An affine aligner inhabits only part of [`CigarOp`].** The enum is production's
    /// *parser* vocabulary and must represent everything a BAM can contain, including
    /// `Skip`, `HardClip` and `Padding`, which describe an input record rather than a
    /// placement this module computed. Exactly which subset the affine aligner emits is
    /// **not settled here** — its algorithm lands in Milestone E and is currently gated,
    /// and arch §6 still has the aligner's output shape open. The obligation is recorded
    /// rather than guessed: the step that writes the producer states the subset, and only
    /// then may a consumer treat the remaining variants as unreachable.
    pub cigar: Vec<CigarOp>,
}

// A note on derives: `Hash` is not derived because [`CigarOp`] does not derive it, and
// adding it there is a production edit. If an `Alignment` ever needs to be a map key,
// that is the blocker to clear first.

/// Where a read's repeat lies, in read coordinates, **and what held it there**.
///
/// A span is only a *measurement* when both flanks anchored. Otherwise the read ran off its
/// own end somewhere inside the repeat, and what the span records is a **lower bound** on
/// the repeat's length, not the length. Those are different facts, and a caller that
/// confuses them counts a truncated read as a short allele.
///
/// This is where ng goes past production. Production's `Delimited`
/// ([alignment.rs](../../../ssr/pileup/alignment.rs)) has a single `BorderOffEnd` variant
/// covering every unanchored case, so it is **side-blind**: it cannot say *which* flank was
/// missing, and the STR read preparer cannot otherwise tell a measurement from a bound
/// (arch §2.1, §5; spec §8).
///
/// The offsets are read coordinates — `u64` and 0-based, indexing the read slice the
/// aligner was handed, exactly like [`RepeatInterval`] and for the same reason: ng speaks
/// one width, so nothing narrows and no off-by-width bug is possible.
///
/// A [`Range`] rather than ng's existing [`RegionSpan`], deliberately: `RegionSpan` is the
/// region seam's *reference*-coordinate vocabulary, and this module knows no coordinates
/// and no callers. The cost is that this type is not `Copy`.
///
/// **Producer's invariant:** every span is well-formed, `start <= end`. Nothing enforces
/// it — the aligner is the only producer — and an inverted range would read as a
/// zero-length tract, which is the same short-allele failure the enum exists to prevent.
/// [`Self::length_lower_bound`] saturates rather than wrapping so a violation cannot turn
/// into an enormous length.
///
/// [`RepeatInterval`]: crate::ng::tandem_repeat::RepeatInterval
/// [`RegionSpan`]: crate::ng::tandem_repeat::RegionSpan
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RepeatSpan {
    /// Both flanks anchored, so the span **pins the repeat's length**. The only variant
    /// that is a measurement.
    Between(Range<u64>),
    /// Only the left flank anchored: the read ends before the right flank is reached, so
    /// the span runs to the end of the read and the true repeat is **at least** this long.
    FromLeft(Range<u64>),
    /// Only the right flank anchored: the read begins after the left flank, so the span
    /// starts at the read's start and is likewise a lower bound.
    ///
    /// (Described by read start and end rather than 5′/3′ on purpose: reads are held
    /// reference-forward, so strand vocabulary would invert for a reverse-strand read and
    /// mean the opposite of what it says.)
    FromRight(Range<u64>),
    /// Neither flank anchored — the read lies wholly inside the repeat. It carries **no
    /// per-read fact** about this repeat's length beyond "at least the whole read", which
    /// is why this variant has no span to report.
    Unanchored,
}

impl RepeatSpan {
    /// The repeat's length in the read — **only when the read actually measured it**.
    ///
    /// This is the accessor a genotype path wants, and the reason it exists separately
    /// from [`Self::observed_span`]: taking a length off the raw span is a one-line
    /// mistake that leaves no trace. `observed_span().map(|s| s.end - s.start)` returns a
    /// perfectly plausible number for a read that ran off its own end, and that number is
    /// a **short allele that was never observed**. The correct version differs from the
    /// wrong one by an `if`, so the safe form is given a name and the unsafe one is not
    /// reachable by accident.
    #[must_use]
    pub fn measured_length(&self) -> Option<u64> {
        match self {
            Self::Between(span) => Some(span.end.saturating_sub(span.start)),
            Self::FromLeft(_) | Self::FromRight(_) | Self::Unanchored => None,
        }
    }

    /// The largest repeat length this read rules out — a length the repeat is **at least**.
    ///
    /// Defined for every case, because even an unanchored read says something: the repeat
    /// is at least as long as the stretch of it the read shows. For [`Self::Between`] the
    /// bound is the measurement itself.
    ///
    /// **The [`Self::Unanchored`] answer presumes the read was aligned against a real
    /// locus.** That variant is also what an *empty* reference produces — there was no
    /// locus to delimit — and there the `read_len` fallback would assert a repeat at least
    /// as long as the read at a locus with no bases at all. The caller knows which case it
    /// is in; this accessor cannot.
    #[must_use]
    pub fn length_lower_bound(&self, read_len: u64) -> u64 {
        match self {
            Self::Between(span) | Self::FromLeft(span) | Self::FromRight(span) => {
                span.end.saturating_sub(span.start)
            }
            // The read lies wholly inside the repeat, so all of it is repeat.
            Self::Unanchored => read_len,
        }
    }

    /// The stretch of the read the aligner read out as repeat, **whether or not it
    /// measures the repeat's length**.
    ///
    /// For extracting the observed bases — which is valid for a lower bound too, and is
    /// what an interrupted repeat needs, since its sequence must come out verbatim. Do
    /// **not** take a length from this; use [`Self::measured_length`] or
    /// [`Self::length_lower_bound`], which cannot lie about which one they are.
    #[must_use]
    pub fn observed_span(&self) -> Option<&Range<u64>> {
        match self {
            Self::Between(span) | Self::FromLeft(span) | Self::FromRight(span) => Some(span),
            Self::Unanchored => None,
        }
    }

    /// Whether this span **measures** the repeat, rather than bounding it below.
    ///
    /// Written as an exhaustive `match` rather than `matches!`, deliberately. `matches!`
    /// compiles an implicit `_ => false`, so a fifth variant that *is* a measurement —
    /// an interrupted-repeat case, say — would silently read as a lower bound and yield
    /// systematically short alleles, with the compiler saying nothing. This is the one
    /// place in the type where a wildcard would absorb exactly the change that matters.
    #[must_use]
    pub fn is_measurement(&self) -> bool {
        match self {
            Self::Between(_) => true,
            Self::FromLeft(_) | Self::FromRight(_) | Self::Unanchored => false,
        }
    }
}

/// Where the repeat sits inside the reference stretch the aligner is given, measured in
/// bases from that stretch's start.
///
/// **Both flank lengths are measured, never assumed equal.** A repeat near a contig end has
/// a short flank or none at all, and code that halves a total or reuses one flank's length
/// for the other silently mis-locates the tract exactly where the evidence is thinnest
/// (arch §2.2).
///
/// This travels **per call, not as constructor state**: it changes at every locus, so
/// holding it would mean building a new aligner per locus across millions of loci. That is
/// production's own shape — a stateless model plus a per-call context (arch §2.2, §4).
///
/// (Arch §2.2 adds "and `Motif` is an allocation" as a second reason. That part is not true
/// of ng's [`Motif`], which is `Copy` over an inline `[u8; 6]` and allocates nothing. The
/// per-locus cost stands on its own.)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepeatGeometry {
    /// Bases of flank before the repeat starts. May be zero at a contig start.
    pub left_flank_len: Bp,
    /// Bases of flank after the repeat ends. May be zero at a contig end.
    pub right_flank_len: Bp,
    /// The repeat unit. Needed only by the algorithms that price whole-unit slips
    /// (algorithm 4); the flat-gap delimiter needs the tract boundaries alone.
    pub motif: Motif,
}

impl RepeatGeometry {
    /// Whether this geometry fits a reference stretch of `reference_len` bases — that is,
    /// whether the two flanks leave room for the repeat between them.
    ///
    /// **This is the caller's precondition, given a name.** Arch §3 makes it the caller's
    /// to uphold because nothing on the hot path can check it cheaply, and production is
    /// safe only because its locus type enforces it upstream: production computes the tract
    /// boundary by subtracting the right flank from the reference length, which would
    /// **underflow** if this were violated. Implementations assert this in debug builds;
    /// having it as a function means a caller that wants to check can, without recomputing
    /// the reasoning.
    ///
    /// Defined as "[`Self::repeat_len`] has an answer" rather than by its own arithmetic,
    /// so the two cannot drift into disagreeing about the boundary case.
    #[must_use]
    pub fn fits_reference(&self, reference_len: Bp) -> bool {
        self.repeat_len(reference_len).is_some()
    }

    /// How many reference bases the repeat itself spans, inside a stretch of
    /// `reference_len` bases: the stretch minus both flanks.
    ///
    /// `None` when the geometry does not fit, rather than a saturated zero — a zero would
    /// be indistinguishable from a genuinely empty tract, and this is the arithmetic whose
    /// unchecked form underflows.
    #[must_use]
    pub fn repeat_len(&self, reference_len: Bp) -> Option<Bp> {
        reference_len
            .get()
            .checked_sub(self.left_flank_len.get())?
            .checked_sub(self.right_flank_len.get())
            .map(Bp)
    }
}

/// What varies from one repeat to the next — the per-call argument every repeat-aware
/// algorithm is handed.
///
/// Both parts change at every locus, which is why they travel **per call rather than as
/// constructor state**: holding them would mean building a new aligner per locus across
/// millions of loci. This mirrors production, whose read-likelihood model is a stateless
/// value taking a per-call scoring context (arch §2.2, §4).
///
/// Borrowed rather than owned, so a caller that already holds a geometry and a model for a
/// locus pays nothing to pass them to each of that locus's reads.
///
/// *(The plan lists this in step A2 with [`RepeatGeometry`]; it landed in A3, because it
/// names a [`StutterModel`] and that type arrives there.)*
#[derive(Debug, Clone, Copy)]
pub struct RepeatContext<'a> {
    /// Where the repeat sits in the reference stretch.
    pub geometry: &'a RepeatGeometry,
    /// How the repeat slips. **Only the algorithms that price whole-unit slips read this**
    /// — algorithm 4. The flat-gap delimiter (algorithm 3) ignores it, which is the point:
    /// its measurement is content-agnostic and must not be pulled toward tidy lengths.
    pub stutter: &'a StutterModel,
}

/// A read's bases and their per-base qualities, held together because they are only
/// meaningful together.
///
/// **This type exists to dissolve a precondition rather than document one.** The two slices
/// have to be the same length — quality *i* belongs to base *i* — and the aligner trait used
/// to take them as two separate `&[u8]` arguments with that stated as a caller obligation
/// plus a `debug_assert!`. Two problems with that. The arguments were positionally
/// interchangeable with each other *and* with the reference, so a transposition was a silent
/// wrong answer; and the assertion compiled out of the release build the project runs, so
/// the obligation was unenforced exactly where it mattered.
///
/// Checking once at construction costs nothing on the hot path — this is built per read,
/// while the thing arch §3's no-`Result` rule defends is per *matrix cell* — and it makes the
/// invariant true by construction everywhere downstream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadBases<'a> {
    bases: &'a [u8],
    qualities: &'a [u8],
}

impl<'a> ReadBases<'a> {
    /// Pair a read's bases with its qualities.
    ///
    /// # Errors
    ///
    /// [`DomainError::ReadQualityLengthMismatch`] if the two slices differ in length.
    pub fn try_new(bases: &'a [u8], qualities: &'a [u8]) -> Result<Self, DomainError> {
        if bases.len() != qualities.len() {
            return Err(DomainError::ReadQualityLengthMismatch {
                bases: bases.len(),
                qualities: qualities.len(),
            });
        }
        Ok(Self { bases, qualities })
    }

    /// The read's bases.
    #[must_use]
    pub fn bases(&self) -> &'a [u8] {
        self.bases
    }

    /// The per-base qualities, the same length as [`Self::bases`].
    #[must_use]
    pub fn qualities(&self) -> &'a [u8] {
        self.qualities
    }

    /// How many bases the read has.
    #[must_use]
    pub fn len(&self) -> usize {
        self.bases.len()
    }

    /// Whether the read has no bases. A legal, if useless, read.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    /// The quality of base `index`, as the domain newtype.
    #[inline]
    #[must_use]
    pub fn quality_at(&self, index: usize) -> BaseQual {
        BaseQual(self.qualities[index])
    }
}

/// Line a read up against a reference sequence the single most probable way.
///
/// # What varies between implementations
///
/// The three associated types are what let one trait cover both families. `Output` is
/// [`Alignment`] for the affine aligner and [`RepeatSpan`] for the repeat-aware ones —
/// they answer different questions, so a single output type would force one of them to
/// lie. `Context` is `()` for the affine aligner, which needs nothing locus-specific, and
/// [`RepeatContext`] for the others. `Scratch` is per-algorithm, never one shared buffer.
///
/// **`Context` carries a lifetime, and is passed by value.** The architecture sketch has a
/// plain `type Context` taken by reference, but the repeat context *is* two references —
/// the geometry and the stutter model both travel per call and are borrowed, so that they
/// need not be rebuilt per locus across millions of loci. A plain associated type could
/// only name that as `RepeatContext<'static>`, which would force every caller to own its
/// locus data forever. So `Context<'a>` is generic over the lifetime, and since the value
/// is two pointers and `Copy`, taking a reference *to* it would only add indirection.
///
/// # Contract
///
/// **Pure per call.** The output is a function of the arguments and the implementation's
/// own constructor state, with no hidden mutation beyond `Scratch` — so results never
/// depend on call order or thread count, which the cohort's byte-identity guarantee rests
/// on. `Scratch` is caller-owned and reused; an implementation that allocates per call is a
/// defect, not a slow path.
///
/// **Determinism, and why it is a correctness rule rather than a nicety.** Two reads of the
/// same molecule must measure the same repeat. Where several line-ups score equally, an
/// arbitrary choice between them lets identical input produce different repeat lengths, so
/// two rules settle every tie (spec §4.2):
///
/// - prefer **match, then deletion, then insertion**, both within the alignment and when
///   choosing which state the alignment ends in;
/// - a gap landing exactly on a flank/repeat junction belongs to the block on its **5′
///   side** — an insertion at the left junction joins the flank, one at the right junction
///   joins the repeat.
///
/// **No `Result`, deliberately.** A fallible alignment would push error handling onto the
/// caller's hottest path for cases that are answers, not failures: an empty reference gives
/// [`RepeatSpan::Unanchored`], not an error. Two conditions are instead the **caller's** to
/// uphold, because nothing here can check them cheaply:
///
/// - the geometry must fit the reference it is used with — [`RepeatGeometry::fits_reference`].
///
/// It is a documented precondition plus a debug assertion, and it is worth being blunt about
/// what that means: **a debug assertion compiles out of the release build this project
/// runs**, so a violated precondition is a wrong answer rather than a crash. If it turns out
/// to be reachable from untrusted input, it becomes a checked constructor on the context
/// type — not a `Result` on this call (arch §3).
///
/// Arch §3 named a *second* such condition — that the read and its qualities be the same
/// length. That one is no longer the caller's to remember: [`ReadBases`] makes it true by
/// construction, which is the same remedy arch §3 prescribes, applied one step earlier.
///
/// The `Sized` supertrait is deliberate: implementations are selected by generic type
/// parameter, never behind `dyn`, because this runs per read across a whole cohort (arch
/// §4). Stating it here makes that a compile error rather than a convention.
pub trait BestPathAligner: Sized {
    /// Reused matrices and traceback buffers — allocated per worker, never per read.
    ///
    /// **`Default` must decide nothing that changes a result.** The bound is the only way
    /// a caller can build one, so whatever `Default` picks is invisible at every call
    /// site: a scratch that defaulted a band half-width or a matrix cap would be choosing
    /// the answer, and a bake-off comparing two aligners could not see or report the
    /// value it was comparing at. Scratch is *buffers only* — empty, grown on demand, and
    /// inert with respect to output. An aligner that needs a result-affecting parameter
    /// puts it in its own constructor state, where it is visible, or in the per-call
    /// context. If one ever needs a scratch it must configure, the escape hatch is a
    /// `fn scratch(&self) -> Self::Scratch` on this trait, not a cleverer `Default`.
    type Scratch: Default;
    /// [`Alignment`] for the affine aligner; [`RepeatSpan`] for the repeat-aware ones.
    type Output;
    /// `()` for the affine aligner; [`RepeatContext`] for the repeat-aware ones.
    ///
    /// `Copy` because the by-value signature rests on it: a caller aligning many reads
    /// against one locus has to hand the same context to each, and without the bound a
    /// generic helper doing that would not compile. Both implementations are two pointers
    /// or nothing at all; an implementation with a large owned context passes `&'a Big`.
    type Context<'a>: Copy;

    /// Align `read` against `reference`.
    ///
    /// `reference` is a bare stretch of bases, not a genome region — this module knows no
    /// coordinates and no callers. The read carries its qualities with it ([`ReadBases`]),
    /// so the two cannot be mismatched or transposed here.
    fn align(
        &self,
        read: ReadBases<'_>,
        reference: &[u8],
        context: Self::Context<'_>,
        scratch: &mut Self::Scratch,
    ) -> Self::Output;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A0 lands no algorithm, so this test's job is to be a **compile-time anchor**, and
    /// it is written to fail if the things it anchors change: that the type keeps both
    /// fields, derives structural equality over *both* of them, and still builds its
    /// operations out of production's [`CigarOp`].
    ///
    /// The whole-struct comparisons are the part that carries weight — a hand-written
    /// `PartialEq` that ignored either field would pass field-by-field assertions and
    /// fail these.
    ///
    /// TODO(Milestone E): the two contract properties stated on the fields have no test
    /// and cannot have one until a producer exists. The step that writes the affine
    /// aligner owes `align_returns_offset_relative_to_the_reference_stretch` (the
    /// stretch-versus-genome confusion is a wrong answer, not a panic) and
    /// `align_returns_operations_in_read_order`.
    #[test]
    fn alignment_distinguishes_offset_and_operation_order() {
        let alignment = Alignment {
            reference_offset: 7,
            cigar: vec![CigarOp::Match(10), CigarOp::Insertion(2), CigarOp::Match(5)],
        };

        assert_eq!(alignment.reference_offset, 7);
        assert_eq!(
            alignment.cigar,
            vec![CigarOp::Match(10), CigarOp::Insertion(2), CigarOp::Match(5)]
        );

        // Equality must see the offset: same operations, different placement.
        assert_ne!(
            alignment,
            Alignment {
                reference_offset: 8,
                cigar: alignment.cigar.clone(),
            }
        );

        // Equality must see operation order: same operations, permuted.
        assert_ne!(
            alignment,
            Alignment {
                reference_offset: 7,
                cigar: vec![CigarOp::Insertion(2), CigarOp::Match(10), CigarOp::Match(5)],
            }
        );
    }

    fn geometry(left_flank_len: u64, right_flank_len: u64) -> RepeatGeometry {
        RepeatGeometry {
            left_flank_len: Bp(left_flank_len),
            right_flank_len: Bp(right_flank_len),
            motif: Motif::new(b"CAG").expect("CAG is a valid period-3 motif"),
        }
    }

    /// **The mismatch that used to be a documented precondition.** It was the caller's to
    /// remember, backed by a `debug_assert!` that compiled out of the release build — so the
    /// one build that matters did not check it. Now it cannot be built wrong.
    #[test]
    fn read_bases_rejects_qualities_that_do_not_match_the_bases() {
        assert!(ReadBases::try_new(b"ACGT", &[30, 30, 30, 30]).is_ok());
        assert!(ReadBases::try_new(b"", &[]).is_ok());

        assert!(matches!(
            ReadBases::try_new(b"ACGT", &[30, 30, 30]),
            Err(DomainError::ReadQualityLengthMismatch {
                bases: 4,
                qualities: 3
            })
        ));
        assert!(ReadBases::try_new(b"ACGT", &[]).is_err());
        assert!(ReadBases::try_new(b"", &[30]).is_err());
    }

    #[test]
    fn read_bases_hands_back_what_it_was_given() {
        let read = ReadBases::try_new(b"ACGT", &[10, 20, 30, 40]).expect("matched lengths");
        assert_eq!(read.bases(), b"ACGT");
        assert_eq!(read.qualities(), &[10, 20, 30, 40]);
        assert_eq!(read.len(), 4);
        assert!(!read.is_empty());
        assert_eq!(read.quality_at(2), BaseQual(30));
        assert!(
            ReadBases::try_new(b"", &[])
                .expect("empty is legal")
                .is_empty()
        );
    }

    /// **The distinction the type exists to carry.** Only `Between` measures a repeat;
    /// every other case is a lower bound, and a consumer that treats one as the other
    /// counts a truncated read as a short allele — a wrong genotype, with nothing in the
    /// output to show for it.
    #[test]
    fn only_a_two_flank_span_is_a_measurement() {
        assert!(RepeatSpan::Between(4..10).is_measurement());
        assert!(!RepeatSpan::FromLeft(4..10).is_measurement());
        assert!(!RepeatSpan::FromRight(4..10).is_measurement());
        assert!(!RepeatSpan::Unanchored.is_measurement());
    }

    /// **A length is only available where one was measured.** This is the accessor split
    /// that keeps the wrong answer out of reach: all three anchored variants carry an
    /// identical span, so anything that reads a length off the span alone cannot tell a
    /// measurement from a truncated read.
    #[test]
    fn only_a_two_flank_span_yields_a_measured_length() {
        assert_eq!(RepeatSpan::Between(4..10).measured_length(), Some(6));
        assert_eq!(RepeatSpan::FromLeft(4..10).measured_length(), None);
        assert_eq!(RepeatSpan::FromRight(4..10).measured_length(), None);
        assert_eq!(RepeatSpan::Unanchored.measured_length(), None);

        // The trap, spelled out: the raw span cannot distinguish them.
        assert_eq!(
            RepeatSpan::Between(4..10).observed_span(),
            RepeatSpan::FromLeft(4..10).observed_span()
        );
    }

    /// Every case bounds the length below, including the unanchored one — a read lying
    /// wholly inside a repeat still proves the repeat is at least as long as the read.
    #[test]
    fn every_case_bounds_the_length_below() {
        let read_len = 100;
        assert_eq!(RepeatSpan::Between(4..10).length_lower_bound(read_len), 6);
        assert_eq!(RepeatSpan::FromLeft(4..10).length_lower_bound(read_len), 6);
        assert_eq!(RepeatSpan::FromRight(4..10).length_lower_bound(read_len), 6);
        assert_eq!(
            RepeatSpan::Unanchored.length_lower_bound(read_len),
            read_len
        );
    }

    /// Degenerate and extreme payloads. The inverted range is the producer-invariant
    /// violation the type documents: it must saturate to zero rather than wrap to an
    /// enormous length, since a huge allele is a far louder wrong answer than a short one
    /// but neither should be manufactured by arithmetic.
    #[test]
    fn span_lengths_handle_empty_inverted_and_extreme_payloads() {
        assert_eq!(RepeatSpan::Between(7..7).measured_length(), Some(0));
        // Built field-by-field: a literal `10..4` is rejected outright by clippy's
        // `reversed_empty_ranges`, which is itself a good sign — the invariant is one the
        // language's own tooling treats as a mistake.
        let inverted = Range { start: 10, end: 4 };
        assert_eq!(
            RepeatSpan::Between(inverted.clone()).measured_length(),
            Some(0)
        );
        assert_eq!(RepeatSpan::Between(inverted).length_lower_bound(50), 0);
        assert_eq!(
            RepeatSpan::Between(0..u64::MAX).measured_length(),
            Some(u64::MAX)
        );
    }

    /// All four cases are distinguishable, including the two one-flank cases from each
    /// other — a mis-assigned side is a wrong observation class that still looks like a
    /// perfectly good read, which is why the sides are separate variants rather than one
    /// side-blind case as in production's `Delimited`.
    #[test]
    fn the_two_one_flank_cases_are_not_interchangeable() {
        assert_ne!(RepeatSpan::FromLeft(4..10), RepeatSpan::FromRight(4..10));
        assert_ne!(RepeatSpan::Between(4..10), RepeatSpan::FromLeft(4..10));
        assert_ne!(RepeatSpan::Between(4..10), RepeatSpan::Unanchored);
    }

    #[test]
    fn an_observed_span_is_reported_for_every_case_but_the_unanchored_one() {
        assert_eq!(RepeatSpan::Between(4..10).observed_span(), Some(&(4..10)));
        assert_eq!(RepeatSpan::FromLeft(4..10).observed_span(), Some(&(4..10)));
        assert_eq!(RepeatSpan::FromRight(4..10).observed_span(), Some(&(4..10)));
        assert_eq!(RepeatSpan::Unanchored.observed_span(), None);
    }

    /// The flanks are **measured, never assumed equal** — so an asymmetric geometry, which
    /// is what a repeat near a contig end produces, must round-trip intact rather than
    /// being normalised to a common flank length.
    #[test]
    fn geometry_keeps_the_two_flanks_apart() {
        let asymmetric = geometry(12, 3);
        assert_eq!(asymmetric.left_flank_len, Bp(12));
        assert_eq!(asymmetric.right_flank_len, Bp(3));
        assert_ne!(asymmetric, geometry(3, 12));
    }

    /// A repeat at a contig edge has no flank on that side. Zero is a legal length, not a
    /// missing value — the aligner reports which flanks *anchored*, and a flank that does
    /// not exist cannot.
    #[test]
    fn a_geometry_may_have_no_flank_on_either_side() {
        assert!(geometry(0, 5).fits_reference(Bp(20)));
        assert!(geometry(5, 0).fits_reference(Bp(20)));
        assert_eq!(geometry(0, 0).repeat_len(Bp(20)), Some(Bp(20)));
    }

    /// The empty reference — the degenerate input the trait contract names, whose defined
    /// answer is `Unanchored` rather than an error. With no flanks it is a fit with an
    /// empty repeat; with any flank at all it is not.
    #[test]
    fn an_empty_reference_fits_only_a_geometry_with_no_flanks() {
        assert!(geometry(0, 0).fits_reference(Bp(0)));
        assert_eq!(geometry(0, 0).repeat_len(Bp(0)), Some(Bp(0)));
        assert!(!geometry(1, 0).fits_reference(Bp(0)));
        assert!(!geometry(0, 1).fits_reference(Bp(0)));
    }

    #[test]
    fn geometry_fits_a_reference_exactly_large_enough_for_its_flanks() {
        let flanks = geometry(8, 6);
        assert!(flanks.fits_reference(Bp(20)));
        // Exactly the two flanks and an empty repeat between them: still a fit.
        assert!(flanks.fits_reference(Bp(14)));
        assert_eq!(flanks.repeat_len(Bp(14)), Some(Bp(0)));
        assert!(!flanks.fits_reference(Bp(13)));
    }

    /// **The underflow this precondition exists to prevent.** Production computes the tract
    /// boundary by subtracting the right flank from the reference length; with a geometry
    /// that does not fit, the unchecked form wraps to an enormous number instead of failing.
    /// `repeat_len` returns `None` rather than a saturated zero, because a zero would be
    /// indistinguishable from a genuinely empty tract.
    ///
    /// **Both subtractions are exercised.** A left flank that alone overruns short-circuits
    /// at the first one; the mirror case — left flank fits, right flank overruns — is the
    /// only way to reach the second, so a regression there would otherwise pass.
    #[test]
    fn repeat_len_reports_no_answer_rather_than_underflowing() {
        let flanks = geometry(10, 10);
        assert!(!flanks.fits_reference(Bp(19)));
        assert_eq!(flanks.repeat_len(Bp(19)), None);
        assert_eq!(flanks.repeat_len(Bp(0)), None);

        // First subtraction: the left flank alone overruns.
        let left_overruns = geometry(u64::MAX, 1);
        assert!(!left_overruns.fits_reference(Bp(u64::MAX)));
        assert_eq!(left_overruns.repeat_len(Bp(u64::MAX)), None);

        // Second subtraction: the left flank fits, the right one does not.
        let right_overruns = geometry(1, u64::MAX);
        assert!(!right_overruns.fits_reference(Bp(u64::MAX)));
        assert_eq!(right_overruns.repeat_len(Bp(u64::MAX)), None);
        // ...and it fits exactly when the reference is one base longer than both flanks
        // would need — which `u64` cannot express here, so check the small mirror instead.
        let small = geometry(1, 5);
        assert_eq!(small.repeat_len(Bp(6)), Some(Bp(0)));
        assert_eq!(small.repeat_len(Bp(5)), None);
    }

    /// The two accessors must never disagree about the boundary, since `fits_reference` is
    /// the named precondition and `repeat_len` is the arithmetic implementations actually
    /// perform. Swept across every combination in a small window, including the exact-fit
    /// and one-short cases where a `<=`/`<` slip would live.
    #[test]
    fn fitting_is_exactly_having_a_repeat_length() {
        for left in 0..8u64 {
            for right in 0..8u64 {
                let geometry = geometry(left, right);
                for reference_len in 0..20u64 {
                    let reference_len = Bp(reference_len);
                    assert_eq!(
                        geometry.fits_reference(reference_len),
                        geometry.repeat_len(reference_len).is_some(),
                        "disagreement at flanks {left}/{right}, reference {reference_len:?}"
                    );
                    // When it fits, the three parts must reconstruct the whole.
                    if let Some(repeat) = geometry.repeat_len(reference_len) {
                        assert_eq!(left + right + repeat.get(), reference_len.get());
                    }
                }
            }
        }
    }

    // TODO(Milestone B): the `BestPathAligner` contract states four properties that no
    // test can reach until an implementation exists. The step that writes algorithm 3 owes
    // `align_breaks_ties_toward_match_then_deletion_then_insertion`,
    // `align_assigns_a_junction_gap_to_the_block_on_its_five_prime_side`,
    // `align_returns_unanchored_for_an_empty_reference`, and a debug-build test that the
    // two documented preconditions are asserted.
}
