//! ng's STR admission policy — which detected repeats become STR loci.
//!
//! **This began as a port of [`crate::ssr::catalog::postprocess`], not a rewrite.**
//! The rule set (period scope, score gate, compound-motif drop, bundle drop,
//! minimal trim, copy floor, purity floor, contig-edge drop) is a working, tested
//! implementation — itself a faithful port of GangSTR's `minimal_trim.py` /
//! `remove_bundles.py` — and re-deriving it would be daft (spec §5). So the
//! **decisions are transcribed unchanged**; the shape is ng's.
//!
//! **One rule is deliberately not ported: the flank embed.** Production's
//! `build_loci` stores each tract's bases plus a flank each side, because its
//! consumer genotypes loci without a FASTA open. This module's job is to say what
//! each stretch of the genome *is*, and a payload nobody here reads is not part of
//! that (owner, 2026-07-17 — see [`Locus`]). The **flank requirement stays**: a
//! tract without clean sequence either side is not an STR, which is a question
//! about distances, not bases.
//!
//! ## Why a copy and not a call
//!
//! Step 3 needs admission windowed, 1-based/`u64`, `RepeatInterval`-driven,
//! all-knobs, and handing bundle members back rather than dropping them. That is
//! five changes to `build_loci`, which is not "a small tweak" — and production is
//! frozen (spec Revision 2026-07-16, owner): `src/ssr/` stays exactly as it is,
//! so it remains an **independent yardstick** for the experiments ng exists to
//! run. A production catalog we had rebased to suit ng would be a mirror, not an
//! oracle.
//!
//! The cost is two copies of one policy, which can silently diverge. The
//! mitigation is that divergence is **tested for** rather than prevented:
//! the `differential_vs_production_build_loci_*` tests below drive both
//! implementations from the same intervals and assert identical loci (spec
//! §8.0). What sharing one function used to guarantee by construction is now a
//! test — a weaker guarantee, and the price of a production tree an experiment
//! cannot break.
//!
//! **The differential's blind spot, stated because it is not obvious.** Both
//! sides run the *same* inputs through *transcribed* logic, so an input class no
//! test supplies is not merely uncovered — it is **invisible**: both
//! implementations would be wrong together and the test would stay green. That is
//! why the suite deliberately drives the gates that a single fixed configuration
//! never fires (`min_score`, `min_purity`), a tract `minimal_trim` actually trims,
//! and soft-masked sequence. Each was verified by mutation: break the code, the
//! test fails.
//!
//! ## Divergences from production, and they are the whole list
//!
//! 1. **Coordinates are 1-based inclusive** ([`Locus`]), where production's are
//!    0-based half-open (spec §4). ng is 1-based end to end; production is not
//!    asked to move.
//! 2. **`u64`** coordinates and lengths, where production is `u32` (spec §4).
//! 3. **Input is [`RepeatInterval`]** (ng's scanner) rather than `TrfRecord`
//!    (trf-mod's parse shape). `build_loci` only ever reads `start`/`end`/
//!    `period`/`score` — exactly `RepeatInterval`'s fields — so this costs
//!    nothing. **ng does not depend on trf-mod** (spec Revision).
//! 4. **The pre-filter lives here**, beside the policy it is inseparable from,
//!    rather than in a test file (spec §5.1).
//! 5. **[`Motif`] is ng's too**, against spec §4's expectation — see below.
//!
//! Milestone A2 adds the knobs (period scope, copy floors) and A3 the windowing;
//! both keep the differential green. **Nothing else differs on purpose** — and
//! anything that differs by accident is what the differential exists to catch.
//!
//! ## Order (unchanged from production; `postprocess`'s own numbering)
//!
//! 1. **Period 2..=6** + an early `score` floor. Period-1 homopolymers are
//!    dropped before bundling, so a poly-A run cannot bundle-drop a real STR.
//! 2. **Drop compound-motif loci** — a motif itself internally periodic
//!    (`ATAT` = `(AT)²`).
//! 3. **Drop bundles** — any repeat within `flank_bp` of another goes
//!    with its whole cluster. (ng keeps the *selection* and rejects the
//!    *disposal* in A3 — spec §2.4 — but A1 transcribes it as-is.)
//! 4. **End-trim** to whole-motif boundaries, then the per-period copy floor.
//! 5. **Recompute purity** from the trimmed tract, then the purity floor.
//! 6. **The flank check**: clamp `tract ± flank_bp` at the *contig's* ends and drop
//!    any locus whose flank clamped to zero — a tract with no anchor. Arithmetic
//!    only; production embeds the bases here, and ng does not (above).

use std::fmt;

use crate::ng::tandem_repeat::{PeriodRange, RepeatInterval};
use crate::ng::types::{Bp, Position};

// ---------------------------------------------------------------------
// The motif — ported, though the spec expected it reused
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
/// `ssr::types::Motif` is `pub(crate)`. ng's [`Locus`] is `pub` (the ng-sibling
/// convention) and returns a motif, so reusing it trips rustc's
/// `private_interfaces` lint — a `pub` item leaking a `pub(crate)` type. The
/// three ways out: widen `Motif` in `src/ssr/types.rs` (**touching production —
/// forbidden**); demote ng's whole admission surface to `pub(crate)` (bends ng's
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

// ---------------------------------------------------------------------
// The locus — ng's own, born 1-based
// ---------------------------------------------------------------------

/// A [`Locus`] could not be built because its inputs broke a documented invariant.
///
/// `PartialEq` but not `Eq`: `BadPurity` carries the offending `f32` verbatim
/// (including the `NaN` that may have caused the rejection), and `NaN != NaN`.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum LocusError {
    /// Coordinates are not ordered `1 <= start <= end`. The `1 <=` is ng's:
    /// these are 1-based, so `0` is not a position.
    #[error("locus coordinates out of order: require 1 <= start ({start}) <= end ({end})")]
    BadCoordinates { start: u64, end: u64 },
    /// `purity_fraction` is not a finite value in `[0.0, 1.0]`.
    #[error("purity fraction {purity_fraction} is not finite in [0.0, 1.0]")]
    BadPurity { purity_fraction: f32 },
}

/// One STR locus — where a repeat is, and what it is.
///
/// **It carries no reference bases, and that is the module boundary** (owner,
/// 2026-07-17). This module's job is to say what each stretch of the genome *is*;
/// the bases are already in the reference, and whoever wants them has it open. An
/// earlier version embedded the tract plus a flank each side — production's
/// catalog does, because its consumer genotypes loci without a FASTA in hand — and
/// nothing in ng ever read them. What that payload cost is in [`finish_locus`] and
/// in the walk: a second reference fetch per window, a second margin on top of the
/// scan's, and two panics for the tracts that fell off its edge.
///
/// The **flank is still a typing criterion** — a tract needs clean sequence either
/// side or it is not an STR (spec §2.4, and [`RejectionReason::FlankClamped`]) —
/// but that is a question about *distances*, which coordinates answer. Only the
/// payload is gone.
///
/// ng's port of [`crate::ssr::types::Locus`], now divergent by more than the two
/// differences spec §4 records (**1-based inclusive** and **`u64`**, against
/// production's 0-based/`u32`).
///
/// `chrom` is the contig **name**, matching the project's name-based contig model.
/// The fields are private and the invariant `1 <= start <= end` (plus
/// `purity_fraction` finite ∈ `[0.0, 1.0]`) is enforced by [`Locus::new`], so the
/// accessors are infallible by construction.
#[derive(Clone, Debug, PartialEq)]
pub struct Locus {
    /// Contig name.
    chrom: Box<str>,
    /// Tract start (1-based, inclusive).
    start: u64,
    /// Tract end (1-based, inclusive).
    end: u64,
    /// The repeat unit (verbatim, phase-faithful; `period = motif.period()`).
    motif: Motif,
    /// Fraction of the tract matching a perfect motif tiling, in `[0.0, 1.0]`.
    /// A *degree*, not a flag — see [`Self::is_perfect`].
    purity_fraction: f32,
}

impl Locus {
    /// Build a locus, validating its coordinate and purity invariants.
    ///
    /// The checks are not ceremony: [`admit`] derives these coordinates by
    /// arithmetic on detector output, and an off-by-one here is a *wrong locus,
    /// not a crash*.
    pub fn new(
        chrom: Box<str>,
        start: u64,
        end: u64,
        motif: Motif,
        purity_fraction: f32,
    ) -> Result<Self, LocusError> {
        if !(1 <= start && start <= end) {
            return Err(LocusError::BadCoordinates { start, end });
        }
        if !purity_fraction.is_finite() || !(0.0..=1.0).contains(&purity_fraction) {
            return Err(LocusError::BadPurity { purity_fraction });
        }
        Ok(Self {
            chrom,
            start,
            end,
            motif,
            purity_fraction,
        })
    }

    /// Contig name.
    #[inline]
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    /// Tract start (1-based, inclusive).
    #[inline]
    pub fn start(&self) -> u64 {
        self.start
    }

    /// Tract end (1-based, inclusive).
    #[inline]
    pub fn end(&self) -> u64 {
        self.end
    }

    /// The repeat unit (verbatim, phase-faithful).
    #[inline]
    pub fn motif(&self) -> Motif {
        self.motif
    }

    /// Fraction of the tract matching a perfect motif tiling, in `[0.0, 1.0]`.
    #[inline]
    pub fn purity_fraction(&self) -> f32 {
        self.purity_fraction
    }

    /// The period (motif length, in bases).
    #[inline]
    pub fn period(&self) -> usize {
        self.motif.period()
    }

    /// The tract's length in bases. **Inclusive**, so this is where the `+ 1`
    /// the rebase introduces lives — the one place the arithmetic does not
    /// cancel (spec §4).
    #[inline]
    pub fn tract_len(&self) -> u64 {
        self.end - self.start + 1
    }

    /// Whether the tract is a perfect (uninterrupted) motif tiling.
    #[inline]
    pub fn is_perfect(&self) -> bool {
        self.purity_fraction == 1.0
    }
}

// ---------------------------------------------------------------------
// The knobs
// ---------------------------------------------------------------------

/// Purity floor: a tract matching less than this fraction of a perfect motif
/// tiling is rejected. The catalog's value.
pub const DEFAULT_MIN_PURITY: f32 = 0.8;

/// Score floor: zero, which **admits everything the scanner can emit** — see
/// [`SsrAdmissionParams::min_score`] for why that is deliberate and what it
/// costs. The catalog's value.
pub const DEFAULT_MIN_SCORE: i32 = 0;

/// Flank margin (bp) embedded each side of a tract. The catalog's value, and
/// (spec §2.4) also the bundle radius.
pub const DEFAULT_FLANK_BP: u64 = 50;

/// The narrowest period admitted by default: **2**, excluding period-1
/// homopolymers. The standard GangSTR/HipSTR drop — error-prone for STR
/// genotyping, and not the di/tri/tetra-nucleotide target. Dropping them
/// *before* bundling is load-bearing: it stops a long poly-A/T run from
/// bundle-dropping an adjacent real STR.
pub const DEFAULT_MIN_PERIOD: u8 = 2;

/// The widest period admitted by default: **6**, the microsatellite ceiling.
pub const DEFAULT_MAX_PERIOD: u8 = 6;

/// A period ceiling above [`MAX_MOTIF_LEN`] cannot work — [`Motif`] could not
/// hold the motif. The default must satisfy it; a *configured*
/// [`SsrAdmissionParams::periods`] is checked at [`admit`] instead, since A2
/// makes it a knob (see [`SsrAdmissionParams::periods`]).
const _: () = assert!(
    DEFAULT_MAX_PERIOD as usize <= MAX_MOTIF_LEN,
    "the default period ceiling must fit in a Motif"
);
const _: () = assert!(DEFAULT_MIN_PERIOD <= DEFAULT_MAX_PERIOD);
/// `PeriodRange::new` rejects a zero floor, so `SsrAdmissionParams::default`'s
/// `expect` rests on this. Not idle: §10's period experiment is precisely about
/// moving this floor toward 1, and 0 is one step further.
const _: () = assert!(
    DEFAULT_MIN_PERIOD >= 1,
    "PeriodRange::new rejects a zero period floor"
);

/// The fewest motif copies a tract must have, per period, to be admitted.
///
/// **The "length" half of the period × length routing frontier** the experiments
/// exist to measure (spec §5.2, §10) — requiring *n* copies at period *p* is
/// requiring a tract of *n·p* bases. It is a knob for that reason, and it is what
/// A2 exists for: production hardcodes it in a `const fn`, so the question could
/// not even be *expressed* against that code.
///
/// **The name.** This is the per-period, stricter sibling of
/// [`ScanParams::min_copies`](crate::ng::tandem_repeat::ScanParams::min_copies) —
/// whose own doc calls itself "a permissive floor" and defers to exactly this. So
/// it takes that vocabulary rather than coining a third word for one quantity.
/// It is deliberately **not** called a "copy number": in this crate that phrase
/// means CNV (`relative_copy_number`, `carrier_copy_numbers` in `paralog/`), so
/// for the geneticist reader it would point at the wrong concept entirely.
///
/// **This is the single table.** A1 carried two copies — `copy_number_floor`
/// (admission's) and a nested `copy_floor` inside [`prefilter`] — because
/// production has two: the pre-filter's cleanup has to apply the floor *before*
/// the bundle drop (which runs before admission's own floor), so the number was
/// written twice. They were **duplicates, not rivals**: identical for every
/// period either could reach. A2 folds them into this one value, which both
/// [`prefilter`] and [`admit`] read — so a swept floor moves both, which is the
/// whole point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MinCopies {
    /// `by_period[p - 1]` is the minimum for period `p`, so the array is exactly
    /// GangSTR's `{1:10, 2:5, …}` values in order, with no dead slot to explain.
    by_period: [u32; MAX_MOTIF_LEN],
    /// The minimum for a period this table does not cover.
    ///
    /// Faithful to production's `_ => 3` arm, and it earns its keep: [`prefilter`]
    /// sees **raw detector output**, so an interval wider than [`MAX_MOTIF_LEN`]
    /// (the detector's range is the caller's choice, not admission's) must be
    /// gated on its copy count exactly as production gates it — not indexed out
    /// of bounds, and not silently admitted.
    for_wider_periods: u32,
}

impl MinCopies {
    /// A per-period table: `by_period[p - 1]` is the minimum for period `p`.
    /// `for_wider_periods` covers periods above [`MAX_MOTIF_LEN`].
    pub fn new(by_period: [u32; MAX_MOTIF_LEN], for_wider_periods: u32) -> Self {
        Self {
            by_period,
            for_wider_periods,
        }
    }

    /// The same minimum at every period — the shape for a sweep that varies one
    /// dimension at a time (spec §5.2).
    pub fn uniform(copies: u32) -> Self {
        Self::new([copies; MAX_MOTIF_LEN], copies)
    }

    /// The minimum for `period`.
    ///
    /// Period `0` and periods above [`MAX_MOTIF_LEN`] both fall to
    /// `for_wider_periods` — which is what production's two `_ => 3` arms do for
    /// them, so the total function agrees with the originals everywhere, not just
    /// on the reachable range. (Period 0 is unreachable in practice:
    /// `PeriodRange::new` rejects a zero floor.)
    #[inline]
    pub fn for_period(&self, period: u8) -> u32 {
        period
            .checked_sub(1)
            .and_then(|i| self.by_period.get(usize::from(i)))
            .copied()
            .unwrap_or(self.for_wider_periods)
    }
}

impl Default for MinCopies {
    /// GangSTR's table (`minimal_trim.py`: `thresholds = {1:10, 2:5, 3:4, 4:3,
    /// 5:3, 6:3}`, default 3), which production hardcodes in two places and this
    /// folds into one. Period 1's entry is unreachable at the default period scope
    /// — [`DEFAULT_MIN_PERIOD`] excludes it — and is kept for parity with the
    /// source table, and because §10's experiment may put period 1 back.
    fn default() -> Self {
        //         period:  1  2  3  4  5  6
        Self::new([10, 5, 4, 3, 3, 3], 3)
    }
}

/// Admission's parameters — ng's copy of `ssr::catalog::CatalogParams`, **with
/// every rule a knob** (spec §5).
///
/// A1 transcribed it field-for-field; **A2 makes it a config that can express the
/// question ng exists to ask.** Production's `CatalogParams` cannot: the two
/// dimensions spec §5.2 wants swept — **period** and **length** — are precisely
/// the two it hardcodes (`MIN_PERIOD`/`MAX_PERIOD` consts, a `copy_number_floor`
/// `const fn`), so the routing experiment could not be run against that code at
/// all. Here they are [`periods`](Self::periods) and
/// [`min_copies`](Self::min_copies). A config that cannot express
/// the question is not a config.
///
/// **`bundle_threshold` is gone — it and `flank_bp` were always one number**
/// (spec §2.4). They are two histories, not two designs: `bundle_threshold` is
/// GangSTR's `THRESH=50`, ported from a panel builder with **no flank concept at
/// all**, so over there the number related to nothing; `flank_bp` is ours, the
/// clean sequence a locus must have either side. Both landed on 50 independently, and
/// production's documented `bundle_threshold >= flank_bp` invariant records that
/// coincidence rather than resolving it. **The flank requirement is the
/// primitive and bundle-ness is derived from it** — a repeat is bundled exactly
/// when another sits too close for it to have a clean flank — so one number says
/// it, and the §10 experiment on flank size moves the bundle definition with it
/// for free. `default_matches_the_frozen_catalog_params` pins that the collapse
/// is safe: production ships both at the same value.
///
/// `Default` is the catalog's own values — for §8's comparability **only**, not
/// an endorsement (spec §5.2).
///
/// The fields are `pub` and unvalidated, as production's `CatalogParams` is;
/// [`admit`] `debug_assert`s the contracts that are otherwise only prose.
#[derive(Debug, Clone, PartialEq)]
pub struct SsrAdmissionParams {
    /// The period range admitted, inclusive. Default:
    /// [`DEFAULT_MIN_PERIOD`]`..=`[`DEFAULT_MAX_PERIOD`] (2..=6) — production's
    /// hardcoded consts, now a knob (spec §5).
    ///
    /// **The "period" half of the routing frontier** (spec §5.2, §10). Widening
    /// the floor to 1 puts homopolymers on the STR path; narrowing the ceiling
    /// takes penta/hexamers off it. Both are questions ng exists to answer, and
    /// neither could be asked of production's code.
    ///
    /// **`max` must not exceed [`MAX_MOTIF_LEN`].** [`Motif`] cannot hold a
    /// longer unit, so a wider tract would clear this gate and then be discarded
    /// at `Motif::new` — the knob would *appear* to move and would not. A1 caught
    /// that as a static assert; a knob cannot be checked at compile time, so
    /// [`admit`] `debug_assert`s it and [`finish_locus`] fails loudly rather than
    /// dropping the locus quietly.
    pub periods: PeriodRange,
    /// The fewest motif copies a tract needs, per period — the "length" half of
    /// the routing frontier. Default: [`MinCopies::default`]. Read by **both**
    /// [`prefilter`] and [`admit`], so a swept minimum moves both (spec §5.2).
    pub min_copies: MinCopies,
    /// Purity floor applied after recomputation (a degeneracy cutoff in
    /// `[0, 1]`); imperfect-but-above-floor loci are kept. Default:
    /// [`DEFAULT_MIN_PURITY`].
    pub min_purity: f32,
    /// Early accept-gate on the detector's `score`: records **below** it are
    /// dropped (`score >= min_score` is admitted).
    ///
    /// **At [`DEFAULT_MIN_SCORE`] (`0`) this gate never fires, deliberately — and
    /// that is a trap worth knowing** (spec §5c). The catalog can afford a `0`
    /// floor because the external `trf-mod` binary applies its own `-s 30`
    /// upstream. **ng has no trf-mod** (spec Revision), and
    /// `RepeatInterval::score` is a Ruzzo–Tompa segment total, **not a TRF
    /// score** — the two are not on the same scale, so a threshold carried across
    /// from TRF would not mean there what it meant here.
    ///
    /// Precisely: Ruzzo–Tompa only emits **positive-scoring** segments (`score =
    /// r - l > 0`, `tandem_repeat.rs`), so `score >= 0` cannot reject anything the
    /// scanner produces. It is a no-op *for scanner output* — not literally "no
    /// gate": a negative score would be dropped, and nothing emits one. The
    /// copy-number and purity floors are the real gates, and the 16/16 scanner
    /// parity ran exactly this way. Inherited knowingly, not by accident;
    /// `admit_at_default_never_rejects_a_score_the_scanner_can_emit` pins it.
    pub min_score: i32,
    /// How much clean sequence (bp) a tract must have either side to be a locus —
    /// **and, since A2, the bundle radius too** (spec §2.4; see the type's docs).
    /// Default: [`DEFAULT_FLANK_BP`].
    ///
    /// One number, two jobs, because they are the same requirement seen twice: a
    /// locus needs `flank_bp` of clean sequence each side to anchor reads, so a
    /// repeat with another repeat inside `flank_bp` cannot have one — which is
    /// what makes it a bundle member rather than a locus. The contig's own end is
    /// the third face of it ([`RejectionReason::FlankClamped`]).
    ///
    /// **It is a distance, and nothing here reads the bases it measures** — the
    /// flank is a test a tract passes, not something a [`Locus`] carries (spec
    /// §1.2).
    pub flank_bp: u64,
}

impl Default for SsrAdmissionParams {
    /// The catalog's pinned Stage-0 defaults, field for field.
    ///
    /// **This must equal `CatalogParams::default()`**, and that is not decoration
    /// — spec §8.1's `.cat` parity oracle compares ng's walk against a catalog
    /// built at those values, and the comparison is meaningless if the two drift.
    /// Production is frozen, so it can only drift by an ng-side edit; that edit
    /// should be deliberate, so `default_matches_the_frozen_catalog_params` pins
    /// every field.
    ///
    /// Note `min_score`'s value is a no-op gate — see the field's own doc.
    fn default() -> Self {
        Self {
            // Infallible at the defaults: the const asserts above prove
            // `1 <= DEFAULT_MIN_PERIOD <= DEFAULT_MAX_PERIOD`, which is exactly
            // what `PeriodRange::new` rejects.
            periods: PeriodRange::new(DEFAULT_MIN_PERIOD, DEFAULT_MAX_PERIOD)
                .expect("the default period range is const-asserted valid"),
            min_copies: MinCopies::default(),
            min_purity: DEFAULT_MIN_PURITY,
            min_score: DEFAULT_MIN_SCORE,
            flank_bp: DEFAULT_FLANK_BP,
        }
    }
}

// ---------------------------------------------------------------------
// The candidate — retired at B2
// ---------------------------------------------------------------------
//
// A1 introduced a private `Candidate` — "exactly the four fields `build_loci`
// reads off a `TrfRecord`, widened to ng's `u64`" — because `RepeatInterval` was
// `u32` and the port needed one widening site rather than a `u64::from` scattered
// through the policy. Its doc said B2 would make the conversion an identity but
// that the type would "still earn its keep as the 'what admission actually needs'
// statement".
//
// **It did not.** With `RepeatInterval` at `u64` (spec §4), `Candidate` was
// field-for-field identical to it, reachable only through two `From` impls that
// were both no-ops — and one of them carried a `try_into().expect()` and a
// `fallible_impl_from` exemption for a narrowing that can no longer narrow. That
// is ceremony asserting a distinction the types no longer have, so the policy now
// works on `RepeatInterval` directly. `Admitted::bundled` hands the caller's own
// intervals back with no conversion at all, which is what "handed back verbatim"
// should have meant all along.
//
// **Coordinates are 0-based half-open throughout the policy** — the detector's and
// the slice's own space. The conversion to ng's 1-based inclusive happens exactly
// once, at [`Locus`] construction in [`finish_locus`], which is the one place the
// arithmetic does not cancel (spec §4).

// ---------------------------------------------------------------------
// The pre-filter
// ---------------------------------------------------------------------

/// Clean raw scanner intervals before [`admit`] — **not optional** (spec §5b).
///
/// A validated finding, not a guess (`scanner_parity`, Milestone D of the
/// scanner plan). trf-mod hands `build_loci` a *clean* candidate set:
/// statistically significant repeats, redundancy already eliminated. Our scanner
/// is deliberately permissive (`min_copies = 2`), so it also emits low-copy noise
/// in aperiodic sequence **and every period-multiple of every real tract**. Fed
/// straight in, that noise trips the bundle drop — which runs *before* the copy
/// floor — and cascades the real loci away (the first parity run scored 0/16).
/// With this filter the scanner reproduces the golden catalog at **16/16**.
///
/// Two cleanups, mirroring what trf-mod bakes in: the per-period **copy floor**,
/// then **period-multiple redundancy elimination** (trf-mod's `IsRedundant` — a
/// higher-period interval overlapping a divisor-period one is the same tract).
///
/// Ported from the test-only `catalog_prefilter`
/// (`ng::scanner_parity` — `#[cfg(test)]`, so not linkable), which had no
/// production home; here it sits
/// beside the policy whose ordering makes it necessary (spec §5.1).
///
/// **It reads the same [`MinCopies`] as [`admit`] — since A2, there is one
/// table.** A1 carried production's two copies verbatim (see
/// [`MinCopies`]); folding them is what makes a swept floor actually move
/// the whole policy, rather than half of it. `p.periods.min()` likewise replaces
/// the bare `2` this was ported with.
///
/// **The period *ceiling* is deliberately not applied here — only the floor is —
/// and the reason is transcription fidelity, nothing cleverer.** Production's
/// pre-filter gates `period >= 2` and leaves the ceiling to `build_loci`, so this
/// does the same. Since no test can compare it (above), matching the original by
/// inspection is the whole of the evidence, and staying close to it is the point.
///
/// *Two arguments an earlier draft gave here were wrong, recorded so they are not
/// re-invented: (1) that applying the ceiling early "would change which survivors
/// eliminate which" — it cannot, because `floored` sorts ascending by period, so
/// an out-of-scope interval only ever sorts after in-scope ones and can only be a
/// victim of redundancy elimination, never an eliminator; every interval it could
/// eliminate is its own multiple, which [`admit`]'s ceiling drops regardless. And
/// (2) that "the differential would catch it otherwise" — it cannot, and the
/// paragraph above says so.*
///
/// **The two floors run at different points, and that is the whole of this
/// function's subtlety** (fixed 2026-07-17; see
/// `a_homopolymer_does_not_survive_as_a_period_two_repeat`). The **copy** floor runs
/// first, the **period** floor runs *last* — after redundancy elimination — because
/// **period 1 divides every period**, which cuts both ways:
///
/// - A homopolymer tiles under every motif length: `AAAA…` is a perfect period-2
///   `AA`, period-3 `AAA`, period-5 `AAAAA` repeat, and the scanner emits all of
///   them. Only the period-1 interval can eliminate those aliases. Drop period 1
///   *before* redundancy — as this did until 2026-07-17 — and the aliases outlive
///   their only eliminator: a 30 bp poly-A came out as **three** "repeats" (periods
///   2, 3 and 5, identical spans), which then fed coverage, the satellite cap and
///   the rejection counts as if three tandem repeats were there. **The period floor
///   filtered out the period-1 *label*, not the homopolymer**, so
///   `periods.min() == 2` did not mean what it said.
/// - The reason it was ordered that way is real, and the **copy** floor is what
///   answers it: a surviving period-1 interval eliminates every real STR it
///   overlaps. Aperiodic sequence is full of 2 bp period-1 specks (`GG`, `TT`), and
///   those must never be eliminators. [`MinCopies`]'s period-1 entry is **10** — so
///   they are dropped here, by the copy floor, and only a genuine homopolymer run
///   survives to eliminate its own aliases. That entry was dead weight until this
///   ordering (nothing reached it — see
///   `the_two_copy_floor_tables_agree_on_every_reachable_period`); it is now exactly
///   what separates a poly-A from a speck.
///
/// **The scanner must therefore scan period 1 even when nothing admits it**
/// (`TypedRegionConfig::periods`, which is why it is a separate knob from
/// `SsrAdmissionParams::periods`): an eliminator that was never detected cannot
/// eliminate. Not seeing a homopolymer and dropping one are *not* the same, and
/// this ordering is what makes that sentence true.
///
/// The result is that **the period range means what it says**: at `2..=6` a
/// homopolymer contributes nothing at all, and at `1..=6` it is admitted as a
/// period-1 tract under its own motif. Which of those is wanted is the caller's
/// parameter, not this function's opinion.
///
/// **Malformed intervals are dropped, not trusted.** `RepeatInterval`'s fields are
/// `pub` and carry no ordering invariant, and this fn is `pub` — unlike the
/// private `#[cfg(test)]` helper it was ported from, which only ever saw its own
/// scanner's output. An `end < start` interval would underflow `end - start`:
/// debug panic, or release wrap-to-huge — and a wrapped copy count sails straight
/// through the very floor this filter exists to apply. So the ordering is checked
/// rather than assumed. Well-formed intervals are unaffected, so the port stays
/// faithful.
pub fn prefilter(intervals: &[RepeatInterval], params: &SsrAdmissionParams) -> Vec<RepeatInterval> {
    // The COPY floor only. The period floor is deliberately not here — it runs
    // last, once redundancy elimination has had the period-1 intervals it needs
    // (see the fn docs). The copy floor is what stops a 2 bp `GG` speck being an
    // eliminator: `MinCopies`' period-1 entry is 10.
    let mut floored: Vec<RepeatInterval> = intervals
        .iter()
        .copied()
        .filter(|iv| {
            // `iv.end > iv.start` first: it guards the subtraction below (see the
            // fn docs), and it is what `admit` independently requires anyway.
            iv.end > iv.start
                && (iv.end - iv.start) / u64::from(iv.period)
                    >= u64::from(params.min_copies.for_period(iv.period))
        })
        .collect();
    // Process low periods first so a fundamental tract is kept and its multiples
    // dropped. A homopolymer sorts first of all, which is what lets it take its own
    // period-2/3/4/5/6 aliases with it.
    floored.sort_by_key(|iv| (iv.period, iv.start));
    let mut kept: Vec<RepeatInterval> = Vec::new();
    for iv in floored {
        let redundant = kept.iter().any(|a| {
            a.period < iv.period
                && iv.period % a.period == 0
                && iv.start < a.end
                && a.start < iv.end
        });
        if !redundant {
            kept.push(iv);
        }
    }
    // The PERIOD floor, last: an out-of-scope interval has now done its
    // eliminating. Dropping a homopolymer here drops the homopolymer, not merely
    // its period-1 label — its aliases went with it above.
    kept.retain(|iv| iv.period >= params.periods.min());
    kept
}

// ---------------------------------------------------------------------
// Admission
// ---------------------------------------------------------------------

/// What [`admit`] made of one window's candidates.
///
/// **`bundled` is the A3 change of substance** (spec §2.4, §6a). Production's
/// `build_loci` *deletes* the members of a too-close cluster: they are neither
/// loci nor anything else, and the bases they cover become a hole nobody accounts
/// for. ng keeps the **selection** and rejects the **disposal** — the same records
/// are set aside, but handed back so the walk can route them as an `SsrBundle`
/// region and a consumer can decide later, with the evidence in hand, rather than
/// having the decision taken here by deletion.
///
/// Keeping selection and disposal separable is also what preserves spec §8's
/// oracle: the *selection* is unchanged, so ng's locus set stays comparable with
/// production's (a bundle member is not a locus in either system).
#[derive(Debug, Clone, PartialEq)]
pub struct Admitted {
    /// The admitted STR loci, start-sorted. Exactly what `build_loci` returns.
    pub loci: Vec<Locus>,
    /// The cluster members `build_loci` silently drops — every repeat that
    /// cleared the period/score/compound gates but sits within `flank_bp` of
    /// another, so no clean flank can be built around it (spec §2.4).
    ///
    /// **Coordinates are the input's**, i.e. offsets into `bases`, not genomic —
    /// these are the caller's own `RepeatInterval`s handed straight back. The walk
    /// re-bases them when it emits the region's hull, the same way it re-bases
    /// everything else it reads out of a window.
    pub bundled: Vec<RepeatInterval>,
    /// Every repeat admission turned down, **with the reason** — the breakdown spec §3.1
    /// needs and that one total cannot give (a purity rejection and a copy-floor one are
    /// exactly the distinction spec §10's routing question turns on).
    ///
    /// **Per record, not a tally, and that is what makes the walk's count
    /// window-invariant.** `admit` runs over a window's whole fetched slice, margins
    /// included, so the same repeat is rejected again by every window that can see it. An
    /// aggregate here would be counted once per window and the total would depend on
    /// `window_bp` — a memory knob changing a reported number. Handing the records back
    /// lets the walk attribute each to exactly one core, the same way it attributes
    /// everything else.
    ///
    /// Coordinates are the input's, like [`Self::bundled`]'s.
    pub rejected: Vec<(RepeatInterval, RejectionReason)>,
}

/// Why admission turned a repeat down (spec §3.1).
///
/// **No `str_bundle`**: since spec §2.4 that is a *route*, not a rejection —
/// [`Admitted::bundled`] carries those, and they are not here.
///
/// **The scope and score gates are not here either**, and that is a real gap rather than
/// an oversight: a repeat outside `periods` or under `min_score` is not *rejected*, it is
/// out of the question being asked — and `prefilter` has already dropped most of it before
/// `admit` sees it (spec §5b). The five below are the gates that turn down a repeat that
/// was, on the face of it, a candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RejectionReason {
    /// Fewer whole motif copies than the per-period floor ([`MinCopies`]).
    CopyFloor,
    /// Below the purity floor once re-measured against a perfect motif tiling.
    Purity,
    /// The motif is itself a repeat (`ATAT` is `AT` twice), so the period is a lie.
    Compound,
    /// No whole-motif boundaries to trim to — nothing here tiles.
    NoCleanTrim,
    /// A flank clamped to nothing at a **contig** end: the tract abuts base 1 or the last
    /// base, so there is no unique sequence to anchor reads against (spec §2.6).
    FlankClamped,
}

/// Repeat bp turned down, by reason — spec §3.1's breakdown of
/// `TypedRegionCounts::repeat_bp_with_no_locus`.
///
/// **In bp, not per repeat** (spec §3.1): admission trims every survivor, so a repeat that
/// admits one locus and sheds 200 bp contributes nothing to a per-repeat counter. The bp
/// charged is the repeat's **detected** length, before trimming.
///
/// **It does not partition `repeat_bp_with_no_locus`, and must not be read as if it did.**
/// Two rejected repeats can overlap — the scanner is deliberately permissive — so their bp
/// are both charged; and bases with no locus for a reason that is not a *rejection*
/// (bundled, capped as a satellite, out of scope) are in the total and not here. It is a
/// diagnosis of admission's gates, not an account of the genome.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RejectionCounts {
    pub copy_floor: u64,
    pub purity: u64,
    pub compound: u64,
    pub no_clean_trim: u64,
    pub flank_clamped: u64,
}

impl RejectionCounts {
    /// Charge `bp` to `reason`.
    pub fn add(&mut self, reason: RejectionReason, bp: u64) {
        let slot = match reason {
            RejectionReason::CopyFloor => &mut self.copy_floor,
            RejectionReason::Purity => &mut self.purity,
            RejectionReason::Compound => &mut self.compound,
            RejectionReason::NoCleanTrim => &mut self.no_clean_trim,
            RejectionReason::FlankClamped => &mut self.flank_clamped,
        };
        *slot += bp;
    }

    /// Total bp turned down, over every reason.
    pub fn total(&self) -> u64 {
        self.copy_floor + self.purity + self.compound + self.no_clean_trim + self.flank_clamped
    }
}

/// Admit detected repeats into start-sorted STR loci — ng's `build_loci`,
/// **windowed**.
///
/// `chrom` is the contig name. `bases` is a slice of that contig (any case — the
/// tract and motif are upper-cased here for case-stable identity); `recs` are
/// offsets into `bases`. **No base outside a tract is read**, so the slice needs no
/// margin of its own — hand it exactly what was scanned (spec §2.6). Feed it
/// [`prefilter`]ed intervals, never raw scanner output (spec §5b).
///
/// # The two facts a bare slice cannot tell you (spec §2.6)
///
/// - **`bases_start`** — the 1-based genomic coordinate of `bases[0]`.
/// - **`contig_len`** — how long the contig *really* is.
///
/// Production takes neither, because it only ever sees whole contigs, so
/// `contig_seq.len()` silently plays two different roles: the bound on what can be
/// read, **and** the position of the chromosome's end. Hand that code a 100 kb
/// window and it believes the chromosome ends at your window edge — so it clamps
/// every flank there and **throws away every locus within `flank_bp` of every
/// boundary, a different set for every `window_bp`**. It makes no noise at all
/// when it does. Splitting the two roles is the fix, and it is the whole reason
/// this signature exists.
///
/// **Whole-contig is the degenerate case**: `bases_start = 1`,
/// `contig_len = bases.len()`. Then the two roles coincide again and the output is
/// production's, unchanged — which is what keeps
/// `differential_vs_production_build_loci_*` a valid oracle for the windowed code
/// (spec §5a). Every window-only path is therefore proved *by* the degenerate
/// case staying green, plus `admit_keeps_a_locus_the_window_edge_would_have_eaten`.
///
/// # Panics
///
/// If `bases` cannot supply a survivor's flanks — i.e. the caller windowed without
/// leaving a margin (spec §2.6 fetches core ± `max_repeat_len`, and 1 kb is far
/// more than the 50 bp of flank). Loudly, because the alternative is dropping the
/// locus quietly, which is the very bug this signature exists to kill.
///
/// Output coordinates are **1-based inclusive genomic** ([`Locus`]); the input's
/// are 0-based half-open offsets into `bases`.
pub fn admit(
    recs: Vec<RepeatInterval>,
    chrom: &str,
    bases: &[u8],
    bases_start: Position,
    contig_len: Bp,
    p: &SsrAdmissionParams,
) -> Admitted {
    // **Newtypes, not two bare `u64`s** (A3 review). They were adjacent, same-typed,
    // and silently transposable — and the guard below catches a transposition only
    // sometimes. `Position` and `Bp` are different types, so the compiler catches it
    // always. This is the cheap half of that finding; the other half is *provenance*
    // (below), which no signature can fix.
    let bases_start = bases_start.get();
    let contig_len = contig_len.get();
    assert!(
        bases_start >= 1,
        "bases_start is 1-based; got {bases_start}"
    );
    // **A real `assert!`, for the same reason the period ceiling below is one:
    // sweeps run in release.** An understated `contig_len` clamps `ref_end` early,
    // and the `ref_end == tract_end` test still passes — so a locus comes out with
    // a **silently truncated right flank**, or vanishes. That is spec §2.6's bug
    // reappearing inside the code written to kill it, and §10's flank and routing
    // sweeps drive this path in release. Once per `admit` call.
    //
    // **It is inherently incomplete, and that is worth stating rather than
    // trusting.** `<=` admits equality, so a caller passing *the window's own end*
    // as `contig_len` — production's exact mistake — is arithmetically legal and
    // no check here can catch it. The other half needs **provenance**: a
    // `contig_len` read from the reference's contig table, not derived from
    // whatever slice is in hand.
    //
    // `WindowedRefSeq::contigs()` (B3) is that table, and reading from it is the
    // **walk's** obligation (Milestone D) — it is the only party that holds both
    // the reference and the window. Nothing `admit` can be handed proves where the
    // number came from, so this is a contract stated here and discharged there,
    // not a check that was forgotten.
    assert!(
        bases_start + bases.len() as u64 - 1 <= contig_len,
        "the window [{bases_start}, {}] runs past the contig end ({contig_len}): \
         contig_len must be the CONTIG's length, never the window's",
        bases_start + bases.len() as u64 - 1
    );
    // `flank_bp = 0` admits nothing at all, from any input: every tract's own flank
    // test (`ref_start == tract_start`) fires, so every locus is dropped and
    // nothing is even bundled — an empty result, no error, from a knob spec §10
    // plans to sweep. The M7 shape again: the knob appears to move, and instead
    // switches the whole step off.
    assert!(
        p.flank_bp >= 1,
        "flank_bp must be at least 1: at 0 every locus fails its own flank test \
         and admit silently returns nothing"
    );
    // A NaN `min_purity` bites quietly: every `purity < p.min_purity` comparison
    // below is then false, so the purity gate passes **everything** rather than
    // failing. `assert!` for the same reason as its neighbours — the purity floor
    // is a swept knob (spec §10), sweeps run in release, and a release sweep that
    // silently admitted every tract would report a finding about the data that is
    // really a fact about a `NaN`.
    assert!(
        p.min_purity.is_finite() && (0.0..=1.0).contains(&p.min_purity),
        "min_purity must be finite in [0, 1], got {}",
        p.min_purity
    );
    // A1 could assert this at compile time; A2 made the ceiling a knob, so it
    // moves here. Above `MAX_MOTIF_LEN` a tract clears the period gate and is then
    // discarded at `Motif::new` — the knob appears to move and does not.
    //
    // **A real `assert!`, not a `debug_assert!`** — `PeriodRange::new` bounds
    // neither end, so `PeriodRange::new(2, 7)` is legal and this is reachable. The
    // point of the knob is to be **swept by the §5.2 routing experiment, and
    // sweeps run in release**: a debug-only check would let that experiment record
    // "period 7 admits nothing" when the code never tried, which is a wrong
    // *result*, not a missing panic. Once per `admit` call — free next to the scan
    // it guards.
    assert!(
        p.periods.max() as usize <= MAX_MOTIF_LEN,
        "period ceiling {} exceeds MAX_MOTIF_LEN ({MAX_MOTIF_LEN}): wider tracts would pass \
         this gate and then be silently dropped at Motif::new",
        p.periods.max()
    );

    // 1. scope + score gate, then 2. compound-motif drop. Both need the tract,
    //    so we slice the (upper-cased) prefix motif here and filter on it.
    //
    //    `r.end <= bases.len()` stays a SLICE bound, not a contig bound: `r` is an
    //    offset into `bases`, and this guards the slicing two lines down. The
    //    contig end is a different question, and it is asked in `finish_locus`.
    //
    //    Two `filter`s until E1e; one loop now, because a rejection has to be *recorded*
    //    and not merely fall through. Same gates, same order, same outcome — pinned by
    //    the differential against production's `build_loci`.
    let mut rejected: Vec<(RepeatInterval, RejectionReason)> = Vec::new();
    let mut kept: Vec<RepeatInterval> = Vec::new();
    for r in recs {
        // Scope, score, and malformed. **Not recorded as rejections** — a repeat outside
        // the period range or under the score floor is out of the question being asked,
        // not turned down by it (`RejectionReason`).
        if r.period < p.periods.min()
            || r.period > p.periods.max()
            || r.score < p.min_score
            || r.end <= r.start
            || (r.end as usize) > bases.len()
        {
            continue;
        }
        let period = r.period as usize;
        let tract = &bases[r.start as usize..r.end as usize];
        // Need at least one full motif to form the prefix. Charged to the copy floor:
        // a tract that cannot hold even one whole copy is that floor at its extreme, and
        // it is the same fact the floor exists to measure.
        if tract.len() < period {
            rejected.push((r, RejectionReason::CopyFloor));
            continue;
        }
        let motif = upper(&tract[..period]);
        if is_compound(&motif) {
            rejected.push((r, RejectionReason::Compound));
            continue;
        }
        kept.push(r);
    }

    // 3. drop bundles (on the raw, pre-trim coordinates). Records must be
    //    start-sorted for the streaming clustering.
    kept.sort_by_key(|r| (r.start, r.end));
    // `flank_bp` IS the bundle radius (spec §2.4) — see `SsrAdmissionParams`.
    // A3: the cluster members come back rather than being deleted (`Admitted`).
    let (kept, bundled) = split_bundles(kept, p.flank_bp);

    // 4-5. per record: end-trim + copy floor, recompute purity + floor, embed.
    let mut loci = Vec::with_capacity(kept.len());
    for r in &kept {
        match finish_locus(r, chrom, bases, bases_start, contig_len, p) {
            Ok(locus) => loci.push(locus),
            Err(reason) => rejected.push((*r, reason)),
        }
    }
    Admitted {
        loci,
        // No conversion: since B2 these ARE the caller's own intervals, which is
        // what `Admitted::bundled`'s "handed back verbatim" always claimed.
        bundled,
        rejected,
    }
}

/// Steps 4-6 for one record: end-trim, copy-number floor, motif, purity floor, the
/// flank check. `Err(reason)` if the record fails any gate.
///
/// **Every `None` became an `Err(reason)` at E1e**, and nothing else changed. The gates
/// and their order are exactly as before; what is new is that the caller can now tell a
/// purity rejection from a copy-floor one, which is the distinction spec §10's routing
/// question turns on and which one `None` could never carry (spec §3.1).
///
/// **The one place the coordinate base changes**, and (since A3) the one place
/// the *contig* end is asked about. Everything above works in the detector's
/// 0-based half-open offsets into `bases` — which is also the slice's own space,
/// so the arithmetic stays production's unchanged. The 1-based inclusive genomic
/// [`Locus`] is built at the bottom, once. Spec §4 predicted exactly this: "the
/// arithmetic cancels … exactly one site does not".
fn finish_locus(
    r: &RepeatInterval,
    chrom: &str,
    bases: &[u8],
    bases_start: u64,
    contig_len: u64,
    p: &SsrAdmissionParams,
) -> Result<Locus, RejectionReason> {
    let period = r.period as usize;
    let raw_tract = upper(&bases[r.start as usize..r.end as usize]);
    // A tract too short to hold one motif: the copy floor at its extreme (see `admit`,
    // which charges the same case the same way). Unreachable from `admit` — its own gate
    // above rejects `tract.len() < period` first — but this fn is not `admit`'s alone to
    // reason about.
    let motif_bytes = raw_tract
        .get(..period)
        .ok_or(RejectionReason::CopyFloor)?
        .to_vec();

    // End-trim to clean whole-motif boundaries (GangSTR minimal_trim).
    let (st, en) = minimal_trim(&raw_tract, &motif_bytes).ok_or(RejectionReason::NoCleanTrim)?;
    let new_start = r.start + st as u64;
    let new_end = r.start + en as u64;
    let trimmed = &raw_tract[st..en];

    // Copy-number floor — GangSTR computes copies from the ORIGINAL span
    // (integer division), as an accept-gate after trimming.
    let ref_copy = (r.end - r.start) / u64::from(r.period);
    if ref_copy < u64::from(p.min_copies.for_period(r.period)) {
        return Err(RejectionReason::CopyFloor);
    }

    // After minimal_trim the tract starts on a motif boundary, so
    // `trimmed[..period] == motif_bytes`; use it as the phase-faithful motif.
    //
    // A1 wrote `.ok()?` here, faithfully to production, where it was unreachable:
    // the period ceiling was a `const` a static assert held at `<= MAX_MOTIF_LEN`.
    // **A2 made that ceiling a knob and the arm reachable** — set
    // `periods.max() > MAX_MOTIF_LEN` and every widest-period tract clears the
    // gate, fails here, and vanishes through the same rejection as a legitimate
    // policy one. That is M7's silent failure exactly: the knob appears to
    // move and does not. So it fails loudly instead (release behaviour unchanged —
    // still a rejection; `admit` asserts the same contract up front).
    //
    // E1e: charged to the copy floor when it does happen, for the reason the short-tract
    // case above is — a motif that cannot be formed is a tract with no whole copies. It
    // is **unreachable** in any case: `admit` asserts `periods.max() <= MAX_MOTIF_LEN`
    // (a real assert, in release) and gates `period >= periods.min() >= 1`, so
    // `motif_bytes.len()` is in `1..=MAX_MOTIF_LEN` and `BadLength` is the only way
    // `Motif::new` fails. The `debug_assert` is the message that matters if it ever fires.
    let motif = match Motif::new(&motif_bytes) {
        Ok(m) => m,
        Err(e) => {
            debug_assert!(
                false,
                "period {} admitted but has no valid Motif ({e}) — is \
                 SsrAdmissionParams::periods.max() above MAX_MOTIF_LEN ({MAX_MOTIF_LEN})?",
                r.period
            );
            return Err(RejectionReason::CopyFloor);
        }
    };

    // Recompute purity from the trimmed tract vs a perfect motif tiling.
    let purity = recompute_purity(trimmed, &motif_bytes);
    if purity < p.min_purity {
        return Err(RejectionReason::Purity);
    }

    // **The window-boundary fix (spec §2.6), and the whole of it.**
    //
    // Move to genomic 1-based inclusive first, because the flank clamp is a
    // question about the *contig*, and `bases` cannot answer it. `[s, e)` at
    // offset `s` is `[bases_start + s, bases_start + e - 1]`.
    let tract_start = bases_start + new_start;
    let tract_end = bases_start + new_end - 1;

    // **The flank is a typing criterion, and it is a question about distances.**
    // A tract needs clean sequence either side or it is not an STR (spec §2.4) —
    // and whether that sequence exists is answered by coordinates and the contig's
    // length, not by reading any of it. Clamp at the CONTIG's ends — 1 and
    // `contig_len` — never at the slice's. Production writes
    // `.min(contig_seq.len())` because for it the two are the same thing; for a
    // window they are not, and believing the chromosome stops at the window edge is
    // what silently eats every locus near every boundary.
    let ref_start = tract_start.saturating_sub(p.flank_bp).max(1);
    let ref_end = (tract_end + p.flank_bp).min(contig_len);

    // Drop a locus whose flank clamped to nothing on either side — a tract
    // abutting base 1 of the contig (empty left flank) or ending on its last base
    // (empty right flank). Nothing can be anchored against a zero-length flank, so
    // the tract is not an STR. These are genuinely about the contig, so a window
    // edge no longer impersonates one.
    if ref_start == tract_start || ref_end == tract_end {
        return Err(RejectionReason::FlankClamped);
    }

    // **No bases are read past the tract**, which is what lets the caller hand over
    // exactly the slice it scanned. The flank test above is arithmetic, and the
    // tract's own bytes are in the slice by construction — the interval came from
    // scanning it. Until 2026-07-17 this read `[tract - flank_bp, tract + flank_bp]`
    // to embed in the locus, and a tract at the slice's own edge has no such margin:
    // hence two panics here, a second margin on the caller's fetch, and a second
    // reference read per window. All of it served a payload nothing consumed.

    match Locus::new(
        chrom.to_string().into_boxed_str(),
        tract_start,
        tract_end,
        motif,
        purity,
    ) {
        Ok(locus) => Ok(locus),
        // Unreachable: every gate above guarantees the invariants (`minimal_trim`
        // gives `st < en`, so `tract_start <= tract_end`; `recompute_purity` returns
        // `[0, 1]`).
        //
        // Which is exactly why the verdict is worth keeping: this error can now only
        // fire if **the arithmetic is wrong**, and that is otherwise silent (a bad
        // locus would leave through the same rejection as a routine policy one).
        // `debug_assert` makes the differential — which runs in debug — fail loudly
        // instead of quietly disagreeing.
        //
        // The reason it carries in release is `NoCleanTrim`: of the five, that is the
        // one that means "this tract did not turn out to be a well-formed repeat", which
        // is the closest true thing to say about a locus whose invariants did not hold.
        // A count is not the point here; the `debug_assert` is.
        Err(e) => {
            debug_assert!(
                false,
                "admit built an invalid locus, so the arithmetic is wrong: {e} \
                 (tract [{tract_start}, {tract_end}], window starts at {bases_start}, \
                  contig is {contig_len} long)"
            );
            Err(RejectionReason::NoCleanTrim)
        }
    }
}

/// ASCII upper-case copy of `bytes`.
fn upper(bytes: &[u8]) -> Vec<u8> {
    bytes.iter().map(|b| b.to_ascii_uppercase()).collect()
}

/// Count non-overlapping greedy occurrences of `motif` in `repeat` (GangSTR
/// `count_motif`): scan left-to-right, advancing by `motif.len()` on a match and
/// by 1 otherwise.
fn count_motif(repeat: &[u8], motif: &[u8]) -> usize {
    let m = motif.len();
    let mut c = 0;
    let mut s = 0;
    while s < repeat.len() {
        if s + m <= repeat.len() && &repeat[s..s + m] == motif {
            c += 1;
            s += m;
        } else {
            s += 1;
        }
    }
    c
}

/// `true` if `motif` is itself internally periodic — a non-fundamental period
/// such as `ATAT = (AT)²` (GangSTR `is_compound`, threshold 0.8). A motif whose
/// shorter prefix `sub` tiles more than 80% of it is compound.
fn is_compound(motif: &[u8]) -> bool {
    let l = motif.len();
    const THRESHOLD: f64 = 0.8;
    for i in 1..=(l / 2) {
        let sub = &motif[..i];
        if (count_motif(motif, sub) * i) as f64 > l as f64 * THRESHOLD {
            return true;
        }
    }
    false
}

/// End-trim a tract to clean whole-motif boundaries (GangSTR `minimal_trim`):
/// find the smallest `start_offset` where two motif copies (`motif*2`) appear,
/// and the largest `end_offset` (exclusive) where they appear before it. Returns
/// `(start_offset, end_offset)` into `rep`, or `None` if the tract has no clean
/// motif boundary within `min(3·|motif|, |rep|/2)` of each end.
///
/// Offsets, not coordinates — no rebase applies.
fn minimal_trim(rep: &[u8], motif: &[u8]) -> Option<(usize, usize)> {
    let m = motif.len();
    let ll = m * 2;
    let max_trim_len = (m * 3).min(rep.len() / 2);

    let two_motifs = |slice: &[u8]| -> bool {
        slice.len() == ll && &slice[..m] == motif && &slice[m..] == motif
    };

    // Smallest start offset whose next 2·m bytes are motif*2. GangSTR checks the
    // bound first and bails (not continues) — replicate.
    let mut start = None;
    for so in 0..=max_trim_len {
        if so + ll >= rep.len() {
            return None;
        }
        if two_motifs(&rep[so..so + ll]) {
            start = Some(so);
            break;
        }
    }
    let st = start?;

    // Largest end offset (exclusive) whose preceding 2·m bytes are motif*2.
    let mut end = None;
    for eo in (rep.len().saturating_sub(max_trim_len)..=rep.len()).rev() {
        if eo < ll {
            return None;
        }
        if two_motifs(&rep[eo - ll..eo]) {
            end = Some(eo);
            break;
        }
    }
    let en = end?;
    if st >= en {
        return None;
    }
    Some((st, en))
}

/// Fraction of `tract` matching a perfect tiling of `motif` from phase 0. With
/// `motif == tract[..period]`, the first repeat always matches, so a perfect
/// tract scores 1.0 and interruptions lower it proportionally.
fn recompute_purity(tract: &[u8], motif: &[u8]) -> f32 {
    if tract.is_empty() || motif.is_empty() {
        return 0.0;
    }
    let m = motif.len();
    let matches = tract
        .iter()
        .enumerate()
        .filter(|&(i, &b)| b == motif[i % m])
        .count();
    matches as f32 / tract.len() as f32
}

/// Two records are "close" (bundle candidates) if any of their start/end
/// coordinates are within `thresh` bp (GangSTR `is_close`, `check_motif=False`;
/// chrom equality is implicit — these are one contig's records).
fn is_close(a: &RepeatInterval, b: &RepeatInterval, thresh: u64) -> bool {
    a.start.abs_diff(b.start) < thresh
        || a.start.abs_diff(b.end) < thresh
        || b.start.abs_diff(a.end) < thresh
        || b.end.abs_diff(a.end) < thresh
}

/// Group [`Admitted::bundled`] back into its clusters — one `Vec` per bundle.
///
/// [`admit`] returns the members flattened (spec §5a's shape), but the walk emits
/// **one region per cluster**, so it needs the grouping back. Re-deriving it is
/// exact rather than approximate: `bundled` is the concatenation of the clusters
/// in coordinate order, and two *adjacent* clusters are by definition not close
/// (or they would be one cluster) — so re-running the same [`is_close`] chaining
/// reproduces exactly the grouping [`split_bundles`] found.
///
/// It lives here, beside the flank test, rather than in the walk: **the cluster
/// rule is admission's**, and a second copy of `is_close` is how the two would
/// drift.
///
/// `flank_bp` must be the same radius [`admit`] was given, or the grouping will
/// not match the split.
///
/// # A singleton is a bug, and the walk is what keeps it that way
///
/// Every cluster has ≥ 2 members by definition (spec §2.4), so a singleton means `bundled`
/// and the split that produced it disagree. The check earned its keep at D3, catching a
/// truncated member that could not re-chain.
///
/// E2 briefly made it conditional: a walk restricted to *part* of a contig saw `bundled`
/// cut at the scan edge, and a member whose partner lay outside arrived alone — not a
/// disagreement, just the split minus what was never scanned. That was a symptom of
/// scanning less than a contig, and it is gone: the scan set is now **whole contigs**
/// (spec §2.5, owner 2026-07-17), so the only edges are the contig's, beyond which there
/// is nothing to have missed. The rule is unconditional again, and the assert with it.
pub fn bundle_clusters(bundled: &[RepeatInterval], flank_bp: u64) -> Vec<Vec<RepeatInterval>> {
    let mut out: Vec<Vec<RepeatInterval>> = Vec::new();
    for &iv in bundled {
        match out.last_mut() {
            Some(cluster) if is_close(cluster.last().expect("non-empty"), &iv, flank_bp) => {
                cluster.push(iv);
            }
            _ => out.push(vec![iv]),
        }
    }
    debug_assert!(
        out.iter().all(|c| c.len() >= 2),
        "a bundle has >= 2 members by definition (spec §2.4); a singleton means \
         the grouping disagrees with the split that produced it"
    );
    out
}

/// Split every maximal run of consecutive close records off from the isolated
/// ones (GangSTR `remove_bundles.py`'s selection). `recs` must be start-sorted.
/// Returns `(isolated, bundled)`.
///
/// **The selection is production's, verbatim; only the disposal differs** (spec
/// §2.4). `remove_bundles.py` — and A1, transcribing it — *deletes* the cluster;
/// this hands it back. The partition is the same one production computes, member
/// for member, which is what keeps spec §8's oracle valid: a bundle member is not
/// a locus in either system, so the two locus sets stay comparable. What changes
/// is only that ng can now *route* the cluster (as an `SsrBundle` region) instead
/// of leaving its bases an unaccounted hole.
///
/// No source comment of GangSTR's says why bundles are dropped at all (spec §10);
/// keeping selection and disposal separable is what lets that question be asked
/// later, with the evidence in hand.
fn split_bundles(
    recs: Vec<RepeatInterval>,
    thresh: u64,
) -> (Vec<RepeatInterval>, Vec<RepeatInterval>) {
    let mut isolated = Vec::new();
    let mut bundled = Vec::new();
    let n = recs.len();
    let mut i = 0;
    while i < n {
        if i + 1 < n && is_close(&recs[i], &recs[i + 1], thresh) {
            // Advance through the whole close cluster; set all of it aside.
            let mut j = i;
            while j + 1 < n && is_close(&recs[j], &recs[j + 1], thresh) {
                j += 1;
            }
            bundled.extend_from_slice(&recs[i..=j]);
            i = j + 1;
        } else {
            isolated.push(recs[i]);
            i += 1;
        }
    }
    (isolated, bundled)
}

impl fmt::Display for Locus {
    /// `chrom:start-end motif` — a readable identity for test diffs.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}:{}-{} {}",
            self.chrom,
            self.start,
            self.end,
            String::from_utf8_lossy(self.motif.as_bytes())
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The catalog's settings with a 5 bp flank, so the fixtures can be small.
    fn params() -> SsrAdmissionParams {
        matched_params(0.8, 0, 5).0
    }

    /// `admit` over a whole contig — spec §5a's **degenerate case**:
    /// `bases_start = 1`, `contig_len = bases.len()`, so the window IS the contig
    /// and the two roles `contig_seq.len()` used to play coincide again.
    ///
    /// This is what makes the A1 differential a valid oracle for A3's windowed
    /// code: production only ever ran whole-contig, so it is the only
    /// configuration the two can be compared at — and every window-only path is
    /// proved by this staying green plus the direct edge tests.
    fn admit_whole_contig(
        recs: Vec<RepeatInterval>,
        chrom: &str,
        contig: &[u8],
        p: &SsrAdmissionParams,
    ) -> Vec<Locus> {
        admit(recs, chrom, contig, Position(1), Bp(contig.len() as u64), p).loci
    }

    /// A 0-based half-open interval, the way the scanner emits one.
    fn iv(start: u64, end: u64, period: u8, score: i32) -> RepeatInterval {
        RepeatInterval {
            start,
            end,
            period,
            score,
        }
    }

    // ---- the ported motif ---------------------------------------------------

    #[test]
    fn motif_keeps_bytes_verbatim_and_reports_its_period() {
        let m = Motif::new(b"CAG").expect("valid period-3 motif");
        assert_eq!(
            m.as_bytes(),
            b"CAG",
            "verbatim: not rotated to a canonical form"
        );
        assert_eq!(m.period(), 3);
    }

    #[test]
    fn motif_rejects_lengths_outside_the_ssr_range() {
        assert!(matches!(
            Motif::new(b"").expect_err("empty"),
            MotifError::BadLength { len: 0 }
        ));
        assert!(matches!(
            Motif::new(b"ACGTACG").expect_err("period 7 > MAX_MOTIF_LEN"),
            MotifError::BadLength { len: 7 }
        ));
        assert!(Motif::new(b"ACGTAC").is_ok(), "period 6 is the ceiling");
    }

    /// The invariant the derived `Eq`/`Hash` rest on: the unused tail is zeroed,
    /// so two equal motifs built from different-length buffers still compare
    /// equal and no stale byte can make them differ.
    #[test]
    fn motif_equality_ignores_the_unused_tail() {
        assert_eq!(Motif::new(b"AT").unwrap(), Motif::new(b"AT").unwrap());
        assert_ne!(Motif::new(b"AT").unwrap(), Motif::new(b"ATA").unwrap());
        // A short motif built after a long one must not inherit its tail.
        let long = Motif::new(b"ACGTAC").unwrap();
        let short = Motif::new(b"AC").unwrap();
        assert_ne!(long, short);
        assert_eq!(short.as_bytes(), b"AC");
    }

    /// The port must not drift from the type it was copied from. Cheap to check:
    /// both are in this crate.
    #[test]
    fn motif_matches_productions_on_the_same_bytes() {
        for bytes in [b"AT".as_ref(), b"CAG".as_ref(), b"ACGTAC".as_ref()] {
            let ours = Motif::new(bytes).expect("ng motif");
            let theirs = crate::ssr::types::Motif::new(bytes).expect("production motif");
            assert_eq!(ours.as_bytes(), theirs.as_bytes());
            assert_eq!(ours.period(), theirs.period());
        }
        // And they agree on what is invalid.
        assert!(Motif::new(b"ACGTACG").is_err());
        assert!(crate::ssr::types::Motif::new(b"ACGTACG").is_err());
    }

    // ---- the locus type ----------------------------------------------------

    #[test]
    fn locus_reports_its_span_period_and_purity() {
        // 1-based inclusive: [6, 21] is 16 bases, 8 copies of AT.
        let l =
            Locus::new("chr1".into(), 6, 21, Motif::new(b"AT").unwrap(), 1.0).expect("valid locus");
        assert_eq!(l.start(), 6);
        assert_eq!(l.end(), 21);
        assert_eq!(l.tract_len(), 16, "inclusive: end - start + 1");
        assert_eq!(l.period(), 2);
        assert!(l.is_perfect());
    }

    #[test]
    fn locus_rejects_position_zero_because_it_is_one_based() {
        let err = Locus::new("chr1".into(), 0, 5, Motif::new(b"AT").unwrap(), 1.0)
            .expect_err("0 is not a 1-based position");
        assert!(matches!(
            err,
            LocusError::BadCoordinates { start: 0, end: 5 }
        ));
    }

    #[test]
    fn locus_rejects_an_inverted_span() {
        let err = Locus::new("chr1".into(), 9, 4, Motif::new(b"AT").unwrap(), 1.0)
            .expect_err("end precedes start");
        assert!(matches!(
            err,
            LocusError::BadCoordinates { start: 9, end: 4 }
        ));
    }

    #[test]
    fn locus_rejects_non_finite_purity() {
        let err = Locus::new("chr1".into(), 1, 6, Motif::new(b"AT").unwrap(), f32::NAN)
            .expect_err("NaN purity");
        assert!(matches!(err, LocusError::BadPurity { .. }));
    }

    // ---- the ported helpers (production's own cases) ------------------------

    #[test]
    fn is_compound_flags_atat_but_not_at_or_atc() {
        assert!(is_compound(b"ATAT"), "ATAT = (AT)^2 is compound");
        assert!(is_compound(b"ATATAT"), "ATATAT = (AT)^3 is compound");
        assert!(!is_compound(b"AT"), "fundamental AT is not compound");
        assert!(!is_compound(b"ATC"), "ATC is not internally periodic");
        assert!(!is_compound(b"ATCG"), "ATCG is not internally periodic");
    }

    #[test]
    fn purity_is_one_for_perfect_tiling() {
        assert_eq!(recompute_purity(b"ATATATAT", b"AT"), 1.0);
        assert_eq!(recompute_purity(b"CAGCAGCAG", b"CAG"), 1.0);
    }

    #[test]
    fn purity_drops_with_one_interruption() {
        assert_eq!(recompute_purity(b"ATATCTAT", b"AT"), 7.0 / 8.0);
    }

    #[test]
    fn minimal_trim_strips_partial_motifs_at_both_ends() {
        let clean = b"GATATATATATATC";
        let (st, en) = minimal_trim(clean, b"AT").expect("clean AT run trims");
        assert_eq!(&clean[st..en], b"ATATATATATAT");
    }

    #[test]
    fn minimal_trim_returns_none_without_two_clean_copies() {
        assert!(minimal_trim(b"ACGT", b"AT").is_none());
    }

    #[test]
    fn split_bundles_sets_aside_whole_clusters_and_keeps_isolated() {
        let recs = vec![
            (iv(100, 130, 2, 100)),
            (iv(150, 180, 3, 100)),
            (iv(200, 230, 2, 100)),
            (iv(5000, 5030, 2, 100)),
        ];
        let (isolated, bundled) = split_bundles(recs, 50);
        assert_eq!(isolated.len(), 1, "only the isolated D survives");
        assert_eq!(isolated[0].start, 5000);
        // A3: the cluster comes back rather than vanishing (spec §2.4).
        assert_eq!(bundled.len(), 3, "A, B, C are handed back, not deleted");
        assert_eq!(
            bundled.iter().map(|c| c.start).collect::<Vec<_>>(),
            vec![100, 150, 200]
        );
    }

    #[test]
    fn split_bundles_keeps_all_when_none_close() {
        let recs = vec![
            (iv(100, 130, 2, 100)),
            (iv(1000, 1030, 2, 100)),
            (iv(2000, 2030, 2, 100)),
        ];
        let (isolated, bundled) = split_bundles(recs, 50);
        assert_eq!(isolated.len(), 3);
        assert!(bundled.is_empty(), "nothing close, nothing set aside");
    }

    // ---- admit end-to-end, in ng's 1-based coordinates ----------------------

    /// The rebase's headline case: production asserts `start() == 5` on these
    /// exact bytes; ng must say **6**, and the tract bytes must be identical.
    #[test]
    fn admit_emits_a_clean_locus_one_based() {
        let left = b"CGCGC"; // 5 bp left flank
        let tract = b"ATATATATATATATAT"; // 16 bp = 8 copies of AT
        let right = b"GCGCG"; // 5 bp right flank
        let contig = [left.as_ref(), tract.as_ref(), right.as_ref()].concat();
        let loci = admit_whole_contig(vec![iv(5, 21, 2, 100)], "chr1", &contig, &params());

        assert_eq!(loci.len(), 1);
        let l = &loci[0];
        assert_eq!(l.chrom(), "chr1");
        assert_eq!(l.start(), 6, "1-based: production's 0-based 5, plus one");
        assert_eq!(l.end(), 21, "1-based inclusive == 0-based exclusive");
        assert_eq!(l.tract_len(), 16);
        assert_eq!(l.motif().as_bytes(), b"AT");
        assert_eq!(l.purity_fraction(), 1.0);
        // The flanks are why this locus was admitted, and they are not in it: a
        // locus says where a repeat is and what it is. The bases are in the
        // reference, which the caller has open (see `Locus`).
    }

    #[test]
    fn admit_drops_period_over_six_and_compound() {
        let contig = b"AAAACGATATATATATATATCCCC";
        let too_long = iv(0, 14, 7, 100); // period 7 > MAX_PERIOD
        let compound = iv(6, 20, 4, 100); // ATAT = (AT)^2
        let loci = admit_whole_contig(vec![too_long, compound], "chr1", contig, &params());
        assert!(loci.is_empty(), "both records are dropped, got {loci:?}");
    }

    #[test]
    fn admit_drops_below_copy_number_floor() {
        // period 2 needs >= 5 copies; (AT)*3 is only 3 → dropped.
        let contig = b"CCCCCATATATCCCCC";
        assert!(admit_whole_contig(vec![iv(5, 11, 2, 100)], "chr1", contig, &params()).is_empty());
    }

    #[test]
    fn admit_keeps_a_locus_whose_flank_clamps_to_one_base() {
        // Tract (AT)*8 starts at 0-based 1, so only 1 bp of the 50 bp left flank
        // exists. A 1 bp flank is still a (weak) anchor, so the locus is kept — the
        // gate is "clamped to NOTHING", not "clamped".
        let contig = b"GATATATATATATATATCGCGCGCGCG";
        let loci = admit_whole_contig(vec![iv(1, 17, 2, 100)], "chr1", contig, &params());
        assert_eq!(loci.len(), 1);
        assert_eq!(
            loci[0].start(),
            2,
            "the tract begins at the contig's 2nd base"
        );
    }

    /// **Every rejection reason, charged to the gate that fires it** (E1e, spec §3.1).
    ///
    /// One total cannot separate a purity rejection from a copy-floor one, and that is the
    /// distinction spec §10's routing question turns on — so each column has to be right,
    /// not just their sum. Each case below is built to trip **one** gate and to reach it:
    /// the gates run in a fixed order, so a fixture that trips two tells you nothing about
    /// the second.
    ///
    /// Driven through `admit` directly, without `prefilter`, and that is not a shortcut —
    /// it is the only place three of these are reachable at all (see
    /// `the_pre_filter_makes_three_of_admissions_gates_unreachable_from_the_walk`).
    #[test]
    fn each_rejection_reason_names_the_gate_that_fired() {
        let reasons = |recs: Vec<RepeatInterval>, contig: &[u8]| -> Vec<RejectionReason> {
            admit(
                recs,
                "chr1",
                contig,
                Position(1),
                Bp(contig.len() as u64),
                &params(),
            )
            .rejected
            .into_iter()
            .map(|(_, reason)| reason)
            .collect()
        };

        // CopyFloor: (AT)*3, under the period-2 floor of 5. Flanks either side, so the
        // flank gate cannot be what fires.
        let contig = b"CGCGCGCGCGCGCGCGCGCGATATATCGCGCGCGCGCGCGCGCGCG";
        assert_eq!(
            reasons(vec![iv(20, 26, 2, 100)], contig),
            vec![RejectionReason::CopyFloor],
            "three copies is under the floor of five"
        );

        // Compound: the motif is itself a repeat (`ATAT` is `AT` twice), so period 4 is a
        // lie about this tract. Rejected before any of the later gates.
        let contig = b"CGCGCGCGCGCGCGCGCGCGATATATATATATATATATATCGCGCGCGCGCGCGCGCGCG";
        assert_eq!(
            reasons(vec![iv(20, 40, 4, 100)], contig),
            vec![RejectionReason::Compound],
            "ATAT is AT twice: the period is a lie"
        );

        // NoCleanTrim: the prefix motif is `AT`, but there is no `ATAT` within the trim
        // window, so nothing here tiles and there are no whole-motif boundaries to trim
        // to. This gate runs BEFORE the copy floor, which is why it is what fires.
        let contig = b"CGCGCGCGCGCGCGCGCGCGATGGGGGGGGGATATATATCGCGCGCGCGCGCGCGCGCG";
        assert_eq!(
            reasons(vec![iv(20, 39, 2, 100)], contig),
            vec![RejectionReason::NoCleanTrim],
            "no two consecutive motifs within the trim window"
        );

        // Purity: clean ends (so it trims) and enough copies (so the floor passes), but
        // interrupted enough that a re-measure against a perfect tiling falls under 0.8.
        let contig = b"CGCGCGCGCGCGCGCGCGCGATATATATGGGGGGGGGGATATATATCGCGCGCGCGCGCGCGCGCG";
        assert_eq!(
            reasons(vec![iv(20, 46, 2, 100)], contig),
            vec![RejectionReason::Purity],
            "trims and clears the floor, but is not pure enough"
        );

        // FlankClamped: a tract abutting base 1 — a fact about the CONTIG's end (spec
        // §2.6), not about the tract.
        let contig = b"ATATATATATATATATCGCGC";
        assert_eq!(
            reasons(vec![iv(0, 16, 2, 100)], contig),
            vec![RejectionReason::FlankClamped],
            "no left flank to anchor against"
        );
    }

    /// **The walk reaches only ONE of admission's five gates** — `FlankClamped`. The other
    /// four columns of `TypedRegionCounts::rejected_by_reason` are structurally zero, and a
    /// reader has to know it: "no impure tracts in this genome" is a wrong thing to read
    /// from a zero the *scanner* caused.
    ///
    /// # This count has been wrong twice; here is the archaeology
    ///
    /// E1e wrote that `Compound` was unreachable and `Purity` live. The walk-level fixtures
    /// (spec §8's homopolymer and impure-tract cases) showed the reverse, and the count
    /// became two. Then the 2026-07-17 pre-filter ordering fix removed `Compound`'s only
    /// customer — the homopolymer alias — and the count is one. The lesson worth keeping is
    /// that **every one of these zeroes is caused by a stage upstream of the gate**, so the
    /// answer moves whenever that stage does; it is not a property of admission at all.
    ///
    /// What is true now, and why:
    ///
    /// - **`CopyFloor` — unreachable.** `prefilter` applies the **same `MinCopies` table**
    ///   with the same arithmetic, so nothing under the floor survives to be turned down
    ///   by it. A2 folded production's two copy-floor tables into one; both *call sites*
    ///   remain, and this is the visible consequence. *(E1e had this right.)*
    /// - **`Purity` — unreachable, and this is the interesting one.** Not because of the
    ///   pre-filter, but because of **Ruzzo–Tompa**: the scanner emits *maximal-scoring*
    ///   segments, and a tract impure enough to fail the 0.80 floor always contains a purer
    ///   sub-segment that scores higher — so the scanner emits *that*, and admission is
    ///   handed the pure core. Measured: a 0.79-purity fixture comes back a **locus**.
    ///   E1e argued from the two thresholds (scanner ≈0.78 vs admission 0.80) that tracts
    ///   would land in the gap; they cannot, because the segmenter never hands the gap over.
    /// - **`NoCleanTrim` — unreachable**, for the same reason: a tract with no motif pair
    ///   in its trim window is a tract the scanner segments away first.
    /// - **`Compound` — unreachable, since 2026-07-17.** E1e's *reasoning* was right all
    ///   along — a compound motif is a period-multiple, and redundancy elimination drops it
    ///   (`ATAT` is eliminated by `AT`) — but it was **false of period 1**, which
    ///   `prefilter` used to drop by the period floor *before* it could eliminate anything.
    ///   That let a poly-A's period-2 alias reach `admit` as motif `AA`, and made this the
    ///   homopolymer gate. The ordering fix gives period 1 its eliminating pass, so the
    ///   alias never survives and E1e's reasoning holds for every period. The gate stays as
    ///   [`admit`]'s own guard — it fires for a caller that skips the pre-filter — but the
    ///   walk no longer reaches it.
    /// - **`FlankClamped` — reachable**: the contig's ends are nobody else's business.
    ///
    /// This test pins the *mechanism* behind each zero, so that a pre-filter or scanner
    /// change fails here — where the explanation is — rather than as a column that
    /// mysteriously grows numbers. It did exactly that on 2026-07-17.
    #[test]
    fn the_walk_reaches_only_one_of_admissions_five_gates() {
        let p = params();

        // CopyFloor: `prefilter` enforces the floor FIRST, with the same table.
        let low_copy = vec![iv(20, 26, 2, 100)];
        assert!(
            prefilter(&low_copy, &p).is_empty(),
            "the pre-filter enforces the copy floor first, so `admit`'s copy-floor gate \
             cannot fire from the walk"
        );

        // Compound, on an STR multiple: redundancy elimination drops it, so THIS compound
        // never reaches the gate. (Both must be present, as the scanner emits them.)
        let fundamental_and_multiple = vec![iv(20, 40, 2, 100), iv(20, 40, 4, 100)];
        assert_eq!(
            prefilter(&fundamental_and_multiple, &p),
            vec![iv(20, 40, 2, 100)],
            "a period-multiple of a surviving tract is eliminated"
        );

        // **And a homopolymer's multiple is eliminated the same way** — since the
        // 2026-07-17 ordering fix, which is what took `Compound`'s only customer away.
        // The period-1 interval is kept as an eliminator, does its work, and is dropped
        // by the period floor afterwards.
        let homopolymer_and_multiple = vec![iv(20, 40, 1, 100), iv(20, 40, 2, 100)];
        assert!(
            prefilter(&homopolymer_and_multiple, &p).is_empty(),
            "the homopolymer eliminates its own period-2 alias, then the period floor \
             drops it: nothing survives, which is what `periods.min() == 2` asks for"
        );

        // Compound is therefore unreachable from the walk: every compound motif is a
        // period-multiple of a shorter one, and redundancy elimination now reaches all of
        // them — period 1 included. It stays as `admit`'s own guard, for a caller handing
        // in intervals the pre-filter never saw.
        let contig = b"CGCGCGCGCGCGCGCGCGCGAAAAAAAAAAAAAAAAAAAACGCGCGCGCGCGCGCGCGCG";
        assert_eq!(
            admit(
                vec![iv(20, 40, 2, 100)],
                "chr1",
                contig,
                Position(1),
                Bp(contig.len() as u64),
                &p,
            )
            .rejected
            .into_iter()
            .map(|(_, r)| r)
            .collect::<Vec<_>>(),
            vec![RejectionReason::Compound],
            "the gate still fires when fed a compound motif directly — `AA` is itself a \
             repeat — but the pre-filter no longer lets one reach it from the walk"
        );
    }

    #[test]
    fn admit_drops_locus_with_empty_left_flank() {
        // Tract abuts 0-based position 0 — no left flank, nothing to anchor.
        let contig = b"ATATATATATATATATCGCGC";
        assert!(
            admit_whole_contig(vec![iv(0, 16, 2, 100)], "chr1", contig, &params()).is_empty(),
            "a tract at contig position 0 has no left flank and is not genotypeable"
        );
    }

    #[test]
    fn admit_drops_locus_with_empty_right_flank() {
        let contig = b"CGCGCATATATATATATATAT";
        assert!(
            admit_whole_contig(vec![iv(5, 21, 2, 100)], "chr1", contig, &params()).is_empty(),
            "a tract ending on the contig's last base has no right flank"
        );
    }

    #[test]
    fn admit_drops_period_one_homopolymer_and_spares_the_neighbour_ssr() {
        let mut contig = vec![b'A'; 20];
        for _ in 0..10 {
            contig.extend_from_slice(b"CAG");
        }
        contig.extend_from_slice(b"TTTTT");
        let loci = admit_whole_contig(
            vec![iv(0, 20, 1, 100), iv(20, 50, 3, 100)],
            "chr1",
            &contig,
            &params(),
        );
        assert_eq!(loci.len(), 1, "the homopolymer is gone, the CAG survives");
        assert_eq!(loci[0].motif().as_bytes(), b"CAG");
        assert_eq!(loci[0].period(), 3);
    }

    // ---- the pre-filter -----------------------------------------------------

    #[test]
    fn prefilter_drops_low_copy_noise_and_period_multiples() {
        // A real (AT)*10 tract, re-detected at its period-4 multiple, plus a
        // 2-copy noise interval. Only the fundamental period-2 tract survives.
        let real = iv(100, 120, 2, 100);
        let multiple = iv(100, 120, 4, 100);
        let noise = iv(500, 504, 2, 10);
        let kept = prefilter(&[real, multiple, noise], &params());
        assert_eq!(
            kept,
            vec![real],
            "fundamental kept; multiple + noise dropped"
        );
    }

    #[test]
    fn prefilter_keeps_a_higher_period_tract_that_is_not_a_multiple() {
        // Period 3 overlapping period 2 is not a period-multiple → both stay.
        let two = iv(100, 120, 2, 100);
        let three = iv(100, 121, 3, 100);
        let kept = prefilter(&[two, three], &params());
        assert_eq!(kept.len(), 2);
    }

    // ---- spec §8.0: the port-fidelity differential --------------------------
    //
    // The oracle the Revision made necessary. ng's `admit` and production's
    // `build_loci` are two copies of one policy; this drives both from the same
    // intervals and asserts they still agree. It is the tripwire on silent
    // drift, and it is expected to be re-pinned the day an experiment
    // deliberately moves ng's admission away from the catalog's rules.

    use crate::ssr::catalog::CatalogParams;
    use crate::ssr::catalog::postprocess::build_loci;
    use crate::ssr::catalog::trf::TrfRecord;
    use crate::ssr::types::Locus as ProdLocus;

    /// One settings pair, built once, for both implementations.
    ///
    /// The differential is only meaningful if the two sides run **the same
    /// settings**, and two hand-maintained literal tables is exactly how they
    /// would quietly stop doing so — the test would keep passing while comparing
    /// two different questions. Deriving production's from ng's makes that
    /// impossible by construction.
    ///
    /// **`flank_bp` maps onto *both* of production's knobs**, which is the A2
    /// collapse stated as code: ng has one number where production has
    /// `flank_bp` and `bundle_threshold`, and this pair is only equivalent
    /// because production ships them equal (spec §2.4;
    /// `default_matches_the_frozen_catalog_params` pins it). Anywhere the two
    /// sides are compared, they must therefore agree — so the differential is
    /// also the evidence the collapse costs nothing.
    fn matched_params(
        min_purity: f32,
        min_score: i32,
        flank_bp: u32,
    ) -> (SsrAdmissionParams, CatalogParams) {
        (
            SsrAdmissionParams {
                min_purity,
                min_score,
                flank_bp: u64::from(flank_bp),
                // `periods` and `min_copies` take the catalog's values —
                // which is what makes them comparable with production's hardcoded
                // ones at all. A2's knobs are exercised by the tests that set them
                // explicitly, never by a silent default drifting in here.
                ..SsrAdmissionParams::default()
            },
            CatalogParams {
                min_purity,
                min_score,
                flank_bp,
                bundle_threshold: flank_bp,
            },
        )
    }

    /// Bridge ng's scanner intervals into production's parse shape. The only
    /// route from a `RepeatInterval` to `build_loci`, and `#[cfg(test)]` — which
    /// is exactly where it belongs (spec §5c).
    fn as_trf(intervals: &[RepeatInterval]) -> Vec<TrfRecord> {
        intervals
            .iter()
            .map(|iv| {
                // ng's `RepeatInterval` is `u64` (spec §4 / B2); `TrfRecord` is
                // production's `u32` parse shape. The differential's whole job is to
                // compare across that seam, so it narrows here — and `expect`s rather
                // than casting, because a truncated coordinate would compare the wrong
                // tract and the test would pass while asserting nothing.
                TrfRecord::for_test(
                    u32::try_from(iv.start).expect("fixture coordinates fit u32"),
                    u32::try_from(iv.end).expect("fixture coordinates fit u32"),
                    u16::from(iv.period),
                    iv.score,
                    b"",
                )
            })
            .collect()
    }

    /// **The conversion, stated once.** ng's 1-based inclusive `[start, end]` is
    /// production's 0-based half-open `[start, end)` shifted by exactly one on
    /// `start`; `end`, `motif` and `purity` do not move. Stating it in one place is
    /// what makes this test pin the arithmetic rather than restate a bug in both
    /// directions.
    #[track_caller]
    fn assert_same_locus(ng: &Locus, prod: &ProdLocus) {
        assert_eq!(ng.chrom(), prod.chrom(), "chrom");
        assert_eq!(
            ng.start(),
            u64::from(prod.start()) + 1,
            "start: 1-based == 0-based + 1 ({ng})"
        );
        assert_eq!(
            ng.end(),
            u64::from(prod.end()),
            "end: inclusive == exclusive ({ng})"
        );
        // **`ref_bytes` and the flanks are not compared, because ng no longer has
        // them** (2026-07-17). Production embeds the tract's bases plus a flank each
        // side; ng deliberately does not — a locus says where a repeat is and what it
        // is, and the bases are in the reference (see `Locus`). This is a divergence,
        // not drift, so the differential does not chase it. What it still pins is
        // every decision production makes: which tracts admit, and at what span, motif,
        // period and purity — which is what "the port is faithful" has to mean once the
        // payload is gone.
        // By bytes: the two `Motif`s are now distinct types (see `Motif`'s docs
        // — production's is `pub(crate)`, ng's is ported). Comparing the bytes is
        // also what would catch the ported type drifting from production's.
        assert_eq!(
            ng.motif().as_bytes(),
            prod.motif().as_bytes(),
            "motif ({ng})"
        );
        assert_eq!(ng.period(), prod.period(), "period ({ng})");
        assert_eq!(
            ng.purity_fraction(),
            prod.purity_fraction(),
            "purity ({ng})"
        );
        // The length identity the rebase turns on: production's `end - start`
        // (half-open) is ng's `end - start + 1` (inclusive) — same bases.
        assert_eq!(
            ng.tract_len(),
            u64::from(prod.end() - prod.start()),
            "tract_len ({ng})"
        );
    }

    /// Drive both implementations from one interval set, at the catalog's
    /// settings, and compare.
    #[track_caller]
    fn assert_agrees(intervals: &[RepeatInterval], chrom: &str, seq: &[u8], case: &str) {
        let (ng_p, prod_p) = matched_params(0.8, 0, 5);
        assert_agrees_at(&ng_p, &prod_p, intervals, chrom, seq, case);
    }

    /// As [`assert_agrees`], at an explicit configuration.
    ///
    /// The differential is a fixed-config regression test (spec §8.0), but a gate
    /// that never *fires* is not compared at all — both sides skip it together and
    /// the test stays green. So the gates whose floor is a no-op at the catalog's
    /// settings (`min_score`) or that no default fixture reaches (`min_purity`) get
    /// driven here at a floor that actually bites.
    #[track_caller]
    fn assert_agrees_at(
        ng_p: &SsrAdmissionParams,
        prod_p: &CatalogParams,
        intervals: &[RepeatInterval],
        chrom: &str,
        seq: &[u8],
        case: &str,
    ) {
        let ours = admit_whole_contig(intervals.to_vec(), chrom, seq, ng_p);
        let theirs = build_loci(as_trf(intervals), chrom, seq, prod_p);

        assert_eq!(
            ours.len(),
            theirs.len(),
            "{case}: locus count differs — ng {:?} vs production {:?}",
            ours.iter().map(Locus::to_string).collect::<Vec<_>>(),
            theirs
                .iter()
                .map(|l| format!("{}:{}-{}", l.chrom(), l.start(), l.end()))
                .collect::<Vec<_>>()
        );
        for (ng, prod) in ours.iter().zip(theirs.iter()) {
            assert_same_locus(ng, prod);
        }
    }

    /// The crafted cases: production's own `build_loci` tests, replayed through
    /// both implementations. Each one exercises a different gate.
    #[test]
    fn differential_vs_production_build_loci_on_crafted_cases() {
        // A clean isolated tract with flanks both sides.
        let contig = [
            b"CGCGC".as_ref(),
            b"ATATATATATATATAT".as_ref(),
            b"GCGCG".as_ref(),
        ]
        .concat();
        assert_agrees(&[iv(5, 21, 2, 100)], "chr1", &contig, "clean locus");

        // Period > MAX_PERIOD, and a compound motif.
        assert_agrees(
            &[iv(0, 14, 7, 100), iv(6, 20, 4, 100)],
            "chr1",
            b"AAAACGATATATATATATATCCCC",
            "period-scope + compound drop",
        );

        // Below the copy-number floor.
        assert_agrees(
            &[iv(5, 11, 2, 100)],
            "chr1",
            b"CCCCCATATATCCCCC",
            "copy floor",
        );

        // Left flank clamped but non-empty — kept.
        assert_agrees(
            &[iv(1, 17, 2, 100)],
            "chr1",
            b"GATATATATATATATATCGCGCGCGCG",
            "clamped left flank",
        );

        // Empty left flank / empty right flank — dropped.
        assert_agrees(
            &[iv(0, 16, 2, 100)],
            "chr1",
            b"ATATATATATATATATCGCGC",
            "empty left flank",
        );
        assert_agrees(
            &[iv(5, 21, 2, 100)],
            "chr1",
            b"CGCGCATATATATATATATAT",
            "empty right flank",
        );

        // Bundle drop: three close tracts + one isolated.
        let mut bundled = b"TTTTT".to_vec();
        bundled.resize(100, b'T');
        bundled.extend_from_slice(b"ATATATATATATATATATAT"); // 100..120
        bundled.resize(150, b'C');
        bundled.extend_from_slice(b"CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"); // 150..180
        bundled.resize(5000, b'G');
        bundled.extend_from_slice(b"ATATATATATATATATATAT"); // 5000..5020
        bundled.resize(5100, b'T');
        // **At a 50 bp radius, not `params()`'s 5.** A2 collapsed `flank_bp` and
        // `bundle_threshold` into one number, and this fixture was built for A1's
        // `flank_bp: 5, bundle_threshold: 50` — its closest pair is 30 bp apart, so
        // at a radius of 5 nothing bundles and the case silently stops testing the
        // bundle drop. It kept passing because *both sides moved together*: the
        // exact "wrong together, still green" blindness the differential has by
        // construction. Pinned by `admit_bundles_at_the_flank_radius` below.
        let (ng_p, prod_p) = matched_params(0.8, 0, 50);
        assert_agrees_at(
            &ng_p,
            &prod_p,
            &[
                iv(100, 120, 2, 100),
                iv(150, 180, 3, 100),
                iv(5000, 5020, 2, 100),
            ],
            "chr1",
            &bundled,
            "bundle drop keeps only the isolated tract",
        );

        // Period-1 homopolymer beside a real STR.
        let mut poly = vec![b'A'; 20];
        for _ in 0..10 {
            poly.extend_from_slice(b"CAG");
        }
        poly.extend_from_slice(b"TTTTT");
        assert_agrees(
            &[iv(0, 20, 1, 100), iv(20, 50, 3, 100)],
            "chr1",
            &poly,
            "homopolymer dropped before bundling",
        );
    }

    /// **B1 — soft-masked (lower-case) reference sequence.**
    ///
    /// The gap this closes: every other fixture is an upper-case literal, and the
    /// synthetic reference's only lower-case bytes are its contig *names*. So the
    /// `upper()` calls in `admit`/`finish_locus` were identity functions in 100% of
    /// the suite — including both differentials — and deleting them would have gone
    /// unnoticed.
    ///
    /// This is not an exotic input. Real references **soft-mask repeats**, which
    /// is precisely the sequence this module exists to process, so lower case is
    /// the mainstream case. And the failure is silent: `Locus` compares by value,
    /// so a lower-case motif reaching Milestone D simply fails to match the
    /// catalog's — a wrong locus, not a crash.
    ///
    /// *Two `upper()` calls now, not three: the third folded `ref_bytes`, which
    /// `Locus` no longer carries. The two that remain are the motif and the tract
    /// purity is computed against, which is the whole of a locus's case-sensitive
    /// surface.*
    #[test]
    fn admit_upper_cases_a_soft_masked_tracts_motif() {
        let left = b"cgcgc";
        let tract = b"atatatatatatatat"; // 8 copies, soft-masked
        let right = b"gcgcg";
        let contig = [left.as_ref(), tract.as_ref(), right.as_ref()].concat();

        let loci = admit_whole_contig(vec![iv(5, 21, 2, 100)], "chr1", &contig, &params());
        assert_eq!(loci.len(), 1, "a soft-masked tract is still a locus");
        let l = &loci[0];
        assert_eq!(l.motif().as_bytes(), b"AT", "motif upper-cased");
        assert_eq!(l.purity_fraction(), 1.0, "purity is case-insensitive");

        // And the case-folding matches production's, bytes included.
        assert_agrees(&[iv(5, 21, 2, 100)], "chr1", &contig, "soft-masked tract");

        // Mixed case, the real soft-masking shape: unique flanks upper, repeat
        // lower. Same locus as the all-upper contig — that is what "case-stable
        // identity" means.
        let mixed = [
            b"CGCGC".as_ref(),
            b"atatatatatatatat".as_ref(),
            b"GCGCG".as_ref(),
        ]
        .concat();
        let mixed_loci = admit_whole_contig(vec![iv(5, 21, 2, 100)], "chr1", &mixed, &params());
        assert_eq!(mixed_loci, loci, "case must not change the locus");
        assert_agrees(&[iv(5, 21, 2, 100)], "chr1", &mixed, "mixed-case contig");
    }

    /// **M3 — a tract `minimal_trim` actually trims.**
    ///
    /// The gap this closes: every other crafted case yields `st = 0`, so the
    /// `+ st` term in `new_start` was undetectable — deleting it passed the whole
    /// suite. The test that calls itself "the rebase's headline case" is itself an
    /// `st = 0` case, so the rebase's most delicate term was pinned by nothing.
    ///
    /// `ATC` prefixes the tract, so `minimal_trim` must skip 3 bases to reach the
    /// first clean `ATAT` boundary: `st = 3`, and the emitted start is the
    /// *trimmed* one, not the detector's.
    #[test]
    fn admit_reports_the_trimmed_start_not_the_detected_one() {
        let left = b"CGCGC"; // 5 bp flank
        let junk = b"ATC"; // trimmed away: no clean motif boundary
        let tract = b"ATATATATATAT"; // 6 clean copies of AT
        let right = b"GCGCG";
        let contig = [left.as_ref(), junk.as_ref(), tract.as_ref(), right.as_ref()].concat();
        // The detector reports the whole junk+tract span, as detectors do.
        let detected = iv(5, 20, 2, 100);

        let loci = admit_whole_contig(vec![detected], "chr1", &contig, &params());
        assert_eq!(loci.len(), 1);
        let l = &loci[0];
        assert_eq!(
            l.start(),
            9,
            "trimmed start: 0-based 5 (detected) + 3 (st) + 1 (1-based), NOT 6"
        );
        assert_eq!(l.end(), 20, "end unmoved");
        assert_eq!(
            l.tract_len(),
            tract.len() as u64,
            "the span covers the clean tract, junk trimmed off"
        );
        assert_eq!(l.motif().as_bytes(), b"AT");

        // And production agrees, which is what pins `+ st` against the original.
        assert_agrees(&[detected], "chr1", &contig, "trimmed start (st > 0)");
    }

    /// **M2 — the purity floor, which no fixture reached.**
    ///
    /// Every other tract is a perfect tiling against a floor of 0.8, so the gate
    /// never rejected on either side. That also left production's **one documented
    /// divergence from GangSTR** — "imperfect single-motif loci are kept"
    /// (`postprocess.rs`) — pinned by nothing at all.
    ///
    /// Both sides of the floor, driven through the differential.
    #[test]
    fn admit_keeps_an_imperfect_tract_above_the_floor_and_drops_one_below() {
        // (AT)*8 with one substitution: 15/16 = 0.9375, above a 0.8 floor.
        let left = b"CGCGC";
        let impure = b"ATATATATATCTATAT"; // one C where a A belongs
        let right = b"GCGCG";
        let contig = [left.as_ref(), impure.as_ref(), right.as_ref()].concat();
        let detected = iv(5, 21, 2, 100);

        let (ng_p, prod_p) = matched_params(0.8, 0, 5);
        let loci = admit_whole_contig(vec![detected], "chr1", &contig, &ng_p);
        assert_eq!(loci.len(), 1, "imperfect but above the floor: KEPT");
        assert!(
            !loci[0].is_perfect() && loci[0].purity_fraction() > 0.8,
            "purity {} should be imperfect and above 0.8",
            loci[0].purity_fraction()
        );
        assert_agrees_at(
            &ng_p,
            &prod_p,
            &[detected],
            "chr1",
            &contig,
            "impure, above floor",
        );

        // Same tract, floor raised above its purity: now the gate fires. Driving
        // the differential here is the point — it is the only way the purity gate
        // is compared against production at all.
        let (ng_hi, prod_hi) = matched_params(0.99, 0, 5);
        assert!(
            admit_whole_contig(vec![detected], "chr1", &contig, &ng_hi).is_empty(),
            "purity below the floor: DROPPED"
        );
        assert_agrees_at(
            &ng_hi,
            &prod_hi,
            &[detected],
            "chr1",
            &contig,
            "impure, below floor",
        );
    }

    /// **M1 — the `min_score` gate, which the differential was blind to.**
    ///
    /// Every fixture ran `min_score: 0` against `score: 100`, so the gate never
    /// fired on either side and **deleting it left the whole suite green**. A2's
    /// acceptance criterion ("differential green at `Default`") is also
    /// `min_score: 0`, so a mis-transcription would have survived A1 *and* A2 and
    /// first bitten an experiment that set a real floor.
    #[test]
    fn admit_applies_the_score_gate_when_the_floor_actually_bites() {
        let left = b"CGCGC";
        let tract = b"ATATATATATATATAT";
        let right = b"GCGCG";
        let contig = [left.as_ref(), tract.as_ref(), right.as_ref()].concat();

        let (ng_p, prod_p) = matched_params(0.8, 50, 5);
        // Above the floor: admitted.
        assert_eq!(
            admit_whole_contig(vec![iv(5, 21, 2, 100)], "chr1", &contig, &ng_p).len(),
            1,
            "score 100 >= floor 50"
        );
        assert_agrees_at(
            &ng_p,
            &prod_p,
            &[iv(5, 21, 2, 100)],
            "chr1",
            &contig,
            "score above floor",
        );

        // Below the floor: dropped. This is the assertion that makes deleting the
        // gate fail.
        assert!(
            admit_whole_contig(vec![iv(5, 21, 2, 49)], "chr1", &contig, &ng_p).is_empty(),
            "score 49 < floor 50 must be dropped"
        );
        assert_agrees_at(
            &ng_p,
            &prod_p,
            &[iv(5, 21, 2, 49)],
            "chr1",
            &contig,
            "score below floor",
        );

        // Exactly at the floor: `>=`, kept. Pins the boundary direction.
        assert_eq!(
            admit_whole_contig(vec![iv(5, 21, 2, 50)], "chr1", &contig, &ng_p).len(),
            1,
            "the gate is `score >= min_score`, so equality is admitted"
        );
    }

    /// **M9 — at `Default`, the score gate never fires, and that is deliberate.**
    ///
    /// Announces the trap the field doc describes: ng has no trf-mod `-s 30`
    /// upstream, so `DEFAULT_MIN_SCORE = 0` gates nothing the scanner can produce.
    /// If someone ever "tidies" the default to a TRF-shaped number, this fails and
    /// points at the scale mismatch (`RepeatInterval::score` is a Ruzzo–Tompa
    /// segment total, not a TRF score).
    ///
    /// **Writing this test corrected the claim it was written to prove.** The doc
    /// said "no score gate at all"; the gate is `score >= 0`, so `i32::MIN` *is*
    /// rejected. The honest statement is narrower and is what is asserted here:
    /// Ruzzo–Tompa emits only positive-scoring segments (`score = r - l > 0`,
    /// `tandem_repeat.rs`), so the floor is unreachable **for scanner output**.
    #[test]
    fn admit_at_default_never_rejects_a_score_the_scanner_can_emit() {
        let mut contig = b"CGCGC".to_vec();
        contig.resize(60, b'C'); // room for the 50 bp default flank
        contig.extend_from_slice(b"ATATATATATATATAT");
        contig.resize(140, b'G');
        let at_default = |score| {
            admit_whole_contig(
                vec![iv(60, 76, 2, score)],
                "chr1",
                &contig,
                &SsrAdmissionParams::default(),
            )
        };

        // 1 is the least a Ruzzo–Tompa segment can score, so this is the whole
        // reachable input range: at Default, nothing the scanner emits is gated.
        assert_eq!(
            at_default(1).len(),
            1,
            "the lowest emittable RT score is admitted"
        );
        assert_eq!(at_default(100).len(), 1);

        // The literal boundary, for the record: `>=` admits 0, and a negative
        // score would be dropped. Neither is reachable from the scanner — which is
        // exactly why "no score gate at all" was the wrong way to say this.
        assert_eq!(at_default(0).len(), 1, "the gate is `>=`");
        assert!(
            at_default(-1).is_empty(),
            "a negative score IS gated — the default is a no-op for scanner output, not no gate"
        );
    }

    /// **M8 — `Default` must equal the frozen catalog's, field for field.**
    ///
    /// Spec §8.1's `.cat` parity oracle compares ng's walk against a catalog built
    /// at these values; the comparison is meaningless if they drift. Production is
    /// frozen, so drift can only come from an ng-side edit — which should be
    /// deliberate, and is now loud. Before this test, `::default()` had **no caller
    /// in the tree** and no field of it was pinned by anything.
    #[test]
    fn default_matches_the_frozen_catalog_params() {
        let ours = SsrAdmissionParams::default();
        let theirs = CatalogParams::default();
        assert_eq!(ours.min_purity, theirs.min_purity, "min_purity");
        assert_eq!(ours.min_score, theirs.min_score, "min_score");
        assert_eq!(ours.flank_bp, u64::from(theirs.flank_bp), "flank_bp");

        // **The A2 collapse, justified rather than asserted.** ng has one number
        // where production has two. That is only lossless because production ships
        // them equal — so this is the evidence for the collapse, and the tripwire
        // if the premise ever stops holding.
        assert_eq!(
            u64::from(theirs.bundle_threshold),
            ours.flank_bp,
            "spec §2.4 collapses bundle_threshold into flank_bp; that is safe ONLY \
             because the catalog ships them at the same value"
        );

        // The knobs A2 hoisted: `Default` must still be the catalog's hardcoded
        // rules, or spec §8's oracles stop comparing like with like.
        //
        // **Against literals, not against our own consts.** Asserting
        // `ours.periods.min() == DEFAULT_MIN_PERIOD` would compare the constant
        // with itself — `Default` is *defined* as `PeriodRange::new(
        // DEFAULT_MIN_PERIOD, DEFAULT_MAX_PERIOD)` — so it could not fail, and
        // the one A2 knob that can silently decalibrate both §8 oracles would be
        // free to drift with this test green. Production's `MIN_PERIOD` /
        // `MAX_PERIOD` are private consts ng cannot read, so they are restated
        // here as the oracle, exactly as
        // `the_folded_copy_floor_table_reproduces_both_of_productions` does for
        // the minimums. If production ever moves, this fails and says so.
        assert_eq!(
            ours.periods.min(),
            2,
            "postprocess.rs's MIN_PERIOD is 2 — period-1 homopolymers are dropped"
        );
        assert_eq!(
            ours.periods.max(),
            6,
            "postprocess.rs's MAX_PERIOD is 6 — the microsatellite ceiling"
        );
        assert_eq!(
            ours.min_copies,
            MinCopies::default(),
            "the folded minimum table; its VALUES are pinned against production by \
             the_folded_copy_floor_table_reproduces_both_of_productions"
        );
    }

    /// **A2 — the folded table reproduces *both* of production's, exactly.**
    ///
    /// A1 carried production's two copy-floor tables verbatim; A2 folds them into
    /// one [`MinCopies`]. This is the fold's proof obligation: the single
    /// table must give what **each** original gave, at every period each could
    /// reach. The two originals live in frozen production and are private
    /// (`postprocess::copy_number_floor`, and `copy_floor` inside
    /// `scanner_parity`'s copy), so they are restated here as the
    /// oracle — which is the point: if production ever changes, this fails and says
    /// so, rather than ng drifting quietly.
    ///
    /// The docs called these tables "disagreeing" (spec §5/§10, and A1's code).
    /// They do not: their one numeric difference is period 1, which neither gate
    /// can reach. That is why A2 is a *deletion*, not a reconciliation.
    #[test]
    fn the_folded_copy_floor_table_reproduces_both_of_productions() {
        // `postprocess::copy_number_floor` — admission's, verbatim.
        fn productions_admission_floor(period: usize) -> u64 {
            match period {
                1 => 10,
                2 => 5,
                3 => 4,
                4 => 3,
                5 => 3,
                6 => 3,
                _ => 3,
            }
        }
        // `scanner_parity::catalog_prefilter`'s nested `copy_floor` — the second
        // copy, verbatim.
        fn productions_prefilter_floor(period: u8) -> u32 {
            match period {
                2 => 5,
                3 => 4,
                _ => 3,
            }
        }

        let folded = MinCopies::default();
        let defaults = SsrAdmissionParams::default();

        for period in defaults.periods.min()..=defaults.periods.max() {
            assert_eq!(
                u64::from(folded.for_period(period)),
                productions_admission_floor(period as usize),
                "folded table must match production's ADMISSION floor at period {period}"
            );
            assert_eq!(
                folded.for_period(period),
                productions_prefilter_floor(period),
                "folded table must match production's PRE-FILTER floor at period {period}"
            );
        }

        // The `_ => 3` arm both originals carry, for a period wider than the table.
        assert_eq!(folded.for_period(7), productions_prefilter_floor(7));
        assert_eq!(
            u64::from(folded.for_period(7)),
            productions_admission_floor(7)
        );

        // Period 1 is the ONLY value the two originals differ on (10 vs 3), and it
        // is unreachable: the default period floor is 2. That unreachability is
        // exactly why "disagreeing" was the wrong word for them.
        assert_ne!(
            productions_admission_floor(1),
            u64::from(productions_prefilter_floor(1)),
            "period 1 is the one value the two originals differ on"
        );
        assert!(
            defaults.periods.min() > 1,
            "if period 1 became admissible by default, the two originals WOULD \
             disagree and the fold would need a decision, not a deletion"
        );
    }

    // ---- A3: windowing, and what admission sets aside ------------------------

    /// A contig with one clean, isolated (AT)*8 tract at a known offset, padded so
    /// a window can be cut around it with room for flanks either side.
    fn contig_with_one_tract_at(tract_offset: usize, total: usize) -> (Vec<u8>, RepeatInterval) {
        let mut contig = vec![b'C'; tract_offset];
        contig.extend_from_slice(b"ATATATATATATATAT"); // 16 bp, 8 copies
        contig.resize(total, b'G');
        let iv = iv(tract_offset as u64, (tract_offset + 16) as u64, 2, 100);
        (contig, iv)
    }

    /// **A3 — the spec §2.6 bug, tested rather than argued.**
    ///
    /// A locus 10 bp from the *window's* right edge, but 10 kb from the *contig's*.
    /// Production's rule clamps the flank at `contig_seq.len()`, so handed this
    /// window it believes the chromosome stops at the edge, clamps the right flank
    /// to nothing, and drops the locus — **and a different set of loci for every
    /// `window_bp` you pick, with no error**. Telling `admit` where the contig
    /// really ends is the entire fix.
    #[test]
    fn admit_keeps_a_locus_the_window_edge_would_have_eaten() {
        // Contig: 10 kb. Tract at 0-based 1000..1016. The window is bases[900..1026]
        // = 1-based [901, 1026], so the tract ends **10 bp from the window's right
        // edge** — far less than the 50 bp flank.
        //
        // That window used to be illegal: `finish_locus` read the tract ± flank_bp to
        // embed it, so it panicked and the caller had to fetch a wider slice. Since
        // the locus stopped carrying bases (2026-07-17) the flank is a question about
        // `contig_len`, which the window cannot contradict — so this is exactly the
        // case that must work, and the tight window is the point of the fixture.
        let (contig, tract) = contig_with_one_tract_at(1000, 10_000);
        let win_start_1 = 901u64;
        let win = &contig[900..1026];
        // Re-base the interval onto the window slice.
        let in_window = iv(100, 116, 2, 100); // 1000-900 .. 1016-900

        let admitted = admit(
            vec![in_window],
            "chr1",
            win,
            Position(win_start_1),
            Bp(contig.len() as u64),
            &SsrAdmissionParams::default(), // flank_bp = 50
        );
        assert_eq!(
            admitted.loci.len(),
            1,
            "the contig continues 9 kb past the window, so the flank is real and the \
             locus must survive — even though the window holds only 10 bp of it"
        );
        let l = &admitted.loci[0];
        // Genomic 1-based, NOT window-relative: the window must not appear in the
        // output (spec §2.3).
        assert_eq!(l.start(), 1001, "genomic, not window-relative");
        assert_eq!(l.end(), 1016);

        // And it is the SAME locus the whole-contig call produces — the window is a
        // memory knob, never an answer (spec §2.3).
        let whole =
            admit_whole_contig(vec![tract], "chr1", &contig, &SsrAdmissionParams::default());
        assert_eq!(whole.len(), 1);
        assert_eq!(
            &admitted.loci, &whole,
            "windowed output must equal whole-contig output, byte for byte"
        );
    }

    /// The other half of §2.6: a tract genuinely at the contig's end still loses
    /// its flank and is still dropped. `contig_len` must be believed when it says
    /// "this really is the end", not only when it says "keep going".
    #[test]
    fn admit_still_drops_a_locus_at_the_real_contig_end() {
        // Tract runs to the last base of a 1016 bp contig, viewed through a window
        // that starts at 901.
        let (contig, _) = contig_with_one_tract_at(1000, 1016);
        let win = &contig[900..1016];
        let admitted = admit(
            vec![iv(100, 116, 2, 100)],
            "chr1",
            win,
            Position(901),
            Bp(contig.len() as u64),
            &SsrAdmissionParams::default(),
        );
        assert!(
            admitted.loci.is_empty(),
            "no right flank exists at the true contig end: not genotypeable"
        );
    }

    /// **A tract flush against the window's own edge is still a locus** — the clamps
    /// are `.max(1)` and `.min(contig_len)`, the *contig's* ends, and the window has
    /// no vote.
    ///
    /// **Mutation is what says these two tests are needed, and the fixture has to be
    /// flush.** The live mutants are `.max(1)` → `.max(bases_start)` and
    /// `.min(contig_len)` → `.min(bases_end)`: the window's ends instead of the
    /// contig's. Whenever a tract sits *anywhere inside* its window, both forms agree
    /// (`tract_start - flank_bp` either way), so no ordinary fixture can tell them
    /// apart. They diverge only when the clamp actually bites — a tract starting
    /// exactly at `bases_start`, or ending exactly at `bases_end` — where the mutant
    /// makes `ref_start == tract_start` and reports `FlankClamped`: **the locus
    /// silently vanishes because of where a memory knob fell**, which is spec §2.6's
    /// bug exactly.
    ///
    /// These two replace a pair of `should_panic` tests. Until 2026-07-17 admission
    /// *read* the flank to embed it, so this window was a caller bug that panicked;
    /// now the flank is arithmetic, the window is legal, and the assertion is about
    /// the answer rather than the crash.
    #[test]
    fn a_tract_at_the_windows_left_edge_keeps_its_flank_from_the_contig() {
        let (contig, _) = contig_with_one_tract_at(1000, 10_000);
        // The window STARTS at the tract: bases[1000..1076] = 1-based [1001, 1076].
        // Not one base of the 50 bp left flank is inside it — and the contig has
        // 1000 bp of it, so the locus stands.
        let win = &contig[1000..1076];
        let admitted = admit(
            vec![iv(0, 16, 2, 100)],
            "chr1",
            win,
            Position(1001),
            Bp(contig.len() as u64),
            &SsrAdmissionParams::default(),
        );
        assert_eq!(
            admitted.loci.len(),
            1,
            "the left flank is 1000 bp of contig the window cannot see: \
             {:?}",
            admitted.rejected
        );
        assert_eq!(admitted.loci[0].start(), 1001);
    }

    #[test]
    fn a_tract_at_the_windows_right_edge_keeps_its_flank_from_the_contig() {
        let (contig, _) = contig_with_one_tract_at(1000, 10_000);
        // The window ENDS at the tract: bases[900..1016] = 1-based [901, 1016].
        let win = &contig[900..1016];
        let admitted = admit(
            vec![iv(100, 116, 2, 100)],
            "chr1",
            win,
            Position(901),
            Bp(contig.len() as u64),
            &SsrAdmissionParams::default(),
        );
        assert_eq!(
            admitted.loci.len(),
            1,
            "the right flank is 9 kb of contig the window cannot see: {:?}",
            admitted.rejected
        );
        assert_eq!(admitted.loci[0].end(), 1016);
    }

    /// **`bundled` is in the INPUT's coordinate space, not genomic.**
    ///
    /// Pinned by nothing until now: every other test that reads `bundled` runs at
    /// `bases_start = 1`, where window offsets and genomic coordinates coincide —
    /// so the convention could silently flip and no test would notice. Milestone D
    /// re-bases these to build each `SsrBundle`'s hull; if they were already
    /// genomic it would **double-shift every hull by `bases_start`**, silently.
    ///
    /// Read this test as the contract D must code against.
    #[test]
    fn bundled_coordinates_are_window_offsets_not_genomic() {
        // Two tracts 30 bp apart, seen through a window starting at genomic 901.
        let mut contig = vec![b'C'; 1000];
        contig.extend_from_slice(b"ATATATATATATATATATAT"); // 1000..1020
        contig.resize(1050, b'C');
        contig.extend_from_slice(b"ATATATATATATATATATAT"); // 1050..1070
        contig.resize(10_000, b'C');

        let win = &contig[900..1200]; // genomic [901, 1200]
        // Window-relative: 1000-900 = 100, 1050-900 = 150.
        let a = iv(100, 120, 2, 100);
        let b = iv(150, 170, 2, 100);

        let admitted = admit(
            vec![a, b],
            "chr1",
            win,
            Position(901),
            Bp(contig.len() as u64),
            &SsrAdmissionParams::default(),
        );
        assert!(
            admitted.loci.is_empty(),
            "30 bp apart: both are bundle members"
        );
        assert_eq!(
            admitted.bundled,
            vec![a, b],
            "handed back verbatim in the caller's own space — NOT genomic \
             (which would read 1001/1051 here)"
        );
        // Said the other way, so the intent survives a careless edit: the genomic
        // positions do NOT appear.
        assert!(
            !admitted
                .bundled
                .iter()
                .any(|iv| iv.start == 1000 || iv.start == 1050),
            "genomic coordinates must not leak into `bundled`"
        );
    }

    /// `bases_start` is 1-based; `0` is not a position (spec §4). Untested until
    /// review caught it.
    #[test]
    #[should_panic(expected = "bases_start is 1-based")]
    fn admit_rejects_a_zero_bases_start() {
        let contig = [
            b"CGCGC".as_ref(),
            b"ATATATATATATATAT".as_ref(),
            b"GCGCG".as_ref(),
        ]
        .concat();
        let _ = admit(
            vec![iv(5, 21, 2, 100)],
            "chr1",
            &contig,
            Position(0),
            Bp(contig.len() as u64),
            &params(),
        );
    }

    /// A `contig_len` that understates the contig is spec §2.6's bug wearing the
    /// caller's clothes: `ref_end` clamps early, the `ref_end == tract_end` test
    /// still passes, and a locus comes out with a **silently truncated flank**.
    /// Caught in release, not just debug — §10's sweeps run the walk in release.
    #[test]
    #[should_panic(expected = "contig_len must be the CONTIG's length")]
    fn admit_rejects_a_contig_len_that_understates_the_contig() {
        let (contig, _) = contig_with_one_tract_at(1000, 10_000);
        let win = &contig[900..1076]; // genomic [901, 1076]
        let _ = admit(
            vec![iv(100, 116, 2, 100)],
            "chr1",
            win,
            Position(901),
            Bp(1000), // a lie: the window itself already reaches 1076
            &SsrAdmissionParams::default(),
        );
    }

    /// `flank_bp = 0` switches the whole step off — every tract fails its own flank
    /// test — and spec §10 plans to sweep `flank_bp`. The M7 shape: a knob that
    /// appears to move and instead silently returns nothing.
    #[test]
    #[should_panic(expected = "flank_bp must be at least 1")]
    fn admit_rejects_a_zero_flank() {
        let contig = [
            b"CGCGC".as_ref(),
            b"ATATATATATATATAT".as_ref(),
            b"GCGCG".as_ref(),
        ]
        .concat();
        let no_flank = SsrAdmissionParams {
            flank_bp: 0,
            ..params()
        };
        let _ = admit_whole_contig(vec![iv(5, 21, 2, 100)], "chr1", &contig, &no_flank);
    }

    /// **A3 — the bundle members come back** (spec §2.4, §6a).
    ///
    /// Production deletes them, so their bases become a hole nobody accounts for.
    /// ng keeps the *selection* and rejects the *disposal*, which is what lets the
    /// walk route them as an `SsrBundle` and lets spec §10's "what are bundles
    /// worth?" question be asked at all — the answer was previously deleted
    /// uncounted.
    #[test]
    fn admit_hands_back_the_bundle_members_it_sets_aside() {
        // Two tracts 30 bp apart (inside a 50 bp flank radius) plus one isolated
        // far away.
        let mut contig = vec![b'C'; 100];
        contig.extend_from_slice(b"ATATATATATATATATATAT"); // 100..120
        contig.resize(150, b'C');
        contig.extend_from_slice(b"ATATATATATATATATATAT"); // 150..170
        contig.resize(5000, b'C');
        contig.extend_from_slice(b"ATATATATATATATATATAT"); // 5000..5020
        contig.resize(5100, b'C');
        let a = iv(100, 120, 2, 100);
        let b = iv(150, 170, 2, 100);
        let lone = iv(5000, 5020, 2, 100);

        let admitted = admit(
            vec![a, b, lone],
            "chr1",
            &contig,
            Position(1),
            Bp(contig.len() as u64),
            &SsrAdmissionParams::default(),
        );

        assert_eq!(admitted.loci.len(), 1, "only the isolated tract is a locus");
        assert_eq!(admitted.loci[0].start(), 5001);
        assert_eq!(
            admitted.bundled,
            vec![a, b],
            "the cluster is handed back, verbatim and in coordinate order — \
             production would have deleted it"
        );

        // The selection is production's: the loci sets still agree, which is what
        // keeps spec §8's oracle valid (a bundle member is not a locus in either
        // system).
        let (ng_p, prod_p) = matched_params(0.8, 0, 50);
        assert_agrees_at(
            &ng_p,
            &prod_p,
            &[a, b, lone],
            "chr1",
            &contig,
            "bundled + isolated",
        );
    }

    /// Nothing is set aside when nothing clusters — `bundled` is not a dumping
    /// ground for every rejection. A repeat admission turns down for any *other*
    /// reason (impure, low-copy, compound) is the walk's `Generic` territory, not
    /// a bundle (spec §2.2), and must not appear here.
    #[test]
    fn admit_bundles_only_what_the_flank_test_rejects() {
        let contig = [
            b"CGCGC".as_ref(),
            b"ATATATATATATATAT".as_ref(),
            b"GCGCG".as_ref(),
        ]
        .concat();

        // A clean isolated tract: a locus, nothing bundled.
        let admitted = admit(
            vec![iv(5, 21, 2, 100)],
            "chr1",
            &contig,
            Position(1),
            Bp(contig.len() as u64),
            &params(),
        );
        assert_eq!(admitted.loci.len(), 1);
        assert!(
            admitted.bundled.is_empty(),
            "nothing close, nothing bundled"
        );

        // A tract rejected by the copy floor: neither a locus NOR a bundle member.
        let short = b"CCCCCATATATCCCCC";
        let admitted = admit(
            vec![iv(5, 11, 2, 100)],
            "chr1",
            short,
            Position(1),
            Bp(short.len() as u64),
            &params(),
        );
        assert!(admitted.loci.is_empty(), "below the copy floor");
        assert!(
            admitted.bundled.is_empty(),
            "a copy-floor rejection is Generic territory, not a bundle (spec §2.2)"
        );
    }

    /// **A2 — the period scope is a knob, and moving it moves the admitted set.**
    ///
    /// Spec §5.2's routing frontier is *"for which (period, tract length) cells does
    /// the STR path beat the generic one?"*. Production cannot be asked: its period
    /// scope is a `const`. This is the test that says ng can.
    #[test]
    fn narrowing_the_period_scope_changes_what_is_admitted() {
        // A (CAG)*10 tract — period 3, admitted at the default 2..=6 scope.
        let mut contig = b"TTTTT".to_vec();
        for _ in 0..10 {
            contig.extend_from_slice(b"CAG");
        }
        contig.extend_from_slice(b"TTTTT");
        let cag = iv(5, 35, 3, 100);

        assert_eq!(
            admit_whole_contig(vec![cag], "chr1", &contig, &params()).len(),
            1,
            "period 3 is inside the default 2..=6 scope"
        );

        // Take period 3 off the STR path: same tract, no locus.
        let dinucs_only = SsrAdmissionParams {
            periods: PeriodRange::new(2, 2).unwrap(),
            ..params()
        };
        assert!(
            admit_whole_contig(vec![cag], "chr1", &contig, &dinucs_only).is_empty(),
            "a 2..=2 scope must exclude the period-3 tract — the knob has to bite"
        );

        // And the floor moves too: period 1 is off by default, on if asked.
        let mut poly = vec![b'T'; 5];
        poly.extend(std::iter::repeat_n(b'A', 30));
        poly.extend_from_slice(b"TTTTT");
        let homopolymer = iv(5, 35, 1, 100);
        assert!(
            admit_whole_contig(vec![homopolymer], "chr1", &poly, &params()).is_empty(),
            "period-1 homopolymers are off the STR path by default"
        );
        let with_homopolymers = SsrAdmissionParams {
            periods: PeriodRange::new(1, 6).unwrap(),
            ..params()
        };
        assert_eq!(
            admit_whole_contig(vec![homopolymer], "chr1", &poly, &with_homopolymers).len(),
            1,
            "a 1..=6 scope admits it — the question spec §10 wants asked"
        );
    }

    /// **A2 — `flank_bp` really is the bundle radius, and moving it moves bundling.**
    ///
    /// The collapse's live wire: `drop_bundles(kept, p.flank_bp)`. Deleting that
    /// argument's link to `flank_bp` — or the bundle drop entirely — must not pass.
    /// Two tracts 30 bp apart: bundled at a 50 bp radius, both admitted at 5.
    #[test]
    fn admit_bundles_at_the_flank_radius() {
        // Two (AT) tracts with a 30 bp gap, each with room for a 50 bp flank.
        let mut contig = vec![b'C'; 100];
        contig.extend_from_slice(b"ATATATATATATATATATAT"); // 100..120
        contig.resize(150, b'C');
        contig.extend_from_slice(b"ATATATATATATATATATAT"); // 150..170
        contig.resize(300, b'C');
        let pair = [iv(100, 120, 2, 100), iv(150, 170, 2, 100)];

        // Radius 50: the 30 bp gap is inside it → one cluster → both dropped.
        let wide = matched_params(0.8, 0, 50).0;
        assert!(
            admit_whole_contig(pair.to_vec(), "chr1", &contig, &wide).is_empty(),
            "30 bp apart is inside a 50 bp flank radius: both are bundle members"
        );

        // Radius 5: the gap clears it → neither is bundled → both admitted.
        let narrow = matched_params(0.8, 0, 5).0;
        assert_eq!(
            admit_whole_contig(pair.to_vec(), "chr1", &contig, &narrow).len(),
            2,
            "30 bp apart clears a 5 bp flank radius: both are loci"
        );

        // And production agrees at both radii — which is the collapse's real
        // claim: ng's one number does what production's two did.
        let (ng_w, prod_w) = matched_params(0.8, 0, 50);
        assert_agrees_at(
            &ng_w,
            &prod_w,
            &pair,
            "chr1",
            &contig,
            "bundled at radius 50",
        );
        let (ng_n, prod_n) = matched_params(0.8, 0, 5);
        assert_agrees_at(
            &ng_n,
            &prod_n,
            &pair,
            "chr1",
            &contig,
            "unbundled at radius 5",
        );
    }

    /// **A2 — the per-period minimum is read *per period*, not as one number.**
    ///
    /// Found by review: `raising_the_copy_floor_changes_what_is_admitted` uses
    /// `uniform(9)`, which by construction makes every entry identical — so it
    /// pins *that* the table is consulted, never *which entry*. A `for_period`
    /// that ignored its argument and always returned `by_period[1]` passed every
    /// test in this file. And **non-uniform is the only shape spec §5.2's frontier
    /// needs**: the experiment varies the minimum *by period*.
    #[test]
    fn the_minimum_is_looked_up_by_the_tracts_own_period() {
        // A period-2 tract (8 copies) and a period-3 tract (10 copies), each
        // isolated and flanked.
        let mut contig = vec![b'C'; 20];
        contig.extend_from_slice(b"ATATATATATATATAT"); // 20..36, 8 copies of AT
        contig.resize(200, b'C');
        for _ in 0..10 {
            contig.extend_from_slice(b"CAG"); // 200..230, 10 copies of CAG
        }
        contig.resize(400, b'C');
        let di = iv(20, 36, 2, 100);
        let tri = iv(200, 230, 3, 100);

        // A table that admits the period-3 tract and rejects the period-2 one —
        // so a lookup ignoring the period cannot produce this answer.
        //          period:  1  2   3  4  5  6
        let picky = SsrAdmissionParams {
            min_copies: MinCopies::new([10, 9, 4, 3, 3, 3], 3),
            ..params()
        };
        let loci = admit_whole_contig(vec![di, tri], "chr1", &contig, &picky);
        assert_eq!(
            loci.len(),
            1,
            "only the period-3 tract clears its own minimum"
        );
        assert_eq!(loci[0].motif().as_bytes(), b"CAG");
        assert_eq!(loci[0].period(), 3);

        // Mirror it: swap which period is picky, and the surviving locus swaps.
        //           period:  1  2  3   4  5  6
        let mirrored = SsrAdmissionParams {
            min_copies: MinCopies::new([10, 5, 11, 3, 3, 3], 3),
            ..params()
        };
        let loci = admit_whole_contig(vec![di, tri], "chr1", &contig, &mirrored);
        assert_eq!(loci.len(), 1, "now only the period-2 tract clears its own");
        assert_eq!(loci[0].motif().as_bytes(), b"AT");

        // The same discrimination must reach `prefilter` — the differential is
        // blind to it, so nothing else would notice.
        assert_eq!(
            prefilter(&[di, tri], &picky),
            vec![tri],
            "pre-filter, by period"
        );
        assert_eq!(
            prefilter(&[di, tri], &mirrored),
            vec![di],
            "pre-filter, mirrored"
        );
    }

    /// `for_wider_periods` must be distinguishable from the last named period —
    /// `Default` and `uniform` both set them equal, so every other test would pass
    /// with a `for_period` that ignored it.
    #[test]
    fn min_copies_uses_the_fallback_only_beyond_the_table() {
        let m = MinCopies::new([10, 5, 4, 3, 3, 3], 99);
        assert_eq!(m.for_period(1), 10);
        assert_eq!(m.for_period(6), 3, "the last named period");
        assert_eq!(m.for_period(7), 99, "the first period past the table");
        assert_eq!(m.for_period(u8::MAX), 99);
        // Period 0 is unreachable (PeriodRange rejects a zero floor) but the fn is
        // total, and it agrees with production's `_ => 3` arm there.
        assert_eq!(m.for_period(0), 99);

        let u = MinCopies::uniform(4);
        assert_eq!(u.for_period(1), 4);
        assert_eq!(u.for_period(6), 4);
        assert_eq!(
            u.for_period(7),
            4,
            "uniform means uniform, fallback included"
        );
    }

    /// **A2 — the period ceiling is enforced in release, not just in debug.**
    ///
    /// `PeriodRange::new` bounds neither end, so `PeriodRange::new(2, 7)` is legal
    /// and this is reachable. It is an `assert!` rather than a `debug_assert!`
    /// because the knob exists to be **swept by the §5.2 experiment, and sweeps
    /// run in release** — a debug-only guard would let the experiment record
    /// "period 7 admits nothing" when the code never tried. That is a wrong
    /// result, not a missing panic.
    #[test]
    #[should_panic(expected = "exceeds MAX_MOTIF_LEN")]
    fn admit_rejects_a_period_ceiling_no_motif_can_hold() {
        let contig = [
            b"CGCGC".as_ref(),
            b"ATATATATATATATAT".as_ref(),
            b"GCGCG".as_ref(),
        ]
        .concat();
        let too_wide = SsrAdmissionParams {
            periods: PeriodRange::new(2, (MAX_MOTIF_LEN + 1) as u8).unwrap(),
            ..params()
        };
        let _ = admit_whole_contig(vec![iv(5, 21, 2, 100)], "chr1", &contig, &too_wide);
    }

    /// The `min_purity` contract `admit` guards: a `NaN` floor would make every
    /// `purity < p.min_purity` false, silently passing everything.
    #[test]
    #[should_panic(expected = "min_purity")]
    fn admit_rejects_a_nan_purity_floor() {
        let contig = [
            b"CGCGC".as_ref(),
            b"ATATATATATATATAT".as_ref(),
            b"GCGCG".as_ref(),
        ]
        .concat();
        let nan = SsrAdmissionParams {
            min_purity: f32::NAN,
            ..params()
        };
        let _ = admit_whole_contig(vec![iv(5, 21, 2, 100)], "chr1", &contig, &nan);
    }

    /// **A2 — the copy floor is a knob, and it is the "length" axis.**
    ///
    /// A floor of *n* copies at period *p* is a minimum tract length of *n·p*
    /// bases, so this knob **is** spec §5.2's length dimension. Production
    /// hardcodes it in a `const fn`, in two places.
    #[test]
    fn raising_the_copy_floor_changes_what_is_admitted() {
        // (AT)*8 = 8 copies at period 2. The default floor for period 2 is 5.
        let contig = [
            b"CGCGC".as_ref(),
            b"ATATATATATATATAT".as_ref(),
            b"GCGCG".as_ref(),
        ]
        .concat();
        let tract = iv(5, 21, 2, 100);

        assert_eq!(
            admit_whole_contig(vec![tract], "chr1", &contig, &params()).len(),
            1,
            "8 copies clears the default floor of 5"
        );

        // Raise the floor above the tract: dropped.
        let strict = SsrAdmissionParams {
            min_copies: MinCopies::uniform(9),
            ..params()
        };
        assert!(
            admit_whole_contig(vec![tract], "chr1", &contig, &strict).is_empty(),
            "8 copies is below a floor of 9"
        );

        // And it moves the PRE-FILTER too, which is the point of folding the two
        // tables: one knob, whole policy. A1's two copies would have moved half.
        assert_eq!(
            prefilter(&[tract], &params()),
            vec![tract],
            "kept at floor 5"
        );
        assert!(
            prefilter(&[tract], &strict).is_empty(),
            "the same knob must move the pre-filter — otherwise a swept floor \
             silently applies to only half the policy"
        );
    }

    /// **Mi6 — `is_close`'s boundary is strict `<`, on every clause that can show it.**
    ///
    /// A `<=` slip otherwise passes the whole suite, differential included: no
    /// other fixture puts two tracts exactly `thresh` apart. Ported from GangSTR's
    /// `is_close`, where the strictness is the definition.
    ///
    /// **`is_close` is four `abs_diff` comparisons, and they mask each other**, so
    /// which fixture you pick decides which clauses you actually pin. Callers hold
    /// `a.start <= b.start` (`admit` sorts before `drop_bundles`), and that
    /// constrains what is reachable:
    ///
    /// - **Disjoint tracts** put only the *gap* clause (`b.start - a.end`) on the
    ///   boundary — the other three sit well outside it and stay `false` either
    ///   way. Verified by mutation: this fixture alone leaves 3 of the 4 mutants
    ///   alive.
    /// - **Overlapping equal-length tracts** put *three* clauses on the boundary
    ///   at once (`b.start - a.start`, `b.start - a.end`, `b.end - a.end` are all
    ///   `thresh`), so flipping any one of them to `<=` flips the verdict.
    /// - The **`a.start`/`b.end` clause is unobservable at its boundary**: it needs
    ///   `b.end == a.start + thresh`, which forces `b.start - a.start < thresh`, so
    ///   the first clause has already fired. An equivalent mutant — no input kills
    ///   it, and that is a property of the ported predicate, not a gap here.
    #[test]
    fn is_close_is_strict_at_the_threshold() {
        // Disjoint: only the gap clause is at the boundary.
        let a = iv(100, 130, 2, 100);
        let gap_exactly = iv(180, 210, 2, 100);
        assert!(
            !is_close(&a, &gap_exactly, 50),
            "a gap of exactly `thresh` is not close (strict <)"
        );
        assert!(
            is_close(&a, &(iv(179, 209, 2, 100)), 50),
            "a gap of thresh - 1 is close"
        );

        // Overlapping, equal length: start-start, gap, and end-end are ALL exactly
        // `thresh`, so this one fixture pins three clauses at once.
        let long = iv(100, 200, 2, 100);
        let shifted = iv(150, 250, 2, 100);
        assert!(
            !is_close(&long, &shifted, 50),
            "three clauses sit exactly at `thresh`; strict < keeps them all false"
        );
        assert!(
            is_close(&long, &(iv(149, 249, 2, 100)), 50),
            "shift one closer and all three fire"
        );
    }

    /// **A homopolymer does not come back as a period-2 repeat** — the ordering
    /// bug fixed 2026-07-17, and the reason the period floor runs last.
    ///
    /// `AAAA…` tiles under `AA`, `AAA`, `AAAAA`, so the scanner emits the *same
    /// span* at every period in scope. Only the period-1 interval divides them all,
    /// so only it can eliminate them. The floor used to run first, and the aliases
    /// outlived their eliminator: periods 2, 3 and 5 all survived (4 and 6 died as
    /// multiples of 2), so one homopolymer entered the cleaned set as **three
    /// repeats** and fed coverage, the satellite cap and the rejection counts as if
    /// it were three. `periods.min() == 2` dropped the period-1 *label*; the
    /// homopolymer sailed on under a wrong one.
    ///
    /// This is what makes the period range mean what it says, which is all the
    /// caller asked for: **out of scope must mean gone, not relabelled.**
    #[test]
    fn a_homopolymer_does_not_survive_as_a_period_two_repeat() {
        // What the scanner really emits for a 30 bp poly-A: every period tiles it.
        // (Verified against `find_tandem_repeats` — this is its output, not a guess.)
        let aliases: Vec<RepeatInterval> = (1..=6).map(|p| iv(100, 130, p, 100)).collect();

        assert!(
            prefilter(&aliases, &params()).is_empty(),
            "at periods 2..=6 a homopolymer contributes NOTHING: the period-1 \
             interval eliminates its own aliases, then the period floor drops it. \
             Before the fix this returned periods 2, 3 and 5 — three 'repeats' \
             where the sequence has one homopolymer"
        );

        // And the range means what it says in the other direction too: put period 1
        // in scope and the homopolymer is a period-1 tract under its own motif —
        // one interval, not six.
        let with_homopolymers = SsrAdmissionParams {
            periods: PeriodRange::new(1, 6).unwrap(),
            ..params()
        };
        assert_eq!(
            prefilter(&aliases, &with_homopolymers),
            vec![iv(100, 130, 1, 100)],
            "at periods 1..=6 the homopolymer is kept — once, at its true period"
        );
    }

    /// **The copy floor is what keeps the period-1 eliminator honest**, and it is
    /// load-bearing precisely because the period floor now runs last.
    ///
    /// Aperiodic sequence is full of 2 bp period-1 specks (`GG`, `TT` — 45 of them
    /// in 300 bp of the test filler). If those reached redundancy elimination they
    /// would eliminate every real STR they overlap, since period 1 divides
    /// everything: the poly-A cascade, which scored 0/16 on the first parity run.
    /// [`MinCopies`]'s period-1 entry of **10** is what stops them — an entry that
    /// was unreachable dead weight until this ordering.
    #[test]
    fn a_two_bp_speck_is_not_an_eliminator_but_a_real_homopolymer_is() {
        let speck = iv(100, 102, 1, 10); // 2 copies at period 1 — under the floor of 10
        let real_str = iv(100, 130, 3, 100); // 10 copies at period 3, overlapping it

        assert_eq!(
            prefilter(&[speck, real_str], &params()),
            vec![real_str],
            "the speck is dropped by period 1's copy floor (10) before it can \
             eliminate anything: it never becomes an eliminator"
        );

        // The floor is a knob, not a constant: drop it and the speck does eliminate
        // the real tract. That is the cascade, and it is what the 10 buys.
        let no_floor = SsrAdmissionParams {
            min_copies: MinCopies::uniform(1),
            ..params()
        };
        assert!(
            prefilter(&[speck, real_str], &no_floor).is_empty(),
            "with the floor at 1 the speck survives, eliminates the period-3 tract \
             as its multiple, and is then dropped by the period floor — the whole \
             locus gone. This is why the copy floor must run FIRST"
        );
    }

    /// **M6 — a malformed interval is dropped, not underflowed.**
    ///
    /// `RepeatInterval`'s fields are `pub` with no ordering invariant, and
    /// `prefilter` is now `pub` — unlike the private test helper it was ported
    /// from. `end < start` would wrap `end - start` to a huge value in release and
    /// sail straight through the copy floor.
    #[test]
    fn prefilter_drops_an_inverted_interval_instead_of_underflowing() {
        let inverted = iv(120, 100, 2, 100);
        let good = iv(100, 120, 2, 100);
        assert_eq!(
            prefilter(&[inverted, good], &params()),
            vec![good],
            "an end < start interval is dropped, not wrapped into a huge copy count"
        );
        assert!(
            prefilter(&[iv(100, 100, 2, 100)], &params()).is_empty(),
            "empty span"
        );
    }

    /// The real case: the scanner's own output over the synthetic STR-diversity
    /// reference the golden catalog was built from. This is the input shape
    /// production has **never seen** (it only ever ate trf-mod's), so it is where
    /// a transcription slip would actually hide.
    #[test]
    fn differential_vs_production_build_loci_on_scanner_output() {
        use crate::ng::tandem_repeat::{PeriodRange, ScanParams, find_tandem_repeats};
        use std::fs::File;
        use std::io::BufReader;
        use std::path::Path;

        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("data")
            .join("tandem_repeat")
            .join("synthetic_ref.fa");
        let file = File::open(&fixture).expect("synthetic reference fixture");
        let mut reader = noodles_fasta::io::Reader::new(BufReader::new(file));

        let mut contigs = 0;
        let mut total_loci = 0;
        for result in reader.records() {
            let rec = result.expect("fasta record");
            let name = String::from_utf8_lossy(rec.name()).into_owned();
            let seq = rec.sequence().as_ref();

            // The catalog scans period 1..=6; the pre-filter and admission both
            // drop period 1 (`scanner_parity`'s configuration).
            let intervals =
                find_tandem_repeats(seq, PeriodRange::new(1, 6).unwrap(), &ScanParams::default());
            let cleaned = prefilter(&intervals, &params());
            assert_agrees(&cleaned, &name, seq, &format!("scanner output on {name}"));

            total_loci += admit_whole_contig(cleaned.clone(), &name, seq, &params()).len();
            contigs += 1;
        }
        assert!(contigs > 0, "the fixture must have contigs");
        assert!(
            total_loci > 0,
            "the fixture must admit loci — a differential over two empty sets proves nothing"
        );
    }
}
