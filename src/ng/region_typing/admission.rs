//! ng's STR admission policy — which detected repeats become STR loci.
//!
//! **This is a port of [`crate::ssr::catalog::postprocess`], not a rewrite.**
//! The rule set (period scope, score gate, compound-motif drop, bundle drop,
//! minimal trim, copy floor, purity floor, flank embed, contig-edge drop) is a
//! working, tested implementation — itself a faithful port of GangSTR's
//! `minimal_trim.py` / `remove_bundles.py` — and re-deriving it would be daft
//! (spec §5). So the **logic is transcribed unchanged**; only the shape is ng's.
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
//! 6. **Embed `ref_bytes`** (tract + flank each side, clamped at contig ends) and
//!    drop any locus whose flank clamped to zero — a tract with no anchor.

use std::fmt;
use std::ops::Range;

use crate::ng::tandem_repeat::{PeriodRange, RepeatInterval};

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
    /// Coordinates are not ordered `1 <= ref_bytes_start <= start <= end`.
    /// The `1 <=` is ng's: these are 1-based, so `0` is not a position.
    #[error(
        "locus coordinates out of order: require 1 <= ref_bytes_start ({ref_bytes_start}) \
         <= start ({start}) <= end ({end})"
    )]
    BadCoordinates {
        ref_bytes_start: u64,
        start: u64,
        end: u64,
    },
    /// The tract `end` runs past the embedded reference span.
    #[error(
        "tract end ({end}) exceeds embedded reference span \
         [{ref_bytes_start}, {ref_bytes_end}]"
    )]
    TractBeyondRefBytes {
        end: u64,
        ref_bytes_start: u64,
        ref_bytes_end: u64,
    },
    /// `ref_bytes` is empty. A 1-based **inclusive** span cannot represent an
    /// empty range (unlike production's half-open one), and a locus with no
    /// embedded reference is meaningless regardless.
    #[error("locus has empty ref_bytes")]
    EmptyRefBytes,
    /// `purity_fraction` is not a finite value in `[0.0, 1.0]`.
    #[error("purity fraction {purity_fraction} is not finite in [0.0, 1.0]")]
    BadPurity { purity_fraction: f32 },
}

/// One STR locus — a single repeat, carrying its own local reference bases.
///
/// ng's port of [`crate::ssr::types::Locus`]. Same fields, same meaning, two
/// differences (spec §4): coordinates are **1-based inclusive** and **`u64`**.
/// Production's stays 0-based/`u32`; neither converts to the other outside
/// [`self`]'s differential test.
///
/// `ref_bytes` spans `[ref_bytes_start, ref_bytes_start + ref_bytes.len() - 1]`
/// — the tract `[start, end]` plus a flank margin each side (clamped at contig
/// ends), upper-cased. `chrom` is the contig **name**, matching the project's
/// name-based contig model.
///
/// The fields are private and the invariant
/// `1 <= ref_bytes_start <= start <= end <= ref_bytes_start + ref_bytes.len() - 1`
/// (plus `purity_fraction` finite ∈ `[0.0, 1.0]`, `ref_bytes` non-empty) is
/// enforced by [`Locus::new`], so the accessors are infallible by construction.
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
    /// Embedded reference bases: the tract plus a flank margin each side,
    /// upper-cased, clamped at contig ends.
    ref_bytes: Box<[u8]>,
    /// Genomic coordinate (1-based) of `ref_bytes[0]`.
    ref_bytes_start: u64,
}

impl Locus {
    /// Build a locus, validating its coordinate and purity invariants.
    ///
    /// Mirrors `ssr::types::Locus::new`, rephrased for 1-based inclusive spans.
    /// The checks are not ceremony: [`admit`] derives these coordinates by
    /// arithmetic on detector output, and an off-by-one here is a *wrong locus,
    /// not a crash*.
    pub fn new(
        chrom: Box<str>,
        start: u64,
        end: u64,
        motif: Motif,
        purity_fraction: f32,
        ref_bytes: Box<[u8]>,
        ref_bytes_start: u64,
    ) -> Result<Self, LocusError> {
        if ref_bytes.is_empty() {
            return Err(LocusError::EmptyRefBytes);
        }
        if !(1 <= ref_bytes_start && ref_bytes_start <= start && start <= end) {
            return Err(LocusError::BadCoordinates {
                ref_bytes_start,
                start,
                end,
            });
        }
        // Inclusive: the last embedded base is at `ref_bytes_start + len - 1`.
        let ref_bytes_end = ref_bytes_start + ref_bytes.len() as u64 - 1;
        if end > ref_bytes_end {
            return Err(LocusError::TractBeyondRefBytes {
                end,
                ref_bytes_start,
                ref_bytes_end,
            });
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
            ref_bytes,
            ref_bytes_start,
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

    /// Embedded reference bases (tract plus flank margins, upper-cased).
    #[inline]
    pub fn ref_bytes(&self) -> &[u8] {
        &self.ref_bytes
    }

    /// Genomic coordinate (1-based) of `ref_bytes()[0]`.
    #[inline]
    pub fn ref_bytes_start(&self) -> u64 {
        self.ref_bytes_start
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

    /// The tract's byte range within `ref_bytes`. Both bounds hold by
    /// construction ([`Locus::new`]).
    #[inline]
    fn tract_range(&self) -> Range<usize> {
        debug_assert!(
            self.ref_bytes_start <= self.start,
            "ref_bytes precedes tract start"
        );
        let offset = (self.start - self.ref_bytes_start) as usize;
        offset..offset + self.tract_len() as usize
    }

    /// The reference tract bytes — the REF allele's sequence.
    pub fn ref_tract(&self) -> &[u8] {
        &self.ref_bytes[self.tract_range()]
    }

    /// The left flank: embedded reference bases before the tract.
    pub fn left_flank(&self) -> &[u8] {
        &self.ref_bytes[..self.tract_range().start]
    }

    /// The right flank: embedded reference bases after the tract.
    pub fn right_flank(&self) -> &[u8] {
        &self.ref_bytes[self.tract_range().end..]
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
/// margin embedded in `ref_bytes`. Both landed on 50 independently, and
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
    /// Flank margin (bp) embedded each side of the tract in `ref_bytes` — **and,
    /// since A2, the bundle radius too** (spec §2.4; see the type's docs).
    /// Default: [`DEFAULT_FLANK_BP`].
    ///
    /// One number, two jobs, because they are the same requirement seen twice: a
    /// locus needs `flank_bp` of clean sequence each side to anchor reads, so a
    /// repeat with another repeat inside `flank_bp` cannot have one — which is
    /// what makes it a bundle member rather than a locus.
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
// The candidate — the four fields the policy reads
// ---------------------------------------------------------------------

/// One candidate tract inside the port: exactly the four fields `build_loci`
/// reads off a `TrfRecord`, widened to ng's `u64`.
///
/// It exists so the transcription stays line-comparable with production's while
/// the input type differs, and so the widening from [`RepeatInterval`]'s `u32`
/// happens **once, at the entry**, rather than scattering `u64::from` through
/// the policy. Milestone B2 widens `RepeatInterval` itself, at which point the
/// conversion becomes an identity — this struct still earns its keep as the
/// "what admission actually needs" statement.
///
/// **Coordinates here are 0-based half-open** — the detector's and the slice's
/// own space. The conversion to ng's 1-based inclusive happens exactly once, at
/// [`Locus`] construction in [`finish_locus`], which is the one place the
/// arithmetic does not cancel (spec §4).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Candidate {
    /// Tract start, inclusive (0-based).
    start: u64,
    /// Tract end, exclusive (0-based).
    end: u64,
    /// Repeat period (motif length, bp). `u8`, like [`RepeatInterval::period`]
    /// and [`PeriodRange`] — A1 chose `u16` to mirror `TrfRecord`'s field, but
    /// `Candidate` is not a `TrfRecord`, and the width had no source and no sink:
    /// it only bought a `u16::from` at each gate and an unchecked `as u8` at the
    /// copy-count lookup.
    period: u8,
    /// Detector segment score.
    score: i32,
}

impl From<RepeatInterval> for Candidate {
    fn from(iv: RepeatInterval) -> Self {
        Self {
            start: u64::from(iv.start),
            end: u64::from(iv.end),
            period: iv.period,
            score: iv.score,
        }
    }
}

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
/// (`ssr::catalog::scanner_parity`), which had no production home; here it sits
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
/// The **floor**, by contrast, is load-bearing here and is tested directly: it
/// runs before redundancy elimination, and **period 1 divides every period**, so a
/// homopolymer that survives it eliminates every real STR it overlaps (see
/// `prefilter_drops_period_one_before_it_can_eliminate_a_real_str`).
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
    let mut floored: Vec<RepeatInterval> = intervals
        .iter()
        .copied()
        .filter(|iv| {
            // `iv.end > iv.start` first: it guards the subtraction below (see the
            // fn docs), and it is what `admit` independently requires anyway.
            iv.period >= params.periods.min()
                && iv.end > iv.start
                && (iv.end - iv.start) / u32::from(iv.period)
                    >= params.min_copies.for_period(iv.period)
        })
        .collect();
    // Process low periods first so a fundamental tract is kept and its multiples
    // dropped.
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
    kept
}

// ---------------------------------------------------------------------
// Admission
// ---------------------------------------------------------------------

/// Admit detected repeats into start-sorted STR loci — ng's `build_loci`.
///
/// `chrom` is the contig name; `contig_seq` is the full contig (any case — the
/// tract, motif, and `ref_bytes` are upper-cased here for case-stable identity).
/// Feed it [`prefilter`]ed intervals, never raw scanner output (spec §5b).
///
/// **A1: whole-contig, as production is.** Milestone A3 adds `bases_start` +
/// `contig_len` so a window can be admitted without the slice being mistaken for
/// the whole chromosome (spec §2.6), and swaps the return for
/// `Admitted { loci, bundled }` (spec §2.4).
///
/// Output coordinates are **1-based inclusive** ([`Locus`]); the input's are
/// 0-based half-open. That conversion is the only intended numeric difference
/// from production's `build_loci`, and the two
/// `differential_vs_production_build_loci_*` tests pin it.
pub fn admit(
    recs: Vec<RepeatInterval>,
    chrom: &str,
    contig_seq: &[u8],
    p: &SsrAdmissionParams,
) -> Vec<Locus> {
    // The `SsrAdmissionParams` contracts that are otherwise only prose. A NaN
    // `min_purity` is the one that bites quietly: every `purity < p.min_purity`
    // comparison below is then false, so the purity gate passes everything rather
    // than failing.
    debug_assert!(
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
    let mut kept: Vec<Candidate> = recs
        .into_iter()
        .map(Candidate::from)
        .filter(|r| {
            r.period >= p.periods.min()
                && r.period <= p.periods.max()
                && r.score >= p.min_score
                && r.end > r.start
                && (r.end as usize) <= contig_seq.len()
        })
        .filter(|r| {
            let period = r.period as usize;
            let tract = &contig_seq[r.start as usize..r.end as usize];
            // Need at least one full motif to form the prefix.
            if tract.len() < period {
                return false;
            }
            let motif = upper(&tract[..period]);
            !is_compound(&motif)
        })
        .collect();

    // 3. drop bundles (on the raw, pre-trim coordinates). Records must be
    //    start-sorted for the streaming clustering.
    kept.sort_by_key(|r| (r.start, r.end));
    // `flank_bp` IS the bundle radius (spec §2.4) — see `SsrAdmissionParams`.
    let kept = drop_bundles(kept, p.flank_bp);

    // 4-5. per record: end-trim + copy floor, recompute purity + floor, embed.
    let mut out = Vec::with_capacity(kept.len());
    for r in &kept {
        if let Some(locus) = finish_locus(r, chrom, contig_seq, p) {
            out.push(locus);
        }
    }
    out
}

/// Steps 4-5 for one record: end-trim, copy-number floor, motif, purity floor,
/// `ref_bytes` embed. `None` if the record fails any gate.
///
/// **The one place the coordinate base changes.** Everything above works in the
/// detector's 0-based half-open space (which is also the slice's, so the
/// arithmetic is production's unchanged); the 1-based inclusive [`Locus`] is
/// built at the bottom, once. Spec §4 predicted exactly this: "the arithmetic
/// cancels … exactly one site does not".
fn finish_locus(
    r: &Candidate,
    chrom: &str,
    contig_seq: &[u8],
    p: &SsrAdmissionParams,
) -> Option<Locus> {
    let period = r.period as usize;
    let raw_tract = upper(&contig_seq[r.start as usize..r.end as usize]);
    let motif_bytes = raw_tract.get(..period)?.to_vec();

    // End-trim to clean whole-motif boundaries (GangSTR minimal_trim).
    let (st, en) = minimal_trim(&raw_tract, &motif_bytes)?;
    let new_start = r.start + st as u64;
    let new_end = r.start + en as u64;
    let trimmed = &raw_tract[st..en];

    // Copy-number floor — GangSTR computes copies from the ORIGINAL span
    // (integer division), as an accept-gate after trimming.
    let ref_copy = (r.end - r.start) / u64::from(r.period);
    if ref_copy < u64::from(p.min_copies.for_period(r.period)) {
        return None;
    }

    // After minimal_trim the tract starts on a motif boundary, so
    // `trimmed[..period] == motif_bytes`; use it as the phase-faithful motif.
    //
    // A1 wrote `.ok()?` here, faithfully to production, where it was unreachable:
    // the period ceiling was a `const` a static assert held at `<= MAX_MOTIF_LEN`.
    // **A2 made that ceiling a knob and the arm reachable** — set
    // `periods.max() > MAX_MOTIF_LEN` and every widest-period tract clears the
    // gate, fails here, and vanishes through the same `None` as a legitimate
    // policy rejection. That is M7's silent failure exactly: the knob appears to
    // move and does not. So it fails loudly instead (release behaviour unchanged —
    // still `None`; `admit` debug_asserts the same contract up front).
    let motif = match Motif::new(&motif_bytes) {
        Ok(m) => m,
        Err(e) => {
            debug_assert!(
                false,
                "period {} admitted but has no valid Motif ({e}) — is \
                 SsrAdmissionParams::periods.max() above MAX_MOTIF_LEN ({MAX_MOTIF_LEN})?",
                r.period
            );
            return None;
        }
    };

    // Recompute purity from the trimmed tract vs a perfect motif tiling.
    let purity = recompute_purity(trimmed, &motif_bytes);
    if purity < p.min_purity {
        return None;
    }

    // Embed ref_bytes: trimmed tract + flank each side, clamped at contig ends.
    let ref_start = new_start.saturating_sub(p.flank_bp);
    let ref_end = (new_end + p.flank_bp).min(contig_seq.len() as u64);

    // Drop a locus whose flank clamped to nothing on either side — a tract
    // abutting position 0 of the contig (empty left flank) or ending on the
    // contig's last base (empty right flank). The delimiter anchors the repeat
    // region on *both* flank junctions; a zero-length flank leaves nothing to
    // anchor against, so the tract is not genotypeable.
    if ref_start == new_start || ref_end == new_end {
        return None;
    }

    let ref_bytes = upper(&contig_seq[ref_start as usize..ref_end as usize]);

    // 0-based half-open -> 1-based inclusive (spec §4). `[s, e)` is `[s+1, e]`:
    // the start shifts by one, the end does not move, and the length is
    // unchanged — which is why every span above could stay in production's
    // arithmetic.
    match Locus::new(
        chrom.to_string().into_boxed_str(),
        new_start + 1,
        new_end,
        motif,
        purity,
        ref_bytes.into_boxed_slice(),
        ref_start + 1,
    ) {
        Ok(locus) => Some(locus),
        // Unreachable: every gate above guarantees the invariants (`minimal_trim`
        // gives `st < en`, the flank-clamp check gives a non-empty `ref_bytes`
        // strictly containing the tract, `recompute_purity` returns `[0, 1]`).
        //
        // Which is exactly why the verdict is worth keeping. Production writes
        // `.ok()` here and so discards it; that is faithful but wrong for a port,
        // because this error can now only fire if **the transcription is wrong** —
        // the one failure A1 exists to catch, and one that is otherwise silent (a
        // bad locus would leave through the same `None` as a routine policy
        // rejection). `debug_assert` keeps release behaviour byte-identical to
        // production while making the differential — which runs in debug — fail
        // loudly instead of quietly disagreeing.
        Err(e) => {
            debug_assert!(
                false,
                "admit built an invalid locus, so the port's arithmetic is wrong: {e} \
                 (tract {new_start}..{new_end} 0-based, ref_bytes {ref_start}..{ref_end})"
            );
            None
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
fn is_close(a: &Candidate, b: &Candidate, thresh: u64) -> bool {
    a.start.abs_diff(b.start) < thresh
        || a.start.abs_diff(b.end) < thresh
        || b.start.abs_diff(a.end) < thresh
        || b.end.abs_diff(a.end) < thresh
}

/// Drop bundles: discard every maximal run of consecutive close records in its
/// entirety (GangSTR `remove_bundles.py`). `recs` must be start-sorted.
///
/// **A3 keeps this selection and rejects its disposal** — the cluster members are
/// handed back as an `SsrBundle` rather than deleted (spec §2.4). A1 transcribes
/// the drop as-is, because the differential can only compare what production
/// also does.
fn drop_bundles(recs: Vec<Candidate>, thresh: u64) -> Vec<Candidate> {
    let mut out = Vec::new();
    let n = recs.len();
    let mut i = 0;
    while i < n {
        if i + 1 < n && is_close(&recs[i], &recs[i + 1], thresh) {
            // Advance through the whole close cluster; drop all of it.
            let mut j = i;
            while j + 1 < n && is_close(&recs[j], &recs[j + 1], thresh) {
                j += 1;
            }
            i = j + 1;
        } else {
            out.push(recs[i]);
            i += 1;
        }
    }
    out
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

    /// A 0-based half-open interval, the way the scanner emits one.
    fn iv(start: u32, end: u32, period: u8, score: i32) -> RepeatInterval {
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
    fn locus_accessors_split_ref_bytes_into_flanks_and_tract() {
        // Built by concatenation so the coordinates cannot drift from the bytes.
        // 1-based: ref_bytes covers [1, 26]; the tract is [6, 21] = 16 bases.
        let left = b"CGCGC";
        let tract = b"ATATATATATATATAT"; // 8 copies of AT
        let right = b"GCGCG";
        let ref_bytes = [left.as_ref(), tract.as_ref(), right.as_ref()].concat();
        let l = Locus::new(
            "chr1".into(),
            left.len() as u64 + 1,
            (left.len() + tract.len()) as u64,
            Motif::new(b"AT").unwrap(),
            1.0,
            ref_bytes.into_boxed_slice(),
            1,
        )
        .expect("valid locus");
        assert_eq!(l.start(), 6);
        assert_eq!(l.end(), 21);
        assert_eq!(l.left_flank(), left);
        assert_eq!(l.ref_tract(), tract);
        assert_eq!(l.right_flank(), right);
        assert_eq!(l.tract_len(), 16, "inclusive: end - start + 1");
        assert_eq!(l.period(), 2);
        assert!(l.is_perfect());
    }

    #[test]
    fn locus_rejects_position_zero_because_it_is_one_based() {
        let err = Locus::new(
            "chr1".into(),
            0,
            5,
            Motif::new(b"AT").unwrap(),
            1.0,
            (*b"ATATAT").into(),
            0,
        )
        .expect_err("0 is not a 1-based position");
        assert!(matches!(err, LocusError::BadCoordinates { .. }));
    }

    #[test]
    fn locus_rejects_empty_ref_bytes() {
        let err = Locus::new(
            "chr1".into(),
            1,
            1,
            Motif::new(b"AT").unwrap(),
            1.0,
            [].into(),
            1,
        )
        .expect_err("an inclusive span cannot be empty");
        assert!(matches!(err, LocusError::EmptyRefBytes));
    }

    #[test]
    fn locus_rejects_tract_past_ref_bytes() {
        // ref_bytes covers [1, 6]; a tract ending at 7 runs past it.
        let err = Locus::new(
            "chr1".into(),
            2,
            7,
            Motif::new(b"AT").unwrap(),
            1.0,
            (*b"ATATAT").into(),
            1,
        )
        .expect_err("tract must fit inside ref_bytes");
        assert!(matches!(err, LocusError::TractBeyondRefBytes { .. }));
    }

    #[test]
    fn locus_rejects_non_finite_purity() {
        let err = Locus::new(
            "chr1".into(),
            1,
            6,
            Motif::new(b"AT").unwrap(),
            f32::NAN,
            (*b"ATATAT").into(),
            1,
        )
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
    fn drop_bundles_discards_whole_clusters_keeps_isolated() {
        let recs = vec![
            Candidate::from(iv(100, 130, 2, 100)),
            Candidate::from(iv(150, 180, 3, 100)),
            Candidate::from(iv(200, 230, 2, 100)),
            Candidate::from(iv(5000, 5030, 2, 100)),
        ];
        let kept = drop_bundles(recs, 50);
        assert_eq!(kept.len(), 1, "only the isolated D survives");
        assert_eq!(kept[0].start, 5000);
    }

    #[test]
    fn drop_bundles_keeps_all_when_none_close() {
        let recs = vec![
            Candidate::from(iv(100, 130, 2, 100)),
            Candidate::from(iv(1000, 1030, 2, 100)),
            Candidate::from(iv(2000, 2030, 2, 100)),
        ];
        assert_eq!(drop_bundles(recs, 50).len(), 3);
    }

    // ---- admit end-to-end, in ng's 1-based coordinates ----------------------

    /// The rebase's headline case: production asserts `start() == 5` on these
    /// exact bytes; ng must say **6**, and the tract bytes must be identical.
    #[test]
    fn admit_emits_a_clean_locus_with_flanks_one_based() {
        let left = b"CGCGC"; // 5 bp left flank
        let tract = b"ATATATATATATATAT"; // 16 bp = 8 copies of AT
        let right = b"GCGCG"; // 5 bp right flank
        let contig = [left.as_ref(), tract.as_ref(), right.as_ref()].concat();
        let loci = admit(vec![iv(5, 21, 2, 100)], "chr1", &contig, &params());

        assert_eq!(loci.len(), 1);
        let l = &loci[0];
        assert_eq!(l.chrom(), "chr1");
        assert_eq!(l.start(), 6, "1-based: production's 0-based 5, plus one");
        assert_eq!(l.end(), 21, "1-based inclusive == 0-based exclusive");
        assert_eq!(l.tract_len(), 16);
        assert_eq!(l.motif().as_bytes(), b"AT");
        assert_eq!(l.purity_fraction(), 1.0);
        assert_eq!(l.ref_tract(), tract);
        assert_eq!(l.left_flank(), left);
        assert_eq!(l.right_flank(), right);
        assert_eq!(l.ref_bytes_start(), 1, "ref_bytes start at contig base 1");
    }

    #[test]
    fn admit_drops_period_over_six_and_compound() {
        let contig = b"AAAACGATATATATATATATCCCC";
        let too_long = iv(0, 14, 7, 100); // period 7 > MAX_PERIOD
        let compound = iv(6, 20, 4, 100); // ATAT = (AT)^2
        let loci = admit(vec![too_long, compound], "chr1", contig, &params());
        assert!(loci.is_empty(), "both records are dropped, got {loci:?}");
    }

    #[test]
    fn admit_drops_below_copy_number_floor() {
        // period 2 needs >= 5 copies; (AT)*3 is only 3 → dropped.
        let contig = b"CCCCCATATATCCCCC";
        assert!(admit(vec![iv(5, 11, 2, 100)], "chr1", contig, &params()).is_empty());
    }

    #[test]
    fn admit_clamps_flanks_at_contig_ends() {
        // Tract (AT)*8 starts at 0-based 1; the 5 bp left flank clamps to 1 bp.
        // A 1 bp flank is still a (weak) anchor, so the locus is kept.
        let contig = b"GATATATATATATATATCGCGCGCGCG";
        let loci = admit(vec![iv(1, 17, 2, 100)], "chr1", contig, &params());
        assert_eq!(loci.len(), 1);
        let l = &loci[0];
        assert_eq!(l.ref_bytes_start(), 1, "clamped to the contig's first base");
        assert_eq!(l.left_flank(), b"G", "only 1 bp of left flank available");
    }

    #[test]
    fn admit_drops_locus_with_empty_left_flank() {
        // Tract abuts 0-based position 0 — no left flank, nothing to anchor.
        let contig = b"ATATATATATATATATCGCGC";
        assert!(
            admit(vec![iv(0, 16, 2, 100)], "chr1", contig, &params()).is_empty(),
            "a tract at contig position 0 has no left flank and is not genotypeable"
        );
    }

    #[test]
    fn admit_drops_locus_with_empty_right_flank() {
        let contig = b"CGCGCATATATATATATATAT";
        assert!(
            admit(vec![iv(5, 21, 2, 100)], "chr1", contig, &params()).is_empty(),
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
        let loci = admit(
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
            .map(|iv| TrfRecord::for_test(iv.start, iv.end, u16::from(iv.period), iv.score, b""))
            .collect()
    }

    /// **The conversion, stated once.** ng's 1-based inclusive `[start, end]` is
    /// production's 0-based half-open `[start, end)` shifted by exactly one on
    /// `start` and `ref_bytes_start`; `end`, `ref_bytes`, `motif`, and `purity`
    /// do not move. Stating it in one place is what makes this test pin the
    /// arithmetic rather than restate a bug in both directions.
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
        assert_eq!(
            ng.ref_bytes_start(),
            u64::from(prod.ref_bytes_start()) + 1,
            "ref_bytes_start ({ng})"
        );
        assert_eq!(ng.ref_bytes(), prod.ref_bytes(), "ref_bytes ({ng})");
        assert_eq!(ng.ref_tract(), prod.ref_tract(), "ref_tract ({ng})");
        assert_eq!(ng.left_flank(), prod.left_flank(), "left_flank ({ng})");
        assert_eq!(ng.right_flank(), prod.right_flank(), "right_flank ({ng})");
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
        let ours = admit(intervals.to_vec(), chrom, seq, ng_p);
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
    /// three `upper()` calls in `admit`/`finish_locus` were identity functions in
    /// 100% of the suite — including both differentials — and deleting them would
    /// have gone unnoticed.
    ///
    /// This is not an exotic input. Real references **soft-mask repeats**, which
    /// is precisely the sequence this module exists to process, so lower case is
    /// the mainstream case. And the failure is silent: `Locus` compares by value,
    /// so a lower-case motif or `ref_bytes` reaching Milestone D simply fails to
    /// match the catalog's — a wrong locus, not a crash.
    #[test]
    fn admit_upper_cases_a_soft_masked_tract_and_its_flanks() {
        let left = b"cgcgc";
        let tract = b"atatatatatatatat"; // 8 copies, soft-masked
        let right = b"gcgcg";
        let contig = [left.as_ref(), tract.as_ref(), right.as_ref()].concat();

        let loci = admit(vec![iv(5, 21, 2, 100)], "chr1", &contig, &params());
        assert_eq!(loci.len(), 1, "a soft-masked tract is still a locus");
        let l = &loci[0];
        assert_eq!(l.motif().as_bytes(), b"AT", "motif upper-cased");
        assert_eq!(l.ref_tract(), b"ATATATATATATATAT", "tract upper-cased");
        assert_eq!(l.left_flank(), b"CGCGC", "left flank upper-cased");
        assert_eq!(l.right_flank(), b"GCGCG", "right flank upper-cased");
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
        let mixed_loci = admit(vec![iv(5, 21, 2, 100)], "chr1", &mixed, &params());
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

        let loci = admit(vec![detected], "chr1", &contig, &params());
        assert_eq!(loci.len(), 1);
        let l = &loci[0];
        assert_eq!(
            l.start(),
            9,
            "trimmed start: 0-based 5 (detected) + 3 (st) + 1 (1-based), NOT 6"
        );
        assert_eq!(l.end(), 20, "end unmoved");
        assert_eq!(l.ref_tract(), tract, "the clean tract, junk trimmed off");
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
        let loci = admit(vec![detected], "chr1", &contig, &ng_p);
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
            admit(vec![detected], "chr1", &contig, &ng_hi).is_empty(),
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
            admit(vec![iv(5, 21, 2, 100)], "chr1", &contig, &ng_p).len(),
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
            admit(vec![iv(5, 21, 2, 49)], "chr1", &contig, &ng_p).is_empty(),
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
            admit(vec![iv(5, 21, 2, 50)], "chr1", &contig, &ng_p).len(),
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
            admit(
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
    /// `scanner_parity`'s `#[cfg(test)]` module), so they are restated here as the
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
            admit(vec![cag], "chr1", &contig, &params()).len(),
            1,
            "period 3 is inside the default 2..=6 scope"
        );

        // Take period 3 off the STR path: same tract, no locus.
        let dinucs_only = SsrAdmissionParams {
            periods: PeriodRange::new(2, 2).unwrap(),
            ..params()
        };
        assert!(
            admit(vec![cag], "chr1", &contig, &dinucs_only).is_empty(),
            "a 2..=2 scope must exclude the period-3 tract — the knob has to bite"
        );

        // And the floor moves too: period 1 is off by default, on if asked.
        let mut poly = vec![b'T'; 5];
        poly.extend(std::iter::repeat_n(b'A', 30));
        poly.extend_from_slice(b"TTTTT");
        let homopolymer = iv(5, 35, 1, 100);
        assert!(
            admit(vec![homopolymer], "chr1", &poly, &params()).is_empty(),
            "period-1 homopolymers are off the STR path by default"
        );
        let with_homopolymers = SsrAdmissionParams {
            periods: PeriodRange::new(1, 6).unwrap(),
            ..params()
        };
        assert_eq!(
            admit(vec![homopolymer], "chr1", &poly, &with_homopolymers).len(),
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
            admit(pair.to_vec(), "chr1", &contig, &wide).is_empty(),
            "30 bp apart is inside a 50 bp flank radius: both are bundle members"
        );

        // Radius 5: the gap clears it → neither is bundled → both admitted.
        let narrow = matched_params(0.8, 0, 5).0;
        assert_eq!(
            admit(pair.to_vec(), "chr1", &contig, &narrow).len(),
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
        let loci = admit(vec![di, tri], "chr1", &contig, &picky);
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
        let loci = admit(vec![di, tri], "chr1", &contig, &mirrored);
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
        let _ = admit(vec![iv(5, 21, 2, 100)], "chr1", &contig, &too_wide);
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
        let _ = admit(vec![iv(5, 21, 2, 100)], "chr1", &contig, &nan);
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
            admit(vec![tract], "chr1", &contig, &params()).len(),
            1,
            "8 copies clears the default floor of 5"
        );

        // Raise the floor above the tract: dropped.
        let strict = SsrAdmissionParams {
            min_copies: MinCopies::uniform(9),
            ..params()
        };
        assert!(
            admit(vec![tract], "chr1", &contig, &strict).is_empty(),
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
        let a = Candidate::from(iv(100, 130, 2, 100));
        let gap_exactly = Candidate::from(iv(180, 210, 2, 100));
        assert!(
            !is_close(&a, &gap_exactly, 50),
            "a gap of exactly `thresh` is not close (strict <)"
        );
        assert!(
            is_close(&a, &Candidate::from(iv(179, 209, 2, 100)), 50),
            "a gap of thresh - 1 is close"
        );

        // Overlapping, equal length: start-start, gap, and end-end are ALL exactly
        // `thresh`, so this one fixture pins three clauses at once.
        let long = Candidate::from(iv(100, 200, 2, 100));
        let shifted = Candidate::from(iv(150, 250, 2, 100));
        assert!(
            !is_close(&long, &shifted, 50),
            "three clauses sit exactly at `thresh`; strict < keeps them all false"
        );
        assert!(
            is_close(&long, &Candidate::from(iv(149, 249, 2, 100)), 50),
            "shift one closer and all three fire"
        );
    }

    /// **The pre-filter's period floor is load-bearing, not a tidy-up.**
    ///
    /// Found by mutation: deleting `iv.period >= p.periods.min()` left the whole
    /// suite green, and **no differential can ever catch it** — `assert_agrees`
    /// runs `prefilter` once and hands the same set to both sides (review finding
    /// M4), so the pre-filter is compared against nothing. Only a direct test
    /// reaches it.
    ///
    /// Why it matters: the floor runs *before* redundancy elimination, and
    /// **period 1 divides every period**. So a homopolymer that survives the floor
    /// stage eliminates every real STR it overlaps, as its "period-multiple". That
    /// is the poly-A cascade production's `iv.period >= 2` exists to prevent —
    /// the same failure the `admit`-side period-1 drop guards from the other
    /// direction, and the reason both are ordered before their bundling/redundancy
    /// step.
    #[test]
    fn prefilter_drops_period_one_before_it_can_eliminate_a_real_str() {
        // A 30 bp poly-A: 30 copies at period 1, which clears period 1's floor of
        // 10 — so the floor alone will not save us.
        let homopolymer = iv(100, 130, 1, 100);
        // A real period-2 tract over the same span: 15 copies, clears its floor.
        // `2 % 1 == 0` and they overlap, so without the period gate the
        // homopolymer eliminates it as a period-multiple.
        let real_str = iv(100, 130, 2, 100);

        assert_eq!(
            prefilter(&[homopolymer, real_str], &params()),
            vec![real_str],
            "the homopolymer must be dropped BEFORE redundancy elimination, or it \
             takes the real STR with it"
        );

        // The gate follows the knob, not a hardcoded 2: admit period 1 and the
        // cascade is back — which is a real cost of that experiment (spec §10),
        // and should be visible rather than surprising.
        let with_homopolymers = SsrAdmissionParams {
            periods: PeriodRange::new(1, 6).unwrap(),
            ..params()
        };
        assert_eq!(
            prefilter(&[homopolymer, real_str], &with_homopolymers),
            vec![homopolymer],
            "at a 1..=6 scope the homopolymer is kept and eliminates the period-2 \
             tract as its multiple — the documented cost of putting period 1 on \
             the STR path"
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

            total_loci += admit(cleaned.clone(), &name, seq, &params()).len();
            contigs += 1;
        }
        assert!(contigs > 0, "the fixture must have contigs");
        assert!(
            total_loci > 0,
            "the fixture must admit loci — a differential over two empty sets proves nothing"
        );
    }
}
