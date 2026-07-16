//! ng step 3 — the typed-region generator: walk the reference and cut it into
//! consecutive typed regions, each a span plus *what the sequence there is*.
//! Design: `doc/devel/ng/spec/typed_regions.md` (spec) and
//! `doc/devel/ng/arch/typed_regions.md` (types & interfaces).
//!
//! **Build status (incremental).** [`admission`] — ng's own copy of the STR
//! admission policy (Milestone A) — plus the walk's types: [`TypedRegion`],
//! [`RegionKind`], [`TypedRegionConfig`], [`TypedRegionCounts`],
//! [`TypedRegionError`] (C2). `GenomeRegions` (C3) and `TypedRegionIterator` plus
//! the walk itself (D–E) are still to come, so **nothing here has logic yet** —
//! these are the shapes the walk will fill.
//!
//! **A folder, not a file, and not because of a bake-off** (there is none —
//! spec §6). The admission port is a second concern with its own dense test
//! suite, so it gets its own module beside the walk.
//!
//! ## Production is frozen; ng owns its copies
//!
//! Step 3 needs an STR admission policy that is windowed, 1-based/`u64`,
//! driven by `RepeatInterval`s, all-knobs, and that hands bundle members back
//! instead of dropping them. `ssr::catalog::postprocess::build_loci` is none of
//! those things, and **reshaping it in place is not on the table** (spec
//! Revision 2026-07-16, owner): production stays exactly as it is, so that it
//! remains an *independent yardstick* for the experiments ng exists to run.
//!
//! So [`admission`] is a **port**: the logic transcribed unchanged, the shape
//! ng's. What sharing one function used to guarantee for free, a test now pins
//! — see [`admission`]'s differential against production (spec §8.0).

pub mod admission;

use std::path::Path;

use crate::ng::ref_seq::RefSeqError;
use crate::ng::tandem_repeat::{PeriodRange, RepeatInterval, ScanParams, find_tandem_repeats};
use crate::ng::types::{Bp, ContigId, GenomeRegion, Position};
use crate::regions::{BedError, ContigBounds, RegionSet};
use admission::{Locus, SsrAdmissionParams};

// ---------------------------------------------------------------------
// What to walk
// ---------------------------------------------------------------------

/// The set of genome regions to walk — sorted, non-overlapping, coalesced,
/// clamped, in genomic order.
///
/// **Wraps production's `RegionSet` read-only; reimplements nothing.** That type
/// already parses BED, coalesces overlapping and adjacent spans, clamps to contig
/// lengths, resolves names against the contig table, and drops zero-length
/// contigs — and it is the same code the production caller's `--regions` runs, so
/// ng and production agree on what a BED *means* by construction rather than by
/// coincidence. `src/regions.rs` is not edited (spec Revision): this wrapper adds
/// ng's width and ng's names, and nothing else.
///
/// **A user BED is not a special case**, which `regions.rs` settled first: *"'Whole
/// genome' is not a special case — it is the region set whose every region covers
/// an entire contig."* [`Self::whole_contigs`] is the default, not a bypass.
///
/// ## The conversion this owns — and it is smaller than the spec expected
///
/// Spec §4 called this "the one conversion seam", *"`GenomeRegions` widening (and
/// **rebasing**) `RegionSet`'s `u32`"*. **There is no rebasing.** `regions::Region`
/// is already **1-based inclusive** — its own doc says so, and its invariant is
/// `1 <= start <= end`. So production and ng already agree on the base, and the
/// only conversion here is `u32` → `u64`, which is lossless and infallible.
///
/// That is worth stating rather than quietly enjoying: the spec anticipated an
/// off-by-one seam at ng's busiest boundary and there is none, because the
/// production author had already made the same call for the same reason.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomeRegions {
    inner: RegionSet,
}

impl GenomeRegions {
    /// One full-length span per contig — **the default** (spec §2.5).
    ///
    /// Zero-length contigs contribute no span, so they never reach the walk
    /// (`RegionSet`'s rule, and the reason spec §2.3 can say "zero-length contigs
    /// never reach the walk" without a guard of its own).
    pub fn whole_contigs(contigs: &[ContigBounds]) -> Self {
        Self {
            inner: RegionSet::whole_contigs(contigs),
        }
    }

    /// Parse a BED and resolve it against the contig table.
    ///
    /// Every failure mode — a short line, non-numeric coordinates, `end <= start`,
    /// an unknown contig name, a span past a contig's end — is `RegionSet`'s to
    /// reject, **up front**. That is what lets `TypedRegionError` have exactly one
    /// variant and `TypedRegionIterator::over_regions` be infallible (spec §8.2):
    /// by the time the walk runs, the only thing left that can fail is reading the
    /// reference.
    pub fn from_bed_path(bed: &Path, contigs: &[ContigBounds]) -> Result<Self, BedError> {
        Ok(Self {
            inner: RegionSet::from_bed_path(bed, contigs)?,
        })
    }

    /// The regions, in genomic order, as ng's [`GenomeRegion`].
    ///
    /// The `u32` → `u64` widening lives here and only here (above).
    pub fn iter(&self) -> impl Iterator<Item = GenomeRegion> + '_ {
        self.inner.iter().map(|r| GenomeRegion {
            contig: ContigId(r.chrom_id),
            start: Position(u64::from(r.start)),
            end: Position(u64::from(r.end)),
        })
    }

    /// How many regions will be walked.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Whether there is nothing to walk.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

// ---------------------------------------------------------------------
// The walk's output
// ---------------------------------------------------------------------

/// A genome region plus **what the sequence there is** — the walk's output.
///
/// `region` is a field, not a per-variant repeat: every typed region has one,
/// *structurally*, and the partition invariant (spec §2.3) reads it off directly.
/// It is also the one place ng's 1-based base is stated for this step (spec §4).
///
/// Rejected: a `region()` accessor over four variants, which makes "every typed
/// region has a region" a convention rather than a fact of the type. The evidence
/// it fails: an earlier draft written against the accessor design stated the
/// invariant in **0-based** phrasing — in the very property the spec calls its
/// spine.
#[derive(Debug, Clone, PartialEq)]
pub struct TypedRegion {
    pub region: GenomeRegion,
    pub kind: RegionKind,
}

/// What the sequence in a region **is** — one of four, and exactly one of them is
/// a *genetic* object (spec §1.1).
///
/// | kind | is | why |
/// |---|---|---|
/// | `SsrLocus` | a **locus** | the only kind the reference alone hands you as a genetic object |
/// | `SsrBundle` | a region | real repeats, none with clean flanks, so no locus can be named |
/// | `Generic` | a region | nothing more specific can be said from the reference alone |
/// | `Satellite` | a region | a tandem array too long to be a microsatellite |
///
/// **This step types regions; it does not decide their fate.** What a consumer
/// then does with each — genotype it, pile it up, mask it, skip it — is a
/// decision downstream (spec §1).
///
/// `Generic` and `Satellite` carry nothing because they *are* just their region.
/// That is not only tidiness: an open generic run costs **two coordinates however
/// many megabases it spans**, which is what makes spec §2.1's emission rule
/// affordable.
#[derive(Debug, Clone, PartialEq)]
pub enum RegionKind {
    /// ng's own [`Locus`] — motif, borders, purity, and the embedded
    /// flank+tract+flank bases. No wrapper: it is 1-based like the rest of ng
    /// (spec §4), and [`TypedRegion`] already carries the region.
    SsrLocus(Locus),
    /// A cluster of repeats none of which has clean flanks (spec §2.4). Carries
    /// the tracts as coordinates — enough to see the structure (each interval has
    /// its period) without this step pre-deciding what it is for. The hull is the
    /// [`TypedRegion`]'s own region.
    ///
    /// `>= 2` members, coordinate-ordered. **This variant is why bundles exist as
    /// a type at all**: production *deletes* these records, so their bases become
    /// a hole nobody accounts for; carrying them lets the decision be taken later,
    /// with the evidence in hand (spec §1, §10).
    SsrBundle { tracts: Box<[RepeatInterval]> },
    /// Nothing more specific can be said from the reference alone — **the
    /// default**, not a leftover. The other three are exceptions carved out of it
    /// (spec §2.2), and a repeat admission turns down for any reason other than
    /// bundling lands back here rather than becoming a hole.
    Generic,
    /// A tandem array longer than `max_repeat_len` — an array, not a
    /// microsatellite. A **typing** claim, and `max_repeat_len` is its parameter,
    /// not a constant of nature (spec §2.1, §10).
    Satellite,
}

// ---------------------------------------------------------------------
// Config and counts
// ---------------------------------------------------------------------

/// The walk's policy. Mirrors `ReadFilterConfig`'s shape (`read_filtering.md`
/// §2.4): defaults as named consts, `Default` = what the lab runs, no dormant
/// knobs.
///
/// `Default` is **the catalog's settings, for spec §8's comparability only** — not
/// an endorsement of them (spec §5.2). The catalog is a yardstick, not an
/// authority.
#[derive(Debug, Clone, PartialEq)]
pub struct TypedRegionConfig {
    /// The period range the **scanner** looks for.
    ///
    /// **Not the same knob as [`SsrAdmissionParams::periods`], and the difference
    /// is load-bearing.** The scanner is deliberately permissive (it scans 1..=6
    /// and emits every period-multiple of every tract); admission is strict (2..=6
    /// by default). Collapsing them would either blind the scanner to the
    /// homopolymers the pre-filter must *see in order to drop them before
    /// redundancy elimination* — the poly-A cascade (`admission::prefilter`) — or
    /// silently widen what gets admitted. Two ranges, two jobs.
    pub periods: PeriodRange,
    /// The scanner's scoring weights.
    pub scan: ScanParams,
    /// The satellite cap **and** the window's detection margin — one field,
    /// because they must be the same number (spec §2.6): the margin exists to
    /// capture whole any repeat that is not a satellite, so it is exactly the
    /// length at which a repeat becomes one.
    pub max_repeat_len: Bp,
    /// The walk's memory unit. **Must not change the output** (spec §2.3) — it is
    /// a memory knob and nothing else, which is why window-invariance is an
    /// acceptance test rather than a nicety.
    pub window_bp: Bp,
    /// Admission's rules — all of them (spec §5).
    pub admission: SsrAdmissionParams,
}

/// The default scan period floor: **1**, wider than admission's.
///
/// The scanner must *see* period-1 homopolymers even though nothing admits them,
/// because `prefilter` drops them **before** redundancy elimination — and period 1
/// divides every period, so a homopolymer that survives that stage eliminates any
/// real STR it overlaps (`admission::prefilter`). Not seeing them is not the same
/// as dropping them.
pub const DEFAULT_SCAN_MIN_PERIOD: u8 = 1;

/// The default scan period ceiling: **6**, the microsatellite ceiling.
pub const DEFAULT_SCAN_MAX_PERIOD: u8 = 6;

/// The satellite cap and detection margin: **1 kb**, comfortably above any real
/// STR locus. Spec §10 asks whether it is right; it is a parameter, so the
/// experiment can answer.
pub const DEFAULT_MAX_REPEAT_LEN: u64 = 1000;

/// The walk's window core: **100 kb**. Memory only — it must never change the
/// output.
pub const DEFAULT_WINDOW_BP: u64 = 100_000;

impl Default for TypedRegionConfig {
    fn default() -> Self {
        Self {
            periods: PeriodRange::new(DEFAULT_SCAN_MIN_PERIOD, DEFAULT_SCAN_MAX_PERIOD)
                .expect("the default scan period range is valid: 1 <= 1 <= 6"),
            scan: ScanParams::default(),
            max_repeat_len: Bp(DEFAULT_MAX_REPEAT_LEN),
            window_bp: Bp(DEFAULT_WINDOW_BP),
            admission: SsrAdmissionParams::default(),
        }
    }
}

/// Why a repeat did not become a locus. The breakdown spec §3.1 needs, because
/// one total cannot separate a purity rejection from a copy-floor one — and that
/// distinction is exactly what spec §10's routing question turns on.
///
/// No `str_bundle` variant: since spec §2.4 that is a **route**, not a rejection.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RejectionCounts {
    pub copy_floor: u64,
    pub purity: u64,
    pub compound: u64,
    pub no_clean_trim: u64,
    pub flank_clamped: u64,
}

/// The walk's running tally — readable mid-walk, complete once exhausted.
///
/// **"No silent caps"**: a base typed away from the STR path must be accounted
/// for. This is the live-caller view of the catalog's measured ~35% STR coverage
/// gap.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TypedRegionCounts {
    /// Requested regions walked.
    pub spans: u64,
    pub ssr_loci: u64,
    pub ssr_bundles: u64,
    /// **The number spec §10's bundle question needs, and has never had** —
    /// because the answer was previously deleted uncounted (production drops
    /// bundle members without recording them).
    pub ssr_bundle_bp: u64,
    pub generic: u64,
    pub satellites: u64,
    pub satellite_bp: u64,
    /// Repeat coverage that yielded no locus, **in bp and broken out by reason**
    /// ([`Self::rejected_by_reason`]).
    ///
    /// In bp, not per repeat, because a per-repeat count answers the wrong
    /// question twice: admission trims every survivor, so a repeat that admits one
    /// locus and sheds 200 bp contributes nothing to a per-repeat counter (spec
    /// §3.1).
    pub repeat_bp_with_no_locus: u64,
    pub rejected_by_reason: RejectionCounts,
}

// ---------------------------------------------------------------------
// The walk — resident (D1)
// ---------------------------------------------------------------------

/// Cut one whole contig into typed regions, holding it entirely in memory.
///
/// **The five steps** (spec §2.1), and the first three are ng's port of the
/// catalog's implementation (spec §5):
///
/// 1. **Detect** — [`find_tandem_repeats`] → raw, overlapping candidates. No policy.
/// 2. **Clean** — [`admission::prefilter`]: per-period copy floor, then
///    period-multiple redundancy. **Not optional** (spec §5b).
/// 3. **Admit** — [`admission::admit`] → the STR loci *and* the tracts it set
///    aside as bundle members.
/// 4. **Cap** — merge the cleaned intervals into coverage runs; a run over
///    `max_repeat_len` is a satellite. Admitted loci inside one are dropped.
/// 5. **Partition** — emit `SsrLocus` at each surviving tract, `Satellite` at each
///    run, `Generic` across everything else.
///
/// # Why this exists when the windowed walk is what ships
///
/// It is **D3's oracle** (impl plan, Milestone D): the windowed walk is proven by
/// *matching this*, which is exactly the window-invariance spec §2.3 demands.
/// Simplest implementation first, as the next one's yardstick — so a windowing bug
/// shows up as a disagreement with a version that has no windows to get wrong.
///
/// # A rejected repeat is generic territory, not a hole (spec §2.2)
///
/// The generic path is the **default**; the other three are exceptions carved out
/// of it. So a repeat admission turns down for being impure, or low-copy, or
/// compound simply stays `Generic` — it is not a bundle, and it is certainly not a
/// hole. Only the *flank test* makes a bundle: a repeat with another repeat within
/// `flank_bp` of it, which is exactly the set [`admission::admit`] hands back.
pub fn partition_resident(
    chrom: &str,
    contig: ContigId,
    bases: &[u8],
    config: &TypedRegionConfig,
) -> Vec<TypedRegion> {
    let contig_len = bases.len() as u64;
    if contig_len == 0 {
        // A zero-length contig has no 1-based position to cover, so it has no
        // regions. `GenomeRegions` drops these before the walk anyway (C3); this
        // keeps the function total for a direct caller.
        return Vec::new();
    }

    // 1. Detect.
    let raw = find_tandem_repeats(bases, config.periods, &config.scan);
    // 2. Clean.
    let cleaned = admission::prefilter(&raw, &config.admission);
    // 3. Admit — whole-contig is the degenerate window (spec §5a).
    let admitted = admission::admit(
        cleaned.clone(),
        chrom,
        bases,
        Position(1),
        Bp(contig_len),
        &config.admission,
    );
    // 4. Cap: coverage runs over the *cleaned* intervals, then the satellite test.
    //
    //    Over the cleaned set, not the raw one (spec §2.4, §8): the raw scanner is
    //    deliberately permissive, so capping its coverage would let detector noise
    //    declare a satellite and silently swallow the real loci underneath it.
    //
    //    **Known untested.** Mutation says `coverage_runs(&raw)` passes the whole
    //    suite, including `.cat` parity — because `raw ⊇ cleaned`, the two differ
    //    only where raw-only coverage runs ≥ `max_repeat_len` *contiguously*, and
    //    no fixture here produces that (the golden's 534 raw intervals on ctg1 do
    //    not union past 1 kb, and the crafted fixtures' noise is adjacent to
    //    nothing). So this ordering rests on inspection and on the spec's argument,
    //    not on a test — stated because an untested claim that reads like a tested
    //    one is worse than either. A discriminating fixture wants ≥ 1 kb of
    //    contiguous low-copy noise abutting a sub-cap array; if the routing
    //    experiments ever build one, pin it here.
    let runs = coverage_runs(&cleaned);
    let max_repeat_len = config.max_repeat_len.get();

    // 5. Partition. Collect the non-generic features in coordinate order, then fill
    //    every gap with `Generic`.
    let mut features: Vec<TypedRegion> = Vec::new();
    for run in &runs {
        if run.len() > max_repeat_len {
            features.push(TypedRegion {
                region: GenomeRegion {
                    contig,
                    start: Position(run.start),
                    end: Position(run.end),
                },
                kind: RegionKind::Satellite,
            });
        }
    }
    // A satellite array is typed as one object, not searched for loci inside it
    // (spec §2.1) — so anything admission produced inside one is dropped, and the
    // one-label-per-base rule needs no tie-break.
    let swallowed = |start: Position| {
        runs.iter()
            .any(|r| r.len() > max_repeat_len && r.contains(start))
    };

    for locus in admitted.loci {
        let span = GenomeRegion {
            contig,
            start: Position(locus.start()),
            end: Position(locus.end()),
        };
        if swallowed(span.start) {
            continue;
        }
        features.push(TypedRegion {
            region: span,
            kind: RegionKind::SsrLocus(locus),
        });
    }

    // Bundles (D2). `admit` set these aside — repeats too close to each other for
    // any of them to have a clean flank (spec §2.4) — and handed them back rather
    // than deleting them, which is the whole point of `Admitted::bundled`. Each
    // cluster becomes **one** region spanning the hull of its tracts: the gaps
    // between members are inside it, and rightly so — they are shorter than a
    // flank, so nothing can be anchored in them either.
    for cluster in admission::bundle_clusters(&admitted.bundled, config.admission.flank_bp) {
        // 0-based half-open → 1-based inclusive, the same one conversion as
        // everywhere else (spec §4).
        let start = Position(cluster.first().expect("non-empty cluster").start + 1);
        let end = Position(cluster.iter().map(|iv| iv.end).max().expect("non-empty"));
        if swallowed(start) {
            continue;
        }
        features.push(TypedRegion {
            region: GenomeRegion { contig, start, end },
            kind: RegionKind::SsrBundle {
                tracts: cluster.into_boxed_slice(),
            },
        });
    }

    features.sort_by_key(|f| f.region.start);

    fill_generic_gaps(features, contig, contig_len)
}

/// A merged run of repeat coverage, 1-based inclusive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CoverageRun {
    start: u64,
    end: u64,
}

impl CoverageRun {
    fn len(self) -> u64 {
        self.end - self.start + 1
    }
    fn contains(self, pos: Position) -> bool {
        self.start <= pos.get() && pos.get() <= self.end
    }
}

/// Union the intervals into maximal runs of covered bases, 1-based inclusive.
///
/// Overlapping **and abutting** runs merge: two tracts that touch cover a
/// contiguous stretch of repeat, and whether that stretch is a satellite is a
/// question about the stretch, not about how the detector happened to split it.
///
/// Input is `RepeatInterval`'s 0-based half-open; output is ng's 1-based
/// inclusive, so `[s, e)` becomes `[s + 1, e]` — the same one conversion `admit`
/// makes (spec §4).
fn coverage_runs(intervals: &[RepeatInterval]) -> Vec<CoverageRun> {
    let mut spans: Vec<CoverageRun> = intervals
        .iter()
        .filter(|iv| iv.end > iv.start)
        .map(|iv| CoverageRun {
            start: iv.start + 1,
            end: iv.end,
        })
        .collect();
    spans.sort_by_key(|s| (s.start, s.end));

    let mut out: Vec<CoverageRun> = Vec::new();
    for s in spans {
        match out.last_mut() {
            // `s.start <= last.end + 1` merges abutting runs as well as
            // overlapping ones.
            Some(last) if s.start <= last.end + 1 => last.end = last.end.max(s.end),
            _ => out.push(s),
        }
    }
    out
}

/// Fill every gap between `features` with `Generic`, so the result tiles
/// `[1, contig_len]` exactly.
///
/// **Maximality is a correctness requirement here, not tidiness** (spec §2.3): a
/// generic region is territory the pileup mints loci *inside*, so its reach is
/// bounded by the region it was handed. Split a run at *p* and an indel spanning
/// *p* is callable by neither half — it never appears, and nothing fails. Hence
/// one `Generic` per gap, however long.
///
/// `features` must be coordinate-ordered and non-overlapping.
fn fill_generic_gaps(
    features: Vec<TypedRegion>,
    contig: ContigId,
    contig_len: u64,
) -> Vec<TypedRegion> {
    let generic = |start: u64, end: u64| TypedRegion {
        region: GenomeRegion {
            contig,
            start: Position(start),
            end: Position(end),
        },
        kind: RegionKind::Generic,
    };

    let mut out = Vec::with_capacity(features.len() * 2 + 1);
    let mut pos = 1u64;
    for f in features {
        if f.region.start.get() > pos {
            out.push(generic(pos, f.region.start.get() - 1));
        }
        pos = f.region.end.get() + 1;
        out.push(f);
    }
    if pos <= contig_len {
        out.push(generic(pos, contig_len));
    }
    out
}

// ---------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------

/// A fatal, walk-level error.
///
/// **One variant**, because `GenomeRegions` owns everything else that could go
/// wrong (unknown contig, bad BED line, a region past a contig's end) and rejects
/// it up front — so by the time the walk runs, the only thing left that can fail
/// is reading the reference.
///
/// `#[non_exhaustive]`: matchers must accept future variants.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum TypedRegionError {
    /// A reference read failed mid-walk. Fatal: every window reads the reference,
    /// and a `None` meaning both "end of the walk" and "a window in chromosome 7
    /// failed" would **silently un-call the rest of the genome** (spec §8.2).
    #[error("reference read failed during the walk")]
    Reference(#[from] RefSeqError),
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `Default` is the catalog's settings, for spec §8's comparability — so the
    /// oracle compares like with like. Pinned against the **literals**, not
    /// against the consts that define it: asserting `periods.min() ==
    /// DEFAULT_SCAN_MIN_PERIOD` would compare a constant with itself and could not
    /// fail (the same tautology the A2 review caught in `SsrAdmissionParams`).
    #[test]
    fn default_config_is_the_catalogs_settings() {
        let c = TypedRegionConfig::default();
        assert_eq!(c.periods.min(), 1, "the scanner scans period 1..");
        assert_eq!(c.periods.max(), 6, "..to 6");
        assert_eq!(c.max_repeat_len, Bp(1000), "the 1 kb satellite cap");
        assert_eq!(c.window_bp, Bp(100_000), "100 kb window");
        assert_eq!(c.scan, ScanParams::default());
        assert_eq!(c.admission, SsrAdmissionParams::default());
    }

    /// **The scanner scans wider than admission admits, and that is deliberate.**
    ///
    /// Two `periods` knobs look redundant until you ask why: `prefilter` must drop
    /// period-1 homopolymers *before* redundancy elimination, and period 1 divides
    /// every period — so a homopolymer that is never *scanned* cannot be dropped,
    /// and one that survives to redundancy elimination takes every real STR it
    /// overlaps with it (the poly-A cascade). Not seeing them is not the same as
    /// dropping them.
    ///
    /// If these two ever coincide, one of the two jobs has been lost.
    #[test]
    fn the_scan_range_is_wider_than_the_admission_range() {
        let c = TypedRegionConfig::default();
        assert!(
            c.periods.min() < c.admission.periods.min(),
            "the scanner must SEE period {} so the pre-filter can drop it before \
             redundancy elimination; admission starts at {}",
            c.periods.min(),
            c.admission.periods.min()
        );
        assert_eq!(
            c.periods.max(),
            c.admission.periods.max(),
            "the ceilings agree: nothing is gained by scanning wider than the \
             widest admissible period"
        );
    }

    #[test]
    fn counts_start_at_zero() {
        let c = TypedRegionCounts::default();
        assert_eq!(c.spans, 0);
        assert_eq!(c.ssr_loci, 0);
        assert_eq!(c.ssr_bundle_bp, 0);
        assert_eq!(c.rejected_by_reason, RejectionCounts::default());
    }

    /// The error is fatal-in-stream and converts from the one thing that can fail
    /// mid-walk. `#[from]` is what lets the walk write `?`.
    #[test]
    fn a_reference_failure_converts_into_the_walk_error() {
        let err: TypedRegionError = RefSeqError::UnknownContig(ContigId(7)).into();
        assert!(matches!(err, TypedRegionError::Reference(_)));
        assert!(err.to_string().contains("reference read failed"));
    }

    // ---- GenomeRegions (C3) ---------------------------------------------

    const CONTIGS: &[ContigBounds] = &[
        ContigBounds {
            name: "chr1",
            length: 100,
        },
        ContigBounds {
            name: "chr2",
            length: 50,
        },
    ];

    /// `whole_contigs` is the default, and the spans are what `regions.rs` calls
    /// "the region set whose every region covers an entire contig" — full-length,
    /// **1-based inclusive**, one per contig, in table order.
    #[test]
    fn whole_contigs_covers_each_contig_end_to_end() {
        let g = GenomeRegions::whole_contigs(CONTIGS);
        let regions: Vec<_> = g.iter().collect();

        assert_eq!(g.len(), 2);
        assert!(!g.is_empty());
        assert_eq!(regions[0].contig, ContigId(0));
        assert_eq!(regions[0].start, Position(1), "1-based: starts at 1, not 0");
        assert_eq!(
            regions[0].end,
            Position(100),
            "inclusive: the last base IS 100"
        );
        assert_eq!(regions[0].len(), 100, "a 100 bp contig walks 100 bases");
        assert_eq!(regions[1].contig, ContigId(1));
        assert_eq!(regions[1].end, Position(50));
        assert_eq!(regions[1].len(), 50);
    }

    /// **No rebasing happens here, and that is the finding.** Spec §4 expected this
    /// seam to widen *and rebase*; `regions::Region` is already 1-based inclusive
    /// (its own invariant is `1 <= start <= end`), so only the width converts.
    ///
    /// This test is the guard on that: if production's base ever moved, the
    /// coordinates below would shift by one and ng's whole 1-based contract would
    /// quietly break at its busiest boundary.
    #[test]
    fn the_seam_widens_but_does_not_rebase() {
        let production = RegionSet::whole_contigs(CONTIGS);
        let ours = GenomeRegions::whole_contigs(CONTIGS);

        for (p, n) in production.iter().zip(ours.iter()) {
            assert_eq!(
                u64::from(p.start),
                n.start.get(),
                "start is carried across verbatim — production is already 1-based"
            );
            assert_eq!(u64::from(p.end), n.end.get(), "end likewise");
            assert_eq!(p.chrom_id, n.contig.get(), "and the id is the same index");
        }
    }

    /// Zero-length contigs contribute no span, so they never reach the walk —
    /// `RegionSet`'s rule, inherited. This is why spec §2.3 can assert "zero-length
    /// contigs never reach the walk" without the walk guarding for it.
    #[test]
    fn a_zero_length_contig_is_dropped_before_the_walk_sees_it() {
        let contigs = &[
            ContigBounds {
                name: "empty",
                length: 0,
            },
            ContigBounds {
                name: "chr1",
                length: 10,
            },
        ];
        let g = GenomeRegions::whole_contigs(contigs);
        let regions: Vec<_> = g.iter().collect();
        assert_eq!(regions.len(), 1, "the empty contig contributes nothing");
        assert_eq!(
            regions[0].contig,
            ContigId(1),
            "and the ids do NOT renumber"
        );
    }

    /// A BED round-trip: ng inherits `RegionSet`'s parsing, its 0-based-BED → 1-based
    /// conversion, and its coalescing — none of which ng reimplements. The
    /// overlapping pair must come back merged.
    #[test]
    fn from_bed_path_parses_converts_and_coalesces() {
        use std::io::Write;
        std::fs::create_dir_all("tmp").unwrap();
        let dir = tempfile::tempdir_in("tmp").unwrap();
        let bed = dir.path().join("r.bed");
        {
            let mut f = std::fs::File::create(&bed).unwrap();
            // BED is 0-based half-open: [0,10) is 1-based [1,10].
            writeln!(f, "chr1\t0\t10").unwrap();
            // Overlaps the first — must coalesce into [1, 20].
            writeln!(f, "chr1\t5\t20").unwrap();
            writeln!(f, "chr2\t0\t5").unwrap();
        }
        let g = GenomeRegions::from_bed_path(&bed, CONTIGS).expect("valid bed");
        let regions: Vec<_> = g.iter().collect();

        assert_eq!(regions.len(), 2, "the chr1 pair coalesced");
        assert_eq!(regions[0].contig, ContigId(0));
        assert_eq!(
            (regions[0].start, regions[0].end),
            (Position(1), Position(20)),
            "BED 0-based [0,20) becomes 1-based inclusive [1,20]"
        );
        assert_eq!(regions[1].contig, ContigId(1));
        assert_eq!(
            (regions[1].start, regions[1].end),
            (Position(1), Position(5))
        );
    }

    // ---- D1: the resident partition -------------------------------------

    /// **The invariant — the acceptance test** (spec §2.3).
    ///
    /// Within a walked region the typed regions are **contiguous**
    /// (`start == prev.end + 1`), **non-overlapping**, **complete** (their union is
    /// the whole span), and **maximal** (no two consecutive share a kind).
    ///
    /// One property: *concatenating the regions reconstructs what was asked for,
    /// exactly.* Every way this design fails shows up as a violation — a rejected
    /// repeat left as a hole breaks completeness; a flank counted as ownership
    /// breaks non-overlap; a generic run split at a window edge breaks maximality.
    #[track_caller]
    fn assert_partitions(regions: &[TypedRegion], contig: ContigId, contig_len: u64, case: &str) {
        assert!(
            !regions.is_empty(),
            "{case}: a non-empty contig has regions"
        );
        let mut expected_start = 1u64;
        let mut prev_kind: Option<std::mem::Discriminant<RegionKind>> = None;
        for r in regions {
            assert_eq!(r.region.contig, contig, "{case}: contig");
            assert_eq!(
                r.region.start.get(),
                expected_start,
                "{case}: gap or overlap at {} (expected {expected_start}); regions: {regions:#?}",
                r.region.start.get()
            );
            assert!(
                r.region.end >= r.region.start,
                "{case}: empty region {:?}",
                r.region
            );
            let kind = std::mem::discriminant(&r.kind);
            assert_ne!(
                Some(kind),
                prev_kind,
                "{case}: two consecutive regions share a kind at {} — MAXIMALITY. \
                 For Generic this is a correctness bug, not untidiness: the pileup \
                 mints loci inside a Generic region, so a split run makes an indel \
                 across the join callable by neither half.",
                r.region.start.get()
            );
            prev_kind = Some(kind);
            expected_start = r.region.end.get() + 1;
        }
        assert_eq!(
            expected_start - 1,
            contig_len,
            "{case}: the partition must cover exactly [1, {contig_len}] — COMPLETENESS"
        );
    }

    /// A contig with one clean isolated (AT)*8 tract: Generic / SsrLocus / Generic.
    /// The smallest case that shows the partition doing its job.
    #[test]
    fn a_lone_tract_partitions_as_generic_locus_generic() {
        // 60 bp of unique sequence either side, so the default 50 bp flanks fit.
        let mut bases = b"CGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCA".to_vec();
        let tract_start_0 = bases.len();
        bases.extend_from_slice(b"ATATATATATATATAT"); // 8 copies
        bases.extend_from_slice(b"CGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCA");
        let len = bases.len() as u64;

        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "lone tract");

        let kinds: Vec<_> = regions
            .iter()
            .map(|r| match &r.kind {
                RegionKind::SsrLocus(_) => "locus",
                RegionKind::SsrBundle { .. } => "bundle",
                RegionKind::Generic => "generic",
                RegionKind::Satellite => "satellite",
            })
            .collect();
        assert_eq!(kinds, vec!["generic", "locus", "generic"]);

        // And the locus is where the tract is — 1-based, so 0-based + 1.
        let RegionKind::SsrLocus(l) = &regions[1].kind else {
            unreachable!()
        };
        assert_eq!(l.start(), tract_start_0 as u64 + 1);
        assert_eq!(l.motif().as_bytes(), b"AT");
        assert_eq!(regions[1].region.start.get(), l.start(), "region == tract");
    }

    /// A 2 kb array is **one** `Satellite`, and the locus inside it is **swallowed**
    /// — typed as one object, not searched for loci inside (spec §2.1).
    ///
    /// **The swallow has to be checked positively**, and this test's first version
    /// did not: it used `vec![b'C'; 60]` "flanks", which are a period-1 homopolymer
    /// — so the scanner found one period-2 tract spanning the *whole contig*, which
    /// starts at base 1, has no left flank, and was dropped. Nothing was swallowed
    /// because nothing was ever admitted, and the assertion passed for entirely the
    /// wrong reason (mutation caught it: "don't drop loci inside a satellite"
    /// survived). The control below is the fix: at the same settings but a cap
    /// above the array, the same bases DO yield a locus.
    #[test]
    fn a_long_array_is_one_satellite_and_swallows_the_locus_inside_it() {
        let mut bases = filler(60);
        for _ in 0..1000 {
            bases.extend_from_slice(b"AT"); // 2 kb, over the 1 kb cap
        }
        bases.extend(filler(60));
        let len = bases.len() as u64;
        let config = TypedRegionConfig::default();

        let regions = partition_resident("chr1", ContigId(0), &bases, &config);
        assert_partitions(&regions, ContigId(0), len, "satellite");

        let satellites: Vec<_> = regions
            .iter()
            .filter(|r| matches!(r.kind, RegionKind::Satellite))
            .collect();
        assert_eq!(satellites.len(), 1, "one array, ONE satellite region");
        assert!(satellites[0].region.len() >= 2000);
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::SsrLocus(_))),
            "the locus is swallowed by the satellite — the §2.4 cost, made visible"
        );

        // **The control.** Raise the cap above the array and the SAME bases admit a
        // locus — so the absence above is the cap doing its job, not admission
        // quietly rejecting the tract for some unrelated reason. This is also what
        // makes `max_repeat_len` a parameter rather than a fact of nature (spec §10).
        let uncapped = TypedRegionConfig {
            max_repeat_len: Bp(10_000),
            ..config
        };
        let regions = partition_resident("chr1", ContigId(0), &bases, &uncapped);
        assert_partitions(&regions, ContigId(0), len, "satellite, uncapped");
        assert_eq!(
            regions
                .iter()
                .filter(|r| matches!(r.kind, RegionKind::SsrLocus(_)))
                .count(),
            1,
            "above the cap the very same tract IS a locus"
        );
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::Satellite)),
            "and nothing exceeds the raised cap, so no satellite"
        );
    }

    /// Coverage runs merge **abutting** tracts, not only overlapping ones: two
    /// tracts that touch cover a contiguous stretch of repeat, and whether *that*
    /// is a satellite is a question about the stretch, not about where the detector
    /// happened to split it.
    ///
    /// Untested until mutation said so — no other fixture puts two runs exactly
    /// end to end.
    #[test]
    fn abutting_coverage_runs_merge_into_one() {
        let iv = |start, end, period| RepeatInterval {
            start,
            end,
            period,
            score: 1,
        };
        // 1-based: [1,10] and [11,20] touch → one run [1,20]. [30,40] is separate.
        let runs = coverage_runs(&[iv(0, 10, 2), iv(10, 20, 3), iv(29, 40, 2)]);
        assert_eq!(
            runs,
            vec![
                CoverageRun { start: 1, end: 20 },
                CoverageRun { start: 30, end: 40 },
            ],
            "abutting runs merge; a gap of even one base does not"
        );
        assert_eq!(runs[0].len(), 20, "inclusive length");

        // The consequence, and the reason the merge exists: two 600 bp tracts that
        // touch are a 1200 bp run — a satellite — though neither tract is.
        let runs = coverage_runs(&[iv(0, 600, 2), iv(600, 1200, 2)]);
        assert_eq!(runs.len(), 1);
        assert_eq!(
            runs[0].len(),
            1200,
            "over the 1 kb cap, though neither tract alone is"
        );
    }

    /// Aperiodic filler — **not** a homopolymer, and that matters.
    ///
    /// A `vec![b'C'; 60]` "flank" is a period-1 tract, and worse: with an `(AT)n`
    /// array between two of them the scanner finds a **single period-2 tract
    /// spanning the whole contig**, which then starts at base 1, has no left flank,
    /// and is dropped — so a test asserting "no locus here" passes for entirely the
    /// wrong reason. This filler has no repeat at any period 1..=6 (see
    /// `a_repeat_free_contig_is_one_generic_region`, which is one `Generic` over
    /// it).
    fn filler(n: usize) -> Vec<u8> {
        b"ACGTTGCAAGCTTGCA"
            .iter()
            .copied()
            .cycle()
            .take(n)
            .collect()
    }

    /// A repeat-free contig is exactly one `Generic` region. Maximality: not many.
    #[test]
    fn a_repeat_free_contig_is_one_generic_region() {
        let bases = b"ACGTTGCAAGCTTGCAACGTTGCAAGCTTGCAACGTTGCAAGCTTGCA".repeat(3);
        let len = bases.len() as u64;
        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());

        assert_partitions(&regions, ContigId(0), len, "repeat-free");
        assert_eq!(regions.len(), 1, "one Generic, not a run of them");
        assert!(matches!(regions[0].kind, RegionKind::Generic));
        assert_eq!(regions[0].region.start, Position(1));
        assert_eq!(regions[0].region.end, Position(len));
    }

    /// **A rejected repeat is generic territory, not a hole** (spec §2.2). A
    /// low-copy tract admission turns down must still be *covered* — completeness
    /// is what the invariant is for.
    #[test]
    fn a_rejected_repeat_leaves_no_hole() {
        // (AT)*3 — below the period-2 copy floor of 5.
        let mut bases = vec![b'C'; 60];
        bases.extend_from_slice(b"ATATAT");
        bases.extend(std::iter::repeat_n(b'G', 60));
        let len = bases.len() as u64;

        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "rejected repeat");
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::SsrLocus(_))),
            "below the copy floor: not a locus"
        );
        assert_eq!(
            regions.len(),
            1,
            "and the whole contig is Generic — no hole"
        );
    }

    /// A tract at position 1 has no left flank, so it is not genotypeable and lands
    /// in `Generic` (spec §8's edge list). The partition still tiles from base 1.
    #[test]
    fn a_tract_at_position_one_is_generic_and_the_partition_still_starts_at_one() {
        let mut bases = b"ATATATATATATATAT".to_vec();
        bases.extend(std::iter::repeat_n(b'C', 80));
        let len = bases.len() as u64;

        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "tract at base 1");
        assert_eq!(regions[0].region.start, Position(1));
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::SsrLocus(_))),
            "no left flank to anchor against: not a locus"
        );
    }

    #[test]
    fn an_empty_contig_has_no_regions() {
        assert!(
            partition_resident("chr1", ContigId(0), b"", &TypedRegionConfig::default()).is_empty()
        );
    }

    // ---- D2: bundle detection --------------------------------------------

    /// Build a contig with `(AT)*10` tracts at the given 0-based offsets, aperiodic
    /// filler between and around them. Each tract is 20 bp.
    fn contig_with_tracts_at(offsets: &[usize], total: usize) -> Vec<u8> {
        let mut bases = filler(total);
        for &off in offsets {
            bases[off..off + 20].copy_from_slice(b"ATATATATATATATATATAT");
        }
        bases
    }

    /// **Two tracts 10 bp apart are ONE `SsrBundle` carrying both** — not two
    /// `Generic` regions, and not a hole. Neither can have a clean flank, so
    /// neither is a locus; but they are real repeats, and production would simply
    /// have deleted them (spec §2.4).
    #[test]
    fn two_close_tracts_become_one_bundle_carrying_both() {
        // Tracts at [60,80) and [90,110): a 10 bp gap, well inside flank_bp = 50.
        let bases = contig_with_tracts_at(&[60, 90], 240);
        let len = bases.len() as u64;
        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "two close tracts");

        let bundles: Vec<_> = regions
            .iter()
            .filter_map(|r| match &r.kind {
                RegionKind::SsrBundle { tracts } => Some((r.region, tracts)),
                _ => None,
            })
            .collect();
        assert_eq!(bundles.len(), 1, "ONE bundle, not two Generic regions");
        assert_eq!(bundles[0].1.len(), 2, "carrying BOTH tracts");
        // **The hull covers both tracts** — asserted as bounds, not as an exact
        // coordinate. The detector decides where a repeat starts and stops, and a
        // ±1–2 bp boundary/phase wobble is expected of it (`scanner_parity`
        // documents the same thing against trf-mod); pinning the exact edge would
        // be testing the detector's phase, not the bundle's hull.
        let (hull, tracts) = (bundles[0].0, bundles[0].1);
        assert!(hull.start <= Position(61), "hull reaches the first tract");
        assert!(hull.end >= Position(110), "hull reaches the last tract");
        assert_eq!(
            hull.start.get(),
            tracts.iter().map(|t| t.start).min().unwrap() + 1,
            "the hull IS the tracts' span, 1-based"
        );
        assert_eq!(hull.end.get(), tracts.iter().map(|t| t.end).max().unwrap());
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::SsrLocus(_))),
            "neither tract has a clean flank, so neither is a locus"
        );
    }

    /// **Three tracts chained 30 bp apart are ONE bundle of three** — emergent
    /// transitivity. There is no separate transitive rule to implement: membership
    /// is the local flank test, and the cluster falls out of it (spec §2.4).
    ///
    /// A–B–C where A and C are 70 bp apart — further than `flank_bp` — still chain,
    /// because B is close to both.
    #[test]
    fn three_chained_tracts_become_one_bundle_of_three() {
        // [60,80), [110,130), [160,180): each gap 30 bp; A→C is 80 bp apart.
        let bases = contig_with_tracts_at(&[60, 110, 160], 300);
        let len = bases.len() as u64;
        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "three chained tracts");

        let bundles: Vec<_> = regions
            .iter()
            .filter_map(|r| match &r.kind {
                RegionKind::SsrBundle { tracts } => Some((r.region, tracts)),
                _ => None,
            })
            .collect();
        assert_eq!(bundles.len(), 1, "the chain is ONE bundle, not two");
        assert_eq!(
            bundles[0].1.len(),
            3,
            "all three, though the outer pair is further apart than flank_bp — \
             transitivity is emergent, not a rule"
        );
        // Bounds, not exact edges — the detector's ±1–2 bp phase wobble is its
        // business, not this test's.
        assert!(bundles[0].0.start <= Position(61));
        assert!(
            bundles[0].0.end >= Position(180),
            "the hull spans the whole chain, outermost tract to outermost tract"
        );
    }

    /// **Bundles do not spread.** A tract far enough from a cluster is admitted
    /// normally — the flank test is local, so a bundle does not contaminate its
    /// neighbourhood.
    #[test]
    fn a_tract_clear_of_a_bundle_is_still_admitted() {
        // A close pair at [60,80)+[90,110), and a lone tract at [400,420) — far
        // from everything, with room for its 50 bp flanks.
        let bases = contig_with_tracts_at(&[60, 90, 400], 600);
        let len = bases.len() as u64;
        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "bundle + lone tract");

        assert_eq!(
            regions
                .iter()
                .filter(|r| matches!(r.kind, RegionKind::SsrBundle { .. }))
                .count(),
            1,
            "the close pair bundles"
        );
        let loci: Vec<_> = regions
            .iter()
            .filter_map(|r| match &r.kind {
                RegionKind::SsrLocus(l) => Some(l),
                _ => None,
            })
            .collect();
        assert_eq!(
            loci.len(),
            1,
            "the lone tract is a locus — bundles don't spread"
        );
        assert_eq!(loci[0].start(), 401);
    }

    /// **A repeat rejected for any reason OTHER than the flank test is `Generic`,
    /// not a bundle** (spec §2.2). Only closeness makes a bundle; being low-copy
    /// makes you ordinary.
    #[test]
    fn a_low_copy_repeat_is_generic_not_a_bundle() {
        let mut bases = filler(240);
        // (AT)*3 — below the period-2 copy floor of 5, and isolated.
        bases[100..106].copy_from_slice(b"ATATAT");
        let len = bases.len() as u64;
        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());

        assert_partitions(&regions, ContigId(0), len, "low-copy repeat");
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::SsrBundle { .. })),
            "a copy-floor rejection is Generic territory, NOT a bundle"
        );
        assert_eq!(regions.len(), 1, "and no hole: the contig is one Generic");
    }

    /// The clustering must reproduce exactly the split that produced it — a
    /// singleton "bundle" would mean the two disagree.
    #[test]
    fn bundle_clusters_regroups_exactly_what_the_split_set_aside() {
        let iv = |start, end| RepeatInterval {
            start,
            end,
            period: 2,
            score: 100,
        };
        // Two clusters: {100,130 / 150,180} and {5000,5030 / 5040,5070}.
        let bundled = [iv(100, 130), iv(150, 180), iv(5000, 5030), iv(5040, 5070)];
        let clusters = admission::bundle_clusters(&bundled, 50);
        assert_eq!(clusters.len(), 2, "two clusters, not one and not four");
        assert_eq!(clusters[0].len(), 2);
        assert_eq!(clusters[1].len(), 2);
        assert_eq!(clusters[0][0].start, 100);
        assert_eq!(clusters[1][0].start, 5000);
    }

    /// **`.cat` parity — the oracle, and D1's real bar** (spec §8.1).
    ///
    /// The walk at the catalog's settings must reproduce the catalog: every golden
    /// locus is present, **or absent *and* inside a satellite run**. A strict
    /// subset, and that shape is *earned* by spec §2.4's ordering — the cap applies
    /// to the *cleaned* coverage, after admission, so the difference can only go
    /// one way. Capping raw coverage would make it bidirectional and untestable.
    ///
    /// **The oracle is the committed trf-mod-built golden catalog** — a different
    /// detector, a different code path, nothing ng touched. Its detector difference
    /// is already characterised (`scanner_parity`: 16/16, 15 exact, one ±1–2 bp
    /// boundary/phase wobble, one genuine scanner-only locus trf-mod's significance
    /// model rejected), so it is a yardstick rather than a confound — hence overlap
    /// matching, inherited from `scanner_parity` for exactly that reason.
    ///
    /// What this proves: the plumbing — the scan, the pre-filter, the admission
    /// call, the coordinate arithmetic. What it does **not** prove: that the
    /// catalog's settings are *right* (spec §5). A fixed-config regression test,
    /// pinned to those settings explicitly rather than to whatever `Default` is, or
    /// it starts failing the first time someone moves a floor and reads a *result*
    /// as a bug.
    #[test]
    fn the_resident_partition_reproduces_the_golden_catalog() {
        use crate::ssr::catalog::io::CatalogReader;
        use std::fs::File;
        use std::io::BufReader;
        use std::path::Path;

        let fixture = |name: &str| {
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("data")
                .join("tandem_repeat")
                .join(name)
        };

        // The golden catalog and the settings it was built at — read, never
        // written (production is frozen).
        let mut golden_reader =
            CatalogReader::new(File::open(fixture("golden.ssr_catalog.bed.gz")).unwrap()).unwrap();
        let cat_params = golden_reader.header().params.clone();
        let golden = golden_reader.read_all().unwrap();
        assert!(!golden.is_empty(), "the golden catalog must have loci");

        // ng's walk at the SAME settings, pinned explicitly.
        let config = TypedRegionConfig {
            admission: SsrAdmissionParams {
                min_purity: cat_params.min_purity,
                min_score: cat_params.min_score,
                flank_bp: u64::from(cat_params.flank_bp),
                ..SsrAdmissionParams::default()
            },
            ..TypedRegionConfig::default()
        };

        let file = File::open(fixture("synthetic_ref.fa")).unwrap();
        let mut reader = noodles_fasta::io::Reader::new(BufReader::new(file));
        let mut ours: Vec<(String, u64, u64)> = Vec::new();
        let mut satellites: Vec<(String, u64, u64)> = Vec::new();
        for (idx, result) in reader.records().enumerate() {
            let rec = result.unwrap();
            let name = String::from_utf8_lossy(rec.name()).into_owned();
            let bases = rec.sequence().as_ref();
            let regions = partition_resident(&name, ContigId(idx as u32), bases, &config);

            // The partition must hold on REAL sequence, not only the crafted
            // fixtures above.
            assert_partitions(
                &regions,
                ContigId(idx as u32),
                bases.len() as u64,
                &format!("golden contig {name}"),
            );

            for r in &regions {
                match &r.kind {
                    RegionKind::SsrLocus(l) => ours.push((name.clone(), l.start(), l.end())),
                    RegionKind::Satellite => {
                        satellites.push((name.clone(), r.region.start.get(), r.region.end.get()))
                    }
                    _ => {}
                }
            }
        }
        assert!(!ours.is_empty(), "the walk must find loci");

        // Overlap match. Production's `Locus` is 0-based half-open, ng's is 1-based
        // inclusive, so production's `[s, e)` is ng's `[s + 1, e]` (spec §4).
        let overlaps =
            |a: &(String, u64, u64), b: &(String, u64, u64)| a.0 == b.0 && a.1 <= b.2 && b.1 <= a.2;
        let mut missed = Vec::new();
        for g in &golden {
            let g1 = (
                g.chrom().to_string(),
                u64::from(g.start()) + 1,
                u64::from(g.end()),
            );
            if !ours.iter().any(|o| overlaps(&g1, o)) {
                // Absent is legal ONLY inside a satellite run — the one expected
                // divergence (the catalog applies no cap; the walk does).
                if !satellites.iter().any(|s| overlaps(&g1, s)) {
                    missed.push(format!("{}:{}-{}", g1.0, g1.1, g1.2));
                }
            }
        }
        assert!(
            missed.is_empty(),
            "every golden locus must be present, or absent AND inside a satellite \
             run. At the catalog's settings a locus missing for any other reason is \
             a machinery bug. Missing: {missed:#?}"
        );
    }

    /// Every BED failure is `RegionSet`'s to reject **up front** — which is what
    /// lets `TypedRegionError` have one variant and the walk's constructor be
    /// infallible (spec §8.2).
    #[test]
    fn a_bad_bed_is_rejected_before_any_walk_exists() {
        use std::io::Write;
        std::fs::create_dir_all("tmp").unwrap();
        let dir = tempfile::tempdir_in("tmp").unwrap();
        let bed = dir.path().join("bad.bed");
        {
            let mut f = std::fs::File::create(&bed).unwrap();
            writeln!(f, "nosuchcontig\t0\t10").unwrap();
        }
        assert!(
            GenomeRegions::from_bed_path(&bed, CONTIGS).is_err(),
            "an unknown contig is caught at construction, not mid-walk"
        );
    }
}
