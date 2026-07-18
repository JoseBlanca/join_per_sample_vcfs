//! ng step 3 — the typed-region generator: walk the reference and cut it into
//! consecutive typed regions, each a span plus *what the sequence there is*.
//! Design: `doc/devel/ng/spec/typed_regions.md` (spec) and
//! `doc/devel/ng/arch/typed_regions.md` (types & interfaces).
//!
//! **Build status (incremental).** [`segment_criteria`] — ng's own copy of the STR
//! classification policy (Milestone A) — plus the walk's types: [`TypedRegion`],
//! [`RegionKind`], [`TypedRegionConfig`], [`TypedRegionCounts`],
//! [`TypedRegionError`] (C2). `GenomeRegions` (C3) and `TypedRegionIterator` plus
//! the walk itself (D–E) are still to come, so **nothing here has logic yet** —
//! these are the shapes the walk will fill.
//!
//! **A folder, not a file, and not because of a bake-off** (there is none —
//! spec §6). The classification port is a second concern with its own dense test
//! suite, so it gets its own module beside the walk.
//!
//! ## Production is frozen; ng owns its copies
//!
//! Step 3 needs an STR classification policy that is windowed, 1-based/`u64`,
//! driven by `RepeatInterval`s, all-knobs, and that hands bundle members back
//! instead of dropping them. `ssr::catalog::postprocess::build_loci` is none of
//! those things, and **reshaping it in place is not on the table** (spec
//! Revision 2026-07-16, owner): production stays exactly as it is, so that it
//! remains an *independent yardstick* for the experiments ng exists to run.
//!
//! So [`segment_criteria`] is a **port**: the logic transcribed unchanged, the shape
//! ng's. What sharing one function used to guarantee for free, a test now pins
//! — see [`segment_criteria`]'s differential against production (spec §8.0).

pub mod segment_criteria;

/// The E3 port anchor: the whole stack, on a real multi-contig FASTA (`anchor.rs`).
#[cfg(test)]
mod anchor;

use std::path::Path;

use crate::ng::ref_seq::{ContigTable, EvictableRefSeq, RawRefSeq, RefSeqError};
use crate::ng::tandem_repeat::{
    PeriodRange, RegionSpan, RepeatInterval, ScanParams, ScannedWindow, SegmentOptions,
    WindowCursor, WindowPlan, find_tandem_repeats, scan_window,
};
use crate::ng::types::{Bp, ContigId, GenomeRegion, Position};
use crate::regions::{BedError, ContigBounds, RegionSet};
use segment_criteria::{RejectionCounts, SsrSegment, SsrSegmentCriteria};

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
/// | `SsrSegment` | a **locus** | the only kind the reference alone hands you as a genetic object |
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
    /// ng's own [`SsrSegment`] — motif, borders, purity. **Coordinates, no bases**: the
    /// bases are in the reference the caller already has open (see [`SsrSegment`]). No
    /// wrapper: it is 1-based like the rest of ng (spec §4), and [`TypedRegion`]
    /// already carries the region.
    SsrSegment(SsrSegment),
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
    /// (spec §2.2), and a repeat classification turns down for any reason other than
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
    /// **Not the same knob as [`SsrSegmentCriteria::periods`], and the floor is where
    /// the difference bites.** The scanner is deliberately permissive (it scans 1..=6
    /// and emits every period-multiple of every tract); classification is strict (2..=6 by
    /// default).
    ///
    /// The scanner must reach **period 1 even though nothing classifies it**, because
    /// `segment_criteria::prefilter` uses the period-1 interval as the *eliminator* that
    /// removes a homopolymer's own aliases: `AAAA…` is a perfect `AA`, `AAA` and
    /// `AAAAA` repeat, and only period 1 divides them all. Scan from 2 and those
    /// aliases have nothing to eliminate them — a poly-A comes back as three
    /// "repeats". **An eliminator that was never detected cannot eliminate**, so this
    /// floor is machinery, not policy; it is not the knob a caller means when it says
    /// which periods to analyse. See `prefilter`'s docs for the ordering this pairs
    /// with (fixed 2026-07-17).
    pub periods: PeriodRange,
    /// The scanner's scoring weights.
    pub scan: ScanParams,
    /// The satellite cap **and** the window's detection margin — one field,
    /// because they must be the same number (spec §2.6): the margin exists to
    /// capture whole any repeat that is not a satellite, so it is exactly the
    /// length at which a repeat becomes one.
    ///
    /// **It is also constrained from below by [`SsrSegmentCriteria::flank_bp`]**, and
    /// [`partition_windowed`] asserts it: the margin is what lets a window see a core
    /// tract's *neighbours*, and the bundle flank test is a `flank_bp`-radius question.
    /// A margin narrower than that radius would classify tracts as clean loci that a
    /// whole-contig walk bundles — silently, and differently per `window_bp`.
    pub max_repeat_len: Bp,
    /// The walk's memory unit. **Must not change the output** (spec §2.3) — it is
    /// a memory knob and nothing else, which is why window-invariance is an
    /// acceptance test rather than a nicety.
    pub window_bp: Bp,
    /// Admission's rules — all of them (spec §5).
    pub criteria: SsrSegmentCriteria,
}

/// The scan period floor: **1**, wider than classification's, and it should stay there
/// whatever classification's floor is set to.
///
/// The scanner must *see* period-1 homopolymers even though nothing classifies them:
/// `prefilter` keeps the period-1 interval through redundancy elimination precisely
/// so it can eliminate the run's own `AA` / `AAA` / `AAAAA` aliases, and drops it by
/// the period floor afterwards. **Not seeing them is not the same as dropping them**
/// — that sentence is true, but not for the reason this doc used to give: the danger
/// is not a surviving homopolymer, it is a surviving *alias* of one, which is what a
/// scanner blind to period 1 leaves behind (`segment_criteria::prefilter`, fixed
/// 2026-07-17).
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
            criteria: SsrSegmentCriteria::default(),
        }
    }
}

/// The walk's running tally — readable mid-walk, complete once exhausted
/// ([`TypedRegionIterator::counts`]).
///
/// **"No silent caps"**: a base typed away from the STR path must be accounted
/// for. This is the live-caller view of the catalog's measured ~35% STR coverage
/// gap.
///
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
    /// Repeat coverage that yielded **no locus**, in bp: every base of cleaned repeat
    /// coverage the walk did not emit as an `SsrSegment` — because it was bundled,
    /// capped as a satellite, or rejected by one of classification's gates.
    ///
    /// In bp, not per repeat, because a per-repeat count answers the wrong
    /// question twice: classification trims every survivor, so a repeat that classifies one
    /// locus and sheds 200 bp contributes nothing to a per-repeat counter (spec
    /// §3.1).
    ///
    /// **Mid-walk this leads**: coverage is counted when a window is scanned, and the
    /// loci that cancel it are subtracted as they are emitted, a block later. It is
    /// exact once the walk is exhausted, which is what the type promises.
    pub repeat_bp_with_no_locus: u64,
    /// [`Self::repeat_bp_with_no_locus`] **broken out by classification's reason** — because
    /// one total cannot separate a purity rejection from a copy-floor one, and that is
    /// exactly the distinction spec §10's routing question turns on.
    ///
    /// **It does not partition the total**, and [`RejectionCounts`] says why: overlapping
    /// rejected repeats are both charged, and bases with no locus for a reason that is not
    /// a *rejection* (bundled, capped as a satellite, out of scope) are in the total and
    /// not here. A diagnosis of classification's gates, not an account of the genome.
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
/// 2. **Clean** — [`segment_criteria::prefilter`]: per-period copy floor, then
///    period-multiple redundancy. **Not optional** (spec §5b).
/// 3. **Admit** — [`segment_criteria::classify`] → the STR loci *and* the tracts it set
///    aside as bundle members.
/// 4. **Cap** — merge the cleaned intervals into coverage runs; a run over
///    `max_repeat_len` is a satellite. Classified loci inside one are dropped.
/// 5. **Partition** — emit `SsrSegment` at each surviving tract, `Satellite` at each
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
/// of it. So a repeat classification turns down for being impure, or low-copy, or
/// compound simply stays `Generic` — it is not a bundle, and it is certainly not a
/// hole. Only the *flank test* makes a bundle: a repeat with another repeat within
/// `flank_bp` of it, which is exactly the set [`segment_criteria::classify`] hands back.
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
    let cleaned = segment_criteria::prefilter(&raw, &config.criteria);
    // 3. Admit — whole-contig is the degenerate window (spec §5a).
    let classified = segment_criteria::classify(
        cleaned.clone(),
        chrom,
        bases,
        Position(1),
        Bp(contig_len),
        &config.criteria,
    );
    // 4. Cap: coverage runs over the *cleaned* intervals, then the satellite test.
    //
    //    Over the cleaned set, not the raw one (spec §2.4, §8): the raw scanner is
    //    deliberately permissive, so capping its coverage would let detector noise
    //    declare a satellite and silently swallow the real loci underneath it.
    //
    //    **Tested since D3** (`the_satellite_cap_applies_to_the_cleaned_coverage_not_
    //    the_raw`) — this comment used to say "known untested", and D1 was right that
    //    nothing then could discriminate the two. D3's 6 kb windowing fixture can: it
    //    carries real scanner noise (1181 raw runs vs 5 cleaned) and one noise interval
    //    abuts the array, so capping the raw set moves the satellite's edge.
    //
    //    D1's prediction of what a discriminating fixture needed — "≥ 1 kb of
    //    contiguous low-copy noise abutting a sub-cap array" — was too strong. Noise
    //    merely *touching* an over-cap run is enough to show the choice; the ≥ 1 kb case
    //    is what would show the *stakes* (noise inventing a satellite and swallowing the
    //    loci under it), and no fixture produces that yet.
    let runs = coverage_runs(&cleaned);

    // 5. Partition — the whole contig as **one block** (below), then fill every gap
    //    with `Generic`.
    let features = resolve_features(&runs, classified.loci, &classified.bundled, contig, config);

    fill_generic_gaps(features, contig, contig_len)
}

/// Step 5 for one **block**: cap the runs, place the loci the surviving runs do not
/// swallow, cluster the bundle members — the non-generic features, coordinate-ordered.
/// Generic is not this function's business: it is whatever is left over, and only the
/// caller knows how far "left over" reaches.
///
/// # A block, and why the windowed walk can work one at a time
///
/// A *block* is a stretch of repeat structure separated from the next by more than
/// `flank_bp` of repeat-free sequence. Every rule here has a **radius**, and the block
/// boundary is wider than all of them: runs merge only when they abut (radius 0),
/// clustering chains members within `flank_bp` (radius `flank_bp`), and swallowing is
/// containment (radius 0). So no input outside a block can change anything inside it,
/// and the partition of a contig is the concatenation of its blocks' partitions.
///
/// That is the whole licence for [`partition_windowed`]: it resolves one block at a
/// time and carries nothing else. [`partition_resident`] calls this once, with the
/// contig as a single block — which is the degenerate case, and why it is the oracle.
///
/// **This is shared, so window-invariance does not test it.** The windowed walk is
/// proven by matching the resident one, and both bottom out here — a bug in this
/// function is invisible to that comparison, by construction. It is covered instead by
/// D1's own bar: the partition invariant and `.cat` parity. What the comparison *does*
/// prove is what only windowing can get wrong: the carries and the attribution.
///
/// # Absorption: a satellite swallows what it touches (spec §2.4a)
///
/// **The rule (owner, 2026-07-16): a microsatellite or a cluster too close to a satellite
/// is swallowed by the satellite, which expands to cover it.**
///
/// A microsatellite 20 bp from a 1 kb array is not genotypeable, and the reason is the
/// array: there is no clean flank on that side, because the flank *is* array. A satellite
/// is already the region type that says *"an array — do not look for loci in here"* (spec
/// §2.1), so it is the one that should say it here too. The alternative — exempting the
/// array from the flank test so the neighbour becomes a clean locus — is wrong for a
/// right-sounding reason (spec §2.4).
///
/// **Both kinds, because both arise, and by different routes.** A cluster reaches here
/// when the array's own tract passes classification's gates and bundles with its neighbour. A
/// bare locus reaches here when it does *not*: a satellite run is built from **cleaned
/// coverage**, while bundling only ever sees tracts that cleared the scope/score/compound
/// gates — so an array rejected by one of those still forms a satellite and never bundles
/// anything. One rule covers both.
///
/// Absorption is iterated to a **fixed point** (a grown satellite reaches further than
/// the run it came from) and it subsumes spec §2.1's swallow: containment and adjacency
/// are the same predicate ([`absorb_into`]). Bundle members are absorbed **before** they
/// are clustered — an ordering that turned out to be load-bearing twice, once for the
/// rule and once for windowing; the reasons are at the loop.
///
/// **What this replaces, and why it was not a nicety.** The old test read a cluster's
/// **start** only, so the answer depended on which *side* of the array the microsat sat:
/// on the left it emitted a bundle **overlapping** the satellite — an invalid partition;
/// on the right the hull's start fell inside the run, so the cluster was dropped whole
/// and the microsat's bases silently became `Generic`. Same situation, two different
/// wrong answers. Probed before believing it, then fixed —
/// `a_microsatellite_beside_a_satellite_is_absorbed_into_it`.
fn resolve_features(
    runs: &[CoverageRun],
    loci: Vec<SsrSegment>,
    bundled: &[RepeatInterval],
    contig: ContigId,
    config: &TypedRegionConfig,
) -> Vec<TypedRegion> {
    let max_repeat_len = config.max_repeat_len.get();
    let flank_bp = config.criteria.flank_bp;
    let mut features: Vec<TypedRegion> = Vec::new();

    // The satellites: over-cap coverage runs — as spans that can still GROW (below).
    let mut satellites: Vec<CoverageRun> = runs
        .iter()
        .copied()
        .filter(|r| r.len() > max_repeat_len)
        .collect();

    let mut loci: Vec<(CoverageRun, SsrSegment)> = loci
        .into_iter()
        .map(|l| {
            (
                CoverageRun {
                    start: l.start(),
                    end: l.end(),
                },
                l,
            )
        })
        .collect();
    // 0-based half-open → 1-based inclusive, the same one conversion as everywhere
    // else (spec §4).
    let mut members: Vec<RepeatInterval> = bundled.to_vec();

    // **Absorption** (spec §2.4a; see the fn docs). Anything a satellite overlaps *or*
    // comes within `flank_bp` of is absorbed into it, and the satellite grows to cover
    // it. Iterated to a fixed point: a satellite that has grown reaches `flank_bp`
    // further than the run it came from, so it can absorb something the run could not.
    //
    // It terminates: each pass either absorbs at least one of a finite set of features
    // or stops. And it cannot reach out of the block — absorption's reach is exactly
    // `BlockWalk::block_barrier`, which is where the block ends.
    //
    // **Bundle members are absorbed one by one, BEFORE they are clustered, and that
    // ordering is load-bearing twice over:**
    //
    // 1. It absorbs whole clusters, never part of one — so the survivors are still
    //    complete clusters and `bundle_clusters` can regroup them (its precondition:
    //    `bundled` is the concatenation of clusters in coordinate order). If a member
    //    is absorbed, every member chained to it is too: `is_close` implies "within
    //    `flank_bp`" (each of its four clauses puts a pair of endpoints inside the
    //    radius), so the grown satellite reaches the next member in the chain, and so
    //    on by induction.
    // 2. It is what makes the **windowed** walk correct at all. A window truncates a
    //    detection at its scanned slice's edge, and only a tract longer than
    //    `max_repeat_len` can be truncated — i.e. an array. A truncated member's `end`
    //    is wrong, so `is_close` cannot re-chain it, and `bundle_clusters` sees a
    //    singleton. Absorbing members first means a truncated one never reaches
    //    clustering: it is over-cap, so it lies inside its own satellite run.
    //    (`bundle_clusters`' `debug_assert` is what caught this, exactly as D2 built it
    //    to.)
    loop {
        let mut absorbed = false;
        members.retain(|iv| {
            let span = CoverageRun {
                start: iv.start + 1,
                end: iv.end,
            };
            !absorb_into(&mut satellites, span, flank_bp, &mut absorbed)
        });
        loci.retain(|(span, _)| !absorb_into(&mut satellites, *span, flank_bp, &mut absorbed));
        if !absorbed {
            break;
        }
    }

    // Bundles (D2). `classify` set these aside — repeats too close to each other for
    // any of them to have a clean flank (spec §2.4) — and handed them back rather
    // than deleting them, which is the whole point of `Classified::bundled`. Each
    // cluster becomes **one** region spanning the hull of its tracts: the gaps
    // between members are inside it, and rightly so — they are shorter than a
    // flank, so nothing can be anchored in them either.
    //
    // A surviving cluster cannot touch a satellite even though only its *members* were
    // tested: its hull adds only the gaps between members, and those are shorter than a
    // flank — far too short to hide an over-cap run — so a hull that reaches a satellite
    // has a member that reaches it.
    let clusters: Vec<(CoverageRun, Vec<RepeatInterval>)> =
        segment_criteria::bundle_clusters(&members, flank_bp)
            .into_iter()
            .map(|cluster| {
                let hull = CoverageRun {
                    start: cluster.first().expect("non-empty cluster").start + 1,
                    end: cluster.iter().map(|iv| iv.end).max().expect("non-empty"),
                };
                (hull, cluster)
            })
            .collect();

    for satellite in &satellites {
        features.push(TypedRegion {
            region: GenomeRegion {
                contig,
                start: Position(satellite.start),
                end: Position(satellite.end),
            },
            kind: RegionKind::Satellite,
        });
    }
    for (span, locus) in loci {
        features.push(TypedRegion {
            region: GenomeRegion {
                contig,
                start: Position(span.start),
                end: Position(span.end),
            },
            kind: RegionKind::SsrSegment(locus),
        });
    }
    for (hull, cluster) in clusters {
        features.push(TypedRegion {
            region: GenomeRegion {
                contig,
                start: Position(hull.start),
                end: Position(hull.end),
            },
            kind: RegionKind::SsrBundle {
                tracts: cluster.into_boxed_slice(),
            },
        });
    }

    features.sort_by_key(|f| f.region.start);
    features
}

/// Absorb `span` into any satellite it overlaps or lies within `flank_bp` of: those
/// satellites and `span` become **one** span covering all of them, gaps included.
/// Returns whether it was absorbed (and sets `absorbed`, the fixed-point loop's flag).
///
/// # Why the gap, and not `segment_criteria::is_close`
///
/// This asks a **flank** question — *are there `flank_bp` clean bases between this
/// feature and the array?* — so the gap is the measure, and `< flank_bp` is strict, as
/// classification's own gap clause is (`is_close_is_strict_at_the_threshold`).
///
/// `is_close` itself is the wrong tool here despite testing the same relation between
/// two *tracts*: it is four `abs_diff` clauses ported from GangSTR, and three of them
/// compare start-to-start and end-to-end — meaningful between two 25 bp tracts, noise
/// against a span that can be 2 Mb long, where a locus 20 bp past the end is millions
/// of bases from the start. The relation is the same; the predicate cannot be.
///
/// The gaps swept in are deliberate: fewer than `flank_bp` bases between two repeats is
/// sequence nothing can be anchored in — the same reasoning that puts a cluster's
/// internal gaps inside its hull (spec §2.4).
fn absorb_into(
    satellites: &mut Vec<CoverageRun>,
    span: CoverageRun,
    flank_bp: u64,
    absorbed: &mut bool,
) -> bool {
    // `span.start <= s.end + flank_bp` is "the gap on this side is < flank_bp", and it
    // is also true whenever the two overlap — so containment (a locus *inside* a
    // satellite, spec §2.1's swallow) and adjacency are one rule, not two. They were
    // two, and reading only the hull's start is what made the answer depend on which
    // side of the array the feature sat.
    let touching: Vec<usize> = satellites
        .iter()
        .enumerate()
        .filter(|(_, s)| span.start <= s.end + flank_bp && s.start <= span.end + flank_bp)
        .map(|(i, _)| i)
        .collect();
    if touching.is_empty() {
        return false;
    }

    let mut union = span;
    // Back to front: removing by index cannot disturb an index still to come.
    for &i in touching.iter().rev() {
        let s = satellites.remove(i);
        union.start = union.start.min(s.start);
        union.end = union.end.max(s.end);
    }
    satellites.push(union);
    merge_runs(satellites);
    *absorbed = true;
    true
}

/// Union overlapping and abutting runs in place, leaving them ascending.
///
/// Shared by the satellite absorption above and by [`coverage_runs`] — the same rule
/// (touching spans are one span), so the same code.
fn merge_runs(runs: &mut Vec<CoverageRun>) {
    runs.sort_by_key(|s| (s.start, s.end));
    let mut merged: Vec<CoverageRun> = Vec::with_capacity(runs.len());
    for s in runs.drain(..) {
        match merged.last_mut() {
            // `s.start <= last.end + 1` merges abutting runs as well as overlapping ones.
            Some(last) if s.start <= last.end + 1 => last.end = last.end.max(s.end),
            _ => merged.push(s),
        }
    }
    *runs = merged;
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
}

/// Union the intervals into maximal runs of covered bases, 1-based inclusive.
///
/// Overlapping **and abutting** runs merge: two tracts that touch cover a
/// contiguous stretch of repeat, and whether that stretch is a satellite is a
/// question about the stretch, not about how the detector happened to split it.
///
/// Input is `RepeatInterval`'s 0-based half-open; output is ng's 1-based
/// inclusive, so `[s, e)` becomes `[s + 1, e]` — the same one conversion `classify`
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
    merge_runs(&mut spans);
    spans
}

// ---------------------------------------------------------------------
// The walk — windowed (D3)
// ---------------------------------------------------------------------

/// Cut one contig into typed regions **without ever holding it in memory** — the walk
/// that ships (spec §2.3, §2.6).
///
/// Identical output to [`partition_resident`] at every `window_bp`: the window is a
/// memory knob and nothing else, which is why window-invariance is an acceptance test
/// (`windowed_matches_the_resident_oracle_*`) and not a nicety.
///
/// # How it can be windowed at all: blocks and three carries
///
/// [`resolve_features`] explains the licence — every rule the partition applies has a
/// radius, and a *block* (repeat structure bounded by more than `flank_bp` of
/// repeat-free sequence) is wider than all of them. So the walk holds **one open
/// block**, resolves it the moment the next feature proves it closed, and carries:
///
/// - the **open coverage run** — a satellite is longer than any window by definition,
///   so its verdict cannot be reached inside one (`CoverageRun`, extended window to
///   window);
/// - the **open bundle cluster** — its members chain across a boundary just as happily
///   as inside one; and
/// - the **open generic run** — one `u64`, [`BlockWalk::emitted_upto`]. This is what
///   makes maximality affordable: a megabase of generic sequence is not accumulated,
///   it is *not yet emitted*, and it costs one coordinate until the next feature ends
///   it.
///
/// # Where each fact comes from
///
/// - **`contig_len` comes from the reference's contig table** ([`ContigTable`]), never
///   from the window. This discharges spec §2.6's provenance obligation, outstanding
///   since A3: `classify`'s guard classifies a caller passing the *window's own end* as
///   `contig_len` — production's exact mistake, and arithmetically legal, so no
///   signature can catch it. Only a caller holding both the reference and the window
///   can, and this is that caller. Get it wrong and every locus within `flank_bp` of
///   every window boundary silently vanishes, a different set for each `window_bp`.
/// - **The window geometry comes from [`scan_windowed`]**, not from arithmetic here:
///   one copy of the core/margin rules, per spec §6.1.
/// - **The bases are raw** ([`RawRefSeq`]), because `SsrSegment` compares by value and
///   canonical bytes would make every IUPAC-carrying locus unequal to the catalog's
///   (spec §6).
///
/// # One fetch per window, and one margin
///
/// Admission runs over **the slice the scan already read**. It needs each tract's own
/// bases (motif, purity), which are in that slice by construction — the interval came
/// from scanning it — and it answers the flank question by arithmetic against
/// `contig_len`, which no slice could answer anyway.
///
/// So there is one margin, spec §2.6's: `max_repeat_len`, what makes a repeat *segment*
/// identically. Until 2026-07-17 there was a second — `flank_bp` on top of it — and a
/// second fetch of every window to get it, because `SsrSegment` embedded `tract ± flank_bp`
/// of bases and a tract at the slice's own edge could not supply them. Removing the
/// payload removed the margin, the re-fetch, and the two panics that policed it.
/// **Collected**, so it holds the contig's regions after all — this is the walk's shape
/// for a test or a small job, and [`TypedRegionIterator`] is the one that ships. It is a
/// wrapper over that iterator, not a second implementation: every window-invariance test
/// below runs through this and therefore through it.
pub fn partition_windowed<R>(
    reference: R,
    contig: ContigId,
    config: &TypedRegionConfig,
) -> Result<Vec<TypedRegion>, TypedRegionError>
where
    R: RawRefSeq + ContigTable + EvictableRefSeq,
{
    TypedRegionIterator::over_contig(reference, contig, config.clone())?.collect()
}

/// Walks the reference region by region, in genomic order, gap-free — **step 3's public
/// surface** (arch §interface).
///
/// Holds one window plus one block, never a contig, let alone the genome (spec §2.6, §7).
/// The regions it yields tile the requested spans exactly: contiguous, non-overlapping,
/// complete, and maximal. Pure function of the reference, the spans, and the config;
/// `window_bp` changes memory and never the output.
///
/// # It owns its inputs
///
/// `reference` and `spans` are taken **by value** so the whole iterator can be moved onto
/// a producer thread (spec §7). That ownership is also what lets it evict: the walk slides
/// forward and never looks back, so it releases the bases it has passed
/// ([`EvictableRefSeq`]) — without which a windowed reference's buffer would grow to hold
/// the contig and the memory bound would be a claim rather than a fact.
///
/// # Generic over the reference, and why that is not the arch's "concrete"
///
/// The arch doc specified a concrete `WindowedRefSeq`, reasoning that the walk needs raw
/// bytes and eviction — "impl capabilities, not trait methods". They are trait methods
/// now ([`RawRefSeq`], [`ContigTable`], [`EvictableRefSeq`]), each added where a capability
/// had to be *required* rather than assumed. Being generic is what lets **the same walk**
/// run over an `InMemoryRefSeq` in tests, which is what makes window-invariance against
/// `partition_resident` testable at all — the strongest check step 3 has.
///
/// # Errors are in-stream and fatal
///
/// Every window reads the reference. `None` meaning both "done" and "a window failed"
/// would silently un-call the rest of the genome, so a read failure is yielded once as
/// `Some(Err(_))` and then the iterator is done — fused, and it says so
/// ([`std::iter::FusedIterator`]) rather than leaving a caller to find out (spec §8.2).
pub struct TypedRegionIterator<R> {
    reference: R,
    config: TypedRegionConfig,
    /// The **scan** spans still to walk, each carrying the requested spans inside it,
    /// **reversed** so the next one is a `pop` (spec §2.5).
    remaining: Vec<ScanSpan>,
    /// The span being walked, if any.
    current: Option<SpanWalk>,
    /// Regions resolved and waiting to be handed out — one block's worth at most.
    queue: std::collections::VecDeque<TypedRegion>,
    /// The window's bases, reused across windows (`typed_regions.md` §7).
    bases: Vec<u8>,
    counts: TypedRegionCounts,
    /// Latched by a fatal error or by exhaustion; the fused contract.
    done: bool,
}

/// One window's contribution to the running tally — the parts that cannot be read off the
/// regions the walk emits, because they are about repeats that never became one.
struct WindowTally {
    /// Cleaned repeat coverage this window owns, in bp.
    repeat_bp: u64,
    /// Repeat bp this window's core owns that classification turned down, by reason.
    rejected: RejectionCounts,
}

/// One span's walk in progress: where the windows are up to, and the block carries.
struct SpanWalk {
    contig: ContigId,
    /// The contig's name — from the reference's table, and the name every `SsrSegment` gets.
    chrom: String,
    /// The contig's **true** length, from the reference's table (spec §2.6's provenance
    /// obligation). Not the span's, and never a window's.
    contig_len: u64,
    cursor: WindowCursor,
    walk: BlockWalk,
    /// What the user asked for inside this scan span — the **emit** set (spec §2.5). The
    /// walk types the whole scan span; only what overlaps these comes out, and territory
    /// is clipped to them ([`emit_into`]).
    requested: Vec<GenomeRegion>,
    /// The last base **emitted**, for keeping the output non-overlapping when an object
    /// is emitted whole past a requested edge. Distinct from `BlockWalk::emitted_upto`,
    /// which tracks the walk's own generic carry across the *scan* span.
    emitted_upto: u64,
}

impl<R: RawRefSeq + ContigTable + EvictableRefSeq> TypedRegionIterator<R> {
    /// Walk `spans`. **Infallible in the arch doc, fallible here** — and the difference is
    /// a real one the arch could not have known: `GenomeRegions` validated the spans
    /// against a contig table, but nothing ties *that* table to the one this reference
    /// carries, so a span naming a contig the reference does not have is still reachable.
    /// Checking it up front, once, beats discovering it mid-walk (spec §8.2's whole
    /// point).
    pub fn over_regions(
        reference: R,
        spans: GenomeRegions,
        config: TypedRegionConfig,
    ) -> Result<Self, TypedRegionError> {
        Self::over_spans(reference, spans.iter().collect(), config)
    }

    /// Walk one contig, end to end — the whole-contig case as a span like any other.
    pub fn over_contig(
        reference: R,
        contig: ContigId,
        config: TypedRegionConfig,
    ) -> Result<Self, TypedRegionError> {
        let len = contig_entry(&reference, contig)?.length;
        // A zero-length contig has no 1-based position to cover, so it contributes no
        // span — `RegionSet`'s rule (C3), restated here for a direct caller.
        let spans = if len == 0 {
            Vec::new()
        } else {
            vec![GenomeRegion {
                contig,
                start: Position(1),
                end: Position(len),
            }]
        };
        Self::over_spans(reference, spans, config)
    }

    fn over_spans(
        reference: R,
        requested: Vec<GenomeRegion>,
        config: TypedRegionConfig,
    ) -> Result<Self, TypedRegionError> {
        // **A swept knob that would silently un-bundle** (spec §10 sweeps both of these).
        // The scan's margin is `max_repeat_len`, so a window sees every repeat within
        // `max_repeat_len` of its core. The bundle flank test needs to see every repeat
        // within `flank_bp` of a core repeat. If the margin were the narrower of the two,
        // a core tract's neighbour could fall outside the window, go unseen, and the tract
        // would be classified as a clean locus instead of bundled — no error, and a
        // different answer for every `window_bp`.
        //
        // **An error, not an `assert!`.** Both knobs are on the command line, so this is
        // reachable from a typo — and user input must not panic. The reason it was an
        // `assert!` still holds and is served better here: A2's rule is that a swept
        // knob's guard must survive `--release`, and a `Result` does, unconditionally.
        if config.max_repeat_len.get() < config.criteria.flank_bp {
            return Err(TypedRegionError::MarginNarrowerThanFlank {
                max_repeat_len: config.max_repeat_len.get(),
                flank_bp: config.criteria.flank_bp,
            });
        }
        // Fail now, not in the middle of chromosome 7 (above).
        for span in &requested {
            contig_entry(&reference, span.contig)?;
        }
        let mut scan = scan_set(&reference, &requested)?;
        // Reversed so the next span is a `pop`.
        scan.reverse();
        Ok(Self {
            reference,
            config,
            remaining: scan,
            current: None,
            queue: std::collections::VecDeque::new(),
            bases: Vec::new(),
            counts: TypedRegionCounts::default(),
            done: false,
        })
    }

    /// The running tally — readable mid-walk, complete once exhausted (spec §3.1).
    pub fn counts(&self) -> &TypedRegionCounts {
        &self.counts
    }

    /// Do one step of work: start a span, scan a window, or finish a span. `Ok(false)`
    /// when there is nothing left at all.
    fn advance(&mut self) -> Result<bool, TypedRegionError> {
        if self.current.is_none() {
            let Some(span) = self.remaining.pop() else {
                return Ok(false);
            };
            let entry = contig_entry(&self.reference, span.scan.contig)?;
            let contig_len = entry.length;
            let chrom = entry.name.clone();
            self.current = Some(SpanWalk {
                contig: span.scan.contig,
                chrom,
                contig_len,
                // Cores tile the **scan** span; margins come from the whole contig, which
                // is not the same number once the span is a BED region (`WindowCursor`).
                cursor: WindowCursor::new(
                    RegionSpan {
                        start: span.scan.start.get() - 1,
                        end: span.scan.end.get(),
                    },
                    contig_len,
                    &self.segment_options(),
                ),
                walk: BlockWalk::new(span.scan.contig, span.scan.end.get()),
                // One span counted per span the **user asked for**, not per span walked:
                // the scan set is this walk's business, and coalescing two requested
                // regions into one scan is not the user losing a region.
                requested: span.requested,
                emitted_upto: 0,
            });
            self.counts.spans += self.current.as_ref().expect("just set").requested.len() as u64;
            return Ok(true);
        }

        let state = self.current.as_mut().expect("just checked");
        let span_done = match state.cursor.next_window() {
            Some(plan) => {
                let tally = Self::scan_and_absorb(
                    &self.reference,
                    &mut self.bases,
                    state,
                    &self.config,
                    &plan,
                )?;
                // Every base of cleaned repeat coverage this window owns. Clipped to the
                // core, and cores tile, so no base is counted twice. The loci that cancel
                // it are subtracted as they are emitted (`tally`).
                self.counts.repeat_bp_with_no_locus += tally.repeat_bp;
                let r = &mut self.counts.rejected_by_reason;
                r.copy_floor += tally.rejected.copy_floor;
                r.purity += tally.rejected.purity;
                r.compound += tally.rejected.compound;
                r.no_clean_trim += tally.rejected.no_clean_trim;
                r.flank_clamped += tally.rejected.flank_clamped;
                // Release what the walk has passed — it never looks back. This is the
                // difference between "holds one window" as a claim and as a fact.
                self.reference.evict_before(plan.fetched.start + 1);
                false
            }
            None => {
                state.walk.finish(&self.config);
                true
            }
        };

        // Drain whatever the walk resolved through the emit policy (spec §2.5): the walk
        // types the whole **scan** span, and only what the user asked for comes out.
        // `&mut self.queue` and `self.current` are disjoint fields, which is why they are
        // borrowed by name rather than through a method.
        let queue = &mut self.queue;
        let state = self.current.as_mut().expect("still walking");
        while let Some(region) = state.walk.out.pop_front() {
            emit_into(queue, region, &state.requested, &mut state.emitted_upto);
        }
        if span_done {
            self.current = None;
        }
        Ok(true)
    }

    /// One window: detect (the cursor planned it), clean, classify, absorb.
    ///
    /// An associated function rather than a method, because it needs `&self.reference`
    /// and `&mut self.current` at once — disjoint fields, which the borrow checker will
    /// grant only when they are named separately.
    ///
    /// Returns this window's share of the tally ([`WindowTally`]).
    fn scan_and_absorb(
        reference: &R,
        bases: &mut Vec<u8>,
        state: &mut SpanWalk,
        config: &TypedRegionConfig,
        plan: &WindowPlan,
    ) -> Result<WindowTally, TypedRegionError> {
        let contig = state.contig;

        // 1. Detect. The scanned slice is the cursor's to decide; this reads exactly it.
        reference.fetch_raw_into(
            contig,
            plan.fetched.start + 1,
            plan.fetched.end - plan.fetched.start,
            bases,
        )?;
        let window = scan_window(plan, bases, config.periods, &config.scan);

        // 2. Clean. **Over the whole fetched slice, margins included**, which is what
        //    makes the verdict the resident one: redundancy elimination is a neighbour
        //    rule, and a core interval's neighbours live in the margin.
        let cleaned = segment_criteria::prefilter(&window.detections, &config.criteria);

        // 3. Admit, over **the slice we already scanned** — no second fetch, no second
        //    margin. Admission reads each tract's own bases (motif, purity) and answers
        //    the flank question by arithmetic against `contig_len`, so the scan slice is
        //    exactly what it needs. Until 2026-07-17 `SsrSegment` embedded tract ± flank_bp of
        //    bases, which a tract at the slice's own edge could not supply; that cost a
        //    wider re-fetch of every window, and the payload had no consumer.
        let bases_start = plan.fetched.start;
        let recs = cleaned
            .iter()
            .map(|iv| RepeatInterval {
                start: iv.start - bases_start,
                end: iv.end - bases_start,
                ..*iv
            })
            .collect();
        let classified = segment_criteria::classify(
            recs,
            &state.chrom,
            bases,
            Position(bases_start + 1),
            // The contig's length, from the reference's table. Never `bases.len()`, never
            // the window's end — spec §2.6's provenance obligation, and the mistake no
            // signature can catch.
            Bp(state.contig_len),
            &config.criteria,
        );

        // The tally's raw material, taken before `absorb` consumes the window.
        let mut tally = WindowTally {
            // This window's own share of the repeat coverage — clipped to the core, so
            // cores tile and nothing is double-counted.
            repeat_bp: window
                .coverage_in_core(&cleaned)
                .into_iter()
                .map(|(s, e)| e - s)
                .sum(),
            rejected: RejectionCounts::default(),
        };
        // **Rejections are attributed to a core, exactly like everything else.** `classify`
        // ran over the whole fetched slice, so every window that can see a repeat rejects
        // it again; counting them all would make the tally depend on `window_bp`, which
        // is a memory knob. The interval's own start decides whose it is.
        for (iv, reason) in &classified.rejected {
            let start = iv.start + bases_start;
            if start >= plan.core.start && start < plan.core.end {
                tally.rejected.add(*reason, iv.end - iv.start);
            }
        }

        // 4-5. Cap and partition, one block at a time.
        state
            .walk
            .absorb(&window, &cleaned, classified, bases_start, config);
        Ok(tally)
    }

    fn segment_options(&self) -> SegmentOptions {
        SegmentOptions {
            // For the scan this is the detection margin only — the satellite cap is the
            // walk's own business (`resolve_features`). Same number, because spec §2.6
            // says it must be: the margin exists to capture whole any repeat that is not
            // a satellite.
            max_repeat_len: self.config.max_repeat_len.get(),
            window_bp: self.config.window_bp.get(),
            // `RegionScanner`'s smoothing knobs; the window rules read neither. Zeroed
            // rather than defaulted so a reader does not go looking for where this walk
            // smooths. It does not.
            merge_gap: 0,
            min_repeat_len: 0,
        }
    }

    fn tally(&mut self, region: &TypedRegion) {
        let bp = region.region.len();
        match &region.kind {
            RegionKind::SsrSegment(_) => {
                self.counts.ssr_loci += 1;
                // This repeat coverage DID yield a locus, so it is not part of the gap.
                // The locus's bases are a subset of the coverage counted when its window
                // was scanned (a locus is a trimmed tract, and a tract is coverage), so
                // this cannot underflow.
                self.counts.repeat_bp_with_no_locus -= bp;
            }
            RegionKind::SsrBundle { .. } => {
                self.counts.ssr_bundles += 1;
                // **The number spec §10's bundle question needs and has never had** —
                // production drops these uncounted.
                self.counts.ssr_bundle_bp += bp;
            }
            RegionKind::Generic => self.counts.generic += 1,
            RegionKind::Satellite => {
                self.counts.satellites += 1;
                self.counts.satellite_bp += bp;
            }
        }
    }
}

impl<R: RawRefSeq + ContigTable + EvictableRefSeq> Iterator for TypedRegionIterator<R> {
    type Item = Result<TypedRegion, TypedRegionError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(region) = self.queue.pop_front() {
                self.tally(&region);
                return Some(Ok(region));
            }
            if self.done {
                return None;
            }
            match self.advance() {
                Ok(true) => continue,
                Ok(false) => {
                    self.done = true;
                    return None;
                }
                Err(e) => {
                    // Fatal: report once, then done. Continuing would scan a hole, and a
                    // hole in a partition is invisible to everything downstream.
                    self.done = true;
                    return Some(Err(e));
                }
            }
        }
    }
}

impl<R: RawRefSeq + ContigTable + EvictableRefSeq> std::iter::FusedIterator
    for TypedRegionIterator<R>
{
}

/// One span to **scan** — a whole contig — and the spans inside it the user actually
/// **asked for** (spec §2.5).
///
/// Two sets, because a BED edge must not cut a decision in half: a repeat at 999 is
/// bundled away by a neighbour at 1030, so a walk that never looks past 1000 classifies a
/// locus the whole-genome run rejects — same reference, different calls, because of
/// `--regions`. [`scan_set`] says why the scan set is the whole contig and not something
/// cheaper.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ScanSpan {
    /// What the walk tiles with cores: **the whole contig** ([`scan_set`]).
    scan: GenomeRegion,
    /// The requested spans inside it, in order — what the walk is allowed to emit.
    /// Non-empty by construction.
    requested: Vec<GenomeRegion>,
}

/// **The scan set is the whole of every contig the request touches** (spec §2.5, owner
/// 2026-07-17).
///
/// # Why the whole contig, and not the requested span grown by a margin
///
/// E2 grew each span by `max_repeat_len` and coalesced, which spec §2.5 chose over
/// whole-contig scans on cost: *"makes a 10 kb region pay for a 90 Mb chromosome"*. **A
/// fixed margin cannot deliver the emit rule.** Every finding that intersects a requested
/// edge is returned **whole** ([`clips_at_a_bed_edge`]) — and "whole" is a fact about the
/// feature, which does not care what the walk scanned:
///
/// - a **satellite** is *by definition* longer than `max_repeat_len`, and the margin **is**
///   `max_repeat_len` — so an array running past it came back cut. Measured: a 3 kb array
///   reported `1001–4001` whole-genome and `1001–2300` under a BED, silently.
/// - a **bundle** chains with no reach at all (spec §2.6: *"A–B–C–D each 30 bp apart runs
///   past any margin you choose"*), so a fixed margin was never going to hold one.
///
/// §2.6 had already settled the principle for *window* edges — *"the data tells you when
/// it's over; you don't have to guess a reach"* — and a grown span is a guessed reach.
/// §2.5 was the odd one out.
///
/// **The cost is time, not memory.** A 10 kb BED on a 90 Mb chromosome now scans 90 Mb —
/// but through the same window, so peak memory is unchanged (spec §7), and the emit set is
/// as narrow as ever. What it buys is that a BED costs the user nothing in *correctness*:
/// every finding is exactly the one a whole-genome run reports, which
/// `a_bed_returns_the_same_findings_the_whole_genome_run_does` asserts as object identity
/// rather than as a resemblance.
///
/// It also **deletes** a rule rather than adding one: with no scan edges but the contig's,
/// a cluster can no longer be cut, so `resolve_features` needs no edge case and
/// `bundle_clusters` gets its "a bundle has ≥ 2 members" assert back at full strength.
fn scan_set<R: ContigTable>(
    reference: &R,
    requested: &[GenomeRegion],
) -> Result<Vec<ScanSpan>, TypedRegionError> {
    let mut out: Vec<ScanSpan> = Vec::new();
    for &span in requested {
        match out.last_mut() {
            // `GenomeRegions` hands them over in genomic order (C3), so a contig's spans
            // are adjacent and one pass groups them.
            Some(last) if last.scan.contig == span.contig => last.requested.push(span),
            _ => {
                let contig_len = contig_entry(reference, span.contig)?.length;
                out.push(ScanSpan {
                    scan: GenomeRegion {
                        contig: span.contig,
                        start: Position(1),
                        end: Position(contig_len),
                    },
                    requested: vec![span],
                });
            }
        }
    }
    Ok(out)
}

/// Whether a typed region is **clipped** to the user's BED edge — and only `Generic` is
/// (spec §2.5).
///
/// # Every *finding* comes back whole; only the default clips
///
/// **Owner's rule, 2026-07-17:** a microsatellite, a bundle *or a satellite* that
/// intersects a BED edge is returned whole. `Generic` alone is clipped, because it is the
/// only kind that is not a finding — *"nothing more specific can be said here"* stays true
/// of any stretch of a generic run, so a clipped one still says exactly what it said.
///
/// Each of the other three is a **claim about its own extent**, and half of it is a
/// different claim:
///
/// - `SsrSegment` — a genetic object. Half a locus is not a locus: its coordinates, its
///   motif and its copy number all describe the whole tract.
/// - `SsrBundle` — carries its member tracts; clip the hull and the members describe bases
///   outside their own region.
/// - `Satellite` — **the case E2 first got wrong.** The label means *"an array too long to
///   be a microsatellite"*, so the extent is the claim and `max_repeat_len` is what makes
///   it. Clip the fixture's 1.2 kb array to a 300 bp request and the result is a
///   `Satellite` region of 300 bp — a span that contradicts the very cap that produced the
///   label. E2 argued the opposite from `RegionKind`'s shape (`Satellite` carries no
///   payload, so nothing could be left misdescribed); that read the *type* correctly and
///   the *meaning* wrongly. What a region carries is not what it claims.
///
/// Spec §2.5 named only loci and bundles and was silent on satellites; it now states this
/// rule.
fn clips_at_a_bed_edge(kind: &RegionKind) -> bool {
    match kind {
        RegionKind::Generic => true,
        RegionKind::SsrSegment(_) | RegionKind::SsrBundle { .. } | RegionKind::Satellite => false,
    }
}

/// Narrow one walked region to what the user asked for, and queue what survives (spec
/// §2.5).
///
/// - **Outside every requested span** → dropped. It was scanned so that the regions inside
///   would be *right*, not to be shown.
/// - **A finding** — locus, bundle or satellite ([`clips_at_a_bed_edge`]) → emitted
///   **whole**, even past the edge: the requested span grows to hold it, and that grown
///   span — the *effective* region — is what the partition invariant is stated over.
/// - **`Generic`** → **clipped** to each requested span it overlaps, which may be more
///   than one: every span on a contig shares that contig's scan, so a generic run across
///   two of them must come back as two regions with the gap dropped, not one region
///   covering ground the user did not ask for.
///
/// `emitted_upto` keeps the output non-overlapping when a finding has just been emitted
/// whole past an edge and the next requested span starts inside it.
fn emit_into(
    queue: &mut std::collections::VecDeque<TypedRegion>,
    region: TypedRegion,
    requested: &[GenomeRegion],
    emitted_upto: &mut u64,
) {
    let overlaps = |r: &GenomeRegion| {
        r.contig == region.region.contig
            && r.start.get() <= region.region.end.get()
            && region.region.start.get() <= r.end.get()
    };

    if !clips_at_a_bed_edge(&region.kind) {
        if requested.iter().any(overlaps) {
            *emitted_upto = region.region.end.get();
            queue.push_back(region);
        }
        return;
    }

    for r in requested.iter().filter(|r| overlaps(r)) {
        let start = region
            .region
            .start
            .get()
            .max(r.start.get())
            .max(*emitted_upto + 1);
        let end = region.region.end.get().min(r.end.get());
        if start > end {
            continue;
        }
        *emitted_upto = end;
        queue.push_back(TypedRegion {
            region: GenomeRegion {
                contig: region.region.contig,
                start: Position(start),
                end: Position(end),
            },
            kind: region.kind.clone(),
        });
    }
}

/// The reference's own table entry for `contig` — **the provenance seam** (spec §2.6).
/// Every contig length and name the walk uses comes through here.
fn contig_entry<R: ContigTable>(
    reference: &R,
    contig: ContigId,
) -> Result<&crate::fasta::ContigEntry, TypedRegionError> {
    reference
        .contigs()
        .entries
        .get(contig.get() as usize)
        .ok_or(TypedRegionError::Reference(RefSeqError::UnknownContig(
            contig,
        )))
}

/// The windowed walk's state: the open block, and the three carries.
///
/// Everything here is bounded by *one block* plus three coordinates — never by the
/// contig. That is the claim `partition_windowed` exists to make good on, and the
/// reason the pieces are named as carries rather than accumulators.
struct BlockWalk {
    contig: ContigId,
    contig_len: u64,
    /// The open block's coverage runs, ascending; the last one is still growing.
    /// **Carry:** a satellite is longer than a window by definition, so a run's verdict
    /// is not reachable inside the window that started it.
    runs: Vec<CoverageRun>,
    /// The open block's classified loci, awaiting their runs' verdict (a locus inside a
    /// satellite is swallowed, spec §2.1).
    loci: Vec<SsrSegment>,
    /// The open block's bundle members, contig coordinates. **Carry:** a cluster chains
    /// across a window boundary as happily as inside one.
    bundled: Vec<RepeatInterval>,
    /// **Carry:** the open generic run — everything up to here is emitted, so the gap
    /// from here to the next feature is generic. One `u64` for a megabase.
    emitted_upto: u64,
    /// Resolved regions waiting to be read out, in genomic order. A queue rather than a
    /// list because `TypedRegionIterator` drains it as it goes: it holds one block's
    /// worth, never a contig's (spec §7).
    out: std::collections::VecDeque<TypedRegion>,
}

/// One thing a window hands the walk, at a 1-based start. Sorted into a single
/// ascending stream per window, because the block boundary is a question about
/// *position* and answering it per kind would ask it three times.
enum WindowItem {
    Coverage(CoverageRun),
    Segment(SsrSegment),
    Bundled(RepeatInterval),
}

impl WindowItem {
    fn start(&self) -> u64 {
        match self {
            Self::Coverage(run) => run.start,
            Self::Segment(l) => l.start(),
            // 0-based half-open → 1-based inclusive (spec §4).
            Self::Bundled(iv) => iv.start + 1,
        }
    }

    /// Coverage sorts first at a tie: a feature's block must be open before the feature
    /// arrives, and it is coverage that opens it.
    fn rank(&self) -> u8 {
        match self {
            Self::Coverage(_) => 0,
            _ => 1,
        }
    }
}

impl BlockWalk {
    fn new(contig: ContigId, contig_len: u64) -> Self {
        Self {
            contig,
            contig_len,
            runs: Vec::new(),
            loci: Vec::new(),
            bundled: Vec::new(),
            emitted_upto: 0,
            out: std::collections::VecDeque::new(),
        }
    }

    /// Fold one window's worth of coverage, loci and bundle members into the open
    /// block, closing it whenever an item proves it over.
    ///
    /// The three attributions — and each is the answer to a different question:
    ///
    /// - **coverage** is *clipped to the core*, so the cores' contributions abut and a
    ///   satellite rejoins across every window it spans;
    /// - **bundle members** are attributed by *the interval's* start, so each is
    ///   contributed once; and
    /// - **loci** by *the locus's* start (post-trim), which is also exactly once —
    ///   every window whose slice holds the tract classifies the same locus (that is what
    ///   the margin buys), and exactly one core holds its start.
    fn absorb(
        &mut self,
        window: &ScannedWindow,
        cleaned: &[RepeatInterval],
        classified: segment_criteria::Classified,
        bases_start: u64,
        config: &TypedRegionConfig,
    ) {
        let core = window.core;
        let in_core = |start_0based: u64| start_0based >= core.start && start_0based < core.end;

        let mut items: Vec<WindowItem> = Vec::new();
        // Coverage of the **cleaned** intervals, not the raw ones (spec §2.4): capping
        // detector noise would let it declare a satellite and swallow the real loci
        // underneath. This is why the window rules are methods over a caller-chosen
        // set — `ScannedWindow::coverage_in_core`.
        items.extend(
            window
                .coverage_in_core(cleaned)
                .into_iter()
                // 0-based half-open → 1-based inclusive.
                .map(|(s, e)| {
                    WindowItem::Coverage(CoverageRun {
                        start: s + 1,
                        end: e,
                    })
                }),
        );
        items.extend(
            classified
                .loci
                .into_iter()
                .filter(|l| in_core(l.start() - 1))
                .map(WindowItem::Segment),
        );
        items.extend(
            classified
                .bundled
                .into_iter()
                // Back to contig coordinates: `classify` hands these back in the offsets
                // it was given, which were offsets into the bases slice.
                .map(|iv| RepeatInterval {
                    start: iv.start + bases_start,
                    end: iv.end + bases_start,
                    ..iv
                })
                .filter(|iv| in_core(iv.start))
                .map(WindowItem::Bundled),
        );
        items.sort_by_key(|i| (i.start(), i.rank()));

        for item in items {
            let start = item.start();
            if self.block_is_open() && start > self.block_barrier(config) {
                self.close_block(config);
            }
            match item {
                WindowItem::Coverage(run) => match self.runs.last_mut() {
                    // Abutting runs merge as well as overlapping ones — the same rule
                    // as `coverage_runs`, and the reason cores tiling is enough for a
                    // satellite to rejoin across windows.
                    Some(last) if run.start <= last.end + 1 => last.end = last.end.max(run.end),
                    _ => self.runs.push(run),
                },
                // A feature's start lies inside a cleaned interval, hence inside a run,
                // hence inside the open block — coverage sorts first, so that run has
                // already arrived. If this ever fires, the attribution is wrong and the
                // partition is about to be built out of order; loudly, in release, is
                // the only useful way to learn that.
                WindowItem::Segment(l) => {
                    assert!(
                        self.block_is_open(),
                        "a locus at {} arrived with no coverage under it",
                        l.start()
                    );
                    self.loci.push(l);
                }
                WindowItem::Bundled(iv) => {
                    assert!(
                        self.block_is_open(),
                        "a bundle member at {} arrived with no coverage under it",
                        iv.start + 1
                    );
                    self.bundled.push(iv);
                }
            }
        }
    }

    fn block_is_open(&self) -> bool {
        !self.runs.is_empty()
    }

    /// The last position that still belongs to the open block: one flank past its
    /// rightmost base. Wider than every radius the partition's rules have — see
    /// [`resolve_features`].
    fn block_barrier(&self, config: &TypedRegionConfig) -> u64 {
        self.runs
            .last()
            .expect("the block is open")
            .end
            .saturating_add(config.criteria.flank_bp)
    }

    /// Resolve the open block and emit it, with the generic run that preceded it.
    fn close_block(&mut self, config: &TypedRegionConfig) {
        let loci = std::mem::take(&mut self.loci);
        let bundled = std::mem::take(&mut self.bundled);
        let features = resolve_features(&self.runs, loci, &bundled, self.contig, config);
        self.runs.clear();

        for f in features {
            let start = f.region.start.get();
            // The open generic run ends here — the third carry, spent. One `Generic`
            // however long it is: maximality is a correctness requirement, not tidiness
            // (spec §2.3, `fill_generic_gaps`).
            if start > self.emitted_upto + 1 {
                self.out
                    .push_back(self.generic(self.emitted_upto + 1, start - 1));
            }
            // Non-overlap, checked rather than assumed. `partition_resident` leaves this
            // as `fill_generic_gaps`'s unchecked precondition; here the features arrive
            // from three carries across many windows, so it is worth one comparison per
            // feature to learn about a broken partition from the walk rather than from a
            // downstream consumer.
            assert!(
                start > self.emitted_upto,
                "features overlap: a region starting at {start} follows one ending at {}",
                self.emitted_upto
            );
            self.emitted_upto = f.region.end.get();
            self.out.push_back(f);
        }
    }

    fn generic(&self, start: u64, end: u64) -> TypedRegion {
        TypedRegion {
            region: GenomeRegion {
                contig: self.contig,
                start: Position(start),
                end: Position(end),
            },
            kind: RegionKind::Generic,
        }
    }

    /// Close the last block and run the generic carry out to the contig's end, leaving
    /// everything in [`Self::out`] for the caller to drain.
    ///
    /// `&mut self` rather than consuming: `TypedRegionIterator` reaches this while the
    /// walk is a field it is holding, and it drains the queue afterwards like any other
    /// step. Idempotent — a second call has no open block and nothing left to reach.
    fn finish(&mut self, config: &TypedRegionConfig) {
        if self.block_is_open() {
            self.close_block(config);
        }
        if self.emitted_upto < self.contig_len {
            self.out
                .push_back(self.generic(self.emitted_upto + 1, self.contig_len));
            self.emitted_upto = self.contig_len;
        }
    }
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
    /// The config's detection margin is narrower than its bundle radius, which
    /// would make the answer depend on `window_bp` (see [`TypedRegionConfig`]).
    ///
    /// Raised by the constructor, before any work: both knobs are on the command
    /// line, so this is what a typo looks like, and a typo must not panic.
    #[error(
        "max_repeat_len ({max_repeat_len}) is the window's detection margin and must not \
         be narrower than flank_bp ({flank_bp}), the bundle radius: a window that cannot \
         see a core tract's neighbours would classify them as loci instead of bundling them"
    )]
    MarginNarrowerThanFlank { max_repeat_len: u64, flank_bp: u64 },
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::ref_seq::{InMemoryRefSeq, RefSeq};

    /// `Default` is the catalog's settings, for spec §8's comparability — so the
    /// oracle compares like with like. Pinned against the **literals**, not
    /// against the consts that define it: asserting `periods.min() ==
    /// DEFAULT_SCAN_MIN_PERIOD` would compare a constant with itself and could not
    /// fail (the same tautology the A2 review caught in `SsrSegmentCriteria`).
    #[test]
    fn default_config_is_the_catalogs_settings() {
        let c = TypedRegionConfig::default();
        assert_eq!(c.periods.min(), 1, "the scanner scans period 1..");
        assert_eq!(c.periods.max(), 6, "..to 6");
        assert_eq!(c.max_repeat_len, Bp(1000), "the 1 kb satellite cap");
        assert_eq!(c.window_bp, Bp(100_000), "100 kb window");
        assert_eq!(c.scan, ScanParams::default());
        assert_eq!(c.criteria, SsrSegmentCriteria::default());
    }

    /// **The scanner scans wider than classification classifies, and that is deliberate.**
    ///
    /// Two `periods` knobs look redundant until you ask why: a homopolymer tiles
    /// under every motif length, so the scanner emits its span at period 2, 3, 4, 5
    /// and 6 as well as 1 — and **only the period-1 interval divides them all**, so
    /// only it can eliminate them in `prefilter`. Scan from 2 and those aliases have
    /// no eliminator: a poly-A enters the partition as three "repeats" (periods 2, 3
    /// and 5). **An eliminator that was never detected cannot eliminate** — which is
    /// what "not seeing them is not the same as dropping them" actually means.
    ///
    /// If these two ever coincide at the floor, the aliasing bug of 2026-07-17 is
    /// back. `a_homopolymer_does_not_survive_as_a_period_two_repeat` is the
    /// behavioural statement; this is the config-level guard.
    #[test]
    fn the_scan_range_is_wider_than_the_criteria_range() {
        let c = TypedRegionConfig::default();
        assert!(
            c.periods.min() < c.criteria.periods.min(),
            "the scanner must SEE period {} so the pre-filter has an eliminator for \
             that run's aliases; classification starts at {}",
            c.periods.min(),
            c.criteria.periods.min()
        );
        assert_eq!(
            c.periods.max(),
            c.criteria.periods.max(),
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
        assert_eq!(c.repeat_bp_with_no_locus, 0);
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

    /// A contig with one clean isolated (AT)*8 tract: Generic / SsrSegment / Generic.
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
                RegionKind::SsrSegment(_) => "locus",
                RegionKind::SsrBundle { .. } => "bundle",
                RegionKind::Generic => "generic",
                RegionKind::Satellite => "satellite",
            })
            .collect();
        assert_eq!(kinds, vec!["generic", "locus", "generic"]);

        // And the locus is where the tract is — 1-based, so 0-based + 1.
        let RegionKind::SsrSegment(l) = &regions[1].kind else {
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
    /// because nothing was ever classified, and the assertion passed for entirely the
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
                .any(|r| matches!(r.kind, RegionKind::SsrSegment(_))),
            "the locus is swallowed by the satellite — the §2.4 cost, made visible"
        );

        // **The control.** Raise the cap above the array and the SAME bases classify a
        // locus — so the absence above is the cap doing its job, not classification
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
                .filter(|r| matches!(r.kind, RegionKind::SsrSegment(_)))
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

    /// A contig with a **3 kb array** — three times the 1 kb margin, so a walk that scanned
    /// only a grown BED span could not see the whole of it.
    fn contig_with_a_long_array() -> Vec<u8> {
        let mut bases = filler(8000);
        for i in 0..1500 {
            bases[1000 + i * 2..1002 + i * 2].copy_from_slice(b"AT");
        }
        bases
    }

    /// The kinds, in order — for asserting a partition's shape without spelling out four
    /// `matches!` arms at every site.
    fn kinds_of(regions: &[TypedRegion]) -> Vec<&'static str> {
        regions
            .iter()
            .map(|r| match &r.kind {
                RegionKind::SsrSegment(_) => "locus",
                RegionKind::SsrBundle { .. } => "bundle",
                RegionKind::Generic => "generic",
                RegionKind::Satellite => "satellite",
            })
            .collect()
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

    /// **An impure tract is `Generic`** (spec §8's fixture list, §2.2) — with the control
    /// that shows it is the interruption's doing.
    ///
    /// # But "impure → Generic" is not one rule, and it took this test to see it
    ///
    /// The scanner decides an impure tract's fate long before classification does, because
    /// Ruzzo–Tompa returns **maximal-scoring** segments. So an interruption has three
    /// possible outcomes, and only the third is this fixture:
    ///
    /// 1. **small** — the surrounding matches pay for it, the tract stays whole, and it is
    ///    pure enough to classify: a **locus**;
    /// 2. **large, with long pieces** — the segment splits, and the pure pieces are close
    ///    together: a **bundle** (or two loci, if far apart);
    /// 3. **large, with short pieces** — the segment splits and each piece falls under the
    ///    copy floor: **`Generic`**, which is this fixture (two 8 bp halves = 4 copies,
    ///    under the floor of 5).
    ///
    /// **Admission's purity gate is not what does any of this.** It is unreachable from the
    /// walk: a tract impure enough to fail the 0.80 floor always contains a purer
    /// sub-segment that scores higher, so Ruzzo–Tompa emits *that* instead — measured, not
    /// argued (a 0.79-purity fixture comes back a **locus**, its pure core). See
    /// `segment_criteria::the_walk_reaches_only_one_of_classifications_five_gates`.
    #[test]
    fn an_impure_tract_is_generic_when_its_pieces_fall_under_the_copy_floor() {
        let mut bases = filler(240);
        bases[100..126].copy_from_slice(b"ATATATATGGGGGGGGGGATATATAT");
        let len = bases.len() as u64;

        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "impure tract");
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::SsrSegment(_))),
            "an impure tract is not a locus: {:?}",
            kinds_of(&regions)
        );
        // Not a hole either — spec §2.2's property, which is what the fixture is for.
        assert_eq!(regions.len(), 1, "the whole contig is one Generic region");

        // **The control.** Remove the interruption and the same span at the same settings
        // IS a locus — so `Generic` above is the interruption's doing, not a tract that was
        // never admissible (D1's trap).
        let mut bases = filler(240);
        bases[100..126].copy_from_slice(b"ATATATATATATATATATATATATAT");
        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_eq!(
            kinds_of(&regions),
            vec!["generic", "locus", "generic"],
            "pure, the same span is a locus"
        );
    }

    /// **A homopolymer is `Generic` at periods 2..=6** (spec §8's fixture list): nothing
    /// classifies period 1, and the bases are still covered.
    ///
    /// `prefilter` is what removes it, and **the whole homopolymer, not just its period-1
    /// label** — the 2026-07-17 ordering fix. A poly-A tiles under `AA`, `AAA`, `AAAAA`, so
    /// the scanner emits the same span at every period in scope; only the period-1 interval
    /// divides them all, so it is kept as an eliminator through redundancy elimination and
    /// dropped by the period floor afterwards. Before that fix, periods 2, 3 and 5 survived
    /// and one homopolymer entered the partition as three "repeats"
    /// (`a_homopolymer_does_not_survive_as_a_period_two_repeat` pins it at the unit level).
    ///
    /// The neighbouring tract is the other half: the copy floor drops period-1 specks
    /// before they can eliminate anything, so a real STR is never taken out by one — the
    /// poly-A cascade.
    #[test]
    fn a_homopolymer_is_generic_and_does_not_take_its_neighbour_with_it() {
        let mut bases = filler(300);
        // 20 bp of poly-A, then a clean (CAG) tract 60 bp away — far enough not to bundle.
        bases[100..120].copy_from_slice(&[b'A'; 20]);
        bases[180..204].copy_from_slice(b"CAGCAGCAGCAGCAGCAGCAGCAG");
        let len = bases.len() as u64;

        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "homopolymer");

        let loci: Vec<_> = regions
            .iter()
            .filter_map(|r| match &r.kind {
                RegionKind::SsrSegment(l) => Some(l),
                _ => None,
            })
            .collect();
        assert_eq!(
            loci.len(),
            1,
            "exactly one locus: the (CAG) tract. The homopolymer is not one, and — the \
             point — it did not eliminate the tract either: {:#?}",
            kinds_of(&regions)
        );
        assert_eq!(loci[0].motif().as_bytes(), b"CAG");
        assert!(
            loci[0].start() >= 181 && loci[0].start() <= 185,
            "and it is the tract at ~181, not something the cascade left behind: {}",
            loci[0].start()
        );
        // The homopolymer's own bases are covered, not a hole (spec §2.2).
        assert!(
            regions
                .iter()
                .any(|r| r.region.contains(Position(110)) && matches!(r.kind, RegionKind::Generic)),
            "the poly-A run is Generic territory"
        );

        // **And it reaches classification not at all** — it is gone by the pre-filter, so it is
        // turned down by no gate and counted under no reason. That is what "out of scope"
        // has to mean: before the 2026-07-17 fix this fixture recorded THREE `Compound`
        // rejections for one homopolymer, because its period-2/3/5 aliases each reached
        // `classify` separately. A rejection count that moves with the number of divisors of
        // a run's length is not measuring anything.
        let counts = counts_over(&bases, 100_000);
        assert_eq!(
            counts.rejected_by_reason,
            RejectionCounts::default(),
            "the homopolymer is filtered out, not rejected — and the (CAG) tract is \
             classified, so nothing here is turned down at all: {:?}",
            counts.rejected_by_reason
        );
    }

    /// **A rejected repeat is generic territory, not a hole** (spec §2.2). A
    /// low-copy tract classification turns down must still be *covered* — completeness
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
                .any(|r| matches!(r.kind, RegionKind::SsrSegment(_))),
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
                .any(|r| matches!(r.kind, RegionKind::SsrSegment(_))),
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
                .any(|r| matches!(r.kind, RegionKind::SsrSegment(_))),
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

    /// **Bundles do not spread.** A tract far enough from a cluster is classified
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
                RegionKind::SsrSegment(l) => Some(l),
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
        let clusters = segment_criteria::bundle_clusters(&bundled, 50);
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
    /// to the *cleaned* coverage, after classification, so the difference can only go
    /// one way. Capping raw coverage would make it bidirectional and untestable.
    ///
    /// **The oracle is the committed trf-mod-built golden catalog** — a different
    /// detector, a different code path, nothing ng touched. Its detector difference
    /// is already characterised (`scanner_parity`: 16/16, 15 exact, one ±1–2 bp
    /// boundary/phase wobble, one genuine scanner-only locus trf-mod's significance
    /// model rejected), so it is a yardstick rather than a confound — hence overlap
    /// matching, inherited from `scanner_parity` for exactly that reason.
    ///
    /// What this proves: the plumbing — the scan, the pre-filter, the classification
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
            criteria: SsrSegmentCriteria {
                min_purity: cat_params.min_purity,
                min_score: cat_params.min_score,
                flank_bp: u64::from(cat_params.flank_bp),
                ..SsrSegmentCriteria::default()
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
                    RegionKind::SsrSegment(l) => ours.push((name.clone(), l.start(), l.end())),
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

    // ---- D3: the windowed walk -------------------------------------------

    /// A contig with real structure, laid out so that **the default 1 kb window's
    /// core edges cut through features** — which is the only way a window-invariance
    /// test proves anything:
    ///
    /// | offset | what | why it is there |
    /// |---|---|---|
    /// | 990 | a lone `(AT)*10` tract | straddles the 1000 core edge: the locus itself crosses a boundary |
    /// | 1990, 2020 | two tracts 10 bp apart | a **bundle cluster** straddling the 2000 edge — its members land in different cores |
    /// | 3500 | a lone tract | an interior control: the same feature nowhere near an edge |
    /// | 4000..5200 | a 1.2 kb `(AT)` array | a **satellite** spanning two whole cores, plus the locus it swallows |
    ///
    /// The rest is aperiodic filler, which leaves generic runs several windows long
    /// (the maximality case). 6 kb at `window_bp = 1000` means a fetched slice of ~3 kb
    /// against a 6 kb contig — genuinely windowed, which a small fixture would not be:
    /// with the 1 kb margin, a contig under ~2 kb is fetched whole every time and every
    /// window-invariance assertion passes for free.
    fn windowing_fixture() -> Vec<u8> {
        let mut bases = filler(6000);
        let mut tract_at = |off: usize| {
            bases[off..off + 20].copy_from_slice(b"ATATATATATATATATATAT");
        };
        tract_at(990);
        tract_at(1990);
        tract_at(2020);
        tract_at(3500);
        for i in 0..600 {
            bases[4000 + i * 2..4002 + i * 2].copy_from_slice(b"AT");
        }
        bases
    }

    fn reference_over(chrom: &str, bases: &[u8]) -> InMemoryRefSeq {
        InMemoryRefSeq::from_named_contigs(vec![(chrom.to_string(), bases.to_vec())])
    }

    /// A **chain**: 20 bp `(AT)` tracts every 50 bp — a 30 bp gap between each, under the
    /// 50 bp flank — running from 1000 for `span_bp`. Every tract is close to the next, so
    /// the whole thing is one cluster however long it is.
    fn chained_cluster(span_bp: usize, total: usize) -> Vec<u8> {
        let mut bases = filler(total);
        let mut at = 1000;
        while at + 20 <= 1000 + span_bp {
            bases[at..at + 20].copy_from_slice(b"ATATATATATATATATATAT");
            at += 50;
        }
        bases
    }

    /// **A BED returns the same findings a whole-genome run does — the same *objects*, not
    /// merely the same kinds** (spec §2.5, owner 2026-07-17). This is what the whole-contig
    /// scan set buys, and it is the strongest form of BED-invariance available.
    ///
    /// # The two fixtures that broke the old design
    ///
    /// E2 grew each requested span by `max_repeat_len` and scanned that. Both of these
    /// defeat a margin, and both are pinned here because they are exactly the shapes a
    /// margin cannot hold:
    ///
    /// - **a 3 kb array** — a satellite is *by definition* longer than `max_repeat_len`,
    ///   and the margin **is** `max_repeat_len`. Under the old design a BED reported it as
    ///   `1001–2300` where the truth is `1001–4001`: a `Satellite` cut to 1300 bp,
    ///   silently. This is the common case, not an exotic one.
    /// - **a 1.5 kb chain** of tracts 30 bp apart — clustering has no reach at all (spec
    ///   §2.6: *"A–B–C–D each 30 bp apart runs past any margin you choose"*). §2.5 hoped
    ///   such a chain would be *"dense-repeat territory heading for `Satellite`"*; it is
    ///   not — a satellite is over-cap **coverage**, and coverage runs merge only where
    ///   they **abut**, so every run in this chain is 20 bp. It stays a bundle, and the old
    ///   design cut it to 24 members of 30.
    #[test]
    fn a_bed_returns_the_same_findings_the_whole_genome_run_does() {
        let config = TypedRegionConfig::default();

        // --- the 3 kb array: a satellite three times the old margin ---
        let bases = contig_with_a_long_array();
        let whole = partition_resident("chr1", ContigId(0), &bases, &config);
        let truth = whole
            .iter()
            .find(|r| matches!(r.kind, RegionKind::Satellite))
            .expect("the fixture has one array");
        assert!(
            truth.region.len() > 3000,
            "the array must be far longer than max_repeat_len: {:?}",
            truth.region
        );

        // Ask about 200 bp of its left end. The satellite comes back WHOLE — all 3 kb,
        // past both requested edges, identical to the whole-genome region.
        let got = walk_bed(&bases, &[(1100, 1300)], 333);
        let sat = got
            .iter()
            .find(|r| matches!(r.kind, RegionKind::Satellite))
            .expect("still a satellite");
        assert_eq!(
            sat.region, truth.region,
            "a BED must not shorten an array: the extent IS the claim, and a margin cannot \
             hold a feature defined as being longer than it"
        );

        // --- the 1.5 kb chain: a bundle no margin could hold ---
        let bases = chained_cluster(1500, 6000);
        let whole = partition_resident("chr1", ContigId(0), &bases, &config);
        let (truth_hull, truth_tracts) = whole
            .iter()
            .find_map(|r| match &r.kind {
                RegionKind::SsrBundle { tracts } => Some((r.region, tracts.clone())),
                _ => None,
            })
            .expect("the chain is one bundle");
        assert!(
            truth_hull.len() > 1400 && truth_tracts.len() >= 25,
            "the chain must run well past max_repeat_len: {truth_hull:?}"
        );
        assert!(
            !whole
                .iter()
                .any(|r| matches!(r.kind, RegionKind::Satellite)),
            "and it is NOT a satellite — §2.5 hoped it would be, but coverage runs merge \
             only where they abut and every run here is 20 bp: {whole:#?}"
        );

        let got = walk_bed(&bases, &[(900, 1200)], 333);
        let (hull, tracts) = got
            .iter()
            .find_map(|r| match &r.kind {
                RegionKind::SsrBundle { tracts } => Some((r.region, tracts.clone())),
                _ => None,
            })
            .expect("still a bundle");
        assert_eq!(hull, truth_hull, "the same hull");
        assert_eq!(
            tracts, truth_tracts,
            "and the same member tracts — all 30, not the 24 inside a grown span"
        );
    }

    /// A `(CAG)*8` microsatellite `gap` bp from a 1.3 kb `(AT)` array, on either side of
    /// it. Different motifs, so the detector cannot join the two into one tract.
    fn micro_near_satellite(micro_left: bool, gap: usize) -> Vec<u8> {
        let mut bases = filler(6000);
        let (micro_at, array_at) = if micro_left {
            (1000, 1000 + 24 + gap)
        } else {
            (1000 + 1300 + gap, 1000)
        };
        // (CAG)*8 — a different motif from the array's, so the detector cannot join them.
        for i in 0..8 {
            bases[micro_at + i * 3..micro_at + i * 3 + 3].copy_from_slice(b"CAG");
        }
        for i in 0..650 {
            bases[array_at + i * 2..array_at + i * 2 + 2].copy_from_slice(b"AT");
        }
        bases
    }

    /// **A microsatellite too close to a satellite is swallowed by it, from either
    /// side** (spec §2.4a — the owner's rule, 2026-07-16).
    ///
    /// The situation is physical and symmetric: a `(CAG)*8` tract 20 bp from a 1.3 kb
    /// array cannot be genotyped, because the flank on that side *is* array. So the
    /// answer must not depend on which side it sits — and before this rule it did, in
    /// two different wrong ways (`resolve_features`' docs).
    ///
    /// The 200 bp arm is the **control**, and it is what makes the 20 bp arm mean
    /// something: at 200 bp the same two features, built by the same code, give a clean
    /// `Generic / locus / Generic / Satellite` — so absorption at 20 bp is the *rule*
    /// firing, not classification quietly rejecting a tract that was never viable.
    #[test]
    fn a_microsatellite_beside_a_satellite_is_absorbed_into_it() {
        let config = TypedRegionConfig::default();
        let kinds = |regions: &[TypedRegion]| -> Vec<&'static str> {
            regions
                .iter()
                .map(|r| match &r.kind {
                    RegionKind::SsrSegment(_) => "locus",
                    RegionKind::SsrBundle { .. } => "bundle",
                    RegionKind::Generic => "generic",
                    RegionKind::Satellite => "satellite",
                })
                .collect()
        };

        for micro_left in [true, false] {
            let side = if micro_left { "left" } else { "right" };

            // --- absorbed: 20 bp of gap, less than the 50 bp flank ---
            let bases = micro_near_satellite(micro_left, 20);
            let regions = partition_resident("chr1", ContigId(0), &bases, &config);
            assert_partitions(
                &regions,
                ContigId(0),
                bases.len() as u64,
                &format!("micro {side} of a satellite, 20 bp"),
            );
            assert_eq!(
                kinds(&regions),
                vec!["generic", "satellite", "generic"],
                "micro {side}: the satellite swallows it — no bundle, no locus, and \
                 (crucially) no second region over the same bases"
            );

            // The satellite **grew to cover the microsatellite**: it is not merely that
            // the tract went missing. The absorbed span reaches past the array's own
            // 1.3 kb, gap included.
            let satellite = regions
                .iter()
                .find(|r| matches!(r.kind, RegionKind::Satellite))
                .expect("one satellite");
            assert!(
                satellite.region.len() >= 1300 + 20 + 24,
                "micro {side}: the satellite must EXPAND over the gap and the tract — \
                 array (1300) + gap (20) + tract (24); got {}",
                satellite.region.len()
            );

            // --- the control: 200 bp of gap, and the flank is clean ---
            let bases = micro_near_satellite(micro_left, 200);
            let regions = partition_resident("chr1", ContigId(0), &bases, &config);
            assert_partitions(
                &regions,
                ContigId(0),
                bases.len() as u64,
                &format!("micro {side} of a satellite, 200 bp"),
            );
            assert!(
                kinds(&regions).contains(&"locus"),
                "micro {side}: 200 bp away the SAME tract is a locus — so the absorption \
                 above is the rule, not a tract that was never admissible"
            );
            assert!(kinds(&regions).contains(&"satellite"));
        }
    }

    /// Absorption must not depend on the window either: the windowed walk agrees with
    /// the oracle on both arms and both sides. (The absorbed feature and the array are
    /// within `flank_bp`, so they are one block — this is what pins that.)
    #[test]
    fn absorption_is_window_invariant() {
        for micro_left in [true, false] {
            for gap in [20usize, 200] {
                let bases = micro_near_satellite(micro_left, gap);
                let resident =
                    partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
                for window_bp in [100u64, 700, 1000] {
                    let config = TypedRegionConfig {
                        window_bp: Bp(window_bp),
                        ..TypedRegionConfig::default()
                    };
                    let windowed =
                        partition_windowed(reference_over("chr1", &bases), ContigId(0), &config)
                            .expect("fetch");
                    assert_eq!(
                        windowed, resident,
                        "micro_left={micro_left} gap={gap} window_bp={window_bp}"
                    );
                }
            }
        }
    }

    /// **The cap applies to the CLEANED coverage, not the raw** (spec §2.4) — pinned at
    /// last, and it took D3's fixture to do it.
    ///
    /// D1 could only assert this by inspection: mutation showed `coverage_runs(&raw)`
    /// passing the entire suite, `.cat` parity included, because `raw ⊇ cleaned` and no
    /// fixture then existed where the extra coverage changed a *satellite*. The
    /// windowing fixture is the first that does — 6 kb of aperiodic filler carries real
    /// scanner noise (1181 raw runs against 5 cleaned), and one noise interval abuts the
    /// array, moving the satellite's edge.
    ///
    /// Written against the two candidate computations rather than a coordinate literal:
    /// it must fail when the cap moves to the raw set, and **not** when the detector's
    /// phase shifts by a base, which is the detector's business (`scanner_parity`).
    ///
    /// **What this does not do.** The difference it catches is a boundary, not the
    /// failure the spec's argument is about — noise unioning past 1 kb, declaring a
    /// satellite of its own, and swallowing the loci underneath. No fixture produces
    /// that. This pins the choice; it does not demonstrate the stakes.
    ///
    /// It also pins the walk's **second** copy of the decision
    /// (`BlockWalk::absorb`'s `coverage_in_core(cleaned)`), but only via
    /// window-invariance, which compares the two walks — so it catches either site
    /// moving, and not both moving together (verified).
    #[test]
    fn the_satellite_cap_applies_to_the_cleaned_coverage_not_the_raw() {
        let bases = windowing_fixture();
        let config = TypedRegionConfig::default();
        let raw = find_tandem_repeats(&bases, config.periods, &config.scan);
        let cleaned = segment_criteria::prefilter(&raw, &config.criteria);

        // The permissive scanner's noise must genuinely be here, or the two sets are the
        // same and everything below passes for free.
        assert!(
            coverage_runs(&raw).len() > 100 * coverage_runs(&cleaned).len(),
            "the fixture must carry raw noise: {} raw runs vs {} cleaned",
            coverage_runs(&raw).len(),
            coverage_runs(&cleaned).len()
        );

        let over_cap = |intervals: &[RepeatInterval]| -> Vec<(u64, u64)> {
            coverage_runs(intervals)
                .into_iter()
                .filter(|r| r.len() > config.max_repeat_len.get())
                .map(|r| (r.start, r.end))
                .collect()
        };
        let from_raw = over_cap(&raw);
        let from_cleaned = over_cap(&cleaned);
        assert_ne!(
            from_raw, from_cleaned,
            "this fixture must DISCRIMINATE the two, or the assertion below cannot fail \
             — which is exactly how the claim went untested through D1 and D2"
        );

        let satellites: Vec<(u64, u64)> = partition_resident("chr1", ContigId(0), &bases, &config)
            .iter()
            .filter(|r| matches!(r.kind, RegionKind::Satellite))
            .map(|r| (r.region.start.get(), r.region.end.get()))
            .collect();
        assert_eq!(
            satellites, from_cleaned,
            "the walk's satellites are the CLEANED coverage's over-cap runs: capping the \
             raw set would let detector noise decide where an array begins — and, with 1 \
             kb of it in a row, that an array exists at all (spec §2.4)"
        );
    }

    /// **The acceptance test: `window_bp` is a memory knob and must not move the
    /// output** (spec §2.3), proven by matching [`partition_resident`] — the simplest
    /// implementation, which has no windows to get wrong (impl plan, Milestone D).
    ///
    /// Several window sizes, because a windowing bug is usually a bug at one specific
    /// alignment of feature to boundary: 1000 cuts the fixture's features (above), 700
    /// and 333 cut it somewhere else and give the carries many more chances to drop or
    /// double something, and 100_000 is the degenerate one-window case.
    #[test]
    fn windowed_matches_the_resident_oracle_at_every_window_size() {
        let bases = windowing_fixture();
        let base = TypedRegionConfig::default();
        let resident = partition_resident("chr1", ContigId(0), &bases, &base);
        assert_partitions(&resident, ContigId(0), bases.len() as u64, "resident");

        // The fixture must contain what it claims to, or every comparison below is
        // "empty == empty". This is the guard on the guard.
        let count = |rs: &[TypedRegion], f: fn(&RegionKind) -> bool| {
            rs.iter().filter(|r| f(&r.kind)).count()
        };
        assert!(
            count(&resident, |k| matches!(k, RegionKind::SsrSegment(_))) >= 2,
            "the fixture must classify loci: {resident:#?}"
        );
        assert_eq!(
            count(&resident, |k| matches!(k, RegionKind::SsrBundle { .. })),
            1,
            "the fixture must bundle the close pair"
        );
        assert_eq!(
            count(&resident, |k| matches!(k, RegionKind::Satellite)),
            1,
            "the fixture must have a satellite"
        );

        for window_bp in [1000u64, 700, 333, 100_000] {
            let config = TypedRegionConfig {
                window_bp: Bp(window_bp),
                ..base.clone()
            };
            let windowed = partition_windowed(reference_over("chr1", &bases), ContigId(0), &config)
                .expect("fetch");
            assert_partitions(
                &windowed,
                ContigId(0),
                bases.len() as u64,
                &format!("windowed at {window_bp}"),
            );
            assert_eq!(
                windowed, resident,
                "window_bp = {window_bp} changed the output — it is a memory knob and \
                 nothing else (spec §2.3)"
            );
        }
    }

    /// The same, on **real sequence** rather than a crafted fixture: every contig of
    /// the golden reference, windowed small, must equal the resident walk that
    /// reproduces the trf-mod-built catalog (`the_resident_partition_reproduces_the_
    /// golden_catalog`). So the windowed walk inherits `.cat` parity — through the
    /// oracle, which is the point of having one.
    ///
    /// Crafted fixtures place features where the author thought to; a real assembly
    /// places them where they are, at densities and spacings nobody chose.
    #[test]
    fn windowed_matches_the_resident_oracle_on_the_golden_reference() {
        use std::fs::File;
        use std::io::BufReader;
        use std::path::Path;

        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("data")
            .join("tandem_repeat")
            .join("synthetic_ref.fa");
        let file = File::open(fixture).unwrap();
        let mut reader = noodles_fasta::io::Reader::new(BufReader::new(file));

        let mut checked = 0;
        for result in reader.records() {
            let rec = result.unwrap();
            let name = String::from_utf8_lossy(rec.name()).into_owned();
            let bases: Vec<u8> = rec.sequence().as_ref().to_vec();
            // One reference per contig, so the id is 0 in each — the walk is
            // per-contig, and `E3` is where a multi-contig reference is walked whole.
            let contig = ContigId(0);

            let resident = partition_resident(&name, contig, &bases, &TypedRegionConfig::default());
            for window_bp in [500u64, 1500] {
                let config = TypedRegionConfig {
                    window_bp: Bp(window_bp),
                    ..TypedRegionConfig::default()
                };
                let windowed = partition_windowed(reference_over(&name, &bases), contig, &config)
                    .expect("fetch");
                assert_eq!(
                    windowed, resident,
                    "{name} at window_bp = {window_bp} diverged from the resident oracle"
                );
            }
            checked += 1;
        }
        assert!(checked > 0, "the golden reference must have contigs");
    }

    /// **Maximality across windows** (spec §2.3): a generic run longer than a window
    /// is **one** region, not one per window.
    ///
    /// This is the open-generic-run carry, and it is a correctness requirement rather
    /// than tidiness: the pileup mints loci *inside* a `Generic` region, so a run split
    /// at a window edge makes an indel spanning that edge callable by neither half — it
    /// never appears, and nothing fails.
    #[test]
    fn a_generic_run_longer_than_a_window_is_one_region() {
        let bases = windowing_fixture();
        let window_bp = 100;
        let config = TypedRegionConfig {
            window_bp: Bp(window_bp),
            ..TypedRegionConfig::default()
        };
        let regions = partition_windowed(reference_over("chr1", &bases), ContigId(0), &config)
            .expect("fetch");
        assert_partitions(&regions, ContigId(0), bases.len() as u64, "maximality");

        // The stretch between the tract at 3500 and the array at 4000 is ~480 bp — five
        // windows — and the one before 990 is ~940. Both must be single regions.
        let longest = regions
            .iter()
            .filter(|r| matches!(r.kind, RegionKind::Generic))
            .map(|r| r.region.len())
            .max()
            .expect("the fixture has generic sequence");
        assert!(
            longest > window_bp * 3,
            "a generic run must span windows whole: the longest is {longest} bp at \
             window_bp = {window_bp}"
        );
    }

    /// **The contig-end clamp guard, at the walk level** (spec §2.6) — and the
    /// **provenance obligation**, discharged and pinned.
    ///
    /// `classify` clamps a locus's flanks at the *contig's* ends and drops a locus whose
    /// flank clamped to nothing. Hand it the **window's** end as the contig's length —
    /// production's exact mistake, and arithmetically legal, so `classify`'s own guard
    /// cannot catch it — and every locus within `flank_bp` of every window boundary
    /// silently vanishes, a different set for each `window_bp`.
    ///
    /// The fixture's tract at 990 straddles the 1000 core edge; the one at 3500 sits in
    /// the middle of a core. Only reading `contig_len` from the reference's contig table
    /// keeps the first one. Mutation-verified: passing `window.core.end` (or
    /// `window.fetched.end`) instead of the table's length drops it and this fails.
    #[test]
    fn a_locus_astride_a_window_edge_is_not_dropped() {
        let bases = windowing_fixture();
        let config = TypedRegionConfig {
            window_bp: Bp(1000),
            ..TypedRegionConfig::default()
        };
        let regions = partition_windowed(reference_over("chr1", &bases), ContigId(0), &config)
            .expect("fetch");

        let loci: Vec<u64> = regions
            .iter()
            .filter_map(|r| match &r.kind {
                RegionKind::SsrSegment(l) => Some(l.start()),
                _ => None,
            })
            .collect();
        // 1-based, and the detector may shift the tract's edge by a base or two — its
        // business, not this test's (`scanner_parity` characterises the same wobble).
        assert!(
            loci.iter().any(|&s| (989..=1000).contains(&s)),
            "the locus astride the 1000 core edge must survive: loci at {loci:?}"
        );
        assert!(
            loci.iter().any(|&s| (3499..=3510).contains(&s)),
            "the interior control locus must be there too: loci at {loci:?}"
        );
    }

    /// A locus is emitted **exactly once**, whichever core holds its start — the
    /// attribution rule. Doubling is the other half of the failure the carries can
    /// cause, and the invariant catches it as an overlap; this says so directly.
    #[test]
    fn every_locus_is_emitted_exactly_once() {
        let bases = windowing_fixture();
        for window_bp in [100u64, 333, 1000] {
            let config = TypedRegionConfig {
                window_bp: Bp(window_bp),
                ..TypedRegionConfig::default()
            };
            let regions = partition_windowed(reference_over("chr1", &bases), ContigId(0), &config)
                .expect("fetch");
            let mut starts: Vec<u64> = regions
                .iter()
                .filter(|r| matches!(r.kind, RegionKind::SsrSegment(_)))
                .map(|r| r.region.start.get())
                .collect();
            let before = starts.len();
            starts.dedup();
            assert_eq!(
                before,
                starts.len(),
                "window_bp = {window_bp}: a locus doubled"
            );
        }
    }

    /// A reference read that fails is fatal, and surfaces as `Err` — the walk never
    /// scans a hole (spec §8.2). An unknown contig is the reachable form of that here.
    #[test]
    fn an_unreadable_contig_is_an_error_not_an_empty_partition() {
        let err = partition_windowed(
            reference_over("chr1", b"ACGT"),
            ContigId(9),
            &TypedRegionConfig::default(),
        )
        .expect_err("contig 9 does not exist");
        assert!(matches!(err, TypedRegionError::Reference(_)));
    }

    #[test]
    fn a_zero_length_contig_windows_to_nothing() {
        let regions = partition_windowed(
            reference_over("empty", b""),
            ContigId(0),
            &TypedRegionConfig::default(),
        )
        .unwrap();
        assert!(regions.is_empty());
    }

    /// **A margin narrower than the bundle radius silently un-bundles**, so the walk
    /// refuses it — with an **error**, and before any work. Without the guard, a core
    /// tract whose neighbour fell outside the window would be classified as a clean locus:
    /// no error, and a different answer for every `window_bp`.
    ///
    /// It is a `Result`, not an `assert!`, because both knobs are command-line flags
    /// (`typed_regions_cli.md` §2.1) and a typo must not panic. A2's rule — a swept knob's
    /// guard must survive `--release` — is what a `debug_assert` would have broken, and a
    /// `Result` keeps unconditionally.
    #[test]
    fn a_detection_margin_narrower_than_the_bundle_radius_is_refused() {
        let config = TypedRegionConfig {
            max_repeat_len: Bp(10),
            ..TypedRegionConfig::default()
        };
        let err = partition_windowed(
            reference_over("chr1", &windowing_fixture()),
            ContigId(0),
            &config,
        )
        .expect_err("a margin narrower than flank_bp is refused");
        assert!(
            matches!(
                err,
                TypedRegionError::MarginNarrowerThanFlank {
                    max_repeat_len: 10,
                    flank_bp: 50
                }
            ),
            "and it names both numbers, so the message says which flag to move: {err}"
        );
        assert!(
            err.to_string()
                .contains("must not be narrower than flank_bp"),
            "the operator-facing message survives: {err}"
        );
    }

    // ---- E1: the iterator surface ----------------------------------------

    /// `over_regions` walks **every** span, in genomic order, gap-free across contigs —
    /// and each contig's regions are exactly what walking that contig alone gives.
    ///
    /// The multi-contig case is the one the per-contig tests cannot reach: the walk has
    /// to put down one contig and pick up the next without carrying anything across
    /// (`emitted_upto` and the block carries are per span, and a leak would show as a
    /// gap or an overlap at the seam).
    #[test]
    fn over_regions_walks_every_contig_in_order() {
        let a = windowing_fixture();
        let b = micro_near_satellite(true, 20);
        let contigs = &[
            ContigBounds {
                name: "chrA",
                length: a.len() as u32,
            },
            ContigBounds {
                name: "chrB",
                length: b.len() as u32,
            },
        ];
        let reference = InMemoryRefSeq::from_named_contigs(vec![
            ("chrA".to_string(), a.clone()),
            ("chrB".to_string(), b.clone()),
        ]);
        let config = TypedRegionConfig {
            window_bp: Bp(700),
            ..TypedRegionConfig::default()
        };

        let regions: Vec<TypedRegion> = TypedRegionIterator::over_regions(
            reference,
            GenomeRegions::whole_contigs(contigs),
            config.clone(),
        )
        .expect("valid spans")
        .collect::<Result<_, _>>()
        .expect("no read fails");

        // Each contig's slice of the output partitions that contig...
        let (from_a, from_b): (Vec<_>, Vec<_>) = regions
            .iter()
            .cloned()
            .partition(|r| r.region.contig == ContigId(0));
        assert_partitions(
            &from_a,
            ContigId(0),
            a.len() as u64,
            "chrA via over_regions",
        );
        assert_partitions(
            &from_b,
            ContigId(1),
            b.len() as u64,
            "chrB via over_regions",
        );

        // ...and is identical to walking it on its own: nothing carries across the seam.
        assert_eq!(
            from_a,
            partition_resident("chrA", ContigId(0), &a, &config),
            "chrA"
        );
        assert_eq!(
            from_b,
            partition_resident("chrB", ContigId(1), &b, &config),
            "chrB"
        );

        // Genomic order, contigs in table order.
        assert!(
            regions
                .windows(2)
                .all(|w| (w[0].region.contig, w[0].region.start)
                    <= (w[1].region.contig, w[1].region.start)),
            "regions must come out in genomic order"
        );
    }

    /// The **running tally**: readable mid-walk, complete once exhausted (spec §3.1).
    ///
    /// Checked against the regions actually yielded rather than against literals — the
    /// counts are a claim *about the output*, so anything else would be two independent
    /// guesses at the fixture.
    #[test]
    fn counts_tally_the_regions_yielded() {
        let bases = windowing_fixture();
        let reference = reference_over("chr1", &bases);
        let mut iter =
            TypedRegionIterator::over_contig(reference, ContigId(0), TypedRegionConfig::default())
                .expect("valid contig");

        assert_eq!(
            *iter.counts(),
            TypedRegionCounts::default(),
            "nothing walked, nothing counted"
        );

        let mut yielded: Vec<TypedRegion> = Vec::new();
        // Mid-walk the tally must already describe what has come out so far.
        for r in iter.by_ref().take(3) {
            yielded.push(r.expect("no read fails"));
        }
        assert_eq!(
            iter.counts().ssr_loci + iter.counts().generic + iter.counts().satellites,
            3,
            "the tally counts what has been HANDED OUT, not what has been scanned"
        );

        for r in iter.by_ref() {
            yielded.push(r.expect("no read fails"));
        }
        let counts = iter.counts();
        let kind_count =
            |f: fn(&RegionKind) -> bool| yielded.iter().filter(|r| f(&r.kind)).count() as u64;
        assert_eq!(counts.spans, 1);
        assert_eq!(
            counts.ssr_loci,
            kind_count(|k| matches!(k, RegionKind::SsrSegment(_)))
        );
        assert_eq!(
            counts.ssr_bundles,
            kind_count(|k| matches!(k, RegionKind::SsrBundle { .. }))
        );
        assert_eq!(
            counts.generic,
            kind_count(|k| matches!(k, RegionKind::Generic))
        );
        assert_eq!(
            counts.satellites,
            kind_count(|k| matches!(k, RegionKind::Satellite))
        );
        assert!(counts.satellites > 0 && counts.ssr_loci > 0 && counts.ssr_bundles > 0);

        let bp_of = |f: fn(&RegionKind) -> bool| -> u64 {
            yielded
                .iter()
                .filter(|r| f(&r.kind))
                .map(|r| r.region.len())
                .sum()
        };
        assert_eq!(
            counts.satellite_bp,
            bp_of(|k| matches!(k, RegionKind::Satellite))
        );
        assert_eq!(
            counts.ssr_bundle_bp,
            bp_of(|k| matches!(k, RegionKind::SsrBundle { .. }))
        );

        // **The gap `repeat_bp_with_no_locus` measures**, pinned exactly rather than
        // bounded: it is the contig's cleaned repeat coverage minus the bases that came
        // out as loci. Computed here from the resident path — an independent route to the
        // same number, so this fails if either the accumulate or the subtract is wrong.
        // (A `>=` bound was the first version, and it left the subtraction untested:
        // dropping it makes the number bigger, and bigger still satisfies `>=`.)
        let config = TypedRegionConfig::default();
        let raw = find_tandem_repeats(&bases, config.periods, &config.scan);
        let cleaned = segment_criteria::prefilter(&raw, &config.criteria);
        let coverage_bp: u64 = coverage_runs(&cleaned).iter().map(|r| r.len()).sum();
        let locus_bp = bp_of(|k| matches!(k, RegionKind::SsrSegment(_)));
        assert!(
            coverage_bp > 0 && locus_bp > 0,
            "the fixture must have both"
        );
        assert_eq!(
            counts.repeat_bp_with_no_locus,
            coverage_bp - locus_bp,
            "repeat coverage ({coverage_bp} bp) that yielded no locus ({locus_bp} bp of \
             it did)"
        );
        assert!(
            counts.repeat_bp_with_no_locus >= counts.ssr_bundle_bp,
            "and the bundled tracts are part of it"
        );
    }

    /// Walk a fixture and hand back the finished tally.
    fn counts_over(bases: &[u8], window_bp: u64) -> TypedRegionCounts {
        let config = TypedRegionConfig {
            window_bp: Bp(window_bp),
            ..TypedRegionConfig::default()
        };
        let mut iter =
            TypedRegionIterator::over_contig(reference_over("chr1", bases), ContigId(0), config)
                .expect("valid contig");
        for r in iter.by_ref() {
            r.expect("no read fails");
        }
        *iter.counts()
    }

    /// **The rejection breakdown must not depend on `window_bp`** (E1e).
    ///
    /// This is the whole reason `classify` hands rejections back **per record** instead of
    /// tallying them itself. It runs over a window's entire fetched slice, margins and
    /// all, so every window that can see a repeat rejects it again; a tally taken inside
    /// would count one repeat once per window and the number would move with a **memory
    /// knob**. Attributing each record to the core holding its start is what makes the
    /// count a fact about the genome instead of about the walk.
    ///
    /// The whole tally is compared, not just the breakdown: every count here is a claim
    /// about the reference, and `window_bp` may not touch any of them (spec §2.3).
    #[test]
    fn the_tally_does_not_depend_on_the_window() {
        // The windowing fixture, plus a tract abutting base 1 — which classification turns
        // down for having no left flank (spec §2.6). That rejection is what gives the
        // breakdown something to count: of classification's five gates, the pre-filter makes
        // three unreachable from the walk and this is the reachable one that a fixture
        // can be sure of (`segment_criteria::the_pre_filter_makes_three_of_classifications_gates_
        // unreachable_from_the_walk`).
        let mut bases = windowing_fixture();
        bases[0..20].copy_from_slice(b"ATATATATATATATATATAT");
        let baseline = counts_over(&bases, 100_000);

        // Or this compares zeroes — and it would, without the tract above.
        assert!(
            baseline.rejected_by_reason.total() > 0,
            "the fixture must reach a gate: {:?}",
            baseline.rejected_by_reason
        );
        assert!(baseline.ssr_loci > 0 && baseline.satellites > 0);

        for window_bp in [100u64, 333, 700, 1000] {
            assert_eq!(
                counts_over(&bases, window_bp),
                baseline,
                "window_bp = {window_bp} moved the tally — `window_bp` is a memory knob \
                 (spec §2.3). A rejection counted once per window that could SEE it is \
                 exactly how that happens: at window_bp = 100 the tract at base 1 is in \
                 the fetched slice of the first ten windows, and `classify` turns it down in \
                 every one of them"
            );
        }
    }

    /// A reference-read failure mid-walk is **fatal and in-stream**: `Some(Err(_))` once,
    /// then `None` forever (spec §8.2). `None` meaning both "done" and "chromosome 7
    /// failed" would silently un-call the rest of the genome.
    ///
    /// The failure is injected by a reference that reads a prefix and then refuses, so
    /// the error lands **mid-walk** — a constructor-time failure would prove nothing
    /// about the stream.
    #[test]
    fn a_read_failure_mid_walk_is_yielded_once_and_then_the_iterator_is_done() {
        /// Fails every read starting past `fail_from` (1-based).
        struct FailsLate {
            inner: InMemoryRefSeq,
            fail_from: u64,
        }
        impl RefSeq for FailsLate {
            fn fetch_into(
                &self,
                contig: ContigId,
                start: u64,
                len: u64,
                dst: &mut Vec<u8>,
            ) -> Result<(), RefSeqError> {
                self.inner.fetch_into(contig, start, len, dst)
            }
        }
        impl RawRefSeq for FailsLate {
            fn fetch_raw_into(
                &self,
                contig: ContigId,
                start: u64,
                len: u64,
                dst: &mut Vec<u8>,
            ) -> Result<(), RefSeqError> {
                if start > self.fail_from {
                    return Err(RefSeqError::InvalidStart);
                }
                self.inner.fetch_raw_into(contig, start, len, dst)
            }
        }
        impl ContigTable for FailsLate {
            fn contigs(&self) -> &crate::fasta::ContigList {
                self.inner.contigs()
            }
        }
        impl EvictableRefSeq for FailsLate {
            fn evict_before(&mut self, _pos: u64) {}
        }

        let bases = windowing_fixture();
        let reference = FailsLate {
            inner: reference_over("chr1", &bases),
            fail_from: 2000,
        };
        let config = TypedRegionConfig {
            window_bp: Bp(500),
            ..TypedRegionConfig::default()
        };
        let mut iter = TypedRegionIterator::over_contig(reference, ContigId(0), config)
            .expect("the contig table is readable; only the BASES fail");

        let mut errors = 0;
        let mut items = 0;
        for r in iter.by_ref() {
            match r {
                Ok(_) => items += 1,
                Err(e) => {
                    assert!(matches!(e, TypedRegionError::Reference(_)));
                    errors += 1;
                }
            }
        }
        assert_eq!(errors, 1, "the failure is reported exactly once");
        assert!(items > 0, "the windows before the failure still yielded");

        // Fused: done means done. A caller polling on is not handed a partial walk that
        // looks complete.
        assert!(iter.next().is_none(), "fused: nothing after the error");
        assert!(iter.next().is_none());
    }

    /// **The walk releases the bases it has passed** — the difference between "holds one
    /// window, never a contig" (spec §7) as a claim and as a fact.
    ///
    /// Nothing else can catch this. The reference impls the other tests use hold their
    /// bytes outright, so their eviction is an honest no-op and a walk that never evicted
    /// would pass every test in this file — while a real `WindowedRefSeq`'s buffer grew to
    /// the whole contig, silently, which is exactly the memory profile the windowed walk
    /// exists to avoid. So the reference here records the asks.
    #[test]
    fn the_walk_evicts_the_bases_it_has_passed() {
        struct EvictionSpy {
            inner: InMemoryRefSeq,
            evictions: std::rc::Rc<std::cell::RefCell<Vec<u64>>>,
        }
        impl RefSeq for EvictionSpy {
            fn fetch_into(
                &self,
                contig: ContigId,
                start: u64,
                len: u64,
                dst: &mut Vec<u8>,
            ) -> Result<(), RefSeqError> {
                self.inner.fetch_into(contig, start, len, dst)
            }
        }
        impl RawRefSeq for EvictionSpy {
            fn fetch_raw_into(
                &self,
                contig: ContigId,
                start: u64,
                len: u64,
                dst: &mut Vec<u8>,
            ) -> Result<(), RefSeqError> {
                self.inner.fetch_raw_into(contig, start, len, dst)
            }
        }
        impl ContigTable for EvictionSpy {
            fn contigs(&self) -> &crate::fasta::ContigList {
                self.inner.contigs()
            }
        }
        impl EvictableRefSeq for EvictionSpy {
            fn evict_before(&mut self, pos: u64) {
                self.evictions.borrow_mut().push(pos);
            }
        }

        let bases = windowing_fixture();
        let evictions = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));
        let reference = EvictionSpy {
            inner: reference_over("chr1", &bases),
            evictions: std::rc::Rc::clone(&evictions),
        };
        let config = TypedRegionConfig {
            window_bp: Bp(1000),
            ..TypedRegionConfig::default()
        };
        let regions: Vec<TypedRegion> =
            TypedRegionIterator::over_contig(reference, ContigId(0), config)
                .expect("valid contig")
                .collect::<Result<_, _>>()
                .expect("no read fails");
        assert!(!regions.is_empty());

        let evictions = evictions.borrow();
        assert_eq!(
            evictions.len(),
            6,
            "one per window: a 6 kb contig at window_bp = 1000. Evicting once at the end \
             would bound nothing"
        );
        assert!(
            evictions.windows(2).all(|w| w[0] <= w[1]),
            "eviction only ever moves forward — the walk never looks back: {evictions:?}"
        );
        assert!(
            *evictions.last().unwrap() > 4000,
            "and it follows the walk to the end of the contig: {evictions:?}"
        );
    }

    /// A span naming a contig the reference does not have fails **at construction**, not
    /// mid-walk. The arch doc expected this constructor to be infallible because
    /// `GenomeRegions` validates; it validates against *a* contig table, and nothing ties
    /// that one to the reference's.
    #[test]
    fn a_span_the_reference_does_not_have_is_rejected_up_front() {
        let contigs = &[
            ContigBounds {
                name: "chr1",
                length: 100,
            },
            ContigBounds {
                name: "chr2",
                length: 100,
            },
        ];
        // The spans know two contigs; the reference has one.
        let reference = reference_over("chr1", &filler(100));
        let result = TypedRegionIterator::over_regions(
            reference,
            GenomeRegions::whole_contigs(contigs),
            TypedRegionConfig::default(),
        );
        assert!(
            matches!(result, Err(TypedRegionError::Reference(_))),
            "chr2 is not in this reference, and that must be caught before any walk"
        );
    }

    // ---- E2: scan wider than you emit ------------------------------------

    /// Walk `bases` restricted to the 1-based inclusive spans `want`.
    fn walk_bed(bases: &[u8], want: &[(u64, u64)], window_bp: u64) -> Vec<TypedRegion> {
        let requested: Vec<GenomeRegion> = want
            .iter()
            .map(|&(s, e)| GenomeRegion {
                contig: ContigId(0),
                start: Position(s),
                end: Position(e),
            })
            .collect();
        let config = TypedRegionConfig {
            window_bp: Bp(window_bp),
            ..TypedRegionConfig::default()
        };
        TypedRegionIterator::over_spans(reference_over("chr1", bases), requested, config)
            .expect("valid spans")
            .collect::<Result<_, _>>()
            .expect("no read fails")
    }

    /// **BED-invariance — the acceptance test** (spec §2.5): *"whether a base is STR /
    /// bundle / generic / satellite must not depend on the BED. The BED chooses what you
    /// are shown, never what things are."*
    ///
    /// Every requested span is compared against the **whole-genome** run restricted to the
    /// same coordinates. The spans are chosen to cut through everything the fixture has —
    /// a BED edge 5 bp inside the bundle's cluster, one mid-satellite, one splitting a
    /// generic run — because an edge in the middle of nowhere tests nothing.
    ///
    /// This is what "scan wider than you emit" buys, and it is not free: the naive walk
    /// (scan only what was asked for) gets the *bundle* wrong, because the neighbour that
    /// makes it a bundle sits outside the span.
    #[test]
    fn a_bed_does_not_change_what_things_are() {
        let bases = windowing_fixture();
        let whole = partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());

        // Spans that cut through real structure: across the tract at 990; across the
        // bundle at 1990..2040; inside the satellite at 4000..5200; and a plain interior.
        for want in [
            (900u64, 1100u64),
            (1995, 2100),
            (4500, 4800),
            (3000, 3600),
            (1, 6000),
        ] {
            for window_bp in [333u64, 1000] {
                let got = walk_bed(&bases, &[want], window_bp);

                // What the whole-genome run says about those same bases.
                for r in &got {
                    for pos in [r.region.start.get(), r.region.end.get()] {
                        // Only ask about bases the user actually requested: an object
                        // emitted whole reaches past the edge, which is the point.
                        if pos < want.0 || pos > want.1 {
                            continue;
                        }
                        let truth = whole
                            .iter()
                            .find(|w| w.region.contains(Position(pos)))
                            .unwrap_or_else(|| panic!("the whole-genome run covers {pos}"));
                        assert_eq!(
                            std::mem::discriminant(&r.kind),
                            std::mem::discriminant(&truth.kind),
                            "BED {want:?} at window_bp = {window_bp}: base {pos} is {:?} \
                             with the BED and {:?} without it — the BED changed what a \
                             base IS (spec §2.5)",
                            r.kind,
                            truth.kind
                        );
                    }
                }

                // And the loci are the same objects, not merely the same kind.
                let loci_in = |rs: &[TypedRegion]| -> Vec<(u64, u64)> {
                    rs.iter()
                        .filter(|r| matches!(r.kind, RegionKind::SsrSegment(_)))
                        .map(|r| (r.region.start.get(), r.region.end.get()))
                        .filter(|(s, _)| *s >= want.0 && *s <= want.1)
                        .collect()
                };
                assert_eq!(
                    loci_in(&got),
                    loci_in(&whole),
                    "BED {want:?} at window_bp = {window_bp}: the loci starting inside the \
                     span must be exactly the whole-genome run's"
                );
            }
        }
    }

    /// The requested span is **tiled exactly** — the partition invariant, restated where
    /// it belongs once a BED is involved: over the *effective* region, which is what the
    /// user asked for grown to hold a straddling object whole (spec §2.5).
    #[test]
    fn a_bed_span_is_tiled_by_what_comes_back() {
        let bases = windowing_fixture();
        for want in [(900u64, 1100u64), (1995, 2100), (3000, 3600)] {
            let got = walk_bed(&bases, &[want], 333);
            assert!(!got.is_empty());

            // Contiguous and non-overlapping, and covering every requested base. The
            // effective region may start before / end after the request (an object
            // straddling the edge), so the *bounds* are asserted as reaching at least the
            // request, not equalling it.
            assert!(got[0].region.start.get() <= want.0);
            assert!(got.last().unwrap().region.end.get() >= want.1);
            for pair in got.windows(2) {
                assert_eq!(
                    pair[1].region.start.get(),
                    pair[0].region.end.get() + 1,
                    "gap or overlap inside a BED span: {want:?}"
                );
                assert_ne!(
                    std::mem::discriminant(&pair[0].kind),
                    std::mem::discriminant(&pair[1].kind),
                    "two consecutive regions share a kind — MAXIMALITY (spec §2.3)"
                );
            }
        }
    }

    /// **`Generic` is clipped to the user's edge; every *finding* straddling it comes back
    /// whole** — locus, bundle **and satellite** (spec §2.5, owner 2026-07-17;
    /// [`clips_at_a_bed_edge`]).
    #[test]
    fn generic_clips_at_the_edge_and_findings_come_back_whole() {
        let bases = windowing_fixture();

        // A span whose edges fall in plain generic sequence: nothing straddles, so the
        // partition starts and ends exactly where the user asked.
        let got = walk_bed(&bases, &[(3000, 3600)], 333);
        assert_eq!(
            got[0].region.start,
            Position(3000),
            "Generic is clipped to the user's edge, not grown to the scan span"
        );
        assert!(matches!(got[0].kind, RegionKind::Generic));
        assert_eq!(got.last().unwrap().region.end, Position(3600));

        // A span wholly INSIDE the 1.2 kb array: the satellite comes back whole, reaching
        // past both requested edges.
        //
        // **Clipping it would produce a `Satellite` region of 301 bp** — a span that
        // contradicts the `max_repeat_len` (1 kb) test that produced the label. The extent
        // is the claim: "an array too long to be a microsatellite". That is what E2 got
        // wrong by reasoning from the type (`Satellite` carries no payload, so nothing
        // could be left misdescribed) instead of from the meaning.
        let got = walk_bed(&bases, &[(4500, 4800)], 333);
        assert_eq!(got.len(), 1, "the whole span is inside the array: {got:#?}");
        assert!(matches!(got[0].kind, RegionKind::Satellite));
        assert!(
            got[0].region.start.get() < 4500 && got[0].region.end.get() > 4800,
            "the Satellite is emitted WHOLE, past both edges: {:?}",
            got[0].region
        );
        assert!(
            got[0].region.len() > 1000,
            "and its span is over the cap that made it a satellite — which a clipped one \
             ({} bp of request) could not be",
            4800 - 4500 + 1
        );
        // The same region the whole-genome run reports, not a version of it.
        let whole = partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        let truth = whole
            .iter()
            .find(|r| matches!(r.kind, RegionKind::Satellite))
            .expect("the fixture has one");
        assert_eq!(got[0].region, truth.region, "the SAME satellite");

        // A span cutting into the bundle's cluster: the bundle carries its member tracts,
        // so clipping it would leave them outside their own region. It comes back whole,
        // and the effective region grows to hold it.
        let got = walk_bed(&bases, &[(1995, 2100)], 333);
        let bundle = got
            .iter()
            .find(|r| matches!(r.kind, RegionKind::SsrBundle { .. }))
            .expect("the cluster at 1990..2040 straddles this edge");
        assert!(
            bundle.region.start.get() < 1995,
            "the bundle is emitted WHOLE, reaching back past the requested edge to {} \
             (the effective region grows to hold it)",
            bundle.region.start.get()
        );
        let RegionKind::SsrBundle { tracts } = &bundle.kind else {
            unreachable!()
        };
        assert!(
            tracts
                .iter()
                .all(|t| t.start + 1 >= bundle.region.start.get()
                    && t.end <= bundle.region.end.get()),
            "every member tract is inside the hull — which is what clipping would break"
        );
    }

    /// A locus straddling the edge likewise comes back **whole**, with its bases intact —
    /// half a locus is not a locus.
    #[test]
    fn a_locus_straddling_the_bed_edge_is_emitted_whole() {
        let bases = windowing_fixture();
        // The tract at 990..1010 (1-based 991..1010); ask for a span ending inside it.
        let got = walk_bed(&bases, &[(900, 1000)], 333);
        let locus = got
            .iter()
            .find_map(|r| match &r.kind {
                RegionKind::SsrSegment(l) => Some(l),
                _ => None,
            })
            .expect("the tract at ~991 straddles the requested edge at 1000");
        assert!(
            locus.end() > 1000,
            "the locus reaches past the user's edge ({}), whole",
            locus.end()
        );
        // The same object the whole-genome run builds, bases and all.
        let whole = partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        let truth = whole
            .iter()
            .find_map(|r| match &r.kind {
                RegionKind::SsrSegment(l) if l.start() == locus.start() => Some(l),
                _ => None,
            })
            .expect("the same locus exists in the whole-genome run");
        assert_eq!(locus, truth, "and it is the SAME locus, not a clipped one");
    }

    /// **Two requested spans on one contig share its scan**: the walk scans the contig
    /// once, and the ground between them — which the user did not ask for — must not come
    /// back (spec §2.5).
    ///
    /// This is the case that makes `Generic` clip against *each* requested span rather than
    /// against the scan span: one generic run covers both, and it has to come back as two
    /// regions with the gap dropped. Since the scan set became whole contigs it is no
    /// longer an edge case at all — **every** pair of spans on a contig is this case, which
    /// is a good reason for the rule to be the general one.
    #[test]
    fn two_spans_sharing_a_scan_span_do_not_leak_the_gap_between_them() {
        let bases = windowing_fixture();
        let got = walk_bed(&bases, &[(3000, 3200), (3400, 3600)], 333);

        assert!(
            got.iter()
                .all(|r| (3000..=3200).contains(&r.region.start.get())
                    || (3400..=3600).contains(&r.region.start.get())),
            "nothing may start in the gap the user did not ask for: {got:#?}"
        );
        let covered: u64 = got.iter().map(|r| r.region.len()).sum();
        assert_eq!(
            covered,
            201 + 201,
            "exactly the two requested spans come back — the gap between them is scanned \
             (they share a scan span) and not emitted"
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
