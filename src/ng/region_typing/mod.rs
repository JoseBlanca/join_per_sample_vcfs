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
use crate::ng::tandem_repeat::{PeriodRange, RepeatInterval, ScanParams};
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
