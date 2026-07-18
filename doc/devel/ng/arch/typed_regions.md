# ng step 3 — typed-region generator: types & interfaces

*Architecture draft (2026-07-16), the code-facing companion to
[`../spec/typed_regions.md`](../spec/typed_regions.md). It gives the types and signatures as they
appear in code; **every *why* points back to the spec** — this doc does not re-argue them. It
supersedes the step-3 sketch (`LocusRouter`/`LocusSource`) in the shared
[`ng_step_interfaces.md`](ng_step_interfaces.md) §3, and the `locus_router/` entry in
[`module_layout.md`](module_layout.md); the edits to those two docs land alongside this one.
Naming: **STR** in prose, `ssr` in code.*

---

## Revision — 2026-07-16: ng owns its copies

Tracks the spec's Revision of the same date: **production (`src/ssr/`, `src/regions.rs`) is frozen**;
ng copies what it needs rather than reshaping production in place, and ng must not depend on trf-mod.
What changed in this doc: the module is a folder; `SsrSegment` holds **ng's** `SsrSegment`, not
`ssr::types::Locus`; `TypedRegionConfig.catalog: CatalogParams` becomes
`criteria: SsrSegmentCriteria` (ng's own); and the reconciliation table's "use directly / reuse,
widened" rows become "copy into ng". The interface, the iterator contract, and every decision below
them are unaffected — none of them ever depended on where the code lived.

## Module home

`src/ng/region_typing/` — a folder, though **not** because of a bake-off (there is none; spec §6):
it carries the ~500-line port of the classification policy alongside the walk, and those are two concerns
with two test suites.

- `region_typing/mod.rs` — `TypedRegion`, `RegionKind`, `TypedRegionConfig`, `TypedRegionCounts`,
  `TypedRegionError`, `GenomeRegions`, `TypedRegionIterator`, the walk.
- `region_typing/segment_criteria.rs` — ng's port (spec §5): `Motif`, `SsrSegment`, `SsrSegmentCriteria`,
  `Classified`, the pre-filter, `classify`, and the differential test against production's `build_loci`
  (spec §8.0).

It is a **step** module, so it sits in the step layer of the `src/ng/` tree, replacing the
`locus_router/` folder the layout doc reserved. `GenomeRegion` and `Position` are shared vocabulary
and land in `ng::types`, not here.

## The types

The final shapes. Contracts are stated below each; the *why* is the cross-ref.

```rust
/// A genome region plus what the sequence there IS. The walk's output.
/// `region` is a field, not a per-variant repeat: every typed region has one, structurally, and it
/// is the one place ng's 1-based base is stated. (spec §3, §4)
pub struct TypedRegion { pub region: GenomeRegion, pub kind: RegionKind }

/// All four kinds are physical regions — this module produces no genetic locus (spec §1.1).
/// An `SsrSegment` is a stretch carrying an STR; it becomes a locus only if variable, downstream.
/// `Generic`/`Satellite` carry nothing — they *are* just their region.
pub enum RegionKind {
    SsrSegment(SsrSegment),                               // ng's own SsrSegment (segment_criteria.rs) — see below
    SsrBundle { tracts: Box<[RepeatInterval]> },   // >= 2, coordinate-ordered; hull = TypedRegion.region
    Generic,
    Satellite,
}
```

**`SsrSegment` wraps nothing — it holds ng's own `SsrSegment` directly.** `SsrSegment` is a port of
`ssr::types::Locus` (spec §5), born 1-based/`u64` where production's is 0-based/`u32`; production's is
left alone (spec Revision). With no coordinate mismatch to hide, the earlier `ng::SsrSegment` wrapper
has no job. It carries motif, borders and purity — **coordinates, no bases** (spec §1.2): production
embeds tract+flanks because its consumer genotypes without a FASTA open, and this module types
regions. It is `ReadPreparer::Locus` (closes `read_preparation_ssr.md` §8), and that consumer fetches
the tract from the reference it already opens.

```rust
/// ng's port of the catalog's `Locus`, minus its `ref_bytes` — 1-based inclusive, u64. Private
/// fields + a validating `new`, whose invariant is now just `1 <= start <= end` plus a finite
/// purity. (spec §4, §5, §1.2)
pub struct SsrSegment { /* chrom, start, end, motif: Motif, purity_fraction */ }

/// Ported too, though spec §4 expected `ssr::types::Motif` reused: production's is `pub(crate)`, so
/// a `pub` ng `SsrSegment` returning one trips `private_interfaces` — and widening it is a production
/// edit. Coordinate-free, so the port is 40 lines. (spec §4, corrected at the A1 review)
pub struct Motif { /* buf: [u8; MAX_MOTIF_LEN], len: u8 */ }

/// ng's classification policy — the port of `postprocess::build_loci`, windowed and all-knobs. (spec §5a)
pub fn classify(recs: Vec<RepeatInterval>, chrom: &str,
             bases: &[u8], bases_start: u64, contig_len: u64,
             p: &SsrSegmentCriteria) -> Classified;

pub struct Classified {
    pub loci: Vec<SsrSegment>,                 // what build_loci returns today
    pub bundled: Vec<RepeatInterval>,     // what it silently drops today (spec §2.4, §6a)
}
```

```rust
/// A set of genome regions to walk — sorted, non-overlapping, coalesced, clamped, genomic order.
/// Wraps production's `RegionSet` (parses BED, coalesces, clamps, converts BED's 0-based text once,
/// drops zero-length contigs). The wrapper owns the one u32 → u64 widening at this seam. (spec §2.5, §4)
pub struct GenomeRegions { /* RegionSet */ }
impl GenomeRegions {
    pub fn whole_contigs(contigs: &[ContigBounds]) -> Self;                       // the default
    pub fn from_bed_path(bed: &Path, contigs: &[ContigBounds]) -> Result<Self, BedError>;
}

/// The walk's policy. `Default` = the catalog's settings, for oracle comparability only — not an
/// endorsement of them (spec §5, §5.2). Every field is a parameter, including the period scope and
/// copy floors production hardcodes — ng's copy makes them knobs (spec §5).
pub struct TypedRegionConfig {
    pub periods: PeriodRange,
    pub scan: ScanParams,
    pub max_repeat_len: Bp,       // the satellite cap AND the detection margin — one field (spec §2.6)
    pub window_bp: Bp,            // memory knob; must not change the output (spec §2.3)
    pub criteria: SsrSegmentCriteria,  // ng's own (segment_criteria.rs), NOT ssr's CatalogParams: same
                                  // defaults + the hardcoded period scope and copy floors as real
                                  // knobs; no separate `bundle_threshold` (spec §2.4, §5)
}

/// Running tally — "no silent caps": a base typed away from the STR path is accounted for. (spec §3.1)
pub struct TypedRegionCounts {
    pub spans: u64,
    pub ssr_loci: u64,
    pub ssr_bundles: u64,   pub ssr_bundle_bp: u64,   // the number for spec §10's bundle question
    pub generic: u64,
    pub satellites: u64,    pub satellite_bp: u64,
    pub repeat_bp_with_no_locus: u64,
    pub rejected_by_reason: RejectionCounts,          // { copy_floor, purity, compound, no_clean_trim, flank_clamped }
}

/// Fatal, walk-level. `#[non_exhaustive]`. One variant, because `GenomeRegions` has already
/// rejected everything else that could go wrong (unknown contig, bad BED line, out-of-range span). (spec §8.2)
pub enum TypedRegionError { Reference(RefSeqError) }
```

## The interface — a concrete iterator, no trait

Step 3 has **no swappable-trait bake-off** (spec §6): its variation is config knobs, and its arch
doc's imagined competitor — active-region detection — is data-driven and lives inside the pileup.
So it is a concrete iterator, not a trait.

```rust
/// Walks the reference: region by region in genomic order, gap-free. Holds one ~102 kb window plus
/// three unfinished coordinates — never a contig, let alone the genome (spec §2.6, §7).
/// Concrete over the accessor (needs raw bytes + eviction, impl capabilities not trait methods) and
/// **owns its inputs**, so it can be moved onto a producer thread (spec §7).
pub struct TypedRegionIterator { /* … */ }

impl Iterator for TypedRegionIterator {
    /// Every window reads the reference, so a read can fail mid-walk. `None` meaning both "done" and
    /// "a window failed" would silently un-call the rest of the genome, so the error is in-stream.
    /// Fused: `Some(Err(_))` once, then `None`. (spec §8.2)
    type Item = Result<TypedRegion, TypedRegionError>;
}

impl TypedRegionIterator {
    /// The only constructor. Infallible — `GenomeRegions` validated everything already. `reference`
    /// and `spans` are taken **by value** (ownership, above).
    pub fn over_regions(reference: WindowedRefSeq, spans: GenomeRegions,
                        config: TypedRegionConfig) -> Self;
    /// The running tally — readable mid-walk, complete once exhausted.
    pub fn counts(&self) -> &TypedRegionCounts;
}
```

**Contract.** The items tile the requested regions exactly: contiguous, non-overlapping, complete,
and **maximal** (no two consecutive share a kind — for `Generic` this is a correctness requirement,
spec §2.3). Pure function of the reference + regions + config; `window_bp` changes memory, never the
output. The queue that fans work out to analysis workers is the pipeline's, not this iterator's
(spec §7).

## Decisions — one line each, *why* in the spec

- **`SsrSegment(SsrSegment)`, no wrapper** — ng's `SsrSegment` is born 1-based, so there is nothing to wrap
  (spec §4).
- **ng copies; production is frozen** — the classification policy, the pre-filter, `SsrSegment`, `Motif`, and
  the classification knobs are ng's own (spec Revision, §5). Only `RegionSet` is reused as-is, because
  wrapping it read-only costs production nothing. `Motif` was expected to reuse and does not: it is
  `pub(crate)`, so a `pub` ng type cannot return it (spec §4).
- **Concrete iterator, no trait** — no bake-off; `LocusRouter`/`LocusSource`/`RefWindow` retire
  (spec §6). `LocusKind` does **not** — it stays as the per-*locus* core's input (`Ssr`/`Generic`),
  which is a different type from this step's per-*region* `RegionKind` (`ng_step_interfaces.md` §3).
- **`Item = Result<_, _>`, one error variant** — fatal errors in-stream so a failure can't look like
  end-of-walk; `GenomeRegions` owns all the input validation (spec §8.2).
- **Owns its inputs, no lifetime** — so it moves onto a producer thread (spec §7).
- **`GenomeRegions` wraps `RegionSet`** — read-only reuse, the one conversion seam (spec §2.5, §4).
- `OPEN:` none at this layer — the open questions (period × length frontier, bundle disposal,
  satellite cap, flank size) are all *parameter values*, not interface shape (spec §10).

## Reconciliation with existing code

Every row read at the cited line (2026-07-16). Convergence, not new types.

| ng name | existing code | action |
|---|---|---|
| `GenomeRegion` | `regions::Region` ([regions.rs:38](../../../../src/regions.rs)), `ContigInterval` ([bam/alignment_input.rs:544](../../../../src/bam/alignment_input.rs)) | the consolidation `ng_step_interfaces.md` §6 reserved — this step is its first use (1-based, `u64`) |
| `GenomeRegions` | `RegionSet` ([regions.rs:74](../../../../src/regions.rs)) | **wrap read-only**, not reimplement — it already sorts/coalesces/clamps/converts-BED; wrapper adds ng's width + base + name. `regions.rs` does not move |
| `SsrSegment`'s payload | `ssr::types::Locus` ([ssr/types.rs:136](../../../../src/ssr/types.rs)) | **copy into ng**, born 1-based/`u64` (spec §4, §5); production's untouched |
| ng `SsrSegment`'s motif | `ssr::types::Motif` ([ssr/types.rs:36](../../../../src/ssr/types.rs)) | **copy into ng** — no coordinates/width to rebase, but production's is `pub(crate)` and leaks through ng's `pub` `SsrSegment` (spec §4, corrected at the A1 review) |
| `RegionKind::SsrBundle` tracts | `RepeatInterval` ([ng/tandem_repeat.rs:186](../../../../src/ng/tandem_repeat.rs)) | reuse; widened to `u64` (ng's own code, spec §4) |
| `TypedRegionConfig.criteria` | `CatalogParams` ([ssr/catalog/mod.rs:42](../../../../src/ssr/catalog/mod.rs)) | **copy into ng** as `SsrSegmentCriteria`: same defaults, no `bundle_threshold`, gains the hardcoded period/floor knobs (spec §5) |
| `segment_criteria::classify` | `build_loci` ([ssr/catalog/postprocess.rs:69](../../../../src/ssr/catalog/postprocess.rs)) | **copy into ng** — logic transcribed unchanged; takes `RepeatInterval`, windowed, returns `Classified { loci, bundled }` (spec §5a). Production's `build_loci` and `TrfRecord` stay |
| the pre-filter | `catalog_prefilter` ([ng/scanner_parity.rs](../../../../src/ng/scanner_parity.rs), moved out of `src/ssr/` at B2) | **copy into ng**, beside `classify` (spec §5.1) |
| the windowed scan | `collect_windowed` ([ng/tandem_repeat.rs:530](../../../../src/ng/tandem_repeat.rs)) | promote (private → crate) and stream (spec §6.1) |
| reference access | `WindowedRefSeq` ([ng/ref_seq.rs](../../../../src/ng/ref_seq.rs)) | + a raw-bytes path (`ref_seq.md` YAGNI now spent) + a `contigs()` accessor |
| `TypedRegionError` | new (`thiserror`, `#[non_exhaustive]`) | one variant; mirrors `ReadFilterError`'s in-stream-fatal shape |

**Two oracles, both in spec §8:**

- **Port fidelity** (spec §8.0) — ng's `classify` vs production's `build_loci` on the same intervals,
  whole-contig, at the catalog's settings; `SsrSegment` sets identical modulo the coordinate base. The
  bridge is `TrfRecord::for_test` (`#[cfg(test)] pub(crate)`, same crate). This is what the Revision
  bought: what was once true by construction is now a test.
- **`.cat` parity** (spec §8.1) — the committed **trf-mod-built** golden catalog on the same
  reference; the `SsrSegment` set is a strict subset (differs only by the satellite cap), with
  `scanner_parity`'s overlap tolerance for detector wobble.

## Test & bench shape

Tests beside the code in `region_typing/` — the walk's in `mod.rs`, the port's (including spec §8.0's
differential) in `segment_criteria.rs`. No `bench/` — no competing impls (spec §6). The regression anchors
are the two oracles above and the invariant tests (partition, window-invariance, BED-invariance)
detailed in spec §8; the arch doc adds no test surface of its own.
