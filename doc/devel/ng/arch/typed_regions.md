# ng step 3 — typed-region generator: types & interfaces

*Architecture draft (2026-07-16), the code-facing companion to
[`../spec/typed_regions.md`](../spec/typed_regions.md). It gives the types and signatures as they
appear in code; **every *why* points back to the spec** — this doc does not re-argue them. It
supersedes the step-3 sketch (`LocusRouter`/`LocusSource`) in the shared
[`ng_step_interfaces.md`](ng_step_interfaces.md) §3, and the `locus_router/` entry in
[`module_layout.md`](module_layout.md); the edits to those two docs land alongside this one.
Naming: **STR** in prose, `ssr` in code.*

---

## Module home

`src/ng/region_typing.rs` — a single file (a step with no bake-off is a file, not a folder; spec
§6). It is a **step** module, so it sits in the step layer of the `src/ng/` tree, replacing the
`locus_router/` folder the layout doc reserved. `GenomeRegion` and `Position` are shared vocabulary
and land in `ng::types`, not here.

Promoted to `region_typing/` only if it outgrows one file — the `tandem_repeat.rs` precedent.

## The types

The final shapes. Contracts are stated below each; the *why* is the cross-ref.

```rust
/// A genome region plus what the sequence there IS. The walk's output.
/// `region` is a field, not a per-variant repeat: every typed region has one, structurally, and it
/// is the one place ng's 1-based base is stated. (spec §3, §4)
pub struct TypedRegion { pub region: GenomeRegion, pub kind: RegionKind }

/// Exactly one kind is a genetic object (a locus); the other three are physical regions.
/// `Generic`/`Satellite` carry nothing — they *are* just their region. (spec §1.1)
pub enum RegionKind {
    SsrLocus(Locus),                               // ssr::types::Locus, used directly (see below)
    SsrBundle { tracts: Box<[RepeatInterval]> },   // >= 2, coordinate-ordered; hull = TypedRegion.region
    Generic,
    Satellite,
}
```

**`SsrLocus` wraps nothing — it holds the catalog's `Locus` directly.** This is the payoff of
rebasing `Locus` to 1-based (spec §4, §9): with no coordinate mismatch to hide, the earlier
`ng::SsrLocus` wrapper has no job. `Locus` carries motif, borders, purity, and the embedded
flank+tract+flank bases, and is `ReadPreparer::Locus` (closes `read_preparation_ssr.md` §8).

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
/// copy floors that are hardcoded consts today (spec §9 moves them into `CatalogParams`).
pub struct TypedRegionConfig {
    pub periods: PeriodRange,
    pub scan: ScanParams,
    pub max_repeat_len: Bp,       // the satellite cap AND the detection margin — one field (spec §2.6)
    pub window_bp: Bp,            // memory knob; must not change the output (spec §2.3)
    pub catalog: CatalogParams,   // all admission knobs; `flank_bp` is now also the bundle threshold
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

- **`SsrLocus(Locus)`, no wrapper** — `Locus` is rebased to 1-based, so there is nothing to wrap
  (spec §4).
- **Concrete iterator, no trait** — no bake-off; `LocusRouter`/`LocusSource`/`RefWindow` retire
  (spec §6). `LocusKind` does **not** — it stays as the per-*locus* core's input (`Ssr`/`Generic`),
  which is a different type from this step's per-*region* `RegionKind` (`ng_step_interfaces.md` §3).
- **`Item = Result<_, _>`, one error variant** — fatal errors in-stream so a failure can't look like
  end-of-walk; `GenomeRegions` owns all the input validation (spec §8.2).
- **Owns its inputs, no lifetime** — so it moves onto a producer thread (spec §7).
- **`GenomeRegions` wraps `RegionSet`; `SsrLocus` holds `Locus`; `TypedRegionConfig.catalog` is
  `CatalogParams`** — reuse production types directly, widened/rebased rather than duplicated (recon
  below).
- `OPEN:` none at this layer — the open questions (period × length frontier, bundle disposal,
  satellite cap, flank size) are all *parameter values*, not interface shape (spec §10).

## Reconciliation with existing code

Every row read at the cited line (2026-07-16). Convergence, not new types.

| ng name | existing code | action |
|---|---|---|
| `GenomeRegion` | `regions::Region` ([regions.rs:38](../../../../src/regions.rs)), `ContigInterval` ([bam/alignment_input.rs:544](../../../../src/bam/alignment_input.rs)) | the consolidation `ng_step_interfaces.md` §6 reserved — this step is its first use (1-based, `u64`) |
| `GenomeRegions` | `RegionSet` ([regions.rs:74](../../../../src/regions.rs)) | **wrap**, not reimplement — it already sorts/coalesces/clamps/converts-BED; wrapper adds ng's width + name |
| `SsrLocus`'s payload | `ssr::types::Locus` ([ssr/types.rs:136](../../../../src/ssr/types.rs)) | **use directly**, rebased to 1-based + widened to `u64` (spec §9); no wrapper |
| `RegionKind::SsrBundle` tracts | `RepeatInterval` ([ng/tandem_repeat.rs:186](../../../../src/ng/tandem_repeat.rs)) | reuse as-is |
| `TypedRegionConfig.catalog` | `CatalogParams` ([ssr/catalog/mod.rs:42](../../../../src/ssr/catalog/mod.rs)) | reuse; loses `bundle_threshold`, gains the hardcoded period/floor knobs (spec §5, §9) |
| the admission call | `build_loci` ([ssr/catalog/postprocess.rs:69](../../../../src/ssr/catalog/postprocess.rs)) | windowed + returns `Admitted { loci, bundled }`; `#[cfg(test)]` bridge must open (spec §5) |
| the windowed scan | `collect_windowed` ([ng/tandem_repeat.rs:530](../../../../src/ng/tandem_repeat.rs)) | promote (private → crate) and stream (spec §6.1) |
| reference access | `WindowedRefSeq` ([ng/ref_seq.rs](../../../../src/ng/ref_seq.rs)) | + a raw-bytes path (`ref_seq.md` YAGNI now spent) + a `contigs()` accessor |
| `TypedRegionError` | new (`thiserror`, `#[non_exhaustive]`) | one variant; mirrors `ReadFilterError`'s in-stream-fatal shape |

**Parity oracle:** an `ssr-catalog` `.cat` on the same reference, at the catalog's settings — the
`SsrLocus` set is a strict subset (differs only by the satellite cap). Full contract in spec §8.

## Test & bench shape

Tests beside the code in `region_typing.rs`. No `bench/` — no competing impls (spec §6). The
regression anchors are the parity oracle and the invariant tests (partition, window-invariance,
BED-invariance) detailed in spec §8; the arch doc adds no test surface of its own.
