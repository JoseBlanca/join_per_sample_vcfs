# ng locus generation (the shared shape) — implementation plan

**Status:** draft, 2026-07-22. The build order for the **shared locus-generation shape**: the
`locus_generation/` module, the locus type every generator produces
(`SampleLocusObservations`), the `LocusGenerator<S>` contract, the count-only `NoLoci` generator,
the dispatcher, and the `SampleLocusObservationsIterator` public surface. Design is settled in
[`../spec/locus_generation.md`](../spec/locus_generation.md) (spec) and
[`../arch/locus_generation.md`](../arch/locus_generation.md) (types & interfaces), under the shared
arch docs ([step interfaces](../arch/ng_step_interfaces.md),
[module layout](../arch/module_layout.md)). This turns that design into build order; it is **not**
a place for new design — the shape's only remaining open item is a *type choice* deferred to the
pileup generator (spec §12), not an interface question here.

This plan ships **no real generator** — every region kind routes to `NoLoci`, exactly what the
shape spec ships (spec §2). The first real generator (STR) is
[`locus_generation_ssr.md`](locus_generation_ssr.md), its own plan.

---

## Scope

**In:** `src/ng/locus_generation/mod.rs`; the locus types (`SampleLocusObservations`,
`ObservedSequence`, `ReadCoverage`) and the shared `LocusKind` / `SsrDetail`; the derived
`num_obs_along_locus()` / `complete_observations()`; the `LocusGenerator<S>` trait; `NoLoci` +
`UnhandledReason`; the dispatcher over `RegionKind`; `LocusConfig` / `LocusCounts`;
`LocusGenerationError`; the `SampleLocusObservationsIterator`.

**Out (later plans):**

- **The STR locus generator** — `SsrSegment` → one locus (fetch + align + tally). Specced in
  [`../spec/locus_generation_ssr.md`](../spec/locus_generation_ssr.md); its own plan.
- **The pileup / generic generator** — `Generic` region → many loci from the data
  (`locus_generation/pileup/`, spec §11); its own spec + plan.
- **The `SsrBundle` generator** — depth-per-position at least (spec §11). Routed to
  `NoLoci { NotImplemented }`, counted, meanwhile.
- **The cohort merge, the artifact, and the paralog windowed statistics** — downstream of a
  sample's loci (spec §3, §11).
- **Parallelism** — deferred whole (spec §9).

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The dispatch + accounting is the heart; prove it with `NoLoci` before the iterator plumbs
  it.** What this shape must get right is "every kind assigned, every base accounted for, the two
  kinds of nothing kept apart" — all provable with the count-only generator alone (spec §13), so
  Milestone B builds and proves it before Milestone C wraps it in a stream.
- **Reuse over rewrite.** `ObservedSequence` mirrors production's per-allele shape
  (`AlleleObservation` + `AlleleSupportStats`), `ChainId` is reused as-is, the iterator mirrors
  `TypedRegionIterator`'s lazy/counts/fatal-in-stream shape, and the error wraps upstream errors —
  nothing is re-derived (arch §8).
- **Verify against ground truth.** The north-star is the shape's acceptance test (spec §13):
  exhaustive dispatch, reconciled accounting, preserved order — asserted over a multi-kind
  fixture, not self-consistency.
- **Incremental, with pauses.** One milestone, stop for review, next.
- **Ungated / container builds.** `ng` stays a plain module; all `cargo` via `./scripts/dev.sh`
  (CLAUDE.md); native host build at completion. **No `bench/`** — the dispatcher has no bake-off
  (the generators do), so its measure is the fixture test, not a frontier (arch §Test & bench).

## Preconditions (already in place)

- **Step 3 is complete** — `region_typing/` gives `TypedRegion`, `RegionKind` (`SsrSegment` /
  `SsrBundle` / `Generic` / `Satellite`), `TypedRegionIterator`, `TypedRegionCounts`,
  `TypedRegionError` ([region_typing/mod.rs](../../../../src/ng/region_typing/mod.rs); typed_regions
  plan Checkpoint E). `Motif` and `SsrSegment` live in `region_typing/segment_criteria.rs`.
- **Read ingestion is merged into this branch** (via `main`, commit 3ec3f97) — `SampleReads`,
  `reads_in_region`, `IngestError` ([read/input/mod.rs](../../../../src/ng/read/input/mod.rs)),
  `MappedRead` ([bam/alignment_input.rs](../../../../src/bam/alignment_input.rs)).
- **Foundations are done** — `ng::types` (`GenomeRegion`, `Position`, `ContigId`, `Bp`),
  `RefSeqError`, and the `RefSeq` / `RawRefSeq` / `EvictableRefSeq` accessors.
- **The reuse target for `ObservedSequence`** — `AlleleObservation` + `AlleleSupportStats` +
  `ChainId` ([pileup_record.rs](../../../../src/pileup_record.rs)), and `CohortPileupRecord`
  ([var_calling/types.rs](../../../../src/var_calling/types.rs)) as the cohort-composition model —
  are readable (arch §8).

---

## The steps

### Milestone A — the shared types (types, then the one derived method)

**✅ A1. Scaffold the `locus_generation/` module.** *(landed with A2, one commit.)*
`src/ng/locus_generation/mod.rs` with its `#[cfg(test)]` block; wire `pub mod locus_generation;`
into `ng/mod.rs`. A **folder**, not a file — a step whose implementations compete (module_layout
principle 1); `ssr.rs` and `pileup/` land in later plans. *Source:* arch §Module home.

**✅ A2. The locus types.**
`SampleLocusObservations` (owned, no lifetimes), `ObservedSequence`, `ReadCoverage`
(`Complete` / `PartialLeft(u16)` / `PartialRight(u16)`), and the shared `LocusKind`
(`Generic` / `Ssr(SsrDetail)` / `SsrBundle`, `#[non_exhaustive]`) + `SsrDetail` (`motif: Motif`
from `region_typing`, `left_flank` / `right_flank`). Declared here in the shared module because
`SampleLocusObservations.kind` names `LocusKind`; the STR generator *fills* `SsrDetail` (plan 2).
No logic. *Depends:* A1. *Source:* spec §3, arch §1; ng_step_interfaces §2.

**✅ A3. `num_obs_along_locus()` + `complete_observations()`.**
The derived per-position depth (a `Complete` counts at every position, a `PartialLeft(n)` at the
leftmost `n`, `PartialRight(n)` at the rightmost `n`; length = `region.len()`) and the
complete-only iterator guard. **Own step** — an off-by-one in the partial-coverage arithmetic is a
silently-wrong depth, not a crash, so it lands with its own unit tests green (spec §13 tests 5, 6):
depth derives correctly over a fixture mixing complete and partial observations, and a `Complete`
vs `PartialLeft` of identical `bases` stay distinct entries. *Depends:* A2. *Source:* spec §3, §13.

> **Checkpoint A:** the locus types compile; the depth derivation and the complete/partial
> distinction are unit-tested. Pause for review.

### Milestone B — the contract, `NoLoci`, and the dispatcher (the heart)

**✅ B1. `LocusGenerator<S>` + `NoLoci` + `UnhandledReason`.** *(also `LocusGenerationError`, pulled forward from C1 — the trait signature needs it.)*
The trait (`begin_segment`, `next_locus`, `S` a trait parameter not an associated type — arch §2);
`NoLoci { reason: UnhandledReason }` with `UnhandledReason::{OutOfScope, NotImplemented}`, and
`impl<S> LocusGenerator<S> for NoLoci` (counts the segment and its bases, emits nothing). *Depends:*
A. *Source:* spec §4, §5, arch §2, §3. Verified: `NoLoci` over a segment yields `None` and tallies
the bases; the reason is carried, not fixed.

**✅ B2. The dispatcher + `LocusCounts`.** *(GeneratorSet with trait-object slots collapses the arch `<T,G>` to concrete; `LocusConfig` deferred — an empty struct is a dormant lever.)*
The exhaustive `match` over `RegionKind` routing each branch's payload to the generator supplied
for it, with `Satellite → NoLoci { OutOfScope }` and every unfilled kind → `NoLoci { NotImplemented }`
(spec §5); `LocusCounts` (regions/loci + the two `unhandled_*` counters and their `_bp`);
`LocusConfig` (minimal — dispatcher-only knobs, generators own theirs, spec §7); the concrete
`GeneratorSet` shape (one slot per kind — the impl-time confirmation arch §3 leaves open, pinned
here, not new design). *Depends:* B1. *Source:* spec §5, §7, arch §3, §5. Verified: the `match` is
exhaustive (compiler); **each kind reaches the generator it should** and the two kinds of nothing
land in **different** counters (spec §13 tests 1, 3), over a fixture mixing a satellite and a
generic region.

> **Checkpoint B:** the dispatch is exhaustive and every kind is accounted for — the shape's
> correctness, proven with `NoLoci` alone, before any stream wraps it. Pause for review.

### Milestone C — the public iterator (plumbing) + the acceptance anchor

**✅ C1. `LocusGenerationError`.** *(landed in B1.)*
`#[non_exhaustive]` wrapping `TypedRegionError` / `IngestError` / `RefSeqError` — this step mints
no error of its own (arch §4). *Depends:* B. *Source:* spec §6, arch §4.

**✅ C2. `SampleLocusObservationsIterator<T>`.** *(concrete GeneratorSet; no config param.)*
`new(regions, reads, generators, config)`; `Iterator` with
`Item = Result<SampleLocusObservations, LocusGenerationError>` — pull a region, `begin_segment`,
call `next_locus` until `None`, then the next region; **no buffer of loci**; `counts()` running
tally; fused, fatal-in-stream (`Some(Err(_))` once then `None`). Mirrors `TypedRegionIterator`.
*Depends:* C1. *Source:* spec §6, arch §4. Verified: iterate a small multi-region fixture; a
seeded upstream error surfaces as `Some(Err(_))` once then `None`.

**✅ C3. The acceptance anchor.**
Over a multi-region, multi-kind fixture (all routed to `NoLoci`): `regions_in` reconciles against
loci emitted plus `unhandled_not_implemented` plus `unhandled_out_of_scope`, and the two `_bp`
counters sum to the bases those regions cover; loci come out in **coordinate order** across kinds;
the two kinds of nothing stay separate. The shape's north-star (spec §13 tests 1–4). *Depends:* C2.
*Source:* spec §13.

> **Checkpoint C:** locus generation runs end-to-end over a typed-region stream; the accounting
> reconciles and order is preserved. **The shared shape is complete** (every real generator plugs
> in via a later plan). Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | type-level tests; **depth derivation** from mixed complete/partial `read_coverage`, and complete≠partial identity (spec §13.5, §13.6) |
| B | exhaustive `match` (compiler); **each kind routes** + **two kinds of nothing separate**, via `NoLoci` (spec §13.1, §13.3) |
| C | iterator + running counts + fatal-in-stream; **acceptance fixture** — full accounting reconciles and **coordinate order** holds across a multi-kind stream (spec §13.1–13.4) |

## Out of scope (next plans)

- **The STR locus generator** — [`locus_generation_ssr.md`](locus_generation_ssr.md).
- **The pileup / generic generator** and **the `SsrBundle` generator** — `locus_generation/`
  siblings, their own specs + plans (spec §11).
- **The cohort merge, the artifact, the paralog windowed statistics, parallelism** — downstream of
  a sample's loci (spec §3, §9, §11).
