# ng locus generation — Milestones B & C (the contract, dispatcher, iterator), 2026-07-22

Plan-driven ([impl plan](../../ng/impl_plan/locus_generation.md), spec
[locus_generation.md](../../ng/spec/locus_generation.md), arch
[locus_generation.md](../../ng/arch/locus_generation.md)); one implement → review → apply → commit
loop per step. Milestones B and C complete **the shared locus-generation shape** — the contract
every generator plugs into, the dispatcher that routes and accounts, and the public iterator. Ships
**no real generator**: every kind resolves to `NoLoci`, which is exactly what the shape spec ships.

## Plan

B (the heart, proven with `NoLoci` alone): the `LocusGenerator<S>` contract, `NoLoci`, and the
dispatcher + accounting. C (the plumbing): the iterator wrapping it, and the acceptance anchor.

## Assumptions / deviations surfaced (recorded, not silent)

- **`LocusGenerationError` was pulled forward from C1 into B1** — the trait's `next_locus` signature
  names it, so it cannot land later. C1 is therefore a no-op.
- **`GeneratorSlot` uses `Box<dyn LocusGenerator<S>>` trait objects**, per `ng_step_interfaces.md`
  §4's "`Box<dyn _>` for the lab" choice. This collapses the arch's iterator `<T, G>` (generic over
  the generator set) into a **concrete `GeneratorSet`** — the per-run swap lives in the boxes, not a
  type parameter. Recorded; the arch §4 / spec §6 surface is left as the proposal (design docs are
  not edited during implementation).
- **`LocusConfig` was deferred, not shipped as an empty struct.** An empty config the iterator
  stores and never reads is a *dormant lever*, which the repo's "no dormant knobs" discipline
  (read_filtering) forbids. `new(regions, reads, generators)` drops the arch's `config` param; it
  lands when it has a field.
- **The over-reach clamp in `num_obs_along_locus` (A3) is documented, not `debug_assert`ed** — this
  repo has recorded debug-only guards biting twice (they compile out of release).

## Changes

All in [`src/ng/locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs):

- **B1:** `LocusGenerator<S>` (`begin_segment` / `next_locus`, `S` a trait parameter);
  `NoLoci { reason }` + `UnhandledReason`; `impl<S> LocusGenerator<S> for NoLoci`;
  `LocusGenerationError` (`#[non_exhaustive]` thiserror, wraps `TypedRegionError` / `IngestError` /
  `RefSeqError`).
- **B2:** `LocusCounts` (regions_in, **regions_handled**, loci_emitted, the two unhandled counters +
  `_bp`) with `record_unhandled` as the one chokepoint; `GeneratorSlot<S>` (Generator / Unfilled);
  `GeneratorSet` — the stateful dispatcher (`begin_region` / `next_locus` / `counts` /
  `all_unimplemented` / `new`), one region resident, no locus buffer. `Satellite` has no slot.
- **C2:** `SampleLocusObservationsIterator<T>` — `new` / `counts` / `Iterator` / `FusedIterator`.

## Tests (16 in the module)

Contract: `NoLoci` emits nothing and carries its reason. Dispatch: every kind accounted with the two
nothings separate; **each kind routes to its own slot** (distinguishable generators 2/3/5, so a
mis-route shows a wrong count); a filled slot's region is *handled*, not unhandled. Iterator: a
multi-kind stream drains and accounts; loci in coordinate order; **a region yielding several loci
streams them all before advancing** (spec §13.4); both the stream-error and the generator-error
paths are fatal and fuse. Plus the A-milestone locus-type and depth tests.

## Validation

Container (`./scripts/dev.sh`): `cargo fmt --check` clean; `cargo clippy --lib --tests --all-features
-- -D warnings` clean; `cargo test --lib` — **2109 passed, 0 failed**. **Native host** (macOS):
`cargo test --lib ng::locus_generation` — 16 passed. (A pre-existing example, `dhat_ng_merge`, has
unrelated `-D warnings` issues that predate this work and do not touch `locus_generation`.)

## Tradeoffs / follow-ups

- The **real-generator `next_locus` dispatch path** is written and unit-tested via fakes
  (`FixedCountGenerator` / `EchoGenerator`), but no *production* generator exercises it until the STR
  generator lands (its plan is blocked on the ng STR aligner).
- **`GeneratorSet` payload types** for `Generic` and `SsrBundle` are `()` for now — they refine when
  those generators are designed (the arch's deferred "confirm when the second generator lands").
- `Box<dyn LocusGenerator<S>>` carries **no `Send` bound** (single-threaded v1); documented as a
  conscious omission for the deferred producer-thread parallelism (arch §9).
- **The shared shape is complete.** Next locus-generation work is a real generator, each its own
  plan.
