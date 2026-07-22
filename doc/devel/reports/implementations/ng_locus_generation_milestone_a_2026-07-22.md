# ng locus generation — Milestone A (the shared types), 2026-07-22

Plan-driven ([impl plan](../../ng/impl_plan/locus_generation.md), spec
[locus_generation.md](../../ng/spec/locus_generation.md), arch
[locus_generation.md](../../ng/arch/locus_generation.md)); one implement → review → apply →
commit loop per step. Milestone A of the **shared locus-generation shape** — the locus type every
generator will produce, and the one derived method on it. No generator, no dispatcher, no iterator
yet (Milestones B, C).

## Plan

The types-first milestone: scaffold the module, define the locus types, then the derived
per-position depth. A1 (pure scaffold) paired with A2 in one commit per the plan's granularity note.

## Assumptions (silent choices surfaced)

- **`LocusKind` / `SsrDetail` are declared in the shared `locus_generation/mod.rs`**, not in
  `ssr.rs`, because `SampleLocusObservations.kind` names `LocusKind` — the type graph forces the
  declaration into the shared module. The STR generator (a later plan) only *fills* `SsrDetail`.
  This is the placement the arch doc left slightly loose (arch §1); resolved the only clean way, no
  new design.
- **Over-reach of a partial's covered extent (`n > region.len()`) is clamped, not asserted.** It is
  a producer invariant (enforced where `ReadCoverage` is minted, the STR generator); the consumer
  clamps defensively so the derivation is panic-free on any input. No `debug_assert` — this repo has
  recorded being burned twice by debug-only guards compiled out of release (PROJECT_STATUS, the
  step-3 A3 review and the CLI E-review). A deliberate deviation from the review's "add a
  debug_assert" suggestion, with the reason recorded at the call site.

## Changes

- New [`src/ng/locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs): the folder
  module (a step folder — implementations compete), wired into
  [`src/ng/mod.rs`](../../../../src/ng/mod.rs).
- Types: `SampleLocusObservations` (owned, no lifetimes), `ObservedSequence`, `ReadCoverage`,
  `LocusKind` (`#[non_exhaustive]`, derives `Eq`), `SsrDetail`.
- `impl SampleLocusObservations`: `num_obs_along_locus()` (derived per-position observation depth)
  and `complete_observations()` (the complete-only guard).
- Reuse: `ChainId` (`pileup_record`), `Motif` (`region_typing`), `GenomeRegion` (`ng::types`); the
  moments modelled on `AlleleObservation` / `AlleleSupportStats`, minus the anchor-relative
  read-position-bias fields (the pileup generator's).

## Tests (9)

Build a locus of each `LocusKind` incl the zero-coverage SSR case; complete≠partial of identical
bases (the dedup-key property); depth derivation over a mixed complete + left-partial + right-partial
fixture whose expected vector `[5,5,5,5,3,3,3,8,8,8]` is asymmetric so a left↔right transposition
cannot alias; single-base locus; empty (zero-coverage) locus; both-ends over-reach clamp;
complete-only filtering pinned to yield *all* completes.

## Validation

Container (`./scripts/dev.sh`): `cargo fmt --check` clean, `cargo clippy --lib -- -D warnings`
clean, `cargo test --lib ng::locus_generation` — 7 module tests pass (plus the 2 A2 smoke tests).

## Tradeoffs / follow-ups

- Carried to **B/A-later**: the observation dedup (when the STR generator folds observations) must
  use an explicit `(bases, read_coverage)` key — never derived `==`, which compares every field.
- `num_obs_along_locus()` is *observation* depth: it omits covering-but-unobserved reads and
  double-counts across overlapping generic loci. Both caveats are the paralog filter's (spec §11),
  documented on the method.
- Next: Milestone B — the `LocusGenerator<S>` contract, `NoLoci`, and the dispatcher (the heart,
  proven with `NoLoci` alone).
