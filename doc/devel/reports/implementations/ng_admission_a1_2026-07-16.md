# ng typed-region generator — Milestone A1: the admission port

**Date:** 2026-07-16 · **Plan:** [typed_regions.md](../../ng/impl_plan/typed_regions.md) (Milestone A1)
· **Design:** [spec](../../ng/spec/typed_regions.md) §4, §5, §5.1, §8.0 · [arch](../../ng/arch/typed_regions.md)

## What landed

`src/ng/region_typing/` — ng's own copy of the STR admission policy, the rules deciding which
detected tandem repeats become genotypeable STR loci.

- `mod.rs` — module scaffold (the walk's own types land in C2).
- `admission.rs` — `Motif`, `Locus`, `SsrAdmissionParams` (+ `DEFAULT_*` consts), `prefilter`,
  `admit`, and the ported helpers (`upper`, `count_motif`, `is_compound`, `minimal_trim`,
  `recompute_purity`, `is_close`, `drop_bundles`), plus 36 tests including the spec §8.0 differential.
- `src/ng/mod.rs` — wires it in, and records the frozen-production rule.

**Source:** `ssr::catalog::postprocess::build_loci` (itself a port of GangSTR's `minimal_trim.py` /
`remove_bundles.py`), plus `catalog_prefilter` from `ssr::catalog::scanner_parity`. The **logic is
transcribed unchanged**; only the shape is ng's.

## Why a copy at all

Owner decision, 2026-07-16 (spec Revision): *"don't touch `src/ssr` … copy the code to `ng/` and
modify there … leave production as is and create a fresh ng caller from scratch."* ng needs admission
windowed, 1-based/`u64`, `RepeatInterval`-driven, all-knobs, and handing bundle members back — five
changes to one function, which is not "a small tweak", which is exactly when the rule says copy.
Production stays an **independent yardstick** for the experiments ng exists to run.

**This plan touches no file outside `src/ng/`.** Verified by `git status` at every commit.

## Divergences from production (the whole list)

1. **1-based inclusive** coordinates (production: 0-based half-open) — spec §4.
2. **`u64`** (production: `u32`) — spec §4.
3. **Input is `RepeatInterval`** via an internal `Candidate`, not trf-mod's `TrfRecord`. `build_loci`
   only ever reads `start`/`end`/`period`/`score` — exactly `RepeatInterval`'s fields. **ng does not
   depend on trf-mod.**
4. **The pre-filter lives beside the policy** it is inseparable from (it exists *because* the bundle
   drop runs before the copy floor), not in a test file — spec §5.1.
5. **`Motif` is ported too** — see deviations.

The rebase is applied at **exactly one site** (`finish_locus`'s `Locus::new` call): `[s, e)` → `[s+1,
e]`, so `start` shifts by one, `end` does not move, the length is unchanged, and every interior
computation stays in production's 0-based arithmetic. That is what makes the port diffable against its
source, and three review agents confirmed it independently by mechanical diff.

## How it is verified

**Spec §8.0's differential** — `admit` and production's real `build_loci` driven from the same
intervals, compared field-by-field modulo the coordinate base. Available *without touching
production* because `TrfRecord::for_test` is `#[cfg(test)] pub(crate)` in the same crate. Two cases:
production's own seven `build_loci` tests replayed, and the scanner's output over the synthetic
reference (**17 real loci across 2 contigs**; 534 raw intervals → 19 after the pre-filter on ctg1 —
confirmed non-vacuous by an instrumented probe).

This is a **stronger** check than the plan's original Milestone A had. That plan proved a production
rebase by "the golden fixture is still green" — one implementation against itself. This compares two
implementations directly.

**Its limit, learned at review and now documented in the module:** both sides run the same inputs
through transcribed logic, so an input class no test supplies is not merely uncovered — it is
*invisible*, because both sides are wrong together and the test stays green. That is why the suite
deliberately drives the gates a single fixed configuration never fires. **Every coverage claim is
mutation-verified**; see the [fixes report](../reviews/fixes_applied_2026-07-16_ng_admission_a1.md).

## Deviations from the plan (recorded, not escalated)

- **`Motif` is ported, not reused.** Spec §4 said reuse it — "no coordinates, no width, nothing to
  rebase, the Revision's 'costs production nothing' case exactly". The premise is true; **the
  conclusion was wrong**. `ssr::types::Motif` is `pub(crate)`, and ng's `pub` `Locus::motif()` returns
  one, tripping `private_interfaces`. The escapes were: widen it in production (**forbidden**), demote
  ng's whole admission surface to `pub(crate)` (bends the ng-sibling convention and buys `dead_code`
  warnings until D), or port 40 coordinate-free lines. Ported. **Dividend:** `region_typing/` now
  names nothing from `src/ssr/` outside its `#[cfg(test)]` differential. Spec §4, the arch, and the
  plan are corrected in this commit. *General lesson: "costs production nothing" must be checked
  against **visibility**, not just coordinates and width.*
- **`region_typing/` is a folder**, where the arch specified a file — the port is a second concern
  with its own test suite. Not a bake-off; `module_layout.md` principle 1 is unbent (spec §6 + arch
  updated).
- **`Candidate`** (internal) stands in for `TrfRecord`, so the widening from `RepeatInterval`'s `u32`
  happens once at the entry rather than scattering `u64::from` through the policy.

## Two documented claims this milestone disproved

Both were in the spec I wrote hours earlier; both are corrected in code and spec.

1. **"At `Default`, ng ships no score gate at all"** (§5c). The test written to prove it failed: the
   gate is `score >= 0`, so `i32::MIN` *is* rejected. True claim: Ruzzo–Tompa emits only
   positive-scoring segments, so the floor is unreachable **for scanner output** — a no-op in
   practice, not an absent gate.
2. **"The pre-filter has a second, *disagreeing* copy-floor table"** (§5/§10). They agree on every
   reachable period; the only difference (period 1) is unreachable behind two gates. The duplication
   is structural, not behavioural — so **A2 deletes a copy rather than resolving a conflict**, a
   smaller job than the spec described.

## Validation (container)

`cargo fmt --check` clean · `cargo clippy --lib --tests` clean · `ng::region_typing` **36 passed** ·
`cargo test --lib` **1802 passed, 0 failed** · `cargo doc --no-deps` 0 errors in scope.

Pre-existing and out of scope (all reproduce on a stashed clean tree): `examples/ssr_psp_seqdump.rs:41`
clippy error; `benches/psp_writer_perf.rs:386` panic; 8 `cargo doc` unresolved links elsewhere.

## Next

**A2** — hoist `MIN_PERIOD`/`MAX_PERIOD`/`copy_number_floor` into `SsrAdmissionParams`, fold the
pre-filter's second copy-floor table into it, collapse `bundle_threshold` into `flank_bp`. The
differential must stay green at `Default`. **A3** — windowing (`bases_start` + `contig_len`) and
`Admitted { loci, bundled }`. Then **Checkpoint A** (a planned pause for review).

Carried into A2/A3 from the review: the `Vec`-by-value signature, the `PrefilteredIntervals` newtype,
the test-fixture case table, and a `proptest` property over `assert_agrees` (recommended as its own
step — the strongest remaining lever on the port).
