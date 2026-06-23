# Code Review: ssr-call genotyping+pre-pass — Milestone A
**Date:** 2026-06-23
**Reviewer:** rust-code-review skill (orchestrator, focused inline pass)
**Scope:** Milestone A diff (commits `2229b16` A1 types, `4e74b78` A2 simulator)
**Status:** Approve-with-changes

---

## 1. Scope
- **Reviewed:** the Milestone A diff on branch `ssr-cohort`.
- **In-scope files:** [param_estimation.rs](../../../../src/ssr/cohort/param_estimation.rs),
  [candidate_set.rs](../../../../src/ssr/cohort/candidate_set.rs),
  [rung_ladder.rs](../../../../src/ssr/cohort/rung_ladder.rs),
  [stutter.rs](../../../../src/ssr/cohort/stutter.rs),
  [sim.rs](../../../../src/ssr/cohort/sim.rs), [mod.rs](../../../../src/ssr/cohort/mod.rs).
- **Out of scope:** the reading layer (`types`/`reader`/`merge`/`driver`) — unchanged;
  the bench harnesses.
- **Categories applied:** reliability, errors, naming, idiomatic, refactor_safety,
  module_structure, smells, extras (determinism/stable-output). Skipped:
  unsafe_concurrency (no `unsafe`/threads/atomics), tooling (no `Cargo.toml` change).

## 2. Verdict
**Approve-with-changes.** Foundational types + a deterministic simulator, cleanly
factored and well-tested. No Blockers, no Majors. The findings are one defensive
hardening of the determinism primitive and a handful of test/naming refinements.

## 3. Execution status
- `cargo fmt --check` → pass.
- `cargo clippy --all-targets --all-features -- -D warnings` → pass (clean).
- `cargo test --all-features` → **1181 lib pass**, integration + doctests green.
- `cargo test --all-targets` → pre-existing unrelated failure in
  `benches/psp_writer_perf.rs:386` (index OOB; untouched by this diff).
- Findings "Needs verification": 0.

## 4. Open questions and assumptions
1. None blocking. The `FixedPointAccum` input-magnitude bound (Mi1) is documented in
   the parameters arch (`reads × O(few nats)`), so the guard is defensive, not a
   correctness fix yet — but it should land before D2 wires the `ℓ_pen` reduce.

## 5. Top 3 priorities
1. **Mi1** — `FixedPointAccum::add` silently swallows non-finite inputs and saturates
   extreme ones; add a debug assert + doc the bound, so a broken `ℓ_pen` surfaces
   rather than vanishing to `0`.
2. **Mi2** — add the missing simulator tests (separated-het bimodality; per-group
   *shape* divergence), since A2's stated purpose is to make those exact properties
   testable downstream.
3. **Nit** — `below()` reads as a predicate but returns an index; rename.

## 6. Findings

### Minor

- `src/ssr/cohort/param_estimation.rs:157` — **[Minor]** `FixedPointAccum::add` silently swallows non-finite and extreme inputs
- **Confidence:** High (behavior); Medium (impact, no callers yet)
- **Problem:** `self.acc += (x * SCALE).round() as i128;` — a `NaN` input casts to `0`
  (Rust saturating float→int), and a value beyond `i128::MAX / SCALE` saturates to
  `i128::MAX`, both silently. The accumulator's reason to exist is making the `ℓ_pen`
  plateau a *sharp diagnostic* ("non-monotone ⇒ a bug"); a `NaN` from a broken E-step
  would be absorbed as `0` and hide exactly that.
- **Why it matters:** A determinism/diagnostic primitive that silently eats `NaN`
  defeats its own purpose once the D2 reduce starts feeding it real log-likelihoods.
- **Suggested fix:** assert finiteness in debug and document the magnitude contract.
  ```rust
  pub(crate) fn add(&mut self, x: f64) {
      debug_assert!(x.is_finite(), "FixedPointAccum::add got a non-finite value: {x}");
      self.acc += (x * SCALE).round() as i128;
  }
  ```
  Add to the doc comment: inputs are expected to be finite and bounded well within
  `i128::MAX / 2⁴⁰ ≈ 1.5e26` (per-locus log-normalizers are `reads × O(few nats)`).

- `src/ssr/cohort/sim.rs` (tests) — **[Minor]** Missing tests for two A2-purpose properties
- **Confidence:** High
- **Problem:** A2 exists to make "recover what we put in" testable, and the plan calls
  out **separated-het** and **per-group shape divergence** as required support. The
  current tests cover the homozygote mode and per-group *level* divergence, but not
  (a) that a separated het deposits support at *both* allele lengths, nor (b) that two
  groups differing only in stutter *shape* (`decay`) produce different magnitude
  distributions. A regression that collapsed a het to one mode, or ignored per-group
  shape, would pass today.
- **Why it matters:** These are the exact properties D1/D3/M3 will lean on; catching a
  generator bug here is far cheaper than debugging a recovery failure later.
- **Suggested fix:** see `## Missing tests`.

### Nits

- `src/ssr/cohort/sim.rs:72` — `below(n)` returns a uniform index but reads like a
  boolean predicate (cf. `chance`). Rename to `index_below` (or `uniform_index`).
- `src/ssr/cohort/sim.rs:~290` — the record builder computes the ref-tract length via
  `build_tract(&locus.motif, locus.ref_units).len()`, allocating a `Vec` only to read
  its length. Use `locus.motif.period() * locus.ref_units as usize`.

## 7. Out of scope observations
- `benches/psp_writer_perf.rs:386` — pre-existing index-out-of-bounds panic in the
  bench harness (`--all-targets`). Separate follow-up; not introduced here.

## 8. Missing tests to add now

Grouped under `sim`:

- `separated_het_deposits_support_at_both_allele_lengths`
  - **Input class:** a diploid `(4, 9)`-unit genotype at low `ε`, moderate level.
  - **Bug it catches:** a generator that drew all reads from one allele copy, or
    mis-tiled, would show one mode instead of two.
  - **Body:** simulate one such sample/locus; assert both `build_tract(motif, 4)` and
    `build_tract(motif, 9)` appear in `seq_counts` with non-trivial support.

- `differing_group_shape_changes_the_slip_magnitude_distribution`
  - **Input class:** two chemistries identical except `shape.decay` (e.g. `0.05` vs
    `0.5`), `ε = 0`, fixed allele, high depth.
  - **Bug it catches:** a forward model that ignored `decay` (per-group shape) when
    drawing magnitudes — the M3 axis A2 promises to support.
  - **Body:** call `simulate_cell` for each; compare the mean `|Δ|` of non-faithful
    reads (higher `decay` ⇒ larger mean magnitude).

- `fixed_point_accum_rejects_non_finite_in_debug` (paired with Mi1)
  - **Input class:** `NaN` / `+Inf`.
  - **Bug it catches:** silent absorption of a non-finite contribution.
  - **Body:** `#[should_panic]` debug test calling `add(f64::NAN)`.

## 9. What's good
- `FixedPointAccum`'s order/grouping-independence is proven by *bit-identity*
  (`to_bits()`), not approximate equality — the right bar for the determinism claim
  ([param_estimation.rs](../../../../src/ssr/cohort/param_estimation.rs)).
- The simulator drives the **real** `CohortMerger` reader path in its round-trip test
  rather than a shortcut, so it exercises decode + merge ([sim.rs](../../../../src/ssr/cohort/sim.rs)).
- A `BTreeMap<Box<[u8]>, _>` is used to get the "ascending by bytes" `observed`
  contract for free, matching the Stage-1 writer ([sim.rs](../../../../src/ssr/cohort/sim.rs)).
- Types seeded in their destination files with "shape, not final" notes where the
  shape genuinely depends on later milestones — honest forward-compatibility.

## 10. Commands to re-verify
- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-features`

### Author response convention
Address each finding by ID (Mi1, Mi2, Nits) with `fixed in <commit>` / `disputed` /
`deferred`.
