# Code review — hidden-paralog filter P1 (`callable_positions`)

## 1. Scope
- **Reviewed:** the working-tree diff implementing plan step **P1** (store the
  per-sample callable-position total).
- **Against:** branch `tomato2-paralog-filter`, uncommitted diff.
- **In-scope files:** `src/sample_summary/mod.rs`, `src/sample_summary/coverage.rs`,
  `src/sample_summary/het.rs`, `examples/dump_sample_summary.rs`,
  `src/vcf/writer.rs` (one unrelated pre-existing doc-lint fix).

## 2. Method
Two parallel review sub-agents per `ai/skills/rust-code-review`:
reliability + refactor_safety; naming + errors + idiomatic + defaults + smells +
module_structure. Verification (dev container): `cargo test --lib` 1395 passed;
`pileup_cli_integration` + `psp_to_pileup_integration` 11 + 4 passed;
`cargo clippy --lib --tests --example dump_sample_summary -D warnings` clean;
`cargo fmt --check` clean. `cargo clippy --all-targets` fails **only** in
out-of-scope benches/examples (`min_mapq_diff_t`/`no_mapq_diff_filter`/
`DEFAULT_MIN_MAPQ_DIFF_T` — a removed MAPQ-diff CLI API, broken at HEAD).

## 3. Verdict
**Sound.** No Blocker or Major findings. Two Minor (both doc-consistency), a
handful of missing boundary/regression tests, two nits. All refactor-safety
checks clean (every `CoverageByGcHistogram` construction site updated; no
`..Default` catch-all hiding a missed field).

## 4. Findings

- **Mi1 (naming, `mod.rs`, HetCounts doc).** The pre-existing `HetCounts` doc
  still asserted `Hobs = n_het / (n_het + n_hom_alt)`, which P1's own new
  `callable_positions` doc and the edited example explicitly contradict. One
  file carried two definitions of `Hobs`. → **Fixed**: doc rewritten to define
  `Hobs = n_het_sites / callable_positions` and label the confident ratio a
  reference-divergence proxy.
- **Mi2 (refactor_safety, `mod.rs`, version doc).** A stale on-disk v1 `.psp`
  fails as a missing-field `ParseToml` *before* the version guard runs (the
  required field has no `#[serde(default)]`), so the documented
  `UnsupportedVersion { got: 1 }` path is unreachable for real v1 files.
  Fail-closed, but the diagnostic is misleading. → **Fixed (option a)**: version
  doc-block now states v1 is rejected at the TOML layer, guard catches only
  future/zeroed versions; remedy is re-run `pileup`.
- **Missing tests.** The `callable_positions >= n_tiles` invariant, its
  `== n_tiles` accept boundary, and the multi-tile `finish`-drain sum were
  unpinned. → **Fixed**: added `validate_rejects_callable_positions_below_n_tiles`,
  `validate_accepts_callable_positions_equal_to_n_tiles`,
  `callable_positions_counts_every_covered_position_across_many_tiles`.
- **Nit — example het-rate label** ("rate/Mb" terse). → **Fixed**:
  `"Hobs (het/Mb)  (n_het_sites per Mb callable)"`.
- **Nit — `callable_pos` abbreviation.** → **Skipped**: the abbreviation
  preserves the dump's aligned label column; cosmetic.
- **Cross-category — `saturating_add` cap has no test.** → **Skipped with
  rationale**: testing the `u64::MAX` cap would need a `#[cfg(test)]` mutator on
  a private accumulator field for a practically-unreachable path; the sibling
  `depth_sum` saturating_add is untested for the same reason, so adding surface
  only here would be inconsistent. Behaviour is documented.
- **Cross-category — zero-callable guard.** `validate()` admits
  `callable_positions == 0` (all-`N` sample, `n_tiles == 0`). Correct; the
  downstream `Hobs` consumers (built in S1/Q4) must guard the zero denominator.
  Carried forward to those steps.

## 5. What's good
- The field is on the right type (`CoverageByGcHistogram` — it sums the coverage
  stage's per-tile `covered`, which the het path never sees; no back-reference).
- Version bump is a documented `pub const` enumerating each version's additions.
- Typed validation on both serialise and parse paths; `saturating_add` matches
  existing accumulator discipline.

Audit trail: `tmp/review_2026-07-01_paralog-p1/{reliability_refactor,naming_etc}.md`.
