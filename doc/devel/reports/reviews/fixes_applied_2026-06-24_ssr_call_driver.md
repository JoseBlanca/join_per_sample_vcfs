# Fix Application Report: ssr_call_driver_2026-06-24.md

**Date:** 2026-06-24
**Source review:** `doc/devel/reports/reviews/ssr_call_driver_2026-06-24.md`
**Source state reviewed against:** branch `ssr-cohort`, HEAD `ce91077` (post ia→ai/doc move)
**Execution mode:** interactive
**Overall status:** In progress

---

## 1. Executive summary

### Review totals
- Blockers: 1 (B1)
- Majors: 6 (M1–M6)
- Minors: 16 (Mi1–Mi16)
- Nits: grouped

### Outcome totals (running)
- Applied: 0
- Deferred: 0
- (updated incrementally below)

### Validation summary
- Per-finding: `cargo test --lib --all-features ssr::cohort` + `cargo fmt --check` + `cargo clippy --lib`.
- Final full gate: `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --all-targets --all-features`, `cargo doc --no-deps` (recorded at end).
- Perf check: **skipped** — no bench under `benches/` references `ssr`/`cohort::` (confirmed `grep -rl "ssr\|cohort::" benches/` empty), so no Apply touches a bench-covered hot path.

### Unresolved high-priority findings
- (tracked at end)

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| B1 | Blocker | present-order VCF columns | Apply (test-first) | Applied | driver.rs, vcf_out.rs, inbreeding.rs | Pass |
| M1 | Major | `"?"` contig fallback | Apply | — | — | — |
| M2 | Major | `sample_chemistry` silent defaults | Apply | — | — | — |
| M3 | Major | duplicate `G₀` fallback const | Apply | — | — | — |
| M4 | Major | `partial_cmp().unwrap()` NaN-argmax | Apply | — | — | — |
| M5 | Major | unguarded contig/sample names | Apply (test-first) | — | — | — |
| M6 | Major | untested filtered-locus emit | Apply (test-only) | — | — | — |
| Mi1 | Minor | duplicated attribution helper | Defer | — | — | — |
| Mi2 | Minor | `LocusModel` bundle (11-arg) | Defer | — | — | — |
| Mi3 | Minor | per-round alloc churn | Defer | — | — | — |
| Mi4 | Minor | per-locus `f_present` alloc | Defer | — | — | — |
| Mi5 | Minor | tuple bins primitive obsession | Defer | — | — | — |
| Mi6 | Minor | bare `#[allow(too_many_arguments)]` | Apply | — | — | — |
| Mi7 | Minor | `..Default::default()` in test literals | Apply | — | — | — |
| Mi8 | Minor | `level_mult`→`level_multiplier` | Apply | — | — | — |
| Mi9 | Minor | `FrozenParams.params`→`chemistry` | Apply | — | — | — |
| Mi10 | Minor | stale "sweep is serial" comment | Apply | — | — | — |
| Mi11 | Minor | `EmCfg.inbreeding_f` test-only | Apply | — | — | — |
| Mi12 | Minor | unsurfaced threads/queue_depth coercion | Apply | — | — | — |
| Mi13 | Minor | `from_utf8_lossy` alleles | Defer | — | — | — |
| Mi14 | Minor | `QUAL=.` for variable-but-zero locus | Defer | — | — | — |
| Mi15 | Minor | refit non-convergence untested | Apply (test-only) | — | — | — |
| Mi16 | Minor | once-per-run `level_per_group.clone()` | Defer | — | — | — |

## 4. Per-finding log

### B1 — present-order VCF columns
- **Severity:** Blocker
- **Initial decision:** Apply (test-first)
- **Final status:** Applied
- **Reasoning:** Verified against current code: `format_vcf_record` iterated `call.calls` (present-only) producing `present_count` columns; the cursor yields `Ok(None)` for an absent sample and the merger sparse-omits only all-absent loci, so partial-coverage loci reach the formatter. Confirmed reachable by a test that demonstrably failed (15 cols vs 21).
- **Implementation summary:** `format_vcf_record` now takes `n_samples` and builds `sample_fields` dense at width `n_samples`, placing each present call at `locus.present[k]` and leaving absent samples as the `./.:.:.` placeholder. `n_samples` threaded from `run` → `write_genotyped_chunk` → `genotype_locus` → `format_vcf_record`. `genotype_locus` gained a justified `#[allow(clippy::too_many_arguments)]` (8 args; pre-empts Mi6 for this fn).
- **Review suggestion used verbatim?:** No (adapted — review sketched a `Vec<Option<&SampleCall>>`; used a `Vec<String>` placeholder fill keyed by `present`, fewer allocations).
- **Adaptation:** placeholder-string fill rather than an intermediate `Option` vector.
- **Verification performed:** new test `run_emits_dense_sample_columns_for_a_partial_coverage_locus` failed pre-fix (15 vs 21 cols), passes post-fix; full `ssr::cohort` suite green.
- **Files changed:** `src/ssr/cohort/vcf_out.rs`, `src/ssr/cohort/driver.rs`, `src/ssr/cohort/inbreeding.rs` (test call-site arg).
- **Tests added or modified:** `run_emits_dense_sample_columns_for_a_partial_coverage_locus` (new); `format_vcf_record` call-sites in vcf_out/inbreeding tests updated for the new arg.
- **Validation:**
  - `cargo test --lib --all-features ssr::cohort` → 0, `143 passed; 0 failed`
  - `cargo fmt --check` → 0
  - `cargo clippy --lib --all-features -- -D warnings` → 0
- **Follow-up:** None.
- **Residual risk:** None — `present[k] < n_samples` holds by construction (present indices are cohort sample ids `< n_samples`).

## 12. Notes

- Initial deferrals: Mi1/Mi2 (cross-file / signature refactors needing a design choice), Mi3/Mi4/Mi16 (allocation levers the review itself says to bench-gate, and no SSR bench exists yet), Mi5 (multi-site tuple→struct refactor), Mi13 (allele-byte validation boundary — design choice), Mi14 (QUAL `.`-vs-`0.0` is a policy choice the review flags as "either/or"). These remain open follow-ups.
