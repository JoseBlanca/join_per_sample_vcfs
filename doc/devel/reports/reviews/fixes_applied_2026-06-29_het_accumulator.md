# Fix Application Report: het_accumulator_2026-06-29.md

**Date:** 2026-06-29
**Source review:** `doc/devel/reports/reviews/het_accumulator_2026-06-29.md`
**Source state:** branch `tomato2-paralog-filter`, B3 working diff
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers 0; Majors 3 (M1–M3); Minors 4 (Mi1–Mi4); Nits 1.

### Outcome totals
- Applied: 8 (M1, M2, M3, Mi1, Mi2, Mi3, Mi4, Nit)
- Deferred / Disputed: 0

### Validation summary
- `cargo fmt --check` → exit 0
- `cargo clippy --lib -- -D warnings` → clean on scope
- `cargo test --lib` → 1392 passed, 0 failed (+13 het; +5 vs pre-fix 8)
- `cargo doc --no-deps` → pre-existing `ClassicStutterModel` failure only
- `cargo audit` → not run (no dependency change)

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files |
|---|---|---|---|---|---|
| M1 | Major | `k > n` silent miscount (debug-assert only) | Apply | Applied | het.rs |
| M2 | Major | no exact-margin boundary test | Apply | Applied | het.rs |
| M3 | Major | no `k==n`/`k==0` test near gate | Apply | Applied | het.rs |
| Mi1 | Minor | min_depth inclusion untested | Apply | Applied | het.rs |
| Mi2 | Minor | no accumulator identity proptest | Apply | Applied | het.rs |
| Mi3 | Minor | unchecked sum vs validator | Apply | Applied | het.rs |
| Mi4 | Minor | NaN-freedom undocumented | Apply | Applied | het.rs |
| Nit | Nit | ln(10) literal repeated | Apply | Applied | het.rs |

## 3. Questions asked and answers
None.

## 4. Per-finding log (key items)

### M1 — `k > n` silent miscount
Applied via design change: `observe_site` now takes
`SiteCounts { ref_obs, alt_obs }` and derives `n_total = ref + alt`
(saturating). `alt ≤ total` is therefore a type invariant — a `k > n`
site is unrepresentable, the `debug_assert` is gone, and the
argument-transposition risk (two same-typed `u64`s) is eliminated. This
single fix closes the errors, reliability, and idiomatic angles of the
finding. The Stage-1 caller (C2) naturally has `ref_obs` (REF allele obs)
and `alt_obs` (Σ non-REF obs).

### M2 — exact-margin boundaries
Applied. `site_at_exact_gate_margin_is_rejected` (balanced site, margin =
its gate value → gate `> M` false → rejected) and
`site_at_exact_split_margin_is_ambiguous` (alt-biased site whose gate
exceeds the margin while `logLR == +M` → ambiguous, not het). Both
compute the margin from the LLs inline, pinning the strict-comparator
contract.

### M3 / Mi1 / Mi2
- M3: `pure_hom_alt_at_min_depth_is_hom_alt` (`ref_obs == 0` at `n == min_depth`).
- Mi1: `site_at_exactly_min_depth_is_scored` (`n_total == min_depth` admitted).
- Mi2: `accumulator_maintains_count_identity` proptest over arbitrary
  `(ref, alt)` sequences asserting `n_het + n_hom_alt + n_ambiguous == n_variant`.

### Mi3 / Mi4 / Nit
- Mi3: comment on the counters / `finish()` sum naming the genome-site
  bound (`<< u64::MAX`), matching the validator's `checked_add` intent.
- Mi4: comment tying the LLs' finiteness (and hence NaN-free `max`/`<`/`>`)
  to the constructor asserts + by-construction `alt ≤ total`.
- Nit: tests use `std::f64::consts::LN_10` instead of `2.302_585`.

## 5–8. Deferred / Disputed / Failed / Blocked
None.

## 9. Performance check
Skipped — accumulator not yet wired into the pileup (that is C2; cost
check in D3). The classifier is a few multiplies per site.

## 10. Commands run
- `./scripts/dev.sh cargo test --lib sample_summary::het`
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings [-A clippy::doc_lazy_continuation]`

## 11. Command results
- `cargo test --lib` → 1392 passed, 0 failed
- `cargo fmt --check` → exit 0
- `cargo clippy` → 5 pre-existing vcf/writer.rs errors only

## 12. Notes
- This step also amended the spec (§3) and architecture (Premise 1b/2/4)
  to the three-way binomial-LR het classification, and the B1 `HetCounts`
  type from two counts to four — all settled with the user before coding.
