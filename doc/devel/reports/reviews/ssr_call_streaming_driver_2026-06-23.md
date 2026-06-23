# Code Review: ssr-call streaming driver (Step H4)

**Date:** 2026-06-23
**Reviewer:** rust-code-review skill (orchestrator; inline synthesis)
**Scope:** commit `d17f524` on branch `ssr-cohort` — the two-pass streaming `run()` (burn-in → freeze → sweep → VCF), the task's definition of done
**Status:** Approve-with-changes

---

## 1. Scope

- **Reviewed:** PR diff, commit `d17f524` (Step H4).
- **In-scope files:** [driver.rs](../../../../src/ssr/cohort/driver.rs) — the rewritten `run`, `FrozenParams`, `collect_burn_in_subset`, `check_unique_sample_names`, `genotype_locus`, the `PLOIDY`/`BURN_IN_MAX_LOCI` constants, the `DuplicateSampleName`/`ThreadPool` error variants, the removal of the TSV `write_dump`/`format_locus` path, and the 3 integration tests.
- **Out of scope:** the called genotyping/pre-pass functions (reviewed in their milestones); `build_param_set` (H1).
- **Categories considered:** reliability, errors, naming, idiomatic, defaults, smells, refactor_safety, extras. `unsafe_concurrency`: the only concurrency is `ThreadPoolBuilder` + the rayon-internal `par_iter` in the burn-in (no `unsafe`/shared mutable state in this diff); covered under reliability.

## 2. Verdict

**Approve-with-changes.** The two-pass pipeline is correct and the end-to-end integration test proves the DoD (PASS variant emitted with the right genotypes, monomorphic locus dropped, decision-E hard error). No Blocker/Major. The Minors are about the **burn-in subset selection** (positional bias, and its coupling to the decision-E check), a zero-covered-loci error-clarity edge, and an end-to-end test gap (filtered records, cross-thread determinism).

## 3. Execution status

- `cargo fmt --check` clean · `cargo clippy --lib -D warnings` clean · `cargo test --lib` = **1270 passed, 0 failed, 2 ignored** (−4 TSV tests, +3 integration tests vs H3).
- Sub-agent fan-out not used (repeated transient overload on prior steps); reviewed inline. "Needs verification": 0.

## 4. Open questions and assumptions

1. **The decision-E check runs on the *subset*, not the full cohort.** `build_param_set` is called inside the burn-in on the bounded subset, so a sample whose confident loci all fall *outside* the first `BURN_IN_MAX_LOCI` would be flagged `UnresolvedSamples` even though it resolves elsewhere in the full cohort. With the positional (first-`cap`) selection this is a real (if unlikely at `cap = 20_000`) interaction. See M-Mi1. The fix is representative subset selection (the deferred calibration item) — does it also warrant moving the resolution check to a full-cohort pass? Flag for calibration.
2. **H1-Mi1 resolved as "accept the documented boundary."** H4 does not pass a genotyped-loci period set into `build_param_set`; the `shape_by_period` backfill universe stands (benign while `fallback_p == EM default == 0.5`). Recorded as the H4 decision.

## 5. Top 3 priorities

1. **Mi1** — document (and flag for calibration) that the positional subset both biases chemistry and can spuriously trip the decision-E check; representative selection is the fix.
2. **Mi2** — add cross-thread end-to-end determinism + a filtered-record emission test (coverage gaps).
3. **Mi3** — a zero-covered-loci cohort surfaces as `UnresolvedSamples`; consider a clearer error or document it.

## 6. Findings

### Minor

**Mi1: [driver.rs](../../../../src/ssr/cohort/driver.rs#L146-L151) — positional first-`cap` subset both biases the burn-in and couples to decision E.**
`collect_burn_in_subset` takes the first `BURN_IN_MAX_LOCI` loci in catalog order. Beyond the documented chemistry bias, the decision-E resolution check (`build_param_set`) runs on this subset, so a sample resolvable only at later loci would be wrongly rejected. Confidence: High (mechanical). At `cap = 20_000` over a genome it is unlikely to bite, and the integration cohort is far under the cap, so no test fails — but it is a latent correctness coupling, not just a statistical nicety.
*Fix:* none required for Step 1 beyond the doc note added here; the real fix is representative (reservoir/stratified) selection, already the deferred calibration item (reading Q-R5). Strengthen the `BURN_IN_MAX_LOCI` doc to name the decision-E coupling, and carry it as a calibration open item.

**Mi2: end-to-end test gaps — cross-thread determinism and filtered-record emission.**
The integration test runs at `threads = 2` only and asserts a PASS + a dropped locus; it does not assert (a) the same VCF at `threads = 1` vs `> 1` (the headline byte-identity property, end to end) nor (b) that a *filtered* locus is emitted with its reason (only the formatting is unit-tested in `vcf_out`). Confidence: High.
*Fix:* add a cross-thread determinism test (cheap — run twice, compare the VCF strings). A filtered-locus e2e is harder to construct deterministically; acceptable to defer to calibration, but note the gap. See section 8.

**Mi3: [driver.rs](../../../../src/ssr/cohort/driver.rs) — a cohort with no covered loci surfaces as `UnresolvedSamples`.**
If no sample covers any catalog locus, the subset is empty, `grouped` is empty, and every sample is reported `UnresolvedSamples` — technically true but it hides the root cause (no overlapping loci). Confidence: Medium (edge case).
*Fix:* optional — either document this in the `UnresolvedSamples` error, or add a distinct early error when the subset is empty. Low priority; defer unless it confuses real runs.

### Nits

- [driver.rs](../../../../src/ssr/cohort/driver.rs) `genotype_locus` takes 7 parameters (under clippy's threshold, so no warning); if it grows, bundle the four `*Cfg`s into a `GenotypeCfg` struct. No action now.
- The single-threaded sweep preserves catalog order by construction, so the merger's `_seq` is unused; that is correct for Step 1 (the writer-reorder is Milestone J). A one-line comment noting "order preserved by the serial sweep; reorder is J" would preempt the question.

## 7. Out of scope observations
None.

## 8. Missing tests to add now
Under `driver.rs`:
- `run_is_byte_identical_across_thread_counts` — write one cohort, `run()` with `threads = 1` and `threads = 4` to two outputs, assert the VCF bytes are identical (proves the burn-in determinism flows through to the VCF end to end).
- (Deferred) `run_emits_a_filtered_locus_with_its_reason` — needs a locus that admits to a non-PASS verdict while still present; constructing it deterministically is non-trivial, so defer to calibration. The FILTER formatting is already covered by `vcf_out::formats_a_no_call_filtered_record`.

## 9. What's good
- The two passes are cleanly separated and the subset is `drop`ped before the sweep, so peak memory is `max(subset, one streamed locus)` — the cohort-scaling property the design promised. [driver.rs run](../../../../src/ssr/cohort/driver.rs)
- The burn-in runs inside an explicit `ThreadPoolBuilder` scope returning a `Result`, so `--threads` is honoured and a pool-build failure is a typed error rather than a panic.
- `genotype_locus` returns `Option<String>` — the emit/drop decision is expressed in the type, and the single `candidates.admit == Pass && !is_variable` guard is the whole policy in one readable line.
- The integration test builds a *real* on-disk cohort (catalog + `.ssr.psp` per sample) and drives the actual `run()`, so it exercises the file I/O, the two-pass re-open, the header, and the emit policy together — a genuine DoD check, not a stub.

## 10. Commands to re-verify
- `cargo fmt --check` · `cargo clippy --lib --all-features -- -D warnings` · `cargo test --lib ssr::cohort::driver`

### Author response convention
Address Mi1–Mi3 + Nits by id; answer open questions 1–2 first.
