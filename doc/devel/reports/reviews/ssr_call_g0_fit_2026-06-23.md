# Code Review: ssr-call G₀ decay fit (Step G1)

**Date:** 2026-06-23
**Reviewer:** rust-code-review skill (orchestrator; inline synthesis — see Execution status)
**Scope:** commit `640a1e6` on branch `ssr-cohort` — fit the per-period `G₀` pseudocount-decay `p` in the SSR pre-pass
**Status:** Approve-with-changes

---

## 1. Scope

- **Reviewed:** PR diff, commit `640a1e6` (Step G1 of the `ssr-call` driver plan).
- **Reviewed against:** `640a1e6`, branch `ssr-cohort`.
- **In-scope files:**
  - [param_estimation.rs](../../../../src/ssr/cohort/param_estimation.rs) — `G0FitCfg`, `AlleleSpreadAccum`.
  - [prepass.rs](../../../../src/ssr/cohort/prepass.rs) — `PrepassStats.allele_spread_by_period` + its merge; the G₀ tally in `accumulate_locus` + `add_allele_copies`; `fit_g0_decay`; `estimate`/`run_prepass` gain `&G0FitCfg`; 5 new tests.
- **Out of scope:** [sample_groups.rs](../../../../src/ssr/cohort/sample_groups.rs), [inbreeding.rs](../../../../src/ssr/cohort/inbreeding.rs) (signature-only caller updates, verified correct); the rest of `src/ssr/cohort/`.
- **Categories dispatched:** reliability, errors, naming, idiomatic, defaults, smells, refactor_safety, extras (numerical kernel on the pre-pass path, config-bearing, diff-vs-intent). `unsafe_concurrency` skipped (no new primitives; the new reduce is integer-sum, order-independent). `tooling`/`module_structure` skipped (no `Cargo.toml` change; coherent two-file addition).

## 2. Verdict

**Approve-with-changes.** The fit is correct, matches the architecture (`ssr_call_driver.md` §9), and is well-tested for the diploid contract. No Blocker or Major findings. Three Minor findings (a fallback-consistency gap, a forward-compat divisibility assumption, two test gaps) plus Nits.

## 3. Execution status

- `cargo fmt --check` — **clean** (exit 0).
- `cargo clippy --lib --all-features -- -D warnings` — **clean** (exit 0).
- `cargo test --lib` — **1258 passed, 0 failed, 2 ignored.**
- **Sub-agent fan-out unavailable:** the parallel per-category dispatch failed on transient server overload (HTTP 500/529, zero work performed), so the orchestrator performed the category checks inline. Findings below were derived by direct reading of the diff and the downstream consumer (`em.rs::period_decay`); no category was skipped, but the redundant multi-agent cross-check did not run. Re-running the fan-out later would add confidence but is not expected to change the verdict.
- "Needs verification" findings: 0.

## 4. Open questions and assumptions

1. **Who owns "a `p` for every period present"?** Arch §9 says the pre-pass emits a decay for every period present. G1 emits one only for periods that have ≥1 *variable* locus (periods with only monomorphic loci produce no `allele_spread_by_period` entry, hence no `g0_by_period` entry). The downstream EM (`em.rs:115`) already covers the gap with its own hardcoded `p=0.5` fallback, so there is no crash — but two fallback values now exist. Resolve in H1 (see M-scoped Mi1): either the pre-pass tracks all observed periods and emits `G0FitCfg.fallback_p`, or we document the EM fallback as authoritative for unseen periods. Affects Mi1.

## 5. Top 3 priorities

1. **Mi1** — reconcile the two `G₀` fallbacks (`G0FitCfg.fallback_p` vs `em.rs` hardcoded `0.5`) so a period with no variable loci can't silently bypass the configured fallback.
2. **Mi2** — document/assert the diploid assumption behind the integer multiplicity `ploidy / peaks.len()` before polyploidy lands.
3. **Mi3** — add the missing het-multiplicity and asymmetric-spread tests.

## 6. Findings

### Minor

**Mi1: [prepass.rs](../../../../src/ssr/cohort/prepass.rs#L235-L239) / [em.rs](../../../../src/ssr/cohort/em.rs#L115-L122) — two independent `G₀` fallbacks; periods with no variable loci bypass `G0FitCfg.fallback_p`.**
`estimate` builds `g0_by_period` only from `allele_spread_by_period`, which is populated only by **variable** loci. A period whose loci are all monomorphic gets **no** entry, so at call time `period_decay` falls through to its own `const FALLBACK = { p: 0.5 }`. Today both fallbacks are `0.5`, so behaviour is correct — but they are two unsynchronised sources of truth: changing `G0FitCfg.fallback_p` would silently *not* affect unseen periods. Confidence: High.
*Fix (pick one, ideally in H1 where `build_param_set` consumes this):* (a) have the pre-pass record every period it encountered (variable or not) and emit `fallback_p` for those lacking sufficient spread — fulfils arch §9 literally; or (b) drop `em.rs`'s magic `0.5` and make `period_decay` take an explicit fallback sourced from `G0FitCfg`, documenting the EM path as the single owner. Note the choice in the H1 plan item.

**Mi2: [prepass.rs](../../../../src/ssr/cohort/prepass.rs#L158) — integer multiplicity `ploidy as u64 / peaks.len()` assumes divisibility.**
For diploid (the only supported ploidy — asserted in `run_locus_em`) this is exact: 1 peak → 2 copies, 2 peaks → 1 each. For a future polyploid with `peaks.len()` not dividing `ploidy` (e.g. tetraploid, 3 peaks → `4/3 = 1`), it silently under-counts copies and biases the spread. Confidence: High (latent, not reachable today).
*Fix:* add `debug_assert!(ploidy as usize % peaks.len() == 0 ..)` or a one-line comment tying the exactness to the diploid contract, so the assumption is visible when polyploidy work begins.

**Mi3: missing tests — het multiplicity and asymmetric-spread mode reference.**
`g0_spread_counts_only_variable_loci` and the recovery tests exercise the homozygote path only; the 2-peak (het → 1 copy each) branch of the multiplicity calc and the `rungs.modal_length()` mode reference under an asymmetric spread are uncovered. Confidence: High.
*Fix:* see section 8.

### Nits

- [param_estimation.rs](../../../../src/ssr/cohort/param_estimation.rs#L83) — `AlleleSpreadAccum` diverges from the module's `*Profile` / `*Stats` accumulator naming (`SlipProfile`, `SampleStutterStats`); `AlleleSpreadProfile` would be more consistent. Defensible as-is (it is a pair of sums, not a profile).
- [prepass.rs](../../../../src/ssr/cohort/prepass.rs#L158) — `peaks.len().max(1)` is a redundant guard: a `Confident(Peaks)` always has ≥1 peak (the `.expect("a confident genotype has ≥1 peak")` a few lines below relies on it). Harmless; could drop the `.max(1)` or keep it as defensive.
- [prepass.rs](../../../../src/ssr/cohort/prepass.rs#L214-L219) — `fit_g0_decay` relies on `k_bar > 0.0` being false for `NaN` to absorb a hypothetical `0/0` (only reachable if `min_copies == 0`, which never co-occurs with a real entry since entries carry ≥2 copies). Safe, but the NaN-absorption is implicit; a one-line comment or an explicit `total_copies == 0` guard would make it self-evident.
- [prepass.rs](../../../../src/ssr/cohort/prepass.rs#L171-L187) — the post-loop G₀ block lengthens `accumulate_locus`; extracting `accumulate_allele_spread(stats, period, &locus_alleles, mode)` would keep the function single-purpose. Optional.

## 7. Out of scope observations

None. The caller updates in `sample_groups.rs` / `inbreeding.rs` are correct signature ripples.

## 8. Missing tests to add now

Grouped under `prepass.rs`:

- `g0_spread_counts_het_alleles_as_one_copy_each` — a cohort of separated hets (e.g. 6/10, reusing `het_spec`) at a variable locus; assert `allele_spread_by_period[&2].total_copies` equals (n_samples × 2) and that distances are measured from `modal_length()` (covers the 2-peak multiplicity branch — currently only the homozygote `ploidy/1` path is tested).
- `g0_spread_measures_distance_from_the_modal_length` — a variable locus with an asymmetric allele distribution (mode at one length, a minority allele several units away); assert `distance_weighted / total_copies` matches the hand-computed K̄ from the modal length, pinning the mode-reference choice.
- (If Mi1 is fixed toward option (a)) `g0_emits_fallback_for_a_period_without_variable_loci` — a cohort whose only loci are monomorphic; assert `g0_by_period` carries the period with `p == fallback_p` rather than omitting it.

## 9. What's good

- The integer two-sum accumulator (`AlleleSpreadAccum`) is the right call — exact sufficient statistic for K̄, order-independent reduce, and the doc comment explicitly justifies sums-over-histogram (no fixed distance bound). [param_estimation.rs:74-89](../../../../src/ssr/cohort/param_estimation.rs#L74-L89)
- Determinism is preserved for free: the new reduce is integer addition and the existing `prepass_is_byte_identical_across_thread_counts` test compares full `EstimatedParams`, so it now also guards `g0_by_period` without a new test.
- The "variable loci only" filter is implemented exactly as the design requires and is directly tested (`g0_spread_counts_only_variable_loci`), with the modal-mode reference aligned to call-time `g0_pseudocounts`.
- `G0FitCfg` follows the module's established `dev_default()` + "pinned in F2" convention, keeping the provisional constants visible and per-field documented.

## 10. Commands to re-verify

- `cargo fmt --check`
- `cargo clippy --lib --all-features -- -D warnings`
- `cargo test --lib ssr::cohort`
- New tests from section 8 once added.

### Author response convention
Address each finding (Mi1–Mi3, Nits) with `fixed in <commit>` / `disputed because …` / `deferred to <step>`, answering open question 1 first.
