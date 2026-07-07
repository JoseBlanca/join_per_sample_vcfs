# Code Review + Fixes: genotype-EM core Phase 2.2

**Date:** 2026-07-07
**Reviewer:** rust-code-review skill (orchestrator, 4 category sub-agents)
**Scope:** working-tree diff of `src/var_calling/posterior_engine.rs` — introduce the
`GenotypeEmModel` trait + zero-sized `SnpModel`, route the M-step + convergence delta
through it (`run_em_loop` generic over `Model`). Byte-identical extraction.
**Status:** Approve-with-changes → fixes applied.

*(Consolidated review + fix-application report: this is a ~90-line byte-identical
extraction, so the two skill reports are merged into one file rather than split.)*

---

## 1. Execution status

- `cargo test -p pop_var_caller --lib` → 1616 passed, 0 failed, 4 ignored, no warnings.
- Byte-identity: tomato1 63-sample cohort VCF (189 938 records) `diff` = 0 vs the
  pre-Phase-2 (Phase-0) binary. Determinism floor already established at 0.
- `cargo clippy -p pop_var_caller --lib` → only the pre-existing `src/vcf/writer.rs:384`
  needless-lifetime warning (out of scope).

## 2. Categories dispatched

`refactor_safety`, `idiomatic`, `naming`, `reliability` (scoped to the diff; the
change is a byte-identical extraction with no error-handling, unsafe, or module-move
surface yet). Audit trail: `tmp/review_2026-07-07_genotype-em-2.2/`.

## 3. Findings and disposition

| ID | Sev | Categories | Finding | Disposition |
|---|---|---|---|---|
| M1 | Major | reliability | `ref_pseudocount_early_stop_regression_cohort` asserts only `base.iterations == heavy.iterations`; a *uniform* spurious early-stop survives | **Applied** — pinned absolute count `== 32` + `> 1` floor |
| Mi1 | Minor | naming, reliability | `convergence_delta` is a noun name that mutates `expected_counts_prev` (side effect undocumented) — misuse hazard (convergent) | **Applied** — documented the side effect on the trait method; kept the noun for consistency with `m_step` (naming agent's recommendation) |
| Mi2 | Minor | reliability | No test isolates the cross-record scratch-reuse guard or the single-sample no-copy branch | **Applied** — added two tests + `engine_for_records` helper |
| N1 | Nit | naming | `M` (MathBackend) sits beside new `Model`, momentary misread | **Won't fix** — `M` rename out of scope; noted for a future sweep |

**Adaptation on Mi2:** the reliability agent proposed a single-sample test asserting
`iterations == 1`. Verified empirically that single-sample converges at **iteration 2**
(the E-step is `p̂`-independent so the fixed point is reached at iter 1, but the
zero-delta convergence is only *detected* at iter 2). Pinned `== 2` instead; the
proposed `== 1` would have been wrong.

## 4. Tests added / changed

- `ref_pseudocount_early_stop_regression_cohort` — added `iterations > 1` floor and
  `iterations == 32` pin (M1).
- `single_sample_em_detects_convergence_on_the_second_iteration` — pins the
  single-sample no-copy branch (`== 2`, verified).
- `cohort_em_reusing_scratch_ignores_stale_expected_counts_prev` — drives two cohort
  records through one engine (new `engine_for_records` helper) and asserts the second
  runs past iteration 1, guarding the cross-record `expected_counts_prev` reuse.

## 5. What's good

- The extraction is a verbatim move: `m_step` → `convergence_delta` → `mem::swap`
  order, the `p_hat`/`p_hat_next` double-buffer, and the `expected_counts_prev` copy
  timing are all preserved (refactor_safety + reliability both confirmed).
- `EmContext` is `Copy`, so passing `ctx` by value into the trait methods is
  semantically inert.
- Static-dispatch generic (`Model: GenotypeEmModel`, no `dyn`) keeps the SNP copy at
  today's machine code, as the perf design requires.
