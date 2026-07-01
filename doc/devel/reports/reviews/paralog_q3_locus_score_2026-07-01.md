# Code review — hidden-paralog filter Q3 (`score_locus_for_paralogy`)

## 1. Scope
- **Reviewed:** the new `src/paralog/locus_score.rs` (the per-locus H1-vs-H2
  marginal likelihood ratio) + re-export lines in `src/paralog/mod.rs`.
- **Against:** branch `tomato2-paralog-filter`, uncommitted diff.

## 2. Method
Two parallel review sub-agents per `ai/skills/rust-code-review`: (a) reliability
+ refactor_safety with a **numerical-correctness / prototype-faithfulness**
emphasis; (b) naming + errors + idiomatic + defaults + smells + module_structure.
Verification (dev container): `cargo test --lib paralog` 44 passed (was 37; +7
challenge/parity/property tests); `cargo clippy --lib --tests -D warnings` clean;
`cargo fmt --check` clean.

## 3. Verdict
**Numerically faithful and sound after fixes.** The numerical reviewer traced
the port line-by-line against `build_paralog_lr.py` and **confirmed** the H1
p-marginal (`LSE_num − LSE_norm` ≡ the prototype's normalised-SFS marginal), the
carrier-config enumeration (`(3,1),(4,1),(4,2),(6,1),(6,3),(8,1),(8,4)`
exactly), the H2 flat `LSE − ln(count)` marginal, the `LogSumExp` accumulator,
`ln_normal`/`log_add_exp`/Wright priors, the winsorise + veto thresholds, and
`N = cohort size`. No faithfulness bug. All findings were release-mode
safety/misuse on the pure public API.

## 4. Findings (all fixed)

- **M — `total_reads − alt_reads` u32 underflow** (both likelihoods). A
  malformed record with `alt_reads > total_reads` wraps to ~4e9 in release,
  silently poisoning the LR. → **Fixed:** clamp `alt_reads = min(alt_reads,
  total_reads)` once at the source (the usable-sample collection loop), so both
  likelihoods see a coherent pair.
- **M — zero/negative σ₀ → silent `NaN` LR.** `ln_normal` divides by σ₀ with no
  guard on the public entry point. → **Fixed:** drop samples whose σ₀ is
  non-finite or `≤ 0` in the collection loop (an unfit coverage model rather
  than a poisoned score).
- **M — observations/σ₀ length mismatch silently truncated the cohort** (was a
  `debug_assert` only; `zip` truncates in release while `1/2N` still counts the
  dropped sample). → **Fixed:** a hard runtime precondition — a length mismatch
  returns `ParalogScore::neutral()`; removed the now-redundant `debug_assert`.
- **Mi — H2 `PROB_FLOOR` on `P(non-carrier)`** (prototype floors only the
  carrier branch). Inert on the default grid (`min P0 = 0.16`). → **Fixed:**
  documented the intentional guard (protects non-default grids reaching `q→1`).
- **Mi — `PROB_FLOOR` doc said "add" but the code clamps.** → **Fixed:** reworded
  (clamp is equivalent to the prototype's `+1e-300` to f64 precision).
- **Mi — `Vec<(SampleObservation, f64)>` tuple primitive obsession.** →
  **Fixed:** named `UsableSample { observation, single_copy_depth_sd }`, threaded
  through both likelihood helpers.
- **Mi — flat `carrier_branch[s*configs.len()+j]` hand-rolled 2-D index.** →
  **Fixed (light):** named the stride (`n_configs`) + documented the
  `[sample][config]` row-major layout.
- **Nit — dead `total_reads > 0` veto guard.** → **Fixed:** documented it guards
  `homalt_min_depth == 0`.
- **Deferred (efficiency) — 5 per-locus `Vec` allocations on the genome-wide hot
  path.** The pure-function signature is the settled arch interface (Premise 2);
  reusable scratch belongs in the var-calling wiring (S2) where the per-locus
  loop lives and can be measured. Noted for S2.

## 5. Tests added (7)
`alt > total` clamp stays finite; degenerate-σ₀ sample dropped (not NaN);
length-mismatch → neutral; **an absolute-value parity anchor** (single sample on
collapsed grids: hand-computed `logL1 ≈ −7.1965`, `logL2 ≈ −11.3157` — pins the
magnitude the sign-only tests can't); odd/small-`T` config enumeration;
degenerate-grid midpoint; Wright priors sum to 1.

## 6. What's good
The pure-function contract (plain `ParalogScore`, neutral-on-empty, no `Result`)
is the right shape for a statistical kernel over validated data; named/documented
math constants; clean module placement (`super::ParalogModelParams` only, no
`var_calling` dependency, arch Premise 0).

Audit trail: `tmp/review_2026-07-01_paralog-q3/{reliability,naming_etc}.md`.
