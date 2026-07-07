# EM convergence criterion: test the driver, not the readout

**Status:** Phase 0 implemented, 2026-07-07, branch `em-convergence-criterion`
(the fix below is in `run_em_loop`; the previously-failing proptest passes and a
deterministic regression test is pinned). Model-intent spec — *what* quantity the
per-record EM should test for convergence, and *why* the previous choice made the
emitted allele frequencies and genotypes depend on the Dirichlet pseudocounts in
a way that was never intended. Grown from investigating a pre-existing proptest
failure (`var_calling::posterior_engine::larger_ref_pseudocount_cannot_increase_p_alt`)
turned up during the SSR work.

This is **Phase 0 of a larger arc**: unifying the SNP and SSR EM machines onto
one Dirichlet-multinomial core (design doc to follow). The fix here is
self-contained and mergeable on its own, and it moves the SNP engine onto the
same principle the unified core will use — "test the quantity the iteration
actually drives."

**Scope.** This is about the **stopping rule** of the cohort EM in
`src/var_calling/posterior_engine.rs` — the `while` loop in `run_em_loop`. It
does **not** change the E-step, the genotype prior, or the M-step formulas. It
changes *when* the loop decides it has converged, so that the stopping point
stops leaking the pseudocount magnitude into the emitted `p̂`. It is a
**call-changing** fix for the cohort path (not byte-identical); single-sample
records are untouched (§2, §4) and stay byte-identical. It lands behind the usual
re-baselining gate.

**The SSR EM is unaffected**, and for an instructive reason. The SSR π-loop
(`run_pi_em`, `src/ssr/cohort/em.rs`) uses a **plug-in HWE(π)** prior: the
frequency estimate π *is* fed back into the next E-step, and the G₀ pseudocount
is folded into π. There the tested quantity is the driver, so no leak. The SNP
engine has this bug precisely because it moved to the θ̂ site-frequency prior
(the `sfs_genotype_prior` work) which removed p̂ from the E-step — leaving p̂ a
dangling readout the stopping test still watched. The long-term unification
direction is SSR → the SNP-style marginalized core, **not** SNP → plug-in.

---

## 0. Vocabulary

- **EM (expectation–maximisation).** The per-record iteration the posterior
  engine runs: an **E-step** computes each sample's genotype posterior given the
  current parameters, an **M-step** re-estimates the parameters (`p̂`, `f̂_C`)
  from those posteriors, and the two alternate until the parameters stop moving.
- **`p̂` — the allele-frequency estimate.** The M-step's Dirichlet posterior mean
  over the allele simplex; `p̂[alt]` is the reported `AF` and feeds QUAL. This is
  a **readout**: a value the engine reports, derived from the E-step's output.
- **`expected_counts` — the driver.** The posterior-weighted allele copy numbers
  `Σ_{sample,genotype} posterior · (copies of allele in genotype)`. This is the
  quantity the cohort E-step actually feeds back into the next iteration's prior
  (§2). It is what *drives* the iteration; `p̂` is computed *from* it but does not
  feed back.
- **Dirichlet pseudocount.** The prior weight the M-step adds to the observed
  counts before normalising: `ref_pseudocount = 10`, `snp_alt_pseudocount = 0.01`,
  etc. Encodes "reference is usually the common allele."
- **Leave-one-out (LOO) prior.** In a cohort, each sample's genotype prior is
  informed by the *other* samples' allele counts, so a sample is not allowed to
  vote for its own prior. This is the coupling that makes the cohort EM a genuine
  multi-iteration process (§2).

---

## 1. The problem, in one sentence

The EM tests convergence on `p̂`, but `p̂` is a pseudocount-scaled readout of the
real iteration variable — so a larger pseudocount damps the per-iteration change
in `p̂`, trips the stopping threshold **earlier**, and halts the EM at a
different point on an otherwise-identical trajectory. The emitted allele
frequency then depends on the pseudocount through the stopping point, not just
through the intended prior.

### The evidence (the failing proptest, 3-sample cohort)

The invariant "a larger REF pseudocount cannot raise the alt frequency" is
sound at the EM's fixed point but is violated at the default stopping threshold:

| REF pseudocount | iterations to converge | emitted `p̂[alt]` |
|---:|---:|---:|
| 10 (default) | 30 | 0.00200 |
| 100 (×10) | **8** | **0.01913** ← *higher*, wrong direction |

Run both to a tight threshold (`1e-9`, cap 500) and the invariant is restored —
`p̂[alt]` = 0.000128 at pseudocount 100 versus 0.000848 at pseudocount 10, the
expected monotone direction. The violation is **purely early stopping**, not a
math error in the prior or the M-step.

---

## 2. Root cause — the stopping test reads a variable the iteration does not drive

Three facts about the cohort EM, together, are the whole bug:

1. **The E-step is pseudocount-independent.** The genotype prior is the
   θ̂-derived Dirichlet-multinomial site-frequency prior. The cohort leave-one-out
   step (`e_step_cohort_loo`) feeds the *previous* iteration's `expected_counts`
   into the next prior via `α'_s = α_species + (cohort counts − own counts)`. The
   pseudocounts appear **nowhere** in this path. So the trajectory of posteriors
   and `expected_counts` across iterations is **identical for every pseudocount
   setting** — same starting point, same E-step, same feedback variable.

2. **The pseudocount enters only the M-step readout.** `m_step_p_hat` computes
   `p̂[a] = (expected_counts[a] + pseudocount[a]) / denominator`. A larger
   `ref_pseudocount` inflates the shared denominator, which **compresses the
   per-iteration change** `|p̂_new − p̂_old|` — the same underlying motion in
   `expected_counts`, divided by a bigger number.

3. **Convergence is tested on `p̂`** — `run_em_loop` stops when
   `max_a |p̂_new[a] − p̂_old[a]| < convergence_threshold` (default `1e-3`).

Put together: two runs that share an identical `expected_counts` trajectory get
halted at **different iterations** purely because the pseudocount rescales the
delta being thresholded. The larger pseudocount trips the threshold sooner —
before the leave-one-out prior has finished pulling the alt frequency down toward
its fixed point — so it emits a *higher* alt frequency despite encoding a
*stronger* reference prior. The readout the loop watches is not the variable the
loop is actually moving.

**Cohort-specific.** A single sample uses `e_step` (no leave-one-out feedback),
whose posteriors are static across iterations, so `p̂` reaches its fixed point at
once and there is no trajectory to stop early on. The defect only bites at
`n_samples ≥ 2`, where the LOO prior makes the EM a genuine multi-iteration
process.

---

## 3. The fix — test convergence on `expected_counts`, the pseudocount-free driver

**Implemented (Phase 0).** The cohort stopping rule now watches the quantity the
E-step actually feeds back — `expected_counts` — instead of the derived,
pseudocount-scaled `p̂`. In `run_em_loop`:

- A new `expected_counts_prev` scratch buffer holds the previous iteration's
  counts. Each cohort iteration computes
  `last_delta = max_abs_diff(expected_counts_prev, expected_counts) / (ploidy · n_samples)`
  (option (a) below — normalise to a per-chromosome frequency), then copies
  `expected_counts → expected_counts_prev`.
- `p̂` and `f̂_C` are still computed each iteration for the emitted record; they
  are just no longer the cohort convergence signal.
- **Single-sample records keep the historical `p̂`-delta test verbatim.** They
  have no leak (§2) and stay byte-identical; the metric branches on
  `n_samples > 1`.
- The iteration-1 "never converge on the flat prior" guard is unchanged, which
  also means iteration 1's delta is discarded — so a reused-scratch
  `expected_counts_prev` needs no explicit reset.

Because `expected_counts` is pseudocount-independent, two configs that differ
only in pseudocount now **stop at the same iteration** — the same point on the
shared trajectory — and their emitted `p̂` differs *only* through the intended
M-step term. The monotonicity invariant is restored by construction, not by
tightening a threshold. The previously-failing proptest passes, and
`ref_pseudocount_early_stop_regression_cohort` pins both the invariant and the
equal-stopping-iteration property on the original counterexample.

### Why this is the right variable

The stopping test should ask "has the thing the iteration is moving settled?" The
thing the cohort EM moves is `expected_counts` (through the LOO prior). `p̂` is a
report computed *from* it and re-scaled by a prior constant; convergence of `p̂`
and convergence of `expected_counts` coincide *in the limit*, but the finite-
threshold crossing points differ by the pseudocount scaling. Testing the driver
removes the leak.

### Threshold re-interpretation

`convergence_threshold` currently means "max change in an allele *frequency*" — a
dimensionless number in `[0, 1]`, default `1e-3`. On `expected_counts` it becomes
"max change in a posterior-weighted allele *count*", whose scale is
`~ploidy · n_samples`. The threshold must be re-interpreted (and its default
re-derived) accordingly. Two clean options, to settle in §5:

- **(a) Normalise the counts** by `ploidy · n_samples` before differencing, so the
  delta stays a per-chromosome fraction on `[0, 1]` and the existing `1e-3` scale
  and its validation range `(0, 0.1]` carry over unchanged in meaning.
- **(b) Difference the posteriors directly** (each in `[0, 1]`), taking the max
  cell-wise change across the whole `n_samples × n_genotypes` posterior matrix.
  Also dimensionless, and arguably the most literal "have the posteriors stopped
  moving" test.

**Chosen: (a)** — it keeps the reported `final_max_delta_p` diagnostic and the
threshold on the same familiar frequency scale, so downstream tuning and the
`EMNoConv` filter semantics are undisturbed. Concretely the tested quantity is
`q = expected_counts / (ploidy · n_samples)`, the **pseudocount-free empirical
allele frequency** — the same thing `p̂` reports but without the prior nudge.

---

## 4. Scope of the behaviour change

This changes the stopping point of the cohort EM, so it **changes calls** — but,
as measured, only at the margin. Records that previously stopped early (with the
default `ref_pseudocount = 10`, damped p̂ deltas) now run closer to their fixed
point, so a few marginal calls settle differently. The dominant transition is
`0/1 → 0/0`: letting the leave-one-out prior converge nudges borderline hets to
hom-ref, and a handful of now-monomorphic sites drop out.

**Measured impact (old binary vs new, same `.psp` inputs; the diff is exact —
old-vs-old is byte-identical, so there is no thread-nondeterminism confound):**

| benchmark | shape | sites | Δsites | sites w/ GT change | GT cells changed | max ΔQUAL | wall |
|---|---|---:|---:|---:|---:|---:|---:|
| GIAB per_sample HG002 30× | 1 sample | 759 | 0 | 0 | 0 / 759 | 0 | — |
| GIAB mendelian trio | 3-sample cohort | 776 | −1 | 0 | 0 / 2 325 | 0 | — |
| tomato1 | 63-sample cohort | 189 990 | −52 (−0.027%) | 32 (0.017%) | 40 / 11.96 M (0.0003%) | 15.5 | ~66 s (was ~65 s) |

- **Single-sample is byte-identical**, empirically confirmed on GIAB per_sample
  (0 differences of any kind) — the fix does not touch the per-sample path (§2).
- **Cohort blast radius is tiny:** on the 63-sample tomato1 cohort only 0.017% of
  sites change any genotype and 0.0003% of genotype cells move; the GIAB trio
  moved a single site.
- **No perf regression:** the added iterations are on the small population of
  slow-converging records, so cohort wall time is unchanged (~66 s vs ~65 s).

---

## 5. What is decided, and what is open

**Decided and done (Phase 0):**

- The convergence test keys off the **pseudocount-independent driver**
  (`expected_counts`), not the pseudocount-scaled readout `p̂`. Root-cause fix;
  tightening the threshold or relaxing the proptest were rejected as band-aids
  that leave the hidden coupling in place.
- **Threshold basis: option (a)** — normalise the count delta to a per-chromosome
  fraction `q = expected_counts / (ploidy · n_samples)`, keeping the existing
  `1e-3` default and `(0, 0.1]` validation range unchanged in meaning. No default
  re-derivation needed.
- Single-sample behaviour is already at its fixed point and stays
  **byte-identical** (kept on the `p̂` test explicitly).
- The full test suite (1613 lib + integration + doctests) passes with the change,
  including the previously-failing proptest and the new deterministic regression.
- **Empirical impact measured and cleared (§4):** GIAB per_sample byte-identical;
  GIAB trio −1 site; tomato1 63-sample cohort 0.017% of sites change a genotype,
  0.0003% of cells, no perf regression. The old-vs-old determinism floor is exact
  (0 differences), so the diff is purely the fix. This was the remaining
  pre-merge gate.

**Open (for follow-up):**

- **`f̂_C` convergence** — `f_hat_compound` was *already* excluded from the
  stopping test before this change (the old test was `p̂`-only), so this fix does
  not regress it. Still worth confirming the driver-based test captures
  compound-allele convergence, or folding compound expected-counts into the delta.
- **Interaction with `EMNoConv`** — verify switching the tested variable does not
  push a materially different population of records into the non-convergence
  branch at the iteration cap (part of the cohort-benchmark check above).
