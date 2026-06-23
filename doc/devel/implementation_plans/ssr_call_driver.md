# SSR Stage 2 — wiring the genotyper into the `ssr-call` driver (implementation plan)

**Status:** ready to implement, 2026-06-23, branch `ssr-cohort`. Turns the settled
architecture [doc/devel/architecture/ssr_call_driver.md](../architecture/ssr_call_driver.md)
(decisions A–E) into build order. Read the arch doc first; it owns the *why*, this plan
owns the *order* and the *process*.

**Goal.** Replace the driver's Phase-1 TSV dump with a real **VCF**: a bounded
**two-pass streaming** run (decision A) — a pre-pass/burn-in that freezes the
cross-locus-pooled parameters, then a single streaming genotyping sweep of independent
per-locus EMs — emitting a GangSTR/TRtools-style SSR VCF.

---

## The implementation loop (follow this for every step)

Each numbered **step** below is one pass through this loop. Do **not** batch steps.

1. **Implement** the step using
   [ia/skills/feature_implementation_skill.md](../../../ia/skills/feature_implementation_skill.md)
   (types first, then logic, then tests; fmt + clippy `-D warnings` + tests green).
   → **commit** (the implementation substep).
2. **Review** the step's diff using
   [ia/skills/code_review_skill.md](../../../ia/skills/code_review_skill.md). Save the
   review report under `doc/devel/reports/reviews/`.
   → **commit** (the review substep).
3. **Fix** the review findings using
   [ia/skills/code_review_fixes.md](../../../ia/skills/code_review_fixes.md). Save the
   fixes-applied report.
   → **commit** (the fixes substep; one or more commits if the fixes are large).
4. **Report & gate.** Write a summary in the conversation of how the step went
   (what was built, test results, review verdict, what was fixed, any deferrals), then
   **ask for permission to continue** before starting the next step.

**Commit hygiene.** Commit after **each substep** (implement / review / fixes), not just
at the end of a step. Use the project's commit-message convention
(`feat(ssr):` / `docs(ssr):` / `fix(ssr):` …) and the trailer required by the repo.

**Validation gate (every implement/fixes commit).** `cargo fmt --check`,
`cargo clippy --all-targets --all-features -- -D warnings`, `cargo test`. Report real
outcomes; never fabricate. Work in the dev container (`scripts/dev.sh`).

**Determinism gate.** Every step that touches the maths keeps Stage 2 byte-identical
across thread counts (the existing pre-pass / EM property) — add or extend a
thread-count test where the step introduces new reduces or parallelism.

---

## Milestones at a glance

| Milestone | What it delivers | Satisfies |
|---|---|---|
| **G** | `G₀` decay fit in the pre-pass (decision B) | removes a provisional constant; feeds `build_param_set` |
| **H** | **Step 1** — two-pass streaming driver + SSR VCF (freeze-everything) | the **original task's definition of done** (valid VCF end-to-end) |
| **I** | **Step 2** — per-locus stutter adaptation (local `θ_locus` + rate refit, shrunk) | removes the Step-1 frozen-param bias where data support it |
| **J** | Parallelism — the streaming worker pool + byte-identity (parallelism last) | throughput; `--threads` / `--queue-depth` live on the sweep |

Build order is **G → H → I → J**. G is a self-contained pre-pass feature the driver
consumes, so it lands first; H reaches a runnable VCF (the walking skeleton); I deepens
accuracy; J adds throughput. (Order principles inherited from the roadmap: types first,
walking skeleton early, simulator-driven tests, determinism throughout, parallelism last.)

---

## Milestone G — fit `G₀` in the pre-pass (decision B)

### Step G1 — per-period `G₀` decay from confident germline allele spread  ☐
*Arch:* [ssr_call_driver.md §9](../architecture/ssr_call_driver.md). *Source:* spec §4.3/§5.5.

- **Types:** add `g0_by_period: HashMap<u8, G0PseudocountDecay>` to `EstimatedParams`
  ([prepass.rs](../../../src/ssr/cohort/prepass.rs)); a per-period folded-distance
  accumulator on `PrepassStats` (integer counts, like `SlipProfile`); a small
  `G0FitCfg` (min included-allele-copy count, fallback `p`, steepness clamp floor).
- **Accumulate** (in `accumulate_locus`): for each locus with **≥2 distinct confident
  alleles cohort-wide** (variable-loci only), tally chromosome counts per allele length,
  take the modal length (tie-break on the smaller), and add each allele's count-weighted
  `k = |length − mode|` into the per-period histogram. Monomorphic loci skipped.
- **Fit** (in `estimate`): `K̄ = Σ k·C[k] / Σ C[k]`; `p = (√(1+K̄²) − 1)/K̄`; below the
  min-count fall back to the coded default `p`; clamp against over-tightening. Emit a `p`
  for **every** period present.
- **Alignment check:** define distance-from-mode the same way `g0_pseudocounts` does at
  call time (`rungs.modal_length()`).
- **Tests:** recovery (simulate germline alleles around a known-spread mode → `p` ≈
  formula); thin-period falls back to default; **variable-only** (adding monomorphic loci
  does not move `p`); byte-identical across thread counts.
- *Depends:* nothing new (extends D2/D3). *Note:* `build_param_set` consumes
  `g0_by_period` in **H1**; until then G1 is exercised through its own tests.

---

## Milestone H — Step 1: two-pass streaming driver + SSR VCF

Freezes **all** params from the pre-pass (ε, `θ_period`, per-group level line, `G₀`,
`F`) and genotypes the full cohort in one streaming sweep. Reuses the existing per-locus
EM (`run_locus_em` runs on fixed params today — `em.rs:7`), so no per-locus refinement
yet (that is Milestone I).

> **Burn-in reuse.** The burn-in is `run_cohort_em` run on the **bounded subset** (it
> already returns converged `f_per_sample` + `level_per_group`); we freeze those. So
> `run_cohort_em` is **repurposed as the subset burn-in**, not retired — what is retired
> is using it over the *whole* cohort: the full genotyping pass calls `run_locus_em`
> per streamed locus on the frozen params.

### Step H1 — `build_param_set` (pre-pass → frozen `ParamSet`)  ☐
*Arch:* [ssr_call_driver.md §3, §7](../architecture/ssr_call_driver.md) (decision E).

- `build_param_set(&EstimatedParams, &GroupedParams, n_samples) -> Result<ParamSet, SsrCallError>`:
  the field-by-field assembly currently inlined in the test, reading `g0_by_period` (G1)
  into `pseudocount_decay_per_loci_group`.
- **Decision E:** a cohort sample absent from `grouped.group_of_sample` (no confident
  genotype anywhere) is a **hard error** naming the offending sample(s) — new
  `SsrCallError::UnresolvedSamples` variant. No group-0 default; no m2(a) fallback.
- **Carry-over from the G1 review (Mi1):** `est.g0_by_period` only carries periods with
  ≥1 *variable* locus; a period whose loci are all monomorphic is absent (the EM then
  uses its own hardcoded `p=0.5` — two unsynchronised fallbacks). In `build_param_set`,
  fill `pseudocount_decay_per_loci_group` with `G0FitCfg.fallback_p` for **every period
  present** in the genotyped loci, so the configured fallback is the single source of
  truth (arch §9 "a `p` for every period present"). See
  [review](../reports/reviews/ssr_call_g0_fit_2026-06-23.md) Mi1.
- **Tests:** happy path (round-trips every field, dense `group_of_sample`); the
  unresolved-sample hard error; `G₀` carried through for every period present.
- *Depends:* G1.

### Step H2 — merger accessors + second-pass re-open  ☐
*Arch:* [ssr_call_driver.md §5](../architecture/ssr_call_driver.md) (decisions C/C-1).

- On `CohortMerger`: `chromosomes() -> Vec<ParsedChromosome>` (name + **length** + md5,
  global-id order — lengths come from the `.ssr.psp` headers, **not** the catalog) and
  `sample_names() -> Vec<String>` (basename without `.ssr.psp`, from the private
  `labels`).
- The driver runs **two passes**, so it re-opens the merger for the genotyping sweep
  (it holds the catalog + psp paths). Confirm `open` is cheap to call twice (cursors
  re-seek from the block index; reading §8 *re-read, don't cache*).
- **Tests:** `chromosomes()` carries lengths/order; `sample_names()` basename derivation;
  re-open yields the same locus stream.

### Step H3 — the SSR VCF header writer  ☐
*Arch:* [ssr_call_driver.md §5, §7](../architecture/ssr_call_driver.md) (decision C).

- In [vcf_out.rs](../../../src/ssr/cohort/vcf_out.rs):
  `write_vcf_header(out, &chromosomes, &sample_names, &warnings)` emitting plain text:
  `##fileformat=VCFv4.4`; one `##contig=<ID,length>` each; `##INFO=<ID=PERIOD,…>`;
  `##FORMAT` for `GT`/`GQ`/`REPCN`; `##FILTER` for `notPeriodic`/`tooManyAlleles`/
  `lowDepth`; the apparent-`F_IS` `##` warning when present; then the
  `#CHROM … FORMAT <sample_1…N>` line in merger order.
- Dedicated SSR writer — **not** the SNP `src/vcf/writer.rs` (settled; the two diverge,
  and `format_vcf_record` already emits text data lines).
- **Tests:** header lines present & correct; contig lengths; sample columns in order;
  warning line appears iff supplied; output parses as a VCF header.

### Step H4 — the streaming driver (burn-in → freeze → sweep → VCF)  ☐  ← **task definition of done**
*Arch:* [ssr_call_driver.md §2, §4, §6, §8](../architecture/ssr_call_driver.md).

- **Pre-pass / burn-in:** stream the merger once, collect a **bounded stratified subset**
  of `CohortLocus` into a `Vec` (simple cap + per-period stratification; selection
  strategy is calibratable, reading Q-R5); `run_prepass_stats` → `estimate`
  → `group_samples` → `build_param_set` (H1); then `run_cohort_em` **on the subset** →
  freeze `f_per_sample` + `level_per_group`.
- **Genotyping sweep:** re-open the merger (H2); for each streamed `CohortLocus`:
  `build_rungs` → `assemble_candidates` → `run_locus_em_with(frozen F, frozen level)` →
  `apply_fp_control` → emit policy (emit-iff-variable for PASS; filtered always emitted
  with its reason; monomorphic PASS dropped — §6) → `format_vcf_record`.
- **Header:** `write_vcf_header` (H3) before the records; `f_is_warning` folded into the
  warnings.
- **Threads:** map `config.threads` to the rayon pool (the pre-pass is already parallel;
  the sweep stays **single-threaded** in Step 1 — parallelism is Milestone J).
  `config.queue_depth` documented as the future sweep knob.
- **Integration test (DoD):** the **real** `run()` over a simulated-cohort `.ssr.psp`
  set produces a valid VCF — header + records, genotypes correct at high depth,
  a PASS variant with an ALT, a filtered record kept with its reason, a monomorphic
  locus dropped; plus the decision-E hard error when a sample never resolves.
- *Depends:* H1, H2, H3, G1. **Completes the original task's definition of done.**

---

## Milestone I — Step 2: per-locus stutter adaptation (local, shrunk)

Builds the deferred per-locus refit (`em.rs:12` — "θ_locus M-step deferred until D wires
the slip accumulators"). Local to each locus ⇒ stays a single streaming sweep; shrinkage
toward the frozen priors is the anti-oscillation knob.

### Step I1 — per-locus slip attribution + `θ_locus` shape M-step  ☐
*Arch:* [ssr_call_driver.md §4](../architecture/ssr_call_driver.md).

- In the per-locus EM, attribute each read to its genotype-weighted parent allele →
  per-locus slip profile → refine `θ_locus` via `refine_theta_locus(profile, frozen
  θ_period, shape_shrink_strength)`.
- **Tests:** a high-depth locus adapts toward its own MLE; a thin/ambiguous locus stays
  anchored at `θ_period` (no oscillation); byte-identical across thread counts; the
  Step-1 VCF integration test still passes (calls only sharpen).

### Step I2 — shrunk per-locus stutter **rate** refit  ☐
*Arch:* [ssr_call_driver.md §4](../architecture/ssr_call_driver.md).

- From the same slip attribution, refine this locus's stutter **rate**, shrunk toward
  the frozen per-group level line (the group line supplies the length-dependence; the
  locus nudges the overall rate). Remove the Step-1 frozen-rate reliance where superseded.
- **Tests:** rate adapts at high depth, anchored when thin; byte-identical; recovery on a
  simulated locus with a deliberately off-group rate.

---

## Milestone J — parallelism on the genotyping sweep (parallelism last)

### Step J1 — bounded-queue worker pool + writer reorder + byte-identity  ☐
*Arch:* [ssr_call_driver.md §4, §8](../architecture/ssr_call_driver.md);
[ssr_call_reading.md §5](../architecture/ssr_call_reading.md).

- Wire the genotyping sweep onto the reading layer's producer → bounded-queue → worker
  pool → writer (reorder by locus seq) topology; `config.threads` drives the pool,
  `config.queue_depth` the queue (the peak-resident-loci RSS knob).
- **Tests:** byte-identical VCF across thread counts (each locus is already a pure
  independent function — §4); back-pressure bounds resident loci.
- *Note:* optional for the task's DoD; do it once the maths (I) is settled.

---

## Cross-cutting notes

- **Reports.** Save an implementation report per step under
  `doc/devel/reports/implementations/` and review/fixes reports under
  `doc/devel/reports/reviews/`, per the skills.
- **PROJECT_STATUS.md.** Update the Stage-2 `ssr-call` block at the **end of each step**
  (status, report links, closed/opened items, refresh *Last completed task*), per the
  feature-implementation skill's project-status protocol. Do not touch other blocks or
  the About paragraph.
- **Out of scope (whole plan):** calibrating the `dev_default` constants (incl. the `G₀`
  fallback `p` / clamp floor and the subset size); the exact-AF QUAL kernel;
  beta-binomial; polyploidy (ploidy hard-coded 2, the non-diploid panic is the contract —
  decision D). No fabricated calibration numbers.

---

## Step checklist

- [ ] **G1** — `G₀` decay fit in the pre-pass
- [ ] **H1** — `build_param_set` (+ decision-E hard error)
- [ ] **H2** — merger `chromosomes()` / `sample_names()` + second-pass re-open
- [ ] **H3** — SSR VCF header writer
- [ ] **H4** — streaming driver → VCF *(task definition of done)*
- [ ] **I1** — per-locus `θ_locus` shape refit (shrunk)
- [ ] **I2** — per-locus stutter rate refit (shrunk)
- [ ] **J1** — parallel streaming sweep + byte-identity
