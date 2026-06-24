# SSR stutter / read-likelihood model bake-off (Stage 2 scoring)

**Status:** planned (branch `ssr-cohort`) — to be executed in a fresh conversation.
**Area:** Stage 2 (`ssr-call`) — the read likelihood `Qᵣ(obs | candidate)`.
**Background:** [ssr_delimiter_gap_penalty_2026-06-24.md](../reports/research/ssr_delimiter_gap_penalty_2026-06-24.md)
(the delimiter gap investigation that surfaced this), HipSTR sources under
`../pop_var_caller/HipSTR/`.

## 0. One-paragraph intent

The cohort caller scores each read's observed repeat tract against each candidate allele with
a read likelihood `Qᵣ(obs | cand)` — the probability the observed sequence arose from the
candidate via PCR/sequencing **stutter** (length change) plus **base error** (composition
change). That `Qᵣ` drives the per-locus EM, so its shape determines genotyping accuracy and
calibration. There are three credible ways to model it. Rather than argue them on priors, we
**implement all three behind one interface and bake them off on synthetic data with known
truth**, then pick the production model on measured accuracy / calibration / robustness /
speed. This plan defines the three models, the evaluation harness, the metrics, and the
decision rule. **The Stage-1 delimiter is already fixed (tract-aware gap, `26d21c1`) and is
out of scope here** — this is purely about scoring extracted tracts.

## 1. Scope and non-goals

- **In scope:** the Stage-2 read likelihood `Qᵣ(obs | cand)` ([likelihood.rs](../../../src/ssr/cohort/likelihood.rs),
  [stutter.rs](../../../src/ssr/cohort/stutter.rs), [pair_hmm.rs](../../../src/ssr/cohort/pair_hmm.rs)),
  its parameters, and how they are estimated in the pre-pass
  ([prepass.rs](../../../src/ssr/cohort/prepass.rs)).
- **Out of scope:** the Stage-1 delimiter (done); the genotype EM / FP-control / VCF emit
  (unchanged — they consume `Qᵣ`); real-data calibration of absolute parameter values (a
  later effort — this bake-off is about model *shape*, on synthetic data).
- **Important framing:** `obs` and `cand` here are **tract-only** sequences (the delimiter
  already stripped flanks). So "flank vs tract" does not arise inside `Qᵣ`; the models differ
  in how they treat **in-frame (whole-unit) vs out-of-frame (non-unit) length change** and
  base composition.

## 2. The three candidate models

All three compute `Qᵣ(obs | cand, motif, params)` and must be swappable behind a common
trait (see §3). Each also specifies how its parameters are estimated in the pre-pass.

### Model A — HipSTR (explicit in-frame + out-of-frame stutter)

Mirror HipSTR's `StutterModel::log_stutter_pmf` (`../pop_var_caller/HipSTR/src/stutter_model.cpp`):
marginalize over a length change `Δ_bp`, scoring it by a PMF that **splits in-frame from
out-of-frame**:

```
Qᵣ = Σ_{Δ_bp}  P_stutter(Δ_bp | period)  ·  align_subst(obs | resize(cand, Δ_bp))
P_stutter(Δ_bp):  in-frame  (Δ_bp = k·period):  up/down · geom_in ·(1−geom_in)^(|k|−1)
                  out-frame (Δ_bp ≠ k·period):  up/down · geom_out·(1−geom_out)^(|eff|−1)
```

- Two **distinct** geometric sub-models: in-frame in *units*, out-of-frame in *bp* (rarer).
- `align_subst` is substitution-only (length already matched by `resize`); placements for
  impure tracts as today (`reach_variants`).
- **New parameters:** `geom_out`, `out_up`, `out_down` (out-of-frame), estimated alongside the
  existing in-frame shape in the pre-pass — out-of-frame reads binned separately, not floored
  into the in-frame profile.

### Model B — Fix our current code

Keep the current decomposition ([likelihood.rs](../../../src/ssr/cohort/likelihood.rs)):

```
Qᵣ = Σ_{Δ units} S_θ(Δ) · Σ_v (1/|v|) · align_subst(obs | cand ⊕ Δ)
```

- `S_θ(Δ)`: in-frame stutter geometric in units (the only stutter term).
- `align_subst`: substitution in-tract + `FLANK_SLOP`-bounded boundary gaps; out-of-frame
  residual absorbed as boundary slop, larger → the `λ·(1/D)` outlier floor.
- **The fix:** stop the pre-pass flooring out-of-frame reads into the in-frame slip profile
  (`read_units = obs.len()/period` in [prepass.rs](../../../src/ssr/cohort/prepass.rs)) — drop
  or separately account out-of-frame reads so they don't bias the in-frame shape estimate.
- Smallest change; out-of-frame stays an approximation by design.

### Model C — Two-penalty unified pair-HMM (the user's idea)

A single forward pair-HMM `Qᵣ = forward(obs | cand)` with **two gap penalties applied
uniformly along the whole alignment** (no flank/tract or in-frame/out-of-frame *marginalization*
— the distinction lives in the gap transitions themselves):

- match / mismatch emission under `ε`.
- **unit gap:** insert or delete a whole motif unit (`period` bases in one transition) at the
  cheap **stutter** cost `g_unit` — allowed wherever the local motif phase permits.
- **base gap:** insert or delete a single base at the stiffer **error/out-of-frame** cost
  `g_base`.

Stutter is folded *into* the alignment as a unit-jump transition rather than marginalized
outside it. **Two parameters:** `g_unit`, `g_base`, estimated from the pre-pass slip/error
stats. This is the model we have not built; it needs a period-aware forward DP (a unit-jump
transition over `period` columns), the genuinely new piece.

### Model summary

| | stutter (in-frame) | out-of-frame | composition | new machinery |
|---|---|---|---|---|
| A HipSTR | geom in units | **explicit** geom in bp | subst-only | out-of-frame estimator + resize/placement |
| B fix-current | geom in units (`S_θ`) | boundary-slop gap + outlier | subst + flank-slop | pre-pass leak fix only |
| C two-penalty | unit-jump gap `g_unit` | single-base gap `g_base` | match/mismatch | period-aware unit-jump forward DP |

## 3. Common interface (swappable for A/B testing)

Introduce a trait so the EM and the harness call any model interchangeably:

```rust
pub(crate) trait ReadLikelihoodModel {
    /// Qᵣ(obs | cand) under this model's stutter/error treatment.
    fn q_r(&self, obs: &[u8], cand: &[u8], motif: &Motif, params: &Self::Params,
           scratch: &mut Self::Scratch) -> f64;
}
```

- The current `read_likelihood` becomes Model B's impl behind the trait (behaviour-preserving
  refactor first, gated by the existing likelihood tests).
- Each model owns its `Params` (estimated in the pre-pass) and reusable `Scratch` (hot path —
  `Qᵣ` is ~75% of pileup self-time per the perf review; keep per-model scratch).
- The cohort driver picks the model via config so the harness can sweep all three without
  forking the pipeline.

## 4. Evaluation harness (synthetic, known truth)

Reuse / extend the simulator [sim.rs](../../../src/ssr/cohort/sim.rs) (forward stutter draw +
substitutions, deterministic SplitMix64). The harness:

1. **Generate** a cohort of reads from **known genotypes** under a chosen *generative* model
   (see §5 fairness) at controlled depth / stutter rate / ε / motif / allele lengths.
2. **Score + genotype** the same reads through the cohort EM under each of Models A/B/C.
3. **Compare** the called genotypes to truth and record the metrics (§6).

Deliver it as a committed evaluation harness (a `#[cfg(test)]`/`examples/` driver, not a unit
test): it emits a results table (model × scenario × metric) so the decision is reproducible.

## 5. Fairness — cross-generative testing

A scoring model trivially wins on data its own assumptions generated. So generate under
**multiple generative models** and score each dataset with **all three** scoring models:

- **G1 — in-frame-only geometric** (the current `sim.rs` forward model).
- **G2 — HipSTR-style** in-frame + out-of-frame geometric (add out-of-frame draws to the sim).
- **G3 — "messy realistic"**: mixed in-frame + occasional out-of-frame + impure/interrupted
  tracts + a heavier substitution tail. The model that holds up here matters most.

Report the full 3×3 (scoring × generative) matrix per metric. Production choice weights G3
(robustness to model-mismatch) heavily — we do not know the true generative process.

## 6. Metrics

- **Genotype concordance** — fraction of samples whose called allele-length genotype equals
  truth (primary).
- **Allele-length error** — distribution of `|called − true|` in repeat units (captures
  near-misses vs the delimiter's old hard collapse).
- **Calibration** — reliability of GQ / posterior vs empirical correctness (a model can be
  accurate but overconfident).
- **Robustness sweeps** — concordance vs depth (5–100×), stutter rate (low→high), out-of-frame
  fraction, impurity, and allele length (incl. long alleles, the original failure).
- **Runtime** — per-`Q�r` cost and per-locus EM cost (the hot path); Model C's unit-jump DP and
  Model A's out-of-frame marginalization both add work — measure it.
- **Likelihood calibration (optional)** — does `Qᵣ(obs|cand)` track the true generative
  `P(obs|cand)` on held-out draws (KL / rank correlation)?

## 7. Decision rule

Pick the production model on, in priority order: (1) genotype concordance + allele-length
error on **G3 (messy)** and across the robustness sweeps; (2) calibration; (3) runtime within
the EM's budget; (4) implementation/maintenance cost. A model that wins only on data matching
its own assumptions (its diagonal) does not qualify. Record the rationale; the losing models'
code is removed (keep the trait if it aids future swaps, else inline the winner).

## 8. Order of work (for the fresh conversation)

1. **Trait + refactor:** introduce `ReadLikelihoodModel`; make the current code Model B behind
   it (behaviour-preserving; existing likelihood/EM tests green).
2. **Harness + metrics + G1:** the simulate→score→genotype→compare driver with the metrics
   table, on the existing generative model; Models A/B/C plug in as they land.
3. **Model C** (the user's idea — the genuinely new DP): period-aware two-penalty forward
   pair-HMM, with unit tests (faithful read, ±1/±k unit slip, single-base out-of-frame,
   impure tract, long allele).
4. **Model A** (HipSTR): explicit out-of-frame geometric + out-of-frame estimator + resize/
   placement; unit tests mirroring HipSTR's `log_stutter_pmf`.
5. **Generative G2 + G3**, then run the full 3×3 × scenarios bake-off.
6. **Decide + productionize** the winner; update the spec / params docs; remove the losers.

## 9. Risks / open questions

- **Parameter estimation parity:** each model must estimate its params from the *same*
  pre-pass stats, or the bake-off compares estimators, not models. Fix the estimator per model
  and document it.
- **Normalization:** `align_subst` already notes its unequal-length path slightly exceeds a
  proper transition normalization; the EM renormalizes per read, but a bake-off comparing
  likelihoods across models should note where each is/isn't a proper distribution.
- **Determinism:** all models must stay byte-identical across thread counts (the cohort
  contract) — pure per-read functions, as today.
- **Model C phase handling:** the unit-jump transition needs the motif phase at each column
  (cf. HipSTR's `num_upstream_matches`); get this right or impure tracts misbehave.

## 10. References

- Ours: [likelihood.rs](../../../src/ssr/cohort/likelihood.rs), [stutter.rs](../../../src/ssr/cohort/stutter.rs),
  [pair_hmm.rs](../../../src/ssr/cohort/pair_hmm.rs), [param_estimation.rs](../../../src/ssr/cohort/param_estimation.rs),
  [prepass.rs](../../../src/ssr/cohort/prepass.rs), [sim.rs](../../../src/ssr/cohort/sim.rs).
- HipSTR: `../pop_var_caller/HipSTR/src/stutter_model.cpp` (`log_stutter_pmf`),
  `src/SeqAlignment/{StutterAlignerClass,RepeatStutterInfo,AlignmentModel}.{h,cpp}`.
- GangSTR: `../pop_var_caller/GangSTR/src/` (read-class stutter likelihood — a third datapoint).
- Delimiter fix (separate, done): [ssr_delimiter_tract_aware_gap_2026-06-24.md](../reports/implementations/ssr_delimiter_tract_aware_gap_2026-06-24.md).
