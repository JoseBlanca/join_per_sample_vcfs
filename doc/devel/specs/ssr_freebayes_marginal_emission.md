# SSR emission — freebayes-style joint marginal-likelihood polymorphism test

*Status: proposal, 2026-07-08, branch `ssr-freebayes-marginal` (off `main`). One of
two candidate replacements for the current heuristic emission gates in `ssr-call`
(Stage 2); the sibling is the BIC model-selection test
([`ssr_bic_confident_genotype.md`](ssr_bic_confident_genotype.md)), built separately.
This doc builds **only** the freebayes-style approach, behind a default-off toggle,
and reuses the existing cohort machinery — it introduces no new read model, no new EM,
and one small new prior. Where this and [`ssr_cohort_mark2.md`](ssr_cohort_mark2.md)
disagree on intent, Mark-2 wins.*

> **Merged-interface note (2026-07-08).** This model was merged into `ssr-bic-emission`
> and unified with the BIC model behind one selector: **`PVC_SSR_EMIT_MODEL=freebayes`**
> (replacing the standalone `PVC_SSR_FREEBAYES_EMIT` toggle described in §4.1). Its knobs
> are `PVC_SSR_FREEBAYES_MIN_QUAL` (the emit QUAL floor — the precision/recall knob,
> analogue of the BIC margin) and `PVC_SSR_FREEBAYES_THETA` (the SFS `θ`, default 0.01).
> Both emission models decide from the **same** one `data_ll` read-likelihood pass; the
> freebayes marginal now rides in `EmissionEvidence.freebayes_ln_p_mono` rather than a
> `LocusCall` field. The math below is unchanged.

---

## 0. Glossary

- **Emission decision** — the Stage-2 choice, per locus, of whether to write the locus
  to the VCF as a *variable* (polymorphic) call or drop it as monomorphic. Today this
  is `is_variable` + `apply_fp_control` in [`vcf_out.rs`](../../../src/ssr/cohort/vcf_out.rs);
  this doc replaces it (when the toggle is on) with a principled cohort-level test.
- **Segregation** — the cohort-level signal that separates a real allele from stutter:
  a real allele is read-*dominant* in some plants and *absent* in others; systematic
  stutter is a consistent low shoulder present in *every* plant. A per-plant genotyper
  cannot see this; a joint site-level marginal can.
- **`data_ll[s][g]`** — the per-sample, per-genotype read data log-likelihood
  `ln P(reads_s | G_s = g)` the EM already computes from the `Qᵣ` stutter read model
  ([`compute_data_ll`](../../../src/ssr/cohort/em.rs), em.rs:553). This is the *only*
  read-level input the new test consumes; per-chemistry, per-locus stutter is already
  baked into it via [`refine_theta_locus`](../../../src/ssr/cohort/stutter.rs).
- **SFS (site-frequency spectrum) prior** — the neutral-coalescent prior on how many
  chromosomes in the cohort carry each allele. freebayes uses **Ewens' Sampling
  Formula**; its `θ/k` shape penalises rare alleles (singletons, i.e. stutter
  shoulders) and is the key ingredient here.
- **`F` (inbreeding coefficient)** — the frozen per-sample excess-homozygosity term
  (tomato selfer, mean `F_IS ≈ 0.82`), threaded as `f_present` and already used in the
  Wright genotype prior [`genotype_prior`](../../../src/ssr/cohort/em.rs) (em.rs:156).
- **`P(monomorphic | data)`** — the posterior probability the locus is really fixed:
  every sample homozygous for one shared allele. `QUAL = −10·log10 P(monomorphic)`;
  `P(polymorphic) = 1 − P(monomorphic)`.

---

## 1. Motivation — why a joint marginal, not per-plant gates

The SSR caller genotypes a cohort of shallow samples (median ~3 reads/plant/locus).
Its false positives are all the **same failure**: systematic per-chemistry stutter
mis-attributed to a real allele. At 3 reads a stutter shoulder looks like a carrier;
a few coincidences plus frequency strength-borrowing amplify it into a cohort-wide
"variant". The current pipeline genotypes each plant *then* applies heuristic gates
(`is_variable` = "any non-ref MAP allele anywhere" + an allele-balance no-call in
`apply_fp_control`). Those gates are a crude proxy for the real discriminator —
**cross-sample segregation** — and they leak these FPs.

freebayes already solved this shape of problem for SNPs: it does not ask "did any
sample look non-ref?"; it asks **"is the whole cohort's read evidence better explained
by a polymorphic site than by a fixed one?"**, integrating a neutral SFS prior that
makes a rare "allele" pay for itself. We graft that structure onto our stutter read
likelihoods.

The current [marginalized-DM prior](../reports/ssr_marginalized_prior_benchmark_2026-07-07.md)
(`PVC_SSR_MARGINALIZED_PRIOR`) is a *genotype* prior inside the EM — it shapes each
sample's genotype call. This is different and complementary: a *site-level emission
test*. It does not touch the per-sample genotypes; it replaces the emit/QUAL decision.

---

## 2. The freebayes structure we port (verbatim source in `freebayes/`)

freebayes forms, per genotype-combination `C` across all samples, a log posterior
([`Genotype.cpp:1499`](../../../freebayes/src/Genotype.cpp)):

```
posteriorProb(C) = priorProbAf(C)     // Ewens SFS prior on the combo's allele partition
                 + priorProbG_Af(C)   // multinomial genotype-combo sampling prob
                 + Σ_s GL_s(G_s)       // data likelihood (our data_ll)
```

Then ([`freebayes.cpp:526`](../../../freebayes/src/freebayes.cpp)):

```
Z      = logsumexp_C posteriorProb(C)                     // marginal over the SFS
pHom   = Σ_{C fixed}  exp(posteriorProb(C) − Z)           // P(monomorphic | data)
QUAL   = −10·log10(pHom)                                  // ResultData.cpp:65
```

### 2.1 The Ewens SFS prior (the `θ/k` term)

[`Ewens.cpp:4-22`](../../../freebayes/src/Ewens.cpp), for an allele partition where
`a_j` distinct alleles each occur on `j` chromosomes and `M = Σ_j j·a_j` chromosomes total:

```
                    M!                        θ^{a_j}
P(partition) = ───────────────────  ·  ∏_j  ─────────────
               θ·∏_{h=1}^{M-1}(θ+h)          j^{a_j}·a_j!
```

- `θ^{a_j}/j^{a_j}` → each allele of count `j` contributes `θ/j`: **singletons undiscounted,
  common alleles penalised `1/j`** — the neutral-coalescent shape. A lone stutter
  "allele" (count 1 across the cohort) is exactly what this makes expensive relative
  to the fixed-site hypothesis.
- `θ·∏_{h=1}^{M-1}(θ+h)` = the Ewens rising-factorial normaliser (`θ^{(M↑)}`).
- `θ` = population-scaled diversity (freebayes default 0.01). Fixed constant here
  (§4.4), **not** the stutter `θ` (a different quantity) and **not** a user knob.

### 2.2 The genotype-combo multinomial term

[`Genotype.cpp:1391`](../../../freebayes/src/Genotype.cpp): `permutationsln − multinomialCoefficientLn(n, counts)`
— the number of ways to realise this per-sample genotype assignment given only the
pooled allele counts (het permutations in the numerator, population multinomial in the
denominator). **freebayes has no inbreeding term** in this factor (or anywhere); we
add `F` here (§3.2), because a selfer's excess homozygosity is real and load-bearing.

---

## 3. Our port — an exact biallelic cohort marginal

Enumerating genotype combos across 51 samples is exponential; freebayes controls it
with a banded search. We do not need that machinery, because the emission decision is
a **biallelic contrast** and factorises exactly.

### 3.1 Reduction to the two principal alleles

The FP mechanism is one stutter allele one motif-unit from the cohort mode; the real
signal is one line's allele segregating against another's. Both are **biallelic**. So
per locus we take the **two principal alleles** `M` and `A` — the top two by the EM's
cohort frequency estimate `call.pi` (`M` = argmax, `A` = runner-up) — and run the exact
biallelic marginal between them. Remaining candidate alleles (rare → the stutter set)
are outside the contrast; the SFS prior would suppress them anyway, and a genuinely
multi-allelic locus still segregates through its strongest pair. `k < 2` candidates →
trivially monomorphic (`QUAL = 0`, not variable).

*(General multiallelic extension — an exact dynamic program over the allele-count
simplex for the top-`Kmax` alleles — is noted in §7 as future work; it is unnecessary
for the tomato FP problem and would only add cost.)*

### 3.2 The exact biallelic marginal (what we compute)

Let `n = 2N` be the number of cohort chromosomes present at the locus (`N` present
samples, diploid), and let `κ ∈ {0,…,n}` be the number of `A` chromosomes, i.e. allele
frequency `p_A = κ/n`, `p_M = 1 − p_A`. For each present sample `s` with frozen
inbreeding `F_s`, the genotype-marginal read likelihood at frequency `p_A` is

```
L_s(κ) =   P(M/M | p, F_s)·exp(data_ll[s][MM])
         + P(M/A | p, F_s)·exp(data_ll[s][MA])
         + P(A/A | p, F_s)·exp(data_ll[s][AA])
```

with the **Wright genotype prior reused verbatim** from
[`genotype_prior`](../../../src/ssr/cohort/em.rs) (em.rs:156):

```
P(a/a | p,F) = F·p_a + (1−F)·p_a²        P(a/b | p,F) = (1−F)·2·p_a·p_b
```

The site marginal, in ln-space (order-stable → thread-deterministic):

```
term(κ)   = lnEwensBiallelic(κ; n, θ) + Σ_s ln L_s(κ)      // κ = 0…n
Z         = logsumexp_{κ=0}^{n} term(κ)
P(mono)   = [ exp(term(0)) + exp(term(n)) ] / exp(Z)        // fixed-M ∪ fixed-A corners
QUAL      = clamp(−10·log10 P(mono), 0, qual_cap)
```

`term(0)` is "all `M`", `term(n)` is "all `A`" — the two ways the locus is *fixed for
one allele*, matching the generalised monomorphic definition (not just fixed-ref).
`lnEwensBiallelic(κ)` is the §2.1 formula for the partition `{n−κ, κ}` (both classes
count 1; corners collapse to one class of count `n`), reusing `factorialln` / `gammaln`
from [`genetics.rs`](../../../src/ssr/cohort/genetics.rs).

### 3.3 What is reused vs new

| piece | source | new? |
|---|---|---|
| `data_ll[s][g]` (Qᵣ read likelihoods) | `compute_data_ll` (em.rs:553) | reuse |
| genotype indexing `Genotype{i,j}` | `enumerate_diploid_genotypes` (em.rs:145) | reuse |
| Wright genotype prior w/ `F` | `genotype_prior` (em.rs:156) | reuse (make `pub(crate)`) |
| per-sample frozen `F` | `f_present` (inbreeding.rs) | reuse |
| candidate alleles / `ref_idx` / `pi` | `CandidateSet`, `LocusCall` | reuse |
| `ln Γ` / `ln n!` | `genetics.rs` | reuse |
| **Ewens biallelic SFS ln-prior** | — | **new (~20 lines)** |
| **biallelic site marginal + QUAL** | — | **new (~40 lines)** |

The entire novel surface is one module, `freebayes_emit.rs` (the Ewens prior + the
marginal loop). No read model, no EM, no genotype-prior duplication.

---

## 4. Integration

### 4.1 Toggle (default off → byte-identical)

`EmCfg` gains `freebayes_emit: bool` (default `false`) and `sfs_theta: f64`, read once in
[`driver::run`](../../../src/ssr/cohort/driver.rs) with the established pattern:

```rust
freebayes_emit: std::env::var("PVC_SSR_FREEBAYES_EMIT").is_ok_and(|v| v == "1"),
```

Off ⇒ not a single new branch is taken ⇒ VCF byte-identical to today's
`cohort.ssr.vcf` (verified by `diff`).

### 4.2 Where it branches in

The marginal needs `data_ll` and `f_present`, both live inside
[`run_locus_em_with`](../../../src/ssr/cohort/em.rs) (em.rs:387). So we compute the
freebayes QUAL there, right after `genotype_pass` produces `call`, and store it on
`LocusCall` as `freebayes_qual: Option<f64>`. Then in
[`genotype_locus`](../../../src/ssr/cohort/driver.rs) (driver.rs:442-448), when the
toggle is on:

- **skip** `apply_fp_control` (the heuristic allele-balance no-call — this is what the
  principled test replaces);
- keep the record-emit gate `admit == Pass && is_variable(call)` (a locus whose MAP is
  all-ref carries no non-ref allele to report and is genuinely monomorphic);
- set `QUAL = call.freebayes_qual` instead of `site_qual(...)`.

Per-sample `GT:GQ:REPCN` columns are the **unchanged EM MAP calls**
(`final_calls`) — the task's "genotypes = the MAP configuration". The discrimination
lives entirely in the site `QUAL`: a stutter FP locus may still have a MAP het (so it is
emitted as a record) but its `QUAL` is low (`P(monomorphic)` high) and it falls below
any reasonable threshold; a real segregating locus earns a high `QUAL`. This is why the
evaluation sweeps `QUAL` thresholds (§6).

### 4.3 Determinism

Per-locus, no cross-locus reduction; the sample sum and the `κ` logsumexp are
fixed-order. Thread count cannot change the result — same guarantee the rest of the SSR
path already holds. A unit test asserts byte-identity across `--threads 1` vs `4`.

### 4.4 Fixed parameters (no tuning knobs)

Per the "principled replacement, no corroboration/threshold knobs" constraint, the only
free quantity is the SFS `θ`, fixed at a documented constant (start at freebayes'
`0.01`; the report records sensitivity). `qual_cap` reuses `FpControlCfg::qual_cap`
(200). Nothing here is exposed for per-run tuning.

---

## 5. Unit tests

1. **A clean segregating locus emits with high QUAL** — a synthetic cohort where ~half
   the plants are read-dominant for `A` and half absent (all `M`): `P(monomorphic)`
   tiny, `QUAL` large, emitted.
2. **A systematic-stutter locus does not** — every plant a strong `M` mode with a thin
   one-unit `A` shoulder (constant minority fraction): `P(monomorphic)` ≈ 1, `QUAL`
   near 0, filtered at any sane threshold.
3. **Determinism across threads** — byte-identical VCF at `--threads 1` and `4`.
4. **Default-off byte-identity** — a locus routed through both paths with the toggle off
   yields the current QUAL/emit exactly (guarded end-to-end by the `diff` check).
5. **Ewens prior sanity** — `lnEwensBiallelic` is a proper distribution over `κ` (the
   `θ/k` singleton penalty present), corners dominate under weak data.

---

## 6. Evaluation (identical to how BIC will be scored)

Data (reuse existing pileups; no re-pileup):

- pileups: `benchmarks/ssr_tomato1/results_rerun_20260708/ours/cohort/psp/*.ssr.psp` (51)
- catalog: `benchmarks/ssr_tomato1/results/ours/ssr_tomato1.ssr.catalog`
- HipSTR: `benchmarks/ssr_tomato1/results_ssr15k/hipstr/cohort.str.vcf.gz`
- read dump: `benchmarks/ssr_tomato1/results_rerun_20260708/our_reads.tsv`

**Silver standard** (reference impl = §4 of
[`ssr_error_signals_dashboard.py`](../../../benchmarks/ssr_tomato1/scripts/ssr_error_signals_dashboard.py)).
Two read-grounded label sets, then the **confident core** where both agree:

- `true100`: a non-ref tract length read-dominant (≥80% of a plant's reads, plant
  depth ≥4) in ≥2 independent plants (recurrence), or a corroborated balanced het —
  built once from our pileup reads and once from HipSTR's per-read fields.
- `false100`: no plant read-dominated by any non-ref length (all stutter).
- confident core ≈ **561 true / 8,850 false** loci.

Score across several `QUAL` thresholds, over the confident core:

- **recall** = fraction of the 561 `true100` loci emitted as variable (PASS,
  MAP-variable, `QUAL ≥ T`);
- **false-positive rate** = fraction of the 8,850 `false100` loci emitted as variable.
- also: total emission count, and genotype concordance vs HipSTR on comparable loci.

A standalone scorer, `silver_recall_fp_curve.py`, mirrors the dashboard's §4 classifier
exactly and is validated by reproducing the **current caller's 81% recall / 0.16% FP**
before it is trusted on the freebayes VCF.

**Comparison targets** (prior runs, same core):

| approach | recall | FP |
|---|---:|---:|
| current heuristic gates | 81% (108/561 missed) | 0.16% (14/8850) |
| BIC offline prototype (toy stutter, will differ in-caller) | 90/83/80/72% | 3.9/1.0/0.38/0.08% |

The freebayes recall-vs-FP curve across `QUAL` thresholds is reported against these.
Winner (freebayes vs BIC) becomes the caller's emission model.

---

## 7. Honest caveats & future work

- The silver standard is read-grounded, so it partly co-defines "segregation"; the
  non-circular check is the orthogonal truth set (in prep). Noted, not blocking.
- Median 3 reads/plant is a real information floor. The aim is to use the cohort
  optimally (segregation is a cohort-level signal), not to beat physics.
- **Multiallelic exactness** — the biallelic reduction (§3.1) is exact for the ≤2-allele
  contrast that dominates a selfer and targets the stutter FP directly. An exact
  count-simplex DP over the top-`Kmax` alleles (Ewens partition on the full vector,
  `F` as here) is the drop-in generalisation if an outbred multi-allelic cohort ever
  needs it; it changes only the prior/marginal, not the integration.
- **`θ` estimation** — a fixed `θ` is a deliberate simplification; a cohort-estimated
  Watterson `θ` (or per-period `θ`) is a later refinement, recorded as sensitivity now.
</content>
</invoke>
