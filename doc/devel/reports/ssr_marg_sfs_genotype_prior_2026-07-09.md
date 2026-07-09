# Does the Ewens `θ/k` SFS frequency prior genotype SSRs better than the marginalized Dirichlet-multinomial? — a bounded probe

*Status: experiment + verdict, 2026-07-09, branch `ssr-marg-sfs` (off `main`, which
already carries all of `ssr-bic-emission`). Companion to
[`ssr_calling_models_and_optimization.md`](../specs/ssr_calling_models_and_optimization.md)
§3.2/§3.3 and [`ssr_freebayes_marginal_emission.md`](../specs/ssr_freebayes_marginal_emission.md).*

## 1. The question, and the decision it informs

We genotype SSRs with an expectation-maximisation (EM) loop and borrow the freebayes
site-frequency-spectrum (SFS) marginal only for the **emission** decision (emit the
locus, or drop it as monomorphic). A full freebayes-style **genotyper** — one that
replaces the EM and marginalises the allele frequency jointly across the whole cohort
to produce genotypes — is memory-costly at cohort scale and fights this project's
large-cohort thesis, so we decided not to build it.

That decision left one cheap question open. The genotyping EM's **per-sample genotype
prior** builds each plant's prior from a cohort **allele-frequency model**. Today that
model is a symmetric **Dirichlet-multinomial (DM)** with leave-one-out cohort borrowing
(`PVC_SSR_MARGINALIZED_PRIOR`, the "MARG" prior). Would swapping its frequency model for
the **Ewens `θ/k` SFS** — the same neutral spectrum that gives the freebayes *emission*
its edge — produce **better genotypes**, specifically on the contested lone-carrier
heterozygotes where the recall lives?

- **If yes:** worth knowing, and a hint the full genotyper might pay off.
- **If no:** the full genotyper is *definitively* not worth building — even the cheapest
  possible graft of its frequency model yields nothing.

This probe answers that without the memory cost.

## 2. What "the SFS genotype prior" is — and why it is a cheap graft (the CHEAP-ONLY guard)

The guard on this probe was: if bringing the SFS into genotyping needs cohort-joint
genotype-configuration structures, a second read-likelihood pass, or memory that scales
with cohort size, **stop** — that is the rejected joint genotyper. It does not.

The freebayes emission marginal (`freebayes_emit::ln_ewens_reduced`, `θ^K·∏1/n_a`) is a
**cohort-joint configuration** object: it enumerates integer allele-count vectors `n`
over all `2N` cohort chromosomes and weights each by `Ewens(n)·∏ₛ Lₛ(n/2N)`. The
Hardy-Weinberg genotype-frequency structure lives in the per-sample factors `Lₛ(p)`, not
in the Ewens weight. Reusing that object *literally* for genotyping means computing each
sample's genotype posterior from the joint count-vector enumeration (leave-one-out) —
which is exactly the joint genotyper the guard forbids.

The only cheap, per-sample, memory-neutral form is **mean-field**: summarise the rest of
the cohort by its expected allele-count spectrum (leave-one-out) and put an SFS-shaped
prior on the per-sample frequency `p`. And the mathematically faithful finite-`k` SFS
prior on `p` **is a symmetric Dirichlet(`θ/k`)** — that is literally what "Ewens `θ/k`"
names, `θ` diversity spread over `k` alleles. So:

> The cheap SFS genotype prior = the **existing** marginalized-DM machinery
> (`marginalized_genotype_log_priors`, `leave_one_out_alpha`), with its **base measure**
> changed from the mode-centred `G₀` pseudocounts to a flat symmetric `θ/k`.

Only the frequency model's base measure changes. The EM, the leave-one-out borrowing,
the DM integration, the read model, and the emission models are all untouched. It is
`O(k)` per sample, per-locus, allocates one length-`k` vector — memory-neutral. The
"one genuine design question" — how leave-one-out works under the SFS — resolves itself:
the base measure is additive to the same `E[cohort] − E[own]` leave-one-out counts the
DM path already computes (`α'ₛ = θ/k + max(0, E[cohort] − E[own])`).

This is also *why* we do not call `ln_ewens_reduced` directly: per-sample and sound, it
would have to be the finite-`k` Dirichlet form anyway (using it raw omits the HWE
structure); joint and literal, it trips the guard. The DM marginal with a `θ/k` base is
both.

## 3. Implementation

A new `EmCfg` field `sfs_base` (env `PVC_SSR_MARG_SFS=1`, a sibling of
`PVC_SSR_MARGINALIZED_PRIOR`), default off. It implies the marginalized prior and swaps
its base measure to `sfs_theta / k` per allele (`sfs_theta` = `PVC_SSR_FREEBAYES_THETA`,
default 0.01, shared with the freebayes emission). One helper, `marginalized_base`,
returns a borrow of `G₀` when off (byte-identical) and the flat `θ/k` vector when on.

`G₀` is mode-centred (`p^|Δ|`, the modal allele gets 1.0, off-modal alleles decay), total
mass ≈ 2; the SFS base is flat and sparse (`θ/k ≈ 0.0025`, total mass = `θ` = 0.01). So
the SFS base is ~200× weaker and allele-agnostic: it favours the monomorphic corners of
the frequency simplex (the neutral "rare alleles are unlikely" prior) but does not favour
*which* allele — the empirical leave-one-out counts do that. A lone off-modal carrier,
with no cohort support, is therefore held to the sparse neutral prior and pulled toward
homozygous; once an allele recurs, the leave-one-out counts grow and it is admitted.

**Files:** `src/ssr/cohort/em.rs` (`sfs_base` field, `marginalized_base`, `MIN_SFS_BASE`,
two unit tests), `src/ssr/cohort/driver.rs` (env wiring).

### Verification

- **Byte-identical when off:** default VCF is bit-for-bit identical to the pre-change
  `heuristic.marg0` baseline (`diff`, headers included).
- **Thread-deterministic:** `PVC_SSR_MARG_SFS=1` + freebayes at `--threads 1` vs `8` →
  byte-identical VCFs.
- **Tests:** 260 `ssr::cohort` tests pass (debug). Two new: `marginalized_base` returns
  `G₀` off / flat `θ/k` on; and the SFS base raises a lone off-modal carrier's posterior
  homozygosity relative to `G₀` (the predicted extra conservatism). `cargo fmt` /
  `clippy` clean.

## 4. Results (silver confident core: 561 true100 / 8850 false100)

Frequency-prior forms compared holding the emission model fixed:
**{plug-in, DM-marg (current), SFS-marg (new)} × {freebayes, bic, heuristic}**.

### 4.1 Freebayes emission — recall / FP frontier (QUAL swept post-hoc)

This is the clean comparison: freebayes marginalises the frequency in the *emit*
decision, so it does not confound the genotype prior with the emission's own frequency
handling. At matched emission, a recall gap at matched FP **is** a genotyping gap.

| QUAL ≥ | plug-in | **DM-marg** | **SFS-marg** |
|---|---|---|---|
| 0  | 88.9% / 1.24% | 92.9% / 1.48% | 91.3% / **0.52%** |
| 10 | 82.9% / 0.27% | 85.9% / 0.27% | 85.6% / 0.26% |
| 20 | 80.7% / 0.21% | 84.1% / 0.19% | 84.0% / 0.18% |
| 30 | 78.8% / 0.19% | 82.4% / 0.19% | 81.6% / 0.18% |
| 40 | 77.2% / 0.19% | 81.3% / 0.16% | 80.6% / 0.15% |
| 60 | 75.2% / 0.11% | 77.4% / 0.11% | 77.4% / 0.10% |

(DM-marg reproduces the known baseline exactly — Q≥20 = 84.1%/0.19% — validating the
pipeline.)

**At every QUAL floor from 10 upward, DM-marg and SFS-marg are indistinguishable:** recall
within 0.1–0.7 points, FP within 0.01–0.03%. Both clearly beat plug-in (~+3–4 recall pts);
**neither beats the other.** They lie on the same Pareto frontier.

The one place they differ is the **unfiltered** point (Q≥0): SFS-marg is ~3× cleaner on FP
(0.52% vs 1.48%) for −1.6 recall pts. Concretely, of the confident core it suppresses **85
of DM's 131 false-positive lone-carrier calls, at the cost of 9 true ones** — exactly the
extra conservatism the sparse `θ/k` base predicts. **But that gain is redundant with the
QUAL knob:** a QUAL≥10 floor prunes the same low-confidence loci, and the two priors
collapse onto one frontier. The SFS prior buys nothing you cannot get by raising the QUAL
threshold on the DM prior.

### 4.2 BIC (margin 0) and heuristic emission — single points

| run | prior | recall | FP% | | run | prior | recall | FP% |
|---|---|---|---|---|---|---|---|---|
| bic.marg0 | plug-in | 83.4% | 0.33% | | heuristic.marg0 | plug-in | 80.7% | 0.16% |
| bic.marg1 | DM-marg | 87.3% | 0.35% | | heuristic.marg1 | DM-marg | 56.7% | 0.11% |
| bic.margsfs | SFS-marg | 87.0% | 0.34% | | heuristic.margsfs | SFS-marg | 56.1% | 0.09% |

Under BIC, SFS tracks DM within 0.3 recall pts / 0.01% FP. Under the heuristic, both MARG
priors collapse identically to ~56% — the known MARG × per-sample-no-call mutual
exclusivity (the no-call removes exactly the hets MARG recovered); SFS inherits it
unchanged. No emission model separates SFS from DM.

### 4.3 HipSTR concordance (agreement, full emitted set)

| emission | DM-marg | SFS-marg |
|---|---|---|
| freebayes | 87.6% | 87.4% |
| bic | 87.0% | 87.1% |
| heuristic | 96.1% | 95.6% |

Identical within noise everywhere.

## 5. Verdict

**The Ewens `θ/k` SFS frequency prior does not genotype SSRs better than the marginalized
Dirichlet-multinomial — not meaningfully, and mostly not at all.** On the clean freebayes
frontier the two are the same instrument at any usable QUAL floor; under BIC and the
heuristic they track within noise; concordance is identical. The only real difference — a
cleaner unfiltered false-positive rate from the sparser base measure — is fully redundant
with the QUAL threshold the caller already exposes.

This is the expected and, honestly, the likely outcome. It matches the direct evidence we
already had: the freebayes *emission* beat BIC only modestly, which said the SFS prior's
edge over our EM-based frequency handling is small. That edge lives in **marginalising the
frequency in the emit decision**, not in a better **genotype** frequency model — and this
probe confirms it by isolating the frequency model inside genotyping and finding it inert.

**Consequence for the rejected full genotyper.** The cheapest possible graft of a freebayes
genotyper's frequency model — the SFS prior dropped straight into the per-sample genotype
prior — yields no genotyping improvement over the DM prior we already have. So the
memory-costly joint freebayes genotyper, which offers the same frequency model at cohort
scale, is **definitively not worth building.** This closes the question.

Depth (~3 reads/plant) remains the wall, per §8 of the models map: the last recall lives in
lone-carrier heterozygotes that no frequency prior can resolve at this coverage. The gold
standard, not a better prior, is the next move.

## 6. Artifacts

- Toggle: `PVC_SSR_MARG_SFS=1` (default off, byte-identical). `θ` via `PVC_SSR_FREEBAYES_THETA`.
- VCFs: `benchmarks/ssr_tomato1/results_rerun_20260708/ours/cohort/marg_sfs/`
  (`fb.{marg0,marg1,margsfs}.vcf`, `bic.*`, `heuristic.*`, `detT{1,8}.vcf`).
- Scorer: reuses `scripts/silver_standard.py` (`build_core`) + `lib/ssr_concordance.py`.
