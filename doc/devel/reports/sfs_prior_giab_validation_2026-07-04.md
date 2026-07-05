# SFS genotype prior — GIAB validation

**Date:** 2026-07-04
**Scope:** end-to-end GIAB per-sample validation of the SFS-marginalized genotype
prior (spec/arch/plan `sfs_genotype_prior`), wired into the posterior engine
(commit on branch `sfs-genotype-prior`). This is step 2 of the revised plan
(implement → check GIAB → remove old code if improved).
**Verdict:** clear improvement, **zero precision/recall cost** → proceed to code removal.

---

## Setup

`PRESET=high-recall benchmarks/giab/src/run_ours_per_sample.sh {5,10,15,30}x`
(HG002/3/4, single-sample), new prior on by default. Compared against freebayes
and against the committed pre-change baseline (`freebayes_comparison.tsv`).
Analysis: `tmp/gt_confusion_new.py` (GT concordance), `tmp/prec_recall.py`
(precision/recall/FP).

> Note on running the benchmark on macOS: `run_ours_per_sample.sh` execs the
> binary directly on the host, so it uses the **host** `target/release` build,
> not the container `target-container/release` build. Rebuild on the host
> (`cargo build --release`) before running, or the script silently uses a stale
> binary.

## SNP genotype concordance (= 1 − GT-mismatch / TP)

| coverage | ours OLD | ours NEW | freebayes | `1/1→0/1` overcalls OLD→NEW |
|---|---|---|---|---|
| 5× | 83.6% | **94.6%** | 93.4% | 214 → **8** |
| 10× | 97.8% | **98.9%** | 98.3% | 33 → 6 |
| 15× | 99.6% | **99.7%** | 99.2% | 7 → 4 |
| 30× | 99.8% | **99.9%** | 99.4% | 5 → 2 |

The single-sample low-coverage het over-call is fixed: the `1/1→0/1` transition
collapses from 214 to 8 at 5×, matching the prototype's prediction
(83.6 → 94.6 vs the prototype's 94.5). Ours now beats freebayes at every depth.
The residual 5× mismatch is now `0/1→1/1` (70) — true hets sequenced all-alt at
low depth, the same sampling-limited error freebayes has (75); it is not a prior
artefact.

## Precision / recall / FP — unchanged

The genotype prior moves genotypes among *already-emitted* variants (het ↔
hom-alt, both emit the ALT), so the allele-level metrics are identical to the old
prior:

| cov | metric | ours OLD | ours NEW |
|---|---|---|---|
| 5× | TP / FP / FN | 1455 / 13 / 606 | 1455 / 13 / 606 |
| 30× | TP / FP / FN | 2033 / 14 / 28 | 2034 / 14 / 27 |

Precision, recall, and FP are byte-identical (±1 at 30× from pre-existing QUAL
non-determinism). Ours' SNP **F1 beats freebayes at every depth** (5×: 0.825 vs
0.738; 30×: 0.990 vs 0.962) on much higher recall at ~equal precision.

## Indels — unchanged at depth (as diagnosed)

| cov | ours OLD | ours NEW | freebayes |
|---|---|---|---|
| 5× | 24.6% mism | 18.9% | 12.4% |
| 30× | 22.5% | 21.5% | 3.0% |

Indels are also biallelic-diploid, so the prior touches them — it helps mildly at
5× (few reads → prior matters) but the depth-independent `1/1→0/1` plateau is
untouched, because it is driven by inflated REF read support at repeats (a
genuine ~20 % ref VAF the likelihood reads as het), which no genotype prior can
override. This confirms the two-root-cause split from the original investigation:
the indel fix is upstream allele-assignment, separate work.

## Conclusion

The new prior is a strict improvement on SNP genotyping (5× concordance
+11 points, beats freebayes everywhere) at **zero precision/recall/FP cost**, and
does no harm to indels. Proceed to step 3 — remove the now-superseded HWE(p̂)
genotype-prior path. (Caveat for step 3: the HWE path still serves
multiallelic / non-diploid records, so a full removal requires generalizing the
SFS marginalization to arbitrary `(ploidy, n_alleles)` — arch §7.3.)

Reproduce: `tmp/gt_confusion_new.py`, `tmp/prec_recall.py`.
