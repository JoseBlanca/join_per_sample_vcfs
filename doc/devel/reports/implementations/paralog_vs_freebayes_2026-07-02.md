# Hidden-paralog filter vs freebayes — GIAB safety + tomato2 FP-class comparison

**Date:** 2026-07-02. **Branch:** `tomato2-paralog-filter`. Follows the T1/T2
report ([paralog_t1_t2_2026-07-02.md](paralog_t1_t2_2026-07-02.md)). Answers
"how does the new cohort paralog filter behave relative to freebayes?" on the two
benchmarks that have a freebayes track: the human GIAB per-sample set and the
tomato2 cohort.

## 1. GIAB (human) — the filter is safe (a structural no-op)

The GIAB freebayes-comparison dashboard evaluates **single-sample** callsets
(HG002/3/4, each on its own confident BED vs its own truth). The paralog filter
is a **cohort** filter: it needs cross-sample allele frequencies (Hexp) and
scores no locus with fewer than `min_samples` (= 20) usable samples. So on a
single-sample GIAB run it scores **0 loci and drops nothing**.

Confirmed empirically — the new binary run on the 50× GIAB `.psp`, filter on vs
`--no-paralog-filter`, identical input:

| sample | records (off) | records (on) | record diffs | samples fit |
|---|---|---|---|---|
| HG002 | 900 | 900 | **0** | 0 |
| HG003 | 725 | 725 | **0** | 0 |
| HG004 | 725 | 725 | **0** | 0 |

The callsets are byte-for-byte identical (variant records), so **every
precision/recall/F1 number vs freebayes is unchanged** — the paralog-filter
column in the dashboard equals the existing "ours" column. Two independent
reasons it does nothing here, both expected: (a) single-sample eval can't run a
cohort filter, and (b) GIAB's confident regions **exclude segmental
duplications**, so the collapsed-paralog class the filter targets is absent from
GIAB by construction (the earlier "validated on tomato, refuted on GIAB"
coverage-signal result). **Takeaway: the on-by-default filter does no harm to the
clean human benchmark.**

> Side note (not a regression): the HG003 **single-pass** run (`--no-paralog-filter`)
> hit a stack overflow at the default thread-stack size on this specific region
> set; it completes with `RUST_MIN_STACK=64M`, and the **two-pass** (filter-on)
> path completed at the default size. This is in the pre-existing single-pass
> writer path, untouched by S6c — flagged for a separate look.

## 2. tomato2 — freebayes emits the FP class the filter removes

The paralog class lives in tomato (seg-dup-rich genome), so this is where the
comparison is meaningful. freebayes v1.3.10 was run on the same 59-sample cohort,
same 160 × 200 kb regions (per-region, `--bam-list` of the 59 CRAMs, concatenated
+ `bcftools norm -m -any`, biallelic SNPs) → **450 045 SNP loci**. Compared
against our filter's decision (from the T1 run: 19 773 dropped, 249 021 kept
biallelic SNPs):

| our verdict | n | freebayes also calls | ALT match | freebayes QUAL ≥ 30 | fb QUAL median / mean |
|---|---|---|---|---|---|
| **dropped (paralog)** | 19 773 | **79.4 %** (15 699) | 99.2 % | **55 % of all drops** | 144 / 750 |
| kept (real) | 249 021 | 86.7 % | 99.5 % | 66 % of kept | 159 / 585 |

**freebayes leaves the paralog class in.** ~79 % of the loci our filter drops are
also called by freebayes, at the same ALT (99 %), and **>half of our entire drop
set are freebayes calls that pass its default QUAL ≥ 30** — i.e. they would appear
in a standard freebayes callset as confident SNPs. Their freebayes QUAL even
*skews high* (mean 750 vs 585 for real variants): the depth inflation that makes
them paralog artifacts also inflates their QUAL, the "confidently wrong"
signature. The ~21 % freebayes doesn't call are mostly its own low-signal /
representation differences.

The figure [`paralog_vs_freebayes.png`](../../../../benchmarks/tomato2/results/paralog_vs_freebayes.png)
shows it in the tomato2 dashboard's coverage × het space:

- **A** — our dropped loci (red) form a distinct **high-coverage + het-excess
  cluster** above the kept density (grey); the paralog signature.
- **B** — freebayes calls 79 % of our drops (55 % at QUAL ≥ 30) — it does not
  filter this class.
- **C** — freebayes QUAL at our drops has a heavy confident tail, not a
  low-QUAL-only distribution.

**Takeaway:** the filter removes a real false-positive class that a standard
caller (freebayes at its default QUAL gate) emits confidently — the class R1
identified and T1 confirmed carries coverage + het excess. There is no truth set
for tomato2, so this is a cross-caller *agreement-on-existence* argument (freebayes
emits them; our coverage model flags them as collapsed paralogs), not a
precision/recall measurement.

## Reproduction

```
# GIAB safety (filter on vs off on existing 50x .psp)
./scripts/dev.sh bash tmp/paralog_t1t2/giab_check.sh

# freebayes on tomato2 cohort (host binary; 160 regions, 8-way)
awk '{print NR" "$1":"$2+1"-"$3}' benchmarks/tomato2/regions_n160_200kb.bed \
  | xargs -P8 -n2 bash tmp/paralog_t1t2/fb_one.sh
# concat+norm -> fb_all.vcf, then:
uv run tmp/paralog_t1t2/fb_compare.py tmp/paralog_t1t2/sum_off.vcf tmp/paralog_t1t2/sum_on.vcf tmp/paralog_t1t2/fb_all.vcf
uv run tmp/paralog_t1t2/fb_figure.py  sum_off.vcf sum_on.vcf fb_all.vcf benchmarks/tomato2/results/paralog_vs_freebayes.png
```
(Scratch scripts under `tmp/` — gitignored; the figure is committed under
`benchmarks/tomato2/results/`.)
