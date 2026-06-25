# Benchmark: SSR genotyping — ours vs HipSTR (tomato, concordance)

**Date:** 2026-06-25
**Branch:** `main`
**Harness:** `benchmarks/ssr_tomato1/` + `benchmarks/lib/{run_hipstr.sh,ssr_concordance.py}`

First end-to-end comparison of our SSR caller against **HipSTR v0.6.2** on a
tomato cohort. No STR truth set exists for tomato, so this measures
**agreement, not accuracy**.

## Method

Both callers genotype the **same locus set**: HipSTR's `--regions` BED is
projected from our `ssr-catalog` by `scripts/catalog_to_hipstr_bed.py`
(coordinate mapping verified against HipSTR `region.cpp` + our VCF), restricted
to the sliced region set. Sample names are forced to match (HipSTR
`--bam-samps` = the CRAM filename base). HipSTR reads CRAM natively; it is
single-threaded (no `--threads`), our caller runs the cohort pre-pass + sweep
on 4 threads.

**Comparison unit:** per-allele base-pair difference from REF (== HipSTR's GB),
as a sorted multiset. This cancels representation offsets and collapses
same-length SNP alleles to a length-0 difference.

**The two callers define an STR allele differently**, and the report classifies
HipSTR loci accordingly so the gap is reported, not penalised:
- **length-poly** — ≥1 ALT differs from REF by a whole motif unit → the
  *comparable* set (what our REPCN model genotypes).
- **non-unit** — a length change that is not a whole motif unit (e.g. 1 bp in a
  5 bp repeat) → our model flags `notPeriodic`.
- **seq-only** — every ALT is REF-length (a SNP inside the repeat) → outside our
  model; our caller is correctly silent (not a miss).

## Cohort

The SNP `tomato1` regions (80 × 100 kb, ~774 incidental SSRs) were too sparse.
`scripts/pick_ssr_regions.py` samples **15,000 catalog SSRs** (seed 42,
per-period balanced 3,000×p2–p6, ch00 excluded), expands each ±500 bp, merges →
`ssr_regions.bed` (14,455 windows, 15.0 Mb). The whole-genome CRAMs were sliced
to that BED on `rick` (`lib/slice_crams_on_rick.sh`) → **51 samples** (one
truncated CRAM removed). Same ~data budget as the SNP slice, ~20× the SSRs.

## Result

| | |
|---|---|
| HipSTR loci genotyped | 11,352 (541 s, 1 thread) |
| Ours loci emitted | 74,408 (51 samples, 222 s, 4 threads) |
| **Comparable set** (HipSTR whole-unit length-poly) | **1,645** |
| ours also emits | 775 (47.1%); 870 absent |
| **Concordance** (cells called by both) | **28,121 / 29,186 = 96.4%** |

HipSTR locus classes: 8,298 mono (we drop, correct), **1,645 length-poly**,
1,117 seq-only SNP (outside our model), 292 non-unit indel (`notPeriodic`).
So only ~14% of HipSTR's calls are repeat-length polymorphisms; the rest is the
allele-definition difference, not a sensitivity gap.

Discordances (3.6%) are overwhelmingly **±1 repeat unit**, skewed long at
dinucleotides: dominant bucket **+2 bp ×509** (ours one dinucleotide unit longer
than HipSTR) vs −2 bp ×74 — a systematic short-period length bias. Everything
else is small whole-unit stutter disagreement.

## Findings

1. **Pre-pass robustness is coverage-bound, not a bug.** On the sparse SNP
   regions, `ssr-call` hard-failed (10/63 samples had no confident pre-pass
   genotype). On the dense 15k-SSR slice it completes on all 51 — the failure
   was too few confident genotypes per sample, now resolved.
2. **Allele-definition mismatch dominates the apparent "miss" rate.** HipSTR
   calls sequence haplotypes (in-repeat SNPs = alleles); we call repeat length.
   A length-based concordance is the fair common ground.
3. **`F_IS ≈ 0.83`** across runs — partly real (tomato is selfing/inbred), but
   the magnitude still warrants a stutter-estimation check.

## Open questions (next)

- **The 870 comparable loci HipSTR calls but we don't emit (53%).** Under-call
  (no-call / filter / length-mono) vs HipSTR over-call on stutter — sample and
  compare evidence.
- **The +2 bp dinucleotide length bias** — systematic, worth a root-cause look.

## Reproduce

```bash
# 1. region set (needs the catalog; ch00 excluded, seeded)
uv run benchmarks/ssr_tomato1/scripts/pick_ssr_regions.py \
    --catalog <catalog> --n 15000 --flank 500 --seed 42 --per-period \
    --exclude-chrom SL4.0ch00 --out benchmarks/ssr_tomato1/ssr_regions.bed
# 2. slice the whole-genome CRAMs on rick with that BED (lib/slice_crams_on_rick.sh),
#    copy *.bench.cram + *.crai into benchmarks/ssr_tomato1/crams/
# 3. run both callers (reuse the cached catalog via CATALOG=...)
OUT_ROOT=.../results_ssr15k CATALOG=<catalog> \
    benchmarks/lib/run_ours_ssr.sh benchmarks/ssr_tomato1/bench.config.sh cohort
OUT_ROOT=.../results_ssr15k CATALOG=<catalog> \
    benchmarks/lib/run_hipstr.sh   benchmarks/ssr_tomato1/bench.config.sh cohort
# 4. concordance
uv run benchmarks/lib/ssr_concordance.py \
    --ours    .../results_ssr15k/ours/cohort/cohort.ssr.vcf \
    --hipstr  .../results_ssr15k/hipstr/cohort.str.vcf.gz

# Diagnostic: dump per-locus .ssr.psp QC (depth + gate drops + observations)
cargo run --release --example ssr_psp_dump -- <file.ssr.psp> [start_min] [start_max]
```

HipSTR build (macOS, one-time): `cd HipSTR && make HipSTR
CPPFLAGS=-I$(brew --prefix)/include LDFLAGS=-L$(brew --prefix)/lib`.
