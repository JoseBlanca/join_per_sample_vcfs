# Hidden-paralog filter — R1 data-first validation on tomato2

**Date:** 2026-07-01. **Branch:** `tomato2-paralog-filter`. **Plan step:** R1
(the checkpoint after the pure statistics core, Q1–Q5).

## What R1 does

Runs the pure statistics pieces end-to-end on the **real tomato2 data** and
checks them against the data itself (not against the Python prototype
bit-for-bit — the production model deliberately diverges: per-sample σ₀,
folded-SFS `p` prior, het-rate `F`; arch Premise 5):

- **Rust harness** `examples/paralog_score_parity.rs` — reads the 59 per-sample
  `.psp` files (coverage histogram + het counts re-derived from the body, so
  `callable_positions` from P1 is present) and the cohort VCF (`AD` + `GT`);
  fits Q2's `SingleCopyCoverageModel` per sample; forms Q4's `F` from the het
  rate and the cohort `Hexp`; measures each locus's per-sample relative copy
  number from the `.psp` window depth; scores with Q3; calibrates π + the FDR
  curve with Q5; dumps per-locus `(chrom, pos, LR, post, qval)` and per-sample
  `(σ₀, obs_het, F)` to TSVs.
- **Python companion** `benchmarks/tomato2/src/paralog_score_parity.py` (uv) —
  does the checks Rust can't (no parquet reader): the loose Python↔Rust LR
  correlation against the prototype's `paralog_lr.parquet`, the F Spearman
  against the prototype's cohort F (recomputed from the VCF), and the
  flagged-set paralog profile against `snp_features.parquet`.

## Results — every check met

| Check | Target | Measured |
|---|---|---|
| Coverage σ₀ (median across samples) | ≈ 0.28 | **0.282** (range 0.208–0.329) |
| Mode/median guard on real samples | passes | **58/59 fit**; 1 rejected (SRS1839376, ratio 0.261 — a genuinely degenerate sample) |
| Single-copy peak anchor | `relative_copy_number = 1.0` | 1.0 by construction (mode anchor) |
| Python↔Rust LR correlation | ≥ 0.98 (loose) | **Pearson 0.9931**, Spearman 0.987, clip(−40,60) Pearson 0.988, on 268 537 shared loci |
| Prior π (paralog fraction) | ≈ 9 % | **8.88 %** (prototype's own π on its LRs: 11.07 %) |
| F Spearman vs prototype cohort F | ≈ 0.86 | **0.851** (Rust F range [0.00, 0.91] vs prototype [0.00, 0.97]) |
| FDR-flagged set = paralog profile | elevated het + coverage excess | **flagged** obs_het 0.115, mean_rel_cov **1.31×**, fis −0.00; **kept** obs_het 0.043, mean_rel_cov 1.04×, fis +0.27 |

Locus sets match the prototype almost exactly (268 537 Rust vs 268 546
prototype), confirming the biallelic-SNP + `≥ 20`-usable-sample filtering lines
up.

## One real bug the data caught: the `Hexp` scale

The first run gave F ≈ 0.99 for *every* sample (F Spearman 0.34) and π 14.45 %.
Cause: `Hexp` was computed as the mean of `2pq` over **variant sites** (≈ 0.092),
while `Hobs = n_het / callable_positions` is a rate over **all callable
positions** (≈ 0.0005) — different scales, so `F = 1 − Hobs/Hexp` saturated.
Fix (matching arch Premise 3): `Hexp` is the per-callable-position expected het,
`Σ2pq / callable` (median callable across the cohort as the reference scalar),
= 24 798 / 30.5 M = **0.000812**. With that, F varies correctly and π drops to
the expected 8.88 %. This is exactly the kind of scale error the data-first
validation exists to catch — the isolated unit tests could not, because they
never see the two quantities' real magnitudes together. (The production S1
wiring must form `Hexp` on the same per-position scale.)

## Verdict

The maths is trustworthy on real data before any var-calling wiring: the
coverage model anchors at 1.0 with σ₀ ≈ 0.28, `F` tracks the prototype
(Spearman 0.85), the LR correlates with the prototype at Pearson 0.99, π lands
at ≈ 9 %, and the FDR-flagged set carries the hidden-paralog signature (het
excess + coverage excess). Milestone R checkpoint reached.

## Reproduce

```text
# Rust harness (container):
./scripts/dev.sh cargo run --release --example paralog_score_parity -- \
    --psp-dir benchmarks/tomato2/psp_files \
    --vcf     benchmarks/tomato2/results/cohort.vcf.gz \
    --out-lr      tmp/paralog_parity_lr.tsv \
    --out-samples tmp/paralog_parity_samples.tsv
# Python companion (host, uv):
uv run benchmarks/tomato2/src/paralog_score_parity.py
```
