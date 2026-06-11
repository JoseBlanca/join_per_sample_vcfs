# Variant QUAL distribution vs sequencing depth — ours / GATK / freebayes

**Date:** 2026-06-10 · **Branch:** `qual-analysis` · **Benchmark:**
`benchmarks/human_genome_bottle` (single-sample GIAB HG002, ~5 Mb of
GRCh38 high-confidence regions, native CRAM ≈ 301x).

## Questions

1. The QUAL distribution of false positives vs true positives, for each
   of the three SNP callers.
2. How that TP-vs-FP QUAL separation changes with sequencing depth,
   within each caller.

## Method

The native CRAM is ~301x — far above realistic WGS — so depth was
**created by subsampling** (`samtools view -s`, seed 42, fraction =
target/301.35) to a ladder of global depths: **5, 10, 15, 20, 30, 50,
100, and full (≈301x)**. Each subsample lands within 0.1x of its target
mean over the BED. At every depth all three callers were re-run with the
benchmark's fair-comparison flags (every emission gate = 0; `--no-baq`
and `--no-complexity-filter` for ours; `-stand-call-conf 0` for GATK; no
post-call QUAL cap for freebayes), then each call was scored TP/FP
against the GIAB truth set with `compare_to_truth.sh` (allele
concordance on POS+REF+ALT after left-align + biallelic-split).

Orchestrated by `benchmarks/lib/run_depth_sweep.sh`. ours + freebayes
run on the host; GATK runs in the dev container (`/opt/gatk/gatk`), which
needs CRAM 3.0 — the subsamples are written `--output-fmt cram,version=3.0`
because samtools' default 3.1 is rejected by the container's htsjdk.
`compare_to_truth.sh` was extended to emit per-record locus **DP**
alongside QUAL. Outputs:

- `results/depth_sweep/qual_dist_by_depth.tsv` — one row per call
  (`depth, caller, class, status, qual, dp`), 249,873 rows.
- `results/depth_sweep/accuracy_by_depth.tsv` — precision/recall/F1 per
  (depth, caller, class).

Dashboard: `benchmarks/lib/qual_depth_dashboard.py` (marimo).

## Finding 1 — TP/FP QUAL separation at native depth (SNPs)

Median QUAL and *gateability* = share of FP calls whose QUAL sits below
the median TP QUAL (1.0 ⇒ a QUAL cutoff cleanly removes every typical FP):

| caller | TP med | FP med | nFP | gateability |
|--------|-------:|-------:|----:|------------:|
| freebayes | 4261 | **0** | 8650 | 1.00 |
| ours | 5862 | 150 | 1403 | 0.97 |
| gatk | 5879 | 2461 | 123 | 0.66 |

Two distinct error profiles:

- **freebayes & ours** emit *many* FPs but pile them at low QUAL — they
  are noise the QUAL score correctly distrusts, so a gate recovers
  precision almost for free (freebayes' 8650 FPs are essentially all at
  QUAL≈0; ours' 1403 FPs are mostly < 627).
- **GATK** emits *few* FPs (123) but they are **confident-wrong** —
  high-QUAL, overlapping the TP distribution. A QUAL gate cannot remove
  them; only ~⅔ sit below the typical TP. GATK's precision advantage is
  structural, not gate-recoverable.

## Finding 2 — separation vs depth (SNPs)

Gateability is remarkably **flat across depth** for each caller, even as
the absolute QUAL scale slides ~50× from 5x to 301x:

| depth | freebayes gate | ours gate | gatk gate |
|------:|---------------:|----------:|----------:|
| 5x   | 0.95 | 0.91 | 0.61 |
| 10x  | 0.99 | 0.96 | 0.78 |
| 20x  | 1.00 | 0.97 | 0.77 |
| 30x  | 1.00 | 0.97 | 0.79 |
| 100x | 1.00 | 0.97 | 0.68 |
| 301x | 1.00 | 0.97 | 0.66 |

Takeaways:

- **The QUAL scale is strongly depth-dependent.** TP-median SNP QUAL
  rises roughly linearly with depth (ours 86→5862, gatk 124→5879 from
  5x→301x). **A fixed QUAL threshold does not transfer across depths** —
  any recommended `--min-qual` must be depth-aware (or expressed
  per-read-depth, e.g. QUAL/DP).
- **freebayes' FP flood is depth-driven.** FP count climbs 697 (5x) →
  ~13k (30–50x) while FP QUAL stays ≈0; its low-depth precision is
  actually its *best* (fewer spurious sites to call). Its FPs remain
  perfectly gateable at every depth.
- **ours is the most depth-stable separator** (gate ≈ 0.97 from 10x up),
  with a small, low-QUAL FP set (~1.3k) the score reliably flags.
- **GATK's few FPs stay hard to gate at every depth** (0.61–0.81); at low
  depth its precision lead over ours/freebayes narrows because its TP and
  FP QUAL distributions sit closest together.

Indels (smaller counts, noisier): freebayes/GATK FPs are also mostly
low-QUAL (gate 0.83–0.99), while **ours has the least-gateable indel FPs
at low depth** (gate 0.50 at 5x, 0.77–0.81 elsewhere) — our indel FPs are
relatively higher-QUAL than the other callers'. Worth a closer look if
indel precision at low depth matters.

## Dashboard

`uvx marimo edit --sandbox benchmarks/lib/qual_depth_dashboard.py`

Four views (class + clip/bins/density/log-y controls): TP-vs-FP QUAL
histogram **grid** (depth rows × caller cols); **separation vs depth**
(median TP/FP QUAL with IQR bands + the gateability score); **accuracy
vs depth** (precision/recall/F1); and **best QUAL cutoff vs depth** (the
F1-optimal threshold and the F1 it achieves — shows the depth-dependence
of the optimal gate directly). The companion `comparison_dashboard.py`
covers the single-depth TP/FP QUAL + ROC view.

## Reproduce

```
bash benchmarks/lib/run_depth_sweep.sh benchmarks/human_genome_bottle/bench.config.sh host
DEV_EXTRA_MOUNT=$HOME/genomes ./scripts/dev.sh \
    bash benchmarks/lib/run_depth_sweep.sh benchmarks/human_genome_bottle/bench.config.sh gatk
bash benchmarks/lib/run_depth_sweep.sh benchmarks/human_genome_bottle/bench.config.sh compare
bash benchmarks/lib/run_depth_sweep.sh benchmarks/human_genome_bottle/bench.config.sh merge
```

(`full` depth reuses the existing canonical full-depth VCFs under
`results/{ours,gatk,freebayes}/`.)
