# M7 inline-scoring — tomato2 validation (M8)

**Date:** 2026-07-02
**Branch:** `tomato2-paralog-filter`
**Scope:** validate the M7 change (score the hidden-paralog LR once in the caller
worker → single write pass) on the 59-sample tomato2 cohort — drop-profile,
wall, RSS, byte-identity vs the pre-M7 flow (M6) — and record the real σ₀ on
regenerated sliding-window `.psp`.

## What was run

All in the dev container (Apple `container`, **8 cpu / 16 GB**, reference
genome bind-mounted read-only). Regenerated all 59 `.psp` with the current
binary (`pileup`, 200 kb × 160 random regions = 32 Mb BED) — the on-disk `.psp`
predated M5/M6 (no windowed columns, v2 summary), so a fresh build was required
and is the correct current-format artefact going forward.

Three cohort `var-calling` runs on the same 59 `.psp`, all with the
research-retention flags (`--min-qual 0 --no-allele-balance-filter`, DUST on),
`--threads 16`:

| config | binary | wall | peak RSS | data records |
|---|---|---|---|---|
| **M7 filter-on** | `2827eb3` (score-in-worker, single write pass) | **44.9 s** | 350 MB | 262,539 |
| M7 filter-off | `2827eb3` `--no-paralog-filter` | 28.1 s | 400 MB | 283,348 |
| **M6 filter-on** | `dded5fd` (serial two-pass rescoring) | **178.4 s** | 399 MB | 262,539 |

Peak RSS = `/proc/<pid>/status` `VmHWM` (monotonic) polled at 0.5 s; wall =
timestamps. **Caveat:** container-capped and 16 threads on 8 physical cores
(oversubscribed) — the absolute wall is **not** the prod 32-core target; the
**M6:M7 ratio** (same hardware, same inputs) is the trustworthy figure.

## Findings

### Byte-identity (the correctness gate) — PASS

`diff <(zcat m7_on) <(zcat m6_on)` differs on **one line only**: the
`##commandline` header, because the recorded argv embeds the binary path
(`target-container/release/pop_var_caller` vs the copied `pvc_m6`). Every one of
the **262,539 data records** is byte-identical, and the paralog provenance is
identical to six digits:

```
##paralogFilter=target_fdr=0.0100;pi=0.100521;lr_cut=3.8500;em_converged=true   # M7 and M6, verbatim
```

So M7 drops the exact same loci with the exact same π and cut as the pre-M7
flow. The byte-identity argument (worker scores the same record it spills, via
the same pure scorer, from bit-identical inputs; order-independent integer-bin
histogram) holds on real 59-sample data, not just in the unit tests.

### Drop-profile — unchanged

20,809 records dropped = **7.34 %** of the 283,348 unfiltered calls
(283,348 − 262,539). π = **0.1005** (EM converged), LR cut 3.85, F = 0.0000
(default outbred). Consistent with the single-individual-reformulation
validation (~6.6–7 % on tomato2); M7 is a pure perf change and leaves the
callset identical (previous section).

### Wall — the M7 payoff

- End-to-end **3.97×** faster: 178.4 s (M6) → 44.9 s (M7).
- The filter's **marginal** cost (over the shared 28.1 s caller baseline) fell
  from **+150.3 s** (M6, the serial calibrate-score + write-recompute) to
  **+16.8 s** (M7, inline scoring on the parallel workers) — an **~8.9×**
  reduction in the added scoring cost. In ratio terms the filter went from
  ~6.3× the baseline wall to ~1.6×.

This is the direct measurement of what M7 set out to fix: the two serial
transcendental scoring passes (the +11.6× regression flagged in the T2 perf
review) are gone; scoring now overlaps the decode-bound producer on the workers
that already own the cores.

### RSS — flat, slightly better

Peak RSS 350 MB (M7 on) vs 399 MB (M6 on) vs 400 MB (filter off): no
regression; M7 is marginally lower (one fewer spill read + no serial re-score
scratch). Memory-flat in variant and sample count as designed. (Absolute values
are small because the run covers 32 Mb of regions, not the whole genome.)

### σ₀ on the regenerated sliding-window `.psp` (owed from M4)

`examples/tomato2_sigma0.rs` fits `SingleCopyCoverageModel` from each `.psp`
summary and prints σ₀ (`single_copy_depth_sd`):

- **59 / 59 samples fit.**
- σ₀ (relative single-copy depth SD): **median 0.270**, mean 0.270, range
  **0.208 – 0.329**.
- single-copy scale (absolute mode depth over the 500 bp windows): ~5–7×.

Full per-sample table: `benchmarks/tomato2/tmp_m8/sigma0.tsv` (regenerable).

## Conclusion

M7 is validated on tomato2: the callset is byte-identical to the pre-M7 flow
(same 20,809 drops, same π/cut), memory is flat, and the two serial scoring
passes are replaced by one parallel scoring for a **3.97× end-to-end speedup**
(filter marginal cost −8.9×) even on the 8-core container. The prod 32-core
number should be re-measured on the Linux dev box, but the ratio establishes the
win. The regenerated 59 `.psp` are now current-format (windowed columns + v3
summary) and σ₀ is recorded.

### Reproduce

```
# regenerate .psp (container, ref mounted):
DEV_EXTRA_MOUNT=$HOME/genomes ./scripts/dev.sh bash -c \
  'POP_VAR_CALLER_BIN=target-container/release/pop_var_caller JOBS=8 \
   bash benchmarks/tomato2/src/build_psp.sh'
# timed comparison (M7 on/off + M6):
DEV_EXTRA_MOUNT=$HOME/genomes ./scripts/dev.sh bash benchmarks/tomato2/tmp_m8/run_m8.sh
# sigma0:
./scripts/dev.sh target-container/release/examples/tomato2_sigma0 benchmarks/tomato2/psp_files/*.psp
```
