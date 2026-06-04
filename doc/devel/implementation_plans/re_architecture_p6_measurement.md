# Phase 6 — measurement vs `main` (results)

Companion to the [execution plan](re_architecture_execution_plan.md) §P6. Run on
the host (macOS, `/usr/bin/time -l`), release build, **T=1 single-threaded**
(the new pipeline has no parallel topology yet — crossbeam deferred), whole
genome (no `--regions`), 50 real tomato `.psp` samples
(`/Users/jose/devel/pop_var_caller/tmp/aligned_psp/SRR*.p1.psp`).

The new pipeline is driven through the same CLI via a temporary `PVC_NEW_PIPELINE`
env hook in `run_var_calling` (removed at the P7 swap).

## Peak RSS + wall (T=1, whole genome)

| N  | old RSS | new RSS | old wall | new wall | calls identical¹ | records |
|----|---------|---------|----------|----------|------------------|---------|
| 2  | 118 MB  | 238 MB  | 2.7 s    | 2.4 s    | ✅ yes            | 33 463  |
| 8  | 372 MB  | 446 MB  | 6.6 s    | 6.5 s    | ✅ yes            | 64 287  |
| 24 | 1041 MB | 1065 MB | 17.4 s   | 18.5 s   | ✅ yes            | 133 399 |
| 50 | 2200 MB | **1857 MB** | 38.8 s | 53.9 s | ✅ yes            | 186 259 |

¹ "calls identical" = GT/GQ/AD/AF/AC/FILTER byte-identical (QUAL excluded, the
hard constraint §6) — verified with `--min-qual 0` so the QUAL filter can't
perturb the record set (see below). At N=24 the diff is **0** lines over 152 394
records; at N=50 the diff is **0** over 204 118 records.

## Conclusions

1. **Calls are byte-identical to `main` at full scale.** N=50 whole genome,
   204 118 records, zero differences on the contract columns. The
   re-architecture is faithful on real data, not just the oracle fixtures.

2. **Memory is competitive — P5 is NOT needed.** The new/old RSS ratio *shrinks*
   with N (2.0× → 1.2× → 1.02× → **0.84×**): the small-N overhead is fixed and
   amortizes, and at N=50 the new pipeline uses *less* RSS than `main` (1857 vs
   2200 MB). The hoped column-selective-decode memory lever (P5) is therefore
   **dropped** — the record-streaming architecture already meets the memory goal
   (it preserves the AC-pushdown: only variable positions are held cohort-wide).

3. **Wall is comparable at T=1 up to N=24; ~1.4× slower at N=50** (single-thread).
   This is the expected cost of having no parallel topology yet. `main` ships at
   T=4/T=8; the new pipeline must gain the **bounded-crossbeam producer→caller→
   writer topology** to be wall-competitive in production. This is now the
   justified next step (it was "optional perf" pre-measurement).

## The one divergence: borderline-QUAL records under the `--min-qual` filter

With the **default** `--min-qual 30`, N=24 and N=50 differ by a *handful* of
records (1 at N=24). Every such record sits at QUAL ≈ 30: e.g.
`SL4.0ch10:15343880 T→C` has QUAL **30.42 in `main`** but **29.24 in the new
pipeline** — opposite sides of the threshold, so `main` keeps it and the new
pipeline drops it.

This is **the accepted QUAL difference, not a calling bug**:

- QUAL is a *product over per-sample hom-ref posteriors* and is explicitly
  **excluded** from byte-identity (§6; non-deterministic since `03e2221`).
- The underlying calls are identical: with `--min-qual 0` the two VCFs match
  exactly (the record is present in both, only its QUAL differs).
- The per-sample posteriors differ only by float-ordering noise between the
  **row EM** the new caller uses (`run_em_for_record`) and the **columnar EM**
  `main` uses (`run_em_columnar`) — too small to move the integer GQ or the
  argmax genotype, but the QUAL product over 24–50 samples accumulates it to
  ~1 phred, enough to cross the filter for a record already at the boundary.

**Consequence for P7 (swap):** after the swap the new pipeline *is* production,
so its QUAL is self-consistent; there is no ongoing A/B. The migration artifact
is that a few records within ~1 phred of `--min-qual` may appear/disappear vs
`main`. Options: (a) accept (calls unchanged; QUAL already non-deterministic),
(b) pursue bit-exact QUAL parity with the columnar EM (not required by §6, and
would mean reviving the columnar EM input path we deliberately dropped). **(a)
is recommended** — it is within the letter and spirit of the §6 contract.
(User decision: accepted.)

## Addendum — parallel topology (crossbeam), N=50 whole genome

After P6, the bounded-crossbeam producer→caller→writer topology was added
(producer on the main thread; W caller threads share `&VariantCaller`; writer
thread reorders by `chunk_order`). Re-measured at production thread counts
(`--min-qual 0` to compare calls):

| T | old wall | new wall | old RSS | new RSS | calls identical |
|---|----------|----------|---------|---------|-----------------|
| 4 | 11.4 s   | 27.2 s   | 2982 MB | 2099 MB | ✅ 204 118       |
| 8 | 7.8 s    | 27.3 s   | 2998 MB | 2184 MB | ✅ 204 118       |

- **Byte-identity holds under parallelism** (the writer's `chunk_order` reorder
  makes the output independent of worker completion order).
- **Memory now beats `main`** — ~2.1 GB vs ~3.0 GB (≈ −30 %). `main`'s peak grows
  with threads (peak ∝ in-flight blocks ∝ W); the new pipeline's is bounded by
  the queue cap + the single producer.
- **Wall does NOT scale with threads** — flat ~27 s at T=4 and T=8, ~3.5× slower
  than `main` at T=8. **The producer is the serial bottleneck**, not the callers
  (parallelizing the callers left wall unchanged). In `main` the producer's
  per-sample decode is rayon-parallel across samples and the *workers* are the
  bottleneck; here the producer is fully serial — serial per-sample segment
  decode + the record-building allocation churn the design (§2.2) flagged as
  "the cost to watch". Closing the gap needs **intra-producer parallelism**
  (rayon across samples for the per-sample read/decode), which is the next perf
  decision.
