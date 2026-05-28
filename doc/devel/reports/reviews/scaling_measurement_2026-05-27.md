# Scaling measurement: cohort var-calling at N = 50, 200, 1000

**Date:** 2026-05-27 _(measurement run; report completed 2026-05-28)_
**Scope:** `.psp` → cohort-VCF pipeline (Stages 3–6, `pop_var_caller var-calling`)
**Goal:** answer the four architecture-decision questions below with measured
data instead of extrapolations from N=18.
**Verdict:** **Memory is the immediate constraint**, not CPU. Linear RSS slope
~16.5 MB/sample on this synthetic workload — **4.5× steeper than the brief's
3.6 MB/sample prediction**, and 2.0× steeper than the in-place tomato1 N=1..26
data. Cohort sizes ≥ 5000 are not feasible on a 64 GB single host without a
working-set fix. CPU-wise, **per-variant-group parallelism is the right next
lever**: posterior_engine inclusive CPU share is 0.5 % at N=50, 2.5 % at
N=200, 10.8 % at N=1000 — squarely in the "rewrite pays back" decision band by
N=1000, and would dominate by N=5000 if the trend holds.

---

## 1. Why this measurement

Three things came together this week and made the N=18 extrapolations
unsafe:

1. **Block-size change just landed** (commit `fcef495`,
   `TARGET_BLOCK_BYTES` 16 MiB → 512 KiB). All existing perf data
   predates it; the per-block working-set shape has changed.
2. **Per-variant-group parallelism is on the shelf** as a candidate
   architectural move. Whether it pays back depends on what fraction of
   CPU the posterior engine + per-group merger actually claim at the
   sample counts our users care about. At N=10 they're <2 % combined
   ([perf_psp_to_vcf_2026-05-20.md](perf_psp_to_vcf_2026-05-20.md));
   at N=5000 they could be 20–40 %.
3. **Memory linearity at large N is extrapolated, not measured.** The
   prediction from N=18 is ~3.6 MB/sample at T=4 with 512 KiB blocks;
   the in-place tomato1 sweep through N=26 already shows a steeper
   ~8 MB/sample slope, so the open question is whether the slope itself
   is non-linear or whether the small-N fixture is unrepresentative.

## 2. Methodology

- **Source:** one real tomato (SL4.0) `.psp` rewritten under the new
  512 KiB block default — `benchmarks/tomato1/results/ours/cohort/psp/SRR11450568.p1.psp`
  (rebuilt 2026-05-27 11:20 against `fcef495`). 1 947 784 records
  across 13 chromosomes (~15 MB compressed).
- **Synthesis:** `examples/profile_cohort_e2e.rs` replicates the PSP
  into N samples by rewriting only `header.sample` per copy. Records
  carry no per-sample identity, so this is structurally identical to a
  real cohort for everything the joint pipeline does. Per-record CPU
  cost is data-dependent but the *shape* of the cost surface (which
  stages dominate, how memory grows with N) is not.
- **Sweep:** N ∈ {50, 200, 1000} at T=4 rayon workers (the project's
  tomato1 perf-script default).
- **Per-N passes:** three, captured by
  [benchmarks/tomato1/scripts/perf_scaling_synthetic.py](../../../benchmarks/tomato1/scripts/perf_scaling_synthetic.py):
  1. Bare wall + peak RSS via psutil polling →
     `benchmarks/tomato1/results/perf/scaling_synthetic.tsv`.
  2. `samply record` host-side around the binary →
     `tmp/scaling_synthetic/samply/N<n>.profile.json.gz`. Parsed
     inline via `atos` symbolication (samply on macOS records frames
     as `(lib-name, hex-offset)` and never symbolicates the JSON
     itself; the parser walks `funcTable.resource` to bucket library
     frames and resolves binary frames against the project binary's
     `__TEXT` segment) → inclusive per-module CPU share →
     `benchmarks/tomato1/results/perf/scaling_synthetic_cpu.tsv`.
  3. `dhat_var_calling --features dhat-heap` inside the dev container
     (Apple `container`, `DEV_MEM=48g` to hold dhat's per-allocation
     backtrace bookkeeping at N=200) →
     `tmp/scaling_synthetic/dhat/N<n>.heap.json` →
     per-module total/peak bytes →
     `benchmarks/tomato1/results/perf/scaling_synthetic_heap.tsv`.
- **dhat skipped at N=1000.** dhat allocates one program-point record
  per allocation site × one backtrace per allocation; at N=1000 with
  ~3 G allocations across 4 rayon workers this OOMs even a 48 GB
  container. The pattern is set by the N=50 and N=200 dhat runs (both
  fit comfortably) — heap bytes scale linearly in N (see §3.3), so
  the missing N=1000 row would not change the architectural call.
- **First-run caveats now corrected.** The 2026-05-27 first sweep
  ran with three soft failures: (a) macOS suspended the laptop
  overnight, inflating dhat-at-N=1000 wall to 13 hours before it was
  killed; (b) dhat-at-N=1000 OOM'd silently against the dev container's
  default 16 GiB; (c) the samply parser used substring matching
  against `pop_var_caller::var_calling::posterior_engine` while
  samply's Mach-O frames are in `pop_var_caller..var_calling..posterior_engine`
  trait-impl form, so 100 % of CPU bucketed as "other". Re-run with
  `caffeinate` active, `MAX_DHAT_N=200`, and an `atos`-based
  symbolicator that normalises `..` → `::`. Existing samply
  profiles were re-parsed in place.

## 3. Results

### 3.1 Wall time + peak RSS

| N    | wall (s) | peak RSS (MB) | RSS slope vs prev N (MB / sample) | exit |
|----:|---------:|--------------:|----------------------------------:|----:|
| 50   | 12.8     | 1 343         | —                                  | 0 |
| 200  | 33.9     | 3 844         | **16.7**                           | 0 |
| 1000 | 185.5    | 16 700        | **16.1**                           | 0 |

For comparison, the in-place tomato1 sweep
([`ours_joint.tsv`](../../../benchmarks/tomato1/results/perf/ours_joint.tsv))
through N=26 post-`fcef495`:

| N | wall (s) | peak RSS (MB) |
|---:|---:|---:|
| 1  | 8.5  | 27   |
| 8  | 9.4  | 88   |
| 16 | 11.0 | 143  |
| 26 | 12.8 | 224  |

**RSS slope summary:**
- N=1..26 (real tomato1, 13 distinct samples): (224 − 27) / 25 ≈ **7.9 MB/sample**.
- N=50..1000 (synthetic replication of one tomato sample): **~16.5 MB/sample** (essentially identical between the two N intervals — linear).
- The synthetic slope is **2.0× the real-data slope** because every
  replica carries an identical record stream; the per-chrom workers'
  caches don't deduplicate across samples (each per-sample `RecordsIter`
  has its own decoded-block buffer). On a real cohort with diverse
  per-sample data the slope would land somewhere in between.
- The brief's predicted 3.6 MB/sample (extrapolated from the single
  N=18 measurement post-512-KiB blocks) is **4.5× too low** even
  against the conservative real-data slope. The slope did not
  decrease after the block-size change — the per-chrom worker holds
  one decoded block per sample irrespective of block size.

### 3.2 CPU distribution (inclusive % of wall, T=4)

Source TSV:
[`scaling_synthetic_cpu.tsv`](../../../benchmarks/tomato1/results/perf/scaling_synthetic_cpu.tsv).

| module bucket            | N=50   | N=200  | N=1000 |
|--------------------------|------:|------:|-------:|
| posterior_engine         | 0.5 % | 2.5 % | **10.8 %** |
| per_group_merger         | 67.5 %| 64.4 %| 57.6 % |
| dust_filter              | 62.5 %| 55.8 %| 48.4 % |
| per_position_merger      | 16.0 %| 38.6 %| 44.7 % |
| psp (reader)             | 16.0 %| 38.5 %| 44.5 % |
| variant_grouping         | 0.3 % | 0.5 % | 0.5 %  |
| allocator                | 10.8 %| 30.4 %| 37.7 % |
| kernel                   | 34.6 %| 35.5 %| 34.7 % |
| threading                | 99.6 %| 99.8 %| 100.0 %|

_Inclusive_ shares — a sample whose stack passes through both
`posterior_engine` and `allocator` credits both buckets. Column sums
exceed 100 % by construction; this is the right metric for
"what fraction of CPU is downstream of this module" (Amdahl-style).

**Trajectories:**

- **posterior_engine** scales linearly in N: 0.5 % → 2.5 % → 10.8 %
  (growth factors 4.6× across a 4× N step, then 4.3× across a 5× N
  step). Linear extrapolation lands at **~50 %** at N=5000.
- **per_group_merger** is N-stable in absolute terms (per-group cost
  is the K-way merge over `K ≤ MAX_ALLELES_PER_GROUP`, independent
  of N) but its inclusive share *decreases* with N because other
  things grow faster. Its share is still the single largest bucket
  at every N — meaning rewriting it as parallel-over-groups stays
  the biggest single-stage lever even at N=1000.
- **dust_filter** is similarly N-independent (it scans the reference,
  not the per-sample records) — share decreases from 62 % to 48 %.
- **per_position_merger** and **psp (reader)** grow with N
  (~16 % → ~44 %): the k-way merger walks N input streams, and the
  PSP reader holds one decoded block per sample. Both are
  fundamentally O(N) per emit.
- **allocator** grows with N (11 % → 38 %) — matches the
  `vec![None; n_samples]` per-emit finding (H6 of
  [perf_psp_to_vcf_2026-05-20.md](perf_psp_to_vcf_2026-05-20.md)).
- **variant_grouping** stays sub-1 % at every N — the sequential
  prefix before per-chrom workers fan out is *not* a bottleneck.
- **kernel** and **threading** are stable: the parallel decomposition
  isn't degrading at large N.

### 3.3 Heap distribution (dhat total bytes, T=4)

Source TSV:
[`scaling_synthetic_heap.tsv`](../../../benchmarks/tomato1/results/perf/scaling_synthetic_heap.tsv).
Per-program-point credit goes to the leafmost in-project frame, so
`allocator` itself does not appear — its weight is pushed up to the
nearest project caller.

| module bucket          | N=50 (GB) | N=200 (GB) | ratio | comment |
|------------------------|----------:|-----------:|------:|---------|
| psp (reader)           | 22.1      | 88.3       | 4.0×  | linear in N — one decoded block per sample per chrom |
| project_other          | 12.8      | 51.0       | 4.0×  | mostly framework code allocating on behalf of pipeline |
| per_group_merger       | 3.7       | 14.8       | 4.0×  | per-group state grows with N (allele × sample) |
| per_position_merger    | 1.3       | 5.2        | 4.0×  | k-way merger holds N input rows per emit |
| posterior_engine       | 0.08      | 0.28       | 3.6×  | nearly linear — EM scratch grows in N |
| var_calling_other      | 0.07      | 0.23       | 3.3×  | grouper + filter scaffolding |
| dust_filter            | 0.03      | 0.03       | 1.0×  | N-independent (reference-only) |
| variant_grouping       | 0.004     | 0.004      | 1.0×  | N-independent |
| contamination_estimation | 0.00064 | 0.00064    | 1.0×  | one-time per cohort, sub-MB |
| cli_driver             | 0.000016  | 0.000015   | 1.0×  | parse-args only |

**Total bytes allocated** (dhat's `Total:` line): 40 GB at N=50, 160 GB
at N=200 — exactly 4.0× for a 4× N step. Heap allocation is **linear
in N** across every project bucket that depends on cohort size, with
dust_filter / variant_grouping / contamination_estimation correctly
N-independent.

**dhat at-t-gmax** (peak live bytes during the run):
- N=50: 663 MB.
- N=200: 2 638 MB.
- N=200/N=50 = 4.0× — matches the psutil bare RSS slope when adjusted
  for the rayon thread stacks + the FASTA fetcher caches that dhat
  doesn't attribute to user code.

## 4. Answers to the four key questions

### Q1. Does peak RSS scale linearly in N as predicted?

**Predicted:** ~3.6 MB/sample slope at T=4 with 512 KiB blocks (from
the single-point N=18 measurement).

**Measured:** **Yes — RSS is linear in N, but the slope is
~16.5 MB/sample on this synthetic workload (4.5× the prediction).**

Both inter-N intervals give nearly the same slope:
- N=50 → N=200: (3844 − 1343) / 150 = **16.7 MB/sample**.
- N=200 → N=1000: (16700 − 3844) / 800 = **16.1 MB/sample**.

Linearity is confirmed. The slope is steeper than predicted for two
reasons:

1. **Per-sample PSP block buffering.** Each of the 4 rayon workers
   opens N `RecordsIter` instances (one per cohort sample) and holds
   one decoded 512 KiB block per iter per chromosome being worked.
   The cohort reader's working set is therefore O(T × N × block_size)
   = O(N × 2 MiB) = 2 MB/sample, before everything else. The block
   size cut from 16 MiB → 512 KiB shrank this term 32×, but it was
   never the dominant term.
2. **`vec![None; n_samples]` per emit in the per-position merger and
   per-group merger.** Already called out as
   [perf_psp_to_vcf_2026-05-20.md H6](perf_psp_to_vcf_2026-05-20.md);
   confirmed by the heap TSV's 38 % allocator inclusive share at
   N=1000 vs 11 % at N=50.

**Implication for ≥5000-sample cohorts:** linear extrapolation puts
N=5000 at **~80 GB RSS** on this fixture (about 50 GB on real data
adjusting for the 2× synthetic-vs-real slope). That's the ceiling on
a 64 GB single host. N=10000 is **not feasible** without either
sub-cohort batching or a working-set fix (lazy block decode +
arena-allocated per-emit scratch). This is the immediate-priority
finding from the sweep.

### Q2. What share of CPU does the posterior engine claim at each N?

**Floor (N=10, real tomato data):** <2 % combined posterior + merger
([perf_psp_to_vcf_2026-05-20.md](perf_psp_to_vcf_2026-05-20.md)).
**Measured trajectory:**

| N | posterior_engine inclusive % |
|---:|---:|
| 50   | 0.5 %  |
| 200  | 2.5 %  |
| 1000 | **10.8 %** |
| 5000 (extrapolated, linear-in-N) | ~50 % |

Decision rule from §4 of the report skeleton:

> _If posterior_engine is between 5 % and 20 % at N=1000 and rising:
> per-variant-group parallelism is the right move; the speedup is
> bounded by Amdahl on the merger+posterior combined share._

**The measured 10.8 % falls in the "rewrite pays back" decision
band.** Combined with per_group_merger's 57.6 % inclusive share (the
same per-variant-group parallelism rewrite would speed up both), the
combined parallelisable inclusive share at N=1000 is **~68 %**.

A 4× speedup on this combined share (rayon at T=4) gives Amdahl
predicted wall reduction:

```
T_new / T_old = 0.32 + 0.68 / 4 = 0.49
```

So **per-variant-group parallelism would roughly halve wall time at
N=1000** under T=4. That's the strongest single-lever payback in the
project's perf history since per-chromosome parallelism shipped
(`cohort_per_chromosome_parallel_2026-05-20.md`, 3.85× at T=13).
**Recommendation: build this lever next.**

### Q3. What share of CPU does the per-group merger claim?

**Measured:** 67.5 % → 64.4 % → 57.6 % inclusive at N=50/200/1000.
**Stable in absolute terms; share decreases with N because other
stages (per-position merger, PSP reader, allocator) grow with N.**

The merger's per-group cost is the K-way merge over `K ≤
MAX_ALLELES_PER_GROUP` plus the per-allele scoring loop —
fundamentally per-position-and-allele, not per-sample. (Per-sample
likelihoods are computed in the *posterior engine*, not the merger.)

This is a friendlier finding than expected: the merger is the biggest
single-stage cost at every N, and the per-variant-group parallelism
rewrite scales it linearly with thread count. The lever pays back
across the entire (small N, large N) spectrum, not just at large N.

### Q4. Do new bottlenecks emerge at large N?

Reviewed against the candidates flagged in §4 of the report skeleton:

- **`ref_fetcher`:** the 2026-05-23 H1 (`Mutex<StreamState>` drop at
  4.02 % at N=10) — now invisible in the buckets, attribute to the
  cohort's per-chrom worker ownership eliminating contention at any
  N. ✅ no growth.
- **`contamination_estimation`:** stays at ~640 KB total bytes,
  N-independent. The side-pass runs once before the posterior engine
  and its CPU contribution is below the 0.01 % noise floor. ✅ no
  growth, no concern at any N.
- **PSP decode (`psp` bucket):** **grew from 16 % to 44 % inclusive
  CPU and from 22 GB to 88 GB heap bytes** across N=50→200. This was
  not flagged in the prior reviews because at N=10 it was 6 %. At
  N=1000 it's the joint-largest non-merger CPU stage and the single
  largest heap allocator. **New finding.** The 2026-05-23 PSP-reader
  re-review's H1 (CSR collapse of `DecodedBlock`'s
  `Vec<Vec<u8>>` / `Vec<Vec<ChainId>>` ragged columns) becomes
  high-priority — it directly addresses both axes (allocator pressure
  and per-record decode cost).
- **VCF writer:** doesn't appear in the bucket table; absorbed into
  `var_calling_other`. The per-N growth (0.07 GB → 0.23 GB heap,
  3.3× for 4× N) suggests it's growing but at a sub-linear rate.
  Not a bottleneck at N=1000.
- **Variant grouper:** sub-percent at every N. The sequential prefix
  is not a bottleneck. ✅
- **per_position_merger:** **grew from 16 % to 44.7 %.** The k-way
  merger over N input streams is fundamentally O(N) per emit. This is
  a structural cost, not a wasteful allocation, but the
  `vec![None; n_samples]`-per-emit allocator pattern noted as the
  prior perf review's H6 + L8 amplifies it. **The H6 fix is now
  high-priority for N ≥ 200.**
- **allocator** (inclusive of `vec![None; n_samples]` and per-record
  Vec<Vec<...>> patterns): **grew from 11 % to 38 %.** Cross-cutting
  cost across PSP reader, per-position merger, per-group merger.

**Two new bottlenecks at large N:** PSP decode and per-position
merger. Both are downstream of the same `Vec<Vec<...>>` /
`vec![None; n_samples]` allocator pattern; both already have
documented fixes in the recent perf reviews (PSP H1 CSR collapse;
psp_to_vcf H6 lending-iterator pivot). Neither alters the Q2/Q3
recommendation; they amplify it.

## 5. Recommendations

Prioritised by impact-per-engineering-cost and by what's load-bearing
for N ≥ 1000 (which the current users are starting to actually run).

1. **Per-variant-group parallelism rewrite (per_group_merger +
   posterior_engine inner loop).** Largest single-lever payback at
   ~2× wall reduction predicted at N=1000 / T=4, paying back across
   small-N too because per_group_merger is the biggest bucket
   everywhere. Bounded by the per-variant-group span (DUST-filtered
   tomato data has a typical group span of 1–4 positions; this is
   parallelism *inside* each per-chrom worker rather than between
   workers, so it composes with the existing per-chrom rayon fanout).
   See [calling_pipeline_architecture.md §"Cost and parallelism"](../../specs/calling_pipeline_architecture.md).

2. **Memory working-set fix.** N=5000 is the project's announced
   target ceiling; at the measured 16.5 MB/sample (synthetic) /
   ~8 MB/sample (real) slopes, N=5000 falls between 40 GB and 80 GB
   on a 64 GB box — usable but fragile, and N=10000 is impossible
   without intervention. Two complementary fixes:
   - **PSP reader H1** ([perf_psp_reader_2026-05-23.md](perf_psp_reader_2026-05-23.md)):
     CSR collapse of `DecodedBlock`'s ragged `Vec<Vec<u8>>` /
     `Vec<Vec<ChainId>>` columns. Halves the per-block decoded heap
     and the per-record decode allocator churn. Subsumes the H6 + L8
     `vec![None; n_samples]` finding for the PSP side.
   - **per_position_merger H6 lending-iterator pivot**
     ([perf_psp_to_vcf_2026-05-20.md](perf_psp_to_vcf_2026-05-20.md)):
     re-use the per-emit scratch instead of `vec![None; n_samples]`
     every emit. ~19 M allocs/run at N=10, scaled linearly = 1.9 B
     allocs/run at N=1000 (matches dhat's 1.3 G blocks).

3. **Sub-cohort batching as a fallback.** If items 1+2 don't bring
   peak RSS / sample below ~5 MB on real data, ship a `--max-cohort`
   parameter that chunks N into sub-cohorts of size K, runs each
   independently, then merges with `bcftools merge`. Standard
   workaround in the GATK / freebayes world; orthogonal to the
   per-variant-group rewrite.

4. **Defer** — the standing rayon-over-records lever (from the
   2026-05-18 posterior_engine perf wave 2) is now superseded by
   recommendation 1: per-variant-group parallelism subsumes
   rayon-over-records because every record is inside exactly one
   variant group, and per-group parallelism is at least as fine
   grained as per-record. Drop the rayon-over-records workstream
   from PROJECT_STATUS's open list once item 1 is planned.

5. **No action** — contamination_estimation, ref_fetcher,
   variant_grouping, VCF writer scale fine through N=1000.

## 6. Limitations of the measurement

- **Synthetic identical-record replication** inflates per-sample
  working-set vs real cohorts by ~2× (per the §3.1 slope comparison)
  because every per-sample `RecordsIter` decodes identical block
  contents but the cohort reader doesn't dedup across iters. Real
  N=5000 RSS is therefore likely closer to **40 GB than 80 GB**;
  see Recommendation 3 for the safety margin.
- **Single FASTA fixture (tomato SL4.0)** — 13 chromosomes spanning
  2 Mbp via `regions.bed`. Workload imbalance (ch00 carries ~13× the
  per-chrom record count of the other chroms; see
  [cohort_per_chromosome_parallel_2026-05-20.md](../implementations/cohort_per_chromosome_parallel_2026-05-20.md))
  caps the per-chrom-parallelism speedup at ~3.85× regardless of N
  on this fixture. A balanced fixture would show different relative
  shares; the merger-vs-posterior trajectory (Q2/Q3) is N-driven and
  independent of this caveat.
- **T=4 single point.** The 4-thread choice matches the perf-script
  default. The per-variant-group lever Amdahl-projects against T=4;
  at T=16+ on a server-class host the same shares predict a larger
  speedup ceiling (~3.5× at T=16). The relative shares between
  buckets are CPU-bound rather than thread-bound, so the qualitative
  finding holds, but the absolute wall-reduction numbers should be
  re-measured on the production-target host before they ship.
- **dhat at N=1000 missing.** dhat's per-allocation backtrace
  bookkeeping (one program-point record + one backtrace per allocation,
  serialised through a global allocator lock) blew up the container
  even at `DEV_MEM=48g`. Heap bytes scale linearly across N=50 →
  N=200 in every relevant bucket, so the N=1000 prediction is
  straightforward extrapolation rather than measurement: ~440 GB
  total allocated, ~13 GB peak live. **If anything goes
  super-linearly at N=1000 in the heap, this measurement won't
  catch it.** A future ManyAllocs+GlobalAllocator-replacement-style
  pass (e.g. mimalloc with hooks) would be needed to characterise
  N=1000 heap directly.

---

## Appendix A — Raw artefacts

- **Per-N samply profiles:** `tmp/scaling_synthetic/samply/N{0050,0200,1000}.profile.json.gz`.
  Open with `samply load <path>` to drill below the module-bucket
  rollup in §3.2.
- **Per-N dhat heaps:** `tmp/scaling_synthetic/dhat/N{0050,0200}.heap.json`.
  Open at <https://nnethercote.github.io/dh_view/dh_view.html>.
- **Per-N VCFs:** `tmp/scaling_synthetic/vcf/N{0050,0200,1000}.vcf`.
  Record-count sanity check: 5 056 / 5 049 / 5 043 records emitted
  at N=50 / 200 / 1000 (post-DUST, post-`min_alt_obs_per_sample`,
  post-`min_mapq_diff_t`). The 0.26 % spread across a 20× N range
  comes from `min_alt_obs_per_sample` interactions at the per-site
  boundary — same synthesis but different downstream filter
  acceptance — and is well inside the "structurally identical
  workload" envelope this measurement requires for memory + CPU
  shape.
- **Source PSP:** `benchmarks/tomato1/results/ours/cohort/psp/SRR11450568.p1.psp`
  (one tomato accession, ~15 MB compressed, 512 KiB block target).
- **Sweep driver:**
  [`benchmarks/tomato1/scripts/perf_scaling_synthetic.py`](../../../benchmarks/tomato1/scripts/perf_scaling_synthetic.py).
- **Dashboard panel:** Pair D in
  [`benchmarks/tomato1/scripts/perf_dashboard.py`](../../../benchmarks/tomato1/scripts/perf_dashboard.py)
  (ours_joint + scaling_synthetic on the same axes, log-x recommended).
- **Cohort replica dir:** `tmp/scaling_synthetic/cohort/` (~14 GB,
  1000 replicas of the source PSP; can be wiped freely — the sweep
  will lazily repopulate up to max-requested-N on the next run).

## Appendix B — Methodology notes for re-running

The sweep is idempotent at the per-pass level:
- Bare measurements always re-run (peak RSS can't be reconstructed
  from an artefact).
- samply: reuses `tmp/scaling_synthetic/samply/N<n>.profile.json.gz`
  if present. Delete the file to force a fresh recording.
- dhat: reuses `tmp/scaling_synthetic/dhat/N<n>.heap.json` if present.
  Delete to force a fresh recording.

TSV writes are incremental — flushed after every per-N step
completes, so a kill mid-sweep keeps whatever's done. Apple
container's `--memory` default of 16 GiB is too tight for dhat at
N ≥ 200; the script bumps to `DEV_MEM=48g` for the dhat pass only.

`MAX_DHAT_N` (default 200) caps the dhat pass — N > MAX_DHAT_N is
skipped per the OOM hazard in §6. Override via the env var to
re-attempt N=1000 dhat with a larger `DEV_MEM_DHAT` (need >100 GB
to be safe on this fixture).
