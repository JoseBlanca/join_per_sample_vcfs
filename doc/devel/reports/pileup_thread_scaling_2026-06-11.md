# Pileup (Stage 1) thread-scaling analysis — 2026-06-11

**Branch:** `pileup-perf`. **Question:** does `pop_var_caller pileup`
(BAM/CRAM → `.psp`) use `--threads` effectively *within a single
sample*?

**Answer: no.** Single-sample pileup peaks at ~1.8× on an 18-core host
and then *regresses* — it is structurally limited to a thin
parallel-only stage (BAQ) wrapped in a serial pipeline with no
stage overlap.

## Measurement

- Host: macOS (Apple silicon), 18 logical cores, native release build
  (`target/release/pop_var_caller`, not the container binary — the
  Apple-container VM caps CPUs and would distort the curve).
- Input: one tomato sample,
  `benchmarks/tomato1/crams/SRR7279540.p1.bench.cram` (~80 MB CRAM),
  reference `S_lycopersicum_chromosomes.4.00.fa`. Whole-genome (no
  `--regions`). BAQ on (default).
- One run per thread count (~30–60 s each — long enough to be stable;
  variance dominated by the serial floor, not noise).

| Threads | real | user | speedup | parallel efficiency | maxRSS |
|--------:|-----:|-----:|--------:|--------------------:|-------:|
| 1  | 59.2 s | 55.2 s | 1.00× | 100% | 538 MB |
| 2  | 41.4 s | 57.1 s | 1.43× |  72% | 565 MB |
| 4  | 32.7 s | 59.5 s | **1.81×** | 45% | 695 MB |
| 8  | 34.4 s | 67.4 s | 1.72× |  22% | 644 MB |
| 18 | 41.8 s | 76.8 s | 1.42× |   8% | — |

Two things stand out:

1. **Peak is 1.81× at T=4, then it regresses.** T=18 (1.42×) is no
   better than T=2. On 18 cores we leave ~16 of them idle or wasted.
2. **`user` time climbs 55 s → 77 s** as threads increase. Extra
   threads burn *more* total CPU while *worsening* wall time — the
   signature of synchronization/overhead with no useful parallel work
   to absorb it.

Correctness is unaffected: `.psp` bodies are **byte-identical across
all thread counts** (the only differing bytes are the
`created = …Z` timestamp in the TOML header). Threading does not
change results.

## Root cause (from the code)

Single-sample Stage 1 is a 4-stage chain, only one stage of which is
parallel, and the stages do not overlap:

```
reader (BAM/CRAM decode)  →  BaqStream (BAQ-cap)  →  walker  →  PspWriter
       SERIAL                  PARALLEL (rayon)       SERIAL      SERIAL
```

1. **Only BAQ is parallel.** Decode, the walker (pileup accumulation,
   indel left-alignment, chain-ids, open-record closure), and the
   writer all run single-threaded.
   `src/pop_var_caller/stage1_pipeline.rs`. With `--no-baq` the whole
   command is single-threaded (passthrough, no rayon at all).

2. **The parallel stage is bulk-synchronous, not pipelined.**
   `src/pileup/per_sample/baq_stream.rs` `refill_batch` (≈ lines
   190–268): serially pull ~1024 reads from the decoder → BAQ them in
   parallel (`par_drain` / `map_init`) → walker drains the batch
   serially. While BAQ runs, decode and the walker sit idle; while the
   walker drains, all BAQ workers sit idle. The three phases are
   strictly sequential per batch.

3. **Regions are processed serially.** `src/pop_var_caller/cli.rs`
   `run_pileup` loops `for region in region_set.iter()` into one shared
   writer — no cross-region parallelism.

Amdahl: solving the T=4 point (1.81×) gives a parallel fraction of
~0.60, i.e. ~40% of wall time is the serial decode+walker+writer floor.
That alone caps speedup near ~1.8×; the bulk-synchronous stalls then
make T>4 *worse* rather than merely flat.

## Context: is this even the intended usage?

The project's current benchmark stance runs pileup **single-threaded
per sample and parallelizes across samples** via a FIFO process pool —
see `benchmarks/tomato1/scripts/perf_ours_pileup_build.py`, which calls
`pileup … --threads 1`. Under that model, intra-sample `--threads` is
effectively a BAQ-only knob and the poor scaling above is expected, not
a regression.

So the finding splits into two questions for the user:

- **(a)** Is single-sample throughput a goal at all? If the deployment
  model is always "many samples, one core each," intra-sample threading
  may not be worth investing in — the cores are already saturated across
  samples. The cross-sample path already scales well by construction.
- **(b)** If single-sample wall time *does* matter (few samples, deep
  coverage, interactive runs, or a per-sample latency SLA), then the
  serial floor + bulk-synchronous BAQ is the thing to fix.

## Candidate fixes (if (b) matters) — not yet implemented

Ordered roughly by leverage vs. effort:

1. **Pipeline the stages (decode ∥ BAQ ∥ walker) instead of
   bulk-syncing.** Replace the batch barrier with bounded channels so
   decode, BAQ workers, and the walker run concurrently. Ideal wall
   ≈ max(decode, BAQ/N, walker, writer) instead of their sum. This
   directly attacks the rising-user-time regression and likely recovers
   most of the lost cores. Highest leverage. (Mirrors the
   producer-staged-pipeline work already done on the cohort side.)

2. **Region/contig-parallel walker.** The walker is independent per
   disjoint region (the region-driven path already clamps and
   concatenates — `drive_region_into_writer`). Running N regions
   through N independent walker+writer lanes parallelizes the *serial*
   40%, which pipelining alone cannot. Needs an ordered merge into the
   single output `.psp` (or parallel shards + concat). Biggest win for
   the serial floor, more involved (writer ordering, chain-id reach
   across region seams — bounded by `MAX_RECORD_SPAN`).

3. **Tune the BAQ batch size.** 1024 is a guess; the barrier cost
   scales with batch granularity. Cheap experiment, small ceiling — at
   best it shifts the plateau, doesn't remove it.

## Repro

```sh
REF=/Users/jose/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa
CRAM=benchmarks/tomato1/crams/SRR7279540.p1.bench.cram
for T in 1 2 4 8 18; do
  /usr/bin/time -l ./target/release/pop_var_caller pileup \
    --reference $REF --output tmp/pileup-perf/t$T.psp --threads $T $CRAM
done
```
