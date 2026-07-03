# Cohort var-calling performance session — 2026-07-03

**Branch:** `cohort-varcalling-perf`. **Scope:** the cohort SNP-calling stage
(`.psp` → VCF, `var-calling`) — *not* pileup. **Goal:** profile it and make it
faster without changing the callset.

**Result: tomato2 (59 real samples), T16 → ~33s to ~22.5s, ≈ 32% faster wall.**
Five commits; every change verified byte-identical (or call-set-identical where a
calibration change was intended), all unit + integration tests green.

## Measurement setup

- Host: macOS / Apple Silicon, **6 performance cores** (+ efficiency cores, 18
  logical), native release build (`target/release/pop_var_caller`).
- Input: 59 real tomato `.psp` (`benchmarks/tomato2/psp_files/`, ~12 GB),
  `--regions regions_n160_200kb.bed --min-qual 0 --no-allele-balance-filter`
  (the research callset — the hidden-paralog filter is **on**, which is the
  default and the interesting case).
- Profiling: `samply record --save-only` on a `--profile profiling` build
  (full debug info, no fat-LTO frame collapsing), symbolicated with `atos`.
  Phase-split by sample timestamp (main pass vs write pass); self-time and
  inclusive-time computed separately (the distinction turned out to matter a lot).
- Guidance: `ai/skills/rust-performance-review` — profile-first, size before
  rewriting, one change per measurement, correctness first.

## Baseline and the shape of the problem

The first profile overturned the working assumption. The EM (which we expected to
be hot) is ~3%. The costs were:

- **Hidden-paralog filter ≈ 30% of wall** (ON 32.9s vs OFF 23.0s).
- **Producer decode ≈ 30%** (zstd + varint + fold).

But an SFS-grid sweep proved the paralog *scoring* is **off the critical path**:
cutting its H1 transcendentals **8×** moved wall ~2%. It runs on the caller
threads, which have slack. So its 30% CPU share was real but not wall-relevant —
the wall is paced by the **producer** ("the producer is the wall floor", as the
pipeline's own docstring says) and, for paralog-ON, by the **serial write pass**.

## What landed

| Commit | Change | Wall effect |
|---|---|---|
| `b43191b` (→ `main`) | Report real thread budget, not the stale "capped by N chromosomes" | reporting only; byte-identical |
| `df67853` | Skip scoring hom-ref loci; calibrate on variable loci only | correctness + ~5% (noise-level) |
| `291c12a` | Stop spilling `window_coverage` (never re-read) | byte-identical; −168 MB spill |
| `fce6d59` | **Parallelise the write-pass VCF formatting** | **~14%** |
| `eb2f34d` | **Incremental producer fold** | **~21%** |

### 1. Chromosome-cap report fix (`b43191b`)
The stderr summary printed `effective_threads = min(pool, n_chromosomes)` with a
"capped by 13 chromosomes" note — a leftover from the retired per-chromosome
architecture. The pipeline distributes an interval schedule across a caller pool;
there is no per-chromosome thread cap (which is why T16 > T12 with only 13
contigs). Now reports `producer_threads + caller_threads`. Landed on `main`.

### 2. Hom-ref skip (`df67853`)
The caller scored a paralog LR for *every* post-EM record, including the ~25% that
carry no ALT in any sample and are dropped downstream as hom-ref. Worse, those
non-variant LRs were folded into the FDR **calibration histogram** — a non-variant
population setting the cut for the real candidates. Return the `NaN` unscored
sentinel for `!is_variant_call` loci: the histogram drops non-finite LRs, so
calibration is now built from variable loci only (the statistically correct
population). Correctness win; ~5% wall (scoring is off the critical path, so the
gain is small and near the noise floor). Output changed as intended (the shifted
cut moves ~4.8k borderline loci between drop buckets; emitted variants unchanged).

### 3. Drop `window_coverage` from the spill (`291c12a`)
The two-pass paralog filter spills every record to disk, then re-reads it to apply
the FDR cut. `ParalogSpillRecord.window_coverage` (~470 B/record × ~350k ≈ 168 MB)
was written but **never read** in production `run_write_pass` — scoring already
consumed it inline in the main pass. Made it a transient local. Byte-identical,
wall-neutral (the spill is dominated by the full 59-sample `PosteriorRecord`, ~3.67
GB; `window_coverage` was ~4%), but honest dead-weight removal.

### 4. Parallel write-pass formatting (`fce6d59`) — ~14%
Phase timing split the paralog-ON run into **main pass ~25s (parallel)** + **write
pass ~7.5s (serial)**. The write pass re-read the 3.67 GB spill, applied the cut,
and wrote the VCF **single-threaded**. Decomposed it further:

| write-pass sub-phase | time |
|---|---|
| spill decode | ~0.7s |
| **VCF text formatting** (`encode` + noodles serialise) | **~5.5s** |
| bgzf compress (serial sink) | ~1.6s |

The formatting is pure per record, so it parallelises. Rebuilt `run_write_pass` as
a reader → format-workers → sink pipeline (dedicated scoped threads): a reader
thread decodes + applies the cut, `W` workers run the downstream filters + format
each survivor to its VCF line off the writer thread (new `DetachedFormatter`,
Arc-shared header/contigs, per-worker genotype-table cache — byte-identical to
`write_record`), and the sink commits pre-formatted bytes to the bgzf sink in
order via `CohortVcfWriter::write_preformatted`. Only the final bgzf sink stays
serial. Write pass 7.8s → ~2.7s. Verified byte-identical VCF **and** identical
run-summary counters.

### 5. Incremental producer fold (`eb2f34d`) — ~21%, the big one
The fold/plan stage's `rebuild_fold` re-folded the **entire growing window**
`[chunk_start, watermark]` on every read-ahead iteration of `fill_to_target` —
O(iterations × window) per chunk. The code flagged it as "the producer's serial
floor". Though only **7.3% of CPU**, it is the **serial spine that paces the whole
producer** (decode and compact run ahead but wait on it), so its redundancy capped
chunk-production throughput.

Replaced it with `extend_fold(w)`: fold only the new watermark window
`(fold_upto, w]` and **append** it, resetting at each cut. Byte-identical because
the fold is per-position integer `max` + position union (associative +
commutative), and each increment's positions are all strictly greater than
everything already folded — so it appends (a subset of the order-independent
union) to the same final fold the rebuild produced. Correctness rests on
`w ≤ watermark`: every sample's data in the window is already buffered, so no later
read can add a position `≤ fold_upto`.

**Wall ~28.5s → ~22.5s (~21%)** — far more than the fold's 7.3% CPU share, because
the limiter was *serialization*, not CPU. Verified byte-identical at fixed config;
calls + genotypes identical across thread count and chunk size (the QUAL FP wobble
across chunk boundaries is the pre-existing non-determinism, calls unchanged).

## The rewrites we did *not* do (measured, refuted)

Profiling **self-time** (not inclusive time) killed three tempting rewrites before
they were built — the single most valuable discipline of the session:

- **SIMD varint decode.** `decode_varint_column` is only **1.0% self-time**; the
  "19%" that motivated it was the *inclusive* decode subtree, which is mostly zstd.
  LEB128 already has a single-byte fast path and most varints are single-byte.
- **A better fold algorithm** (direct-address / bucket, exploiting dense integer
  keys to drop the merge's `log N`). Correct and elegant, but the fold is only
  **2.8% of CPU** after the incremental change — the redundancy, not the
  per-fold algorithm, was the cost. <1% wall available.
- **A faster codec (lz4 vs zstd).** A micro-bench on *real* tomato2 column
  payloads: lz4 is only **1.27× faster to decode for 1.79× bigger files**. The
  data is so compressible (zstd → 14% of raw) that zstd already decodes it at 2.3
  GB/s; lz4's usual 3–4× edge only appears on less-compressible data. At ~9%
  decode self-time that is ~2% wall for +79% disk — a poor trade.

## Can we give decode more threads?

No. The producer **saturates at ~6 threads = the P-core count**. A split sweep
showed more producer threads (P10–14) never beats the default, and rayon already
`par_iter`s the decode across the whole producer pool automatically. zstd
decompression is compute-heavy and wants a *fast* core each; beyond 6 threads the
extras land on efficiency cores (slow) or oversubscribe. The 32% `[kernel/park]`
is the *caller* threads idling (paralog off-path) while the producer is P-core
bound. A producer-lighter split *might* shave a few %, but it is within this
machine's ~10–20% run-to-run noise. The ceiling here is **6 fast cores × zstd
throughput** — a hardware limit, not a threading one. The code already scales with
P-cores, so the 32-core production target would gain further for free.

## Where it stands

The main pass is now a well-balanced, producer-bound stage (fast zstd decode +
cheap incremental fold + compact) where the easy wins are exhausted. Remaining
levers are all hard: a decode/format redesign, or attacking the pipeline's
structural parallelism (off-CPU concurrency analysis via `hotpath` would be the
next diagnostic). This branch is 4 perf commits ahead of `main`.

### Net
**~33s → ~22.5s (~32%)** on tomato2/59 at T16, byte-identical callset, standout
being the incremental fold. Lesson: rank work by **self-time**, and remember that
a small-CPU stage can still be the wall floor if it is the serial gate.
