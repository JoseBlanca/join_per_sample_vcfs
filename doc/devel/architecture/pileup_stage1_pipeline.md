# Stage 1 (pileup) pipeline — current architecture & perf opportunities

**As of `pileup-perf` @ `5e02144` (2026-06-12).** Living doc: the
topology after the de-barrier (commit 2). Use it to reason about where
the remaining wall-time goes and which levers are worth pulling.

`pop_var_caller pileup` turns one sample's BAM/CRAM(s) into a `.psp`
artefact. At `--threads >= 4` it runs the de-barriered pipeline below;
below 4 it falls back to the inline bulk-synchronous `BaqStream` (cheaper
coordination when there aren't enough cores for the three roles to
overlap).

## Current topology (`--threads >= 4`)

```
            ┌──────────────────────────┐
            │  SOURCE  (1 thread)       │   serial, per input file
            │  AlignmentMergedReader    │
   BAM/CRAM │   • noodles decode (bgzf/ │
   ───────▶ │     cram codec)           │   ← decompress is single-threaded
            │   • k-way merge (by coord)│
            │   • at-locus dedup        │
            │   • order check           │
            │   • min_read_length       │
            │   • RecordBuf→MappedRead   │   ← copies seq/qual/cigar (alloc-heavy)
            │   → single-contig packets  │
            └────────────┬─────────────┘
                         │  RawPacket{seq, reads}        bounded crossbeam
                         │  (≈1024 reads, one contig)    cap = depth·M
                         ▼
        ┌───────────────────────────────────────┐
        │  READ-PROC POOL  (M = n_threads)        │   EMBARRASSINGLY PARALLEL
        │  run_worker × M, process_read per read: │
        │    • G2  cigar_is_bad                    │
        │    • F3  left_align_indels  (raw ref)    │   raw-byte RawContigRefCache
        │    • F1  mismatch_fraction  (raw ref)    │   (per worker, per contig)
        │    • BAQ HMM                (upper ref)   │   ManualEvictChromRefFetcher
        │  → ResultPacket{seq, survivors, counts}  │   (per worker, per contig)
        └───────────────────────┬─────────────────┘
                         │  ResultPacket (out of order)   bounded crossbeam
                         ▼
            ┌──────────────────────────┐
            │  CONSUMER  (1 thread)     │   serial
            │  ReorderReadIter          │   ← BTreeMap reorder by packet seq
            │   → restores coord order  │
            │  PileupWalker             │   ← STATE MACHINE, serial within contig
            │   • active-read set        │     (open-record widen, mate pairing,
            │   • open-record table      │      chain-ids, depth cap, closure)
            │  PspWriter                │   ← block build + zstd compress + IO
            │   → .psp                  │     (shares this thread with walker)
            └──────────────────────────┘
```

**Back-pressure** is the bounded channel depth (`QUEUE_DEPTH_PER_WORKER ·
n_workers`): the source blocks when workers fall behind, workers block
when the walker falls behind. **Byte-identity** comes from the reorder
buffer restoring global coordinate order before the serial walker — the
output is independent of worker count.

## Hardware note (the measurement machine)

The dev laptop is **6 performance cores + 12 efficiency cores** (reports
as "18 logical"). The E-cores are much slower, so effective parallelism
is bounded by ~6 fast cores, *not* 18. Two consequences the perf numbers
below all reflect:

- The all-cores default (`--threads` omitted → 18) is the **worst**
  setting: it spawns 18 BAQ workers that spill onto E-cores and contend.
  Worker-count sweep at `--threads 18`: **W=4 → 17.5 s (best)**, W=6 →
  17.7, W=8 → 18.0, W=18 → 18.5 (user-time climbs 71 s→84 s = wasted
  E-core work). Optimum ≈ the fast-core count.
- The pipeline is **CPU-bound** here: at T4, user≈63 s on 4 cores →
  63/4 ≈ 16 s hard floor, and we sit at ~17 s. Rearranging threads
  cannot beat total-work ÷ fast-cores.

## Where the wall-time goes

Steady-state wall ≈ **max(SOURCE, WORKERS/M, WALKER, WRITER)** — a *max*,
not a sum, because the roles run on different threads and overlap. With
enough workers (`M ≳ 4` saturates) the `WORKERS/M` term drops below the
serial stages, so the floor is:

```
  floor ≈ max( SOURCE , WALKER+WRITER-on-one-thread )
```

Measured (tomato SRR7279540, 18-core native release):

| threads | wall  | notes                                            |
|--------:|------:|--------------------------------------------------|
| 1       | 61.0s | inline path                                      |
| 2       | 42.5s | inline path                                      |
| 4       | 17.1s | pipeline — **3.46×**, ~hits the floor            |
| 8       | 17.8s | pipeline — flat (floor-bound, not worker-bound)  |
| 18      | 18.7s | pipeline — slight oversubscription (M=18)        |

The flat T4→T18 curve is the tell: **past ~4 workers the bottleneck is
the serial stages (SOURCE and WALKER), not the parallel fold.** The
~17 s floor is `max(source, walker+writer)`.

### Measured stage breakdown (#0, 2026-06-12, T4)

Env-gated per-stage timers (`PVC_STAGE_TIMERS=1`), full-genome tomato,
wall 18.06 s with the timer overhead:

```
SOURCE:   read(decode+merge+dedup+convert) = 6.7s   send-blocked = 8.2s
CONSUMER: walker = 9.7s   writer = 5.4s   recv-blocked = 1.0s
```

**The CONSUMER thread binds.** It is busy **15.1 s** (walker 9.7 s +
writer 5.4 s, run *serially* on one thread) and blocked only 1.0 s. The
SOURCE is busy only 6.7 s and spends 8.2 s **blocked waiting to hand
off** — decode outruns the consumer. So:

- The **writer is ~36 % of the consumer thread**, serialized behind the
  walker → **splitting it onto its own thread (#2) overlaps 5.4 s with
  the walker's 9.7 s**, dropping the consumer floor from 15.1 s → ~9.7 s.
  Predicted wall **~17 s → ~10 s** (another ~1.7×).
- After #2 the **walker (9.7 s)** is the floor — only contig-parallel
  (#4) breaks it.
- Decode is **not** the bottleneck (half-idle), which is exactly why #1
  (conversion→workers) changed nothing.

## Opportunity map

Ordered by leverage × (1/effort). "Floor" = does it lower the
`max(...)` floor, i.e. go below 17 s.

| # | Lever | Where | Floor? | Effort | Risk / notes |
|--:|-------|-------|:------:|:------:|--------------|
| 0 | ~~Measure stage breakdown~~ | probe | — | XS | **DONE 2026-06-12** (see above): the **consumer binds** (walker 9.7 s + writer 5.4 s serial); source is half-idle. Points straight at #2, then #4. |
| 1 | ~~Conversion → workers~~ | SOURCE→pool | **no — TRIED** | S | **Tried 2026-06-12 and reverted: byte-identical but zero wall-time change (T4 17.1→17.25).** Confirmed by #0: the source is half-idle (8.2 s send-blocked), so speeding it up can't help. Dead. |
| 2 | ~~Split writer onto its own thread~~ | CONSUMER | **no — TRIED twice** | S–M | **Reverted.** No gain BAQ-on (T4 17.4→17.7). Re-tested under `--no-baq` to *free* the fast cores (workers go cheap) — still no gain (16.03→16.15). So it would **not** pay off on bigger hardware either, as built: the per-record channel (~7.6 M `PileupRecord` hand-offs, writer flushing in zstd bursts) costs about what the overlap saves. Also learned: `--no-baq` is only 1.4 s faster than BAQ-on → **BAQ is fully hidden behind the workers; the wall is bound by the serial WALKER**, not BAQ or the writer. |
| 3a | **Cap worker count** | pipeline | small (~6%) | XS | Measured: `n_workers = n_threads` oversubscribes (18 workers on 6 fast cores). Capping to ~fast-core-count gives the best wall (W=4 17.5 s vs W=18 18.5 s) and fixes the "all-cores default = worst" footgun. Modest but the cheapest real win. **Don't hardcode a laptop-tuned 4** — needs a portable heuristic (or just a saner default than all-logical-cores). |
| 3b | ~~Thread-count `2n`→`n`~~ | pipeline | thread hygiene | XS | **DONE 2026-06-12.** `configure_rayon_pool(n)` spawned an idle n-thread rayon pool the pipeline never uses → `--threads 4` ran **10 threads** (4 workers + 4 idle rayon + source + main). Now rayon is sized only for the inline path; `--threads 4` runs **5** (= n + c). Byte-identical. |
| 3 | **Tuning sweep** | pipeline knobs | partial | S | `n_workers` cap (T18 18.7 vs T4 17.1 = oversubscription); inline↔pipeline threshold (T2/T3); packet size (1024); queue depth (2). Modest, cheap. |
| 4 | **Contig-parallel lanes** | whole pipeline | **YES (big)** | L | Run each contig through its own source+pool+walker lane; parallelizes decode *and* the walker across contigs → floor drops from `max(source,walker)` toward `total / largest-contig` (~5–6× for 12 tomato chromosomes). Needs an **ordered-merge writer** (reorder by contig, assign byte offsets at write — the `.psp` sequential-write concern) and a **chain-id prefix-sum** so ids match the serial numbering (byte-identity), or accept semantic-but-not-byte-identical ids. |
| 5 | **Multi-file parallel decode** | SOURCE | yes (multi-file only) | M | One decode thread per input BAM/CRAM, k-way merged from per-file channels. Helps lane/replicate merges; **does nothing for a single CRAM** (the benchmark). |
| 6 | **Within-file decode parallelism** | SOURCE | yes | L | CRAM slices / bgzf blocks are independently decodable, but noodles doesn't expose this cleanly. Largely subsumed by #4 (per-contig indexed queries decode in parallel). Low priority. |
| 7 | **Walker hot-path micro-opt** | WALKER | yes (small) | M | Fold accumulation / cursor queries; already optimized (lazy CIGAR, fold cache). Only worth it if #0 shows the walker binds *and* profiling finds a specific hotspot. |
| 8 | **Memory (RSS) check** | whole pipeline | — | S | Not yet measured for the pipeline (bounded in-flight packets + per-worker fetchers likely raise RSS vs main). Measure before merge if memory matters; trade queue depth / packet size if needed. |

## Recommended sequencing

1. **#0 measure** the source/walker/writer split (decides everything).
2. Cheap squeezes on the single-lane pipeline that the measurement
   points at: **#1 (conversion→workers)** and/or **#2 (writer thread)**,
   plus **#3 tuning**. These keep the simple topology and stay
   byte-identical with low risk.
3. **#4 contig-parallel** is the only path *well* below 17 s, but it is
   the most work and carries the `.psp` ordered-merge + chain-id
   complexity. Worth it only if single-sample latency is a real goal
   (few/deep samples, interactive, per-sample SLA) rather than the
   cross-sample pool model (which already saturates cores by running
   pileup `--threads 1` across samples).
