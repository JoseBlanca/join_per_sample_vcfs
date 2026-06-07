# Re-architecting the cohort `.psp` → VCF pipeline

**Branch:** `re-architect` (off `main` @ `a71a078`).
**Status:** design largely settled (the §2 architecture, naming, and build
strategy are decided; see the appendix gap list). A few items remain `[OPEN]`
in §8. The ordered, gated phasing is the
[execution plan](re_architecture_execution_plan.md).

## Why we're here

**The primary objective is a clearer architecture.** The cohort
var-calling pipeline's `.psp` → VCF path is correct and fast, but it is
hard to follow. The five *logical* sections it performs (decode → find
cohort-variant positions → compact → genotype → write) do **not** map
onto its runtime structure: they collapse into three fused components
wired by hand-rolled queues, with names ("producer" / "worker" /
"collector") that say nothing about *what* each does. We want a structure
where the runtime stages map to the logical sections, each named for its
job, with explicit hand-offs — so the pipeline reads the way it works.

**Memory or performance gains would be a welcome bonus, not the goal.**
The re-architecture opens the door to several such wins (column-selective
decode, budget-bounded queues), and we'll take them where they fall out
cleanly. But the bar this work is graded against is clarity, not a memory
or wall-time target.

**Method: re-architect first, then measure against `main`.** Once the
restructured pipeline is in place (and producing byte-identical calls —
see the hard constraints, §6), we measure **wall-time and peak RSS**
against `main` at
representative scale to confirm we haven't regressed, and to capture
whatever bonus gains the new structure delivers. The numbers are a
guardrail and a scorecard, not the design driver.

For reference, the `main` baseline we measure against (50 real tomato
samples, whole genome, release `/usr/bin/time -l`; see `tmp/lowmem_measure/`
on the `low-memory-mode` worktree and [[project_worker_bucket_profile]]):

| config | T=4 | T=8 |
|---|---|---|
| default path | 2976 MB | 3140 MB |
| `--low-memory` | 1797 MB | 2127 MB |

Throughput already scales acceptably (T4→T8 ≈ 1.4×; workers are the
bottleneck and the single producer keeps up), so this is explicitly **not**
a speed project. Peak RSS grows with both `N` and `--threads`; a dhat
profile (N=8) attributes ~58% of peak-live bytes to the read/decode side
and ~30% to buffered output records. Those are the levers a clearer
structure can also pull on — described below as bonuses, gated on the same
measurement.

This document covers: (1) a clearer **vocabulary** for the stages
(the clarity core), (2) **column-selective decode** (the largest bonus
lever, on the read side), (3) **budget-bounded queues** (decouple peak
from `threads`), and (4) **crossbeam-channel** to make the wiring
ergonomic. It is deliberately incremental and measurement-gated.

---

## 1. The architecture as it exists today

The pipeline is best read on two levels: the **logical components** — the
distinct jobs the pipeline performs, in data-flow order — and the
**functional components** — the runtime threads/pools those jobs actually
run on. The mismatch between the two (five logical jobs, three functional
components, with the first three fused) is the clarity gap this work
targets.

### 1.1 Logical components (what the pipeline does)

Five logical components, in data-flow order:

1. **decode** — read one sample's `.psp` and decode the columns for a span;
2. **find-variant** — decide which positions are variant **across the
   cohort** (the fan-in / join over all N samples);
3. **compact** — build a block keeping only the variant positions;
4. **genotype** — EM / posterior genotype calling;
5. **write** — VCF writing.

### 1.2 Functional components (how it runs)

These five do **not** map to five runtime stages. They collapse into
**three concurrent functional components**, plus a look-ahead helper,
wired by two bounded hand-offs (all in
[src/var_calling/driver.rs](../../src/var_calling/driver.rs)):

```
              [DUST-ahead pool]  (precomputes low-complexity masks, runs ahead)
                      │ masks (one per covered interval, in genomic order)
                      ▼
 PRODUCER (1 dedicated thread)        WORKERS (W threads)        COLLECTOR (1 thread)
 §1 decode  (rayon par across N)        §4 genotype                §5 write
 §2 find-variant (fold → is_kept)   ▲   (EM/posterior)          ▲  + reorder (BTreeMap by seq_idx)
 §3 compact (keep kept)             │                          │
        │ ReadyBlock                 │ BlockResult              │
        └── BlockQueue (cap = 2W) ───┘  result channel ─────────┘
            (Mutex+Condvar, bounded)    (std::mpsc, unbounded today)
```

The mapping is the heart of the clarity problem: **logical §1+§2+§3 are all
fused into the single `PRODUCER` component**, while §4 and §5 each get their
own. Three functional components carry five logical jobs.

- **`PRODUCER`** (logical §1+§2+§3, **fused**) = `BlockIterator`
  ([`next_block`](../../src/var_calling/driver.rs#L1746) /
  [`produce_block`](../../src/var_calling/driver.rs#L1813)). One dedicated
  thread. §1 decode is parallel across samples (rayon `par_iter`), but §2
  fold and §3 compact are serial. In `--low-memory` it splits into pass-1
  [`summarize_interval`](../../src/var_calling/driver.rs#L1930) (§1+§2,
  streaming per-sample, drop after fold) and pass-2 `materialize_next_chunk`
  (§3, re-read + keep).
- **`WORKERS`** (logical §4) = `process_block` → `run_window`. The only
  embarrassingly-parallel stage.
- **`COLLECTOR`** (logical §5) = the thread running `drive_blocks_parallel`'s
  receive loop — emits in `seq_idx` order via a `BTreeMap` reorder buffer.
- **`DUST-ahead pool`** (look-ahead helper, feeds §2) = a separate worker
  pool precomputing complexity masks ahead of the producer. Not one of the
  three core components — a side helper delivering masks in genomic order.

Vocabulary today: a **chunk** is a `.psp` read unit; a **block**
(`ReadyBlock`) is the materialised columnar output of §3.
See [[project_cohort_block_iterator_design]].

### 1.3 Structures that flow between the logical components

What each logical component hands to the next (the data, not the threads):

```
 .psp on disk
     │  decoded span, one sample at a time
     ▼
 §1 decode ──────────────► §2 find-variant ──────────► §3 compact ──────────► §4 genotype ──────────► §5 write
   SampleColumns             WindowSummary (cohort        ReadyBlock {            BlockResult {           VCF records
   (per sample, per span)    fold) → is_kept: Vec<bool>     seq_idx,               seq_idx,              (bgzf on disk)
                             + chunk_cuts: Vec<u32>          chunk: Materialised      records:
                                  ▲                          Chunk,                  Vec<PosteriorRecord>,
   IntervalDustMask ──────────────┘                          partition,             stats,
   (DUST helper → §2)                                        pre_fetched_ref_bytes  }
                                                           }
```

| edge | structure | what it carries |
|---|---|---|
| disk → §1 | `.psp` blocks (column manifest + compressed columns) | on-disk columnar pileup for one sample |
| §1 → §2 | [`SampleColumns`](../../src/var_calling/columns.rs#L289) | one sample's decoded columns for a span; folded one sample at a time, not retained |
| §2 → §3 | [`WindowSummary`](../../src/var_calling/two_pass.rs#L47) → `is_kept: Vec<bool>` + `chunk_cuts: Vec<u32>` | the cohort fold across all N samples, reduced to a per-position keep-mask and the chunk's partition cuts |
| §3 → §4 | [`ReadyBlock`](../../src/var_calling/driver.rs#L1406) `{ seq_idx, chunk: `[`MaterialisedChunk`](../../src/var_calling/columns.rs#L654)`, partition: `[`WindowPartition`](../../src/var_calling/partition.rs#L92)`, pre_fetched_ref_bytes }` | N × `SampleColumns` for kept positions only, plus the variant-group partition and pre-fetched REF bytes |
| §4 → §5 | [`BlockResult`](../../src/var_calling/driver.rs#L460) `{ seq_idx, records: Vec<`[`PosteriorRecord`](../../src/var_calling/posterior_engine.rs#L690)`>, stats }` | the genotype calls for one block, tagged with `seq_idx` for the reorder |
| §5 → disk | VCF records | bgzf-compressed VCF in `seq_idx` (genomic) order |
| helper → §2 | [`IntervalDustMask`](../../src/var_calling/driver.rs#L954) | low-complexity mask per covered interval, in genomic order |

Two things to notice. The §1→§2 hand-off is **per-sample and ephemeral** —
`SampleColumns` is folded into the `WindowSummary` and dropped, never all
held at once; this is the join. The §3→§4 hand-off, `ReadyBlock`, is the
**heavy** structure: it owns a full N-sample `MaterialisedChunk` and is the
unit the queues bound. The two collapse into one component (`PRODUCER`)
precisely because §2's keep-mask is needed before §3 can decide what to
materialise.

### Where the bytes live

Every in-flight `ReadyBlock` owns a full **N-sample** `MaterialisedChunk`
(all heavy per-allele columns for the span). The number of in-flight blocks
is `2W` (queue) + `W` (workers) + the result channel. So:

```
peak RSS  ≈  (in-flight blocks)  ×  N  ×  (bytes per position, incl. heavy columns)
                    │                │            │
            scales with threads   scales with   §1 decodes ALL columns up front,
            (2W+W+channel)        sample count   even for the ~96% of positions
                                                 that §3 then drops
```

That product is the whole story: peak ∝ `threads × N × heavy-columns`.

---

## 2. The new architecture

> Bird's-eye module/type/function map (and the settled naming) lives in the
> companion appendix
> [re_architecture_module_outline.md](re_architecture_module_outline.md).

A three-component topology (producer → caller → writer), but the **clarity
win is in the class boundaries** — per-sample work is isolated into objects,
and the cohort-wide logic orchestrates those objects instead of being fused
with them.

### 2.1 Classes (per-sample → cohort → output)

**`SamplePspReader`** — one per sample. Created from `(path, region/bed
segments)` — **no dust** (the producer is the sole dust consumer, §2.3).
Owns the `Read + Seek` handle (today's `ColumnSpanReader`, over the
format-level `psp::PspReader`) and hands out `SamplePspChunk`s as the cohort
advances genomically. Stays `!Send` (per-thread ownership preserved).

**`SamplePspChunk`** — one sample's columns for **one psp segment** (the read
unit). *All per-sample responsibility lives here:*

- holds the segment's compressed columns + manifest; **which columns
  load when is encoded in the methods** — `new()` decodes the fold columns
  (`positions`, `nonref_obs`, `ref_span`) eagerly and caches them; the getters
  decode their column(s) lazily on call;
- data crosses to the producer via **typed getters** —
  `take_seq(keep) -> PerAlleleSeq`, `take_chain_ids(keep) -> PerAlleleChainIds`,
  `take_fixed(keep) -> PerAlleleFixed` — each decoding the `keep` rows only,
  **moving the result out** and freeing the underlying compressed bytes
  (fetched at most once). No generic `Column` type: both ends are statically
  typed, so the producer owns the typing; the producer's per-position merge
  then assembles these per-sample pieces into `CohortPileupRecord`s (the
  typed-getter choice over a generic enum/byte view);
- the §4.1 "decode-some, skip-the-rest" lever rides on the psp format's
  **existing** per-block manifest tags — no parallel `ColumnTag` enum needed;
- **[deferred] low-memory mode** (chosen at `new`): instead of holding
  compressed columns resident, re-read each column from the `.psp` on demand.
  Only materialised when asked for; nothing held speculatively. Built only if
  the new architecture's RAM — **measured against `main` first** — calls for
  it; the count-cap + hold-compressed-segment design may already suffice.

**`CohortPileupRecord`** — one variable cohort position with **every
sample's pileup data at it** (alleles, obs counts, per-allele stats). The
unit the producer emits and the grouper consumes. (This is the old
per-position cohort record — cf. `PerPositionPileups` — revived.)

**`PileupCohortChunk`** — the producer's **work-unit** handed to a caller:
the variable cohort positions of one **safe-gap-bounded span, as a list of
`CohortPileupRecord`s** (+ its `RefSpan`, tagged `chunk_order`). **Not a
columnar block** — the producer builds the records directly, so the caller
groups them without first splitting a block (replaces today's columnar
`MaterialisedChunk`). It's a "chunk" only in the sense of *one ordered piece
of parallel work*; its contents are records.

**`OverlappingPileupRecords`** — the overlapping variant groups after
grouping, **record-based**. Built **and consumed inside the caller** — never
crosses a queue (§2.4). The row shape is **in-tree** as the `#[cfg(test)]`
oracle (`OverlappingVariantGroup` + `MergedRecord`) that the current columnar
path is validated against, so we **promote it to production** — byte-identical
by construction, and grouping isn't perf-critical. *This reverts grouping
only:* the **EM keeps its per-record SoA layout + SIMD backend** (appendix
§D.2) — we do **not** revert the EM math (per-group merge may be row-shape as
long as it feeds the EM SoA likelihood buffers). (Today's grouping is
column-native via `WindowPartition`; this reverts grouping, not the EM.)

**`Variant`** — the final record the VCF writer emits (today's
`PosteriorRecord`).

### 2.2 The cohort variant decision (per-sample data isolated, cohort fold)

"Is this position variant in the cohort" is **intrinsically a fold across
all N samples**, not a per-sample test: positions are unioned, `nonref_obs`
and `ref_span` taken as the cross-sample **max** per position, then
reach-clustered, and a cluster is kept iff its summed max-obs ≥
`min_alt_obs` ([two_pass.rs:94-182](../../src/var_calling/two_pass.rs#L94-L182)).
A per-sample boolean cannot reproduce this byte-identically. So:

- per-sample data stays **encapsulated** in each `SamplePspChunk` (its light
  columns, reached via its methods);
- the cohort producer **holds all N `SamplePspChunk`s** and runs the fold
  over their light columns → the keep-mask (`Vec<bool>`) + the chunk cuts;
- it then gathers the heavy columns for the variable positions via the typed
  getters (`take_seq`/`take_chain_ids`/`take_fixed`) and assembles them into
  the `CohortPileupRecord`s.

This is the split agreed in discussion: per-sample concerns (decode,
light/heavy column access) belong to the object; cohort concerns (the fold,
the keep decision, materialisation) belong to the producer that owns the
objects. **Column-selective decode falls out for free** — the fold touches
only the light columns; heavy columns are decoded only for the ~4 % of
positions that survive.

**Producer algorithm** — it **streams `CohortPileupRecord`s**, sliced into
chunks at safe gaps. Per step:

1. **Read one psp segment per sample, light columns only** (positions, obs
   count, ref-span). Across samples this runs in parallel (read + light-decode
   per sample).
2. **Merge across samples by position** to decide which cohort positions are
   variable, applying the AC / `min_alt_obs` rule — the per-position cohort
   join (the old `per_position_merger`, revived). Then **apply dust**: mark
   dust-masked positions non-variable (after the obs decision, matching
   `main`). The dropped positions leave gaps the **caller's grouper bridges**
   (a group whose reach already spans a dropped position continues past it —
   the byte-identity invariant,
   [partition.rs:33-39](../../src/var_calling/partition.rs#L33-L39)).
3. **Build a `CohortPileupRecord` per variable position.** Decode the heavy
   columns for that position via the readers' `take_*` getters — only the
   ~4 % of positions that vary, never the rest — and assemble each sample's
   data at the position into the record. Each segment is decoded **once**;
   drop it once its records are built.
4. **Accumulate records; cut at a safe gap** (a variant-free stretch wider
   than `max_group_span`, found from the records' reach via
   [`find_block_cut`](../../src/var_calling/loader.rs#L743)) — the chunk
   boundary. Because the cut lands in a gap *between* whole records, nothing
   is ever split and no segment is re-decoded.
5. **Fetch the chunk's REF span** and ship the chunk — its
   `CohortPileupRecord`s + `RefSpan`, tagged `chunk_order` — to a caller.
   **Grouping happens in the caller** (§2.4).

**Memory — the key property (must be preserved).** At the cohort-wide level
the *only* thing held is the **consolidated chunks: variable positions only,
with the AC / `min_alt_obs` filter already applied**. The ~96 % of positions
that carry no cohort variant are **never materialised cohort-wide** — they're
dropped at the per-position merge (step 2), before any record is built. So
cohort RAM scales with *variable* positions × samples × in-flight chunks, not
with all positions. This is where the memory savings live (it preserves
today's AC-pushdown win); the implementation must keep it this way — never
hold a full-coverage cohort structure. Each psp segment is decoded once.
**Column-selective decode falls out for free:** step 1 reads only the light
columns; step 3 decodes the heavy ones only for variable positions.

**Parallelism:** the per-sample read + light-decode (step 1) runs in parallel
across samples; the per-position merge that builds the records (steps 2–3) is
sequential by nature — it *is* the cohort join — on the single producer
thread. This matches today's "parallel fold, serial assemble" shape.

**Why records, not a columnar block** (your call): the next stage groups by
position anyway, so emitting per-position records directly avoids building a
columnar block only to split it — and it sidesteps the mid-segment-cut
problem (records are whole; cuts fall in gaps between them). The cost to watch
is allocation churn from many small records, mitigated by only building the
variable ~4 % and reusing buffers if a profile later calls for it.

### 2.3 DUST sharing

The DUST machinery is fine as it stands (the DUST-ahead pool). Masks are
reference-derived (sample-independent), so the **producer is the sole dust
consumer**: it owns the DUST-ahead pool and applies the masks itself at
step 2 of the producer algorithm (§2.2), downstream of the per-sample fold.

The `SamplePspReader`s **do not need dust at all** — they get the **bed / region
segments** (to know which spans to read), not the mask. This keeps the dust
machinery off the per-sample objects entirely. (Earlier we considered
lending an `Arc<IntervalDustMask>` to each reader; the record-based,
dust-at-the-producer design removes that need.)

### 2.4 Functional sections (the three components)

Grouping is **bundled into the caller**, not the producer (revised from the
earlier "keep grouping in the producer" call — bundling gives fewer queue
crossings and keeps the grouped records worker-local).

- **`PileupCohortChunk` producer** (a `CohortChunkIntegrator` class, or a
  `read_psp_cohort` function exposing its output queue) — owns the N
  `SamplePspReader`s and the DUST-ahead pool; runs the per-position merge
  (§2.2) to **build the `CohortPileupRecord`s** for variable positions, cuts
  chunks at **safe gaps** (wider than `max_group_span`, via
  [`find_block_cut`](../../src/var_calling/loader.rs#L743), so each chunk's
  groups are self-contained), and fetches the chunk's **one contiguous REF
  span** to ship alongside. Does **not** group. Feeds the caller queue.
- **`VariantCaller` workers** (W threads) — the only embarrassingly-parallel
  stage, and the bundle of *several modules*: **group** the chunk's records
  into `OverlappingPileupRecords` + **per-group merge** + **EM/posterior** →
  `Variant`s. REF is sliced per group from the chunk's REF span. The grouped
  records are created and consumed **here** — never crossing a queue.
- **VCF writer** — `Variant` queue → VCF, emitted in genomic (`chunk_order`)
  order via a reorder buffer.

Maps onto the logical components: producer = §1+§2+§3; caller = group + §4;
writer = §5. Note this moves grouping from §3 (producer-side, today) into the
caller — the only departure from today's thread shape, and a clean one (just
a chunk of `CohortPileupRecord`s crosses the queue).

### 2.5 Wiring (crossbeam)

Two bounded `crossbeam-channel` hand-offs:

- **producer → caller**: bounded MPMC (callers pop) — carries the chunk: its
  list of `CohortPileupRecord`s + its `RefSpan`. Replaces the hand-rolled
  `BlockQueue` (Mutex + Condvar). The *grouped* records never travel here —
  grouping happens inside the caller.
- **caller → writer**: bounded, single-consumer — carries `Variant`s.
  Replaces today's unbounded `std::mpsc`.

Bounding both sides gives real back-pressure, so peak tracks queue capacity
rather than scheduling. **Ordering**: callers finish out of order, so the
writer reorders by **`chunk_order`** (today's opaque `seq_idx`, renamed) — a
monotonic per-chunk counter the producer stamps in genomic order, drained via
a `BTreeMap`. The determinism contract is unchanged. **Gapless invariant:**
the writer advances on the *next expected* `chunk_order`, so **every chunk
must yield exactly one `CalledChunk` — even an empty one** (a chunk with kept
positions can still produce zero `Variant`s after grouping/merge); dropping
it would stall the reorder. The caller therefore never silently skips a
chunk. **Decided:** the producer→caller queue uses a simple **count cap**
(crossbeam's built-in bound, like today's `cap = 2W`). A byte-budget
admission rule is a deferred bonus (§4.2) — purely additive (a wrapper around
the same queue), so it can be revisited later with a profile in hand.

---

## 3. Diagnosis: what's intrinsic vs incidental

**Intrinsic (cannot be wired away):** deciding "is this position variant
*in the cohort*" is a **fan-in / join across all N samples**. §2 cannot run
per-sample-independently — there is a barrier where all samples' data for a
span meet. Any design still pays this join. It is why §1–§3 cannot become a
naive per-sample parallel pipeline.

**Incidental (data-model / wiring, fixable):**

- **(A) In-flight memory is bounded by block *count* (`2W`), not bytes.**
  Since a block's size ∝ N, count-bounding makes peak ∝ `N × W`. The cap is
  tied to thread count purely to keep workers fed.
- **(B) §2 only needs *light* columns** (position + obs counts) to compute
  `is_kept`; the *heavy* columns (allele sequences, chain-ids, per-allele
  stats) are only needed in §3 for the ~4% of positions that survive. Today
  the default path decodes **everything** up front, for every sample, then
  throws ~96% away. This is the 58% read bucket.

**Neither throughput nor memory is the objective — clarity is.** T4→T8
gives ~1.4× (workers (§4) are the bottleneck and the single producer keeps
up), so this is not a speed project. The incidental items above are where a
clearer structure can *also* buy memory back — (A) and (B) are bonuses we
take if they fall out cleanly. Any change must hold wall-time and peak RSS
against `main` (the hard constraints, §6); decoding fewer dropped positions
would improve both,
but the design is driven by the cleaner stage mapping, not by these
numbers.

---

## 4. Bonus levers (memory / perf, all measurement-gated)

These are memory/perf wins the §2 architecture *enables* but does not
require — bonuses, not the clarity goal. They compose, can land separately,
and each is behind a measurement gate against `main`. Column-selective
decode (4.1) in particular falls out of the `SamplePspChunk` design almost
for free (§2.2).

### 4.1 Column-selective decode  — the big read lever

**Feasibility (verified):** each `.psp` block carries a per-block **column
manifest** (`ColumnManifestEntry`, tag-ascending, each with a
`compressed_len`), and [`decode_block_payload`](../../src/psp/reader.rs#L1130)
walks columns sequentially from a `Read + Seek` source. A decoder told
*which* columns it needs can **decompress+decode only those and seek past
the rest** (the manifest gives the byte length of each skipped column). The
existing "throwaway sink for unknown columns" path already
consumes-and-discards column bytes; this generalises it to
"skip-known-but-unwanted."

**Light columns** (needed by §2 / `is_kept`): `delta_pos`, `n_alleles`,
`allele_obs_count`, and **`allele_seq_len`** — the latter is what gives the
fold the REF allele's `ref_span` (reach) *without* decoding the sequence
bytes, since the format stores allele **lengths** in a column separate from
the `allele_seq` **bytes** ([writer.rs:861](../../src/psp/writer.rs#L861)).
**Heavy columns** (needed only by §3): the ragged CSR `allele_seq` (bytes) +
`allele_chain_ids` (the 80 MB / 71 MB dhat sites) and the per-allele stat
columns (`q_sum_log`, `fwd`, `placed_left/start`, `mapq_sum`, `mapq_sum_sq`).

This lets §2 run over a fraction of the bytes. **The record-streaming
producer (§2.2) realises it as a single read pass per segment:** decode the
light columns to decide variable positions, then decode the heavy columns
only for those positions (building the records), then drop the segment — each
segment touched once, no re-read from disk (the old `--low-memory` two-pass
re-read is gone). A from-disk re-read survives only as the deferred
low-memory mode (gap #9, §2.1).

### 4.2 Budget-bounded hand-offs — decouple peak from `--threads` *(deferred)*

**Decided: not now.** The hand-offs use a simple **count cap** (§2.5). This
subsection records the future option: replace the count cap with a
**byte-budget** admission rule — hold ≤ `--memory-budget` MB of in-flight
chunks regardless of how many or how large, so peak tracks the budget instead
of `N × threads` and RSS becomes a reproducible knob. It's purely additive (a
wrapper around the same queue), so it can land later without disturbing the
structure. Revisit with a profile in hand — it needs a per-chunk
resident-size estimate and a granularity choice (per-queue vs whole-pipeline),
both better chosen against measured chunk-size distributions than on paper.

### 4.3 crossbeam-channel — ergonomics, not scaling

`crossbeam-channel` gives bounded MPMC + `select!` + better perf than
`std::mpsc`. It would replace the hand-rolled `BlockQueue` (Mutex+Condvar)
and the result channel with uniform bounded channels, and make a
multi-stage / budget-bounded topology cleaner to express. **Explicitly an
ergonomics win** — swapping channels changes no bytes; the result-channel
experiment proved the peak is in the blocks, not the wiring. Worth doing
*because* the re-architecture adds stages and back-pressure edges that are
painful to hand-roll, not as a memory fix in itself.

### 4.4 Naming — settled in §2

The uninformative "producer/worker/collector" names are replaced by the
domain vocabulary defined concretely in §2: `SamplePspReader`,
`SamplePspChunk`, `PileupCohortChunk`, `OverlappingPileupRecords`, `Variant`,
with the producer / `VariantCaller` workers / VCF writer components. The earlier
Decoder/Screener/Assembler/Emitter sketch is superseded by those names.

---

## 5. Topology — decided

A **three-component topology** is the chosen design (§2): one producer
(decode + cohort fold + materialise; no grouping), W `VariantCaller` workers
(group + per-group merge + EM), one VCF writer, two bounded crossbeam
hand-offs. Grouping is bundled into the caller rather than the producer, so
only the compact `PileupCohortChunk` crosses the producer→caller queue and
the grouped records stay worker-local. The full SEDA-style split (separate
Decoder pool / Screener / Assembler / Genotyper / Emitter with a queue
between each) was considered and **rejected for now**: the cohort fold is a
barrier regardless, so extra stages buy cleanliness and cross-interval
overlap but not the dominant memory term. Revisit only if a representative-N
profile shows a stage-overlap stall worth the complexity.

---

## 6. Hard constraints

- **Byte-identity of calls.** GT/GQ/AD/AF/AC/FILTER must be identical to
  `main` for every change (validate with a default-vs-modified VCF diff,
  QUAL excluded — QUAL is non-deterministic since `03e2221`, accepted).
  Column-selective decode must produce the *same* `CohortPileupRecord`s (and
  thus the same calls) for variable positions; budget-bounding and crossbeam
  are pure flow-control.
- **Deterministic emit order** — the `chunk_order` (today's `seq_idx`)
  reorder contract stays, including the gapless "one `CalledChunk` per chunk,
  even if empty" invariant (§2.5).
- **Throughput is measured, not a hard gate.** Wall-time + RSS are tracked at
  N=50 whole-genome, T=4 and T=8, default and `--low-memory`, against `main`
  — but since clarity is the objective, a regression is a *conscious trade*,
  not an automatic block. The main thing to watch is the read-side lifecycle
  change (hold-compressed-chunk vs today's stream-and-drop, §2.2); the EM
  keeps its per-record SoA + SIMD path (appendix §D.2), so it is **not**
  expected to regress. (Byte-identity, above, *is* hard.)
- **Measurement-first** — the N=8 profile *misattributed* the split;
  re-profile at representative N (≥24, ideally 50 on a few whole
  chromosomes) before designing details. Reuse the `[profile.profiling]` +
  `tmp/lowmem_measure/dhat_attrib.py` tooling from the `low-memory-mode`
  worktree.

---

## 7. Build strategy & rollout

**Build in a parallel package, not in place.** The new architecture is built
from the ground up in a fresh **`var_calling_new::`** package (flat module
structure inside, §0 of the appendix), rather than morphing `var_calling::`
incrementally. The split of effort:

- **Rebuild the *structure* from scratch** — the part being re-architected:
  `SamplePspReader` / `SamplePspChunk`, the `CohortChunkIntegrator`, the
  producer→caller→writer wiring, the record-based grouping. This is where the
  clarity win lives.
- **Copy the byte-identity-sensitive numeric kernels verbatim** — the SIMD EM
  (`run_em_columnar`, `compute_log_likelihoods_columnar`, the `posterior_engine`
  / `kernels` math) and the `per_group_merger` numerics. We already chose to
  keep these (appendix gap #5 / §D.2); rewriting proven floating-point math is how
  byte-identity gets lost. Rebuild the skeleton, transplant the organs.
- **Reuse the shared crate directly** — `psp::`, `fasta::`, `bam::`,
  `pileup_record::` are *used*, not duplicated; only `var_calling::`-internal
  code is rebuilt-or-copied.

**The old package is a live byte-identity oracle.** While both exist,
`var_calling::` keeps shipping (the CLI stays on it) and an in-process test
runs *both* pipelines on the same `.psp` cohort and diffs the VCFs (QUAL
excluded). For a rewrite whose success criterion *is* byte-identity (§6),
this continuous live oracle is a far stronger gate than a saved golden file —
**stand it up first.**

**Swap at the end, in one commit:** delete old `var_calling::` + its
now-redundant *unit* tests, **port** the correctness / byte-identity
*integration* tests to the new entry point, rename `var_calling_new::` →
`var_calling::`, repoint the CLI.

**Rollout sketch** (the ordered, gated phasing is the separate
[execution plan](re_architecture_execution_plan.md)):

0. **Re-profile at representative N** + pull the `[profile.profiling]` build +
   `dhat_attrib.py` onto this branch (no production code change).
1. **Skeleton + oracle**: `var_calling_new::` reads (`SamplePspReader` /
   `SamplePspChunk`), `CohortChunkIntegrator` (fold + dust + materialise +
   `RefSpan`), copied kernels, producer/caller/writer wired with crossbeam
   (count cap). Build the old-vs-new VCF-diff oracle test alongside.
2. **Grouping + caller**: record-based `pileup_overlaps` + `em_posterior_calc`
   on the copied kernels; reach byte-identity against the oracle.
3. **Column-selective decode** (4.1) — the read lever, once byte-identical.
4. **Measure vs `main`** (RSS + wall, N=50). *Then* decide the deferred memory
   levers: byte-budget hand-offs (§4.2) and low-memory mode (gap #9) — built
   only if the numbers demand them.
5. **Swap** (the one-commit cutover above).

---

## 8. Open questions

- **[RESOLVED]** Single-pass vs two-pass → **single read pass**: the
  `SamplePspChunk` retains the compressed segment and decodes light-then-heavy
  from it (§2.1 / §4.1). From-disk re-read survives only as deferred
  low-memory (gap #9).
- **[RESOLVED]** Queue admission → simple **count cap** now; byte-budget
  deferred (§4.2), revisit with a profile.
- **[RESOLVED]** Topology Option A vs B → Option A (§5). Stage/type names →
  settled in §2.
- **[RESOLVED]** `SamplePspChunk` granularity / mid-segment cut → **dissolved
  by record-streaming** (§2.2): a `SamplePspChunk` is one segment decoded once
  into records; the chunk cut falls between whole records, nothing is split
  (appendix gap #11).
- **[OPEN]** Does column-selective decode change the DUST coupling (the mask
  is per-position — independent of which columns decode, so likely no, but
  confirm).
- **[DEFERRED — measurement-gated]** Low-memory mode (and whether anything
  like `--low-memory` is even needed) — decide only after measuring the new
  architecture's RAM vs `main`; the count-cap + hold-compressed-segment design
  may already suffice (gap #9, §2.1).

## 9. Non-goals

- No change to the calling math, the EM/posterior engine, or the VCF schema.
- Not chasing the worker/EM scratch bucket as a primary target — it's
  ~7.5% and reused ×W (the N=8 "~40% worker" framing was an artefact).
- Not a speed project; wall-time is a guardrail, not the objective.
