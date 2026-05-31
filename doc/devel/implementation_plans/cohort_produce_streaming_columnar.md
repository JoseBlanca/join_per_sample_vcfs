# Cohort `var-calling` — streaming columnar produce (span-addressable reader + fold/compact)

Implementation plan for reworking the **produce** half of the
`cohort_block` driver — the path from the per-sample `.psp` files to a
materialised block of variable variants. Produce is now the limiting
factor on both wall time (producer-bound past N≈24) and, especially,
memory.

This plan does **not** touch the consume half (per-group EM/posterior +
emit), which already runs in parallel across blocks
([`drive_blocks_parallel`](../../../src/var_calling/cohort_block/driver.rs)).
It also does not touch the PSP→VCF *outputs*: every stage is gated
byte-identical against the current driver.

## Motivating measurement

Real 26-sample tomato cohort, `var-calling` PSP→VCF, branch
`cohort-within-chromosome-parallel` vs `main`, 8 threads:

| caller | N=26 wall | N=26 peak RSS |
|---|---:|---:|
| main (per-chromosome streaming) | ~8.7 s | ~440 MB |
| branch (chunk-batch produce + parallel consume) | ~19.3 s | ~3960 MB |

The branch is ~2× slower and ~9× heavier. Memory grows ~135–155 MB per
added sample (linear in N; pre-existing, documented in
[`examples/dhat_var_calling.rs`](../../../examples/dhat_var_calling.rs)),
and is **invariant** to the parallel look-ahead depth, the per-block
variant target, and the chunk genomic span — ruling those out as the
driver.

Per-phase wall-time breakdown (serial path, N=26, additive):

| phase | time | share | scales with N? | parallel today? |
|---|---:|---:|---|---|
| md5 reference verify (once, upfront) | 1.2 s | 2 % | no | yes |
| **PSP load** (decode + fold + compact + boundary) | 29.7 s | 51 % | yes | across samples |
| **DUST mask** | 9.7 s | 16 % | **no** (≈const ∝ span) | **no — serial on producer** |
| partition + REF prefetch | 0.1 s | 0 % | — | — |
| block math (EM/posterior) | 18.7 s | 32 % | yes (~N^1.8) | **yes** (across blocks) |
| emit (filters + VCF write) | 0.5 s | 1 % | mild | no |

DUST is ≈constant regardless of N (9.70 s at N=8, 9.69 s at N=26): it is
a function of the genomic span analysed, not the sample count, and it
runs **serially on the single producer thread**.

## The two structural problems in today's produce

See [the produce-side code](../../../src/var_calling/cohort_block/loader.rs)
and [`BlockIterator::produce_block`](../../../src/var_calling/cohort_block/driver.rs).

1. **Batch-then-compact = all samples resident over a whole span.**
   `load_chunk_from_iters` pulls a *span* of records for **all N samples
   at once** into `raw_per_sample`, folds the entire accumulation to
   find variable positions, and **only compacts at the very end**. When
   variants are sparse the span grows (by doubling) to hit the
   ~1024-variant target, so the raw working set — held for every sample
   simultaneously — balloons. This is the memory gap vs `main`, which
   k-way-merges the readers **one position at a time** and never holds a
   big all-samples window.

2. **Row round-trip + re-fold waste.**
   - The reader yields per-record `PileupRecord` row objects (each
     allocating an allele `Vec`, sequence bytes, …) which `push_record`
     immediately decomposes back into `SampleColumns`. The on-disk PSP
     is *already columnar* and the destination is columnar — the row
     object is pure overhead: **columnar → row → columnar**.
   - Each span doubling re-runs the fold over the *entire* accumulation
     from scratch (`rebuild_position_union_and_is_kept`,
     O(positions × samples) every time).
   - DUST runs inline on the producer's critical path (~10 s serial).

## Target architecture

Replace batch-then-compact with a **span-addressable columnar reader**
feeding a **streaming fold+compact** loop, converging the producer onto
`main`'s memory profile while still emitting the ~1024-variant blocks the
parallel consumer needs. Three independent, separately-gateable changes:

1. **Span-addressable columnar reader** — a per-sample object that
   decodes PSP blocks into columns and serves an arbitrary requested
   chromosome span, hiding block boundaries from the producer. No
   `PileupRecord` row object on the hot path. This is where the
   cross-sample misalignment problem is encapsulated.
2. **Streaming fold+compact loop** — the producer asks every reader for
   a common span, folds + compacts it, and accumulates to the variant
   target. Raw working set drops from "whole span × all samples" to
   "≈one block per reader + a group-sized straddle".
3. **DUST-ahead queue** — a background thread pre-computes DUST masks
   for the covered intervals (DUST depends only on the reference +
   region) and feeds them to the producer through a bounded queue,
   off the critical path.

The end state of produce is essentially `main`'s k-way streaming merge,
but accumulating into blocks for parallel consumption: `main`'s memory,
the branch's parallel-consumable output.

### The misalignment problem and where it lives

Variability is a **cohort** property: position `p` cannot be classified
until *every* sample has been read past `p`. And a variant **group** can
reach across many positions (an MNP/DEL), so a group cannot be finalised
until the read frontier has passed its entire reach. The wrinkle: PSP
blocks are **not** aligned across samples — sample A's next block may be
1000–2000 while B's is 1500–2500 — so "read one block per sample" leaves
the samples at different positions.

The design choice is to **abstract block boundaries away from the
producer loop** and push them into the per-sample reader. To the
producer, a PSP block is an invisible implementation detail of the
reader; it only ever asks for a **genomic span**. The reader exposes two
operations:

- `peek_next_span()` → the chromosome span it can serve next, computed
  from the block index (each block carries `first_pos`/`last_pos`,
  decoded at open) plus any buffered tail. **Free** — no I/O, no
  decompress.
- `read_span(start, end)` → the **columns** over exactly `[start, end)`,
  for *any* span. The reader decodes whatever blocks overlap and buffers
  whatever it didn't serve. The API does **not** require the span to
  match a block boundary.

Alignment is therefore a **producer convention, not an API constraint**:
the producer *chooses* to ask for the span the reader just reported via
`peek_next_span()`, which is what keeps the readers in step on the happy
path. Nothing breaks if it asks for something else — the reader simply
slices/buffers. The producer drives a simple loop, with the cross-sample
frontier made implicit as the **consensus span**:

1. **Peek** every reader; take `W = min over readers of
   peek_next_span().end` — the furthest position *all* samples can serve
   cheaply.
2. **Read** `read_span(cursor, W)` from every reader, in parallel — by
   construction all samples are now aligned over `[cursor, W)`.
3. **Fold + compact** every variant group whose *entire reach* is `< W`:
   classify its positions (variable + homref-inside-reach kept, pure-REF
   dropped) and append the kept columns to the output block. A group
   reaching past `W` stays in a small loader-side straddle until a later
   span passes its reach. Advance `cursor = W`.
4. **Stop** when the compacted block holds ≥ `target_variants` *and* a
   group has just closed — that point is a clean boundary *by
   construction* (see below). Emit the block; the reader tails + group
   straddle carry forward to the next.

**Happy path:** when all blocks are aligned, every `peek_next_span().end`
is equal, so `W` is that common end and every `read_span` returns one
whole block — zero slicing, zero buffering. Misalignment only triggers
the reader's internal buffer, never the producer.

**Clean cuts come for free — the old safe-cut/carryover machinery goes
away.** Because the loop only ever finalises *complete* groups, a block
boundary always lands on a clean group edge; there is no mid-group tail
to reserve and no failed-boundary case to retry. So this rewrite
**removes** today's `finalise_chunk_boundaries` (safe-end search), the
`NoSafeGap` span-doubling retry, and the explicit `carryover` /
`carryover_snapshot` reserve. The only state crossing block boundaries
is the reader tail buffers + the loader's open-group straddle, both of
which persist naturally across blocks within a covered interval (reset
at interval / chromosome boundaries). This is a deliberate
simplification, not a reuse — we build the streaming loader fresh rather
than fit the old batch machinery.

### Two buffer levels, both small

The abstraction cleanly separates two kinds of buffering that the
batch loader conflated into one giant span:

- **Reader-level (block alignment):** decoded-but-unserved tail, ≤ one
  block per reader. Lives inside each span-addressable reader.
- **Loader-level (group reach):** a group straddling `W` held until its
  reach is passed. Bounded by `max_group_span` — group-sized, *tiny*.

Total resident is `N × one block + one group straddle` instead of the
old "grown span × N samples" — flat in N rather than linear in
(span × N). That is the memory fix made precise.

```
   reader[0] ─ peek/read_span ┐                       ┌─ kept cols ─┐
   reader[1] ─ peek/read_span ┤  cols over [cursor,W) │             │
     ...      (decode + tail  ├─► (all samples aligned ├─ fold+      ├─► ReadyBlock
   reader[N-1]  buffer here)  ┘   by construction)     │  compact    │  (~1024 vars,
                          ▲                            │  (groups<W) │  clean cut)
       DUST-ahead queue ──┘ (mask for [cursor,W))      └─ group straddle (tiny)
```

## Staged implementation

Each stage lands independently and is gated byte-identical (header-
stripped md5) at N=1/4/8/26 on the real tomato cohort, serial and at
4/8 threads, before the next begins. Follow the project's
incremental-step discipline: types → impl → verify → commit, pausing
between stages.

### Stage 1 — span-addressable columnar reader

**What.** A per-sample object wrapping the `PspReader` that serves
chromosome spans as columns: `peek_next_span()` (free, from the block
index + buffered tail) and `read_span(start, end)` (decode the
overlapping blocks, return the `[start, end)` slice as columns, buffer
the unserved tail). No intermediate `PileupRecord` on the hot path; the
cross-sample misalignment is hidden here.

**Why.** Removes a large per-record allocation stream
(`PileupRecord` + allele `Vec`s + seq bytes per record) — the
**columnar → row → columnar** round-trip — *and* localizes the
block-alignment buffering in one cohesive place instead of smearing it
across the producer loop.

**Shape.** The PSP block already decodes to column buffers
([reader.rs](../../../src/psp/reader.rs)); `SampleColumns`
([columns.rs](../../../src/var_calling/cohort_block/columns.rs)) is the same
column family. The work:
- A `SampleColumns`-level **append over a position range** from a decoded
  block's columns. Where the decoded layout lines up with `SampleColumns`
  it's a direct column append (the happy path, no per-record allocation);
  where it doesn't it's a per-column transcode — still no `PileupRecord`
  row object, so still strictly fewer allocations than today's
  round-trip. We don't optimise for the mismatch case beyond "no row
  object".
- A one-block tail buffer + the `peek`/`read_span` API. `peek` reads
  `first_pos`/`last_pos` off the index (no decode); only `read_span`
  decodes, and only the blocks the requested span touches.

**Files.** `src/psp/reader.rs` (columnar span server + tail buffer),
`columns.rs` (append over a range), `loader.rs` (consume the new API in
place of the per-record pull).

**Test infra.** Add a **deliberately block-misaligned** multi-sample PSP
fixture (staggered block boundaries across samples) — the existing
fixtures are aligned, so the reader's slice/buffer path would otherwise
be untested. This fixture also backs the Stage 2 byte-identity gate.

**Gate.** Byte-identical (it still feeds the *same* records over the same
spans into the existing batch fold — this stage does not yet change the
loop). Expect a measurable drop in allocator churn (dhat) and a small
wall improvement on load.

### Stage 2 — streaming fold+compact loop (the memory fix)

**What.** Replace the span-batch loop in `load_chunk_from_iters` with the
consensus-span loop above: peek all readers → `W = min next span end` →
`read_span(cursor, W)` from all (parallel) → fold + compact groups whose
reach is `< W` → advance → repeat until the variant target is met and a
group has just closed → emit. Build this fresh — do **not** retrofit the
batch machinery (per the optimal-code-over-effort call).

**Why.** Collapses the raw working set from "grown span × all samples"
to "≈one block per reader + a group-sized straddle" — the structural
memory fix — and removes the re-fold-on-every-doubling waste (each span
is folded once).

**Shape.**
- The loader-side **group straddle**: positions of a group reaching past
  `W`, held until a later span passes its reach. Bounded by
  `max_group_span`; reused across iterations (scratch discipline).
- An incremental fold over each newly-served `[cursor, W)` slice that
  finalises only groups fully below `W`, producing the identical
  kept-set + group boundaries the batch fold produces.
- **Grouping is still a block-close pass**: once the block's variable
  variants are compacted, run the existing grouping/partition step over
  that small, clean chunk — don't build groups incrementally during the
  fold.

**Removed by this stage.** The streaming loop makes every block boundary
a clean group edge by construction, so the following go away entirely:
`finalise_chunk_boundaries` (safe-end search) and its scratch, the
`NoSafeGap` span-doubling retry in `produce_block`, and the explicit
`carryover` / `carryover_snapshot` reserve + its snapshot/restore. Their
job is subsumed by the persistent reader tails + the group straddle.

**Files.** `loader.rs` (the consensus-span loop + group straddle,
reshaped `ChunkLoadScratch`); `driver.rs` (simplified `produce_block` —
no retry/snapshot); **delete** `chunk_boundaries.rs` (or reduce it to
whatever, if anything, the block-close still needs).

**Gate.** Byte-identical; expect peak RSS to drop toward `main`'s
profile (target: roughly flat-in-N — `N × one block` + the straddle +
the compacted block, not ~135 MB/sample). Re-run the main-vs-branch RSS
table.

**Risk.** This is the substantive change. The consensus-span / group-
straddle bookkeeping and the "never finalise an open group" invariant
are where byte-identity can break — design the buffer lifecycle and the
finalisation predicate on paper first (see Open questions).

### Stage 3 — DUST-ahead queue (serial DUST off the critical path)

**What.** A background thread DUSTs the covered intervals (known up front
from the block indices) in genomic order and pushes `(region → mask)`
into a bounded queue. The producer slices out the mask for its current
span instead of computing it inline.

**Why.** DUST is ~10 s of serial work on the producer thread, constant
in N — after Stage 2 it is the largest serial floor in produce. It
depends only on the reference + region (nothing about PSP data or the
fold), so it precomputes trivially.

**Shape.** Reuse `sdust_mask_for_span`
([dust_filter.rs](../../../src/var_calling/dust_filter.rs))
with its existing halo/floor rules (so masks are byte-identical to the
inline path). A queue (not a per-block future): the producer's cut spans
are dynamic, so pre-DUST the *covered intervals* and have the producer
slice; order-preserving and simple. Bound the queue so the DUST-ahead
thread can't run unboundedly far ahead.

**Files.** `driver.rs` (spawn the DUST-ahead thread + consume the
queue), `dust_filter.rs` (unchanged math).

**Gate.** Byte-identical (DUST is deterministic per region); expect the
producer's serial DUST time to disappear behind the fold/compact.

## Byte-identity discipline

The contract (unchanged from the prior phases): header-stripped md5 of
the VCF body is identical to the current driver, at N=1/4/8/26 on the
real tomato cohort, serial and at 4/8 threads, **after each stage**. The
baselines live in `tmp/` during development. The span-addressable reader
is byte-identical iff `read_span` returns exactly the records the
per-record pull returned for the same span (the slice/clamp is exact).
The streaming fold is byte-identical iff it (a) finalises a group only
when its full reach is below `W`, never finalising an open group, and
(b) classifies positions with the same variable / homref-inside-reach /
pure-REF rule as the batch fold. DUST-ahead is byte-identical iff it
keeps the same halo/floor rules as the inline `sdust_mask_for_span`.

(Note: the branch already diverges from `main`'s output by a known,
pre-existing filter-order difference — that is *out of scope* here and
must not be conflated with these gates, which are against the *branch's
own* current output.)

## Open questions to resolve before Stage 2 code

- **`peek_next_span` semantics with a buffered tail.** A reader holding
  a tail must report the end of its *contiguous* servable data (tail end
  if buffered, else next on-disk block end) so `W = min(...)` stays
  correct. Pin this contract.
- **`W` choice / advance policy.** `min` of the readers' next span ends
  is the natural consensus (keeps each `read_span` to ≈one block of new
  decode). Confirm it can't stall (e.g. one reader at end-of-chromosome
  while others continue) and how end-of-interval is signalled.
- **Reader tail-buffer ownership + recycling.** Where the ≤one-block
  decoded tail lives in the reader and how it reuses its column
  allocations across `read_span` calls (scratch discipline).
- **Incremental fold structure.** Whether to keep the
  position-union / has-variant / max-ref-span arrays incremental, or
  rebuild only over the newly-served `[cursor, W)` slice each step;
  and where the group straddle lives.

(Resolved during planning: the safe-cut / `carryover` / `NoSafeGap`
machinery is **removed**, not reused — clean cuts are by construction;
the reader↔`SampleColumns` append is direct on the happy path and a
transcode otherwise, both without a row object; alignment is a producer
convention, not an API constraint; grouping stays a block-close pass.)

## Out of scope / non-goals

- The consume half (EM/posterior/emit) — already parallel across blocks.
- The pileup stage (BAM/CRAM → PSP).
- Closing the known filter-order output divergence vs `main`.
- Changing the block size target or the parallel-consume look-ahead
  (measured irrelevant to memory).
- Multi-producer / per-interval parallelism — a possible later lever for
  *wall time* once memory is fixed and the producer is no longer the
  bottleneck; not part of this plan.

## Expected outcome

- **Memory:** peak RSS roughly flat in N (straddle + one compacted
  block + the in-flight consume blocks), targeting `main`'s ballpark
  instead of ~135 MB/sample.
- **Wall:** load no longer re-folds on growth (Stage 2) and DUST leaves
  the critical path (Stage 3); produce stops being the serial floor that
  starves the parallel consumer, narrowing or closing the gap to `main`
  at realistic N.
- **All byte-identical** to the current branch output at every stage.
