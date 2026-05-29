# Apply chunk-based loading to `estimate-contamination`

Status: planned. Within-chrom rewrite — the last subcommand that
still uses the streaming `PerPositionMerger` over per-sample PSP
record iterators. After this lands, both production subcommands
(`var-calling` and `estimate-contamination`) share the same
chunk loader + variant filter + (in the future) parallel
dispatch.

See the
[Phase A impl report](../reports/implementations/cohort_within_chromosome_parallel_phase_a_2026-05-28.md)
and the
[variant-bounded chunk loading plan](cohort_within_chromosome_parallel_phase_b1_variant_bounded_chunks.md).

## Domain intent

The contamination side-pass (Stage `estimate-contamination`) walks
the cohort-wide per-position pileups and runs an online EM that
estimates per-sample contamination rates `c_s` and per-batch
background allele frequencies `q_b`. The math is well-isolated in
[`estimate_contamination`](../../../src/var_calling/contamination_estimation.rs)
and consumes any `Iterator<Item = Result<PerPositionPileups, _>>`
upstream.

Today
[`run_estimate_contamination`](../../../src/pop_var_caller/estimate_contamination.rs#L347)
builds that iterator from `PerPositionMerger::new(record_iters,
…).records()` — the streaming per-position merger over per-sample
PSP record iterators. That's the same upstream shape `var-calling`
used before Phase A; the chunk-based rewrite has been applied to
`var-calling` only.

This change applies the same chunk-based rewrite to
`estimate-contamination` so:
- The cohort-wide variant filter drops monomorphic positions
  before they reach the side-pass — the side-pass's step 1a
  filter would reject them anyway.
- Memory is bounded per-chunk (driven by
  `--target-variants-per-chunk`) instead of growing with PSP block
  size.
- Phase B's parallel-window dispatch + Phase C's pipelined chunk
  loading become reachable for `estimate-contamination` too,
  whenever they ship.

## Feature contract

**Inputs.** Unchanged at the CLI level:
- `--psp-files`, `--reference`, `--output`,
  `--batch-assignment`, `--threads`, the stopping-mode knobs,
  the filter thresholds, etc.
- New: `--target-variants-per-chunk <N>` (default `0` =
  BP-only loop; same semantics as on `var-calling`).

**Outputs.** Unchanged: the `.estcontam` artefact at `--output`
must remain byte-identical to today's output for any given
input + config.

**Edge cases.**
- Empty cohort, single chromosome, chromosome with zero variant
  positions: same behaviour as today (the side-pass either
  surfaces `InsufficientSites` / `DidNotConverge` or completes,
  depending on stopping mode).
- Carryover semantics: the side-pass does not group positions
  across windows or chunks — each position is processed
  independently. So no carryover from one chunk to the next is
  required for correctness.
- Chunk-boundary positions: irrelevant to the side-pass for the
  same reason — there is no group spanning chunk boundaries to
  worry about.

**Performance.** Wins come from two places:
- The cohort-wide variant filter drops every position the
  side-pass would reject at step 1a (monomorphic cohort) AND
  every position the side-pass would do work for but ultimately
  reject (single allele type after filtering). Today every PSP
  record's position reaches the side-pass; with the chunk-based
  variant filter, only positions where the cohort signal exists
  do.
- Memory bound on the cohort cross-section. Today
  `PerPositionMerger` buffers per-sample records up to its
  internal high-water mark; with chunks the working set is
  bounded by `target_variants_per_chunk × n_samples` plus the
  chunk's columnar overhead.

**API compatibility.** Internal types only. CLI shape unchanged
beyond the new knob.

**Concurrency.** The side-pass itself stays sequential (online EM
updates serialise). The chunk loader is sequential per
chromosome. Phase B's parallel windows are out of scope here —
they apply to the var-calling worker, not the side-pass.

## Design

### What gets reused

- [`load_chunk_from_iters`](../../../src/var_calling/cohort_block/loader.rs)
  — unchanged.
- [`MaterialisedChunk`](../../../src/var_calling/cohort_block/columns.rs)
  + per-sample
  [`SampleColumns`](../../../src/var_calling/cohort_block/columns.rs)
  — unchanged.
- The variant filter (carried inside `load_chunk_from_iters`) —
  unchanged.
- The per-chromosome PSP-iterator wiring from
  [`drive_cohort_chunked`](../../../src/var_calling/cohort_block/driver.rs)
  — pattern reused.

### What gets skipped

The side-pass does not need:
- The pre-pass (`fix_boundaries`) — it picks `safe_end` for
  carryover semantics that don't apply here.
- The window partition + worker — no parallel processing inside
  a chunk because the online EM is sequential.
- The ref fetcher — it's only used by the worker's per-group
  scalar projection.

### What gets built

A new module **`src/pop_var_caller/contamination_chunked_stream.rs`**
(or a sibling under `cohort_block` — see "Module placement"
below) exposes:

```rust
/// Stream of `PerPositionPileups` driven by the chunk loader.
/// Implements `Iterator<Item = Result<PerPositionPileups,
/// ContaminationStreamError>>` so it slots straight into
/// `estimate_contamination`.
pub struct ChunkedPositionStream<R: Read + Seek> {
    psp_readers: Vec<PspReader<R>>,
    chromosomes: Vec<ParsedChromosome>,
    chrom_idx: usize,
    chrom_cursor: u32,
    chunk_scratch: ChunkLoadScratch,
    chunk: MaterialisedChunk,
    chunk_position_cursor: usize,
    chunk_positions: Vec<u32>,
    n_samples: usize,
    chunk_genomic_span: u32,
    target_variants_per_chunk: u32,
    max_chunk_span_growth: u32,
}

#[derive(Error, Debug)]
pub enum ContaminationStreamError {
    #[error("chunk load: {0}")]
    ChunkLoad(#[from] ChunkLoadError<PspReadError>),
    #[error("psp read: {0}")]
    PspRead(#[from] PspReadError),
}
```

`next()` semantics:

1. If `chunk_position_cursor < chunk_positions.len()`, build a
   `PerPositionPileups` for `chunk_positions[chunk_position_cursor]`
   by walking each sample's columns at that position (binary
   search per sample). Advance the cursor and return the item.
2. If the current chunk is exhausted, advance the chunk:
   - If `chrom_cursor < chrom_length`: load the next chunk on the
     same chromosome.
   - Else: bump `chrom_idx`; if past the end, return `None`.
3. On a fresh chunk: compute `chunk_positions` as the sorted +
   dedup'd union of `chunk.per_sample[s].positions`. Set
   `chunk_position_cursor = 0` and recurse.

Carryover isn't needed for the side-pass; the loader is invoked
with empty `carryover` slots so each chunk is independent.

### Module placement

Two options:

- **(a) Live under `cohort_block`** as a sibling of `driver` /
  `worker`. Pros: cleanly factored next to the other consumers
  of the chunk machinery. Cons: `PerPositionPileups` is from
  `per_position_merger`, which is a streaming-pipeline type —
  exposing it from `cohort_block` creates a cross-direction
  dependency (column-shape module producing row-shape type).
- **(b) Live under `pop_var_caller`** next to
  `estimate_contamination.rs`. Pros: it's a side-pass-specific
  bridge from the chunk machinery to a streaming type. Cons:
  duplicates the "open PSP readers + iterate chromosomes" pattern
  that `driver.rs` already owns.

Option (b) wins on shape — the bridge is only used by
`estimate-contamination`; a future migration of the
`PerPositionPileups` type out of `per_position_merger` (Phase
A.2-style) would re-evaluate this. Module placed at
`src/pop_var_caller/contamination_chunked_stream.rs`.

### What `run_estimate_contamination` looks like after

```rust
// 7. Build the chunk-based position stream (replaces
//    `PerPositionMerger::new`).
let stream = ChunkedPositionStream::new(
    readers,
    chromosomes,
    cfg.chunk_genomic_span,
    cfg.target_variants_per_chunk,
)?;

// 8. Run side-pass — `estimate_contamination` is unchanged.
let estimates = estimate_contamination(
    stream,
    n_samples,
    sample_to_batch.clone(),
    n_batches,
    cfg.clone(),
)?;
```

The `ContaminationEstimationConfig` gains `chunk_genomic_span:
u32` + `target_variants_per_chunk: u32`. The CLI gets a new
`--target-variants-per-chunk` knob.

`estimate_contamination`'s signature stays:

```rust
pub fn estimate_contamination<I>(
    upstream: I,
    n_samples: usize,
    sample_to_batch: Vec<usize>,
    n_batches: usize,
    config: ContaminationEstimationConfig,
) -> Result<ContaminationEstimates, ContaminationEstimationError>
where
    I: Iterator<Item = Result<PerPositionPileups, /* error type */>>,
```

We only need to relax the error type from
`PerPositionMergerError` to a trait bound — easy: either
introduce an `Into<ContaminationEstimationError>` requirement on
the iterator's error or thread a new error variant through
`ContaminationEstimationError`. Pick the second to keep the
signature concrete.

## Steps

Each step ships as its own commit + tests. The byte-identity
oracle across every step is the current streaming pipeline's
output on the same fixture; the
`estimate_contamination_then_var_calling_chain` integration test
already covers this end-to-end.

### Step 1 — `ChunkedPositionStream` skeleton + tests

Build the new module with the iterator impl. Unit tests over
synthetic per-sample record fixtures:
- Single-chrom, dense variants: all positions emitted in order.
- Sparse cohort: only kept positions emitted.
- Multiple chunks within one chromosome: positions emitted
  across chunk boundaries.
- Multi-chromosome: positions emitted across chromosome
  boundaries.
- Empty cohort: iterator yields `None` immediately.
- Upstream PSP error: surfaced as `ContaminationStreamError`.

Module file:
[src/pop_var_caller/contamination_chunked_stream.rs](../../../src/pop_var_caller/contamination_chunked_stream.rs)
(new).

### Step 2 — Loosen `estimate_contamination`'s error generic

Today the function is generic only over an `Iterator<Item =
Result<PerPositionPileups, PerPositionMergerError>>`. Generalise:

```rust
pub fn estimate_contamination<I>(...) -> Result<..., ContaminationEstimationError>
where
    I: Iterator<Item = Result<PerPositionPileups, ContaminationStreamError>>,
```

— OR add a separate entry point that accepts the new error.
The first option is cleaner; the streaming pipeline's
`PerPositionMergerError` becomes a `From` conversion into
`ContaminationStreamError::Upstream`.

Re-run the existing 67 `contamination_estimation` lib tests — they
all construct streams of `Result<PerPositionPileups,
PerPositionMergerError>` so they should still pass after the
`From` impl is added.

### Step 3 — Wire `run_estimate_contamination` to the new stream

Replace the `PerPositionMerger::new(...)` + `merger.records()`
call sites in
[run_estimate_contamination](../../../src/pop_var_caller/estimate_contamination.rs#L347)
with a `ChunkedPositionStream::new(...)` call. Add
`chunk_genomic_span` + `target_variants_per_chunk` to
`ContaminationEstimationConfig`; thread them through the CLI.

End-of-run summary line gains `chunks_loaded=N
avg_variants_per_chunk=N.M` matching the var-calling pattern.

### Step 4 — CLI knob

Add `--target-variants-per-chunk <N>` to the
`estimate-contamination` subcommand. Default `0` (legacy BP-only
behaviour). Help text mirrors the var-calling flag's wording.

### Step 5 — Byte-identity gates

Two test gates:
- Re-run the existing `estimate_contamination_then_var_calling_chain`
  cohort CLI integration test → must pass byte-identically.
- Add a new
  `estimate_contamination_byte_identical_across_target_variants_per_chunk`
  test mirroring step 6 of the variant-bounded chunks plan: same
  fixture run twice with `target = 0` vs `target = 1`, asserting
  the produced `.estcontam` artefact contents match.

## Risks

- **Position ordering across chromosomes.** The streaming
  `PerPositionMerger` emits positions in `(chrom_idx, pos)` order.
  `ChunkedPositionStream` must do the same — chunks loaded in
  the order chromosomes appear in `chromosomes`, positions within
  a chunk sorted ascending. Easy to get wrong if the per-sample
  union is built lazily. **Mitigation:** integration test asserts
  byte-identical artefact contents → any ordering drift surfaces.
- **`PerPositionPileups` shape.** Today
  `PerPositionMerger` builds `per_sample: Vec<Option<PileupRecord>>`
  with `None` for samples that don't have a record at the
  position. The chunk-based stream must mirror this. **Mitigation:**
  unit tests assert `per_sample[s].is_none()` for samples without
  records at the position.
- **Carryover semantics.** Today the chunk loader is invoked with
  carryover from the previous chunk; the side-pass doesn't need
  this. Calling with empty carryover slots on every chunk is the
  right thing. The pre-pass / partition / worker are skipped —
  no `safe_end` is set, so `chunk.safe_end == range.end` and the
  next chunk starts at the current chunk's `range.end`. **Risk:**
  records past the loader's stop point on the current chunk would
  be re-pulled on the next chunk. **Mitigation:** the loader pulls
  records `< attempt_end` and pauses the iterator there; the next
  chunk's loader call passes a higher `range_start`, the iterator
  resumes from where it left off because PSP iterators are
  position-indexed.

  Wait — looking again, `region_records(chrom_id, psp_cursor,
  psp_inclusive_end)` constructs a *fresh* iterator each call. So
  per-call re-seeking via the PSP block index per chunk is
  unavoidable; not a perf concern in practice (block index lookup
  is microseconds), but worth flagging.
- **Toolchain clippy debt.** The 17 pre-existing errors are still
  there; this change should not add to them.

## Out of scope

- Phase B parallel windows for the side-pass — the side-pass's
  online EM is inherently sequential. Per-chunk parallelism is
  not on the table here.
- Phase C pipelined chunk loading — applies to both subcommands
  once it lands; not a contamination-specific concern.
- `PerPositionPileups` columnar refactor — the side-pass still
  consumes the row-shape type. Migrating that to column-shape is
  separate.

## Validation

- `cargo fmt --check`
- `cargo clippy --lib --tests -- -D warnings`
- `cargo test --lib`
- `cargo test --tests`
- `cargo test --test cohort_cli_integration`
- 3-tomato contamination estimation byte-identity check with
  `target = 0` and a non-zero value.

## Estimated effort

~2 sessions across 5 commits, comparable to a single phase of the
variant-bounded chunks plan.
