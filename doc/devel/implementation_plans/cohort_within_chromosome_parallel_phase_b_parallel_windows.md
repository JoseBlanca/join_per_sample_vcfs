# Phase B — parallel within-chunk window dispatch

Status: planned. Final substantive piece of the within-chromosome
rewrite. Variant-bounded chunks
([Phase B prereq](cohort_within_chromosome_parallel_phase_b1_variant_bounded_chunks.md))
gives stable per-chunk workloads; the probe-and-slide algorithm in
[`fix_boundaries`](../../../src/var_calling/cohort_block/pre_pass.rs)
(commit `3687a1f`) already supports `target_window_count > 1`. This
step wires those windows to a rayon dispatch.

## Domain intent

The driver currently runs the windows in `chunk.windows` sequentially
inside the chunk-loop body, one `run_window` call per window. With
the algorithm in place to produce `T` non-overlapping windows per
chunk, the per-window math (layers 1+2+3 + EM) is embarrassingly
parallel — each window operates on a disjoint slice of the chunk's
columnar data and produces its own batch of `PosteriorRecord`s.

The CPU budget for the cohort pipeline is dominated by the EM + the
per-group merger work
([2026-05-27 perf review](../reports/reviews/scaling_measurement_2026-05-27.md):
`per_group_merger` ~68 % at N=1000, `posterior_engine` ~11 % at
N=1000). The within-chunk parallelism is where the headline scaling
win lives.

## Feature contract

**Inputs.** Same chunk loop as today, plus:
- `ChunkDriverParams.target_window_count: usize` — desired number
  of parallel windows per chunk. `1` (default) preserves today's
  sequential behaviour byte-for-byte.
- CLI knob `--worker-windows-per-chunk <N>` on `var-calling`.

**Outputs.** Same `PosteriorRecord` sequence as today, in the same
order. Byte-identity is the gate.

**Edge cases.**
- Sparse chunk: `fix_boundaries`' probe-and-slide may return fewer
  than `T` windows (some slides collide). The driver handles the
  actual `chunk.windows.len()` regardless of the configured `T`.
- Empty windows: a window may contain zero variant groups; the
  worker fast-paths to `Ok(())` without ref-fetcher access.
- `T == 1`: the parallel branch reduces to a single sequential call
  (no rayon overhead worth speaking of).
- Cap firing in the loader: orthogonal to window dispatch.

**Performance.** Expected speedup ≈ `min(T, available_cores)` on the
worker compute portion. End-to-end wall improvement depends on what
fraction of the chunk loop is compute vs I/O (load + fmt + emit).

**API compatibility.** Internal types only. CLI surface gains one
knob.

**Concurrency.**
- Workers run on rayon's thread pool. The existing
  `configure_rayon_pool(args.threads)` setup already builds the
  pool once per process.
- Workers within a chunk run concurrently; chunks themselves stay
  sequential (the chunk loader is single-threaded; Phase C
  pipelines chunk loading separately).
- Output ordering across workers is by window index (rayon's
  `par_iter` + `collect` preserves input order — driver collects
  per-window output buffers and drains them in order).

## Design

### Per-worker state pool

Today the driver owns one `ColumnarPipelineScratch`, one
`PartitionScratch`, and one `WindowPartition` — passed by `&mut`
to the sequential `run_window` call. For parallel dispatch, the
driver owns a **pool** sized to `target_window_count`:

```rust
pub struct WorkerPool {
    target: usize,
    scratches: Vec<ColumnarPipelineScratch>,
    partition_scratches: Vec<PartitionScratch>,
    partitions: Vec<WindowPartition>,
    output_bufs: Vec<Vec<PosteriorRecord>>,
    // Per-chrom ref fetchers — see the "Ref fetcher" section below.
    fetchers: Vec<StreamingChromRefFetcher>,
}
```

`ensure_capacity(n: usize, n_samples: usize)` grows the pool when
`n > self.target` and updates `target`. Buffers in entries
`[0, target)` keep their high-water-mark capacity across chunks.

### Ref bytes — pre-fetch serially, share `&[u8]` across workers

`StreamingChromRefFetcher` is `Send + !Sync` by design — it uses a
sliding buffer over a single open file. Sharing one fetcher across
parallel worker threads would race on the buffer state.

Rather than open `T` fetchers per chromosome (which works but
costs `T` file descriptors and `T` FAI loads per chrom), the
driver pre-fetches every group's REF bytes **sequentially** with
the single per-chrom fetcher *before* the parallel dispatch.
Workers read from the pre-fetched bytes via a `&[u8]` borrow —
read-only state is `Sync` for free.

Concretely, the `WorkerPool` carries one
`pre_fetched_ref_bytes: Vec<Vec<Vec<u8>>>` indexed by
`[window_idx][group_idx]`. Populating it is a sequential walk:

```text
for (window_idx, partition) in partitions[..n_windows].iter().enumerate():
    let buf = &mut pre_fetched_ref_bytes[window_idx]
    buf.clear()
    for group_idx in 0..partition.n_groups():
        let span = partition.group_ends[group_idx] - partition.group_starts[group_idx] + 1
        let mut bytes = mem::take(&mut buf_reuse[...])
        fetcher.fetch_into(partition.group_starts[group_idx], span, &mut bytes)
        buf.push(bytes)
```

The single chrom fetcher stays monotonic because the prefetch
walks windows in index order (= BP order) and groups within each
window in BP order. No fetcher pool, no monotonic-access
invariant across workers, no `T` file descriptors.

The trade-off vs the fetcher-pool design:
- Per-group `Vec<u8>` allocation up front (small — group spans
  are typically < 1 KB).
- Pre-fetch is single-threaded; parallel workers consume
  pre-fetched bytes only. Pre-fetch cost is `O(total REF span)`
  per chunk, which the existing sequential code already pays.
- Workers consume `&[u8]` (Sync, no contention).

### Driver flow

`load_and_run_chunk_with_retry` becomes:

```text
load_chunk_from_iters(...)
pre_pass picks safe_end + emits T windows
pool.ensure_capacity(chunk.windows.len(), n_samples)

# Step 1 — partition all T windows sequentially.
for (window, p_scratch, partition) in zip(chunk.windows, pool.partition_scratches, pool.partitions):
    partition_window(chunk, window, ..., p_scratch, partition)

# Step 2 — run T workers in parallel.
windows.par_iter()
    .zip(partitions.par_iter())            # &WindowPartition
    .zip(scratches.par_iter_mut())         # &mut ColumnarPipelineScratch
    .zip(fetchers.par_iter_mut())          # &mut StreamingChromRefFetcher
    .zip(output_bufs.par_iter_mut())       # &mut Vec<PosteriorRecord>
    .try_for_each(|((((window, partition), scratch), fetcher), output_buf)|
        run_window(chunk, partition, fetcher.clone(), per_group, posterior, scratch, output_buf)
    )?

# Step 3 — drain output_bufs in window order.
for output_buf in &mut pool.output_bufs[..chunk.windows.len()]:
    for record in output_buf.drain(..):
        emit_or_drop(record, ...)
```

The clone in `run_window(... fetcher.clone() ...)` doesn't apply
here — the fetcher is by `&mut` not `Arc`. The signature of
`run_window` will change to take `&dyn ChromRefFetcher` or a
generic `&mut F: ChromRefFetcher` instead of `SharedRefFetcher`.

### `run_window` signature change

Currently:

```rust
pub fn run_window(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    ref_fetcher: SharedRefFetcher,
    per_group_config: PerGroupMergerConfig,
    posterior_config: PosteriorEngineConfig,
    scratch: &mut ColumnarPipelineScratch,
    output: &mut Vec<PosteriorRecord>,
) -> Result<(), PosteriorEngineError>
```

After Phase B:

```rust
pub fn run_window(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    pre_fetched_ref_bytes: &[Vec<u8>],  // one entry per group, in group_idx order
    per_group_config: PerGroupMergerConfig,
    posterior_config: PosteriorEngineConfig,
    scratch: &mut ColumnarPipelineScratch,
    output: &mut Vec<PosteriorRecord>,
) -> Result<(), PosteriorEngineError>
```

The fetcher is gone from `run_window`. Internally,
`build_posterior_record_columnar` reads
`&pre_fetched_ref_bytes[group_idx]` instead of calling
`ref_fetcher.fetch_into`. The `scratch.ref_buf` field stays for
local copies if any layer needs an owned `Vec<u8>` (today
`unify_alleles_columnar` takes `&[u8]`, so no copy is needed).

## Steps

Each step ships as its own commit + tests. The byte-identity oracle
across every step is the current sequential baseline
(`target_window_count = 1`).

### Step 1 — `run_window` takes pre-fetched REF bytes

Switch `run_window`'s signature from `SharedRefFetcher` to
`pre_fetched_ref_bytes: &[Vec<u8>]`. The worker's
`build_posterior_record_columnar` reads from
`&pre_fetched[group_idx]` instead of calling `fetch_into`.

Driver-side: before each `run_window` call, the driver does a
small sequential prefetch loop using the per-chrom fetcher to
populate a `Vec<Vec<u8>>` for the window's groups. Behaviour
unchanged — same fetcher, same fetch order. 3-tomato byte-identity
check must pass.

### Step 2 — `WorkerPool` struct + driver refactor

Introduce `WorkerPool` holding the per-window scratches +
partition_scratches + partitions + output_bufs +
pre_fetched_ref_bytes. Driver-side: `drive_cohort_chunked`
allocates one pool, threads it through
`drive_one_chrom_generic` → `load_and_run_chunk_with_retry`.

For `target_window_count == 1`, the pool has exactly one entry
and the chunk loop runs the same sequential code path as today.
Byte-identity check.

### Step 3 — `ChunkDriverParams.target_window_count` + driver loop honours it

Add the field (default `1`). Pass to `fix_boundaries` (replacing
the hard-coded `1`). The pool resizes per chunk based on the
actual `chunk.windows.len()` from `fix_boundaries` (which may
be less than `target_window_count` when slides collide).

Sequential path still works; just iterates the pool's first
entries in window order. Byte-identity check.

### Step 4 — Parallel dispatch

The for-loop over `chunk.windows` becomes a `par_iter().zip(...)`
chain that calls `run_window` per window concurrently. The
prefetch step stays sequential (one fetcher); only the math
parallelises. Output buffers are drained sequentially after the
parallel scope in window order.

Byte-identity check at `target_window_count = 4`:

```text
records_emitted = expected_for_main
chunks_loaded matches T=1 run
VCF body diff vs T=1 run = 0 lines
```

### Step 5 — CLI knob

Add `--worker-windows-per-chunk <N>` to `var-calling`. Default
`1` (preserves legacy behaviour). Plumbed through
`ChunkDriverParams.target_window_count`. Help text explains the
trade-off.

### Step 6 — Byte-identity gate

New integration test mirroring the variant-count knob's
byte-identity test: same fixture run twice with
`target_window_count = 1` vs `target_window_count = 4`, asserting
the VCF bodies match.

## Risks

- **Pre-fetched bytes memory.** `Vec<Vec<u8>>` per window holds
  the REF bytes for every group. Per-chunk peak ≈ Σ group_span ≈
  the post-filter chunk variant footprint × cohort-typical group
  span. For a 10 K-variant chunk with median 10 bp groups that's
  ~100 KB per chunk total — negligible.
- **EM determinism under parallelism.** The EM is deterministic
  given its inputs; per-window inputs are deterministic given the
  chunk + partition + pre-fetched bytes; the partition is
  deterministic given chunk + window. So byte-identity across
  parallel runs is structurally guaranteed.
- **Rayon pool sharing across `var-calling` and other subcommands.**
  The `configure_rayon_pool(args.threads)` call locks the pool size
  globally. If a future subcommand wants different parallelism,
  it'll need separate pools — out of scope here.
- **Output buffer growth.** Per-window `Vec<PosteriorRecord>`
  grows to ~groups-per-window records and is drained on every
  chunk. Capacity reuse via `Vec::drain` keeps allocations bounded.
- **Pre-fetch is sequential.** Only the worker math parallelises.
  If the pre-fetch becomes the bottleneck (unlikely — it's a
  contiguous sliding-buffer read), a future optimisation could
  switch to a `Sync` fetcher (mmap-based) and parallelise the
  pre-fetch too.

## Validation

- `cargo fmt --check`
- `cargo clippy --lib --tests -- -D warnings`
- `cargo test --lib`
- `cargo test --tests`
- `cargo test --test cohort_cli_integration`
- 3-tomato cohort byte-identity at `target_window_count = 4`.

## Out of scope

- Phase C pipelined chunk loading. Composes cleanly on top of this
  but is its own change.
- SIMD perf review (Phase D).
- Adaptive `target_window_count` (e.g. derived from
  `target_variants_per_chunk / target_variants_per_window`). Once
  Phase D establishes per-window cache budgets, the derivation
  becomes obvious.
- Fetcher-sharing optimisations (mmap-based Sync fetcher, etc.).
  Phase B's per-window fetcher pool is correct and bounded; perf
  review motivates further work.

## Estimated effort

~2 sessions across 6 commits, each byte-identity-gated.
Comparable to the variant-bounded-chunks work — the worker
signature change plus the rayon zip-chain are the substantive
pieces.
