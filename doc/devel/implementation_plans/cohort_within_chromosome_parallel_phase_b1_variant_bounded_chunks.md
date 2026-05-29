# Variant-count-bounded chunk loading

Status: planned. Phase B prerequisite. Decouples worker batch size
from PSP block size + density variation; gives Phase B's parallel
window dispatch a stable workload to split.

See the
[Phase A impl report](../reports/implementations/cohort_within_chromosome_parallel_phase_a_2026-05-28.md)
and the
[Phase B algorithm change](../../../src/var_calling/cohort_block/pre_pass.rs)
that landed in commit `3687a1f`.

## Domain intent

A chunk is the loader → worker handoff unit. Today the chunk's size
is BP-bounded by [`ChunkDriverParams::chunk_genomic_span`](../../../src/var_calling/cohort_block/driver.rs)
(default ≈ a few MB). Variant density varies 100× across the genome
(telomeres vs centromeres vs euchromatin), so a fixed BP span yields
unstable variant counts per chunk — sparse regions hit O(10)
variants, dense ones O(10 000). Two consequences:

- **Phase B parallelism is wasted on small chunks.** When a chunk
  has fewer variants than the configured `target_window_count`,
  windows are degenerate; some workers do nothing.
- **Per-chunk overhead** (allocator churn, partition build,
  ref-fetch setup) doesn't scale with chunk variant count, so the
  proportional overhead spikes on sparse chunks.

This feature switches the chunk-loader loop from "advance by a fixed
BP step" to "advance until accumulated variant count ≥
`target_variants_per_chunk`, capped at a BP safety ceiling".

The user-facing knob is variant count, in line with the perf
review's recommendation: workers are CPU-bound on the EM, so what
matters is the work amount they receive, not the BP range it spans.

## Feature contract

**Inputs**

- `ChunkDriverParams.target_variants_per_chunk: u32` — primary
  knob. Soft lower bound on the variant count per loaded chunk.
- `ChunkDriverParams.chunk_genomic_span: u32` — kept as the initial
  BP nominal span; behaves as a hard ceiling per retry attempt.
- `ChunkDriverParams.max_chunk_span_growth: u32` — already exists
  (`MAX_CHUNK_SPAN_GROWTH = 8`); reused for the under-variant retry.

**Outputs**

- Same `MaterialisedChunk` shape. `chunk.range.end` reflects the
  actual BP end the loader stopped at (= the last attempt's
  `chunk_range_end`).
- A new `ChunkDriverStats.under_variant_retries: u64` counter
  reports how often the driver had to grow a chunk's span to hit
  the variant target.

**Edge cases**

- Chromosome shorter than the first attempt's reach: load whatever
  is there; emit a final chunk regardless of target.
- Dense region: the first attempt already exceeds target; no
  retry needed.
- Sparse region: the first attempt is far under target; retry by
  doubling span up to `max_span = chunk_genomic_span ×
  MAX_CHUNK_SPAN_GROWTH`. If still under target at the cap, emit
  whatever is there (graceful degradation; perf hit but no error).
- Carryover already contains target variants: the loader's
  carryover-drain step accounts for them; no PSP reads beyond the
  carryover may be needed.
- `target_variants_per_chunk = 0`: treated as "no minimum"; behaves
  like today's BP-only loop.

**Performance**

- One extra `usize` count per load attempt. The variant count is
  already a byproduct of the loader's variant filter (`is_kept`).
- Retry reuses the existing carryover-snapshot dance, so wasted
  decompression on under-variant retries is the same shape as the
  existing `NoSafeGap` retries.
- Memory: chunk RSS grows linearly with the actual variant count
  loaded. Bounded by `max_span × density`; in pathological
  conditions the BP ceiling caps it.

**API compatibility**

- Internal `cohort_block` types. Public-facing change is the new
  CLI knob.
- `ChunkDriverParams` gains a field; existing constructors expand.
  No on-disk format change.

**Concurrency**

- Single-threaded loader (unchanged). The variant count is observed
  at the end of each load attempt before deciding whether to retry.

## Design

### Why driver-level retry vs loader-level streaming

Two paths were considered:

- **(A) Driver retry with extended range.** Reuse the existing
  retry-with-extended-range loop in
  [`load_and_run_chunk_with_retry`](../../../src/var_calling/cohort_block/driver.rs)
  (currently triggered by `NoSafeGap`). Add an "under variant
  target" branch that grows `attempt_span` the same way.
- **(B) Loader-level streaming.** Restructure
  [`load_chunk_from_iters`](../../../src/var_calling/cohort_block/loader.rs)
  to pull PSP records incrementally, count variants after each
  batch, stop when target reached. Avoids re-decompression but
  requires interleaving per-sample iterator pulls and re-running
  the position-union / variant-filter passes incrementally.

The plan picks (A). Reasons:

- The existing carryover-snapshot retry already handles the
  re-load case; (A) is `cargo diff`-small.
- The loader's variant-filter pass is straightforward to amortise
  on the BP span we already grow during retry — re-decompression
  on the second attempt re-reads PSP blocks that were just decoded,
  which is a known cost shape we accept for `NoSafeGap` already.
- Loader-level streaming (B) is a substantive rewrite of the
  loader's batched position-union logic; deferred until a perf
  review shows the wasted decompression actually matters.

### Where the variant count comes from

The cohort-wide variant count for one chunk is the number of `true`
entries in
[`ChunkLoadScratch.is_kept`](../../../src/var_calling/cohort_block/loader.rs).
Today the loader produces `is_kept` as the post-filter "keep this
position" predicate; counting `true`s is one extra `.iter().filter(|x|
**x).count()` after the existing pass. The count is then exposed via
the loader's return value or a `ChunkLoadStats` field.

Proposal: `load_chunk_from_iters` returns
`Result<ChunkLoadStats, ChunkLoadError<E>>` instead of `Result<(),
…>`. `ChunkLoadStats` carries `variant_count: u32` (the kept-position
count) and any other counters useful for diagnostics. Drivers that
ignore stats just discard the return value.

### Where the retry sits

The retry loop lives in
[`load_and_run_chunk_with_retry`](../../../src/var_calling/cohort_block/driver.rs)
and currently has one termination condition (`NoSafeGap` resolved or
`attempt_span ≥ max_span`). The new termination condition adds:

```text
attempt_span = nominal_span
loop:
    load with chunk_range_end = chunk_range_start + attempt_span
    stats = load result
    if stats.variant_count < target_variants_per_chunk
       AND attempt_span < max_span:
        carryover_snapshot.restore()
        attempt_span *= 2
        under_variant_retries += 1
        continue
    match pre_pass:
        Ok: break
        NoSafeGap + attempt_span < max_span:
            carryover_snapshot.restore()
            attempt_span *= 2
            no_safe_gap_retries += 1
            continue
        Err(_): return
```

Both retry paths use the same `attempt_span *= 2` growth and the
same `carryover_snapshot` restore. Two diagnostic counters
(`under_variant_retries`, `no_safe_gap_retries`) report the cause
mix at end-of-run for capacity planning.

### What the default looks like

A target of 0 means "no minimum"; the loop behaves like today (pure
BP-bounded). Tests that build a driver with default params should
fall through to today's behaviour byte-identically.

CLI default: under discussion. A reasonable starting point is
something like 10 000 — large enough that Phase B's T=4 / 8 / 16
dispatch yields meaningful per-window work even on sparse
chromosomes; small enough that the chunk fits comfortably in
allocator-friendly RSS even at N=1000 samples.

## Steps

Each step ships as its own commit + tests. The byte-identity oracle
across every step is the current Phase B baseline (`target_variants
= 0`).

### Step 1 — loader reports variant count

[`load_chunk_from_iters`](../../../src/var_calling/cohort_block/loader.rs)
returns
`Result<ChunkLoadStats, ChunkLoadError<E>>` where `ChunkLoadStats`
carries:
- `variant_count: u32` (number of `true` in `is_kept`).

The change is mechanical: count after the existing variant-filter
pass; thread the count through the return. All existing callers
discard the value with `.map(|_| ())` until step 3 wires the
counter through.

Unit tests in [`loader.rs`](../../../src/var_calling/cohort_block/loader.rs)
assert the count matches the expected variant set on fixtures
already in the suite.

### Step 2 — driver params + stats fields

[`ChunkDriverParams`](../../../src/var_calling/cohort_block/driver.rs)
gains `target_variants_per_chunk: u32` (default `0`).
[`ChunkDriverStats`](../../../src/var_calling/cohort_block/driver.rs)
gains `under_variant_retries: u64` + `under_variant_retry_failures:
u64`. No behaviour change yet.

### Step 3 — driver loop honours the target

[`load_and_run_chunk_with_retry`](../../../src/var_calling/cohort_block/driver.rs)
adds the under-variant retry branch sketched above. The retry uses
the existing `carryover_snapshot` restore + `attempt_span *= 2`
growth + `max_span` cap. When `target_variants_per_chunk == 0`, the
retry branch is a no-op and the loop matches today's behaviour
byte-for-byte.

New tests use `cohort_cli_integration`-style fixtures with
artificially low `target_variants_per_chunk` (e.g., target = 100 on
a fixture with sparse + dense regions) and check:
1. The chunk loop grows `attempt_span` in the sparse region.
2. The chunk's variant count after growth is `≥ target` or hits
   the `max_span` cap (whichever comes first).
3. Stats reflect the retry count.

### Step 4 — CLI knob

Add `--target-variants-per-chunk <N>` to the `var-calling`
subcommand under [`src/pop_var_caller/var_calling.rs`](../../../src/pop_var_caller/var_calling.rs).
Default: TBD (likely `10000`); revisit after step-3 perf data.
Help text explains the trade-off (large = bigger chunks, fewer
overheads, more memory; small = smaller chunks, more overhead per
chunk).

Plumb through `ChunkDriverParams.target_variants_per_chunk`.

End-of-run summary line adds `under_variant_retries=N` next to the
existing diagnostic counters.

### Step 5 — byte-identity gate + integration test

The 3-tomato cohort run with `--target-variants-per-chunk 0`
(default off) must remain byte-identical to the prior commit.

A new integration test in
[`tests/cohort_cli_integration.rs`](../../../tests/cohort_cli_integration.rs)
runs the same fixture twice — once with `target_variants_per_chunk =
0`, once with a non-zero value — and asserts equal VCF bodies. This
codifies that the variant-count knob is a *workload-sizing* knob,
not a *correctness* knob: changing it must never change emitted
records.

## Tests

Per the rust-feature-implementation skill's testing policy:

- **Happy path.** Sparse region triggers retry; dense region
  doesn't. Chunk variant counts match expectations.
- **Edge cases.**
  - `target = 0` is a no-op (byte-identical to today).
  - `target = u32::MAX` triggers continuous retries until
    `max_span` cap; loads whatever is there and proceeds.
  - Chromosome shorter than `attempt_span`: clamps to chromosome
    end; no infinite retry.
  - Empty chromosome (no variants at all): loop terminates without
    retry; no extra work.
- **Regression.** Phase A's `NoSafeGap` retry must still fire when
  the gap is the limiting constraint, even if variant count is
  satisfied. The two retry conditions must compose correctly when
  both could fire on the same attempt.
- **Byte-identity.** Cohort CLI integration test asserts that
  emitted VCFs do not depend on `target_variants_per_chunk`.

## Risks

- **Retry overshoot.** Doubling `attempt_span` may overshoot the
  target by up to 2× the previous chunk's count. Mitigation: the
  next chunk's load starts from the previous chunk's `safe_end`
  (= a position past the surplus) and the surplus rides in
  carryover, so no variants are lost. The overshoot is a
  performance hit (more PSP reads than needed) but not a
  correctness issue.
- **`max_span` cap insufficient.** In ultra-sparse regions (some
  pathological centromeres / repetitive arrays), `max_span =
  chunk_genomic_span × 8` may still hit zero variants. Mitigation:
  graceful degradation — emit the chunk anyway, bump the
  `under_variant_retry_failures` counter, log the cap-hit in the
  run summary. Operator response: raise `--chunk-genomic-span` or
  lower `--target-variants-per-chunk` at the CLI.
- **Default knob value picks behaviour for users.** Wrong default
  hurts UX. Mitigation: ship default as `0` (no-op) initially;
  pick a non-zero default only after step-3 perf data.
- **Phase B interaction.** The new variant-count knob composes with
  Phase B's `target_window_count` (the second knob the perf
  discussion identified). The default ratios need calibration after
  Phase B's parallel dispatch lands. For now, Phase B remains
  T=1; the two knobs are orthogonal.

## Out of scope

- **Loader-level streaming** (option B above). Deferred until
  perf review shows the wasted decompression on under-variant
  retries matters at scale.
- **Per-chunk variant-count target derivation from cache size**
  (some kind of L2-aware heuristic). The knob is exposed at the
  CLI; auto-tuning is a future cleanup.
- **Variant-count-bounded window placement.** Phase B currently
  takes a fixed `target_window_count`; deriving T from
  `target_variants_per_window ÷ chunk_variant_count` is the perf
  refinement teased in the
  [Phase A.2 plan](cohort_within_chromosome_parallel_phase_a2_em.md)
  — orthogonal to this work.

## Validation

- `cargo fmt --check`
- `cargo clippy --lib --tests -- -D warnings`
- `cargo test --lib`
- `cargo test --tests`
- `cargo test --test cohort_cli_integration`
- 3-tomato cohort byte-identity check with both `target = 0` and
  the new default value.

## Estimated effort

~2 sessions across 5 commits (one per step), comparable in size to
a Phase A.1 layer.
