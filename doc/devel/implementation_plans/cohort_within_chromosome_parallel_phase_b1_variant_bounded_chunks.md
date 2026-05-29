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

### Loader-level streaming (option B)

The loader owns the "how much to load" decision. Two paths were
considered:

- **(A) Driver-level retry with extended range.** Reuse the
  existing retry-with-extended-range loop in
  [`load_and_run_chunk_with_retry`](../../../src/var_calling/cohort_block/driver.rs)
  (currently triggered by `NoSafeGap`) for the under-variant case
  too. Cheapest diff but the loader's contract becomes "load
  whatever fits, you'll call me again if it wasn't enough", and
  every retry re-reads + re-decompresses the same PSP blocks.
- **(B) Loader-level streaming.** Restructure
  [`load_chunk_from_iters`](../../../src/var_calling/cohort_block/loader.rs)
  to pull PSP records incrementally — wrap each per-sample iterator
  in `Peekable`, pull all records with `pos < attempt_end`, run
  the position-union / variant-filter pass on the accumulated
  data, count variants, grow `attempt_end` and resume the pull
  if we're under target. PSP blocks are decoded once.

The plan picks (B). Cleaner end state:

- The loader's contract becomes "give me a chunk with ≥ N
  variants" — semantically aligned with the perf goal.
- The driver's chunk loop becomes one call per chunk; the
  `attempt_span × 2` + `carryover_snapshot` dance shrinks to the
  rare `NoSafeGap` case only.
- No re-decompression on the common (under-variant) retry path.
  The iterators stay open across the loader's internal extension
  steps and resume from where they paused.

The trade-off vs (A) is a bigger loader rewrite. Mitigated by
splitting the work into byte-identity-gated steps below.

### Streaming pull mechanics

Each per-sample iterator is wrapped in `std::iter::Peekable`. The
inner pull is:

```text
for each per-sample iter:
    while iter.peek().is_ok_and(|r| r.pos < attempt_end):
        record = iter.next()?
        push_into raw_per_sample[s]
```

The peek-and-pull stops cleanly at the first record past
`attempt_end` without consuming it — the next extension iteration
picks up at the same point.

The position-union + `is_kept` filter is recomputed on the full
accumulated data after each extension (CPU work proportional to
`n_positions_so_far`, no PSP I/O). Incremental updates are a perf
refinement deferred until profiling justifies them.

### Where the variant count comes from

The cohort-wide variant count for one chunk is the number of `true`
entries in
[`ChunkLoadScratch.is_kept`](../../../src/var_calling/cohort_block/loader.rs).
Today the loader produces `is_kept` as the post-filter "keep this
position" predicate; counting `true`s is one extra
`.iter().filter(|x| **x).count()` after the existing pass.

The loader returns `Result<ChunkLoadStats, ChunkLoadError<E>>`
where `ChunkLoadStats` carries:

- `variant_count: u32` (the kept-position count),
- `attempt_end: u32` (the actual BP end the loader stopped at —
  the chunk's `range.end`),
- `under_variant_extensions: u32` (how many internal extension
  iterations the loader ran on this chunk).

Drivers that ignore stats discard the return value with `.map(|_|
())`.

### Loader contract (new)

```rust
pub fn load_chunk_from_iters<I, E>(
    scratch: &mut ChunkLoadScratch,
    out: &mut MaterialisedChunk,
    chrom_id: u32,
    range_start: u32,
    initial_span: u32,
    target_variants: u32,    // 0 = no minimum (legacy behaviour)
    max_span: u32,
    per_sample_iters: Vec<I>,
    carryover: &mut [SampleColumns],
) -> Result<ChunkLoadStats, ChunkLoadError<E>>
```

`range_start` + `initial_span` set the first attempt's end;
`max_span` caps growth; `target_variants` is the soft lower bound
on the kept-position count. Caller constructs iterators with end =
chromosome length (or a generous bound) so the loader can pull
beyond `initial_span` when needed.

### Driver simplification

[`load_and_run_chunk_with_retry`](../../../src/var_calling/cohort_block/driver.rs)
keeps the `NoSafeGap` retry (with `carryover_snapshot` restore + a
fresh load), but the inner extension loop is gone. Pseudocode:

```text
loop:
    construct fresh per-sample iterators
    stats = load_chunk_from_iters(target_variants, max_span)
    stats.variant_count and stats.attempt_end now reflect the
        loader's final accumulated state.
    match pre_pass(chunk, target_window_count):
        Ok(()): break
        NoSafeGap and max_span_can_still_grow:
            carryover_snapshot.restore()
            max_span *= 2
            no_safe_gap_retries += 1
            continue
        Err(other): return
```

### What the default looks like

A `target_variants = 0` is a no-op — the loader's extension loop
exits immediately on the first attempt and returns whatever fits
in `initial_span`. Tests that build a driver with default params
fall through to today's behaviour byte-identically.

CLI default: under discussion. A reasonable starting point is
~10 000 — large enough that Phase B's T=4 / 8 / 16 dispatch yields
meaningful per-window work even on sparse chromosomes; small
enough that the chunk fits comfortably in allocator-friendly RSS
even at N=1000 samples.

## Steps

Each step ships as its own commit + tests. The byte-identity oracle
across every step is the current Phase B baseline (`target_variants
= 0`).

### Step 1 — `ChunkLoadStats` return + `variant_count`

[`load_chunk_from_iters`](../../../src/var_calling/cohort_block/loader.rs)
returns `Result<ChunkLoadStats, ChunkLoadError<E>>` instead of
`Result<(), …>`. `ChunkLoadStats` carries `variant_count: u32`
(number of `true` in `is_kept` after the existing filter pass).

Mechanical change: one count after the existing variant-filter pass;
thread the stats through the return. All existing callers discard
the value with `.map(|_| ())` until step 3 wires it through.

Unit tests in
[`loader.rs`](../../../src/var_calling/cohort_block/loader.rs)
assert the count matches the expected variant set on fixtures
already in the suite.

### Step 2 — `Peekable` wrapper + restructured pull (no behaviour change)

Wrap each per-sample iterator in `std::iter::Peekable`. Restructure
the carryover-drain + record-pull blocks into the streaming shape
sketched above — pull only records with `pos < attempt_end` where
`attempt_end` is the input `range.end` for now.

No variant-count check or extension yet. The behaviour is byte-
identical to step 1; this step is a pure refactor that sets up the
extension loop in step 3.

Existing loader tests + the cohort CLI integration test catch any
regression.

### Step 3 — under-variant extension loop in the loader

Add `target_variants: u32` + `max_span: u32` parameters to
`load_chunk_from_iters`. Wrap the existing pull + filter pass in
an `extension` loop:

```text
attempt_end = range_start + initial_span
loop:
    pull records with pos < attempt_end (resumable via Peekable)
    rebuild position_union + is_kept on accumulated data
    variant_count = count true in is_kept
    if variant_count >= target_variants or attempt_end >= max_end:
        break
    attempt_end = min(2 * attempt_end - range_start + range_start, max_end)
    stats.under_variant_extensions += 1
```

`target_variants = 0` makes the loop exit after the first iteration
(the `>= target_variants` check is always true with `0`), preserving
today's behaviour exactly.

Unit tests:
- Sparse fixture: target = 100, initial_span = 50bp. Multiple
  extensions; final variant count `≥ 100` or `max_span` reached.
- Dense fixture: target = 100 already met by initial_span. Zero
  extensions; stats.under_variant_extensions == 0.
- Edge: target = 0 → no extensions ever.
- Edge: target very large + max_span hit → graceful degradation,
  variant_count is whatever was loaded.

### Step 4 — driver wiring

[`ChunkDriverParams`](../../../src/var_calling/cohort_block/driver.rs)
gains `target_variants_per_chunk: u32` (default `0`).
[`ChunkDriverStats`](../../../src/var_calling/cohort_block/driver.rs)
gains a `chunk_extensions: u64` counter aggregated across all loads.

[`load_and_run_chunk_with_retry`](../../../src/var_calling/cohort_block/driver.rs)
constructs each per-sample iterator with end = chromosome length
(or a generous bound), passes `target_variants` + `max_span` to
the loader, and folds returned stats into `ChunkDriverStats`. The
existing `NoSafeGap` retry stays; when triggered, it re-creates
iterators and re-calls the loader with a bigger `max_span` (re-
decompression accepted on this rare path).

The advance condition for the chunk loop becomes
`chunk_range_start = chunk.safe_end` (which it already was) — but
now `chunk.range.end` reflects the loader's actual stop point
rather than the input nominal span.

### Step 5 — CLI knob

Add `--target-variants-per-chunk <N>` to the `var-calling`
subcommand under
[`src/pop_var_caller/var_calling.rs`](../../../src/pop_var_caller/var_calling.rs).
Default: `0` initially (preserves today's behaviour); revisit
after end-to-end perf numbers from step 4.

Help text explains the trade-off: large values = bigger chunks,
fewer per-chunk overheads, more steady-state RSS; small values =
smaller chunks, more overhead per chunk, less RSS.

End-of-run summary line adds `chunk_extensions=N` next to the
existing diagnostic counters.

### Step 6 — byte-identity gate + integration test

The 3-tomato cohort run with `--target-variants-per-chunk 0`
(default off) must remain byte-identical to the prior commit.

A new integration test in
[`tests/cohort_cli_integration.rs`](../../../tests/cohort_cli_integration.rs)
runs the same fixture twice — once with `target_variants_per_chunk =
0`, once with a non-zero value — and asserts equal VCF bodies. The
variant-count knob is a *workload-sizing* knob, not a *correctness*
knob: changing it must never change emitted records.

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
