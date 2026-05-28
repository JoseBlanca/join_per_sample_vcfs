# Cohort within-chromosome chunk-parallel rewrite — Phase A

**Status:** in-flight (this is a living document; updated as steps land)
**Branch:** `cohort-within-chromosome-parallel`
**Plan:** [cohort_within_chromosome_parallel.md](../../implementation_plans/cohort_within_chromosome_parallel.md)
**Motivating measurement:** [scaling_measurement_2026-05-27.md](../reviews/scaling_measurement_2026-05-27.md)

## What this report tracks

Phase A of the architectural rewrite that replaces H1's per-chromosome
`rayon::par_iter` with within-chromosome chunk-based parallelism. The
plan calls Phase A out as: byte-identical VCFs at T=1 against `main` on
the existing integration-test fixtures, with the new chunk loader +
pre-pass + worker pipeline in place but parallel windows deferred to
Phase B. This document accumulates the Phase A landings step by step
and — importantly — every piece of work we deliberately deferred or
simplified so it isn't lost.

## Pre-implementation spike

The plan calls for confirming `posterior_engine`'s M-step is purely
intra-group before committing to per-window EM. Confirmed by reading
[posterior_engine.rs:1723](../../../src/var_calling/posterior_engine.rs#L1723)
(`run_em_for_record`) and the M-step kernels at lines
[2432](../../../src/var_calling/posterior_engine.rs#L2432) and
[2505](../../../src/var_calling/posterior_engine.rs#L2505). The EM
reads only per-record `MergedRecord` fields and the per-record
`RecordScratch`; no cross-record state survives a `next()` call beyond
the scratch buffer (which is reset on every record via `resize_to`).
Safe to run EM per-group inside per-window workers. A regression test
codifying the invariant is queued under "deferred work" below.

## What has landed

### Chunk infrastructure (commit `bc28092`)

New module [`src/var_calling/cohort_block/`](../../../src/var_calling/cohort_block/)
with four submodules:

- [`columns`](../../../src/var_calling/cohort_block/columns.rs) —
  `SampleColumns` (per-sample CSR columnar storage of
  `PileupRecord`-shaped data) + `MaterialisedChunk` (N samples + range
  + safe_end + windows). Scratch-friendly methods: `clear`,
  `truncate`, `drain_rows_from_into`, `push_row_from`, `materialise_record`.
- [`loader`](../../../src/var_calling/cohort_block/loader.rs) —
  `load_chunk_from_iters` reads per-sample record iterators, drains
  carryover from the prior chunk, applies the cohort-wide
  variant-position filter (any sample with a non-REF allele at
  `num_obs > 0`), and compacts survivors column-to-column. Generic
  over the upstream iterator's error type.
- [`pre_pass`](../../../src/var_calling/cohort_block/pre_pass.rs) —
  `fix_boundaries` picks the chunk's `safe_end` (gap >
  `max_group_span` on the post-filter timeline AND no allele below
  reaches into the gap), splits records past `safe_end` into next
  chunk's carryover, partitions `[range.start, safe_end)` into worker
  windows. At T=1 emits one window per chunk; Phase B extends.
- [`partition`](../../../src/var_calling/cohort_block/partition.rs) —
  `partition_window` builds the cohort-wide position view inside one
  window and groups overlapping positions into variant groups. CSR
  columnar output (`WindowPartition`) with `(sample_idx, row_idx)`
  handles per position. Inline DUST mask: caller passes sorted
  half-open masked intervals; masked positions advance without
  touching group state.

60 unit tests across the four submodules + shared `test_helpers`
fixture builders.

### Worker adapter (commit `a094f73`)

[`src/var_calling/cohort_block/worker.rs`](../../../src/var_calling/cohort_block/worker.rs)
— Phase A.0 adapter that builds an `OverlappingVariantGroup` from
each partition group + chunk columns, pipes the resulting iterator
through the existing `PerGroupMerger` and `PosteriorEngine` kernels,
appends `PosteriorRecord`s to the caller's output buffer. Empty
partitions are a fast no-op that never queries the ref fetcher.

4 unit tests cover the adapter shape + the empty-partition fast
path. Kernel correctness is covered by the existing
`PerGroupMerger` and `PosteriorEngine` test suites; end-to-end
byte-identity is the Step 6 validation gate.

## Deliberate simplifications and deferred work

This list is the load-bearing tracker for everything we know about
but chose not to do in Phase A. Each item names what was deferred,
why, and what triggers the follow-up.

### From the rewrite plan

1. **Phase A.1 — rewrite kernels native columnar.**
   Step 4d.0 reuses `PerGroupMerger::process_group` and
   `run_em_for_record` via an adapter that materialises
   `OverlappingVariantGroup`s row-by-row from the chunk columns. The
   plan calls for native-columnar kernels (allele unification +
   likelihood construction + EM) that operate directly on
   `SampleColumns` and `WindowPartition`. Trigger: byte-identity
   confirmed at Step 6; perf review motivates the rewrite.

2. **Phase B — parallel windows.**
   `fix_boundaries` currently only supports `target_window_count = 1`
   (returns `FixBoundariesError::UnsupportedTargetWindowCount` for >
   1). Phase B extends to T-1 internal boundaries placed in safe gaps
   near evenly-spaced positions, with deterministic per-window output
   buffers so byte-identity survives parallel execution.

3. **Phase C — pipelined chunk loading.**
   Background-thread chunk loader: while workers process chunk K, the
   loader prepares chunk K+1. Bounded queue depth 2.

4. **Phase D — SIMD perf-review pass.**
   Per-sample-batch SIMD in `compute_log_likelihoods`, cache-locality
   measurement, allocator A/B, etc. Driven by the
   `rust-performance-review` skill once the architecture is proven.

5. **`GroupedVariantsBatch` columnar refactor.**
   Plan called for refactoring `PerGroupMerger`'s emit to a columnar
   batch accumulator. Deferred — the existing per-emit `MergedRecord`
   allocation is preserved in the worker adapter. Lands with Phase A.1
   or as its own pre-Phase B refactor.

6. **`PileupRecordRef<'a>` borrowed-view threading.**
   Plan called for a borrowed-view type threading through
   `PerPositionMerger` + `VariantGrouper` + `PerGroupMerger`. Phase A
   materialises owned `PileupRecord` rows at the per-sample iterator
   boundary instead (one alloc per record, same cost as `main`'s PSP
   reader). Lands when the kernels go native-columnar (Phase A.1) or
   if SIMD pressure motivates it (Phase D).

7. **`setrlimit(RLIMIT_NOFILE)` raise at driver startup.**
   Plan calls for `max(soft, 4N)` raise with fail-fast on rejection.
   Phase A is single-threaded sequential, so fd budget is N (one
   `region_records` per sample for the current chunk). Equivalent to
   today's per-chrom code. Lands in Phase B/C where the parallel
   chunk loader actually opens > 2N fds.

### From `'static` audit

8. **`SharedRefFetcher` lifetime-clean refactor.**
   `pub type SharedRefFetcher = Arc<dyn ChromRefFetcher + Send>` in
   [per_group_merger.rs:573](../../../src/var_calling/per_group_merger.rs#L573)
   carries an implicit `+ 'static` bound that propagates to every
   consumer (including
   [worker.rs:141](../../../src/var_calling/cohort_block/worker.rs#L141)
   `shared_ref_fetcher`'s `F: ChromRefFetcher + Send + 'static`). The
   fetcher in practice only needs to outlive the cohort-driver call,
   not the program. Refactor candidates: `Arc<dyn ChromRefFetcher +
   Send + 'a>` with the lifetime spelled out, or `&'a dyn
   ChromRefFetcher` if shared-Arc semantics aren't needed at the
   call sites. Own workstream; doesn't gate Phase A correctness.

### From Phase A worker simplifications

The kernels reused via the adapter already implement the full
plant-genome feature set. The simplifications below are about what
the **worker driver** asks the kernels to do, not about what the
kernels are capable of — so they're effectively no-ops in Phase A.0
(the kernels' code paths just aren't exercised) but become real
work items in Phase A.1 when the new kernels are written.

9. **Diploid-only assumption documented but not enforced at the worker level.**
   `PerGroupMergerConfig::ploidy` is still tunable end-to-end. The
   project-default ploidy of 2 is what integration tests exercise.
   Phase A.1 must use `genotype_order(ploidy, n_alleles)` as the
   driver of the genotype enumeration; the diploid fast path is an
   optimisation, not a contract.

10. **Contamination correction kept available via existing kernel.**
    Phase A.0 passes the user's `PosteriorEngineConfig` through
    unchanged, so `contamination = Some(_)` still drives the existing
    mixture-likelihood E-step. Phase A.1 must port the mixture path
    too if `contamination = None` is no longer the only test path.

11. **Compound alleles handled by existing kernel; not by the
    column-side code.**
    `SampleColumns` carries the per-allele `chain_ids` column the
    chain-anchor logic needs (since
    [pileup_record.rs:135](../../../src/pileup_record.rs#L135)
    declares the field), but `partition_window` does not currently
    inspect chain_ids — chain-anchoring happens entirely inside
    `PerGroupMerger`. Phase A.1 needs to either port the
    chain-anchoring logic or surface a column-side helper for it.

12. **SIMD math backend availability.**
    `PosteriorEngine` defaults to `InterpUnivariateSimdMath`; Phase
    A.0 keeps that default. If Phase A.1 rewrites the EM, it must
    decide whether to port the SIMD lanes or start scalar-and-iterate.

### From the variant-filter design

13. **Variant-filter `ChromRefFetcher` integration.**
    The plan said the cohort-wide variant-position filter "uses
    `ChromRefFetcher` for the per-position variant predicate." Phase
    A uses the walker invariant `alleles[0] == REF` to collapse the
    predicate to "any allele at index ≥ 1 with `num_obs > 0`", which
    is the plan's authorised fast path. No ref-base fetch needed.
    Documented here so a future reviewer doesn't reach for the
    fetcher.

### From the pre-pass safe-gap design

14. **No-safe-gap retry-with-extended-range.**
    `fix_boundaries` returns `NoSafeGap` when the chunk has no
    qualifying gap. The plan calls for the driver to extend the
    chunk's load range and retry (cap at 4× nominal size). Driver-level
    concern; lands with Step 5.

### From the kernel-reuse decision

15. **Per-emit allocation pattern in the worker adapter.**
    `build_overlapping_variant_group` allocates a fresh
    `Vec<Option<PileupRecord>>` per cohort position and a fresh
    `Vec<PerPositionPileups>` per group. Matches what the streaming
    pipeline does on `main`. Phase A.1's native-columnar kernels
    eliminate both.

16. **Intra-group EM invariance regression test.**
    Pre-spike confirmed `posterior_engine`'s M-step is purely
    intra-group, but the invariant is not codified by a regression
    test. Adding one is queued: an EM run on a single
    `MergedRecord` whose log-likelihoods carry the same per-sample
    values should yield identical posteriors whether the record is
    fed alone or interleaved with other groups. Lands with Step 6 or
    earlier if a regression risk surfaces.

### Chunk-loop driver (commit `acd2107`)

[`src/var_calling/cohort_block/driver.rs`](../../../src/var_calling/cohort_block/driver.rs)
— `drive_cohort_chunked` + `drive_one_chrom_generic` (`W: Read + Seek`
generic so tests + production share one body). Owns the persistent
chunk-loop scratch across every chunk and every chromosome
(`ChunkLoadScratch`, `FixBoundariesScratch`, `PartitionScratch`,
`MaterialisedChunk`, `WindowPartition`, the per-sample carryover, the
per-window `Vec<PosteriorRecord>`) and the single `CohortVcfWriter`
shared across chromosomes. PSP readers opened once, reused via
`region_records` (block-index seek per call); no per-chrom file open.

Per-chrom DUST mask built once from the whole-chrom reference via
[`sdust_mask_streaming`](../../../src/var_calling/dust_filter.rs)
(translated from 0-based half-open slice coords to 1-based half-open
genomic), reused across every window on the chrom. Downstream filters
(`min_alt_obs_per_sample`, `is_variant_call`, `qual_phred`, MAPQ-diff
t-test in `record_fails_mapq_diff_t`) applied post-EM in
`emit_or_drop`, mirroring [`drive_cohort_pipeline`](../../../src/pop_var_caller/cohort_driver.rs)'s
drop semantics so the counters in `ChunkDriverStats` match
`CohortDriveStats` field for field.

The driver uses one notable trick to keep `fix_boundaries` unchanged
at the chrom end: the last chunk inflates `chunk.range.end` past
`chrom_length + max_group_span` so the pre-pass's case-A gap check
trivially succeeds without `fix_boundaries` needing to be
chromosome-end-aware. The PSP iterator still caps reads at
`chrom_length` (`region_records(chrom_id, psp_cursor,
min(chunk.range.end - 1, chrom_length))`), so no record is ever read
past the chromosome.

Two per-chunk state values track the offset semantics:
- `chunk_range_start`: lower bound of the next chunk's logical range
  (carryover positions live in `[chunk_range_start, psp_cursor)`).
- `psp_cursor`: first position the PSP iterator picks up (records
  before this were already consumed in a previous chunk or are in
  carryover).

These diverge after the first chunk; `chunk_range_start` advances by
`safe_end` per iteration, `psp_cursor` advances by the chunk's actual
`range.end`.

### Rewiring `run_var_calling` (commit `2535ebd`)

[`src/pop_var_caller/var_calling.rs`](../../../src/pop_var_caller/var_calling.rs)
— the H1 per-chromosome `rayon::par_iter` at line 397 is gone.
`process_one_chromosome`, the per-chrom `TempDir` + fragment paths,
the `concat_fragments` call, and `CohortPipelineParams` are all
removed from the var-calling path. A single `drive_cohort_chunked`
call replaces them; `chunk_stats_to_cohort_stats` converts the
driver's [`ChunkDriverStats`](../../../src/var_calling/cohort_block/driver.rs)
into the [`CohortDriveStats`](../../../src/pop_var_caller/cohort_driver.rs)
shape `print_run_summary` consumes.

`var-calling-from-bam` is unaffected — it still uses the streaming
`process_one_chromosome` + bgzf-aware `concat_fragments` (the
`crate::vcf::concat` module stays).

Net diff for the rewire: +54 / -79 in [var_calling.rs](../../../src/pop_var_caller/var_calling.rs).

## Validation status

- `cargo fmt --check` — clean across the branch.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib --features dhat-heap` — **977/977 pass.**
- `cargo test --test cohort_cli_integration --features dhat-heap`
  — **18/18 pass.** These tests run `var-calling` end-to-end through
  the new chunk-loop driver, including the
  `var_calling_emits_deterministic_vcf_across_runs` test that diffs
  two runs of the same input.

**Byte-identity vs `main`** — not yet attempted. The integration suite
runs through the new driver and produces a well-formed cohort VCF on
every fixture, but the existing tests assert presence of header
lines + sample-column names, not per-record VCF body content. The
side-by-side `diff` against `main`'s streaming pipeline output on a
real tomato fixture is the genuine byte-identity check and still owed.

## Decision on perf review for Phase A

**Deferred until Phase A.1 lands.** The current worker
([`worker.rs`](../../../src/var_calling/cohort_block/worker.rs)) is
Phase A.0 — it builds an `OverlappingVariantGroup` row-shape from the
columnar chunk and pipes it through the unchanged `PerGroupMerger` +
`PosteriorEngine` kernels. Measuring wall + peak RSS at this point
would be measuring the chunk-loop overhead *on top of* the
unchanged per-record kernels, not the architecture the rewrite was
sold on. The scaling-measurement report's bottleneck buckets
(`per_position_merger` at 45 % at N=1000, `per_group_merger` at 58 %,
allocator at 38 %) are all things the existing kernels generate —
the new path inherits all of them via the adapter.

A perf review only delivers actionable findings once the kernels
themselves run column-native. That's Phase A.1.

## Recommended next steps

Two concrete pieces in this order:

### 1. Byte-identity diff vs `main` (small, one-shot)

Worth doing before Phase A.1 to lock in a known-good baseline:

- Pick a real tomato fixture (the `SRR11450568.p1.psp × N=10`
  used by the scaling measurement is a good size).
- Run `cargo run --release -- var-calling …` on `main` and on
  this branch with identical args.
- `bcftools view` / `diff` the output VCFs. They should be
  byte-identical for Phase A.0's contract to hold.

If they aren't, that's a Phase A.0 bug to chase before touching
the kernels. If they are, Phase A.1 can begin from a confirmed
baseline.

### 2. Phase A.1 — native columnar kernels (the substantive next step)

The Phase A.0 adapter at
[`worker.rs:101 build_overlapping_variant_group`](../../../src/var_calling/cohort_block/worker.rs#L101)
materialises a row-shaped `OverlappingVariantGroup` per partition
group (with a fresh `Vec<Option<PileupRecord>>` per cohort position),
then feeds it through the streaming kernels. Phase A.1 replaces this
end to end with kernels that walk the columnar chunk + partition
directly.

The work splits naturally into four layers, each portable as its own
step (with the existing kernel as the reference oracle for the
diff at every step):

1. **Native allele unification** — walk `chunk.per_sample[s]`'s
   `allele_seq_*` columns over the group's position list; produce a
   `MergedAlleleColumns` (unified allele list, sorted by `(start,
   end, seq_bytes)`). No `PileupRecord` materialisation.
2. **Native per-(sample, allele) stats gather** — for each sample,
   walk its per-allele scalar columns at the group's positions and
   project onto the unified allele set. Output: a columnar
   `PerGroupSampleStats` (the 7 `AlleleSupportStats` scalars +
   per-allele chain_ids, all CSR-laid-out).
3. **Native log-likelihood computation** — fill a flat
   `log_likelihoods[sample × n_genotypes]` column from the per-sample
   stats. Reuse the standard-likelihood and mixture-likelihood
   formulas from `posterior_engine` (the math is well-defined; only
   the input shape changes).
4. **Native EM** — the M-step on `p̂` and `f̂_C` is closed-form +
   intra-group (spike confirmed). Port the EM loop to read from the
   columnar layout directly; reuse the E-step's HWE-with-F prior and
   the QUAL-via-exact-AF marginalisation. Emit a `PosteriorRecord`
   per group (or, eventually, push into a columnar `PosteriorBatch`
   that the writer drains).

Each layer ships on its own commit + tests against the existing
kernel's output on the same input (byte-identity at every layer
boundary). The first layer (allele unification) is the highest-risk
piece because the chain-anchoring of compound alleles is intricate;
isolating it as its own step keeps the diff reviewable.

**Phase A.1 deliverables:**
- New module `src/var_calling/cohort_block/kernels/` (or similar)
  hosting the four layers as submodules.
- `run_window` switches from the row-shape adapter to the native
  kernel chain.
- The Phase A.0 adapter (`build_overlapping_variant_group` +
  `run_window` row-pipe) is removed once integration tests still pass.
- Perf review (queued behind this).

### 3. Then perf review + PROJECT_STATUS final update

Once Phase A.1 lands and integration tests pass byte-identical, the
perf review delivers actionable signal:

- Re-run [`perf_scaling_synthetic.py`](../../../benchmarks/tomato1/scripts/perf_scaling_synthetic.py)
  on N=50/200/1000 vs `main`.
- Update the impl report with wall + peak RSS deltas.
- Update PROJECT_STATUS's "Last completed task" with final numbers.
