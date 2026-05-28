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

## Validation status (so far)

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib --features dhat-heap` — 977/977 pass.
- Byte-identity vs `main` — not yet attempted; Step 6 gate.

## Next steps in this conversation

The skill workflow's Step 5 (driver) and Step 6 (validation +
PROJECT_STATUS update) remain. Step 5 wires the new chunk-loop driver
into `run_var_calling`, replacing the per-chromosome
`rayon::par_iter`. Step 6 runs the existing
[tests/cohort_cli_integration.rs](../../../tests/cohort_cli_integration.rs)
suite + a side-by-side VCF diff against `main` on a tomato fixture +
the wall + peak RSS sweep via
[perf_scaling_synthetic.py](../../../benchmarks/tomato1/scripts/perf_scaling_synthetic.py)
at T=1 on N=50 / N=200 / N=1000.
