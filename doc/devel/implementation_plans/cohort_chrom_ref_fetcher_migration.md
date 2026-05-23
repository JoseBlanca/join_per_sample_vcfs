# Cohort + from-bam migration to `ChromRefFetcher`

Follow-up to the
[`unified_chrom_ref_fetcher`](unified_chrom_ref_fetcher.md) plan
(steps 0-3 shipped) and the
[BAQ migration to `ManualEvictChromRefFetcher`](#)
(commit `bcf6e94`).

The earlier work introduced two new fetcher abstractions —
`ChromRefFetcher` (single-contig, monotonic-forward streaming) and
`ManualEvictChromRefFetcher` (single-contig, caller-managed
eviction). Production code routes through both via thin adapters
that translate to the legacy `RefSeqFetcher` trait so the consumer
side (`DustFilter`, `PerGroupMerger`, `drive_cohort_pipeline`,
`PileupWalker`) didn't have to migrate.

This plan completes the consumer-side migration for the cohort
var-calling and `var_calling_from_bam` paths. After it lands:

- `DustFilter`, `PerGroupMerger`, and `drive_cohort_pipeline` work
  directly with the `ChromRefFetcher` trait (no `chrom_id`
  parameter; typed `ChromRefFetchError`).
- `SingleChromLegacyAdapter` deletes — its only consumer
  (`process_one_chromosome` in the cohort var-calling path) calls
  the new trait directly.
- `var_calling_from_bam` restructures to process **one chrom at a
  time** (mirroring the cohort var-calling architecture): chunk the
  walker's coordinate-sorted record stream by `chrom_id`, build a
  per-chrom `StreamingChromRefFetcher`, drive the pipeline per
  chrom, concatenate VCF fragments at the end.
- The legacy `RefSeqFetcher` trait and `WalkerLegacyAdapter`
  survive — but their *only* remaining consumer is the Stage 1
  pileup walker. They're scoped to one specific code path, not
  the whole codebase.

## Why this is worth doing

The two-traits coexistence has carried real friction:

- `SingleChromLegacyAdapter` is a pure-bridging type that exists
  only to translate the legacy multi-chrom API into the new
  single-contig API. Every call passes through a `chrom_id`
  validation that always matches by construction.
- `DustFilter` and `PerGroupMerger` carry a `chrom_id` argument
  through every internal call site that's always equal to the
  fetcher's bound contig. Reading the code requires understanding
  that the parameter is vestigial.
- The cohort var-calling tests' stub fetchers implement the
  legacy trait but the production fetcher is on the new trait;
  reviewers have to mentally translate between them.

Removing the bridge isn't huge code-wise (~150 lines net in the
trait migration, ~300 lines in the from-bam restructure) but
it removes a layer of conceptual indirection that's been
accumulating cost on every review pass.

The Stage 1 walker is deliberately out of scope — see
§"Out of scope" below.

## Spec / supporting documents

- [unified_chrom_ref_fetcher.md](unified_chrom_ref_fetcher.md) —
  the original migration plan; this is its continuation.
- [reference_fasta_streaming.md](reference_fasta_streaming.md) —
  Phases A/B/C, where `StreamingChromRefFetcher` was introduced.
- Pipeline architecture spec:
  [calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md).

## Today's state (the starting point)

```
ChromRefFetcher (new trait)        ChromRefFetchError
  │
  └── StreamingChromRefFetcher  ──── via SingleChromLegacyAdapter ───┐
  └── (no other impls today)                                          │
                                                                      ▼
                                                              RefSeqFetcher (legacy)
                                                                      │
                                  ┌───────────────────────────────────┤
                                  │                                   │
                          WalkerLegacyAdapter                      DustFilter
                                  │                                   │
                                  ▼                                   ▼
                          Stage 1 walker              PerGroupMerger / drive_cohort_pipeline
                          var_calling_from_bam                       │
                                                                     ▼
                                                              process_one_chromosome (cohort)
```

The cohort var-calling path goes through *both* traits via the
adapter. After this migration:

```
ChromRefFetcher
  │
  └── StreamingChromRefFetcher ───────────────────────┐
                                                       │
                                  ┌────────────────────┤
                                  │                    │
                          process_one_chromosome   var_calling_from_bam
                                  │                  (per-chrom loop)
                                  ▼                    │
                          DustFilter / PerGroupMerger /
                          drive_cohort_pipeline ◄──────┘


RefSeqFetcher (legacy — surviving, scoped)
  │
  └── via WalkerLegacyAdapter ───► Stage 1 walker only
```

`SingleChromLegacyAdapter` deletes. `WalkerLegacyAdapter` and the
legacy `RefSeqFetcher` survive but only the Stage 1 walker uses
them. Migrating the walker is a separate slice (see §"Out of
scope").

## Scope

### In scope

**A. Trait migration on the cohort pipeline consumers.**

- `DustFilter<I, F>`: bound `F: RefSeqFetcher` → `F: ChromRefFetcher`.
  Drop `chrom_id` from the internal `fetch` / `iter_bases` calls.
  Add a per-record `chrom_id` validation in `next()` — incoming
  records must match the bound contig (defensive; mismatch is an
  upstream bug). Today's `ensure_mask_for(chrom_id)` becomes
  `ensure_mask_loaded()` — called once on the first record, then a
  cached-mask check on every subsequent record. The
  `chromosomes: Vec<ParsedChromosome>` field is no longer needed
  for length lookup (the fetcher knows its length via
  `ChromRefFetcher::length`); the field may shrink to just the
  bound contig's metadata or disappear entirely.
- `PerGroupMerger`: `Arc<dyn RefSeqFetcher + Send + Sync>` →
  `Arc<dyn ChromRefFetcher + Send + Sync>`. Fetch calls drop
  `chrom_id`. Internal logic unchanged.
- `drive_cohort_pipeline`: `fetcher` parameter type changes;
  internal wiring updates. The `SharedRefFetcher` type alias is
  redefined as `Arc<dyn ChromRefFetcher + Send + Sync>` (or
  renamed if "shared" no longer fits).
- `process_one_chromosome` (cohort var-calling): drops the
  `SingleChromLegacyAdapter` wrap; constructs `StreamingChromRefFetcher`
  and `Arc::new`s it as `Arc<dyn ChromRefFetcher + Send + Sync>`
  directly.
- Test mocks: `DustFilter::tests::StubFetcher` and
  `PerGroupMerger`'s test fetchers migrate to the new trait.
  Mechanical: drop the `chrom_id` parameter from `fetch`, return
  `Result<Vec<u8>, ChromRefFetchError>`.

**B. `var_calling_from_bam` per-chrom restructure.**

- New `PerChromRecordsIter` helper in
  [src/pop_var_caller/var_calling_from_bam.rs](../../src/pop_var_caller/var_calling_from_bam.rs)
  (or a new sibling module). Wraps a peekable iterator of
  `Result<PileupRecord, _>` and yields records for one chrom at
  a time. Has a `consume_current_chrom(&mut self) -> impl
  Iterator<Item = Result<...>> + '_` method that yields records
  until the chrom_id changes, then returns `None` so the caller
  can advance to the next chrom.
- `run_cohort_pipeline_for_single_sample` (the from-bam driver)
  restructures from one-pipeline-call into a per-chrom loop:
  ```rust
  let mut iter = PerChromRecordsIter::new(walker_output);
  let frags_dir = TempDir::new_in(...);
  let fragment_paths = ...;
  while let Some(chrom_id) = iter.open_next_chrom() {
      let fetcher = StreamingChromRefFetcher::for_contig(
          reference, &chromosomes[chrom_id as usize].name,
      )?;
      let fetcher: SharedRefFetcher = Arc::new(fetcher);
      let merger = PerPositionMerger::new(
          vec![Box::new(iter.consume_current_chrom())],
          sample_names.clone(),
          chromosomes.clone(),
      )?;
      let fragment_path = frags_dir.path().join(fragment_name(chrom_id));
      let writer_cfg = writer_cfg_template.clone();
      let writer_cfg = WriterConfig { output: fragment_path.clone(), ..writer_cfg };
      drive_cohort_pipeline(
          merger, params.clone(), fetcher,
          &fragment_path, metadata.clone(), writer_cfg,
      )?;
  }
  concat_fragments(&output_tmp, &fragment_paths, kind)?;
  fs::rename(&output_tmp, &args.output)?;
  ```
  The shape mirrors `var_calling.rs::run_var_calling` step-for-step.
  The walker-error stash (`ErrorSheddingAdapter`) lives outside the
  loop and is read after the loop finishes.
- Walker construction stays multi-chrom (the walker itself doesn't
  change); only the consumption pattern downstream of the walker
  changes.

**C. Delete the bridge.**

- `SingleChromLegacyAdapter` struct + impl + tests delete.
- The `legacy_io_error` helper stays — it's still used by
  `WalkerLegacyAdapter`.
- The cohort path's stale doc comments referencing the adapter
  get cleaned up.

### Out of scope (deferred)

- **Stage 1 pileup walker migration.** Migrating `pileup::run`
  (`PileupWalker<I, F: RefSeqFetcher>`) to `ChromRefFetcher` would
  require restructuring the walker itself to be per-chrom. The
  walker carries cross-position state (open-record table, mate
  lookups, active-read tracking) that doesn't naturally chunk at
  chrom boundaries — restructuring it is a separate ~500+ line
  refactor with its own justification (e.g. enabling per-chrom-
  parallel Stage 1, mirroring the cohort var-calling parallelism).
  Tracked as a future slice.
- **Renaming `RefSeqFetcher` to reflect its narrowed scope** (e.g.
  `WalkerRefSeqFetcher`). The trait survives as the multi-chrom
  walker API; renaming is cosmetic and can ride along whenever the
  Stage 1 walker is touched next.
- **Per-chrom parallelism in from-bam.** The new per-chrom loop is
  serial by design — from-bam is single-sample, so each per-chrom
  run is cheap and parallelism would race against the walker's
  multi-chrom output stream. Cohort var-calling has per-chrom
  parallelism because its inputs (PSPs) are random-access by
  chrom; the walker stream isn't.

### Why not break this into two commits

The DUST/PerGroupMerger/drive_cohort_pipeline trait migration and
the from-bam restructure are tangled — they share
`drive_cohort_pipeline`. If we migrated the trait first, from-bam
wouldn't compile until the per-chrom restructure lands too. The
inverse (restructure from-bam without migrating the trait) means
restructuring around a legacy API we're about to delete.

The cleanest shape is **one Step 1 helper commit** (the
`PerChromRecordsIter` adapter with unit tests) followed by **one
Step 2 commit** that flips the trait migration AND the from-bam
restructure together. Both pieces ship green.

## Design

### Step 1 — `PerChromRecordsIter`

A small iterator adapter that chunks a coordinate-sorted record
stream by `chrom_id`. Sketch:

```rust
pub struct PerChromRecordsIter<I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    upstream: std::iter::Peekable<I>,
    /// Chrom currently being consumed. `None` until the first
    /// `open_next_chrom()` call.
    current_chrom: Option<u32>,
    /// Set when the upstream yielded an error during peek. The
    /// error is held here until the next `consume_current_chrom`
    /// call gives it back to the caller.
    pending_error: Option<PspReadError>,
}

impl<I> PerChromRecordsIter<I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    pub fn new(upstream: I) -> Self;

    /// Peek the next upstream record; if there is one, set the
    /// current chrom to its chrom_id and return `Some(chrom_id)`.
    /// Returns `None` when the upstream is exhausted (clean EOF).
    /// A pending error from a previous peek is returned via the
    /// next `consume_current_chrom`'s first `next()` call.
    pub fn open_next_chrom(&mut self) -> Result<Option<u32>, PspReadError>;

    /// Yield records while their `chrom_id` matches the current
    /// chrom. Returns `None` at chrom transition or EOF — the
    /// caller then loops with `open_next_chrom()`.
    pub fn consume_current_chrom(&mut self) -> ChromRecordsIter<'_, I>;
}

pub struct ChromRecordsIter<'a, I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    parent: &'a mut PerChromRecordsIter<I>,
}

impl<'a, I> Iterator for ChromRecordsIter<'a, I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    type Item = Result<PileupRecord, PspReadError>;
    fn next(&mut self) -> Option<Self::Item> { ... }
}
```

Tests:

1. `per_chrom_iter_yields_chunked_records` — synthetic records on
   chroms 0, 0, 0, 1, 1, 2. Two-chrom-loop yields {0: 3 records,
   1: 2 records, 2: 1 record}.
2. `per_chrom_iter_handles_single_chrom` — all records on chrom 0;
   one open_next_chrom returns Some(0), then None.
3. `per_chrom_iter_handles_empty_input` — empty upstream;
   open_next_chrom returns None on first call.
4. `per_chrom_iter_propagates_upstream_error` — error in the
   middle; surfaces through the inner ChromRecordsIter.
5. `per_chrom_iter_chrom_transition_at_error` — error happens
   exactly at a chrom transition; still surfaces correctly.

### Step 2 — Trait migration + from-bam restructure

**DustFilter migration.**

Today's signature:
```rust
pub struct DustFilter<I, F>
where
    I: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
    F: RefSeqFetcher,
{
    upstream: I,
    fetcher: F,
    chromosomes: Vec<ParsedChromosome>,
    config: DustFilterConfig,
    loaded: Option<LoadedChrom>,
    is_finished: bool,
}
```

New shape:
```rust
pub struct DustFilter<I, F>
where
    I: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
    F: ChromRefFetcher,
{
    upstream: I,
    fetcher: F,
    /// `chrom_id` the fetcher is bound to. Validated against
    /// every incoming record's chrom_id; mismatch is an upstream
    /// invariant violation surfaced as `DustFilterError::ChromIdMismatch`.
    bound_chrom_id: u32,
    config: DustFilterConfig,
    /// Cached mask for the bound contig. Built lazily on the first
    /// record (no need to chase `ensure_mask_for(chrom_id)` per
    /// record now that the fetcher is single-contig).
    mask: Option<SdustIntervals>,
    sweep: u32,
    is_finished: bool,
}
```

The `chromosomes` table goes away (length comes from
`fetcher.length()`). The `bound_chrom_id` is passed at
construction. `DustFilterError` grows a `ChromIdMismatch { expected, got }`
variant for the validation.

Tests in `dust_filter.rs::tests` migrate: `StubFetcher` switches
to implementing `ChromRefFetcher`. Construction of `DustFilter`
gets the new `bound_chrom_id` argument.

**PerGroupMerger migration.**

`SharedRefFetcher` type alias redefined:
```rust
pub type SharedRefFetcher = Arc<dyn ChromRefFetcher + Send + Sync>;
```

`PerGroupMerger.ref_fetcher: SharedRefFetcher`. Fetch calls drop
chrom_id:
```rust
// Before:
ref_fetcher.fetch(chrom_id, start, span).map_err(|source| ...)
// After:
ref_fetcher.fetch(start, span).map_err(|source| ChromRefFetchError → io::Error)
```

The fetch result is `Vec<u8>` either way — the API shape of
`ChromRefFetcher::fetch` returns `Result<Vec<u8>, ChromRefFetchError>`.
PerGroupMerger maps the error into its own error enum
(`PerGroupMergerError::RefFetch`).

**drive_cohort_pipeline migration.**

Signature unchanged except for fetcher type. Internal call sites
to DUST + PerGroupMerger update accordingly. The `chromosomes` arg
passed to DUST is no longer needed (DUST gets length from the
fetcher).

**process_one_chromosome simplification.**

Before:
```rust
let streaming = StreamingChromRefFetcher::for_contig(fasta_path, &chrom_entry.name)?;
let adapter = SingleChromLegacyAdapter::new(chrom_id, streaming);
let fetcher: SharedRefFetcher = Arc::new(adapter);
```

After:
```rust
let streaming = StreamingChromRefFetcher::for_contig(fasta_path, &chrom_entry.name)?;
let fetcher: SharedRefFetcher = Arc::new(streaming);
```

The `chrom_id` argument is still passed to `process_one_chromosome`
(used for fragment naming, error reporting), it's just no longer
threaded through the adapter.

**from-bam restructure.**

Today the from-bam driver builds one walker → one merger → one
`drive_cohort_pipeline` call → one VCF file. The new shape:

```rust
fn run_cohort_pipeline_for_single_sample(
    ctx: Stage1PipelineContext<'_>,
    reference: &Path,
    output: &Path,
    command_line: String,
    ...
) -> Result<(), VarCallingFromBamCliError> {
    let chromosomes = contigs_to_parsed(ctx.contigs)?;
    let sample_name = ctx.sample_name.to_string();

    // Walker error stashing stays at this layer — errors surface
    // after the per-chrom loop ends.
    let walker_adapter = ErrorSheddingAdapter::new(ctx.walker);
    let walker_error_handle = walker_adapter.error_handle();
    let mut per_chrom = PerChromRecordsIter::new(walker_adapter.map(Ok));

    // Tempdir for per-chrom VCF fragments + final concat target.
    let output_parent = output.parent().unwrap_or_else(|| Path::new("."));
    let frags_dir = TempDir::new_in(output_parent)?;
    let frag_ext = if path_is_bgzf(output) { "vcf.gz" } else { "vcf" };
    let mut fragment_paths: Vec<PathBuf> = Vec::new();

    let metadata = CohortMetadata { ... };
    let writer_cfg_template = WriterConfig::new(output.to_path_buf())
        .with_emit_gp(emit_gp);

    let pipeline_params = CohortPipelineParams { ... };

    while let Some(chrom_id) = per_chrom.open_next_chrom()? {
        let entry = &chromosomes[chrom_id as usize];
        let fetcher = StreamingChromRefFetcher::for_contig(reference, &entry.name)?;
        let fetcher: SharedRefFetcher = Arc::new(fetcher);

        let fragment_path = frags_dir.path().join(
            format!("chr_{chrom_id:03}_{name}.{frag_ext}", name = entry.name)
        );
        fragment_paths.push(fragment_path.clone());
        let writer_cfg = WriterConfig {
            output: fragment_path.clone(),
            ..writer_cfg_template.clone()
        };

        let records_iter: Box<
            dyn Iterator<Item = Result<PileupRecord, PspReadError>> + '_,
        > = Box::new(per_chrom.consume_current_chrom());
        let merger = PerPositionMerger::new(
            vec![records_iter],
            vec![sample_name.clone()],
            chromosomes.clone(),
        )?;

        drive_cohort_pipeline::<_, VarCallingFromBamCliError>(
            merger,
            pipeline_params.clone(),
            fetcher,
            &fragment_path,
            metadata.clone(),
            writer_cfg,
        )?;
    }

    // Concat per-chrom fragments to a final tmp, then atomic rename.
    let output_tmp = tmp_path_for(output);
    let kind = SinkKind::from_path(output);
    concat_fragments(&output_tmp, &fragment_paths, kind)?;
    std::fs::rename(&output_tmp, output)?;

    if let Some(e) = walker_error_handle.take() {
        let _ = std::fs::remove_file(output);
        return Err(VarCallingFromBamCliError::Walker(e));
    }
    Ok(())
}
```

`CohortPipelineParams` needs `Clone` (it already does, for the
per-chrom-parallel cohort path).

The walker is built once (multi-chrom CRAM input → multi-chrom
PileupRecord stream); the per-chrom loop consumes it.

**`SingleChromLegacyAdapter` deletion.**

Once the cohort path migrates, the adapter has no consumers.
Delete struct + impl + 4 unit tests.

## Test strategy

- **Step 1**: 5 new unit tests for `PerChromRecordsIter` (listed
  above).
- **Step 2**: existing tests carry the load.
  - DUST's 38 unit tests migrate (the StubFetcher swap).
  - PerGroupMerger's tests migrate similarly.
  - `cohort_cli_integration.rs::var_calling_*` tests: byte-identical
    VCF output expected; the trait migration is pure plumbing.
  - `pileup_cli_integration.rs` (Stage 1): unchanged behaviour;
    Stage 1 still uses the legacy trait.
  - **var_calling_from_bam integration test** (one in
    `cohort_cli_integration.rs`): byte-identical output expected.
    The per-chrom restructure produces the same records in the
    same order; concat reassembles them in contig-table order
    which matches today's continuous-write order on
    coordinate-sorted input.
- Real-data validation: `pop_var_caller var-calling-from-bam`
  on a small CRAM sample after the migration — should produce
  the same VCF as before. Quick diff check, not a perf bench.

## Migration order (commit-by-commit)

1. **PerChromRecordsIter helper + 5 unit tests.** Standalone,
   no production wiring. The new struct exists in
   `var_calling_from_bam.rs` (or a sibling module) but no caller
   uses it yet.

2. **Big migration commit:**
   - DUST: `F: RefSeqFetcher` → `F: ChromRefFetcher`. Drop
     `chromosomes` field. Add `bound_chrom_id` field. Add
     `DustFilterError::ChromIdMismatch`. Update internal logic.
     Migrate `StubFetcher` in tests.
   - PerGroupMerger: similar bound + signature changes. Migrate
     test fetchers.
   - `SharedRefFetcher` type alias redefined.
   - `drive_cohort_pipeline`: signature update; internal wiring.
   - `process_one_chromosome`: drop adapter wrap.
   - `var_calling_from_bam`: per-chrom loop via the helper from
     step 1.
   - Delete `SingleChromLegacyAdapter` (struct + impl + tests).

This is one big diff (~600-800 lines) but it has to ship as one
commit — partial migration leaves the code uncompileable. The
step 1 commit is the only reversibly-shippable preamble.

## Risks and discovery items

1. **DUST's chrom_id mismatch check.** The new defensive check
   could surface bugs in upstream code that today silently emit
   multi-chrom records into a "single-chrom" worker context.
   Treatment: if a test fires the new error, treat it as
   discovery — fix the upstream emission, not the validation.

2. **from-bam VCF output ordering.** The per-chrom restructure
   produces fragments per chrom and concats in contig-table
   order. Today's single-pass write produces records in walker
   order (which is contig-table order for coordinate-sorted
   CRAM). These should match; if a fixture has reads on contigs
   in non-table order, the migrated output differs. Mitigation:
   the contig table is built from the CRAM header
   (coordinate-sorted), so iteration order matches by
   construction. The integration test asserts byte-identical
   output.

3. **walker_error_handle timing.** The ErrorSheddingAdapter sits
   between the walker and the per-chrom iter. Errors during
   `walker.next()` go into the stash and surface after the per-
   chrom loop completes. If an error fires mid-chrom, the loop
   exits early (current chrom's iter returns None), the next
   `open_next_chrom()` returns the stashed error. Verified by
   the existing walker-error integration test (which still uses
   the legacy single-pass shape, but the error semantics are
   identical).

4. **`PerPositionMerger` lifetime.** The merger takes
   `Vec<Box<dyn Iterator + 'static>>` (or similar). The inner
   `consume_current_chrom()` returns a borrowed iter with a
   lifetime tied to the `PerChromRecordsIter`. Box-erasing it
   into a `Box<dyn Iterator + '_>` is straightforward; the merger
   probably needs its lifetime parameter widened or its inputs
   re-typed. To be resolved at implementation time.

5. **Mock fetchers in tests.** DUST tests' `StubFetcher` impls
   `RefSeqFetcher` today. Migration: drop the `chrom_id` arg,
   change the error type. Mechanical. PerGroupMerger tests
   similar. No semantic change.

## Out-of-scope (explicit non-goals)

- Stage 1 walker migration (see §"Out of scope" above).
- Renaming `RefSeqFetcher` to reflect its scoped use.
- Restructuring drive_cohort_pipeline to take more parameters
  (e.g., chrom_id explicitly) — the chrom_id is captured by the
  fetcher's binding; no separate parameter needed.
- Performance tuning of the from-bam per-chrom loop (e.g.,
  parallelism). Out of scope; the single-sample workload is
  small.
- Renaming `SharedRefFetcher` to reflect the new trait. Cosmetic;
  can ride along later.

## Estimated effort

- Step 1: ~150 lines (helper + 5 unit tests). Half day.
- Step 2: ~600-800 lines net (trait migration across DUST +
  PerGroupMerger + drive_cohort_pipeline + process_one_chromosome
  + from-bam restructure + test fetcher migrations + adapter
  deletion). Day.

Total: ~750-950 lines, two commits.
