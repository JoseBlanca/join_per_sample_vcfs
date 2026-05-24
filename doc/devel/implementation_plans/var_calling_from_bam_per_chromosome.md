# `var-calling-from-bam` — per-chromosome parallelism

Implementation plan for the per-chromosome parallelisation of the
`var-calling-from-bam` subcommand — the "direct" CRAM/BAM → VCF route
that fuses Stage 1 (pileup) and Stages 3–6 (cohort var-calling) without
writing an intermediate `.psp`.

> **Module-path note (2026-05-24).** The plan was originally written
> against the pre-refactor tree. A large module restructure landed on
> 2026-05 (commits `b51d02b`, `99c8e7f`, `5b9e420`, `01de8b7`, …,
> `2660e68`) that flattened `src/per_sample_pileup/` into top-level
> peers (`src/bam/`, `src/baq/`, `src/fasta/`, `src/psp/`, `src/vcf/`,
> `src/pileup/{walker,per_sample}/`, `src/pileup_record.rs`). All
> file/module references below are post-refactor. The design itself
> did not change — only the call-sites.

Today's reality (one sample's CRAM/BAM in, one VCF out, fully serial):

```
[one sample's CRAM(s)]
   │  (single CramMergedReader, k-way merge over inputs)
   ▼
[BAQ chunk pass]   ← already rayon-parallel within a 1024-read chunk
   │
   ▼
[PileupWalker]     ← one state machine over all contigs, serial
   │
   ▼
[PerChromRecordsIter → DUST → grouper → per-group merger → posterior]
   │  ← all serial: one chrom's records leave the walker before the
   │    next chrom's records arrive
   ▼
[CohortVcfWriter → output.vcf.gz]
```

The decomposition lever is the same as cohort `var-calling`'s
[H1](./cohort_per_chromosome_parallel.md): genomes have many
independent contigs; one rayon worker per contig drives an end-to-end
fused pipeline; per-contig VCF fragments concat in contig-table order
at the end. The difference is the *front* of the chain — instead of
sample-parallel `PspReader::region_records`, each worker opens a
per-contig CRAM/BAM query through a new index-based variant of
`CramMergedReader`.

Expected payoff: a 6–10× wall-time reduction on a 13-chromosome input
(theoretical max ~13×; ch00 / decoy-contig read-imbalance — the same
ceiling the cohort H1 hit at 3.85× on tomato — applies here too).
The from-bam path additionally avoids one full PSP round-trip on disk,
so the absolute wall-time win is larger than the H1 multiplier alone
suggests.

## Spec / supporting documents

- Cohort H1 plan (structural template for this work):
  [cohort_per_chromosome_parallel.md](./cohort_per_chromosome_parallel.md).
  Stages 3–6 are identical; the per-chrom worker shape, fragment
  ordering, and concat reuse the cohort-side machinery verbatim.
- Cohort H1 impl report:
  [cohort_per_chromosome_parallel_2026-05-20.md](../reports/implementations/cohort_per_chromosome_parallel_2026-05-20.md).
- Pipeline architecture spec:
  [calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md) —
  Stage 1 contract + Stage 3–6 per-contig independence.
- Existing per-chrom concat (reused as-is):
  [src/vcf/concat.rs](../../src/vcf/concat.rs).
- Existing per-chrom worker (downstream half, reused):
  [`process_one_chromosome`](../../src/pop_var_caller/cohort_driver.rs#L431)
  in `cohort_driver.rs`. The from-bam worker reuses everything from
  `DUST → … → fragment writer` and only differs in how the
  per-position record stream is constructed.

## Scope

### In scope

- New top-level flag `--build-map-file-index` on
  `var-calling-from-bam` (off by default; opt-in to auto-create
  missing `.crai` / `.bai` / `.csi`).
- New pre-flight helper that validates every input CRAM/BAM has an
  alignment index (or builds it when the flag is set) **before** the
  rayon pool starts. Lives in
  [src/bam/](../../src/bam/) as a new sibling module
  (`src/bam/index_preflight.rs`), since alignment-index management is
  a `bam`-domain concern, not a CLI concern.
- New index-based query variant of
  [`CramMergedReader`](../../src/bam/cram_input.rs) —
  takes a contig and returns a per-contig record iterator backed by
  per-input `IndexedReader::query`. The k-way peek-and-scan merge
  logic is unchanged; only the per-input record source changes.
- New per-chrom worker `process_one_chromosome_from_bam` placed in
  [src/pop_var_caller/var_calling_from_bam.rs](../../src/pop_var_caller/var_calling_from_bam.rs)
  next to `run_cohort_pipeline_for_single_sample`. Reuses
  [`drive_cohort_pipeline`](../../src/pop_var_caller/cohort_driver.rs#L198)
  for the DUST → … → writer half.
- Reshape of
  [`run_var_calling_from_bam`](../../src/pop_var_caller/var_calling_from_bam.rs#L176)
  to: (1) pre-flight indexes, (2) read SAM headers + load indexes
  once into `Arc`, (3) `frags_dir + fragment_paths` (mirror of cohort
  `var-calling`), (4) `chromosomes.par_iter().enumerate().map(…)`,
  (5) reuse `concat_fragments` from `src/vcf/concat.rs`,
  (6) atomic rename via the existing writer-side discipline.
- Two new typed-error variants on `VarCallingFromBamCliError`
  (`MissingMapFileIndex` and `IndexBuildFailed`) with rich
  user-facing messages that name both the flag and the
  `samtools index` recipe.
- Integration tests in
  [tests/cohort_cli_integration.rs](../../tests/cohort_cli_integration.rs)
  — golden-comparison serial-vs-parallel, single-chrom degenerate,
  empty-chromosome-in-middle, missing-index-without-flag error,
  missing-index-with-flag builds successfully, pre-flight catches
  before any parallel work.
- Unit tests for the new index pre-flight helper and the new
  `CramMergedReader` query variant.

### Out of scope

- **Path A — per-chrom parallelism for the `pileup` subcommand.**
  Explicitly dropped in the design discussion: the typical PSP
  workflow runs N independent samples and the orchestrator already
  provides the per-sample parallelism. A single `pileup` invocation
  processes exactly one sample and stays serial. No new `psp::concat`
  module is built in this slice.
- **Multi-sample orchestration inside `var-calling-from-bam`.** The
  subcommand keeps its current contract: one sample's CRAMs in, one
  VCF out. Multi-sample workflows compose at the orchestrator layer.
- **Streaming concat / block-level bgzf concat.** Same v2 deferrals
  as the cohort H1 plan; v1 waits for all per-chrom workers to join
  then concats serially via the existing decompress + re-encode
  helper.
- **Sub-chromosome decomposition** (split a giant chromosome / ch00-
  style decoy across workers). Same ch00 ceiling will apply here as
  on the cohort side; address as a separate slice if profiling
  warrants.
- **`RLIMIT_NOFILE` bump.** v1's test fixtures are tomato-scale (≤13
  chromosomes × a handful of CRAMs per sample = well under the 1024
  default). Track as a follow-up note for very-many-CRAM samples.
- **Parallel index building.** If multiple input CRAMs need indexing
  on the same invocation, build them sequentially. Index building is
  I/O-bound; parallel builds would thrash the disk and it is a
  one-time cost per file.
- **`.crai` / `.bai`/`.csi` freshness checks.** If an index file
  exists we trust it. Stale-index detection is the user's problem —
  same convention as samtools / bcftools.
- **`var-calling-from-bam` summary-line changes.** The
  `effective_threads` parenthetical that the cohort H1 added to the
  `var-calling` summary applies here too. Add as a small follow-up
  edit inside the from-bam summary printer; not the focus of this
  slice.

## Design

### 1. CLI surface — `--build-map-file-index`

Added as a **top-level** field on
[`VarCallingFromBamArgs`](../../src/pop_var_caller/var_calling_from_bam.rs#L72)
— **not** on the shared
[`Stage1Args`](../../src/pop_var_caller/cli/shared_args.rs#L69),
because the `pileup` subcommand does not need indexes (Path A is out
of scope) and putting it on `Stage1Args` would expose a flag on
`pileup` that does nothing.

```rust
/// If any input CRAM/BAM lacks its alignment index (.crai for CRAM;
/// .bai or .csi for BAM), build it in place next to the source file
/// before running. Without this flag, missing indexes are a hard
/// error.
///
/// Default: off. Per-chromosome parallelism requires an index for
/// each input; we will not silently create files near user inputs.
#[arg(long)]
pub build_map_file_index: bool,
```

Naming: `--build-map-file-index` was chosen over `--build-missing-index`
to make it unambiguous that this concerns alignment ("map") files
rather than VCF or FASTA indexes (which the project also reads but
never auto-creates). This sets the precedent for any future
filesystem-write side-effect flag in `pop_var_caller`; subsequent
flags should follow the `--build-<artefact-class>-index` pattern.

### 2. Index pre-flight — `src/bam/index_preflight.rs` (new module)

New sibling of `cram_input` and `errors` under
[src/bam/](../../src/bam/):

```rust
/// Pre-flight every input CRAM/BAM to confirm an alignment index is
/// available. For each input:
///   - CRAM: look for `<path>.crai`.
///   - BAM:  look for `<path>.csi`, fall back to `<path>.bai`.
///
/// If an index exists, do nothing (we trust whatever is on disk;
/// stale-index detection is the user's problem).
///
/// If an index is missing:
///   - `build_if_missing = true` → build it in place next to the
///     source file, printing a single progress line to stderr per
///     file ("building index for sample.cram (1 of 3)..."). On write
///     failure (e.g. read-only directory), return
///     `IndexBuildFailed { path, source }`.
///   - `build_if_missing = false` → return
///     `MissingAlignmentIndex { path, expected_index_path }` with an
///     error message that names both the `samtools index` recipe and
///     leaves it to the CLI layer to add the
///     `--build-map-file-index` mention (the CLI knows the flag
///     name; the `bam` module shouldn't).
///
/// Returns once every input has an index (existing or just built).
/// All work happens before any rayon pool starts; if this returns
/// `Err`, no parallel workers have spawned and no fragment files
/// have been created.
pub fn preflight_alignment_indexes(
    inputs: &[PathBuf],
    build_if_missing: bool,
) -> Result<(), AlignmentIndexError>
```

The bam-module-level typed error lives in
[src/bam/errors.rs](../../src/bam/errors.rs) alongside `CramInputError`:

```rust
#[derive(Error, Debug)]
pub enum AlignmentIndexError {
    /// One of the input alignment files has no index.
    #[error(
        "input alignment file '{path}' has no index \
         (looked for '{expected_index_path}')"
    )]
    MissingAlignmentIndex {
        path: PathBuf,
        expected_index_path: PathBuf,
    },

    /// `build_if_missing` was true but index construction failed.
    #[error("failed to build alignment index for '{path}': {source}")]
    BuildFailed {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
}
```

The CLI layer (`VarCallingFromBamCliError`) wraps these into
two user-facing variants whose messages also point at the flag:

```rust
#[error(
    "input alignment file '{path}' has no index (looked for '{expected_index_path}')\n\
     per-chromosome parallelism requires an index for each input.\n\
     either:\n\
       - re-run with --build-map-file-index to have pop_var_caller build it, or\n\
       - run `samtools index {path}` yourself before invoking this command."
)]
MissingMapFileIndex {
    path: PathBuf,
    expected_index_path: PathBuf,
},

#[error("failed to build alignment index for '{path}': {source}")]
IndexBuildFailed {
    path: PathBuf,
    #[source]
    source: std::io::Error,
},
```

Splitting along the bam-module / CLI boundary keeps the `bam` module
free of `--flag` strings (which are CLI vocabulary), while the CLI
controls the user-facing wording. A `From<AlignmentIndexError>` on
the CLI error widens the two `bam::errors` variants into the two
CLI-vocab variants when the pre-flight return value bubbles up.

The progress line uses `eprintln!` — explicit, user-requested action,
single line per build, no recurring chatter. This is consistent with
the project's "no logs" preference (which targets silent-fallback
sites, not interactive progress on explicitly requested work).

**Index-construction call.** Use noodles' built-in builders. Exact
function names need verification against the pinned noodles versions
during commit 1 implementation:

- CRAM → `noodles_cram::crai::fs::write(path)` (or the equivalent
  builder; verify at the current `noodles-cram 0.93` API).
- BAM → `noodles_csi::fs::write(path)` for CSI (preferred over BAI
  for genomes with contigs > 512 Mb; BAI is fine for tomato but CSI
  is the universal answer).

Either way, the index is written via a tempfile-and-rename inside
noodles or by us, so a partial write on (e.g.) ctrl-C doesn't leave
a corrupt index file the next run silently trusts.

### 3. Per-contig CRAM/BAM query — new `CramMergedReader` variant

The existing
[`CramMergedReader::new`](../../src/bam/cram_input.rs)
opens N input CRAMs and performs a streaming k-way merge across them
over the *whole* sample. The new variant queries one contig:

```rust
impl CramMergedReader {
    /// Build a per-contig merged reader over the same input set.
    /// Each input is opened as an `IndexedReader` and queried over
    /// `(contig_name, 1..=contig_length)`. The k-way merge logic is
    /// the same as `new`; only the per-input record source differs.
    ///
    /// The caller passes the shared, already-loaded
    /// `Arc<sam::Header>` per input and `Arc<crai::Index>` per CRAM
    /// (or `Arc<dyn BinningIndex>` per BAM). Index loading is done
    /// once at pre-flight; cloning the `Arc` per worker is one
    /// atomic.
    ///
    /// Returns an iterator with the same record type as `new`. The
    /// returned reader is *not* `Send` (noodles' `Query` iterator
    /// borrows the underlying reader mutably), so each worker
    /// instantiates its own; the shared state above is the only
    /// thing that crosses worker boundaries.
    pub fn query(
        crams: &[PathBuf],
        headers: &[Arc<sam::Header>],
        indexes: &[AlignmentIndex],
        contig: &str,
        config: CramMergedReaderConfig,
    ) -> Result<Self, CramInputError>
}

/// Lives in `src/bam/index_preflight.rs` (or a sibling `index.rs`).
pub enum AlignmentIndex {
    Crai(Arc<noodles_cram::crai::Index>),
    Bai(Arc<noodles_bam::bai::Index>),
    Csi(Arc<noodles_csi::Index>),
}
```

**Why a new constructor, not a parameter on `new`.** The streaming
constructor's k-way merge calls `read_record_into(&mut)` directly on
each underlying CRAM reader; the indexed version goes through
`Reader::query(&header, &index, &region)` which returns a `Query<'_,
'_, '_, R>` iterator that holds a `&mut Reader<R>` for the iteration's
lifetime. The two record-source shapes are different enough that a
single constructor parameterised over both would be uglier than two
constructors with the merge logic factored out into a shared
private helper. Implementation can revisit this if the duplication
turns out to be small.

**Where the indexes are loaded.** The pre-flight step
(§2) does the index *file* existence check / build; it does **not**
load index bytes into memory. The actual `crai::fs::read` /
`bai::fs::read` / `csi::fs::read` calls happen between pre-flight and
the rayon `par_iter`, in `run_var_calling_from_bam`, into a
`Vec<AlignmentIndex>` (one entry per input), each wrapped in `Arc`.
Headers are likewise read once per input into `Arc<sam::Header>`. Both
are shared read-only across all per-chrom workers.

**File-descriptor footprint.** Each per-chrom worker opens its own
`IndexedReader` per input — `n_chroms × n_inputs` open file handles
during the parallel section. For tomato N=1 sample × 1–4 CRAMs × 13
chroms = 13–52 fds. Well under the 1024 default `RLIMIT_NOFILE`.
Multi-CRAM-per-sample workflows at extreme scale (e.g. 256 CRAMs per
sample × 13 chroms = 3328 fds) would exceed the default; same
follow-up note as on the cohort side.

### 4. Per-chrom worker — `process_one_chromosome_from_bam`

New helper in
[src/pop_var_caller/var_calling_from_bam.rs](../../src/pop_var_caller/var_calling_from_bam.rs)
next to `run_cohort_pipeline_for_single_sample`. Signature mirrors
the cohort-side
[`process_one_chromosome`](../../src/pop_var_caller/cohort_driver.rs#L431):

```rust
#[allow(clippy::too_many_arguments)]
fn process_one_chromosome_from_bam(
    chrom_id: u32,
    contig: &ParsedChromosome,
    crams: &[PathBuf],
    headers: &[Arc<sam::Header>],
    indexes: &[AlignmentIndex],
    reference: &Path,
    cram_cfg: CramMergedReaderConfig,
    baq_cfg: BaqConfig,
    walker_cfg: WalkerConfig,
    baq_chunk_size: usize,
    no_baq: bool,
    fragment_path: PathBuf,
    metadata: CohortMetadata,
    writer_cfg_template: WriterConfig,
    pipeline_params: CohortPipelineParams,
) -> Result<(u32 /* chrom_id */, CohortDriveStats), VarCallingFromBamCliError>
```

Body:

1. Open a per-worker `CramMergedReader::query(crams, headers, indexes,
   &contig.name, cram_cfg)` — yields records for *just this chrom* in
   coordinate order.
2. Build the per-worker BAQ stream around it
   (`crate::pileup::per_sample::baq_stream::BaqStream::new`, or
   `prepare_passthrough` from `baq_engine` when `no_baq`). BAQ's own
   `par_drain` uses the same global rayon pool — see §6 on
   nested-rayon below.
3. Build the per-worker `MultiChromStreamingRefFetcher` from
   [`crate::fasta`](../../src/fasta/mod.rs) over the reference FASTA,
   scoped to this one contig.
4. Build the `PileupWalker` from
   [`crate::pileup::walker`](../../src/pileup/walker/mod.rs) over the
   (BAQ-shimmed) record stream, yielding `PileupRecord`s from
   [`crate::pileup_record`](../../src/pileup_record.rs).
5. Build a single-contig `PerPositionMerger` over the walker output.
   With only one sample, the per-position merger sees at most one
   entry per (chrom, pos); the cross-sample tie-break the cohort
   side worries about does not apply here.
6. Build the per-chrom `WriterConfig` from the template with
   `output = fragment_path`.
7. Call the existing
   [`drive_cohort_pipeline`](../../src/pop_var_caller/cohort_driver.rs#L198)
   with the per-chrom merger + the supplied params + the fragment-path
   writer config + the full cohort metadata.
8. Return `(chrom_id, stats)` on success; bubble the typed error
   otherwise.

The walker / BAQ wiring already lives in
[`run_cohort_pipeline_for_single_sample`](../../src/pop_var_caller/var_calling_from_bam.rs#L325)
and the
[`with_stage1_pipeline`](../../src/pop_var_caller/stage1_pipeline.rs)
helper. Most of the per-worker body is a copy of that flow with two
differences: (a) the record source is `CramMergedReader::query`
instead of `CramMergedReader::new`; (b) it writes a fragment, not the
final output. Worth refactoring the shared chunk into a small private
helper rather than two near-copies; the orchestrator-vs-worker split
naturally falls along that line.

### 5. Driver reshape — `run_var_calling_from_bam`

Current shape
([var_calling_from_bam.rs:176-297](../../src/pop_var_caller/var_calling_from_bam.rs#L176)):
configure rayon → build configs → `with_stage1_pipeline` (which opens
one merged reader, one walker) → `run_cohort_pipeline_for_single_sample`
(which drives DUST → … → writer serially over all contigs).

New shape:

```rust
pub fn run_var_calling_from_bam(
    args: &VarCallingFromBamArgs,
) -> Result<(), VarCallingFromBamCliError> {
    // 1. (unchanged) rayon pool.
    // 2. (unchanged) build configs.
    // 3. (new) pre-flight alignment indexes.
    crate::bam::index_preflight::preflight_alignment_indexes(
        &args.crams,
        args.build_map_file_index,
    )
    .map_err(VarCallingFromBamCliError::from)?;  // From<AlignmentIndexError>

    // 4. (new) load every input's SAM header + alignment index once.
    //    Wrap each in Arc for cross-worker sharing.
    let headers: Vec<Arc<sam::Header>> = load_headers(&args.crams)?;
    let indexes: Vec<AlignmentIndex> = load_indexes(&args.crams)?;

    // 5. (new) reconcile contig list across inputs (this currently
    //    happens inside CramMergedReader::new; lift it out so we have
    //    a chromosomes table before opening per-chrom readers).
    let canonical_contigs = reconcile_contigs(&args.crams, &headers)?;
    let chromosomes: Vec<ParsedChromosome> =
        contigs_to_parsed(&canonical_contigs)?;

    // 6. (unchanged) FASTA MD5 verify against the contig table.

    // 7. (new) tempdir + fragment paths — same shape as the cohort
    //    driver. Lives next to args.output so the final atomic rename
    //    stays on one filesystem.
    let frags_dir = TempDir::new_in(output_parent(&args.output))
        .map_err(VarCallingFromBamCliError::Io)?;
    let frag_ext = if path_is_bgzf(&args.output) { "vcf.gz" } else { "vcf" };
    let fragment_paths: Vec<PathBuf> = chromosomes.iter().enumerate()
        .map(|(cid, c)| frags_dir.path().join(
            format!("chr_{cid:03}_{name}.{frag_ext}", name = c.name)))
        .collect();

    // 8. (new) cohort metadata + writer config template (single sample).
    let metadata = CohortMetadata {
        sample_names: vec![sample_name.clone()],
        contigs: chromosomes.clone(),
        tool_string: format!("pop_var_caller {}", env!("CARGO_PKG_VERSION")),
        command_line: current_command_line(),
    };
    let writer_cfg_template =
        WriterConfig::new(args.output.clone()).with_emit_gp(emit_gp);

    // 9. (new) parallel per-chrom drive.
    let per_chrom_results: Vec<_> = chromosomes.par_iter().enumerate()
        .map(|(cid, contig)| process_one_chromosome_from_bam(
            cid as u32,
            contig,
            &args.crams,
            &headers,
            &indexes,
            &args.reference,
            cram_cfg.clone(),
            baq_cfg.clone(),
            walker_cfg.clone(),
            stage1.baq_chunk_size,
            stage1.no_baq,
            fragment_paths[cid].clone(),
            metadata.clone(),
            writer_cfg_template.clone(),
            pipeline_params.clone(),
        ))
        .collect();

    // 10. (new) fail-fast surfacing of the first error, aggregate stats.
    let mut total_stats = CohortDriveStats::default();
    for r in per_chrom_results {
        let (_cid, stats) = r?;
        total_stats.records_written += stats.records_written;
        // ... etc, same as cohort_driver.
    }

    // 11. (new) concat fragments in contig-table order.
    crate::vcf::concat::concat_fragments(&args.output, &fragment_paths)?;
    // frags_dir RAII-drops here.

    // 12. (unchanged-ish) run summary; mirror the cohort `var-calling`
    //     summary's effective_threads parenthetical.
    print_run_summary(/* ... */);
    Ok(())
}
```

The old `with_stage1_pipeline` + `run_cohort_pipeline_for_single_sample`
helpers either become per-chrom or are deleted; pick during
implementation based on whether the from-bam path still has a serial
fallback (it should not — once indexes are required and present, the
parallel path is always used; threads=1 produces one worker that does
all contigs sequentially, which is the correct degenerate case).

### 6. Nested rayon — outer per-chrom vs. inner BAQ chunk

Both layers share the global rayon pool. On the cohort side, the
analogous nested case (outer per-chrom + inner per-group merger
`par_iter`) was resolved by **dropping the inner par_iter** (L1) —
the inner work was <1 % and rayon-bridge overhead was net-negative
on real data.

The from-bam case is different: the inner par_iter is **BAQ chunk
parallelism** in [`baq_stream::BaqStream`](../../src/pileup/per_sample/baq_stream.rs),
which is genuinely expensive (BAQ is a per-read banded-DP that scales
with read length, and a tomato-scale CRAM has millions of reads).
Three options for the v1 design:

- **(a) Keep both layers as-is.** Rayon's work-stealing handles
  nested parallelism correctly in principle — workers from the inner
  par_iter steal from the outer pool when their chunks finish. The
  cost is per-task scheduling overhead at both layers. Simplest,
  ships fastest.
- **(b) Disable BAQ chunk parallelism inside per-chrom workers.**
  Each worker BAQ-caps its own chunks serially. Removes scheduler
  pressure; loses the existing BAQ parallelism for workloads with
  few chromosomes (small viruses, plasmid-only CRAMs) where the
  outer per-chrom layer has nothing to scale on.
- **(c) Disable BAQ chunk parallelism conditionally** — outer rayon
  parallelism present (`n_chroms >= 2`) → inner BAQ serial; outer
  absent (`n_chroms == 1`) → inner BAQ parallel. Best of both
  worlds; tiny added complexity in `BaqStream`'s configuration.

**Decision: (a) for v1, measure, then revisit.** The cohort H1
shipped without changing the inner BAQ parallelism for an equivalent
reason — and the analogous case there (inner per-group `par_iter`)
turned out to be a net loss only after measurement. The same
discipline applies: do not preemptively redesign scheduler shape;
ship the simple thing, measure on the real tomato fixture, and only
revisit (c) if the profile shows the inner par_iter as net-negative
under parallel outer workers.

### 7. Concat — reuse `src/vcf/concat.rs`

No new module. Each per-chrom worker produces a complete
self-contained `.vcf[.gz]` fragment with the same cohort header
(`metadata` is shared verbatim); the existing
[`concat_fragments`](../../src/vcf/concat.rs)
strips the duplicated header from fragments 2..N and copies
the bodies in order. Same atomic-rename + parent-dir-fsync
discipline as the cohort path.

This is the single biggest engineering reuse vs. Path A — the PSP
concat module that Path A would have needed (~300 LOC including
block-index rewriting + XXH3 recomputation) is replaced by zero
new code on this path.

### 8. Determinism

Outputs must be byte-identical between threads=1 and threads=N
(modulo the unused-allocator-side variance in BAQ chunk ordering,
which already exists today). Two failure modes to guard against:

- **Per-contig fragment ordering** — handled by concat walking the
  contig table in declared order. Same invariant as the cohort path.
- **PerPositionMerger insertion order at the same (chrom, pos)** —
  with a single sample, the per-position merger sees at most one
  entry per (chrom, pos), so the cross-sample tie-break the cohort
  side worries about does not apply here. Determinism falls out of
  the per-contig-fragment design.

Add an integration test that runs the same fixture at threads=1 and
threads=4 (or whatever the test host can muster) and diffs the
output VCF bodies.

## Test plan

Integration tests in
[tests/cohort_cli_integration.rs](../../tests/cohort_cli_integration.rs):

1. **`from_bam_parallel_emits_same_records_as_serial`** —
   golden-comparison test at threads=1 vs threads=4. Reuse the
   existing from-bam test fixture; assert record-level byte equality
   of the output VCF bodies (ignore header timestamps + `##source`).
2. **`from_bam_single_chromosome_input_works`** — degenerate case:
   single-contig fixture, one worker, concat with one fragment ≡
   identity.
3. **`from_bam_empty_chromosome_in_middle`** — three-contig fixture
   where the middle contig has zero reads (worker produces a
   header-only fragment).
4. **`from_bam_missing_index_without_flag_errors`** — fixture with a
   `.cram` but no `.crai`; assert `MissingMapFileIndex` error and
   that the error message contains both the flag name and the
   `samtools index` recipe.
5. **`from_bam_missing_index_with_flag_builds_successfully`** —
   same fixture, run with `--build-map-file-index`; assert the
   `.crai` is created next to the source CRAM and that the run
   completes with the same output as a pre-indexed run.
6. **`from_bam_preflight_catches_missing_index_before_parallel_work`**
   — fixture with one indexed and one un-indexed CRAM; assert the
   error fires before any fragment file is created in the tempdir
   (the tempdir should not even be created — pre-flight runs first).
7. **`from_bam_index_build_in_readonly_dir_errors_clearly`** — chmod
   the parent dir to read-only, run with `--build-map-file-index`,
   assert `IndexBuildFailed` with an actionable error message naming
   the path.

Unit tests in `src/bam/index_preflight.rs`:

8. **`preflight_accepts_existing_index`** — `.crai` present, no flag,
   `preflight_alignment_indexes` returns `Ok(())`.
9. **`preflight_builds_when_flag_set_and_missing`** — synthesize a
   tiny CRAM in a tempdir, run pre-flight with the flag, assert the
   `.crai` exists and is non-empty.
10. **`preflight_errors_when_missing_and_flag_unset`** — assert the
    error variant; the bam-module variant carries the path; the CLI-
    layer test (4 above) covers the flag-name wording.
11. **`preflight_idempotent_when_index_exists_and_flag_set`** —
    pre-existing `.crai`, flag set, assert no rebuild (mtime of the
    `.crai` unchanged within the test's timing resolution).

Unit tests in `src/bam/cram_input.rs` for the new query variant:

12. **`query_yields_only_records_for_requested_contig`** — multi-contig
    synthetic CRAM, query one contig, assert no off-target records.
13. **`query_preserves_coordinate_order_within_contig`** —
    same fixture, assert monotonic positions in the output.
14. **`query_k_way_merge_across_inputs_matches_streaming_merge`** —
    build two CRAMs covering the same contig, query and stream both,
    assert the record sequences match.
15. **`query_empty_contig_yields_no_records`** — query a contig with
    no aligned reads, assert the iterator is empty and no error.

Bench / measurement:

- New `examples/profile_from_bam_e2e.rs` mirror of
  [`examples/profile_cohort_e2e.rs`](../../examples/profile_cohort_e2e.rs),
  running the same `RAYON_NUM_THREADS=N` sweep at T=1, 2, 4, 8,
  13, 16 against a real tomato single-sample CRAM. **Acceptance
  threshold: T=8 wall ≥ 4× reduction vs. T=1**, mirroring the
  cohort H1 plan's threshold. Below 2× → the design needs
  re-thinking before merging.
- The standing "CRAM → VCF integration perf bench" from
  PROJECT_STATUS line 802 — landed as part of this slice, mirror of
  [`benches/cohort_e2e_perf.rs`](../../benches/cohort_e2e_perf.rs).
  Two bench families: a core variant that times the parallel
  driver in isolation and a full variant that times
  `run_var_calling_from_bam` end-to-end.

## Sequencing

A single PR. Split into reviewable commits in this order:

1. **Pre-flight + flag + errors.** New flag on
   `VarCallingFromBamArgs`, `index_preflight` module under
   `src/bam/`, `AlignmentIndexError` in `src/bam/errors.rs`, two new
   CLI error variants on `VarCallingFromBamCliError` with the
   `From<AlignmentIndexError>` impl, unit tests 8–11 + integration
   tests 4–7. Lands without touching the parallel path yet — the
   rest of `run_var_calling_from_bam` keeps its current serial body,
   but now requires indexes (or builds them) up front. This is a
   user-visible change in isolation; ships with its own release-note
   bullet.
2. **`CramMergedReader::query` + tests 12–15.** Standalone, no
   callers yet. Pure addition; no existing call path changes.
3. **`process_one_chromosome_from_bam` helper.** Mirrors the cohort
   `process_one_chromosome` shape; reuses `drive_cohort_pipeline`
   from `cohort_driver.rs`. No driver reshape yet.
4. **`run_var_calling_from_bam` reshape + integration tests 1–3.**
   Replaces the serial body with the par_iter + concat. This is
   the slice that delivers the wall-time win.
5. **Bench validation.** Add the `examples/profile_from_bam_e2e.rs`
   driver + the `benches/from_bam_e2e_perf.rs` criterion bench;
   run the T sweep on real tomato CRAMs; record the numbers in the
   impl report.

## Open work / non-goals

1. **`CohortPipelineParams` cloneability.** The cohort H1 added
   `#[derive(Clone)]` here for the same reason; the from-bam path
   uses the same struct via `drive_cohort_pipeline`, so no action is
   needed if Clone is already derived; verify in commit 3.
2. **`SyncRefFetcher` contention (L5).** Same as on the cohort H1 —
   the from-bam workers each build their own
   `MultiChromStreamingRefFetcher` from `args.reference`, so each
   per-chrom worker has its own per-contig sliding buffer with no
   cross-worker sharing. No `Arc<Mutex<…>>` contention by
   construction. If a future shared global FASTA cache is added,
   re-profile.
3. **Per-CRAM index sharing across runs.** Once an index is built
   it sits next to the source forever; subsequent runs see it and
   skip the build. No invalidation logic — stale-index handling is
   the user's problem (same as samtools).
4. **`.bai` vs `.csi` for BAM input.** v1 builds CSI by default
   (universal, contig-size-safe). Detect-and-trust either when
   pre-existing. If a real user request for BAI-preferred output
   appears, add a `--build-bai` opt-in; not in v1.
5. **Index-build cancellation.** If the user ctrl-Cs mid-build,
   noodles' index builders write to a tempfile and rename — partial
   writes won't poison the next run. Verify this assumption during
   implementation (one-line test: kill the process mid-build, run
   again with `--build-map-file-index`, assert success).
6. **`estimate-contamination`-from-bam.** No such subcommand exists
   today, and contamination estimation runs over `.psp` inputs only.
   Out of scope; track as a separate slice if a from-bam contamination
   path is ever wanted.
7. **PSP path (Path A).** Explicitly deferred per the design
   discussion — multi-sample PSP workflows already get parallelism
   from running multiple `pileup` processes. Revisit only if a
   real workload appears where a single sample's `.psp` build
   dominates a pipeline's wall time.

## Forward references

- Once this slice lands, the next ceiling is the same one the cohort
  H1 hit: ch00 / decoy-contig read-imbalance. Sub-chromosome
  decomposition (split a giant contig across multiple workers) is
  the order-of-magnitude lever beyond it. Track as a separate slice
  if profiling on multiple species shows it matters.
- The standing
  [parallel-optimisation integration perf benches](../../PROJECT_STATUS.md)
  open item for the CRAM → VCF arm closes with this slice.

## Estimated effort

- Pre-flight + flag + two error variants (bam) + two CLI error
  variants + unit tests 8–11 + integration tests 4–7: ~280 lines.
- `CramMergedReader::query` + tests 12–15: ~300 lines (the merge
  logic is shared with `new` via a small refactor; mostly new
  query-iterator plumbing).
- `process_one_chromosome_from_bam` helper: ~100 lines net diff
  (most of the body is lifted from `run_cohort_pipeline_for_single_sample`).
- `run_var_calling_from_bam` reshape + integration tests 1–3:
  ~150 lines net diff.
- `examples/profile_from_bam_e2e.rs` + `benches/from_bam_e2e_perf.rs`:
  ~250 lines.

Total: ~1080 net lines. Largest risk surfaces are (a) the
`CramMergedReader::query` k-way merge correctness against real
multi-CRAM samples (cover with unit tests 12–15 + the integration
golden-comparison test), and (b) the pre-flight error wording — get
the message right the first time; it's a user-facing API that's
expensive to change later.
