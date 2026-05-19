# `pop_var_caller` cohort CLI — Stages 3–6

Implementation plan for the second user-visible slice of the multi-sample
SNP caller: three new subcommands that take the pipeline from `.psp`
files (or a BAM, for one-off use) all the way to a multi-sample VCF.

Spec: [calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md)
§"Stage 3 — low-complexity filter" through §"Stage 6 — posterior engine".
Companion: [contamination_estimation.md](../specs/contamination_estimation.md).

This is glue, not new algorithms. Every stage already exists and is
shipped/fixes-applied:

- `.psp` reader: [src/per_sample_pileup/psp/reader.rs](../../src/per_sample_pileup/psp/reader.rs).
- Stage 3 DUST filter: [src/var_calling/dust_filter.rs](../../src/var_calling/dust_filter.rs).
- Stage 4 k-way merger + grouper: [src/var_calling/per_position_merger.rs](../../src/var_calling/per_position_merger.rs), [src/var_calling/variant_grouping.rs](../../src/var_calling/variant_grouping.rs).
- Stage 5 per-group merger: [src/var_calling/per_group_merger.rs](../../src/var_calling/per_group_merger.rs).
- Stage 6 posterior engine: [src/var_calling/posterior_engine.rs](../../src/var_calling/posterior_engine.rs).
- Stage 6 contamination side-pass: [src/var_calling/contamination_estimation.rs](../../src/var_calling/contamination_estimation.rs).
- Stage 6 VCF sink: [src/var_calling/vcf_writer/](../../src/var_calling/vcf_writer/).

The work here is the CLI surface, three orchestrators, the contamination
artefact format, batch-assignment I/O, the engine-side `Config::new ->
Result` validation gap (M3 / Mi12 from the Stage 6 review), and the
end-to-end integration tests that have been waiting on this slice.

## Three subcommands

| Subcommand | Inputs | Output | Contamination? |
|---|---|---|---|
| `var-calling` | one-or-more `.psp` | multi-sample VCF | optional via `--contamination-estimates FILE` |
| `estimate-contamination` | one-or-more `.psp` | contamination TOML | this *is* the contamination subcommand |
| `var-calling-from-bam` | one sample's BAM(s) | single-sample VCF | none (always `c_s = 0`) |

`var-calling` is the standard population-level path. `estimate-contamination`
produces the artefact that `var-calling` consumes. `var-calling-from-bam`
is a one-off convenience: pileup + var-calling fused, no `.psp` on disk.

## Scope

In scope:

- Three new clap-derive subcommands added to
  [PopVarCallerCommand](../../src/pop_var_caller/cli.rs#L55).
- Three new orchestrator modules (`var_calling.rs`,
  `estimate_contamination.rs`, `var_calling_from_bam.rs`) under
  [src/pop_var_caller/](../../src/pop_var_caller/), siblings of
  the existing [cli.rs](../../src/pop_var_caller/cli.rs).
- A shared **contamination artefact** module
  (`src/pop_var_caller/contamination_artifact.rs`): TOML schema,
  serde-derived structs, validator (no two batches with the same id,
  every sample's `batch` field resolves into `[[batches]]`, every
  probability is finite and in `[0, 1]`, q_b row sums to 1.0 within
  tolerance).
- A shared **batch-assignment** module
  (`src/pop_var_caller/batch_assignment.rs`): two-column TSV
  reader.
- **Engine-side validation gap closure** — `Config::new -> Result`
  for `PosteriorEngineConfig` (Mi12), `ContaminationEstimationConfig`
  (currently `with_project_defaults` only, no validation),
  `GrouperConfig`, `PerGroupMergerConfig`, `WriterConfig`.
  `DustFilterConfig::new` already returns `Result` ✓.
- **CLI-side validation** — every `f64`/`u32` knob with a meaningful
  range gets a `value_parser` (mirrors
  [parse_mismatch_fraction](../../src/pop_var_caller/cli.rs#L232)).
  Validation lives in **both** places per the
  [memory: keep settings comprehensive](../../../.claude/projects/-home-jose-devel-join-per-sample-vcfs/memory/feedback_keep_settings_complete.md)-style
  rule of thumb the project applies: defence in depth, library users
  get the same guard CLI users get.
- One small `pileup` refactor: extract the iterator-building part of
  [run_pileup](../../src/pop_var_caller/cli.rs#L276) into a reusable
  helper so `var-calling-from-bam` can drive the Stage 1 pipeline
  without writing a `.psp` (details below).
- Integration tests covering the three paths end-to-end. The
  fixture is a small reference FASTA + a synthetic CRAM (or a couple
  of pre-built `.psp` files) generated at test setup time — same
  approach as the pileup CLI integration test.
- DUST filter CLI bindings (open item from the Stage 3
  fixes-applied report): `--no-complexity-filter`,
  `--complexity-window`, `--complexity-threshold`.

Deferred to follow-up slices (intentionally not in this cut):

- `--regions BED` for restricting calls to a list of regions. Open
  standing item per [TODO.txt](../TODO.txt). The CLI struct stays
  append-only so this slot opens later without breaking the surface.
- `tabix .tbi` index generation alongside `.vcf.gz` outputs.
  Out-of-scope per the Stage 6 VCF-writer report.
- `PL` (phred-scaled likelihoods) `FORMAT` field. Needs Stage 5 →
  `PosteriorRecord` forwarding of `log_likelihoods`; threaded later.
- Per-sample contamination fraction in `INFO`. Wires in once this
  slice lands and is shaken out.
- `bcftools view` / `bcftools stats` manual smoke-test — runs once
  the cohort CLI exists, not part of this plan's deliverable.
- Multi-batch contamination defaults to "all samples in one batch
  named `all_samples`"; user-supplied batch assignments via
  `--batch-assignment` get honoured but no auto-batching heuristic
  (e.g. by sequencing lane) is included.
- Stage 1 BED-skip flag (separate standing item).
- A `Stage 5` implementation report (separate standing item — write
  alongside this slice if natural, otherwise as a follow-up).
- `--per-group-batch-size` exposed as a flag (deliberately kept
  library-default-only per
  [project memory](../../../.claude/projects/-home-jose-devel-join-per-sample-vcfs/memory/project_parallelization_tuning_sequencing.md);
  parallelisation knobs get a coordinated tuning pass after this
  slice lands).
- Cohort-wide rayon parallelism strategy beyond `--threads` sizing
  the global pool. Stage 5 already parallelises per-group via the
  existing rayon usage; this slice does not add a new layer.

Not planned at all:

- **Multi-sample BAM input to `var-calling-from-bam`.** That
  subcommand takes **one sample**'s BAMs (just like `pileup`). For
  multi-sample work the user goes through `.psp`.
- **In-line contamination estimation inside `var-calling`.** The
  user runs `estimate-contamination` first; `var-calling` only ever
  *reads* an estimates artefact.
- **A `--no-contamination` flag.** Absence of
  `--contamination-estimates` already means `c_s = 0`.
- **A `--filter-pass` flag.** v1 always emits `PASS` in FILTER;
  reserved for a future filter slice.
- **External allele-frequency panels** (the
  `--external-allele-frequencies` slot mentioned in the
  contamination spec). v2 follow-up.
- **Auto-batching of samples** (e.g. by sequencing run / lane / flow
  cell from BAM `@RG` tags). User-supplied `--batch-assignment` is
  the only way to define batches.

## Module layout

```
src/pop_var_caller/
  mod.rs                       — modified: declare new modules
  cli.rs                       — modified: add three subcommand variants
  cli/error_bridge.rs          — unchanged
  psp_to_pileup.rs             — unchanged
  var_calling.rs               — NEW: VarCallingArgs, run_var_calling,
                                       VarCallingCliError
  estimate_contamination.rs    — NEW: EstimateContaminationArgs,
                                       run_estimate_contamination,
                                       EstimateContaminationCliError
  var_calling_from_bam.rs      — NEW: VarCallingFromBamArgs,
                                       run_var_calling_from_bam,
                                       VarCallingFromBamCliError
  contamination_artifact.rs    — NEW: TOML schema (serde), reader,
                                       writer, validator
  batch_assignment.rs          — NEW: two-column TSV reader
  stage1_pipeline.rs           — NEW: extracted from run_pileup;
                                       returns a PileupRecord iterator
                                       without a writer
```

Conservation principle: each new file owns one concern. Future
subcommands slot in as further siblings; nothing in
[var_calling/](../../src/var_calling/) (the algorithm-side library
slice) is touched.

## Shared infrastructure

### Contamination artefact — TOML schema

Produced by `estimate-contamination`, consumed by `var-calling`. Single
file, normalised (no redundancy: a sample carries only its `batch`
label and `contamination_fraction`; the per-batch `q_b` lives in
`[[batches]]` keyed by `id`).

```toml
[provenance]
tool         = "pop_var_caller"
version      = "0.1.0"
subcommand   = "estimate-contamination"
created      = 2026-05-19T14:32:11Z

[provenance.inputs]
reference         = "grch38.fa"               # basename only
input_psps        = ["NA12878.psp", "NA12891.psp"]  # basenames
batch_assignment  = "lanes.tsv"               # basename or omitted

[parameters]
# Every effective CLI knob — recorded for reproducibility. Mirrors
# `WriterProvenance.parameters` in the .psp header.
stopping_mode             = "convergence"
block_size                = 1000
min_depth                 = 10
min_major_fraction        = 0.95
min_cohort_minor_count    = 2
min_cohort_minor_fraction = 0.005
min_batch_size            = 5
ref_pseudocount           = 10.0
snp_alt_pseudocount       = 0.01
indel_alt_pseudocount     = 0.00125
c_s_init                  = 0.02
q_b_init_per_class        = 0.25
stability_tolerance       = 0.001
stability_blocks          = 3
threads                   = 8

[[batches]]
id                            = "lane_3"
contaminant_ref_prob          = 0.9989
contaminant_snp_alt_prob      = 0.0010
contaminant_indel_alt_prob    = 0.0001

[[batches]]
id                            = "lane_7"
contaminant_ref_prob          = 0.9991
contaminant_snp_alt_prob      = 0.0008
contaminant_indel_alt_prob    = 0.0001

[[samples]]
name                  = "NA12878"
batch                 = "lane_3"
contamination_fraction = 0.0123

[[samples]]
name                  = "NA12891"
batch                 = "lane_3"
contamination_fraction = 0.0031
```

Validator (runs at load time in `var-calling` and at write time in
`estimate-contamination`):

- Every `[[batches]].id` is unique.
- Every `[[samples]].batch` resolves to an existing
  `[[batches]].id` — unresolved batch reference is a hard error.
- Every probability is finite and in `[0.0, 1.0]`.
- Every batch's three `contaminant_*_prob` fields sum to `1.0`
  within `1e-9` (compound is implicitly 0 per Assumption 2 of the
  contamination impl report).
- Every `contamination_fraction` is finite and in `[0.0, 1.0]`.
- Sample names are unique within `[[samples]]`.

Sample-name reconciliation against the `.psp` inputs (done by
`var-calling` after loading the artefact):

- Every `.psp` sample name must appear in `[[samples]]`. A `.psp`
  sample missing from the artefact is a **hard error**:
  `VarCallingCliError::ContaminationSampleMissing { sample }`.
- Extras in `[[samples]]` (entries that don't correspond to any
  `.psp` input) are **silently ignored**. Rationale: users may run
  contamination over a superset cohort and then call on a subset.

Public Rust API:

```rust
// src/pop_var_caller/contamination_artifact.rs

#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct ContaminationArtefact {
    pub provenance: Provenance,
    pub parameters: BTreeMap<String, toml::Value>,
    pub batches: Vec<BatchEntry>,
    pub samples: Vec<SampleEntry>,
}

pub struct BatchEntry { id: String, contaminant_ref_prob: f64,
                        contaminant_snp_alt_prob: f64,
                        contaminant_indel_alt_prob: f64 }
pub struct SampleEntry { name: String, batch: String,
                         contamination_fraction: f64 }

impl ContaminationArtefact {
    pub fn read(path: &Path) -> Result<Self, ContaminationArtefactError>;
    pub fn write(&self, path: &Path) -> Result<(), ContaminationArtefactError>;
    /// Returns a per-sample map suitable for handing to the posterior engine.
    pub fn to_estimates_for_samples(
        &self, sample_names: &[&str]
    ) -> Result<ContaminationEstimates, ContaminationArtefactError>;
}
```

`ContaminationEstimates` is the existing
[engine-side type](../../src/var_calling/contamination_estimation.rs).
The artefact module is the boundary that converts on-disk TOML to
in-memory typed state.

### Batch assignment — TSV schema

Optional input to `estimate-contamination`. Plain TSV, header row.

```
sample	batch
NA12878	lane_3
NA12891	lane_3
NA12892	lane_7
```

- Two columns exactly: `sample`, `batch`.
- Header row required.
- Tab-separated; no quoting (sample/batch ids may not contain tab or
  newline; this is a hard error at parse time).
- Unknown extra rows (samples not present in the `.psp` inputs) are
  silently ignored. A `.psp` sample missing from the TSV gets
  assigned the default batch id `all_samples` (consistent with the
  no-TSV case).
- A line whose `batch` field is empty is a hard error.

```rust
// src/pop_var_caller/batch_assignment.rs

pub struct BatchAssignment {
    by_sample: HashMap<String, String>,
}

impl BatchAssignment {
    pub fn from_tsv(path: &Path) -> Result<Self, BatchAssignmentError>;
    pub fn empty() -> Self;
    /// Look up a sample's batch; falls back to `"all_samples"`.
    pub fn batch_for(&self, sample: &str) -> &str;
}
```

### Config validation — engine + CLI

The Stage 6 review left M3/Mi12 open: `PosteriorEngineConfig` has
no `Config::new -> Result`. The same pattern is missing for
`ContaminationEstimationConfig`, `GrouperConfig`,
`PerGroupMergerConfig`, and `WriterConfig`. This slice closes all
of them.

Pattern (matches the
[`DustFilterConfig::new`](../../src/var_calling/dust_filter.rs#L501)
shape — that one is already correct):

```rust
impl PosteriorEngineConfig {
    pub fn new(
        convergence_threshold: f64,
        max_iterations: u32,
        ref_pseudocount: f64,
        snp_alt_pseudocount: f64,
        indel_alt_pseudocount: f64,
        compound_alt_pseudocount: f64,
        inbreeding_coefficient: f64,
        max_gq_phred: f64,
    ) -> Result<Self, PosteriorEngineConfigError> {
        // every numeric field validated for finiteness;
        // pseudocounts must be > 0;
        // inbreeding_coefficient in [0.0, 1.0];
        // convergence_threshold in (0.0, 1.0);
        // max_iterations > 0;
        // max_gq_phred in (0.0, 1000.0).
    }
}
```

Validation ranges (used identically in CLI `value_parser` and engine
`new`):

| Field | Range |
|---|---|
| pseudocounts (ref/snp/indel/compound) | finite, `(0.0, 1000.0]` |
| inbreeding_coefficient | finite, `[0.0, 1.0]` |
| convergence_threshold | finite, `(0.0, 0.1]` |
| max_iterations | `1..=500` |
| max_gq_phred | finite, `(10.0, 200.0]` |
| c_s_init | finite, `[0.0, 0.5]` |
| q_b_init_per_class | finite, `[0.0, 1.0]` |
| min_major_fraction | finite, `(0.5, 1.0]` |
| min_cohort_minor_fraction | finite, `[0.0, 0.5)` |
| stability_tolerance | finite, `(0.0, 0.1]` |
| max_alleles_per_var | `2..=16` |
| ploidy | `1..=8` (`2` is the routine path) |

Several ranges above are deliberately **tighter than the mathematical
bound**. The principle: prevent absurd choices when we can recognise
them, rather than accept any value the math doesn't outright reject.
Rationales per knob:

- `convergence_threshold`, `stability_tolerance`: `(0.0, 0.1]`. A
  threshold of `1.0` would terminate EM after one iteration (max
  possible per-allele change is `1.0`). The `0.1` upper bound is
  "10 % per-allele change still counts as converged" — already very
  loose; the default `1e-4` is three orders of magnitude tighter.
- `max_iterations`: `1..=500`. EM converges in 3–5 iterations on the
  GATK reference; the default is `50`. `500` is 10× the default —
  generous headroom for debugging and convergence checks, well below
  the runaway zone.
- `max_gq_phred`: `(10.0, 200.0]`. Phred `200` corresponds to
  `p = 1e-20` — already absurdly confident; Phred `1000` has no
  physical meaning.
- `c_s_init`: `[0.0, 0.5]`. A contamination prior `> 0.5` flips the
  model — the sample's "real" reads become the assumed contaminant.
  Domain-incoherent.
- `min_cohort_minor_fraction`: `[0.0, 0.5)`. A "minor"-allele fraction
  `≥ 0.5` is a contradiction in terms.
- `max_alleles_per_var`: `2..=16`. `1` means no alts ever (REF-only,
  degenerate); `16` is wide enough for the most polyallelic real loci
  (the GATK default is `6`).
- `ploidy`: `1..=8`. Octoploid wheat is the high-water mark for real
  cohorts; `255` is meaningless.
- pseudocounts: `(0.0, 1000.0]`. A pseudocount `> 1000` swamps any
  realistic cohort-level evidence and effectively forces the prior
  unconditionally.

The CLI uses `value_parser = parse_xxx` per knob, with one tiny
helper per range pattern (`parse_finite_unit_interval`,
`parse_finite_positive`, etc.). The same range checks live in the
engine `Config::new`, returning a typed error. Defence in depth:
library users see the same guard CLI users do.

### Error types

One top-level error per subcommand, structured like the existing
[PileupCliError](../../src/pop_var_caller/cli.rs#L247). Each wraps
the underlying module errors via `#[from]`.

```rust
#[derive(thiserror::Error, Debug)]
pub enum VarCallingCliError {
    #[error("psp reader: {0}")]      PspReader(#[from] PspReadError),
    #[error("dust filter: {0}")]     Dust(#[from] DustFilterError),
    #[error("merger: {0}")]          Merger(#[from] PerPositionMergerError),
    #[error("grouper: {0}")]         Grouper(#[from] GrouperError),
    #[error("per-group merger: {0}")] PerGroup(#[from] PerGroupMergerError),
    #[error("posterior engine: {0}")] Posterior(#[from] PosteriorEngineError),
    #[error("vcf writer: {0}")]      Vcf(#[from] VcfWriteError),
    #[error("contamination artefact: {0}")]
    ContamArtefact(#[from] ContaminationArtefactError),
    #[error("sample '{sample}' present in .psp inputs but missing from contamination estimates")]
    ContaminationSampleMissing { sample: String },
    #[error("config: {0}")]          Config(#[from] CohortConfigError),
    #[error("rayon thread pool already initialised — refusing to override")]
    RayonAlreadyConfigured,
    #[error("io: {0}")]              Io(#[from] io::Error),
}
```

`CohortConfigError` aggregates the per-stage `Config::new` errors
behind one variant (`#[from]` chains). The `EstimateContamination`
and `VarCallingFromBam` errors follow the same shape — variants
that don't apply (e.g. `Dust`, `ContamArtefact` in the from-bam
path) are omitted.

## Subcommand: `var-calling`

### CLI surface

Tiered exactly like
[`PileupArgs`](../../src/pop_var_caller/cli.rs#L68): a handful of
common flags visible in `-h`, the rest behind `hide_short_help =
true` and grouped under `help_heading` sections.

```rust
#[derive(clap::Args)]
pub struct VarCallingArgs {
    // ===== Common =============================================
    #[arg(long)] reference: PathBuf,
    #[arg(long)] output: PathBuf,
    #[arg(long)] threads: Option<usize>,
    #[arg(long)] contamination_estimates: Option<PathBuf>,
    #[arg(required = true)] psp_files: Vec<PathBuf>,

    /// Skip the low-complexity (sdust) filter entirely.
    #[arg(long)] no_complexity_filter: bool,

    /// Cohort-wide ploidy.
    #[arg(long, default_value_t = DEFAULT_PLOIDY,
          value_parser = parse_ploidy)]
    ploidy: u8,

    // ===== Advanced — DUST filter ============================
    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_DUST_WINDOW,
          value_parser = parse_dust_window,
          help_heading = "Advanced — DUST filter")]
    complexity_window: u32,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_DUST_THRESHOLD,
          help_heading = "Advanced — DUST filter")]
    complexity_threshold: u32,

    // ===== Advanced — grouper / merger ======================
    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MAX_VARIANT_GROUP_SPAN,
          help_heading = "Advanced — Variant grouping")]
    var_group_max_span: u32,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MAX_ALLELES_PER_RECORD,
          value_parser = parse_max_alleles,
          help_heading = "Advanced — Per-group merger")]
    max_alleles_per_var: usize,

    // ===== Advanced — posterior engine ======================
    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_INBREEDING_COEFFICIENT,
          value_parser = parse_finite_unit_interval,
          help_heading = "Advanced — Posterior engine")]
    inbreeding_coefficient: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_CONVERGENCE_THRESHOLD,
          value_parser = parse_finite_open_unit_interval,
          help_heading = "Advanced — Posterior engine")]
    em_convergence_threshold: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MAX_ITERATIONS,
          help_heading = "Advanced — Posterior engine")]
    em_max_iterations: u32,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_REF_PSEUDOCOUNT,
          value_parser = parse_finite_positive,
          help_heading = "Advanced — Posterior engine")]
    ref_pseudocount: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_SNP_ALT_PSEUDOCOUNT,
          value_parser = parse_finite_positive,
          help_heading = "Advanced — Posterior engine")]
    snp_alt_pseudocount: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_INDEL_ALT_PSEUDOCOUNT,
          value_parser = parse_finite_positive,
          help_heading = "Advanced — Posterior engine")]
    indel_alt_pseudocount: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_COMPOUND_ALT_PSEUDOCOUNT,
          value_parser = parse_finite_positive,
          help_heading = "Advanced — Posterior engine")]
    compound_alt_pseudocount: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = MAX_GQ_PHRED,
          value_parser = parse_max_gq_phred,
          help_heading = "Advanced — Posterior engine")]
    max_gq_phred: f64,

    // ===== Advanced — VCF writer ============================
    #[arg(long, hide_short_help = true,
          help_heading = "Advanced — VCF writer")]
    emit_gp: bool,
}
```

Common flags (in `-h`): `--reference`, `--output`, `--threads`,
`--contamination-estimates`, `--no-complexity-filter`, `--ploidy`,
plus the positional `<psp_files…>`. Everything else is behind
`--help` under three groups: "Advanced — DUST filter",
"Advanced — Variant grouping" / "Per-group merger",
"Advanced — Posterior engine", "Advanced — VCF writer".

### Orchestrator: `run_var_calling`

```rust
pub fn run_var_calling(args: &VarCallingArgs) -> Result<(), VarCallingCliError>;
```

1. Size rayon's global pool from `--threads` (same pattern as
   `run_pileup`).
2. Open every `.psp` with `PspReader` (each opens its own file
   handle; the readers are lazy iterators).
3. Cross-check: every reader's `header.reference` basename equals
   the `--reference` basename. Otherwise error.
   The reference fasta + .fai are opened via the existing
   `SyncRefFetcher`.
4. Build engine configs from CLI args via the validated `new` /
   `with_config` constructors. Any constructor failure surfaces as
   `CohortConfigError`.
5. **Load contamination estimates if supplied.** Read the TOML
   artefact via `ContaminationArtefact::read(path)`. Validate
   sample-name reconciliation: every `.psp` sample must appear in
   `[[samples]]`. Then call `.to_estimates_for_samples(&psp_sample_names)`
   to build the in-memory `ContaminationEstimates` typed value. If
   `--contamination-estimates` is absent, use the engine's "no
   contamination" sentinel (a `ContaminationEstimates::zero()`
   helper — added if it doesn't exist).
6. Wire the pipeline:
   - `Vec<PspReader>` → `PerPositionMerger::new(readers, ...)`.
   - If `!--no-complexity-filter`:
     `DustFilter::new(merger, ref_fetcher, dust_cfg)`.
   - `VariantGrouper::new(filter_or_merger, grouper_cfg)`.
   - `PerGroupMerger::new(grouper, per_group_cfg)`.
   - `PosteriorEngine::new(per_group, posterior_cfg, contamination_estimates)`.
   - `CohortVcfWriter::new(<output>.tmp, header, writer_cfg)`.
7. Drive the pipeline: `for record in posterior_engine { writer.write_record(record)?; }`
   (or whatever the existing seam is).
8. Finalise the writer: flush, fsync, atomic rename
   `<output>.tmp` → `<output>` — same write-tmp-then-rename
   convention as `run_pileup` for the same reasons (partial files
   are silently parseable; we adopt the convention even though
   most bioinformatics CLIs don't).
9. Stderr run-summary: per-stage counts (records filtered by
   DUST, records emitted by grouper, records emitted by writer,
   wall time per stage if cheap).

The structure mirrors `run_pileup` step-for-step; the difference
is which stages are wired and the contamination side-channel.

## Subcommand: `estimate-contamination`

### CLI surface

```rust
#[derive(clap::Args)]
pub struct EstimateContaminationArgs {
    // ===== Common =============================================
    #[arg(long)] reference: PathBuf,
    #[arg(long)] output: PathBuf,           // path for the TOML artefact
    #[arg(long)] threads: Option<usize>,
    #[arg(long)] batch_assignment: Option<PathBuf>,
    #[arg(required = true)] psp_files: Vec<PathBuf>,

    // ===== Advanced — Stopping mode ==========================
    #[arg(long, hide_short_help = true,
          default_value_t = StoppingMode::Convergence,
          help_heading = "Advanced — Stopping mode")]
    stopping_mode: StoppingMode,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_STABILITY_TOLERANCE,
          value_parser = parse_finite_open_unit_interval,
          help_heading = "Advanced — Stopping mode")]
    stability_tolerance: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_STABILITY_BLOCKS,
          help_heading = "Advanced — Stopping mode")]
    stability_blocks: u32,

    // ===== Advanced — Online EM ==============================
    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_BLOCK_SIZE,
          help_heading = "Advanced — Online EM")]
    block_size: u32,

    // ===== Advanced — Informative-site cuts ==================
    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MIN_DEPTH,
          help_heading = "Advanced — Informative-site cuts")]
    min_depth: u32,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MIN_MAJOR_FRACTION,
          value_parser = parse_min_major_fraction,
          help_heading = "Advanced — Informative-site cuts")]
    min_major_fraction: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MIN_COHORT_MINOR_COUNT,
          help_heading = "Advanced — Informative-site cuts")]
    min_cohort_minor_count: u32,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MIN_COHORT_MINOR_FRACTION,
          value_parser = parse_finite_unit_interval,
          help_heading = "Advanced — Informative-site cuts")]
    min_cohort_minor_fraction: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MIN_BATCH_SIZE_FOR_CONTAMINATION,
          help_heading = "Advanced — Informative-site cuts")]
    min_batch_size: u32,

    // ===== Advanced — Priors =================================
    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_REF_PSEUDOCOUNT,
          value_parser = parse_finite_positive,
          help_heading = "Advanced — Priors")]
    ref_pseudocount: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_SNP_ALT_PSEUDOCOUNT,
          value_parser = parse_finite_positive,
          help_heading = "Advanced — Priors")]
    snp_alt_pseudocount: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_INDEL_ALT_PSEUDOCOUNT,
          value_parser = parse_finite_positive,
          help_heading = "Advanced — Priors")]
    indel_alt_pseudocount: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_C_S_INIT,
          value_parser = parse_finite_unit_interval,
          help_heading = "Advanced — Priors")]
    c_s_init: f64,

    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_Q_B_INIT_PER_CLASS,
          value_parser = parse_finite_unit_interval,
          help_heading = "Advanced — Priors")]
    q_b_init_per_class: f64,
}
```

### Orchestrator: `run_estimate_contamination`

```rust
pub fn run_estimate_contamination(
    args: &EstimateContaminationArgs,
) -> Result<(), EstimateContaminationCliError>;
```

1. Size rayon pool from `--threads`.
2. Open every `.psp` with `PspReader`.
3. Cross-check references (same rule as `var-calling`).
4. Load batch assignment if supplied. Default: every sample
   maps to batch `"all_samples"`.
5. Build `ContaminationEstimationConfig` via its validated `new`
   constructor.
6. Call the existing
   [`estimate_contamination`](../../src/var_calling/contamination_estimation.rs)
   library API with the readers, batch map, and config. This is
   a single-pass library call — the orchestrator just hands the
   inputs in.
7. Convert the engine-side `ContaminationEstimates` result back
   into a `ContaminationArtefact` (provenance, parameters,
   batches, samples).
8. Write the artefact to `<output>.tmp` then rename to
   `<output>`. Same write-tmp-then-rename discipline.

Stderr run-summary: per-batch convergence info (blocks consumed,
final stability metric), per-batch sample count, and the number
of samples that fell below `--min-batch-size` and received
`c_s = 0`.

## Subcommand: `var-calling-from-bam`

A convenience subcommand: pileup + var-calling fused, no `.psp` on
disk, no contamination correction. One sample per invocation. Useful
for one-off calls, ad-hoc testing, and the integration-test smoke
path.

### Required refactor: `stage1_pipeline.rs`

The existing `run_pileup` builds the Stage 1 pipeline (CRAM →
BAQ → walker) and drives it into a `.psp` writer. For
`var-calling-from-bam` we need the same pipeline driving into the
Stage 4 merger instead. Extract the pipeline-building part:

```rust
// src/pop_var_caller/stage1_pipeline.rs

/// Build the Stage 1 pipeline (CRAM merged-reader → BAQ → walker)
/// and return its output as an iterator of `PileupRecord` plus
/// the runtime handles needed for end-of-run statistics.
pub struct Stage1Pipeline<'r> {
    pub records: Box<dyn Iterator<Item = PileupRecord> + 'r>,
    pub error_handle: ErrorHandle,
    pub filter_counts: &'r RefCell<FilterCounts>,
    pub baq_skip_counts: Option<&'r RefCell<BaqSkipCounts>>,
}

pub fn build_stage1_pipeline<'r>(
    crams: &[PathBuf],
    reference: &Path,
    cram_cfg: CramMergedReaderConfig,
    baq_cfg: BaqConfig,
    walker_cfg: WalkerConfig,
    baq_chunk_size: usize,
    no_baq: bool,
    /* + lifetimes / arena */
) -> Result<Stage1Pipeline<'r>, PileupCliError>;
```

`run_pileup` is reshaped to call `build_stage1_pipeline(...)` and
feed `records` into the writer.
`run_var_calling_from_bam` calls the same helper and feeds
`records` into a single-sample `PerPositionMerger` (k=1) followed
by the rest of the var-calling pipeline. **No `.psp` on disk
anywhere in this path.**

Lifetime detail: the merged reader and BAQ stream are owned by
the pipeline struct; the iterator borrows from them. The
borrow-juggling pattern in
[`run_pileup`](../../src/pop_var_caller/cli.rs#L329) (the
"borrow order" comment at lines 320-327) needs to translate into
the helper without escaping `&mut` references. A practical
shape: the helper takes a stack-allocated arena built by the
caller (one `CramMergedReader`, one optional `BaqStream`) and
returns an iterator that borrows from it. The arena lives on
the caller's stack for the full pipeline run.

### CLI surface

Composition of `PileupArgs` (the Stage 1 knobs) and the
var-calling knobs, minus contamination, minus `--contamination-estimates`,
plus a single output `.vcf` / `.vcf.gz` instead of a `.psp`.

```rust
#[derive(clap::Args)]
pub struct VarCallingFromBamArgs {
    // ===== Stage 1 knobs (mirror PileupArgs) ==================
    #[arg(long)] reference: PathBuf,
    #[arg(long)] output: PathBuf,
    #[arg(long)] threads: Option<usize>,
    #[arg(required = true)] crams: Vec<PathBuf>,

    #[arg(long, default_value_t = DEFAULT_MIN_MAPQ)]
    min_mapq: u8,
    #[arg(long)] no_baq: bool,

    // (… all the advanced Stage 1 knobs from PileupArgs reused
    //  verbatim — same defaults, same value_parsers, same
    //  help_heading sections; consider sharing the struct via
    //  `#[command(flatten)]` if that works with clap-derive)

    // ===== Stage 3-6 knobs (mirror VarCallingArgs minus contam) ==
    #[arg(long)] no_complexity_filter: bool,
    #[arg(long, default_value_t = DEFAULT_PLOIDY,
          value_parser = parse_ploidy)]
    ploidy: u8,
    // … all the advanced Stage 3–6 knobs from VarCallingArgs ………
}
```

To avoid copy-pasting ~30 fields, factor the shared Stage 1 and
Stage 3–6 advanced-flag groups into `#[derive(Args)]` structs and
`#[command(flatten)]` them into both `VarCallingFromBamArgs` and
the per-subcommand args. clap-derive supports this and the result
is a single source of truth per group.

### Orchestrator: `run_var_calling_from_bam`

```rust
pub fn run_var_calling_from_bam(
    args: &VarCallingFromBamArgs,
) -> Result<(), VarCallingFromBamCliError>;
```

1. Size rayon pool.
2. Open `CramMergedReader` and capture sample name, contigs.
3. Build Stage 1 pipeline via `build_stage1_pipeline(...)` — same
   helper `run_pileup` uses.
4. Wrap the `PileupRecord` iterator into a single-sample input
   for `PerPositionMerger` (k=1; the merger already supports
   this shape — verify in implementation).
5. From there, same wiring as `run_var_calling`: optional DUST,
   grouper, per-group merger, posterior engine
   (with `ContaminationEstimates::zero()`), VCF writer.
6. Finalise output with write-tmp-then-rename.

Stderr run-summary: merged Stage 1 counts (filter counts, baq
counts, walker counts — same as `run_pileup` prints) plus the
Stage 3–6 counts (DUST drops, records emitted by writer).

## Test strategy

### Unit tests (per new module)

- `contamination_artifact.rs`: TOML round-trip on a representative
  artefact; validator hard-errors on every illegal-state class
  (duplicate batch id, dangling sample.batch, q_b that doesn't
  sum to 1, NaN/∞ in any numeric field, contamination_fraction
  out of `[0, 1]`); sample-name reconciliation (missing → error,
  extras → ignored).
- `batch_assignment.rs`: TSV parse happy path; missing header
  error; empty batch field error; extras ignored;
  `batch_for(unknown_sample)` returns `"all_samples"`.
- Each `Config::new`: every documented range bound exercised on
  the boundary (just inside, exactly on, just outside).
- Each `parse_xxx` `value_parser`: same boundary set, plus
  rejection of `NaN`, `±∞`, non-numeric strings.

### Integration tests

`tests/cohort_cli_integration.rs` (one file, multiple
`#[test]` cases). Fixtures live under `tests/data/cohort/` —
small reference FASTA (a few hundred bp) plus three synthetic
CRAMs / `.psp` files generated at test-setup time, similar to
how
[`tests/pileup_cli_integration.rs`](../../tests/pileup_cli_integration.rs)
builds CRAMs from a SAM/FASTA pair.

Cases:

1. **`var-calling` happy path, no contamination.** Three
   pre-built `.psp` files, default config. Assert: VCF is
   produced, opens with `noodles::vcf` cleanly, has the right
   number of samples, ≥ 1 record. Provenance VCF header has
   `##source=pop_var_caller …` with the recorded subcommand.

2. **`var-calling` with contamination estimates.** Same `.psp`
   fixtures, plus a hand-crafted contamination TOML. Assert:
   the engine ran with the supplied `c_s` (visible via the
   `INFO` per-sample contamination once that field exists; for
   v1, assert via an inserted log-line or a small `diagnostics`
   accessor on the engine).

3. **`var-calling` rejects missing-sample artefact.**
   Contamination TOML missing one of the cohort samples; expect
   `VarCallingCliError::ContaminationSampleMissing { sample }`.

4. **`var-calling` accepts extras in artefact.** Contamination
   TOML has one sample the cohort doesn't; runs cleanly.

5. **`var-calling --no-complexity-filter`.** Same fixtures;
   assert: more records emitted than the with-DUST run (the
   fixture is built to include low-complexity regions that DUST
   would otherwise drop).

6. **`estimate-contamination` happy path.** Three `.psp` fixtures
   plus a `batches.tsv` putting two samples in batch `lane_3`
   and one in `lane_7`. Assert: the produced TOML round-trips
   through `ContaminationArtefact::read`, validates, contains
   exactly two `[[batches]]` entries and three `[[samples]]`
   entries.

7. **`estimate-contamination` default batch.** Same `.psp`
   fixtures, no `--batch-assignment`. Assert: one `[[batches]]`
   entry with `id = "all_samples"` and three sample entries
   pointing at it.

8. **End-to-end chain.** Run `estimate-contamination` to a temp
   TOML, then `var-calling --contamination-estimates <that
   TOML>`. Assert: the chain produces a valid VCF.

9. **`var-calling-from-bam` happy path.** One synthetic CRAM,
   small ref. Assert: VCF produced, single sample, ≥ 1 record;
   the per-record output matches what `var-calling` produces
   from the equivalent `.psp` (built via a separate
   `pop_var_caller pileup` call inside the test).

10. **Config validation rejects bad ranges.** Programmatically
    construct args with `inbreeding_coefficient = 1.5`; expect
    `CohortConfigError`. (Bypasses the CLI parser to exercise
    the engine-side `Config::new -> Result`.)

The integration tests do *not* re-cover correctness of each
stage — those are owned by per-module tests. They lock the glue.

### Manual smoke

Before declaring the slice done:

```
$ cargo run --release --bin pop_var_caller -- estimate-contamination \
    --reference ref.fa --output cohort.contam.toml \
    sample1.psp sample2.psp sample3.psp

$ cargo run --release --bin pop_var_caller -- var-calling \
    --reference ref.fa --output cohort.vcf.gz \
    --contamination-estimates cohort.contam.toml \
    sample1.psp sample2.psp sample3.psp

$ bcftools view cohort.vcf.gz | head -20
$ bcftools stats cohort.vcf.gz
```

(The `bcftools` step is the "open item" from the cohort VCF writer
review; it lands with this slice.)

## Risks and trade-offs

- **`var-calling-from-bam` lifetime juggling — and the
  no-intermediate-`.psp` rule.** Stage 1's `run_pileup` already
  does non-trivial `&mut` borrow ordering to thread
  `reader → BaqStream → ErrorSheddingAdapter → walker`
  (see comment at
  [cli.rs:320-327](../../src/pop_var_caller/cli.rs#L320)).
  Extracting `build_stage1_pipeline` to return an iterator
  without escaping borrows is the hardest part of this slice.

  **No temp-`.psp` fallback.** The whole point of
  `var-calling-from-bam` is to avoid the `.psp` intermediate;
  silently writing one to `tmp/` and reading it back would
  defeat the subcommand's reason for existing. If the
  borrow-free `build_stage1_pipeline` extraction proves
  intractable under the current architecture, the correct
  response is to **drop `var-calling-from-bam` from this slice
  entirely** and revisit when (a) the architecture has evolved
  to support it cleanly, or (b) we accept that one-off calls
  always go through `pileup` → `var-calling` as two commands.
  The other two subcommands (`var-calling`, `estimate-contamination`)
  are not at risk — they don't need the Stage 1 extraction.

- **Engine API for "no contamination".** The posterior engine
  currently takes a `ContaminationEstimates`. There is no
  documented `::zero()` sentinel. Two ways to bridge:
  - Add `ContaminationEstimates::zero(sample_names)` that
    populates `c_s = 0`, `q_b` = uniform (or whatever
    no-contamination defaults the engine treats as identity).
  - Use the existing
    [`ContaminationEstimateSource::Mixed`](../../src/var_calling/contamination_estimation.rs)
    slot which is reserved but unused in v1, and let the engine
    interpret missing entries as "no correction".
  Recommended: the first; explicit zero is easier to read and
  test. Add a one-liner constructor.

- **Sample-name basis.** The `.psp` header's `sample` field is
  the authoritative sample name throughout this slice. The
  artefact and the batches TSV reference samples by that name.
  Mismatch (artefact / TSV name doesn't match any `.psp`) is
  the *extras-are-ignored* case; the reverse (a `.psp` sample
  not present in the artefact) is the hard-error case. Per the
  decision in the design discussion.

- **Reference cross-check tolerance.** Different `.psp` files
  may have been produced against the same FASTA at different
  paths (e.g. `grch38.fa` vs `/data/refs/grch38.fa`). The
  cross-check compares basenames only and the per-contig MD5
  set on each PspReader header (which is a real content fingerprint).
  Path mismatches with matching MD5s are fine; MD5 mismatch is
  a hard error.

- **`emit_gp` payload size.** `GP` is `Number=G` (one float
  per genotype). For ploidy=2, `n_alleles=6`, that's
  `(6+1) choose 2 = 21` floats per sample. On a 1000-sample
  cohort with 100K variants, that's 100K × 1000 × 21 × 4 B ≈
  8 GiB. The flag's hidden default-off status protects naive
  invocations; the `--help` text says so explicitly.

- **Atomic rename.** Same write-tmp-then-rename convention as
  `run_pileup`, for the same reasons (partial VCFs are silently
  parseable; we adopt the safer convention even though
  bcftools/GATK don't). The convention is documented once
  in this plan and once in the `--help` epilogue.

- **Threading model.** `--threads` sizes the global rayon pool
  exactly once per process. The Stage 5 per-group merger and the
  posterior engine both use rayon. There is no cohort-wide
  cross-record parallelism in this slice — `rayon-over-records`
  is the explicit follow-up.

## Tasks

This slice is broken into 10 sequential tasks, each with one
deliverable that can be reviewed and merged before the next begins.
Most tasks have a hard dependency on an earlier one; the few that
can run in parallel are noted. Each task should end with `cargo
test` + `cargo clippy -D warnings` + `cargo fmt --check` clean.

### Task 1 — Engine-side `Config::new -> Result` validators

**Deliverable:** Every `Config` struct in the var-calling stack has
a validated `::new(...) -> Result<Self, _>` constructor enforcing
the ranges in the §"Config validation — engine + CLI" table.
`ContaminationEstimates::zero(sample_names)` sentinel constructor
lands alongside.

**Files touched:**
- `src/var_calling/posterior_engine.rs` — `PosteriorEngineConfig::new`,
  closes Mi12. `ContaminationEstimates::zero(sample_names)`.
- `src/var_calling/contamination_estimation.rs` —
  `ContaminationEstimationConfig::new`.
- `src/var_calling/variant_grouping.rs` — `GrouperConfig::new`.
- `src/var_calling/per_group_merger.rs` — `PerGroupMergerConfig::new`.
- `src/var_calling/vcf_writer/mod.rs` — `WriterConfig::new`; drop
  `default_filter_pass` field and `DEFAULT_FILTER_PASS` const.

**Verification:** boundary tests per range (just inside / on / just
outside) for every validated field. Existing call sites continue to
compile and pass — call sites that previously used `Default` or a
plain struct literal switch to `::new(...)?` with documented
defaults.

**Depends on:** —

### Task 2 — CLI value_parser helpers

**Deliverable:** A small `pop_var_caller::cli::parsers` module (or
inline in [cli.rs](../../src/pop_var_caller/cli.rs)) exposing
`parse_finite_positive`, `parse_finite_unit_interval`,
`parse_finite_open_unit_interval`, `parse_ploidy`,
`parse_max_alleles`, `parse_max_iterations`, `parse_max_gq_phred`,
`parse_min_major_fraction`, `parse_c_s_init`,
`parse_min_cohort_minor_fraction`, `parse_dust_window`,
`parse_pseudocount`. Each mirrors the engine-side range from Task 1.

**Files touched:**
- `src/pop_var_caller/cli.rs` (or new `cli/parsers.rs` submodule).

**Verification:** unit tests per parser — accepts in-range values,
rejects out-of-range / `NaN` / `±∞` / non-numeric.

**Depends on:** Task 1 (so the CLI parsers and the engine
constructors share the same numeric ranges by construction).

### Task 3 — Batch-assignment TSV reader

**Deliverable:** `src/pop_var_caller/batch_assignment.rs` with
`BatchAssignment::from_tsv(path)`, `BatchAssignment::empty()`,
`batch_for(sample)` (falls back to `"all_samples"`).

**Files touched:**
- `src/pop_var_caller/batch_assignment.rs` (new).
- `src/pop_var_caller/mod.rs` — module declaration.

**Verification:** unit tests for happy path, missing-header, empty
batch field, extras-ignored, default fallback. ~120 lines total.

**Depends on:** — (can run in parallel with Task 2).

### Task 4 — Contamination artefact TOML schema + I/O + validator

**Deliverable:** `src/pop_var_caller/contamination_artifact.rs`
with `ContaminationArtefact::{read, write,
to_estimates_for_samples}`, validator enforcing every rule in
§"Contamination artefact — TOML schema".

**Files touched:**
- `src/pop_var_caller/contamination_artifact.rs` (new).
- `src/pop_var_caller/mod.rs` — module declaration.

**Verification:** unit tests — TOML round-trip on a representative
artefact; validator rejects each illegal-state class (duplicate
batch id, dangling sample.batch, q_b row that doesn't sum to 1,
NaN/∞ in any numeric field, out-of-range `contamination_fraction`);
sample-name reconciliation (missing → error, extras → ignored).

**Depends on:** Task 1 (uses `ContaminationEstimates` /
`ContaminationEstimates::zero`).

### Task 5 — `estimate-contamination` subcommand

**Deliverable:** Working `pop_var_caller estimate-contamination`
subcommand: `EstimateContaminationArgs`, `run_estimate_contamination`,
`EstimateContaminationCliError`, dispatch from `cli.rs`.

**Files touched:**
- `src/pop_var_caller/estimate_contamination.rs` (new, ~400 lines).
- `src/pop_var_caller/mod.rs` — re-export.
- `src/pop_var_caller/cli.rs` — add `EstimateContamination` variant
  to `PopVarCallerCommand`.
- `src/main.rs` — dispatch the new variant.

**Verification:** unit tests inside the new module for the
args-to-config translation and any helpers. Whole-pipeline behaviour
gets covered by integration tests in Task 9.

**Depends on:** Tasks 1, 2, 3, 4.

### Task 6 — `var-calling` subcommand

**Deliverable:** Working `pop_var_caller var-calling` subcommand:
`VarCallingArgs`, `run_var_calling`, `VarCallingCliError`, dispatch
from `cli.rs`.

**Files touched:**
- `src/pop_var_caller/var_calling.rs` (new, ~500 lines).
- `src/pop_var_caller/mod.rs` — re-export.
- `src/pop_var_caller/cli.rs` — add `VarCalling` variant.
- `src/main.rs` — dispatch the new variant.

**Verification:** unit tests for args-to-config translation,
contamination-load + sample-reconciliation, and the
no-`--contamination-estimates` path producing the right zero
estimate.

**Depends on:** Tasks 1, 2, 4. (Independent of Task 5.)

### Task 7 — Stage 1 pipeline extraction (`stage1_pipeline.rs`)

**Deliverable:** `build_stage1_pipeline(...)` library helper that
returns a `PileupRecord` iterator without writing a `.psp`.
`run_pileup` refactored to call it for code reuse with the from-bam
path. Existing `pop_var_caller pileup` behaviour unchanged.

**Files touched:**
- `src/pop_var_caller/stage1_pipeline.rs` (new, ~150 lines).
- `src/pop_var_caller/cli.rs` — `run_pileup` body refactored;
  CLI surface untouched.
- `src/pop_var_caller/mod.rs` — module declaration.

**Verification:** existing pileup integration tests
(`tests/pileup_cli_integration.rs`) continue to pass — that's the
load-bearing regression gate. New helper-level unit tests for the
borrow-arena shape if non-trivial.

**Depends on:** —. **Risk:** per the "no temp-`.psp` fallback" rule
in §Risks, if this proves intractable, Task 8 is dropped from the
slice. Land Task 7 as a yes/no gate before starting Task 8.

### Task 8 — `var-calling-from-bam` subcommand (or drop)

**Deliverable:** Working `pop_var_caller var-calling-from-bam`
subcommand: `VarCallingFromBamArgs`, `run_var_calling_from_bam`,
`VarCallingFromBamCliError`, dispatch from `cli.rs`. **If Task 7
proves intractable, this task is dropped** and the implementation
report records why.

**Files touched:**
- `src/pop_var_caller/var_calling_from_bam.rs` (new, ~300 lines).
- `src/pop_var_caller/mod.rs` — re-export.
- `src/pop_var_caller/cli.rs` — add `VarCallingFromBam` variant.
- `src/main.rs` — dispatch the new variant.

**Verification:** unit tests for the args-to-config translation;
the k=1 merger wiring is exercised end-to-end in Task 9.

**Depends on:** Task 7, plus Tasks 1, 2, 6 (the from-bam orchestrator
reuses var-calling's Stage 3–6 setup).

### Task 9 — Integration tests + fixture

**Deliverable:** `tests/cohort_cli_integration.rs` covering all 10
cases listed in §"Integration tests", plus the fixture under
`tests/data/cohort/` (small reference FASTA + synthetic CRAM / `.psp`
generators).

**Files touched:**
- `tests/cohort_cli_integration.rs` (new).
- `tests/data/cohort/` (new fixture tree).

**Verification:** all 10 cases pass. The chained
`estimate-contamination → var-calling` end-to-end case is the
load-bearing one. Cases 9 / 10 (from-bam happy path + config
validation hard error) drop if Task 8 was dropped.

**Depends on:** Tasks 5, 6, optionally Task 8.

### Task 10 — Manual smoke + implementation report

**Deliverable:** A recorded smoke-run transcript
(`estimate-contamination → var-calling → bcftools view + stats`)
on a real-ish fixture; the implementation report under
[doc/devel/reports/implementations/pop_var_caller_cohort_cli_<date>.md](../reports/implementations/)
summarising what landed, what was deferred (notably whether Task 8
shipped or was dropped), and the working invocation lines.

**Files touched:**
- `doc/devel/reports/implementations/pop_var_caller_cohort_cli_<date>.md`
  (new).

**Verification:** `bcftools view <output>.vcf.gz | head` parses
without complaint; `bcftools stats` reports non-zero counts; PROJECT_STATUS
"open" items tagged "lands with the cohort CLI" on Stages 3 / 5 /
6 are moved to the appropriate done-state (or kept open with a
note if they slipped).

**Depends on:** Tasks 5, 6, 9, optionally Task 8.

## Validation

The slice is done when:

1. `cargo test` (lib + integration) is clean.
2. `cargo clippy --workspace --all-targets -- -D warnings` is clean.
3. `cargo fmt -- --check` is clean.
4. Manual smoke (above) produces a VCF that `bcftools view`
   parses without complaint and `bcftools stats` reports
   non-zero counts.
5. An implementation report under
   [doc/devel/reports/implementations/](../reports/implementations/)
   summarises what landed, what was deferred, and the smoke-test
   command lines that worked. The report also closes (or moves
   to a "remaining" follow-up list) the Stage 3/5/6 "open"
   items in PROJECT_STATUS that were tagged "lands with the
   cohort CLI".

## File touch list

New:

- `src/pop_var_caller/var_calling.rs` (~500 lines).
- `src/pop_var_caller/estimate_contamination.rs` (~400 lines).
- `src/pop_var_caller/var_calling_from_bam.rs` (~300 lines).
- `src/pop_var_caller/contamination_artifact.rs` (~350 lines
  inc. validator + tests).
- `src/pop_var_caller/batch_assignment.rs` (~120 lines
  inc. tests).
- `src/pop_var_caller/stage1_pipeline.rs` (~150 lines —
  extracted from `cli.rs` plus a thin facade).
- `tests/cohort_cli_integration.rs` (~600 lines for all 10
  cases).
- `tests/data/cohort/` — fixture FASTA + small synthetic CRAM
  / `.psp` generators.
- `doc/devel/reports/implementations/pop_var_caller_cohort_cli_<date>.md` —
  written when the slice lands.

Modified:

- `src/pop_var_caller/mod.rs` — declare new modules; re-export the
  three new `run_*` functions and arg structs.
- `src/pop_var_caller/cli.rs` — add three subcommand variants to
  `PopVarCallerCommand`; `run_pileup` body refactored to call
  `build_stage1_pipeline` for code reuse with the from-bam path.
- `src/var_calling/posterior_engine.rs` — `with_config -> Result`
  (Mi12); `PosteriorEngineConfig::new` validator;
  `ContaminationEstimates::zero(sample_names)` constructor.
- `src/var_calling/contamination_estimation.rs` —
  `ContaminationEstimationConfig::new -> Result` validator.
- `src/var_calling/variant_grouping.rs` —
  `GrouperConfig::new -> Result` validator.
- `src/var_calling/per_group_merger.rs` —
  `PerGroupMergerConfig::new -> Result` validator.
- `src/var_calling/vcf_writer/mod.rs` —
  `WriterConfig::new -> Result` validator. Drop `default_filter_pass`
  field and the `DEFAULT_FILTER_PASS` const (the flag is being
  dropped; v1 always emits PASS).
- `src/main.rs` — dispatch the three new subcommand variants.

Unmodified (verified):

- Everything under
  [src/var_calling/](../../src/var_calling/) except the validator
  additions noted above. Algorithmic code stays as it is.
- All of [src/per_sample_pileup/](../../src/per_sample_pileup/)
  (Stages 1–2). Untouched.

## Open follow-ups for future slices

These are out of scope here but the design accommodates them:

- `--regions BED` for restricting calls (standing item).
- `tabix .tbi` index alongside `.vcf.gz` outputs.
- `PL` `FORMAT` field forwarded from Stage 5.
- Per-sample contamination fraction in `INFO`.
- `bcftools stats` regression baseline checked into
  `tests/data/cohort/` and asserted against in CI (after the
  golden fixtures stabilise).
- `--external-allele-frequencies` for `estimate-contamination`
  (reference-panel `q_b`).
- Cross-record rayon parallelism (the deferred Stage 6 perf
  lever — see
  [posterior_engine_perf_2026-05-18_v2.md](../reports/implementations/posterior_engine_perf_2026-05-18_v2.md)).
- `--per-group-batch-size` exposure (one of the levers in the
  parallelisation-tuning pass that runs after this slice).
- `pop_var_caller pileup` itself extracted into its own sibling
  module (currently lives in [cli.rs](../../src/pop_var_caller/cli.rs);
  this slice already creates the sibling layout for the new
  subcommands, so finishing the symmetry is a small mechanical
  follow-up).
