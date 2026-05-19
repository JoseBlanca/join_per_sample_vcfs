# `pop_var_caller` cohort CLI — review follow-up

Implementation plan for the **deferred** items from the
[cohort CLI review](../reports/reviews/cohort_cli_2026-05-19.md)
and its [fix-application report](../reports/reviews/cohort_cli_2026-05-19_applied.md).
The first fix-application run (db1ec2a) closed 20 findings; this plan
covers the 16 that were Deferred plus the two doc-only follow-ups
left open (M5 and M9).

Spec context unchanged: this is still glue work on the
[cohort CLI plan](pop_var_caller_cohort_cli.md). No new algorithms;
no new subcommands. Every change here improves correctness, API
hygiene, or maintainability of the existing surface.

## What's being deferred from

The cohort CLI slice landed at commits `1523049..147e435` and is
fixes-applied at `db1ec2a`. The review's open items are now:

**Major design-calls (5):**
- **M4** — Post-construction mutation of
  `posterior_cfg.contamination` bypasses any future engine-side
  validation. Pairs with **Mi2** (8 positional `f64` args on
  `PosteriorEngineConfig::new` are silently swappable) and
  **Mi21** (`ContaminationEstimationConfig::validate(&self)` is
  the only Stage-config without a `new(...) -> Result` shape).
- **M8** — `WriterConfig` lacks `#[non_exhaustive]`; the
  `default_filter_pass` drop was an unannounced breaking change
  for any downstream crate using struct-literal construction.
  Pairs with **Mi5** (the on-disk `ContaminationArtefact` /
  `Provenance` / `ProvenanceInputs` / `BatchEntry` /
  `SampleEntry` structs likewise lack `#[non_exhaustive]`) and
  **Mi6** (`provenance.version` is recorded but never validated
  on read).
- **M10** — `VarCallingFromBamArgs` duplicates ~30 fields of
  `VarCallingArgs` + `PileupArgs` verbatim. Three places for the
  same flag-set.
- **M11** — Cohort pipeline wiring (merger → DUST → grouper →
  per-group → posterior → writer) duplicated near-verbatim
  between `var_calling.rs` and the
  `run_cohort_pipeline_for_single_sample` helper.
- **M13** — Three `rayon::ThreadPoolBuilder::build_global()`
  sites with no in-process coordination. Second-call hazard for
  library consumers or a future multi-subcommand test runner.

**Paired Minors (11):**
- **Mi2 / Mi21** — pair with M4 (config-construction shape).
- **Mi5 / Mi6** — pair with M8 (on-disk type hygiene +
  schema-versioning gate).
- **Mi8** — `basename` / `format_md5_hex` / `current_command_line` /
  `rfc3339_now` / `civil_from_days` duplicated across 3–4 files;
  extract `pop_var_caller::common`.
- **Mi13** — file `contamination_artifact.rs` (American) holds
  type `ContaminationArtefact` (British). Pure rename.
- **Mi14** — CLI knob `--min-batch-size` diverges from engine
  field `min_batch_size_for_contamination`. Same concept, three
  names.
- **Mi18** — `build_artefact_from_estimates` linear-scans for a
  "representative sample" per batch. Engine accessor shape
  doesn't fit batch-keyed export; needs a
  `q_b_per_batch()` accessor on `ContaminationEstimates`.
- **Mi19** — Magic `64 * 1024` BufReader/BufWriter capacity
  duplicated at 8+ sites. Promote `DEFAULT_BUFFERED_IO_CAPACITY`.
- **Mi20** — Integration-test fixture-helper duplication between
  `tests/pileup_cli_integration.rs` and
  `tests/cohort_cli_integration.rs`. Extract `tests/common/mod.rs`.
- **Mi23** — Load-bearing
  `estimate-contamination → var-calling` chain integration test
  + the end-to-end CRAM-truncation test for M1/M2's walker-error
  path.

**Doc-only follow-ups (2):**
- **M5 follow-up** — Wire real FASTA→.psp MD5 enforcement. The
  v1 contract (basename-only) was documented in the previous
  fix-application run; this is the long-term answer.
- **M9 follow-up** — Generalise `ErrorSheddingAdapter<I, T, E>`
  so the from-bam walker shim reuses it instead of the
  open-coded inline `.scan()`. The stale doc reference was
  fixed; this is the actual code dedup.

## Five waves

The work splits cleanly along independent dimensions. Each wave is
its own mergeable slice with its own validation pass and its own
fix-application-style report; nothing in one wave blocks another
except where called out explicitly under "Dependencies" below.

| Wave | Theme | Findings | Approx. LOC | Risk |
|---|---|---|---|---|
| 1 | Public-API hygiene | M8, Mi5, Mi6, Mi13 | ~300 | Low — additive `#[non_exhaustive]` + file rename |
| 2 | Config-construction discipline | M4, Mi2, Mi21, Mi14 | ~400 | Medium — touches the engine-side `Config` API |
| 3 | Shared-infrastructure refactor | Mi8, Mi19, M9-followup, M10, M11, Mi18 | ~1500 | High — large mechanical rewrite + clap-derive flatten |
| 4 | Concurrency policy | M13 | ~150 | Low — single new helper |
| 5 | Test infrastructure + missing coverage | Mi20, Mi23, M1/M2-followup, M5-followup | ~600 | Medium — test fixtures + real-MD5 wiring is genuinely new code |

Dependencies (informal):

- **Wave 2 ⇨ Wave 3** — if Wave 2 adopts a builder pattern for
  configs, Wave 3's shared driver wants those builders so the
  per-stage config plumbing inside `drive_cohort_pipeline` reads
  cleanly. Wave 3 can still land if Wave 2 picks `validate(&self)`
  for all configs instead; just a different surface.
- **Wave 3 ⇨ Wave 5** — Mi20 (`tests/common/mod.rs`) wants the
  Wave 3 shared helpers (`basename`, `format_md5_hex`) extracted
  first so the integration tests can reuse them.
- **Waves 1 / 4 are independent** of the others.

Waves can be reordered. The grouping is by theme, not by
dependency order; the table above is what I'd recommend.

---

## Wave 1 — Public-API hygiene

### Scope

Lock the public surface of every type we expect downstream crates
or future slices to consume. Today the slice ships `pub` structs
with all-`pub` fields and no `#[non_exhaustive]`, which means any
future field addition is a breaking change. The dropped
`WriterConfig.default_filter_pass` field on the cohort CLI slice
is the proof-of-life that this matters.

### In scope

- **M8** — Add `#[non_exhaustive]` to
  [WriterConfig](../../src/var_calling/vcf_writer/mod.rs) plus a
  builder-style API:
  ```rust
  #[non_exhaustive]
  #[derive(Debug, Clone)]
  pub struct WriterConfig {
      pub output: PathBuf,
      pub emit_gp: bool,
  }
  impl WriterConfig {
      pub fn new(output: PathBuf) -> Self {
          Self { output, emit_gp: DEFAULT_EMIT_GP }
      }
      pub fn with_emit_gp(mut self, emit_gp: bool) -> Self {
          self.emit_gp = emit_gp;
          self
      }
  }
  ```
  Every in-tree consumer (CLI orchestrators + integration tests)
  migrates from struct-literal to
  `WriterConfig::new(path).with_emit_gp(args.emit_gp)`.

- **Mi5** — Add `#[non_exhaustive]` to the five on-disk artefact
  structs in
  [contamination_artifact.rs](../../src/pop_var_caller/contamination_artifact.rs):
  `ContaminationArtefact`, `Provenance`, `ProvenanceInputs`,
  `BatchEntry`, `SampleEntry`. Test fixtures inside the crate
  continue to compile (struct-literal is fine within the crate);
  the attribute only locks the surface against external callers.

- **Mi6** — Add a `SUPPORTED_VERSIONS: &[&str]` whitelist and a
  read-time gate on `provenance.version`:
  ```rust
  const SUPPORTED_VERSIONS: &[&str] = &["0.1.0"];

  // inside Self::validate():
  if !SUPPORTED_VERSIONS.contains(&self.provenance.version.as_str()) {
      return Err(ContaminationArtefactError::UnsupportedVersion {
          got: self.provenance.version.clone(),
          supported: SUPPORTED_VERSIONS,
      });
  }
  ```
  Add `ContaminationArtefactError::UnsupportedVersion` variant.
  Lenient policy (warn on minor mismatch, error on major)
  considered and rejected — strict allow-list is simpler and
  catches the failure mode the field exists to catch.

- **Mi13** — Rename file
  `src/pop_var_caller/contamination_artifact.rs` to
  `contamination_artefact.rs` so the path uniformly uses the
  British "artefact" (matching the type, the error, and every
  consumer). Mechanical: `git mv` + update `mod.rs:8` +
  update all `use` lines (currently 5 sites:
  `var_calling.rs`, `estimate_contamination.rs`,
  `mod.rs` re-export, plus 2 test files).

### Out of scope (this wave)

- Renaming the type to American spelling (`ContaminationArtifact`).
  That would touch every consumer; Mi13's file rename is the
  smaller direction and the project is already British-default
  in identifiers.
- Other artefact-format changes (new fields, new validators).
  This wave is purely about *guarding* the existing format
  against future drift.

### Tasks

#### Task 1.1 — `#[non_exhaustive]` on `WriterConfig` + builder

Files: `src/var_calling/vcf_writer/mod.rs`; migrate call sites in
`src/pop_var_caller/var_calling.rs`,
`src/pop_var_caller/var_calling_from_bam.rs`,
`tests/cohort_vcf_writer_integration.rs`.

Verification: a new doctest on `WriterConfig::new` showing the
builder pattern; existing 880+ tests stay green.

#### Task 1.2 — `#[non_exhaustive]` on `ContaminationArtefact` and its sub-structs

Files: `src/pop_var_caller/contamination_artifact.rs`.

Verification: existing 19 artefact tests stay green (all fixtures
build via struct literal from inside the crate, which is allowed
on `#[non_exhaustive]` types).

#### Task 1.3 — `provenance.version` validation gate

Files: `src/pop_var_caller/contamination_artifact.rs`.

Verification: new test
`read_rejects_artefact_with_unknown_version` — write a TOML with
`version = "99.0.0"`, expect `Err(UnsupportedVersion)`.

#### Task 1.4 — File rename `contamination_artifact.rs` → `contamination_artefact.rs`

Mechanical. `git mv`, update `mod.rs:8`, sweep `use` lines.

Verification: full `cargo build --tests` clean. No new tests.

**Depends on:** Tasks 1.1–1.3 land first to avoid double-rebasing
the import lines they touch.

---

## Wave 2 — Config-construction discipline

### Scope

Pick **one** validated-construction pattern and apply it to every
Stage 3–6 config. Today the slice has a mix:

- `DustFilterConfig::new -> Result` (positional, 2 args) ✓
- `GrouperConfig::new -> Result` (positional, 1 arg) ✓
- `PerGroupMergerConfig::new -> Result` (positional, 3 args) ✓
- `PosteriorEngineConfig::new -> Result` (positional, **8 args**, mostly `f64`) — Mi2 swap risk
- `PosteriorEngineConfig.contamination` — set by direct field
  mutation after `new`, bypassing any future validation. M4.
- `ContaminationEstimationConfig::validate(&self)` (struct
  literal + `validate()?` at the call site, 14+ fields). Mi21.

The inconsistency is the actual finding. Two reasonable
end-states, both eliminate the silent-swap and
post-construction-mutation hazards:

### Decision recorded — **Option C (hybrid)**

Locked in 2026-05-19 during plan review. Builder for the two
configs with sharp invariants (`WriterConfig`,
`PosteriorEngineConfig`); validate-after-build for the wide
numeric config (`ContaminationEstimationConfig`). The three
options are kept below for context — but the wave implements
Option C.

#### Option A — Validate-after-build everywhere

Every config exposes `pub` fields; construction is via struct
literal (typically anchored on `..Self::with_project_defaults()`);
`config.validate()?` is the single enforcement point.

- ✓ Names every field at the call site (cures Mi2's silent-swap).
- ✓ Friendly for many-field configs (cures the 8/14-arg
  hostility).
- ✓ Easy migration from current `validate(&self)` shape on
  `ContaminationEstimationConfig`.
- ✗ Doesn't help M4: a caller can mutate `cfg.contamination = …`
  after `validate()`. Need to either privatise the field or
  require `validate()` re-call before consumption.

#### Option B — Builder everywhere

Privatise fields; expose `new(required_args)` + per-field
`with_*` setters that **return `Result`**. Final `build()` runs
the full validation.

- ✓ Solves M4 by construction — `cfg.contamination` is
  unreachable except via `with_contamination(estimates)` which
  validates.
- ✓ Stricter than (A); cannot bypass validation.
- ✗ Verbose at the call site for many-field configs (`a().b().c()…`).
- ✗ Larger API surface; more boilerplate.

#### Option C — Hybrid

`WriterConfig` and `PosteriorEngineConfig` adopt the builder
pattern (B) — they're the user-facing surface and they have
fields with non-trivial invariants. `ContaminationEstimationConfig`
(14+ fields, mostly numeric, no cross-field invariants) keeps
the validate-after-build pattern (A) it already has.

This matches the Mi21 reviewer's caveat: "14+ positional args
would be hostile". A builder for a 14-field config is similarly
hostile.

**Chosen:** Option C — builder where the API is sharp,
validate-after-build where the API is wide. Matches the project's
existing instinct.

### In scope

(Assuming Option C; if (A) or (B) is chosen, adapt.)

- **M4** — Privatise
  [PosteriorEngineConfig.contamination](../../src/var_calling/posterior_engine.rs)
  field (and any other field a future validator would gate).
  Expose
  ```rust
  pub fn with_contamination(
      mut self,
      contamination: Option<ContaminationEstimates>,
  ) -> Result<Self, PosteriorEngineConfigError> {
      // future cross-field invariants validated here
      self.contamination = contamination;
      Ok(self)
  }
  ```
  Update both CLI orchestrators (`var_calling.rs:367` and
  `var_calling_from_bam.rs:455`) to use it.

- **Mi2** — Rebuild `PosteriorEngineConfig::new` into a
  small-required-args constructor + per-field `with_*` setters
  (or, if Option A wins the discussion, switch to struct-literal
  + `validate()` like `ContaminationEstimationConfig`). The
  8-`f64`-in-a-row signature has to go.

- **Mi21** — Either leave
  `ContaminationEstimationConfig::validate(&self)` as the
  documented path (Option C), or wrap it in a `build() -> Result`
  pattern (Option B). Either way, document the "construct → call
  X" contract on the type.

- **Mi14** — Rename the CLI flag `--min-batch-size` →
  `--min-batch-size-for-contamination` and the corresponding
  struct field. Mirror the rename in the artefact's parameter-map
  key. **User-visible CLI change** — anyone scripting with the
  old name breaks. Worth doing because the artefact is the
  contract between two subcommands and the same concept should
  carry the same name everywhere.

### Out of scope

- Restructuring `DustFilterConfig`, `GrouperConfig`,
  `PerGroupMergerConfig`. These already have small validated
  `new()` constructors and their argument counts (2 / 1 / 3) are
  not hostile.

### Tasks

#### Task 2.1 — Choose the shape

**Decision:** Option C (hybrid). Locked in 2026-05-19 during plan
review.

#### Task 2.2 — `PosteriorEngineConfig` refactor (M4 + Mi2)

Files: `src/var_calling/posterior_engine.rs`,
`src/pop_var_caller/var_calling.rs`,
`src/pop_var_caller/var_calling_from_bam.rs`.

Verification: existing 14 `config_new_*` tests adapt to the new
shape; the post-construction-mutation seam disappears.

#### Task 2.3 — `ContaminationEstimationConfig` (Mi21)

Files: `src/var_calling/contamination_estimation.rs`. Likely a
doc-only change under Option C; a wrap-in-build()-call under
Option B.

#### Task 2.4 — `--min-batch-size-for-contamination` rename (Mi14)

Files: `src/pop_var_caller/estimate_contamination.rs`,
`src/var_calling/contamination_estimation.rs`,
`src/pop_var_caller/contamination_artifact.rs`
(parameter-map key),
`tests/cohort_cli_integration.rs`.

Verification: parser test (already exists) renames; round-trip
artefact test verifies the new key in `parameters`.

---

## Wave 3 — Shared-infrastructure refactor

### Scope

The biggest wave. Three orchestrators (`var_calling.rs`,
`var_calling_from_bam.rs`, `estimate_contamination.rs`) currently
duplicate:

1. ~5 helper functions verbatim (`basename`,
   `format_md5_hex`, `current_command_line`, `rfc3339_now`,
   `civil_from_days`) across 2–4 files each.
2. The Stage-1 args + Stage 3–6 args, three places in clap-derive
   form (`PileupArgs`, `VarCallingArgs`,
   `VarCallingFromBamArgs`).
3. The cohort pipeline wiring (merger → DUST → grouper →
   per-group → posterior → writer) twice between
   `run_var_calling` and `run_cohort_pipeline_for_single_sample`.

This wave deduplicates all three, plus the engine-side accessor
needed to deduplicate `build_artefact_from_estimates`.

### In scope

- **Mi8** — Create
  [src/pop_var_caller/common.rs](../../src/pop_var_caller/) with
  `pub(crate)` versions of `basename`, `format_md5_hex`,
  `current_command_line`, `rfc3339_now`, `civil_from_days`, plus
  the project-wide buffered-IO capacity constant (Mi19). Every
  current call site `use super::common::*;` and the local
  duplicates disappear (including their identical unit tests).

- **Mi19** — `pub(crate) const DEFAULT_BUFFERED_IO_CAPACITY: usize = 64 * 1024;`
  in `common.rs`. Replace the 3+ in-scope sites + leave the
  orthogonal `BLOCK_HEADER_READ_CAP` (different intent) alone.

- **M9 follow-up** — Generalise `ErrorSheddingAdapter` in
  [error_bridge.rs](../../src/pop_var_caller/cli/error_bridge.rs)
  from `Item = PreparedRead, E = CramInputError` (hard-coded) to
  `ErrorSheddingAdapter<I, T, E>` parametric. The from-bam walker
  shim then reuses the named adapter instead of the open-coded
  inline `.scan()` (and the doc reference at
  `var_calling_from_bam.rs:13` finally points at a real type).

- **M10** — Factor shared clap-derive sub-structs:
  ```rust
  // src/pop_var_caller/cli/shared_args.rs (new)
  #[derive(Debug, Args, Clone)]
  pub struct Stage1Args {
      /* min_mapq, no_baq, BAQ HMM, walker, CRAM-input filters */
  }

  #[derive(Debug, Args, Clone)]
  pub struct CohortPipelineArgs {
      /* dust, var_group_max_span, max_alleles_per_var,
         posterior engine, emit_gp, ploidy */
  }
  ```
  Then:
  ```rust
  pub struct PileupArgs {
      pub reference: PathBuf, pub output: PathBuf,
      pub threads: Option<usize>, pub crams: Vec<PathBuf>,
      #[command(flatten)] pub stage1: Stage1Args,
  }

  pub struct VarCallingArgs {
      pub reference: PathBuf, pub output: PathBuf,
      pub threads: Option<usize>, pub psp_files: Vec<PathBuf>,
      pub contamination_estimates: Option<PathBuf>,
      pub no_complexity_filter: bool,
      #[command(flatten)] pub cohort: CohortPipelineArgs,
  }

  pub struct VarCallingFromBamArgs {
      pub reference: PathBuf, pub output: PathBuf,
      pub threads: Option<usize>, pub crams: Vec<PathBuf>,
      pub no_complexity_filter: bool,
      #[command(flatten)] pub stage1: Stage1Args,
      #[command(flatten)] pub cohort: CohortPipelineArgs,
  }
  ```
  Every `args.foo` access becomes `args.stage1.foo` or
  `args.cohort.foo`; the `*_config_from_args` helpers look one
  field deeper. The `--help` output and CLI surface are
  unchanged (clap-derive `flatten` is transparent at the user
  level).

- **M11** — Extract `drive_cohort_pipeline` in a new
  `src/pop_var_caller/cohort_driver.rs`:
  ```rust
  pub(crate) struct CohortPipelineParams<'a> {
      pub no_complexity_filter: bool,
      pub dust_cfg: DustFilterConfig,
      pub grouper_cfg: GrouperConfig,
      pub per_group_cfg: PerGroupMergerConfig,
      pub posterior_cfg: PosteriorEngineConfig,
      pub fetcher: SharedRefFetcher,
      pub chromosomes: &'a [ParsedChromosome],
  }

  pub(crate) fn drive_cohort_pipeline<I, E>(
      upstream: I,
      sample_names: Vec<String>,
      params: CohortPipelineParams<'_>,
      output: &Path,
      metadata: CohortMetadata,
      writer_cfg: WriterConfig,
  ) -> Result<u64, E>
  where
      I: Iterator<Item = Result<PileupRecord, PspReadError>>,
      E: From<DustFilterError> + From<GrouperError>
          + From<PerGroupMergerError> + From<PosteriorEngineError>
          + From<VcfWriteError>,
  { /* the merger → DUST → grouper → per-group → posterior → writer
       wiring, plus the tmp-cleanup loop M6 added */ }
  ```
  Both `run_var_calling` and
  `run_cohort_pipeline_for_single_sample` call it. The `Box<dyn
  Iterator<Item = Result<_, GrouperError>>>` shape lives in one
  place.

- **Mi18** — Add `pub(crate) fn q_b_per_batch(&self) -> &[[f64; N_ALLELE_CLASSES]]`
  to
  [ContaminationEstimates](../../src/var_calling/contamination_estimation.rs).
  Rewrite `build_artefact_from_estimates` to index it directly;
  drops the linear-scan `position()` workaround and the
  `[0.0; 3]` fall-through (which can become `expect(…)` since
  `build_dense_batches` guarantees every dense batch has ≥1
  sample).

### Out of scope (this wave)

- Performance work. The refactor is byte-for-byte equivalent at
  runtime; the cohort driver is the same iterator chain in one
  place rather than two.
- New tests. The existing integration tests cover the call paths
  end-to-end; the refactor's job is to leave them passing.

### Tasks

#### Task 3.1 — `pop_var_caller::common` module (Mi8 + Mi19)

Files: new `src/pop_var_caller/common.rs` + 4 consumer files. Lift
~80 lines into the shared module; drop ~150 lines of duplication
(helpers + their identical unit tests).

Verification: All `pop_var_caller` lib tests still pass; the
single-set unit tests in `common.rs` cover what was previously
double-tested.

#### Task 3.2 — `ErrorSheddingAdapter` generalisation (M9 follow-up)

Files: `src/pop_var_caller/cli/error_bridge.rs`,
`src/pop_var_caller/var_calling_from_bam.rs`.

The existing CRAM-input adapter call site in
`stage1_pipeline.rs` switches to the parametric form by
inference (no code change needed at the call site beyond the
type annotation).

Verification: All existing tests pass; the walker-shim site at
`var_calling_from_bam.rs:503` (now in the extracted helper) uses
the named adapter.

#### Task 3.3 — Mi18 engine accessor

Files: `src/var_calling/contamination_estimation.rs` (add
accessor), `src/pop_var_caller/estimate_contamination.rs`
(use it in `build_artefact_from_estimates`).

Verification: existing 19 artefact tests pass; the linear-scan
finder is gone.

#### Task 3.4 — Shared `Stage1Args` + `CohortPipelineArgs` (M10)

Files: new `src/pop_var_caller/cli/shared_args.rs`; modify
`cli.rs`, `var_calling.rs`, `var_calling_from_bam.rs` to use the
flatten sub-structs; update the `*_config_from_args` helpers to
look one field deeper.

Risk: this is the largest mechanical sweep in the plan. ~30
fields × 3 args structs. clap-derive's `flatten` is well-trodden
but the type changes propagate to every test fixture that builds
an args struct by literal.

Verification: existing integration tests (which build the args
struct via literal) update mechanically; `--help` output
unchanged.

#### Task 3.5 — `drive_cohort_pipeline` shared driver (M11)

Files: new `src/pop_var_caller/cohort_driver.rs`; modify
`var_calling.rs` and `var_calling_from_bam.rs` to delegate to it.

Verification: existing happy-path integration tests cover both
call sites end-to-end.

**Depends on:** Task 3.4 (the cohort args struct is a natural
parameter for the driver).

---

## Wave 4 — Concurrency policy

### Scope

**M13** — Centralise rayon-pool configuration in a single
idempotent helper:
```rust
// src/pop_var_caller/common.rs (or a new cli/rayon.rs)
use std::sync::OnceLock;

static POOL_CONFIGURED: OnceLock<()> = OnceLock::new();

pub fn configure_rayon_pool(n: Option<usize>) -> Result<(), rayon::ThreadPoolBuildError> {
    let Some(n) = n else { return Ok(()); };
    if POOL_CONFIGURED.get().is_some() {
        // already configured by an earlier call in this process —
        // no-op rather than error
        return Ok(());
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global()?;
    let _ = POOL_CONFIGURED.set(());
    Ok(())
}
```

Each of the three drivers (`run_pileup`, `run_var_calling`,
`run_estimate_contamination`, `run_var_calling_from_bam`) calls
`configure_rayon_pool(args.threads)?`. The drivers' own
`RayonAlreadyConfigured` variants are kept but their `?` site is
the helper, not `build_global` directly.

### Decision recorded — silent no-op

Locked in 2026-05-19 during plan review. Matches how most Rust
CLIs handle rayon's process-global rule and hides the seam from
library consumers and any future multi-subcommand test runner.

### In scope

- **M13** — The helper + migrations of the four call sites + a
  serial integration test (using `serial_test::serial`) that
  calls two `run_*` functions back-to-back with `--threads`
  set, asserts both succeed.

### Out of scope

- Any change to rayon's pool semantics. The helper just gates
  the `build_global` call.

### Tasks

#### Task 4.1 — `configure_rayon_pool` helper

Files: `src/pop_var_caller/common.rs` (or a new submodule).

#### Task 4.2 — Migrate the four call sites

Files: `src/pop_var_caller/cli.rs:301-306`,
`src/pop_var_caller/var_calling.rs:317-323`,
`src/pop_var_caller/var_calling_from_bam.rs:438-443`,
`src/pop_var_caller/estimate_contamination.rs:304-309`.

#### Task 4.3 — Serial integration test

Add `serial_test` to dev-dependencies if not already present.
Add test `rayon_pool_survives_second_call_with_threads` to
`tests/cohort_cli_integration.rs`.

**Depends on:** Task 3.1 (if `configure_rayon_pool` lives in
`common.rs`).

---

## Wave 5 — Test infrastructure + missing coverage

### Scope

The slice closed integration coverage on the happy paths but
deferred several load-bearing tests. This wave adds them, plus
the optional MD5-wiring follow-up (M5).

### In scope

- **Mi20** — `tests/common/mod.rs` with the shared CRAM-builder /
  FASTA-builder / `make_psp_for_sample` fixture helpers.
  `tests/pileup_cli_integration.rs` and
  `tests/cohort_cli_integration.rs` both `mod common;`.

- **Mi23** — The load-bearing chained test:
  `estimate_contamination_then_var_calling_chain`. Build 2–3
  `.psp`s, run `run_estimate_contamination` to a temp TOML, then
  `run_var_calling --contamination-estimates <that>`, assert the
  resulting VCF parses. Plus the deferred unit-level coverage
  the review listed under §8 "Missing tests":
  `estimate_contamination_rejects_reference_mismatch`,
  `var_calling_reports_reference_mismatch_for_psp_with_different_header_basename`,
  `to_estimates_round_trips_floored_batch_with_zero_qb`,
  `to_estimates_orders_dense_batches_by_cohort_first_seen`.

- **M1/M2 follow-up** — End-to-end CRAM-fixture test for the
  walker-error path. Synthesize reads dense enough to trip
  `max_active_reads = 1`; assert `Err(Walker(_))` and that the
  output VCF does not exist. The factored
  `prefer_upstream_or_closure` helper (added in db1ec2a) is
  unit-tested already; this test covers the *wiring*.

- **M5 follow-up — *in scope***. Wire real FASTA→.psp MD5
  enforcement in `run_var_calling` and
  `run_estimate_contamination`. Locked in 2026-05-19 during plan
  review. Either:
  - Read the `.fai` and compute per-contig MD5 by streaming the
    FASTA (correct but adds I/O cost on every run), or
  - Use the in-memory cache that
    [SyncRefFetcher](../../src/per_sample_pileup/ref_fetcher.rs)
    already maintains, hashing each contig as it's fetched the
    first time and comparing against
    `readers[0].header().chromosomes[*].md5`.

  The second is cheaper if `SyncRefFetcher` already exposes the
  fetched bytes; needs an audit. Failure mode: typed error
  `VarCallingCliError::FastaContigMd5Mismatch { contig,
  fasta_md5, psp_md5 }`.

  Worth doing? The plan recommends it because the silent-wrong-
  reference failure mode is exactly the kind of correctness bug
  the project's "no silent intermediates" principle exists to
  prevent. But it's a behaviour change, not a refactor — happy
  to drop if the user prefers the basename-only v1 contract to
  remain the long-term answer.

### Out of scope

- Adding more end-to-end CRAM smoke tests beyond what the review
  listed. The 7 plan-listed integration cases the slice deferred
  are mostly covered in spirit by per-module unit tests; only the
  chain test (Mi23) was load-bearing.

### Tasks

#### Task 5.1 — `tests/common/mod.rs`

Files: new `tests/common/mod.rs`; modify
`tests/pileup_cli_integration.rs` and
`tests/cohort_cli_integration.rs`.

Move `build_fasta`, `build_sam_header`, `build_cram`,
`read_record`, `pileup_args` + `make_psp_for_sample` into the
shared module.

**Depends on:** Task 3.1 (the common helpers in
`pop_var_caller::common` are referenced by `make_psp_for_sample`
via the public `pop_var_caller::PileupArgs` import; no direct
test-side dependency).

#### Task 5.2 — Chained integration test (Mi23)

Files: `tests/cohort_cli_integration.rs`. ~50 LOC.

Add `estimate_contamination_args(reference, output, psp_files)`
factory; assert the chained call produces a parseable VCF.

#### Task 5.3 — Missing-coverage unit tests

Files: `src/pop_var_caller/contamination_artifact.rs`,
`src/pop_var_caller/var_calling.rs`,
`tests/cohort_cli_integration.rs`. Add the four tests named in
§Mi23 above.

#### Task 5.4 — Walker-error CRAM-fixture test (M1/M2 follow-up)

Files: `tests/cohort_cli_integration.rs`. Synthesize a CRAM with
reads dense enough to trip `max_active_reads = 1`; assert
`Err(VarCallingFromBamCliError::Walker(_))` and absence of the
output VCF.

Costs: ~80 LOC including the synthetic-CRAM scaffolding.

#### Task 5.5 — *Optional.* FASTA→.psp MD5 wiring (M5 follow-up)

Audit `SyncRefFetcher`'s exposed surface; pick the
cheap-or-streaming implementation; add the typed error and
matching test. ~100 LOC + a tiny test.

---

## Risks and trade-offs

- **Wave 2 design decision is load-bearing.** If Option A wins,
  Wave 3 is simpler (no builder threading through configs). If
  Option B wins, every config gains a builder and Wave 3's
  driver signature gets shorter. Discussing Wave 2 first makes
  Wave 3's task list less ambiguous.

- **M14's CLI flag rename is the only user-visible behaviour
  change in this whole plan.** Anyone who scripted
  `pop_var_caller estimate-contamination --min-batch-size 5`
  needs to update to `--min-batch-size-for-contamination 5`.
  If we want to soften the migration, we can ship the long form
  with an alias for the short form for one release, then drop
  the alias. clap-derive supports `#[arg(long, alias = "...")]`.

- **Wave 3 is a large mechanical sweep.** The
  `#[command(flatten)]` change (~30 fields × 3 args structs) and
  the `drive_cohort_pipeline` extraction together touch most of
  the cohort CLI surface. Recommend landing as five sequential
  commits (one per Task 3.1–3.5) so each rebases cleanly and
  each green-validation gate is small enough to bisect.

- **Schema versioning (Mi6) is forward-looking only.** Today
  there's only ever been v0.1.0; the gate is dead code on day
  one. The cost is two lines + one test; the value lands when
  v0.2.0 (or whatever) breaks a field.

- **The M5 follow-up (FASTA→.psp MD5 wiring) might prove more
  expensive than it looks** depending on what `SyncRefFetcher`
  exposes. If the cheap path requires a fetcher refactor,
  consider deferring to a separate slice.

## Validation

Each wave ends with the same gate as the original cohort CLI
slice:

1. `cargo fmt --check` clean
2. `cargo clippy --workspace --all-targets -- -D warnings` clean
3. `cargo test --all-targets` clean
4. `cargo doc --no-deps` clean (the pre-existing
   `posterior_engine.rs` `ExactMath` errors remain out of scope
   — that's Stage 6 Mi21)
5. Fix-application-style report under
   `doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_v2.md`
   (or per-wave `_wave1`, `_wave2`, …) summarising what landed.

## File touch list (estimated, by wave)

**Wave 1 (Public-API hygiene):**

- M: `src/var_calling/vcf_writer/mod.rs`,
  `src/pop_var_caller/var_calling.rs`,
  `src/pop_var_caller/var_calling_from_bam.rs`,
  `tests/cohort_vcf_writer_integration.rs`,
  `src/pop_var_caller/contamination_artifact.rs` (then renamed to `_artefact.rs`).
- R: `src/pop_var_caller/mod.rs`, every `use` site of the artefact module.

**Wave 2 (Config discipline):**

- M: `src/var_calling/posterior_engine.rs`,
  `src/var_calling/contamination_estimation.rs`,
  `src/pop_var_caller/var_calling.rs`,
  `src/pop_var_caller/var_calling_from_bam.rs`,
  `src/pop_var_caller/estimate_contamination.rs`,
  `src/pop_var_caller/contamination_artifact.rs`.

**Wave 3 (Refactor):**

- N: `src/pop_var_caller/common.rs`,
  `src/pop_var_caller/cli/shared_args.rs`,
  `src/pop_var_caller/cohort_driver.rs`.
- M: every orchestrator + `error_bridge.rs` + Cargo.toml if a new
  dev-dep is needed for tests.

**Wave 4 (Rayon):**

- M: 4 driver files; possibly a new submodule in
  `src/pop_var_caller/`.
- C: `Cargo.toml` (`serial_test` dev-dep if not already).

**Wave 5 (Tests + MD5 wiring):**

- N: `tests/common/mod.rs`.
- M: `tests/pileup_cli_integration.rs`,
  `tests/cohort_cli_integration.rs`,
  `src/pop_var_caller/var_calling.rs` (M5 path),
  `src/pop_var_caller/estimate_contamination.rs` (M5 path).

## Decisions recorded (2026-05-19 plan review)

1. **Wave 2 shape** — **Option C (hybrid)**. Builder for
   `WriterConfig` and `PosteriorEngineConfig`; validate-after-build
   for `ContaminationEstimationConfig`.
2. **Mi14 flag rename** — **Rename outright** to
   `--min-batch-size-for-contamination` (struct field + artefact
   parameter-map key follow). Hard break for anyone scripting the
   old name; justified because the artefact is the contract
   between two subcommands.
3. **Wave 4 second-call policy** — **Silent no-op**. Match
   conventional Rust CLI handling of rayon's process-global pool.
4. **M5 follow-up scope** — **Wire it in Wave 5**. Real per-contig
   MD5 verification with a typed `FastaContigMd5Mismatch` error.
5. **Wave ordering** — **Accept 1→2→3→4→5**. No parallel waves.
6. **PR/commit shape** — **Per-task commits behind one plan per
   wave**. Same shape as the original cohort CLI slice (~3 commits
   per wave).
7. **Drop candidates** — **Keep everything**. Mi6 (schema gate)
   and Mi13 (file rename) both stay; their cost is small.
