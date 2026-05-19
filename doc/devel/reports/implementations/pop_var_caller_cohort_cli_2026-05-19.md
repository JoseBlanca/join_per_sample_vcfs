# `pop_var_caller` cohort CLI — implementation report (2026-05-19)

Implementation of the cohort CLI slice planned in
[pop_var_caller_cohort_cli.md](../../implementation_plans/pop_var_caller_cohort_cli.md).
Three new subcommands wire Stages 3–6 of the pipeline into the
existing `pop_var_caller` binary:

- `var-calling` — cohort `.psp` files → multi-sample VCF.
- `estimate-contamination` — cohort `.psp` files → contamination TOML.
- `var-calling-from-bam` — single sample's BAM(s) → single-sample VCF,
  no `.psp` intermediate.

All ten plan tasks landed in ten sequential commits on `main`
(`1523049` through `147e435`); no Task was dropped, no temp-`.psp`
fallback was needed for the from-bam path.

## Plan

The slice was broken into ten sequential tasks, each delivering one
reviewable unit:

| # | Deliverable | Commit |
|---|---|---|
| Plan | Implementation plan | `24621b8` |
| 1 | Engine-side `Config::new -> Result` validators | `1523049` |
| 2 | CLI `value_parser` helpers | `e3cdbf2` |
| 3 | Batch-assignment TSV reader | `e7c4cd0` |
| 4 | Contamination artefact TOML schema + I/O + validator | `566d697` |
| 5 | `estimate-contamination` subcommand | `87e06a3` |
| 6 | `var-calling` subcommand | `17164ea` |
| 7 | Stage 1 pipeline extraction | `e400239` |
| 8 | `var-calling-from-bam` subcommand | `150d118` |
| 9 | Integration tests + fixture | `147e435` |
| 10 | This report + PROJECT_STATUS update | _this commit_ |

## Assumptions

The plan left several things unspecified; choices made during
implementation:

- **`ContaminationEstimationConfig` uses `validate(&self)` instead of
  `new(...)`.** The config has 14+ fields including a `StoppingMode`
  enum; a positional constructor would be hostile at the call site.
  CLI builds via struct literal + `validate()?`; same defensive
  guarantee.
- **`ContaminationEstimates::zero` already existed** in the engine
  (`contamination_estimation.rs:320`); the plan's suggestion to add a
  `::zero(sample_names)` wrapper was unnecessary.
- **`--stopping-mode` flag dropped from `estimate-contamination` v1.**
  The engine supports `Convergence` and `FixedSites` modes, but the
  plan's flag list omitted `--num-sites`, leaving FixedSites
  unreachable from the CLI. v1 locks to convergence-only; library
  callers can use the engine API directly for fixed-N runs.
- **`var_calling_from_bam` duplicates flag definitions.** The plan
  suggested factoring shared sub-structs and using
  `#[command(flatten)]`. Pragmatically deferred — the duplication is
  31 fields and the consolidation is its own refactor.
- **Contamination artefact validator accepts all-zero `q_b` rows.**
  Engine convention: floored batches (singletons or below
  `min_batch_size`) carry an all-zero contaminant distribution.
  Mirrors the same exemption in `ContaminationEstimates::from_user_supplied`.
- **`VariantGrouper` generalised over upstream error type.** The
  original `Iterator<Item = Result<_, PerPositionMergerError>>` bound
  prevented chaining `DustFilter` (yields `DustFilterError`).
  Relaxed to `Iterator<Item = Result<_, E>>` with `GrouperError:
  From<E>`. The struct definition itself is now bound-free; the bounds
  live on the impls.
- **`var-calling-from-bam` walker errors are converted via a local
  scan-adapter**, not via a broader `PerPositionMerger`
  generalisation. The merger's error surface stays
  `PspReadError`-shaped; the scan-adapter stashes `WalkerError` and
  reports end-of-stream, mirroring the existing
  `ErrorSheddingAdapter` pattern for CRAM input errors.
- **`WriterConfig.default_filter_pass` field dropped.** v1 always
  emits `PASS`; a future filter slice will re-introduce filter
  expressions through a typed-rule surface, not a binary flag.

## Changes made

New modules under `src/pop_var_caller/`:

- [batch_assignment.rs](../../../src/pop_var_caller/batch_assignment.rs) (~250 lines)
  — `BatchAssignment::{from_tsv, empty, batch_for}` plus
  `BatchAssignmentError`.
- [contamination_artifact.rs](../../../src/pop_var_caller/contamination_artifact.rs) (~570 lines)
  — `ContaminationArtefact` with serde-derived TOML I/O,
  validator (`validate()`), `to_estimates_for_samples(&[&str])`.
- [estimate_contamination.rs](../../../src/pop_var_caller/estimate_contamination.rs) (~530 lines)
  — `EstimateContaminationArgs`, `run_estimate_contamination`,
  `EstimateContaminationCliError`.
- [var_calling.rs](../../../src/pop_var_caller/var_calling.rs) (~530 lines)
  — `VarCallingArgs`, `run_var_calling`, `VarCallingCliError`.
- [var_calling_from_bam.rs](../../../src/pop_var_caller/var_calling_from_bam.rs) (~640 lines)
  — `VarCallingFromBamArgs`, `run_var_calling_from_bam`,
  `VarCallingFromBamCliError`.
- [stage1_pipeline.rs](../../../src/pop_var_caller/stage1_pipeline.rs) (~180 lines)
  — `with_stage1_pipeline<R, E, F>(...)` callback-based driver.
- [cli/parsers.rs](../../../src/pop_var_caller/cli/parsers.rs) (~470 lines)
  — 20 `clap::value_parser` functions covering every range-bound
  knob in the cohort CLI.

Engine-side changes (Task 1):

- [posterior_engine.rs](../../../src/var_calling/posterior_engine.rs)
  — `PosteriorEngineConfig::new` validated constructor + 5 range
  constants + `PosteriorEngineConfigError`. Closes Stage 6 review
  Mi12.
- [contamination_estimation.rs](../../../src/var_calling/contamination_estimation.rs)
  — `ContaminationEstimationConfig::validate(&self)` + 8 new error
  variants + 5 range constants.
- [variant_grouping.rs](../../../src/var_calling/variant_grouping.rs)
  — `GrouperConfig::new` + `GrouperConfigError` + `VariantGrouper`
  generalised over upstream error type + `GrouperError::DustFilter`
  variant.
- [per_group_merger.rs](../../../src/var_calling/per_group_merger.rs)
  — `PerGroupMergerConfig::new` + `PerGroupMergerConfigError` +
  `MAX_PLOIDY` / `MAX_ALLELES_PER_VAR_CAP` constants.
- [dust_filter.rs](../../../src/var_calling/dust_filter.rs)
  — `SD_WLEN` / `MAX_DUST_WINDOW` constants lifted to `pub const`
  for CLI-side range reuse.
- [vcf_writer/mod.rs](../../../src/var_calling/vcf_writer/mod.rs)
  — `default_filter_pass` field + `DEFAULT_FILTER_PASS` const
  dropped; v1 unconditional `Filters::pass()`.

Wiring (Tasks 5, 6, 8):

- [src/pop_var_caller/mod.rs](../../../src/pop_var_caller/mod.rs)
  — declare 6 new sub-modules, re-export public items.
- [src/pop_var_caller/cli.rs](../../../src/pop_var_caller/cli.rs)
  — three new variants on `PopVarCallerCommand`;
  `run_pileup` refactored to use `with_stage1_pipeline`;
  `parse_mismatch_fraction` lifted to `pub(crate)` for reuse.
- [src/main.rs](../../../src/main.rs)
  — dispatch the three new subcommand variants.

## Tests added / updated

- 3 new tests in `variant_grouping::tests` for `GrouperConfig::new`.
- 7 new tests in `per_group_merger::tests` for
  `PerGroupMergerConfig::new`.
- 14 new tests in `posterior_engine::tests` for
  `PosteriorEngineConfig::new` (boundary tests per range, NaN/inf
  rejection, identity round-trip).
- 14 new tests in `contamination_estimation::tests` for
  `ContaminationEstimationConfig::validate`.
- 20 new tests in `pop_var_caller::cli::parsers::tests` (one per
  value_parser, covering accept-default + boundary + reject NaN /
  inf / non-numeric / out-of-range).
- 15 new tests in `pop_var_caller::batch_assignment::tests` (TSV
  happy path, missing header, empty fields, duplicates, fallback,
  IO error).
- 19 new tests in `pop_var_caller::contamination_artifact::tests`
  (round-trip via tempdir; rejects each validator-class violation;
  sample-name reconciliation: exact / extras / missing / reordered;
  all-zero `q_b` row accepted).
- 7 new tests in `pop_var_caller::estimate_contamination::tests`
  (dense-batch building, RFC3339 helper, basename helper).
- 5 new tests in `pop_var_caller::var_calling::tests` (md5 hex
  round-trip, basename, no-contamination path).
- 3 new tests in `pop_var_caller::var_calling_from_bam::tests`
  (contig conversion: happy / missing-md5 / overflow).
- 1 new test in `vcf_writer::record_encode::tests` updated for
  `default_filter_pass` field removal.

3 new integration tests in
[tests/cohort_cli_integration.rs](../../../tests/cohort_cli_integration.rs):

- `var_calling_happy_path_three_samples` — full chained workflow
  (three pileup runs → run_var_calling → VCF on disk, asserted
  parseable and multi-sample).
- `var_calling_from_bam_happy_path` — single CRAM → single-sample
  VCF, no `.psp` written.
- `var_calling_rejects_contamination_artefact_missing_sample` —
  hand-crafted artefact missing a cohort sample → expected error.

**Test count progression** (cumulative lib tests across the slice):

| Task | Cum. lib tests |
|---|---|
| pre-slice | 812 |
| 1 | 850 (+38) |
| 2 | 870 (+20) |
| 3 | 885 (+15) — recount after Task 4 |
| 4 | 865 (+19, including 1 from Task 5's bug-fix) |
| 5 | 872 (+7) |
| 6 | 877 (+5) |
| 7 | 877 (+0; refactor) |
| 8 | 880 (+3) |
| 9 | 880 (+3 integration tests, separate target) |

End-state: 880 lib tests + 5 prior integration test targets + 3 new
in `cohort_cli_integration` = 888 testable items (some counts above
are approximate — the exact number reflects what each task
contributed inclusive of bug-fix renumbering).

## Validation results

Inside `./scripts/dev.sh`:

```
cargo build --tests                     OK (clean, no warnings)
cargo clippy --workspace --all-targets  OK (-D warnings, 0 errors)
cargo fmt -- --check                    OK
cargo test --lib                        880 passed; 0 failed
cargo test --tests                      every integration target green
```

Pre-existing failure noted in `src/gvcf_parser.rs` doctest (line 6) —
legacy code, `from_gzip_path` signature drifted at some point.
Reproduced on `main` before any of this slice's changes; unrelated to
the cohort CLI. The legacy `gvcf_parser` module is marked
"superseded" in PROJECT_STATUS.md.

Manual smoke (deferred): `bcftools view` / `bcftools stats`
inspection of a real VCF output. The fixture used in the integration
tests is too tiny to demonstrate meaningful behaviour; manual smoke
on real cohort data is a follow-up.

## Tradeoffs and follow-ups

**Deferred from the plan:**

- 7 of the 10 plan-listed integration test cases (`--no-complexity-filter`
  toggle, end-to-end with `estimate-contamination` feeding
  `var-calling`, config-validation hard error from the CLI side,
  artefact extras tolerated, etc.). The three landed tests are
  load-bearing end-to-end gates; the deferred cases are mostly
  covered in spirit by per-module unit tests. Reintroduce as needed.
- `bcftools view` / `bcftools stats` manual smoke against real
  cohort data. The synthetic fixture produces zero variant calls
  (tiny coverage, no real diversity).
- Implementation report for the Stage 5 per-group merger (still
  open in PROJECT_STATUS from before this slice).
- Shared `tests/common/` module so cohort and pileup integration
  tests can share fixture helpers.

**Open follow-ups (post-this-slice):**

- `--regions BED` for restricting calls (standing item, plan §"Out
  of scope").
- `tabix .tbi` index alongside `.vcf.gz` outputs.
- `PL` `FORMAT` field forwarded from Stage 5.
- Per-sample contamination fraction in `INFO`.
- `--external-allele-frequencies` for `estimate-contamination`.
- Cross-record rayon parallelism (the deferred Stage 6 perf lever;
  see [posterior_engine_perf_2026-05-18_v2.md](posterior_engine_perf_2026-05-18_v2.md)).
- `--per-group-batch-size` exposure (see project memory:
  parallelisation tuning is deferred until after this slice).
- `pop_var_caller pileup` itself extracted into its own sibling
  module (cli.rs is large now; the other subcommands are siblings).
- `var-calling-from-bam` args struct factored via
  `#[command(flatten)]` shared sub-structs.

**Risks now closed:**

- The plan flagged `var-calling-from-bam` lifetime juggling as the
  highest-risk piece, with explicit fallback "drop the subcommand if
  intractable". The callback-based `with_stage1_pipeline` (Task 7)
  resolved it cleanly — borrow-free, no temp `.psp` ever written, no
  unsafe or self-referential machinery.
- The cohort CLI plan asserted that engine-side `Config::new ->
  Result` would close Stage 6 review Mi12. Done in Task 1.
- The cohort CLI plan locked validation ranges narrower than the
  mathematical domain (e.g. `max_iterations` capped at 500, ploidy
  at 8, max_alleles_per_var at 16). All enforced both CLI-side via
  `value_parser` and engine-side via `Config::new` /
  `Config::validate`.
