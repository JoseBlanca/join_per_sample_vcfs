# Code Review: cohort CLI slice (`pop_var_caller` Stages 3–6 wiring)
**Date:** 2026-05-19
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** PR-shaped diff `03dc622..5bed142` on `main` — 11 commits, ~7600 inserted lines, adding three new subcommands (`var-calling`, `estimate-contamination`, `var-calling-from-bam`) plus shared infrastructure (`stage1_pipeline.rs`, `contamination_artifact.rs`, `batch_assignment.rs`, `cli/parsers.rs`) and engine-side `Config::new -> Result` validators.
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** PR-shaped diff (slice already merged on `main`).
- **Reviewed against:** commit range `03dc622..5bed142` (most-recent first: `5bed142` docs report; `147e435` integration tests; `150d118` `var-calling-from-bam`; `e400239` Stage 1 helper extract; `17164ea` `var-calling`; `87e06a3` `estimate-contamination`; `566d697` contamination artefact; `e7c4cd0` batch-assignment TSV; `e3cdbf2` CLI parsers; `1523049` engine-side validators).
- **In-scope files:**
  - New modules under `src/pop_var_caller/`:
    [var_calling.rs](../../../src/pop_var_caller/var_calling.rs),
    [var_calling_from_bam.rs](../../../src/pop_var_caller/var_calling_from_bam.rs),
    [estimate_contamination.rs](../../../src/pop_var_caller/estimate_contamination.rs),
    [contamination_artifact.rs](../../../src/pop_var_caller/contamination_artifact.rs),
    [batch_assignment.rs](../../../src/pop_var_caller/batch_assignment.rs),
    [stage1_pipeline.rs](../../../src/pop_var_caller/stage1_pipeline.rs),
    [cli/parsers.rs](../../../src/pop_var_caller/cli/parsers.rs).
  - Modified: [cli.rs](../../../src/pop_var_caller/cli.rs) (run_pileup refactor + new subcommand variants), [mod.rs](../../../src/pop_var_caller/mod.rs), [main.rs](../../../src/main.rs) (dispatch).
  - Engine-side surface additions only:
    [variant_grouping.rs](../../../src/var_calling/variant_grouping.rs),
    [contamination_estimation.rs](../../../src/var_calling/contamination_estimation.rs),
    [per_group_merger.rs](../../../src/var_calling/per_group_merger.rs),
    [posterior_engine.rs](../../../src/var_calling/posterior_engine.rs),
    [dust_filter.rs](../../../src/var_calling/dust_filter.rs),
    [vcf_writer/mod.rs](../../../src/var_calling/vcf_writer/mod.rs),
    [vcf_writer/record_encode.rs](../../../src/var_calling/vcf_writer/record_encode.rs).
  - Tests: [cohort_cli_integration.rs](../../../tests/cohort_cli_integration.rs).
- **Deliberately out of scope:**
  - Algorithmic cores under `src/var_calling/` (per-position merger internals, grouper algorithm, per-group merger algorithm, posterior engine algorithm, vcf_writer/* internals beyond the surface additions named above) — pre-existing.
  - The `ebbbab4 docs: architecture diagram` commit — unrelated docs change.
  - Pre-existing failing `src/gvcf_parser.rs` doctest (legacy code, marked superseded).
  - Pre-existing `cargo doc --no-deps` errors at `pileup_to_psp.rs:5-6` and `posterior_engine.rs:725/732/753` — out of scope per task brief.
- **Categories dispatched (10):**
  - `reliability` — always; this is a CLI surface.
  - `errors` — always.
  - `naming` — always.
  - `defaults` — public-API + every clap `default_value_t` knob.
  - `idiomatic` — always.
  - `refactor_safety` — always; the slice closes a "make compiler flag refactors" gap (Mi12).
  - `unsafe_concurrency` — `Arc<dyn _>`, `Rc<RefCell<...>>`, and three `rayon::ThreadPoolBuilder::build_global` sites trigger.
  - `smells` — always.
  - `tooling` — slice ships under `Cargo.toml`.
  - `extras` — slice contains parsers/validators, produces a stable TOML output, accepts untrusted user input, and is PR-shaped (diff-vs-intent applies).

## 2. Verdict

**Request-changes.**

Zero Blockers. Fourteen Major findings — five are silent-failure modes (CRAM-error stash dropped on closure-Err, reference cross-check is basename-only, walker-error path untested end-to-end, post-construction mutation bypasses validation, silent `unwrap_or_default()` on missing MD5), three are public-API hygiene (`WriterConfig` non-exhaustive missing, doc/value mismatch on a user-visible tolerance constant, post-construction config mutation), four are maintainability concerns large enough to deserve fixing in this slice (cohort pipeline duplication, ~95-line closure, three-way `*Args` duplication, helper duplication across 3–4 files), and two are tooling (`cargo doc --no-deps` fails on five new broken intra-doc links; rayon-pool process-singleton hazard with no in-process coordination).

The slice itself is genuinely good code (consistent `#[non_exhaustive]` discipline, exhaustive struct literals, well-shaped error enums, the callback-style `with_stage1_pipeline` cleanly resolves the borrow chain the plan flagged as the highest risk). The findings are mostly *seams that became silent* — places where a future refactor could compile cleanly and produce wrong results without a test catching it.

## 3. Execution status

Commands run (all inside `./scripts/dev.sh`):

| Command | Exit | Result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --workspace --all-targets -- -D warnings` | 0 | clean |
| `cargo test --lib` | 0 | `test result: ok. 880 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 39.30s` |
| `cargo test --tests` | 0 | every integration target green |
| `cargo doc --no-deps` | **101** | 5 new broken intra-doc-link errors + 1 new warning introduced by this slice (see M14). Pre-existing failures in `pileup_to_psp.rs` and `posterior_engine.rs` are noted under "Out of scope observations". |

Findings labeled "Needs verification": 1 (M9, CRLF body rows in `BatchAssignment` — the reviewer reasoned about behaviour from `&str::lines()` semantics but did not run a CRLF fixture).

## 4. Open questions and assumptions

1. **Is `Q_B_SIMPLEX_TOLERANCE = 1e-9` intentionally **stricter** than the engine's `1e-6`, or was `1e-9` a typo for `1e-5`?** Affects **M7**. The doc says "loosened … with slack on top of the engine bound" which contradicts a `1e-9` value.
2. **Should `--threads` on two cohort subcommands in one process succeed, error, or no-op silently?** Affects **M13**. The current behaviour (second call always fails) is a defensible policy but is not exercised by any test and is not surfaced as part of the documented contract.
3. **Is the cohort-CLI's `--reference` FASTA expected to be content-verified against `.psp` per-contig MD5s, or is basename-equality the documented limit?** Affects **M5**. The plan §"Reference cross-check tolerance" says MD5 *is* checked; the code does not. Either the doc or the code needs to move.
4. **Was the `--stopping-mode` flag's omission a v1 deferral or a permanent decision?** Affects no severity-tier finding but the engine retains `FixedSites` machinery, half-reachable. Out of scope for this slice; flag for the plan-vs-impl reconciliation pass.
5. **Should `WriterConfig` field changes be silent-breaking for downstream crate consumers?** Affects **M8**. The slice removed `default_filter_pass` without `#[non_exhaustive]`; if `merge_per_sample_vcfs` is ever published, this is the moment to lock the surface.
6. **The "all-zero `q_b` row signals floored batch" exemption is asymmetric**: the artefact validator accepts it; `to_estimates_for_samples` round-trips it; but no test pins the contract across that boundary. Affects **M11** (challenge test).

## 5. Top 3 priorities

1. **M1 / B-class behaviour without enough evidence to file as Blocker — `stage1_pipeline.rs:167` silently drops a stashed CRAM-input error whenever the closure returns `Err`.** The whole point of the error-shedding adapter is to surface "walker exhausted cleanly but the source stream had an error"; today that signal is lost on the closure-Err path and the user is shown the downstream symptom instead of the root cause. Fix is ~10 lines (prefer the stash over the closure error). See **M1**.
2. **M5 — `var_calling.rs:335-348` reference cross-check is basename-only; the plan's claimed FASTA MD5 enforcement is unimplemented.** A user supplying the wrong `--reference` whose basename matches gets silently-wrong variant calls. Fix is either (a) tighten the docs to match the basename-only reality, or (b) wire per-contig MD5 against the FASTA. The plan promised (b). See **M5**.
3. **M14 — `cargo doc --no-deps` is failing on the slice (5 new broken intra-doc links + 1 warning).** `broken_intra_doc_links = "deny"` is set in `Cargo.toml:20`, so every individual error fails the build. All five are mechanical doc fixes (~5 lines). See **M14**.

## 6. Findings

### Blocker

(none)

### Major

#### M1: [src/pop_var_caller/stage1_pipeline.rs:119-167](../../../src/pop_var_caller/stage1_pipeline.rs#L119-L167) — Stashed CRAM-input error silently dropped when the closure returns Err
- **Categories:** errors, reliability (challenge test)
- **Confidence:** High
- **Problem:** Both branches of `with_stage1_pipeline` snapshot `stashed_upstream_error = error_handle.take()` into a local *before* `let r = result?;` on line 167. When `result` is `Err(_)`, the `?` returns the closure's error and the function exits, dropping `stashed_upstream_error` along with every other local. The `Stage1Outputs { stashed_upstream_error, .. }` struct is never constructed on the error path. So if the closure failed for a downstream reason (writer/posterior) but the **root cause** was an upstream CRAM/BAQ error that the `ErrorSheddingAdapter` translated into end-of-stream (the entire point of the adapter), the caller never sees the CRAM error — they see the downstream symptom.
- **Why it matters:** This is exactly the case the stash exists to surface (the helper's own doc-comment at L62–67 says "walker exhausted cleanly but the source stream had an error"). The user is told "vcf write failed" when the real story is "CRAM record N corrupt". Held at Major (not Blocker) because it requires `f` to return `Err` *and* a stash to be present; not 100%-reproducible.
- **Suggested fix:** Prefer the upstream error over the closure error when both are present:
  ```rust
  // Replace `let r = result?;` (current L167) with:
  let r = match result {
      Ok(v) => v,
      Err(e) => {
          if let Some(upstream) = stashed_upstream_error {
              // Prefer the root cause; conversion exists via E: From<PileupCliError>.
              return Err(PileupCliError::CramInput(upstream).into());
          }
          return Err(e);
      }
  };
  ```
  Add `with_stage1_pipeline_stashes_cram_input_error_when_walker_exhausts_cleanly` and a partner `with_stage1_pipeline_prefers_upstream_error_over_closure_error` test under `#[cfg(test)] mod tests` in `stage1_pipeline.rs`.

#### M2: [src/pop_var_caller/var_calling_from_bam.rs:556-559](../../../src/pop_var_caller/var_calling_from_bam.rs#L556-L559) — Walker-error stash path is end-to-end untested
- **Categories:** reliability
- **Confidence:** High
- **Problem:** The from-bam orchestrator stashes walker errors via `Rc<RefCell<Option<WalkerError>>>` inside an `Iterator::scan`, drives the whole posterior chain to completion, and after `writer.finish()` best-effort deletes the *final* output file. The only integration coverage is the happy path (`var_calling_from_bam_happy_path`); no test demonstrates that (a) a walker mid-stream error surfaces as `Walker(_)` instead of silently terminating the stream, (b) the renamed VCF is removed, (c) the surrounding error chain is unbroken. The scan-adapter is exactly the kind of error-shedding seam where a refactor that drops `walker_error_stash.borrow_mut().take()` produces *wrong, silent success*.
- **Why it matters:** Assumption 5 of the impl report. With no test, "walker error stashed but pipeline ran to completion" becomes the silent success path and a future refactor lands undetected.
- **Suggested fix:** Add integration test `var_calling_from_bam_surfaces_walker_error_and_removes_output`. Synthesize reads dense enough to trip `max_active_reads = 1`; assert `Err(Walker(_))` and that the output VCF does not exist after the call. Pair with a unit test on the scan-adapter pattern in isolation (`walker_error_shim_stashes_and_ends_stream`) so the contract is documented in one place.

#### M3: [src/pop_var_caller/var_calling_from_bam.rs:525-540](../../../src/pop_var_caller/var_calling_from_bam.rs#L525-L540) — Silent `unwrap_or_default()` on missing MD5 and silent truncating cast on contig length in the DUST branch
- **Categories:** errors, reliability
- **Confidence:** High
- **Problem:** Line 477 already calls `contigs_to_parsed(ctx.contigs)?`, which hard-errors with `MissingMd5` / `ContigLengthOverflow`. Five lines later, the DUST branch *re-builds* a parallel `Vec<ParsedChromosome>` directly from `ctx.contigs.entries`, this time with `md5: e.md5.map(format_md5_hex).unwrap_or_default()` and `length: e.length as u32` (silent truncating cast). The first call's invariants are dead-defence against this site today, but the duplication is begging to be silently broken — and the function already owns dedicated typed variants for both failure modes.
- **Why it matters:** Exactly the "silent-fallback / invariant-violation site" pattern the user's `feedback_no_logs_use_errors` memory targets. Promoting to typed errors is the project rule; this site bypasses it.
- **Suggested fix:** Reuse the validated vector rather than rebuilding it:
  ```rust
  // earlier (already exists):
  let chromosomes = contigs_to_parsed(ctx.contigs)?;
  // …
  let merger = PerPositionMerger::new(
      vec![walker_iter],
      vec![sample_name_owned.clone()],
      chromosomes.clone(),
  )?;
  // DUST branch:
  let dust = DustFilter::new(merger, &*fetcher, chromosomes, dust_cfg);
  ```

#### M4: [src/pop_var_caller/var_calling.rs:367](../../../src/pop_var_caller/var_calling.rs#L367) — Post-construction mutation of `posterior_cfg.contamination` bypasses any future validation `PosteriorEngineConfig::new` may enforce on that field
- **Categories:** refactor_safety, defaults (Nit)
- **Confidence:** High
- **Problem:** `PosteriorEngineConfig::new(...)` sets `contamination: None` internally and contains no validation for it. The CLI then mutates the public field directly: `posterior_cfg.contamination = load_contamination(...)?`. A future commit that adds a `Some(_)` invariant inside `new` (e.g. "cohort size of estimates must match cohort size of posterior weights") would still hold for `new`'s internal `None`, then this assignment silently re-introduces an unchecked `Some(_)`. The `var_calling_from_bam.rs:455` site inherits the same hazard (currently sets `None`, which is fine today but the seam is identical).
- **Why it matters:** Validators that look airtight from the engine side are silently bypassed at the CLI seam. The hazard does not show up in the diff (the field is currently unvalidated) but is set up to be silently violated by a future change.
- **Suggested fix:** Either add a validated setter on `PosteriorEngineConfig` (`set_contamination(&mut self, ...) -> Result<(), _>`), or extend `new` to take `contamination: Option<ContaminationEstimates>` as a parameter so the construction path is single-entry. The struct-literal-plus-validate pattern used by `ContaminationEstimationConfig` works here too:
  ```rust
  let posterior_cfg = PosteriorEngineConfig {
      convergence_threshold: args.em_convergence_threshold,
      max_iterations: args.em_max_iterations,
      // ... every field named
      contamination: load_contamination(...)?,
      ..PosteriorEngineConfig::with_project_defaults()
  };
  posterior_cfg.validate()?;
  ```

#### M5: [src/pop_var_caller/var_calling.rs:335-348](../../../src/pop_var_caller/var_calling.rs#L335-L348) — Reference cross-check is basename-only; per-contig MD5 from supplied FASTA is never compared
- **Categories:** errors, extras
- **Confidence:** High
- **Problem:** The cross-check compares `reader.header().reference` against `basename(&args.reference)`. The plan §"Reference cross-check tolerance" claims md5 is also checked via `check_chromosome_agreement`. Verified against [per_position_merger.rs:354-400](../../../src/var_calling/per_position_merger.rs#L354-L400): `check_chromosome_agreement` only compares chromosomes *between readers*, not against the supplied FASTA. `SyncRefFetcher::new` at [ref_fetcher.rs:44-54](../../../src/per_sample_pileup/ref_fetcher.rs#L44-L54) also does no MD5 check. A user supplying a FASTA whose contig names match but whose bytes differ (or whose lengths differ within the FAI) gets silently-wrong variant calls. Same gap in `estimate_contamination.rs:321-330`.
- **Why it matters:** Variant calls would be silently wrong (BAQ/DUST/posterior all key off the FASTA bases). The error message at `var_calling.rs:278-286` reads `"reference mismatch — header has `{psp_ref}`, CLI passed `{supplied_ref}` (basename comparison)"`; the `(basename comparison)` parenthetical is the right concession, but the plan promised more.
- **Suggested fix:** Either (a) tighten the module-level doc-comment (e.g. `var_calling.rs:304`) to state that MD5 enforcement is .psp-vs-.psp only, OR (b) actually verify per-contig MD5 against the supplied FASTA. Option (b) is the long-term answer (load `.fai`, hash each contig, or compute via `noodles_fasta`'s MD5 helper). Option (a) is the pre-merge-cheap answer. **Either way, decide which.** See open question 3.

#### M6: [src/pop_var_caller/var_calling.rs:419-426](../../../src/pop_var_caller/var_calling.rs#L419-L426) — Iterator error leaves `<output>.tmp` orphaned; no cleanup parity with `run_pileup` / `run_var_calling_from_bam`
- **Categories:** extras, defaults (cross-cat), idiomatic (cross-cat), smells (cross-cat)
- **Confidence:** High
- **Problem:** `CohortVcfWriter::new` opens `<output>.tmp` and only `finish()` renames. In `run_var_calling`, any error in the iteration loop (lines 421-426) bubbles via `?` before `writer.finish()`, leaving `<output>.tmp` on disk indefinitely. By contrast `run_pileup` (`cli.rs:353-356`) and `run_var_calling_from_bam` (`var_calling_from_bam.rs:556-560`, `:594-600`) both have explicit best-effort `remove_file` cleanup paths. The plan §"Subcommand: var-calling" step 8 says the tmp-rename convention matches `run_pileup`, but the cleanup half wasn't carried over.
- **Why it matters:** Inconsistent failure-cleanup discipline across the three cohort subcommands. Repeated failed runs accumulate `.tmp` files; users debugging a failed run also have to triage stale tmp paths.
- **Suggested fix:** Wrap the writer loop so failures call `let _ = std::fs::remove_file(tmp_path)` before propagating. Cleaner: give `CohortVcfWriter` a `Drop` impl that removes the tmp file when `finish()` was not called. Either pattern restores parity with `run_pileup`.

#### M7: [src/pop_var_caller/contamination_artifact.rs:73-77](../../../src/pop_var_caller/contamination_artifact.rs#L73-L77) — Doc contradicts the constant value for `Q_B_SIMPLEX_TOLERANCE`
- **Categories:** defaults
- **Confidence:** High
- **Problem:** Constant: `pub const Q_B_SIMPLEX_TOLERANCE: f64 = 1e-9;`. Doc: "engine-side tolerance (`1e-6`) loosened by an additional safety factor — round-trip via `toml::Value::Float` (f64) can shave a few ULP off, so we keep some slack on top of the engine bound." `1e-9` is three orders of magnitude **stricter** than `1e-6`, not looser. The constant is also embedded in the user-visible `BatchProbabilitiesDoNotSum` error message at `:394-398`, so a future "why does my artefact reject?" debug will hit the doc-vs-value mismatch.
- **Why it matters:** Ambiguous default. Either the doc or the value is wrong.
- **Suggested fix:** Pick one — see open question 1. If the validator is intentionally stricter than the engine, rewrite the doc to say so. If `1e-6 + slack` was intended, change the value to ~`1e-5` and keep the doc.

#### M8: [src/var_calling/vcf_writer/mod.rs:67-90](../../../src/var_calling/vcf_writer/mod.rs#L67-L90) — `WriterConfig` lacks `#[non_exhaustive]`; the `default_filter_pass` removal is an unannounced breaking change
- **Categories:** defaults
- **Confidence:** High
- **Problem:** Two `pub` fields, no `#[non_exhaustive]`. Every in-scope call site uses struct literals (`var_calling.rs:384`, `var_calling_from_bam.rs:487`, `tests/cohort_vcf_writer_integration.rs:161/228/261/300`), bypassing the existing `WriterConfig::new(output)` constructor that is the only place `DEFAULT_EMIT_GP` lives. Removing the `default_filter_pass` field is a breaking change for any downstream crate using struct literal form (the rewritten in-tree tests show exactly that shape). Without `#[non_exhaustive]`, the next field addition is silently breaking the same way.
- **Why it matters:** The "single validated/documented constructor that announces every default" rule. Today the constructor pins `emit_gp = DEFAULT_EMIT_GP`, but the CLI ignores the constructor and uses struct literals — the documented default applies only to library callers who happen to use `new()`.
- **Suggested fix:** Add `#[non_exhaustive]` and a builder so future field additions cannot silently break callers and the default for `emit_gp` lives in exactly one place:
  ```rust
  #[derive(Debug, Clone)]
  #[non_exhaustive]
  pub struct WriterConfig {
      pub output: PathBuf,
      pub emit_gp: bool,
  }
  impl WriterConfig {
      pub fn new(output: PathBuf) -> Self { Self { output, emit_gp: DEFAULT_EMIT_GP } }
      pub fn with_emit_gp(mut self, emit_gp: bool) -> Self { self.emit_gp = emit_gp; self }
  }
  ```
  CLI uses `WriterConfig::new(output).with_emit_gp(args.emit_gp)`. Integration tests pivot to the same shape.

#### M9: [src/pop_var_caller/var_calling_from_bam.rs:503](../../../src/pop_var_caller/var_calling_from_bam.rs#L503) — Inline `Iterator::scan` reinvents `ErrorSheddingAdapter`; the doc-comment names a type that doesn't exist
- **Categories:** idiomatic, naming (broken doc link), tooling (M14e)
- **Confidence:** High
- **Problem:** Lines 503–513 build an inline error-shedding adapter from `Rc<RefCell<Option<WalkerError>>>` + `Iterator::scan` to convert `Result<PileupRecord, WalkerError>` into `Result<PileupRecord, PspReadError>`. The exact pattern is already factored as `ErrorSheddingAdapter` in [error_bridge.rs](../../../src/pop_var_caller/cli/error_bridge.rs) — but parameterised only over `Item = PreparedRead, E = CramInputError`. The inline copy doesn't expose a `take()` API (the stash is read via `borrow_mut().take()` four steps away at L557), and the inline shape lies about the iterator's `Item` type (wrapping in `Result<_, PspReadError>` solely to satisfy `PerPositionMerger`'s input bound — the `Err` variant is unreachable). The module doc at L13 even names a `WalkerErrorSheddingAdapter` type that doesn't exist in the crate (one of the broken intra-doc links in M14).
- **Why it matters:** Two copies of a non-obvious self-referential-state pattern drift independently. The stale doc reference suggests the implementer initially planned a named adapter type and never renamed the comment.
- **Suggested fix:** Generalise `ErrorSheddingAdapter` to be parametric over `Item` and `E`, then reuse it. If generalising is too invasive, at minimum factor *this* call site into a dedicated `WalkerErrorSheddingAdapter` in `cli/error_bridge.rs` so the doc-comment becomes truthful. Sketch:
  ```rust
  pub struct ErrorSheddingAdapter<I, T, E> {
      inner: I,
      handle: Rc<RefCell<Option<E>>>,
      done: bool,
      _t: PhantomData<T>,
  }
  // walker site:
  let adapter = ErrorSheddingAdapter::<_, PileupRecord, WalkerError>::new(ctx.walker);
  let handle = adapter.error_handle();
  let walker_iter: Box<dyn Iterator<Item = Result<PileupRecord, PspReadError>>> =
      Box::new(adapter.map(Ok));
  if let Some(e) = handle.take() { return Err(VarCallingFromBamCliError::Walker(e)); }
  ```

#### M10: [src/pop_var_caller/var_calling_from_bam.rs:83-352](../../../src/pop_var_caller/var_calling_from_bam.rs#L83-L352) — `VarCallingFromBamArgs` duplicates ~30 fields of `VarCallingArgs` + `PileupArgs` verbatim
- **Categories:** smells
- **Confidence:** High
- **Problem:** Third copy of every Stage 1 knob (already on `PileupArgs`) and every Stage 3–6 knob (already on `VarCallingArgs`) — same names, defaults, `value_parser`s, `help_heading`s, `hide_short_help`. Plan §"CLI surface" and impl report flag this as a deferred consolidation via `#[command(flatten)]`. With ~25 knobs simultaneously fanning out, default-value drift is the standard failure mode — a CLI knob whose docstring is updated in one place and not in two others lies in `--help`.
- **Why it matters:** Every future tweak to a knob (rename, range change, help-text edit, parser swap) has to be made in three places with no compile guard. The "single source of truth for `-h` text" the plan justifies the CLI side with is already broken on day one.
- **Suggested fix:** Extract `Stage1Args` and `CohortPipelineArgs` sub-structs and use `#[command(flatten)]`. Concrete shape in [smells.md](../../../tmp/review_2026-05-19_cohort-cli/smells.md) — call-site fan-out (`args.min_mapq` → `args.stage1.min_mapq`) is mechanical and the `*_config_from_args` helpers only look one field deeper.

#### M11: [src/pop_var_caller/var_calling.rs:391-426](../../../src/pop_var_caller/var_calling.rs#L391-L426) and [src/pop_var_caller/var_calling_from_bam.rs:515-554](../../../src/pop_var_caller/var_calling_from_bam.rs#L515-L554) — Cohort pipeline wiring (merger → DUST → grouper → per_group → posterior → writer) duplicated near-verbatim
- **Categories:** smells
- **Confidence:** High
- **Problem:** Both files build Stages 3–6 with the same control flow: `PerPositionMerger::new` → conditional `DustFilter`/passthrough wrapped in `Box<dyn Iterator<Item = Result<_, GrouperError>>>` → `VariantGrouper::with_config` → `PerGroupMerger::with_config` → `PosteriorEngine::with_config` → drain loop into `CohortVcfWriter`. ~35 lines each, including the Box-erasure trick and the DUST-on/DUST-off branching. Differences: merger's input shape (k=N vs k=1) and the `chromosomes` source.
- **Why it matters:** The next person to change one of: per-group batch size, ploidy plumbing, `GrouperError` mapping, or to add another optional stage will edit two assemblies. Bonus risk in the `Box<dyn Iterator>` type which is an easy place to introduce a subtle lifetime/generic mismatch.
- **Suggested fix:** Extract a shared `drive_cohort_pipeline` driver parameterised by `From<GrouperError> + From<PerGroupMergerError> + From<PosteriorEngineError> + From<VcfWriteError> + From<DustFilterError>` so both `VarCallingCliError` and `VarCallingFromBamCliError` can adopt it. Sketch in [smells.md](../../../tmp/review_2026-05-19_cohort-cli/smells.md).

#### M12: [src/pop_var_caller/var_calling_from_bam.rs:474-568](../../../src/pop_var_caller/var_calling_from_bam.rs#L474-L568) — ~95-line closure passed to `with_stage1_pipeline`
- **Categories:** smells, idiomatic (Minor in idiomatic)
- **Confidence:** High
- **Problem:** The `|ctx| -> Result<(), VarCallingFromBamCliError>` closure spans 95 lines doing seven distinct things (contigs conversion, metadata, fetcher setup, walker shim, full DUST/grouper/per-group/posterior assembly, VCF drive, walker-error surfacing). The closure captures six values from the outer scope, which makes the capture set the entire interesting parameter surface; promoting to a named function with explicit parameters costs nothing and makes the data flow legible.
- **Why it matters:** Long closures push their binding shape, captures, and return-type annotation offscreen; the named-function rewrite gives the body a doc-comment and a signature. Pairs naturally with the M11 driver extraction — the inner function shrinks to ~20 lines after both refactors land.
- **Suggested fix:** Extract `fn run_cohort_pipeline_for_single_sample(ctx, reference, output, command_line, emit_gp, no_complexity_filter, dust_cfg, grouper_cfg, per_group_cfg, posterior_cfg) -> Result<(), _>`. The `with_stage1_pipeline` call shrinks to one line.

#### M13: [src/pop_var_caller/var_calling.rs:319](../../../src/pop_var_caller/var_calling.rs#L319), [src/pop_var_caller/var_calling_from_bam.rs:427](../../../src/pop_var_caller/var_calling_from_bam.rs#L427), [src/pop_var_caller/estimate_contamination.rs:305](../../../src/pop_var_caller/estimate_contamination.rs#L305) — Three `rayon::ThreadPoolBuilder::build_global()` call sites with no in-process coordination
- **Categories:** unsafe_concurrency, reliability (cross-cat), defaults (cross-cat)
- **Confidence:** High
- **Problem:** Each subcommand driver calls `build_global` from within a `pub` library function, gated only by `if let Some(n) = args.threads`. `build_global` is process-global one-shot. A library consumer or future test that runs e.g. `run_estimate_contamination` then `run_var_calling` in the same process with `--threads` on the second call gets `RayonAlreadyConfigured` — even though the user did nothing wrong. The asymmetry is also surprising: `--threads N` on the first call works regardless of which subcommand; on the second call always fails. Tests dodge the issue by using `threads: None` exclusively, so the `--threads` path has *no coverage* on any of the three subcommands.
- **Why it matters:** Library-API concurrency hazard. The policy "binary calls it at most once" is asserted in doc-comments but not enforced by code structure.
- **Suggested fix:** Centralise pool configuration via a `OnceLock`-gated idempotent helper:
  ```rust
  static POOL_CONFIGURED: OnceLock<()> = OnceLock::new();
  pub fn configure_rayon_pool(n: Option<usize>) -> Result<(), rayon::ThreadPoolBuildError> {
      let Some(n) = n else { return Ok(()); };
      if POOL_CONFIGURED.get().is_some() { return Ok(()); }
      rayon::ThreadPoolBuilder::new().num_threads(n).build_global()?;
      let _ = POOL_CONFIGURED.set(());
      Ok(())
  }
  ```
  Each driver calls `configure_rayon_pool(args.threads)?`. Decide whether mismatched `n` between sequential calls is a silent no-op or a typed error. See open question 2.

#### M14: cargo doc breakage — five new broken intra-doc links + one warning
- **Categories:** tooling, naming (Minor for L14e), reliability (cross-cat)
- **Confidence:** High
- **Problem:** [Cargo.toml:20](../../../Cargo.toml#L20) sets `broken_intra_doc_links = "deny"`. The slice introduces five new errors and one warning:
  - **M14a** — [batch_assignment.rs:30](../../../src/pop_var_caller/batch_assignment.rs#L30): `Self::batch_for` in module-level doc (`Self` not in scope at module level).
  - **M14b** — [contamination_artifact.rs:56](../../../src/pop_var_caller/contamination_artifact.rs#L56): `Self::to_estimates_for_samples` in module-level doc (same reason).
  - **M14c** — [contamination_artifact.rs:277](../../../src/pop_var_caller/contamination_artifact.rs#L277): `ContaminationEstimateSource::UserSupplied` not imported into the module.
  - **M14d** — [stage1_pipeline.rs:15](../../../src/pop_var_caller/stage1_pipeline.rs#L15): `[feedback_no_silent_intermediates]` parsed as a broken intra-doc link (it's a pointer to a user-memory file, opaque to rustdoc).
  - **M14e** — [var_calling_from_bam.rs:13](../../../src/pop_var_caller/var_calling_from_bam.rs#L13): `WalkerErrorSheddingAdapter` named but the symbol doesn't exist in the crate (the actual mechanism is an inline `.scan()` — see M9).
  - **Nit** — [estimate_contamination.rs:296](../../../src/pop_var_caller/estimate_contamination.rs#L296): redundant explicit link target (`[`ContaminationEstimates`](crate::var_calling::contamination_estimation::ContaminationEstimates)` — shortcut form would render identically). Warning only, but trivial to fix in the same pass.
- **Why it matters:** `cargo doc --no-deps` fails with exit 101. Each error fails the build independently. The previous Stage 6 review's Mi21 was deferred specifically pending a clean doc build; this slice extended the debt.
- **Suggested fix:** Five mechanical doc fixes; full text in [tooling.md](../../../tmp/review_2026-05-19_cohort-cli/tooling.md). Quick references:
  - M14a/b: replace `Self::` with the concrete owning type (`BatchAssignment::batch_for`, `ContaminationArtefact::to_estimates_for_samples`).
  - M14c: use fully-qualified path or import the enum.
  - M14d: drop the bracketed link, describe inline ("the "no silent intermediates" principle").
  - M14e: rewrite to name the inline `.scan()` mechanism truthfully, OR extract a named `WalkerErrorSheddingAdapter` type (paired with M9) to make the doc correct.

### Minor

#### Mi1: [src/pop_var_caller/cli/parsers.rs:173-180](../../../src/pop_var_caller/cli/parsers.rs#L173-L180) — `parse_pseudocount` uses a generic flag name in error messages across four distinct flags
- **Categories:** naming, defaults
- **Confidence:** High
- **Problem:** `parse_pseudocount` is wired to four flags (`--ref-pseudocount`, `--snp-alt-pseudocount`, `--indel-alt-pseudocount`, `--compound-alt-pseudocount`). Body uses `name = "pseudocount"` (line 175). The sibling `parse_contam_pseudocount` (`:238-245`) — *deliberately split out per its in-file comment so the flag name in error messages would be specific* — is then *also* wired to three distinct contamination flags, so it has the same problem. The "one parser, four flags" pattern silently inverts the precedent the project just established. Clap usually prepends the flag name in its own error rendering, so user impact is partially masked; but the sub-agent reviewing this category called out that the explicit per-flag rendering precedent makes the inconsistency the surprise.
- **Why it matters:** The same concept-name appears across four distinct knobs in user-facing text. The `parse_contam_pseudocount` precedent shows the project chose verbose specificity; this site chose generic brevity. Pick one convention.
- **Suggested fix:** Split into per-flag entry points mirroring the contamination side:
  ```rust
  fn parse_pseudocount_with_name(s: &str, name: &str) -> Result<f64, String> {
      parse_f64_with(s, name, |v| 0.0 < v && v <= PSEUDOCOUNT_RANGE_MAX, "(0.0, 1000.0]")
  }
  pub fn parse_ref_pseudocount(s: &str) -> Result<f64, String> {
      parse_pseudocount_with_name(s, "ref-pseudocount")
  }
  // …same for snp-alt, indel-alt, compound-alt, and the three contam variants
  ```
  Demoted from Major to Minor at synthesis: clap typically surfaces the flag name, and the existing `parse_contam_pseudocount` precedent is convention-level, not behaviour-level.

#### Mi2: [src/pop_var_caller/var_calling.rs:355](../../../src/pop_var_caller/var_calling.rs#L355) and [src/pop_var_caller/var_calling_from_bam.rs:444-453](../../../src/pop_var_caller/var_calling_from_bam.rs#L444-L453) — `PosteriorEngineConfig::new` takes eight positional, mostly-`f64` arguments; reorder is not type-checked
- **Categories:** refactor_safety (Major M2 in category), idiomatic
- **Confidence:** High
- **Problem:** Eight positional args, six of them `f64`s with default values inside one another's ranges (pseudocounts at 10.0, 0.01, 0.00125, 0.00125; thresholds at 1e-4; inbreeding at 0.0). Swapping any two would compile cleanly and pass engine-side range validators. The companion `ContaminationEstimationConfig` uses the struct-literal-plus-`validate()` pattern explicitly *because* "14+ numeric arguments would be hostile at the call site"; eight identically-typed `f64`s is hostile for the same reason.
- **Why it matters:** Silently-swappable. Use newtypes (`Pseudocount`, `Phred`, `FixationIndex`) for type-level discrimination OR adopt the struct-literal-plus-`validate()` shape so callers name fields. The sub-agent reviewing refactor_safety filed this at Major; demoted to Minor at synthesis because clap surfaces a single `Args` struct, so each call-site `args.x` binding is already named — the silent-swap risk is concentrated in `PosteriorEngineConfig::new`'s 8-argument signature, not at the clap boundary.
- **Suggested fix:** Move `PosteriorEngineConfig` to the validate-after-build pattern (matches contamination_estimation); see M4's fix sketch.

#### Mi3: [src/pop_var_caller/stage1_pipeline.rs:119-164](../../../src/pop_var_caller/stage1_pipeline.rs#L119-L164) — Per-branch local-variable snapshots in `with_stage1_pipeline` are not enforced by a shared struct shape
- **Categories:** refactor_safety
- **Confidence:** Medium
- **Problem:** Both branches end with the same triple of post-closure assignments (`stashed_upstream_error = …; baq_skip_counts = …;` etc.). A future "third counter" added in only one branch would compile silently. `Stage1RunSummary` exists exactly to be the typed snapshot; both branches should *build* one rather than assemble it from independently-set locals.
- **Why it matters:** As fields are added, the struct enforces presence; the current per-branch locals don't.
- **Suggested fix:** Replace the per-branch assignments with a per-branch construction of `(Stage1RunSummary, Option<CramInputError>)`. Then a new field on `Stage1RunSummary` fails to compile in both branches. Full sketch in [refactor_safety.md](../../../tmp/review_2026-05-19_cohort-cli/refactor_safety.md).

#### Mi4: [src/pop_var_caller/estimate_contamination.rs:582](../../../src/pop_var_caller/estimate_contamination.rs#L582) — Silent fallback `_ => 0` on a non-SidePass `ContaminationEstimateSource` reports `sites_processed = 0` to the user
- **Categories:** errors (Major in category), reliability (Minor in category), idiomatic (cross-cat)
- **Confidence:** High
- **Problem:** The match arm's own comment confirms "any other variant would mean the engine reshaped its API and the orchestrator wasn't updated" — an invariant-violation site. Today it silently shows `sites_processed=0` instead of surfacing the API drift. Per `feedback_no_logs_use_errors`, this should be a typed-error variant.
- **Why it matters:** The contamination subcommand is user-facing; a misreported site count masquerading as a successful run is exactly the failure mode the project's "no silent intermediates" + "no logs" memory targets.
- **Suggested fix:** Convert to a hard error:
  ```rust
  let sites_processed = match &estimates.source {
      ContaminationEstimateSource::SidePass { sites_processed, .. } => *sites_processed,
      other => return Err(EstimateContaminationCliError::UnexpectedSource {
          got: format!("{other:?}"),
      }),
  };
  ```
  Add the `UnexpectedSource { got: String }` variant. Demoted from Major (in errors category) to Minor at synthesis because the invariant is engine-internal and the wildcard arm is dead today; promotes to typed error mostly for future-proofing.

#### Mi5: [src/pop_var_caller/contamination_artifact.rs:86-141](../../../src/pop_var_caller/contamination_artifact.rs#L86-L141) — Public on-disk types lack `#[non_exhaustive]`; schema field additions become breaking
- **Categories:** extras
- **Confidence:** High
- **Problem:** `ContaminationArtefact`, `Provenance`, `ProvenanceInputs`, `BatchEntry`, `SampleEntry` are all `pub struct`s with `pub` fields and no `#[non_exhaustive]`. The error enum is correctly `#[non_exhaustive]` (line 351), but the structs that define the on-disk artefact format are not. Combined with the next finding (Mi6 — `provenance.version` is recorded but never validated), the contract surface is locked to the v1 shape.
- **Why it matters:** First breaking change the moment a new schema field lands (e.g. external-allele-frequency support listed in plan §"Open follow-ups").
- **Suggested fix:** Add `#[non_exhaustive]` to all five structs. Test fixtures that build the artefact via struct literal continue to compile because they are inside the crate.

#### Mi6: [src/pop_var_caller/contamination_artifact.rs:103-111](../../../src/pop_var_caller/contamination_artifact.rs#L103-L111), [:195-259](../../../src/pop_var_caller/contamination_artifact.rs#L195-L259) — `provenance.version` is recorded but never validated on read
- **Categories:** extras
- **Confidence:** High
- **Problem:** Schema records `provenance.version = "0.1.0"`, but `ContaminationArtefact::read` / `validate` never look at it. The field exists for forward-compat versioning but the gate isn't wired.
- **Why it matters:** When v2 lands, v1 consumers will either accept the file silently (if structurally compatible) or fail with a generic TOML parse error rather than a clear "this artefact was produced by version X; this binary supports up to Y".
- **Suggested fix:**
  ```rust
  const SUPPORTED_VERSIONS: &[&str] = &["0.1.0"];
  // in validate():
  if !SUPPORTED_VERSIONS.contains(&self.provenance.version.as_str()) {
      return Err(ContaminationArtefactError::UnsupportedVersion {
          got: self.provenance.version.clone(),
          supported: SUPPORTED_VERSIONS,
      });
  }
  ```

#### Mi7: [src/pop_var_caller/batch_assignment.rs:93](../../../src/pop_var_caller/batch_assignment.rs#L93), [:100-108](../../../src/pop_var_caller/batch_assignment.rs#L100-L108) — Header tolerates trailing `\r` via `trim_end()` but body rows do not; CRLF files store batch IDs with embedded `\r` ("Needs verification")
- **Categories:** extras, reliability
- **Confidence:** High (behaviour reasoned from `&str::lines()` semantics; not run against a CRLF fixture)
- **Problem:** Line 93 uses `header.trim_end() != "sample\tbatch"`, tolerating trailing `\r`. Body rows (line 100+) split on `\t` without trimming. A CRLF-terminated body line `NA12878\tlane_3\r` will store `batch = "lane_3\r"`. That label never matches a user-supplied batch reference, so the user gets a "dangling batch" error far from the real cause (line endings). The module's own docs (`:25-27`) say "Trailing whitespace is preserved verbatim — the consumer is expected to feed clean files," but the header's asymmetric `trim_end()` masks the cause.
- **Why it matters:** User-controlled input on a CLI surface should fail in a way that points at the real problem. The current behaviour silently propagates corrupt labels two layers downstream.
- **Suggested fix:** Tighten the header check to exact equality (`header != "sample\tbatch"`) so CRLF fails fast at line 1 with a clear error, OR call `line.trim_end_matches('\r')` on body lines for parity. Either way, add `body_rows_with_trailing_carriage_return_preserve_sample_name` test:
  ```rust
  let m = parse("sample\tbatch\nNA12878\tlane_3\r\n").unwrap();
  assert_eq!(m.batch_for("NA12878"), "lane_3");
  ```
  If the test exposes the corrupt-label behaviour, **promote this finding to Major**.

#### Mi8: helper duplication across 3–4 modules (`basename`, `format_md5_hex`, `current_command_line`, `rfc3339_now`, `civil_from_days`)
- **Categories:** smells, idiomatic, refactor_safety (cross-cat), extras (Nit), reliability (Nit)
- **Confidence:** High
- **Problem:** Pure, stateless byte-shuffling helpers duplicated:
  - `basename` — [cli.rs:494](../../../src/pop_var_caller/cli.rs#L494), [var_calling.rs:497](../../../src/pop_var_caller/var_calling.rs#L497), [estimate_contamination.rs:641](../../../src/pop_var_caller/estimate_contamination.rs#L641) (3×).
  - `format_md5_hex` — [cli.rs:502](../../../src/pop_var_caller/cli.rs#L502), [var_calling_from_bam.rs:684](../../../src/pop_var_caller/var_calling_from_bam.rs#L684) (2×).
  - `current_command_line` — [var_calling.rs:503](../../../src/pop_var_caller/var_calling.rs#L503), [var_calling_from_bam.rs:693](../../../src/pop_var_caller/var_calling_from_bam.rs#L693) (2×).
  - `rfc3339_now` + `civil_from_days` — [cli.rs:633/651](../../../src/pop_var_caller/cli.rs#L633), [estimate_contamination.rs:609/627](../../../src/pop_var_caller/estimate_contamination.rs#L609) (2× each, including the test bodies).
- **Why it matters:** Each duplicate carries a doc-comment justifying it as "self-contained subcommands". That rationale doesn't apply to these specific functions — they take primitives, have no per-subcommand variation, and the subcommands *already* share `cli/parsers`, `contamination_artifact`, `stage1_pipeline` modules. Future fix to any one (e.g. handling shell-quoting in `current_command_line`, or sub-second precision in `rfc3339_now`) has to be made 2–3× and is easy to miss. The `cargo test` evidence shows the tests themselves are textually duplicated (`civil_from_days_matches_known_dates` body identical between `cli.rs:739` and `estimate_contamination.rs:721`).
- **Suggested fix:** Add `src/pop_var_caller/common.rs` (or `cli/util.rs`) with `pub(crate)` versions of all five helpers and a single set of unit tests. Each subcommand `use super::common::{...};`. Visibility stays `pub(crate)`.

#### Mi9: [src/var_calling/posterior_engine.rs:129](../../../src/var_calling/posterior_engine.rs#L129) — `MAX_GQ_PHRED` is named like a hard cap but is the default
- **Categories:** defaults
- **Confidence:** High
- **Problem:** `pub const MAX_GQ_PHRED: f64 = 99.0;` is used as the `default_value_t` for `--max-gq-phred` in both `var_calling.rs:212` and `var_calling_from_bam.rs:338`. The actual validator upper bound is `GQ_PHRED_RANGE_MAX = 200.0` (line 154). A user reading `default_value_t = MAX_GQ_PHRED` believes the engine refuses anything above 99 — but `--max-gq-phred 150` is valid. Sibling constants follow `DEFAULT_*` / `*_RANGE_MAX` discipline; this one doesn't.
- **Why it matters:** "Visible defaults at the call site" — the default constant must announce that it's a default, not a cap.
- **Suggested fix:** Mechanical rename `MAX_GQ_PHRED` → `DEFAULT_MAX_GQ_PHRED`. Update three import sites + the const declaration. Update docstring to "Default GQ cap; the field itself is configurable up to `GQ_PHRED_RANGE_MAX`".

#### Mi10: [src/pop_var_caller/var_calling.rs:218-225](../../../src/pop_var_caller/var_calling.rs#L218-L225) and [src/pop_var_caller/var_calling_from_bam.rs:344-351](../../../src/pop_var_caller/var_calling_from_bam.rs#L344-L351) — `emit_gp` clap field has no `///` doc comment; `--help` shows no rationale
- **Categories:** defaults
- **Confidence:** High
- **Problem:** Every other Advanced flag has a `///` doc comment that clap renders as help text. `emit_gp` does not. `DEFAULT_EMIT_GP = false` has a useful three-line rationale in `vcf_writer/mod.rs:49`, but it's invisible in `--help`.
- **Why it matters:** The CLI is the primary surface. A defaulted flag without help text forces the user into source to discover what `--emit-gp` even does.
- **Suggested fix:** Add a doc comment mirroring the const's rationale:
  ```rust
  /// Emit `GP` (genotype posteriors) `FORMAT` per sample. Off by default
  /// — `GP` is `Number=G`, so the per-sample cell grows as
  /// `(ploidy + n_alleles − 1) choose ploidy` (21 floats at ploidy=2,
  /// n_alleles=6). Opt in when posteriors are wanted on disk.
  ```
  Apply to both files.

#### Mi11: [src/pop_var_caller/estimate_contamination.rs:144-150](../../../src/pop_var_caller/estimate_contamination.rs#L144-L150) — `--min-cohort-minor-count` has no `value_parser` and no engine-side check; `0` is accepted silently
- **Categories:** defaults
- **Confidence:** Medium
- **Assumption:** `min_cohort_minor_count = 0` is an absurd choice rather than a deliberate "no floor" option (neighbouring `--min-depth` and `--min-batch-size` reject it).
- **Problem:** The flag has `default_value_t = DEFAULT_MIN_COHORT_MINOR_COUNT` but no `value_parser`. The engine's `ContaminationEstimationConfig::validate` does not check the field either. A user passing `--min-cohort-minor-count 0` produces a silent semantic shift (every site qualifies on the count axis; the fraction axis becomes the sole filter).
- **Why it matters:** The slice's defence-in-depth intent. This field is the only count knob in the side-pass that escapes both layers.
- **Suggested fix:** Add `parse_min_cohort_minor_count` rejecting `0` to `cli/parsers.rs`, wire it, and add an engine-side check in `validate`.

#### Mi12: [src/pop_var_caller/contamination_artifact.rs:280-342](../../../src/pop_var_caller/contamination_artifact.rs#L280-L342) — `to_estimates_for_samples` does not run `validate()` and partially-defends against invariants `validate()` would catch
- **Categories:** errors, smells (cross-cat — two-pass walk)
- **Confidence:** Medium
- **Problem:** The method assumes `validate()` was called (per the inline comment), but the public-API surface allows construction via all-pub fields without `validate()`. The defensive arm at L306-314 guards only the dangling-sample-batch case; other invariants (simplex tolerance, duplicate ids, non-finite probabilities) are trusted blindly. Inconsistent.
- **Why it matters:** Either trust the invariant unconditionally (`unreachable!` with `// UNREACHABLE:` comment) or run `self.validate()` at the top. The current half-and-half lets a bad row through quietly while loudly catching one specific case.
- **Suggested fix:** Call `self.validate()?` at the top of `to_estimates_for_samples`. Cost: one BTreeMap walk per artefact load (already paid in `read()`); benefit: the defensive arm becomes truly unreachable and can be deleted; the function becomes safe to call from hand-built artefacts. Pairs with simplifying the two-pass walk into a single pass (sketch in [smells.md](../../../tmp/review_2026-05-19_cohort-cli/smells.md)).

#### Mi13: [src/pop_var_caller/contamination_artifact.rs:87](../../../src/pop_var_caller/contamination_artifact.rs#L87) — File path `contamination_artifact.rs` (American) holds type `ContaminationArtefact` (British)
- **Categories:** naming
- **Confidence:** High
- **Problem:** Filename uses "artifact"; type and every consumer use "artefact". Mixed in the same module path (`pop_var_caller::contamination_artifact::ContaminationArtefact`).
- **Why it matters:** Search-and-replace across the file name and the type yields disjoint hit sets. Cosmetic, but the path reads jarringly.
- **Suggested fix:** Rename the file `contamination_artifact.rs` → `contamination_artefact.rs` (and `mod.rs:8` entry) to match the type and every consumer. The type and every reference already use "artefact" so this is the lower-blast-radius direction.

#### Mi14: [src/pop_var_caller/estimate_contamination.rs:172](../../../src/pop_var_caller/estimate_contamination.rs#L172) — CLI knob name `min_batch_size` diverges from the engine field `min_batch_size_for_contamination`
- **Categories:** naming
- **Confidence:** High
- **Problem:** Three names for the same concept: engine field `min_batch_size_for_contamination`, CLI field/flag `min_batch_size` / `--min-batch-size`, artefact parameter-map key `"min_batch_size"`. A future consumer reading the artefact must know that `min_batch_size` *means* `min_batch_size_for_contamination` — undocumented.
- **Why it matters:** Reproducibility lookups: the artefact's parameter map is the contract between subcommands. Same concept, three names, breaks "same name everywhere".
- **Suggested fix:** Pick one. Cheapest: rename the CLI field + flag + parameter-map key to match the engine (`min_batch_size_for_contamination` / `--min-batch-size-for-contamination`). Flag becomes long but lives under `Advanced — Informative-site cuts`, so it's hidden from `-h`. Verbose `--help` text for rarely-touched knobs is acceptable. Alternative (out of slice scope): rename the engine field.

#### Mi15: [src/var_calling/variant_grouping.rs:195-211](../../../src/var_calling/variant_grouping.rs#L195-L211) — `Debug` impl on `VariantGrouper` over-constrains its generics
- **Categories:** idiomatic
- **Confidence:** Medium
- **Problem:** The `Debug` impl carries `where I: Iterator<Item = Result<…, E>>, GrouperError: From<E>` even though the body never touches `E`. The exhaustive destructure prints `upstream`, `config`, `lookahead`, `done` — none typed with `E`. Adding the bound forces every Debug call site to satisfy `From<E>` that has nothing to do with formatting.
- **Why it matters:** Minimum-bound rule. A trait impl's `where` should list only what the impl body needs.
- **Suggested fix:** Drop the bounds:
  ```rust
  impl<I> std::fmt::Debug for VariantGrouper<I> {
      fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
          let Self { upstream: _, config, lookahead, done } = self;
          f.debug_struct("VariantGrouper")
              .field("config", config)
              .field("lookahead_is_some", &lookahead.is_some())
              .field("done", done)
              .finish_non_exhaustive()
      }
  }
  ```

#### Mi16: [src/pop_var_caller/cli/parsers.rs:115](../../../src/pop_var_caller/cli/parsers.rs#L115), :129, :249, :254, :260, :265 — Six range checks against `u32::MAX` are no-ops
- **Categories:** idiomatic
- **Confidence:** High
- **Problem:** `parse_var_group_max_span`, `parse_dust_threshold`, `parse_block_size`, `parse_min_depth`, `parse_min_batch_size`, `parse_stability_blocks` all call `parse_u32_in(..., ..., u32::MAX)`. Every parsed `u32` trivially satisfies `v <= u32::MAX`. Reads as "in `1..=u32::MAX`" — suggesting the upper bound was considered policy when it isn't.
- **Why it matters:** A future maintainer might raise the upper bound under the impression it's a real cap. Error message renders as e.g. "must be in 1..=4294967295" which is silly.
- **Suggested fix:** Introduce two helpers:
  ```rust
  fn parse_u32_min(s: &str, name: &str, lo: u32) -> Result<u32, String> { /* "must be >= lo" */ }
  fn parse_u32_any(s: &str, name: &str) -> Result<u32, String> { /* parse, return */ }
  ```
  Wire the affected sites accordingly.

#### Mi17: [src/pop_var_caller/contamination_artifact.rs:287-301](../../../src/pop_var_caller/contamination_artifact.rs#L287-L301) — Idioms: `BTreeMap` used for lookup-only maps; double-deref `(*cohort_name).to_string()` over `&&str`
- **Categories:** idiomatic
- **Confidence:** Medium
- **Problem:** Two issues at the same site:
  1. `by_name: BTreeMap<&str, &SampleEntry>` and `by_batch: BTreeMap<&str, &BatchEntry>` are built only to support `.get()` / `.contains_key()`. Neither is iterated. `HashMap` expresses lookup-only intent.
  2. `for cohort_name in sample_names` (`sample_names: &[&str]`) binds `cohort_name: &&str`. `(*cohort_name).to_string()` derefs twice; explicit `*` plus method that *also* implicit-derefs reads as noise. Rust idiom: `for &cohort_name in sample_names { ... cohort_name.to_string() ... }`.
- **Why it matters:** Reader signal — `BTreeMap` says "I care about iteration order"; `HashMap` says "lookup table". *(Rust-specific note for the author who's new to Rust: `for &x in slice_of_refs` is pattern-match destructure — strips the outer `&`. Equivalent to `for x in slice_of_refs.iter().copied()` for `Copy` inner types like `&str`. Preferred over explicit `*` derefs when the inner type implements `Copy`.)*
- **Suggested fix:** Change the three `BTreeMap<&str, _>` declarations to `HashMap<&str, _>`; the `Entry::Vacant` branch has an `std::collections::hash_map::Entry::Vacant` equivalent. Change the two `for cohort_name in sample_names` loops to `for &cohort_name in sample_names`.

#### Mi18: [src/pop_var_caller/estimate_contamination.rs:433-498](../../../src/pop_var_caller/estimate_contamination.rs#L433-L498) — `build_artefact_from_estimates` picks "representative sample" per batch via linear `position()`; engine accessor shape doesn't fit the batch-keyed call site
- **Categories:** smells
- **Confidence:** High
- **Problem:** Lines 456-480 build each `BatchEntry` by linearly scanning `estimates.sample_to_batch` for the first sample with that batch index, then calls the sample-keyed `q_b_for_sample(s)`. The code comment admits the awkwardness ("would mean exposing more"). The `match representative_sample { Some(s) => *q_b_for_sample(s), None => [0.0; 3] }` fall-through at L469-472 silently maps "unused batch" to "all-zero floored q_b". The validator accepts all-zero rows, so the downstream `var-calling` reload treats it as floored — that *is* the intended behaviour, but it deserves a doc comment more explicit than the current "Falls back to all-zero (floored signal) if the engine returned an unused batch."
- **Why it matters:** Two distinct smells: (1) the engine's accessor API doesn't fit batch-keyed export; (2) the fall-through buries an invariant ("unused batch should never happen given how `build_dense_batches` constructs the table").
- **Suggested fix:** Add `pub(crate) fn q_b_per_batch(&self) -> &[[f64; N_ALLELE_CLASSES]]` to `ContaminationEstimates` and index directly. Drops the linear scan and the fall-through. If exposing the accessor is undesirable, convert the fall-through into `expect("dense batch must have ≥1 sample by construction")`.

#### Mi19: [src/pop_var_caller/var_calling.rs:329](../../../src/pop_var_caller/var_calling.rs#L329), [src/pop_var_caller/estimate_contamination.rs:315](../../../src/pop_var_caller/estimate_contamination.rs#L315), [src/pop_var_caller/cli.rs:336](../../../src/pop_var_caller/cli.rs#L336) — Magic `64 * 1024` BufReader/BufWriter capacity in 3+ places
- **Categories:** smells
- **Confidence:** High
- **Problem:** The slice uses `64 * 1024` for buffered IO at three new sites; the broader `pop_var_caller` module uses it at five more. `BLOCK_HEADER_READ_CAP: usize = 64 * 1024` already exists in [psp/reader.rs:54](../../../src/per_sample_pileup/psp/reader.rs#L54) — conceptually the same constant.
- **Why it matters:** A future tuning pass that identifies a better buffer size has to be fanned out across 8+ sites.
- **Suggested fix:** Promote `pub(crate) const DEFAULT_BUFFERED_IO_CAPACITY: usize = 64 * 1024;` in `pop_var_caller::common` (or `cli/util.rs`). Replace the three in-scope sites. Leave the orthogonal `BLOCK_HEADER_READ_CAP` alone — different intent.

#### Mi20: [tests/cohort_cli_integration.rs:73-216](../../../tests/cohort_cli_integration.rs#L73-L216) and [tests/pileup_cli_integration.rs:56-149](../../../tests/pileup_cli_integration.rs#L56-L149) — Fixture-helper duplication across integration test files
- **Categories:** smells, reliability (cross-cat)
- **Confidence:** High
- **Problem:** `build_fasta`, `build_sam_header`, `build_cram`, `read_record`, and `pileup_args` are duplicated near-verbatim. CRAM-emitter setup is ~40 lines of noodles plumbing.
- **Why it matters:** Maintenance fan-out concentrated in test-only code; lower stakes than production but the helpers embed format invariants (which `@SQ M5`, which read flags) that need to move in lockstep on a format touch.
- **Suggested fix:** Create `tests/common/mod.rs` with the shared helpers and `mod common;` from both files. Cargo's test discovery treats `tests/common/mod.rs` as a non-test module.

#### Mi21: [src/pop_var_caller/estimate_contamination.rs:348-365](../../../src/pop_var_caller/estimate_contamination.rs#L348-L365) — `ContaminationEstimationConfig::validate(&self)` is the only Stage-config in the slice without `new(...) -> Result`; library callers can hold an un-validated typed value
- **Categories:** extras, refactor_safety (cross-cat)
- **Confidence:** High
- **Problem:** Per the plan §"Engine-side validation", every Stage 3–6 config was supposed to get `Config::new -> Result`. The impl chose `validate(&self)` for `ContaminationEstimationConfig` (14+ fields). The trade-off is documented in the impl report's §Assumptions, but it weakens "library users see the same guard CLI users see": a library caller can construct the struct via literal and pass it to `estimate_contamination()` without ever calling `validate()`.
- **Why it matters:** Asymmetric guarantee across the five Stage configs in the slice. Holding a typed `ContaminationEstimationConfig` value tells you nothing about whether it's valid.
- **Suggested fix:** Either (a) add a `new()` taking a builder-shaped grouping of related sub-configs (StoppingMode, InformativeSiteCuts, Priors), or (b) make the public fields `pub(crate)` so the only public construction path is `with_project_defaults().with_x(...).build()?` where `build()` runs `validate()`. (b) is the smaller surface change.

#### Mi22: [src/pop_var_caller/var_calling_from_bam.rs:425-431](../../../src/pop_var_caller/var_calling_from_bam.rs#L425-L431) — Rayon-already-configured error surfaces as `Stage1(...)` despite the call being subcommand-side
- **Categories:** extras, errors (Nit)
- **Confidence:** Medium
- **Problem:** `run_var_calling_from_bam` initialises the rayon pool itself at L426-430, then `.map_err(|_| PileupCliError::RayonAlreadyConfigured)?`. The `?` lifts via `VarCallingFromBamCliError::Stage1(#[from] PileupCliError)`, so the chain reads `"stage 1: rayon thread pool already initialised"`. But rayon init happened in the from-bam orchestrator *before* any Stage 1 work — it's not a Stage 1 error semantically. The other two subcommands have their own `RayonAlreadyConfigured` variant.
- **Why it matters:** Misleading error attribution on a CLI surface.
- **Suggested fix:** Add a dedicated `RayonAlreadyConfigured` variant to `VarCallingFromBamCliError` and `.map_err(|_| VarCallingFromBamCliError::RayonAlreadyConfigured)?`. Matches the other two subcommands.

#### Mi23: [tests/cohort_cli_integration.rs](../../../tests/cohort_cli_integration.rs) — 3 of 10 plan-listed integration tests landed; the load-bearing chain test (`estimate-contamination → var-calling`) is deferred
- **Categories:** extras, reliability
- **Confidence:** High
- **Problem:** The plan §"Integration tests" listed 10 cases. The impl report acknowledges 3 landed and 7 are "covered in spirit by per-module unit tests". The chained test (`estimate-contamination` produces artefact → `var-calling --contamination-estimates` consumes it) is *not* covered in spirit by per-module tests — only an integration target exercises both subcommands in one process, and that is exactly where bugs at the artefact-handoff boundary would surface.
- **Why it matters:** The artefact handoff is the *whole point* of this slice. A sample-order permutation (very plausible regression target: sort the readers by name) would silently corrupt every downstream contamination correction with no test catching it.
- **Suggested fix:** Land the chained test before the slice is considered done:
  ```rust
  #[test]
  fn estimate_contamination_then_var_calling_chain() {
      let dir = TempDir::new().unwrap();
      let fasta = build_fasta(dir.path());
      let psp_a = make_psp_for_sample(dir.path(), &fasta, "NA00001", &[read_record("r1", 10, b"ACAAA")]);
      let psp_b = make_psp_for_sample(dir.path(), &fasta, "NA00002", &[read_record("r1", 10, b"AAAAA")]);
      let toml = dir.path().join("contam.toml");
      run_estimate_contamination(&estimate_args(fasta.clone(), toml.clone(),
                                                  vec![psp_a.clone(), psp_b.clone()]))
          .expect("estimate-contamination");
      let vcf = dir.path().join("cohort.vcf");
      run_var_calling(&var_calling_args(fasta, vcf.clone(),
                                          vec![psp_a, psp_b], Some(toml)))
          .expect("var-calling");
      assert!(vcf.exists());
  }
  ```
  Document the other 6 deferred cases as `#[ignore]`-marked stubs so they don't get lost.

### Nits

Grouped under one heading per skill convention (low-volume, all clustering around documentation/comment polish or one-line clippy-style nits):

- **Helper-comment style.** [stage1_pipeline.rs:127](../../../src/pop_var_caller/stage1_pipeline.rs#L127) — `expect("ref_id fits u32")` is a genuine panic site; add `// PANIC-FREE: …` invariant comment per the project pattern. Same applies to [var_calling_from_bam.rs:688](../../../src/pop_var_caller/var_calling_from_bam.rs#L688) and [cli.rs:506](../../../src/pop_var_caller/cli.rs#L506) (`expect("writing to a String never fails")`).
- **Best-effort cleanup comments.** [var_calling_from_bam.rs:558](../../../src/pop_var_caller/var_calling_from_bam.rs#L558), [:596](../../../src/pop_var_caller/var_calling_from_bam.rs#L596) — `let _ = std::fs::remove_file(...)` cleanups should carry the same `// best-effort cleanup; rename may not have happened` comment that [cli.rs:354](../../../src/pop_var_caller/cli.rs#L354) does.
- **`#[allow(clippy::too_many_arguments)]` lacks justification.** [stage1_pipeline.rs:87](../../../src/pop_var_caller/stage1_pipeline.rs#L87) — add one-line comment explaining why the helper's eight Stage 1 config args can't reasonably collapse further.
- **Vestigial `let _ = cfg;` placeholder.** [estimate_contamination.rs:566](../../../src/pop_var_caller/estimate_contamination.rs#L566) — `cfg` has already been used to drive the side-pass at L373; the placeholder reads as dead code. Drop it or use it. Pairs with the M11 driver extraction.
- **`DEFAULT_BATCH_ID = "all_samples"`** — the constant name (`_BATCH_ID`, singular) and string value (`all_samples`, plural) pull in different directions; either rename or relabel. Low-priority.
- **`Stage1RunSummary` oversells** — carries `FilterCounts` + `Option<BaqSkipCounts>` only. Consider `stage1_counters`.
- **`Stage1PipelineContext` / `Stage1RunSummary`** repeat the module name. Inside `stage1_pipeline.rs`, `Stage1Context` / `Stage1Counters` would carry the same meaning; the `Stage1` prefix is the load-bearing disambiguator at external use sites, not `Pipeline`.
- **`"all_samples"` literal duplicated** in [estimate_contamination.rs:76](../../../src/pop_var_caller/estimate_contamination.rs#L76) and tests at `:663`, `:710` instead of `DEFAULT_BATCH_ID`. The runtime code uses the constant correctly; only doc + tests drift.
- **Inline absolute-path `#[from]`** at [var_calling.rs:267](../../../src/pop_var_caller/var_calling.rs#L267) and [var_calling_from_bam.rs:397](../../../src/pop_var_caller/var_calling_from_bam.rs#L397) — use `use` import for consistency with the rest of the variants.
- **Redundant turbofish.** [var_calling_from_bam.rs:466](../../../src/pop_var_caller/var_calling_from_bam.rs#L466) — `with_stage1_pipeline::<(), VarCallingFromBamCliError, _>` can rely on inference now that the closure body explicitly annotates `Result<(), VarCallingFromBamCliError>`.
- **Unused `#[from]` variants** — several CLI error variants are unreachable via `?` (`VarCallingCliError::Io`, `::Grouper`, `::PerGroup`; `VarCallingFromBamCliError::Io`, `::Grouper`, `::PerGroup`, `::Psp`; `EstimateContaminationCliError::Io`). Each invites a future `?` site to silently pick the wrong variant. Drop `#[from]` (and the variants) or keep without `#[from]` so future call sites have to choose explicitly.
- **Helper `_` => 0 fallback comment in `print_run_summary`** ([estimate_contamination.rs:580-583](../../../src/pop_var_caller/estimate_contamination.rs#L580-L583)) — exhaustive-match rule. Add `// REVIEW ON UPGRADE:` marker; pairs with Mi4 if that finding is applied (the wildcard arm goes away entirely).

## 7. Out of scope observations

These are pre-existing or out-of-scope per the task brief but worth noting for follow-ups:

- **Pre-existing `cargo doc` errors** outside the slice: [pileup_to_psp.rs:5-6](../../../src/per_sample_pileup/pileup_to_psp.rs#L5-L6) — `pileup::PileupWalker`, `psp::writer::PspWriter` (drift in legacy paths); [posterior_engine.rs:725/732/753](../../../src/var_calling/posterior_engine.rs#L725) — 3× `ExactMath` (the engine refactor PROJECT_STATUS notes left these open under Mi21). Should land in a follow-up doc-cleanup pass alongside M14.
- **Pre-existing failing `src/gvcf_parser.rs` doctest** — legacy code marked superseded; tracked separately.
- **`tool_string` / `WriterConfig` rename follow-ups** from the vcf_writer review (Mi11/Mi12) — touch this slice's two call sites once the writer-side rename lands.
- **`GrouperError::VariantGroupTooWide` message points at the Rust path `GrouperConfig::max_variant_group_span`, not the CLI flag** the user would actually edit. Error-message review item, not in scope here.
- **`command_line` field built from `std::env::args_os()` joined with spaces** ([var_calling.rs:503](../../../src/pop_var_caller/var_calling.rs#L503), [var_calling_from_bam.rs:693](../../../src/pop_var_caller/var_calling_from_bam.rs#L693)) — loses shell-quoting; argv containing a space round-trips as two tokens. Cosmetic; provenance/observability concern.
- **`run_var_calling` is ~115 lines** with ~10 numbered steps but each step is 1–3 lines — flagged in the smells review as a long-function borderline. Not filed because the steps are linear and the function reads as a script.

## 8. Missing tests to add now

Grouped by function under test. Naming follows `function_returns_expected_on_condition`.

**`with_stage1_pipeline` ([stage1_pipeline.rs:88](../../../src/pop_var_caller/stage1_pipeline.rs#L88))** — currently has zero direct tests despite owning the self-referential borrow chain:
- `with_stage1_pipeline_stashes_cram_input_error_when_walker_exhausts_cleanly` — CRAM with corrupted block past first valid block; assert `outs.stashed_upstream_error.is_some()` and `outs.result.is_ok()`. **Pairs with M1 fix.**
- `with_stage1_pipeline_prefers_upstream_error_over_closure_error` — closure returns Err, stash is Some; assert the upstream error is surfaced. **Pairs with M1 fix.**
- `with_stage1_pipeline_baq_on_branch_returns_some_baq_skip_counts` — CRAM with one BAQ-processed read; assert `outs.run_summary.baq_skip_counts.is_some()`.
- `with_stage1_pipeline_no_baq_branch_returns_none_baq_skip_counts` — same CRAM with `no_baq = true`; assert `is_none()`.

**`run_var_calling_from_bam` ([var_calling_from_bam.rs:421](../../../src/pop_var_caller/var_calling_from_bam.rs#L421))** — walker-error path untested:
- `var_calling_from_bam_surfaces_walker_error_and_removes_output` — reads dense enough to trip `max_active_reads = 1`; assert `Err(Walker(_))` and that `output.exists()` is false. **Pairs with M2 fix.**
- `walker_error_shim_stashes_and_ends_stream` — unit test on the scan-adapter pattern in isolation. **Pairs with M9 fix.**

**`run_var_calling` ([var_calling.rs:316](../../../src/pop_var_caller/var_calling.rs#L316))**:
- `var_calling_emits_at_least_one_variant_record` — strengthen `var_calling_happy_path_three_samples` to assert at least one non-comment data line, ideally locating POS=11 with `ALT=C`. (Header-only success currently passes the test.)
- `var_calling_reports_reference_mismatch_for_psp_with_different_header_basename` — build two FASTAs with distinct basenames; run pileup against one, then var-calling pointing at the other. Assert `ReferenceMismatch`.
- `var_calling_removes_tmp_on_iterator_error` — drive `run_var_calling` against a `.psp` that errors mid-stream (e.g. via a `PspReadError::Format` injection). Assert `<output>.tmp` does not exist after. **Pairs with M6 fix.**

**`run_estimate_contamination` ([estimate_contamination.rs:300](../../../src/pop_var_caller/estimate_contamination.rs#L300))** — currently zero integration coverage:
- `estimate_contamination_happy_path_writes_loadable_artefact` — 3 `.psp`s, 2-batch `--batch-assignment` TSV; assert the artefact loads + every cohort sample appears in `[[samples]]` + every batch in `[[batches]]` + `to_estimates_for_samples(cohort_names)` succeeds.
- `estimate_contamination_rejects_reference_mismatch` — .psp with wrong `header.reference`; assert `ReferenceMismatch { .. }`.
- `estimate_contamination_round_trips_into_var_calling` — chain the two subcommands; assert no `SampleMissingFromArtefact`. **Pairs with Mi23 fix.**

**`load_contamination` ([var_calling.rs:439](../../../src/pop_var_caller/var_calling.rs#L439))** — only the no-path branch is covered:
- `load_contamination_returns_estimates_when_path_present_and_cohort_matches` — fixture artefact + cohort match; assert `effective_c_s(0) == 0.0123`.
- `load_contamination_propagates_missing_sample_error` — cohort sample not in artefact; assert `Err(ContamArtefact(SampleMissingFromArtefact { sample: "NA00001" }))`.
- `load_contamination_silently_drops_extras` — artefact superset; runs cleanly.

**`ContaminationEstimationConfig::validate` ([contamination_estimation.rs:282](../../../src/var_calling/contamination_estimation.rs#L282))**:
- `config_validate_rejects_zero_snp_alt_pseudocount` — the loop today only tests `ref_pseudocount = 0.0`; this covers a future field-list drop.
- `config_validate_rejects_zero_indel_alt_pseudocount` — same shape.
- `config_validate_rejects_negative_q_b_init` — only the `> 1.0` case is tested.
- `config_validate_rejects_min_major_above_one` — only the lower-boundary case is tested.
- `config_validate_accepts_fixed_sites_with_positive_n` — ensures the `FixedSites` arm is exercised on the happy path, not just the zero-rejection arm.

**`ContaminationArtefact::to_estimates_for_samples` ([contamination_artifact.rs:280](../../../src/pop_var_caller/contamination_artifact.rs#L280))**:
- `to_estimates_round_trips_floored_batch_with_zero_qb` — all-zero `q_b` row + sample mapped to it; assert `q_b_for_sample(0) == [0.0; 3]`, not `QbNotSimplex` engine error. **Tests Assumption 6's two-sides-agree contract.**
- `to_estimates_orders_dense_batches_by_cohort_first_seen` — cohort spans 2+ batches, ordering distinct from artefact declaration; assert the dense table preserves cohort-first-seen order. **Tests the silent design choice in the dense-batch ordering invariant.**
- `to_estimates_for_samples_returns_engine_estimates_on_empty_cohort` — empty `sample_names`; assert empty result without panic.

**`contigs_to_parsed` / `format_md5_hex` ([var_calling_from_bam.rs:653](../../../src/pop_var_caller/var_calling_from_bam.rs#L653), [:684](../../../src/pop_var_caller/var_calling_from_bam.rs#L684))**:
- `contigs_to_parsed_md5_hex_matches_known_bytes` — md5 = `[0x00, 0x01, …, 0x0f]`; assert `parsed[0].md5 == "000102030405060708090a0b0c0d0e0f"`. (Today the test only checks `len() == 32`.)
- `format_md5_hex_matches_cli_helper_for_same_input` — call both duplicates on the same input; assert they agree. **Cheap drift-detector for Mi8.**

**`BatchAssignment::from_str_with_path` ([batch_assignment.rs:83](../../../src/pop_var_caller/batch_assignment.rs#L83))**:
- `body_rows_with_trailing_carriage_return_preserve_sample_name` — input `"sample\tbatch\nNA12878\tlane_3\r\n"`; assert `batch_for("NA12878") == "lane_3"`. **If this test fails, promote Mi7 to Major.**

**`ContaminationArtefact::write` ([contamination_artifact.rs:169](../../../src/pop_var_caller/contamination_artifact.rs#L169))**:
- `write_to_unwritable_directory_does_not_leave_tmp` — write to `/nonexistent_dir/x.toml`; assert `Err(Io)` and no `.tmp` left in cwd.

**`WriterConfig` shape** (after M8 fix):
- `writer_config_new_pins_default_emit_gp` — assert `WriterConfig::new(path).emit_gp == DEFAULT_EMIT_GP`. Locks the constructor's responsibility.

## 9. What's good

Five concrete patterns worth keeping:

1. **`with_stage1_pipeline`'s callback shape** ([stage1_pipeline.rs:88](../../../src/pop_var_caller/stage1_pipeline.rs#L88)) cleanly resolves the self-referential borrow chain the plan flagged as the highest-risk piece. The `E: From<PileupCliError>` lift on the closure's error type is a smart trick: the helper's internal setup errors lift into the closure's error type without changing the caller's surface.
2. **Engine-side `Config::new -> Result` validators** with range constants pulled out as `pub const *_RANGE_MAX` ([posterior_engine.rs:131-152](../../../src/var_calling/posterior_engine.rs#L131-L152), [per_group_merger.rs:99-105](../../../src/var_calling/per_group_merger.rs#L99-L105)) — both the CLI parser and the constructor reference the same constant, so the two-layer defence stays in sync mechanically.
3. **`VariantGrouper` generalisation over upstream error type** ([variant_grouping.rs:142-203](../../../src/var_calling/variant_grouping.rs#L142-L203)) is a clean rework: the struct itself stays bound-free, the bounds live on the impls, and the chain accepts both `PerPositionMerger` and `DustFilter` without an adapter. (Modulo Mi15 on the `Debug` impl's bound width.)
4. **Atomic tmp-then-rename idiom carried through `ContaminationArtefact::write`** ([contamination_artifact.rs:169-189](../../../src/pop_var_caller/contamination_artifact.rs#L169-L189)) — applied to the new TOML output even though most bioinformatics CLIs don't bother. Consistent with `run_pileup`.
5. **`#[non_exhaustive]` discipline** on every new error enum (`GrouperConfigError`, `PerGroupMergerConfigError`, `PosteriorEngineConfigError`, `ContaminationArtefactError`, `BatchAssignmentError`, all three CLI error enums) plus exhaustive struct destructures in the `Debug` impls — refactor-safe by design.

## 10. Commands to re-verify

Run after applying fixes:

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings
./scripts/dev.sh cargo test --lib                 # was 880 passed
./scripts/dev.sh cargo test --tests
./scripts/dev.sh cargo doc --no-deps              # must drop to 5 errors (pre-existing) once M14 lands
```

New invocations introduced by this review (when fixes land):
- For M1 / M2 / M9 — add `tests/cohort_cli_integration.rs::walker_error_*` cases or `stage1_pipeline.rs::#[cfg(test)] mod tests`.
- For Mi7 — add the CRLF body-row test; if it exposes a corrupt-label bug, promote to Major.
- For Mi23 — land the `estimate_contamination_then_var_calling_chain` end-to-end test.

### Author response convention

Address each finding by its identifier (`B`, `M`, `Mi`, or `Nit` group) with one of:
- `fixed in <commit>`
- `disputed because …`
- `deferred to <issue>`
- `won't fix because …`

Answer open questions in §4 first; they gate several findings.
