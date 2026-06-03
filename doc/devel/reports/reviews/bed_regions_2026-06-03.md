# Code Review: bed-regions
**Date:** 2026-06-03
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** the `bed-regions` feature diff (`--regions <bed>` for `pileup` + `var-calling`), merge-base `4efc627` → `HEAD` on `main`
**Status:** Approve-with-changes

---

## 1. Scope

- **Reviewed:** the merged `bed-regions` feature diff — `git diff 4efc627 HEAD`, source only (1722 insertions / 129 deletions across 17 `.rs` files).
- **Reviewed against:** `main` @ `89b5646` (feature merge `6ecd22c` + reports + tests).
- **In-scope files:** `src/regions.rs` (new), `src/bam/{alignment_input,cram_input,bam_input,errors}.rs`, `src/pop_var_caller/{cli,stage1_pipeline,var_calling,contamination_chunked_stream}.rs`, `src/pileup/per_sample/{pileup_to_psp,baq_stream}.rs`, `src/pileup/walker/driver.rs`, `src/var_calling/driver.rs`, `src/psp/{header,reader,test_fixtures}.rs`, `src/lib.rs`.
- **Out of scope:** the test/example/doc changes (the multi-contig regression tests were added separately); the pre-existing walker chrom-transition tail-drop bug (documented, dormant on this branch — see [bed_regions_perf_and_byte_identity_2026-06-03.md](../bed_regions_perf_and_byte_identity_2026-06-03.md)).
- **Categories dispatched:** reliability, errors, naming, defaults, idiomatic, refactor_safety, module_structure, unsafe_concurrency, smells, extras (10). Skipped: tooling (no `Cargo.toml` change).

## 2. Verdict

**Approve-with-changes.** No Blockers. The feature is cleanly built — correct module placement, lockstep concurrency preserved, idiomatic throughout, all gates green. The Major findings are defensive/maintainability fixes (future-field-safety, a constructible panic, error context, duplication) plus one process item (an unrelated flag folded in). None is an active correctness bug on the supported input path.

## 3. Execution status

All run in the dev container (per `CLAUDE.md`); output verbatim:

- `cargo fmt --check` → exit 0, no diff.
- `cargo clippy --all-targets --all-features -- -D warnings` → exit 0, `Finished` (no warnings).
- `cargo doc --no-deps` → exit 0, `Generated`.
- `cargo test --all-features` → `test result: ok. 1101 passed; 0 failed; 1 ignored`; all integration suites pass (cohort 12, pileup 9, others unchanged).
- `cargo audit` → not run (not required; no dependency changes in the diff).

Findings labeled "Needs verification": 0 (the one most-routed correctness concern, `restrict_intervals_to_regions`, was verified by the orchestrator — see §6 note).

## 4. Open questions and assumptions

1. **Is `--build-map-file-index` intended as part of this feature?** It's index-management, not BED parsing, but region seeking needs an index. Affects **M5**. If yes, no code change — the PR/commit description should say so.
2. **Does any caller need to distinguish pre-flight vs per-input index-load failures?** Affects **Mi3** (the `#[from] AlignmentIndexError` bridge collapses both). If no, the current shape is fine.
3. **Should an empty resolved `RegionSet` (e.g. a comments-only BED) warn instead of silently producing empty output?** Affects **Mi7** — a behaviour decision.
4. **Should the `var-calling` VCF record region provenance** the way the pileup `.psp` header now does? Affects **Mi6**.

## 5. Top 3 priorities

1. **M1** — the three new `merge()` counter-folds sum fields by hand; a future field is silently dropped from the now-load-bearing per-region totals. Exhaustive-destructure + add tests. (Convergent: refactor_safety, idiomatic, reliability.)
2. **M2** — `ContigInterval` has `pub` fields and no validating constructor, so `start == 0` is constructible and aborts the pileup at `query_interval`'s `Position::new(..).expect()`. Add a checked constructor / drop `pub`.
3. **M3** — `BedError::Io` drops the path + line and flattens its `io::Error` into `Display` (no `#[source]`): a mid-stream BED read failure is undebuggable.

## 6. Findings

### Blocker
None.

### Major

**M1: `merge()` counter-folds are not field-addition-safe — [src/pileup/walker/driver.rs:260](../../../../src/pileup/walker/driver.rs#L260), [src/bam/alignment_input.rs:148](../../../../src/bam/alignment_input.rs#L148), [src/pileup/per_sample/baq_stream.rs:71](../../../../src/pileup/per_sample/baq_stream.rs#L71)**
**Categories:** refactor_safety, idiomatic, reliability · **Confidence:** High
`RunSummary::merge` (8 fields, `+` except `active_reads_high_water` which is `max`), `FilterCounts::merge` (10, all `+`), and `BaqSkipCounts::merge` (11, all `+`) sum fields by name. All three currently cover every field, but nothing forces a *new* field to be handled — adding a counter compiles cleanly while `merge` silently omits it from every multi-region pileup's totals (and every whole-genome run now goes region-by-region). `RunSummary` is even `#[non_exhaustive]`. Silent wrong-total, no compile error, no panic, and there is **no unit test** on any of the three.
**Fix:** destructure `other` so a new field is a compile error, keeping the `max` special-case:
```rust
pub fn merge(&mut self, other: &RunSummary) {
    let RunSummary { reads_admitted, records_emitted, record_widen_events,
        mate_overlap_positions, chain_allocations, active_reads_high_water,
        mate_lookup_evictions, column_depth_truncations } = other;
    self.reads_admitted += reads_admitted; /* … each by name … */
    self.active_reads_high_water = self.active_reads_high_water.max(*active_reads_high_water);
}
```
Plus the three tests in §8.

**M2: `ContigInterval` is constructible with `start == 0`, which aborts the pileup — [src/bam/bam_input.rs:288](../../../../src/bam/bam_input.rs#L288), def at [src/bam/alignment_input.rs:540](../../../../src/bam/alignment_input.rs#L540)**
**Categories:** errors (+ module_structure, unsafe_concurrency, defaults cross-notes) · **Confidence:** High
`query_interval` does `Position::new(iv.start as usize).expect("region start >= 1")`. `Position::new` returns `None` only for `0`. The `>= 1` invariant holds for any `RegionSet`-derived interval, but `ContigInterval` is a `pub struct` with `pub start`/`pub end` and no validating constructor — any in-crate (today) or whole-crate (tomorrow) caller can build `ContigInterval { start: 0, .. }` and turn an input-shape problem into a process abort on the hot path. Related: the whole-contig collapse uses `iv.start <= 1` (not `== 1`), tacitly conceding `0` is reachable.
**Fix:** give `ContigInterval` a checked constructor returning `Option`/`Result` and make the fields non-`pub` so `start == 0` is unrepresentable; or at minimum add the `// PANIC-FREE:` marker documenting that this is a *construction* invariant, not a type one.

**M3: `BedError::Io` drops path/line and flattens its source — [src/regions.rs:296](../../../../src/regions.rs#L296)**
**Categories:** errors · **Confidence:** High
`#[error("reading BED file: {0}")] Io(io::Error)` is returned on a mid-file read failure (`line.map_err(BedError::Io)?`) with no path, no line number, and no `#[source]` (a bare tuple field gets neither auto-`#[source]` nor `#[from]`), so `source()` is `None` and the cause is only visible because `{0}` stringifies it into this error's own `Display` — breaking chain rendering and losing provenance.
**Fix:**
```rust
#[error("reading BED file {path}")]
Read { path: PathBuf, line_number: usize, #[source] source: io::Error },
```
Thread `path`/`line_number` from `from_bed_reader`'s loop.

**M4: FASTA-repository construction duplicated three times (and cross-file validation copied) — [src/bam/alignment_input.rs:887](../../../../src/bam/alignment_input.rs#L887)**
**Categories:** smells · **Confidence:** High
The ~12-line "`with_fai_extension` → `MissingFastaIndex` check → `indexed_reader::Builder` → `IndexedReader`/`Repository`" sequence appears verbatim in `new`, `query`, and `load_pileup_inputs`. Separately, the contig/sample cross-file reconciliation (`first_disagreement`→`ContigListMismatch`, `MultipleSampleNames`, the two `reference_input_path.clone().expect(...)`) is a near-verbatim second copy of the block in `AlignmentMergedReader::new` — the `PileupInputs` doc even says invariants are "enforced exactly as `new` does", which the two copies must keep true by hand.
**Fix:** extract `fn open_fasta_repository(fasta: &Path) -> Result<fasta::Repository, _>` and call it from all three sites; extract a small `CrossFileValidator { accept(path, &ExtractedHeader)?; finish() }` fed by both `new` and `load_pileup_inputs`.

**M5: `--build-map-file-index` is an unrelated feature folded into this PR — [src/pop_var_caller/cli.rs:97](../../../../src/pop_var_caller/cli.rs#L97)**
**Categories:** extras (+ naming) · **Confidence:** Medium (depends on Q1)
The pileup diff adds a second new flag with its own `PileupArgs` field, help text, and index pre-flight path. It's adjacent (seek-based region access needs an index) but independently testable and changes default missing-index behaviour. A reviewer scanning for "BED region" changes won't scrutinise the index-build half.
**Fix:** split it out, or state in the commit/PR description that it's part of the regions work so the coupling is explicit. (No code change if the latter.)

**M6: whole-contig range collapse is applied silently with no runtime trace — [src/bam/alignment_input.rs:880](../../../../src/bam/alignment_input.rs#L880)**
**Categories:** defaults · **Confidence:** Medium
`query` normalises a caller range to the unbounded fast path when `iv.start <= 1 && u64::from(iv.end) >= contig_length` — a behaviourally-significant fallback (it changes whether the per-record overlap filter runs) applied with no `debug`-level signal. An operator asking "why did my narrow BED span read the whole contig?" has only the source to consult.
**Fix:** emit a `tracing::debug!` naming the contig + the values that triggered the collapse (downgrade to Minor/doc-only if `tracing` isn't wired on this path).

### Minor

- **Mi1: per-region repeated work undercuts "a sparse BED is cheap."** [src/pop_var_caller/cli.rs:312](../../../../src/pop_var_caller/cli.rs#L312) (smells, idiomatic, reliability, extras). The sequential region loop calls `query` per `Region`, passing `inputs.contigs.clone()` + `inputs.sample_name.clone()` each iteration, and `query` rebuilds the whole FASTA `Repository` every call. The "per-worker, not shareable across threads" justification doesn't apply to this serial reuse. A 10k-region BED pays 10k FASTA-index opens + `ContigList` deep-clones. **Fix:** build the repository once and thread it in; have `query` borrow `&ContigList`/`&str` (or store `Arc` in `PileupInputs`). Verify byte-identical `.psp` after.
- **Mi2: `BedError::Open` double-renders its source.** [src/regions.rs:292](../../../../src/regions.rs#L292) (errors). `#[error("… {path}: {source}")]` on a field named `source` both nests (good) and concatenates the cause into `Display`. Drop the `: {source}` suffix.
- **Mi3: `#[from] AlignmentIndexError` collapses two origin sites.** [src/bam/errors.rs:200](../../../../src/bam/errors.rs#L200) (errors). Pre-flight vs per-input `load_alignment_index` both funnel through one variant; the phase is recoverable only from the inner variant. If a caller needs the phase, drop `#[from]` and map explicitly (Q2).
- **Mi4: `build_map_file_index` is a third name for "alignment file/index".** [src/pop_var_caller/cli.rs:101](../../../../src/pop_var_caller/cli.rs#L101) (naming). The crate uses *alignment file*/*alignment index* everywhere (`AlignmentIndex`, `load_alignment_index`, `alignment_files`); "map file" appears only here, and `load_pileup_inputs`'s param is yet a third spelling (`build_index_if_missing`). Reconcile the field/param to `build_alignment_index`; keep a flag alias if CLI-compat matters.
- **Mi5: `stashed_upstream` drops its noun.** [src/pop_var_caller/cli.rs:301](../../../../src/pop_var_caller/cli.rs#L301) (naming). The long-lived local mirrors `Stage1Outputs::stashed_upstream_error`; name it `stashed_upstream_error`.
- **Mi6: region-provenance asymmetry.** (reliability, refactor_safety cross-notes). The pileup `.psp` header records `regions_bed`/`regions_count`, but `run_var_calling` resolves `--regions` without recording an equivalent marker in the VCF (Q4).
- **Mi7: an empty resolved `RegionSet` silently yields empty output.** [src/pop_var_caller/cli.rs](../../../../src/pop_var_caller/cli.rs) (reliability). A comments-only BED parses to zero spans → empty `.psp`/VCF with a clean exit and no warning (an *unknown contig* still errors loudly). Pin the contract with a test; consider a stderr warning (Q3).
- **Mi8: `command_line` `""` is an undocumented back-compat sentinel.** [src/psp/header.rs:130](../../../../src/psp/header.rs#L130) (defaults). `"" == produced before this field existed` lives only in prose; a future writer leaving it empty would masquerade as legacy. Note in the writer-side docs that new headers must populate it.
- **Mi9: `#[allow(clippy::too_many_arguments)]` without the required justification comment.** [src/pop_var_caller/cli.rs:458](../../../../src/pop_var_caller/cli.rs#L458) (`build_writer_header`, now 8 positional args incl. a run of bare primitives) and [src/pop_var_caller/stage1_pipeline.rs:90](../../../../src/pop_var_caller/stage1_pipeline.rs#L90) (smells). Prefer grouping the layout/region knobs into a small struct so the `#[allow]` can be deleted; else add the lint/why/removal-condition comment.

### Nits
- `load_pileup_inputs`' `.expect(...)` calls ([src/bam/alignment_input.rs](../../../../src/bam/alignment_input.rs)) are sound (set in lockstep on iteration 1, empty-input rejected up front) but use prose instead of the `// PANIC-FREE:` marker; or collapse the two `Option`s into one `Option<(PathBuf, ContigList, String)>` set once.
- `RegionSet::iter` returns `std::slice::Iter<'_, Region>`, exposing the backing container; `impl Iterator<Item = &Region>` hides it. [src/regions.rs:165](../../../../src/regions.rs#L165)
- `var_calling.rs` spells `crate::regions::…` inline; a top-of-file `use` would match `cli.rs`. [src/pop_var_caller/var_calling.rs](../../../../src/pop_var_caller/var_calling.rs)
- `ContigInterval { start: region.start, end: region.end }` (cli.rs) is a manual span copy; `impl From<&Region> for ContigInterval` would name it (single call site, so optional).
- `cram_input.rs` block-scoped `iv` binding could be `interval` for symmetry with the rest of the file.

## 7. Out of scope observations
- The pre-existing walker **chrom-transition tail-drop bug** is dormant on this branch (pileup is single-contig per query) but lives on in the walker; the reliability agent notes a sharper interaction if it ever reactivates — with `--regions`, each region's walker reaches EOF at the region's last read, so a reactivated bug would drop a tail column *per region*, and `drive_region_into_writer`'s clamp can't re-add a column the walker never emitted. Tracked in the perf/byte-identity report; not introduced here.
- `current_command_line()` now lands verbatim in the `.psp` header; if argv can carry secrets, that's a provenance/privacy consideration for whoever owns that policy (out of this diff's lens).

## 8. Missing tests to add now

Grouped by function. (Full bodies in `tmp/review_2026-06-03_bed-regions/{reliability,extras}.md`.)

- `run_summary_merge_sums_counters_but_maxes_high_water` — two summaries where the second's `active_reads_high_water` is *smaller*; catches sum-instead-of-max and dropped additive fields. (`src/pileup/walker/driver.rs` tests)
- `filter_counts_merge_sums_every_field` / `baq_skip_counts_merge_sums_every_field` — distinct non-zero value per field; catch a dropped field. (`alignment_input.rs` / `baq_stream.rs` tests)
- `parse_writer_defaults_command_line_when_key_absent` — deserialize a `[writer]` table with no `command-line` key → `""`; catches a dropped `#[serde(default)]`. Optionally extend `wire_keys_are_pinned` to assert `command-line =` present / `command_line` forbidden. (`psp/header.rs` tests)
- `u64_max_start_errors_without_panicking` + `non_numeric_end_errors` — adversarial BED coordinate (`u64::MAX`) and non-numeric *end* column; pin `BedError::InvalidInterval`/`Parse` rather than a panic on the `bed_start + 1`. (`regions.rs` tests)
- `regions_empty_bed_yields_empty_output` (integration) — comments-only BED → valid empty `.psp`/VCF; locks the Mi7 contract.
- *(belt-and-suspenders)* `restrict_intervals_to_regions` with a covered interval overlapping ≥3 regions and an exact `c.end == region_end_excl` tie — the logic is verified correct, but these pin the boundary.

## 9. What's good
- `src/regions.rs` is a textbook cross-stage interchange peer: **zero `use crate::*`** (only `std`/`thiserror`), takes the neutral `ContigBounds`, consumed by both pipeline stages — correct top-level placement with no back-references.
- The var-calling region restriction preserves the **producer↔DUST-pool lockstep by construction**: narrowing mutates the single shared `chrom_plans` ([src/var_calling/driver.rs:977](../../../../src/var_calling/driver.rs#L977)), so the pool's job list and the producer's walk are projections of one sequence and cannot desync; the `None` path is a literal no-op.
- Whole-genome behaviour is **byte-identical by construction** in var-calling (the intersection is skipped entirely when `region_set` is `None`), not merely by test.
- The new `command_line` wire field is backward-compatible done right: `#[serde(rename = "command-line", default)]` + the derived-`PartialEq` round-trip test catches a forgotten construction site.
- Exhaustive `match` over the `#[non_exhaustive]` `AlignmentIndex` and `is_none_or`/`partition_point`/let-chains throughout — idiomatic and refactor-safe at the enum boundary.

## 10. Commands to re-verify
- `cargo fmt --check` · `cargo clippy --all-targets --all-features -- -D warnings` · `cargo doc --no-deps` · `cargo test --all-features` (all green above).
- After applying M1/§8: `cargo test --lib -- merge` to confirm the new counter-fold tests.

### Author response convention
Address each finding by id (B/M/Mi) with `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer §4 questions first.
