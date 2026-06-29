# Code Review: Stage-1 summary-accumulator wiring (C1+C2)

**Date:** 2026-06-29
**Reviewer:** rust-code-review skill (orchestrator + 1 consolidated 7-category sub-agent)
**Scope:** C1+C2 (merged) — `--gc-window-bp` flag + wiring the coverage/het accumulators into the Stage-1 pileup→psp seam
**Status:** Request-changes (Blocker) → all resolved

---

## 1. Scope
- In-scope: `src/pileup/per_sample/pileup_to_psp.rs` (`SampleSummaryAccumulators`
  bundle + `observe_record` + `drive_region_into_writer` signature),
  `src/pop_var_caller/cli.rs` (`--gc-window-bp`, accumulator construction,
  finalise+attach, `PileupCliError::SampleSummary`),
  `src/sample_summary/mod.rs` (default consts), the two integration tests.
- Out of scope: pre-existing broken examples/benches (`cohort_var_calling_perf`,
  `dhat_var_calling`, `profile_cohort_e2e`) that reference the removed
  `CohortPipelineArgs.min_mapq_diff_t` — unrelated to this diff.
- Categories: reliability, errors, idiomatic, naming, smells, defaults, refactor_safety.

## 2. Verdict
**Request-changes → resolved.** One Blocker (a real correctness bug in the
record→counts shaping), two Majors, two Minors. The error-handling shape,
het input mapping, and finalise ordering were judged correct.

## 3. Execution status
- `cargo test --lib` → 1392 passed (pre-fix). `cargo fmt --check` → exit 0.
- `cargo clippy --lib --tests --bins` → clean on scope. Integration → 22 passed.

## 4. Findings

### Blocker
- **B1 (reliability): `observe_record` double-counted coverage and tripped the
  coordinate-order invariant on multi-base records.** The span-iteration
  attributed a record's depth to *every* position in its REF span, but the
  walker emits a separate per-position record for each interior base, so a
  deletion/MNP anchor (span N at pos P) overlapped the records at P+1..P+N —
  double-counting and feeding the coverage accumulator a backwards position
  (its non-decreasing-`(chrom,pos)` `debug_assert!` panics in debug; silently
  re-enters a left tile in release). Dormant because the only `observe_record`
  test used span-1 `snp_read`. **Fix: anchor-only attribution** — one coverage
  observation per record at its anchor; the per-position records cover the rest.
  *(Applied + regression test `observe_record_handles_multibase_record_in_coordinate_order`.)*

### Major
- **M1 (errors/defaults): `--gc-window-bp 0` panicked in release** (no
  `value_parser`; the accumulator `assert!`s `window_bp >= 1`). *(Applied —
  `value_parser = clap::value_parser!(u32).range(1..)` so it is a clean clap
  error.)*
- **M2 (errors): the stashed-upstream-error path could mask the real error.**
  After the early `break`, the code still ran finish/serialise/attach/finalise
  with `?` on a `.psp` about to be deleted, so a secondary failure shadowed the
  upstream read error. *(Applied — surface the upstream error and remove the
  `.tmp` **before** any summary/finalise work.)*

### Minor
- **Mi1 (defaults): `--gc-window-bp` was the one CLI knob missing from
  `[writer.parameters]`.** *(Applied — threaded through `effective_parameters` /
  `build_writer_header`; the `every-knob` test now asserts it. Bin resolution +
  het ε/margin/min-depth remain recorded in the metadata-section document, noted
  in the doc.)*
- **Mi2 (reliability): `observe_record` index-panicked on an empty-allele
  record.** *(Applied — `record.alleles.first()` guard + `seq.first()` fallback.)*

### Nit
- `offset as u32` cast — **moot** after the anchor-only fix removed the span loop.

## 8. Missing tests added
- `observe_record_handles_multibase_record_in_coordinate_order` (B1).
- `effective_parameters_records_every_knob` now includes `gc_window_bp` (Mi1).
- End-to-end: `pileup_cli_integration::happy_path_default_config` now asserts the
  emitted `.psp` carries a parseable, valid `SampleSummary` with the CLI window.

## 9. What's good
- The seam owns record→counts *shaping*; the accumulators stay pure math
  (split-data-shaping-from-math).
- One accumulator bundle per run, fed exactly the written records, attached via
  A1 before finish — the stored statistics match the `.psp` body.
- Typed `SampleSummaryError` → `PileupCliError::SampleSummary` (`#[from]`).

## 10. Commands to re-verify
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo test --test pileup_cli_integration --test cohort_cli_integration`
- `./scripts/dev.sh cargo clippy --lib --tests --bins --all-features -- -D warnings -A clippy::doc_lazy_continuation`
