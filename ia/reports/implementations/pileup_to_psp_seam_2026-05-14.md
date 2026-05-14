# Stage 1 emission seam: pileup walker → psp writer

Implementation date: 2026-05-14. Follow-on to
[ia/feature_implementation_plans/pileup_pull_iterator.md](../../feature_implementation_plans/pileup_pull_iterator.md)
and its report
[ia/reports/implementations/pileup_pull_iterator_2026-05-14.md](pileup_pull_iterator_2026-05-14.md).

## Plan

Wire the pull-shape `PileupWalker` into `PspWriter::write_record` in a
single small module
[src/per_sample_caller/pileup_to_psp.rs](../../../src/per_sample_caller/pileup_to_psp.rs).
This is the in-tool "Stage 1 emission half" — the half that turns
per-position records into the on-disk `.psp` artefact downstream
stages consume. The CRAM-input → BAQ → walker plumbing is upstream of
this seam and lives in sibling slices.

No new walker/writer logic; just the cross-module driver and a typed
error type that wraps both sides.

## Assumptions / silent choices

- **Function shape: `drive_pileup_to_psp(walker, writer) -> Result<(W, RunSummary), …>`**
  — takes the constructed walker and writer by value, returns the
  writer's inner sink (post-`finish`) and the walker's cumulative
  summary. Caller-managed construction so tests, examples, and the
  future CLI can each build their own header / sink / config without
  the seam taking a position on any of those.
- **Combined error type [`PileupToPspError`]** with `#[from]`
  conversions for both `WalkerError` and `PspWriteError`. Keeps each
  side's context intact (operator triage can still tell *which* side
  faulted) while letting `?` propagate transparently through the
  driver.
- **First-error-wins semantics.** A walker error short-circuits
  before the writer sees the record; a writer error short-circuits
  before the next walker pull. The dropped state (open walker, open
  writer block) is discarded without finalising — there is no
  partial-`.psp` rollback. Documented on the function.
- **Module visibility upgrade for test fixtures.** `pileup::tests`
  and `psp::test_fixtures` were `mod` (private). The seam's test
  module reuses them, so both are now `pub(crate) mod` (still
  gated on `#[cfg(test)]`). This makes shared test fixtures
  available across the crate's `#[cfg(test)]` modules without
  duplication.

## Changes made

### New file: [src/per_sample_caller/pileup_to_psp.rs](../../../src/per_sample_caller/pileup_to_psp.rs)

- `pub enum PileupToPspError { Walker(WalkerError), Psp(PspWriteError) }`
  with thiserror-derived `Display` and `#[from]` conversions.
- `pub fn drive_pileup_to_psp<I, F, W>(walker, writer) -> Result<(W,
  RunSummary), PileupToPspError>` — the seam itself. ~10 lines:
  ```rust
  for item in walker.by_ref() {
      let record = item?;
      writer.write_record(&record)?;
  }
  let summary = walker.summary();
  let sink = writer.finish()?;
  Ok((sink, summary))
  ```
- Four `#[cfg(test)]` tests (see below).

### [src/per_sample_caller/mod.rs](../../../src/per_sample_caller/mod.rs)

- Registered `pub mod pileup_to_psp;` alongside the existing slices.

### [src/per_sample_caller/pileup/mod.rs](../../../src/per_sample_caller/pileup/mod.rs)

- Made `mod tests;` → `pub(crate) mod tests;` so the seam's test
  module can reuse `MockFasta` / `snp_read`.

### [src/per_sample_caller/psp/mod.rs](../../../src/per_sample_caller/psp/mod.rs)

- Made `mod test_fixtures;` → `pub(crate) mod test_fixtures;` so the
  seam's test module can reuse `writer_header`.

## Tests added

All four live in `pileup_to_psp::tests`:

- `drive_pileup_to_psp_passes_records_through_unmodified` — sanity
  check that the seam's iteration side neither drops nor mutates
  records. Compares an independent walker run to confirm parity.
- `writer_rejects_walker_output_when_chain_expires_on_same_record_as_final_reference`
  — **regression pin** for the finding below.
- `walker_error_surfaces_as_walker_variant` — typed error propagation
  from the walker side via `MockFasta` of length 2 against a 4-base
  read (trips `WalkerError::Fasta`).
- `writer_error_surfaces_as_psp_variant` — typed error propagation
  from the writer side via a `chrom_id` mismatch between the walker's
  fasta and the writer's header.

## Validation results

All commands inside the dev container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` —
  clean.
- `cargo test --tests --lib --all-features` — **all green**:
  - 478 unit tests (lib) — pass (+4 vs. previous baseline; the
    seam's four).
  - 88 integration tests across `tests/` — pass.

## Pre-existing finding surfaced by this work

**Title:** Walker emits `expired_chains` on the same record whose
alleles still reference the expiring slot; writer's
`apply_record_to_block` rejects it as
`PhaseChainMarkerInconsistencyKind::AlleleReferencesUnknownSlot`.

**Reproducer:** two SNP reads `r1` and `r2` both covering positions
1..5 of an `ACGTA` reference. Five `PileupRecord`s are emitted; the
last (`pos=5`) carries `expired_chains = [0, 1]` AND `REF.chain_slots
= [0, 1]`.

**Why the walker emits this.** A read contributes at every position
in `[alignment_start, alignment_end]`, including its final position.
The closure rule emits the record at the read's final position when
`walker_pos > alignment_end` — i.e., on the *same tick* that
`expire_passed_reads` releases the read's slot. The
slot-expiration lifecycle marks are then drained and stamped onto
the first record of that closure-step's batch, which is the very
record whose `allele.chain_slots` still references the slot. This is
intrinsic to the walker's tick ordering today and affects essentially
every realistic input.

**Why the writer rejects this.** `apply_record_to_block` applies
`expired_chains` to its running active-slot set *before* validating
each `allele.chain_slots` entry. So a record that simultaneously
expires slot S and references S in its alleles trips the unknown-slot
check immediately. Code:
[src/per_sample_caller/psp/writer.rs:523-574](../../../src/per_sample_caller/psp/writer.rs#L523-L574).

**Why this didn't surface earlier.** The pre-existing walker tests
check only the *aggregate* lifecycle-mark counts (`total_new`,
`total_expired`); they never push records through the writer. The
pre-existing writer tests build records by hand and never exercise
the `expired_chains` + `allele.chain_slots` overlap case. The seam's
roundtrip test is the first integration point where the two contracts
meet.

**Resolution options (out of scope for this seam — flagged for design
discussion):**

1. **Walker stamps expiration onto a *later* synthetic marker record.**
   The walker would emit an extra "marker" record at the position
   immediately after the last contribution, carrying `expired_chains`
   and no allele data. Adds output records; needs a representation
   for a marker-only record (or some other convention).
2. **Walker stamps expiration onto the *first* record of the *next*
   closure step.** Cheap to implement (just defer
   `drain_lifecycle_marks` to the next non-empty close batch), but
   reads that expire on the last tick before end-of-input have no
   "next batch" — they'd need to flow through `flush_chromosome_into`.
3. **Writer applies `expired_chains` *after* allele validation.** The
   spec excerpt in
   [src/per_sample_caller/pileup/mod.rs](../../../src/per_sample_caller/pileup/mod.rs)
   actually says "Stage 2's recv loop applies them before reading
   each record's allele observations" — which is what the writer
   does today. So this option would change the spec, not just the
   code.

I'd recommend option 2 (cheapest, most local) but the choice is
yours.

## Tradeoffs and follow-ups

- **Roundtrip-through-`PspReader` parity test is deferred.** Until
  the walker/writer contract is reconciled, the seam cannot produce
  a `.psp` artefact for any realistic input. The
  `drive_pileup_to_psp_passes_records_through_unmodified` test
  covers the seam's value-passing wiring; the
  `writer_rejects_walker_output_…` test pins the integration
  finding so it can't silently get worse.
- **No production caller yet.** The CLI, the multi-sample driver, and
  the BAQ-in-the-loop wiring all still need to be built. This seam
  is ready for whichever lands first.
