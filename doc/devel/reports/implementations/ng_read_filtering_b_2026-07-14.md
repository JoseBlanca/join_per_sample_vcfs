# Implementation report: ng read filtering — Milestone B (the two-phase cascade)

**Date:** 2026-07-14
**Feature:** ng step 1 — read filtering
**Plan:** [read_filtering.md](../../ng/impl_plan/read_filtering.md) (Milestone B, steps B1–B2)
**Spec / arch:** [spec §3/§5](../../ng/spec/read_filtering.md), [arch §1/§3](../../ng/arch/read_filtering.md)

## 1. Plan

Implement the pure decision half of read filtering: `verdict_pre_decode` (#1–#6,
flag/MAPQ) and `verdict_post_decode` (#7 too-short, #9 bad-CIGAR, #8
high-mismatch), unit-tested in isolation against `InMemoryRefSeq`. B1 and B2 are
the two halves of one cascade and land together.

## 2. Assumptions / decisions (recorded, not silent)

- **Post-decode order = too-short → bad-CIGAR → high-mismatch** (#7 → #9 → #8), a
  deliberate reorder of the spec's #7/#8/#9 *table* so the one
  reference-touching filter runs last. Honors the spec §3 "cheapest-first"
  principle, leaves the keep/drop *set* unchanged (independent filters), and
  charges a both-failing read to the root cause (`BadCigar`) over the symptom
  (`HighMismatchFraction`). Discussed and agreed with the owner before coding.
- **#8 runs on the original (un-left-aligned) CIGAR** — left-alignment is
  deferred to `pileup/` (spec §6); a legal left-shift only re-pairs *equal*
  reference bases, so the mismatch tally is invariant under it.
- **`flag: u16`** for `verdict_pre_decode` (not a `Flags` newtype) — matches
  `MappedRead.flag` and the reused `FLAG_*` constants (arch §6 open item, resolved).
- **Error model — fatal** (owner-confirmed, 2026-07-14): a `RefSeqError` from #8's
  fetch, including `OutOfBounds` (read window past contig end), aborts the run
  rather than degrading to a per-read skip/keep — a validly-aligned read cannot
  cover positions the contig lacks, so it signals a malformed record.
- **`ref_buf: &mut Vec<u8>` scratch parameter** added to `verdict_post_decode`
  (a small deviation from the spec's illustrative signature) so the reference
  read reuses one buffer across reads — no per-read allocation. The iterator (D)
  owns the buffer.

## 3. Changes made

- **`src/ng/read/filtering.rs`** — `verdict_pre_decode(flag, mapq, config)`
  (mirrors production `classify_pre_decode` order exactly) and
  `verdict_post_decode(read, reference, config, ref_buf) -> Result<_, RefSeqError>`.
  Both reuse the production predicates (`cigar_is_bad`, `cigar_ref_span`,
  `read_exceeds_mismatch_fraction`) and `FLAG_*` constants unchanged; #8 reads
  raw reference bytes via `RawRefSeq::fetch_raw_into`. Marked
  `#[cfg_attr(not(test), allow(dead_code))]` until the Milestone-D iterator wires
  them. Guarded `u32` conversions for the reference coordinates (fail loudly, not
  silently truncate).

## 4. Tests added (12)

- Pre-decode: clean-keep; low-MAPQ boundary (kept at threshold, dropped one
  below, unavailable→0); no-minimum; each flag bit → its bucket; toggle-off
  behavior; first-firing attribution; **full cascade order**.
- Post-decode (against `InMemoryRefSeq`): too-short boundary; both bad-CIGAR
  shapes; high-mismatch boundary (exclusive); BQ-floor exclusion; #8-disabled
  makes no reference access (and enabling proves a fetch is attempted);
  **`OutOfBounds` is fatal**; #9-before-#8 and #7-before-#9 attribution;
  zero-ref-span keep.

## 5. Validation

Dev container: `cargo fmt -- --check` (ng clean), `cargo clippy --lib` (clean),
`cargo test --lib -- ng::read::filtering` → **18 tests pass**.

## 6. Tradeoffs and follow-ups

- **Deferred to Milestone D:** the `DropReason`↔`ReadFilterCounts` exhaustive-match
  enforcement (belongs at the tally site) and the reused-buffer no-stale
  assertion (plan D3).
- The `verdict_*` fns are unused until D wires them into the `ReadFilter`
  iterator (dead-code-allowed in the meantime).
