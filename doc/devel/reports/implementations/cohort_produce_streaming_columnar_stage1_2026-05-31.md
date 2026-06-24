# Stage 1 — span-addressable columnar PSP reader (2026-05-31)

Implementation report for Stage 1 of
[cohort_produce_streaming_columnar.md](../../../doc/devel/implementation_plans/cohort_produce_streaming_columnar.md):
the span-addressable, columnar reader that lets the cohort producer ask
a `.psp` for *genomic spans* (not blocks) and receive `SampleColumns`
directly, with no row-shape `PileupRecord` round-trip.

**Additive only** — nothing is wired into the production loader yet
(that is Stage 2). The reader ships as a unit-tested building block.

## Plan

1. **psp side:** expose a columnar block decode — a borrowed
   `BlockColumns` view + a `BlockColumnReader` cursor that walks/decodes
   blocks on demand (reusing the existing `decode_block_payload` and
   block-index helpers), with `peek_block` (free, index-only) and
   `load_current`/`advance`/`columns`.
2. **cohort_block side:** `SampleColumns::append_block_window` — the
   columnar→columnar transcode (positions prefix-sum, bulk per-allele
   fixed-column copies, CSR slices for the ragged columns) replacing the
   `push_record` row path.
3. **cohort_block side:** `ColumnSpanReader` — wraps `BlockColumnReader`,
   tracks a per-region cursor, and serves `peek_next_span()` +
   `read_span(end, &mut SampleColumns)`, buffering any partially-served
   block via the held cursor (the "tail").
4. Tests: equivalence vs the row path over whole / split / clamped /
   empty spans; a misaligned-PSP fixture driving the consensus-span
   loop two samples at a time.

## Assumptions

- The decoded block's per-allele column order matches `PerAlleleFixed`'s
  fields 1:1 (verified by reading both: `allele_obs_count`→`num_obs`,
  `allele_q_sum_log`→`q_sum`, …), so the fixed-column append is a
  straight per-column `extend_from_slice`.
- Stage 1 carries **no production behaviour change** (not wired in), so
  correctness is proven by the equivalence unit tests rather than an
  end-to-end VCF md5. The end-to-end md5 *was* re-checked anyway because
  unrelated clippy cleanups touched `produce_block` (see Validation).
- Position delta-varint overflow is surfaced with the same typed error
  `materialise_next_record` raises (`ColumnElementDecode` /
  `VarintOverflow`), decoded by the caller (`ColumnSpanReader`) during
  its windowing scan rather than inside the column append.

## Changes made

- **[src/psp/reader.rs](../../../src/psp/reader.rs):** new public
  `BlockColumns<'a>` borrowed view + `BlockColumnReader<'r, R>` cursor
  (`index`/`cur_block_idx`/`seek_to`/`peek_block`/`load_current`/
  `advance`/`columns`), and `PspReader::column_blocks()`. Reuses the
  existing decode internals; `RecordsIter` (the row path) is untouched.
- **[src/psp/mod.rs](../../../src/psp/mod.rs):** export `BlockColumns`,
  `BlockColumnReader`.
- **[columns.rs](../../../src/var_calling/cohort_block/columns.rs):**
  `SampleColumns::append_block_window(cols, record_lo, record_hi,
  allele_lo, abs_positions)` — the columnar transcode (bulk position +
  per-allele-fixed copies, CSR pushes for seq/chain-ids, running
  per-record offsets), no `PileupRecord`.
- **[column_span_reader.rs](../../../src/var_calling/cohort_block/column_span_reader.rs)**
  (new): `ColumnSpanReader` — `new`/`reset`/`peek_next_span`/`read_span`.
  Hides PSP block boundaries; a `read_span` ending inside a block
  resumes from the held cursor next call. `+ mod.rs` registration/export.
- **[driver.rs](../../../src/var_calling/cohort_block/driver.rs):**
  clippy cleanups in code introduced by the earlier Stage 1/2 commits,
  surfaced now: removed the vestigial single-pass `loop` in `next_block`
  (`never_loop`, a deny-by-default error), collapsed the DUST `if let …
  { if … }` into a let-chain, reworded a doc comment whose wrapped `+`
  line read as a markdown list. Behaviour-preserving (byte-identical,
  see Validation).

## Tests added

All in `column_span_reader`'s test module:

- `whole_region_one_span_matches_record_path` — multi-block file, one
  unbounded `read_span` equals `region_records` → `push_record`.
- `split_spans_including_midblock_match_record_path` — same, served
  across several `read_span` calls whose ends land inside blocks.
- `clamped_subwindow_matches_record_path` — a `[start, end]` sub-window
  matches the row path and really drops out-of-window records.
- `empty_region_yields_no_records` — a gap between records serves
  nothing and decodes nothing.
- `peek_reports_servable_ends_and_terminates` — `peek_next_span` yields
  strictly-ascending block ends, terminates, and the block-by-block
  drive reconstructs the file.
- `misaligned_blocks_served_at_consensus_span_match` — the **misaligned
  fixture**: the same records written with two different block targets
  (different boundaries) both reconstruct exactly when driven at a
  shared `W = min(peek)` — the Stage 2 cross-sample scenario.

## Validation results

- `cargo test --lib` (container): **1055 passed; 0 failed; 1 ignored**
  (was 1049; +6 new). `cargo test --test cohort_cli_integration`: 20/20.
- `cargo fmt`: applied.
- `cargo clippy --lib` / `--tests`: clean for all touched code. One
  pre-existing warning remains — `needless_range_loop` at
  [dust_filter.rs:626](../../../src/var_calling/dust_filter.rs#L626),
  unrelated to this change (the loop variable is also used in genomic-
  position arithmetic, so it's a borderline false positive); left
  untouched to avoid churning a byte-identity-critical DUST hot loop.
- **End-to-end byte-identity** (header-stripped md5, real 26-sample
  tomato cohort, 8 threads): N=8 and N=26 both identical to the
  committed baseline — confirming the `driver.rs` clippy cleanups did
  not change output.

## Tradeoffs and follow-ups

- **Layering:** the columnar block decode is exposed *publicly* from
  `psp` (`BlockColumns`/`BlockColumnReader`) so that `ColumnSpanReader`
  (which must produce `cohort_block`'s `SampleColumns`) can live in
  `cohort_block` without `psp` back-referencing it. The plan glossed
  this; the resolution is a clean `cohort_block → psp` dependency via a
  borrowed view.
- **Not the happy-path memcpy:** positions are delta-encoded and the
  ragged columns are CSR, so the append is a cheap transcode, not a raw
  `memcpy`. It still eliminates every per-record/per-allele `Vec` of the
  row path — the measured waste.
- **Deferred to Stage 2:** wiring the loader to consume `ColumnSpanReader`
  via the consensus-span streaming loop (the memory fix), which also
  deletes the safe-cut / `carryover` / `NoSafeGap` machinery.
- **Deferred to Stage 3:** DUST-ahead queue.
