# `.psp` container generalization — step 3 (interval-overlap region query)

**Date:** 2026-06-15
**Branch:** `ssr-architecture`
**Plan:** [psp_container_generalization.md](../../../doc/devel/implementation_plans/psp_container_generalization.md)
(architecture §10.5 of [ssr_genotyping_architecture.md](../../../doc/devel/architecture/ssr_genotyping_architecture.md))
**Builds on:** [1a](psp_container_generalization_step1a_2026-06-15.md),
[1b](psp_container_generalization_step1b_2026-06-15.md),
[2](psp_container_generalization_step2_2026-06-15.md)

## Goal

Generalize the region-query overlap to **intervals** (§10.5), with SNP as a
degenerate point. Behaviour- *and* byte-preserving for SNP.

## What was already done (and so didn't change)

A read of the overlap machinery showed two of the three pieces were already
interval-shaped:

- **Block-index overlap** — `find_first_overlapping_block` already tests the
  block range `[first_pos, last_pos]` against the query `[start, end]`
  (`last_pos >= start && first_pos <= end`). Interval-correct for any schema's
  `last_pos`.
- **`last_pos` fill** — `BlockAccumulator::last_pos` is already the schema-owned
  "block max end". SNP records are points, so SNP's `last_pos` is the last
  record's start (= its degenerate end) — unchanged.

So the writer and the block index needed **no code change** (the produced `.psp`
is byte-identical to step 2). The one piece still point-shaped was the
**per-record region clamp** (the `record_coord` hook added in 1b).

## What changed

### [src/psp/kind.rs](../../../src/psp/kind.rs)
- `PspKind::record_coord(&Record) -> (u32, u32)` → `record_interval(&Record) ->
  (u32, u32, u32)` (`chrom_id, start, end`), `end` exclusive.

### [src/psp/writer.rs](../../../src/psp/writer.rs)
- `SnpKind::record_interval` returns `(chrom_id, pos, pos + 1)` — the degenerate
  point `[pos, pos + 1)` (read-side metadata only; does not touch written bytes).

### [src/psp/reader.rs](../../../src/psp/reader.rs)
- `RangeClamp` predicates are now interval-based: `interval_past_window(chrom,
  rec_start)` (`rec_start > q_end`, terminate) and `interval_before_window(rec_end)`
  (`rec_end <= q_start`, skip). `RecordsIter::next` clamps via
  `S::record_interval`.

### [src/psp/index.rs](../../../src/psp/index.rs)
- Doc-only: pinned the `BlockIndexEntry::last_pos` semantics ("block max
  inclusive end; SNP = last record start, SSR = `max(record.end)`").

## Why SNP behaviour is identical

For a SNP point at `pos`, `record_interval = (chrom, pos, pos + 1)`:
- past window: `rec_start > q_end` ⟺ `pos > q_end` (was `pos > end`). Same.
- before window: `rec_end <= q_start` ⟺ `pos + 1 <= q_start` ⟺ `pos < q_start`
  (was `pos < start`). Same.

So the set of records `region_records` yields is unchanged — confirmed by the
reader's `region_records` unit tests and the `region_flag_clamps_to_one_record`
e2e test.

## Gates (dev container)

- `cargo build --lib` / `clippy --all-targets -- -D warnings` / `fmt --check` — clean
- `cargo test --lib` — **1136 passed, 0 failed, 1 ignored**
- SNP e2e: `psp_to_pileup_integration` (incl. region clamp), `pileup_cli_integration` — green
- Produced `.psp` is byte-identical to step 2 (no writer/index encoder change).

## Next (remaining §10.6 steps)

- **4** — `registry_ssr` + `SsrLocusRecord` (all-CSR: `amb_*` CSR + QC scalars;
  `SsrBlock`/`SsrDecoder`/`SsrKind`). This is where `record_interval` returns a
  real locus interval and `last_pos` becomes `max(record.end)`.
- **5** — synthetic `.ssr.psp` round-trip test.
