# SSR Stage 0 (`ssr-catalog`) — trf-mod spawn + BED parse (`catalog/trf.rs`)

**Date:** 2026-06-15
**Branch:** `ssr-architecture`
**Plan:** [ssr_catalog.md](../../../doc/devel/implementation_plans/ssr_catalog.md) §3;
architecture [ssr_catalog.md](../../../doc/devel/architecture/ssr_catalog.md) §2-3.
**Scope:** the third Stage-0 increment — locate, version, run, and parse the
external `trf-mod` detector. Now testable end-to-end: `trf-mod` was installed in
the dev container in the preceding `build(container)` commit.

---

## What landed (in [`catalog/trf.rs`](../../../src/ssr/catalog/trf.rs))

- `locate_trf_mod(override)` — architecture §2.4 layered discovery: explicit
  override → a copy beside our own executable (`current_exe().parent()`) → `PATH`
  → [`CatalogError::TrfModNotFound`].
- `version(bin)` — runs `trf-mod -v` and extracts the `"… Version 4.10.0"` line
  (searches both stdout and stderr) for the catalog header's `trf_mod_version`.
- `run_on_contig(bin, name, seq, temp_root)` — per-contig run via a unique
  `tempfile` dir (collision-free across parallel contigs): writes the contig to
  `input.fa`, runs `trf-mod input.fa` with stdout redirected to `output.bed` and
  stderr discarded, waits, and treats a non-zero exit as
  [`CatalogError::TrfRun`]. **No stdin/stdout pipes** (architecture §3) → no
  deadlock, clean exit-status check. Parses every non-blank BED line.
- `parse_bed_line(line, expect_ctg, line_no)` — asserts the 10-column layout and
  that the contig column matches, then extracts `start`/`end` (0-based
  half-open), `period`, `fracMatch`, `score`, `pattern`; `copyNum`/`fracGap`/
  `entropy` are dropped. A malformed line is a hard [`CatalogError::TrfParse`]
  (catches an upstream format change, §3).
- New `CatalogError` variants: `TrfModNotFound`, `TrfSpawn`, `TrfRun`,
  `TrfVersion`, `TrfParse`.

## Format pinned against the source

The 10-column BED and coordinate convention were read from trf-mod's
`trfrun.h::trf_print_bed`:
`fprintf("%.*s\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%d\t%.2f\t%s", name, first-1, last,
period, copies, matches·.01, indels·.01, score, entropy, pattern)`.
So column 1 = `first - 1` (0-based start), column 2 = `last` (exclusive end),
and `[start, end)` is 0-based half-open — matching `TrfRecord` and the plan.
Version verified: `trf-mod -v` → "Tandem Repeats Finder, Version 4.10.0"
(lh3 fork); `--version` is unrecognised.

## Tests (4, all green)

- `parse_bed_line_extracts_fields` — a real row from
  `trf-mod TRF-mod/t/small_test.fasta` (`0  0  35  7  5.00  1.00  0.00  70  1.95
  TCATCGG`).
- `parse_bed_line_rejects_wrong_column_count` / `..._contig_name_mismatch` —
  the hard-assert paths.
- `run_on_contig_detects_a_synthetic_repeat` — **integration**: runs the
  installed `trf-mod` on 60×`CAG` and asserts it returns a period-3 tract. Runs
  in the dev container (confirmed not skipped); skips gracefully with a printed
  note if `trf-mod` is absent (e.g. a host `cargo test`).

## Deferred (the last Stage-0 piece)

- `run()` orchestrator + the `ssr-catalog` CLI: per-contig fan-out/collect
  (rayon, order-preserving) reading the reference (noodles-fasta), calling
  `run_on_contig` → `postprocess::build_loci` → `CatalogWriter`, computing
  `reference_md5`, building the `CatalogHeader`, and writing the CSI index. Pins
  the `flank_bp` / `min_score` / `bundle_threshold` defaults.

## Validation (dev container)

- `cargo fmt --check` → clean.
- `cargo clippy --all-targets --all-features -- -D warnings` → clean.
- `cargo test --lib` → all green (+4 trf tests).
- `cargo doc --no-deps` → clean.
