# SSR Stage 0 (`ssr-catalog`) — format I/O layer (`catalog/io.rs`)

**Date:** 2026-06-15
**Branch:** `ssr-architecture`
**Plan:** [ssr_catalog.md](../../../doc/devel/implementation_plans/ssr_catalog.md) §5/§7;
architecture [ssr_catalog.md](../../../doc/devel/architecture/ssr_catalog.md) §7.
**Scope:** the first increment of Stage 0 — the catalog **wire format** (the
cross-stage contract Stage 1's `fetch_reads` reader consumes). The detection
front-end (`trf`, `postprocess`) and the `run()` orchestrator + CLI are deferred
to later increments.

---

## What landed

New module `src/ssr/catalog/`:

- [`mod.rs`](../../../src/ssr/catalog/mod.rs) — module doc (incremental build
  status) + `CatalogError` (`#[non_exhaustive]`; today: `Io`, `HeaderParse`,
  `RowParse`, `InvalidMotif`, `InvalidLocus` — the `trf-mod`-spawn / FASTA
  variants come with `trf.rs`).
- [`io.rs`](../../../src/ssr/catalog/io.rs) — the format:
  - `CatalogParams` (recorded build knobs: `flank_bp`, `min_purity`,
    `min_score`, `bundle_threshold`) and `CatalogHeader` (reference +
    `reference_md5`, `trf_mod_version`, params, `tool_version`, `date`).
  - `CatalogWriter<W>` — wraps the sink in a **bgzip** frame
    (`noodles-bgzf`), emits the `##` metadata block + the `#`-column header on
    `new`, one row per `write_locus`, finalises on `finish`.
  - `CatalogReader<R>` — decodes the bgzip frame, parses the header on `new`,
    yields loci via `read_locus` / `read_all`.
  - private `locus_to_row` / `row_to_locus` — the single shared row codec
    (`chrom  start  end  motif  purity_fraction  ref_seq_start  ref_seq`), so
    write and read can't drift. `Motif::new` + `Locus::new` enforce the
    invariants on parse.

Wired `pub mod catalog;` into [`src/ssr/mod.rs`](../../../src/ssr/mod.rs) (the
module-level `#![allow(dead_code)]` covers the not-yet-consumed surface).

## Design notes

- **Format matches the architecture doc verbatim** (§7): self-describing
  bgzip TSV, `##` metadata + `#` column header (both skipped by tabix/CSI),
  7 columns, everything else derived. The column-header line is written
  verbatim and **validated on read**, so an upstream column reshuffle fails
  loudly.
- **No clock read** — `date` is a caller-supplied field, so the writer is a
  pure function of its inputs (deterministic catalog bytes).
- **f32 purity round-trips bit-for-bit** — Rust's shortest-round-trip `f32`
  formatting + `parse::<f32>()`; pinned by a dedicated test.
- **Visibility:** the whole catalog API is `pub(crate)`, matching the
  `Locus`/`Motif` visibility and the SSR-module convention (no cross-crate
  consumer; avoids the `private_interfaces` lint).

## Deferred within Stage 0 (next increments)

- **`write_index` / `CatalogReader::query`** — the CSI coordinate index for the
  `--regions` path (lean `.csi`, no new dep). The sequential `read_locus` path
  is here; region query is the next I/O piece.
- **`trf.rs`** — locate + spawn `trf-mod` per contig + BED parse.
  **Blocked on tooling:** the `trf-mod` binary is not installed on the host or
  in the dev container (only the `TRF-mod/` source is vendored), so this
  increment cannot be built/tested until trf-mod is compiled into the
  container image.
- **`postprocess.rs`** — the GangSTR-port pipeline (period≤6 → drop-compound →
  drop-bundle → end-trim → recompute-purity → embed-`ref_seq`). Pure logic, no
  external dep — buildable next regardless of the trf-mod blocker.
- **`run()` orchestrator + `ssr-catalog` CLI subcommand.**

## Tests (6, all green)

In `io.rs`:
- `locus_row_round_trips` — `row_to_locus(locus_to_row(l)) == l`; pins the 7-column shape.
- `imperfect_purity_round_trips_exactly` — f32 purity survives the text round-trip.
- `writer_reader_round_trip_through_bgzf` — multi-locus write→bgzf→read equality + header equality.
- `header_missing_required_key_is_rejected` — a dropped `reference_md5` → `HeaderParse`.
- `row_with_wrong_column_count_is_rejected` — `RowParse`.
- `row_with_bad_motif_is_rejected` — empty motif → `InvalidMotif`.

## Validation (dev container)

- `cargo fmt --check` → clean.
- `cargo clippy --all-targets --all-features -- -D warnings` → clean (0 warnings).
- `cargo test --lib` → `1148 passed; 0 failed; 1 ignored` (+6 catalog tests over the prior 1142).
- `cargo doc --no-deps` → clean.
- `cargo audit` → not run (cargo-audit unavailable in-container; no deps added).
