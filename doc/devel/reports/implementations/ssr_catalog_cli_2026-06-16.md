# SSR Stage 0 (`ssr-catalog`) — the `ssr-catalog` CLI subcommand

**Date:** 2026-06-16
**Branch:** `ssr-architecture`
**Plan:** [ssr_catalog.md](../../../doc/devel/implementation_plans/ssr_catalog.md) §2.
**Scope:** the fifth Stage-0 increment — wire the `run()` orchestrator behind a
`pop_var_caller ssr-catalog` subcommand. Makes Stage 0 runnable end-to-end from
the command line.

---

## What landed

- New [`src/pop_var_caller/ssr_catalog.rs`](../../../src/pop_var_caller/ssr_catalog.rs):
  - `SsrCatalogArgs` (clap) — `--reference`, `--output`, `--trf-mod-path`,
    `--temp-dir` (default `ssr-catalog-tmp`, CWD-relative + disk-backed), and the
    Advanced knobs `--flank-bp` / `--bundle-threshold` / `--min-purity` /
    `--min-score` (defaults pinned from `CatalogParams::default`: 50 / 50 / 0.8 / 0).
  - `run_ssr_catalog(args)` — validates `bundle_threshold ≥ flank_bp`, stamps the
    build date (`rfc3339_now()` date part — the library reads no clock), maps the
    args onto `CatalogConfig`, and calls `catalog::run`.
  - `SsrCatalogCliError` (`#[non_exhaustive]`): `BundleThresholdTooSmall` and
    `Catalog { source: Box<dyn Error + Send + Sync> }` — the catalog cause is
    **boxed** so this public error does not name the crate-internal `CatalogError`
    (avoids the `private_interfaces` lint while preserving the `source()` chain
    that `main::format_error_chain` walks).
- Wiring: `SsrCatalog(SsrCatalogArgs)` variant on `PopVarCallerCommand`
  ([cli.rs](../../../src/pop_var_caller/cli.rs)); the dispatch arm in
  [main.rs](../../../src/main.rs); the re-export in
  [pop_var_caller/mod.rs](../../../src/pop_var_caller/mod.rs).

## End-to-end verification (in the dev container)

`pop_var_caller ssr-catalog --help` renders the full option list. Running it on a
synthetic `chr1` (poly-T flanks + `CAG×30` + poly-T + non-repetitive + `AT×25`)
produced a correct catalog:

```
## tool: ssr-catalog
## tool_version: 0.1.0
## reference_md5: cb228f6e8c0fc3171d7e1fce7fd3ab42
## trf_mod_version: Tandem Repeats Finder, Version 4.10.0
## flank_bp: 50  ## min_purity: 0.8  ## min_score: 0  ## bundle_threshold: 50
## date: 2026-06-16
#chrom  start  end  motif  purity_fraction  ref_seq_start  ref_seq
chr1    80   170  CAG  1  30   …
chr1    291  341  AT   1  241  …
```

— both tracts at exact coordinates, purity 1.0, 50 bp flanks embedded, and **no
period-1 homopolymer loci** (the poly-T runs were dropped before bundling, so
they did not bundle-drop the CAG/AT tracts).

## Tests

`rejects_bundle_threshold_below_flank` (the CLI invariant fires before any I/O)
and `parses_through_the_top_level_cli` (the subcommand + defaults parse through
`Cli`). The orchestrator's own end-to-end build is covered by the lib-level
`run_builds_a_catalog_for_a_synthetic_reference` (catalog `mod.rs`).

## Stage-0 status after this increment

The `ssr-catalog` tool is **usable end-to-end** (sequential). Remaining Stage-0
follow-ups, both optimizations/extensions rather than blockers:
- **Per-contig rayon fan-out** (architecture §8) — bounded pool sized by a
  `--num-chroms-in-parallel` arg; ordered collector (contig order = sorted).
- **`write_index` (CSI)** — the coordinate index for the catalog reader's
  `--regions` query path (consumed later by Stage 1's `fetch_reads`).

## Validation (dev container)

- `cargo fmt --check` → clean.
- `cargo clippy --all-targets --all-features -- -D warnings` → clean.
- `cargo test --lib` → all green (+2 CLI tests).
- `cargo doc --no-deps` → clean.
- Manual: `ssr-catalog --help` + a synthetic-reference build, output verified.
