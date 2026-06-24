# SSR Stage 0 (`ssr-catalog`) — the `run()` orchestrator (sequential)

**Date:** 2026-06-15
**Branch:** `ssr-architecture`
**Plan:** [ssr_catalog.md](../../../doc/devel/implementation_plans/ssr_catalog.md) §2;
architecture [ssr_catalog.md](../../../doc/devel/architecture/ssr_catalog.md) §2/§5/§8.
**Scope:** the fourth Stage-0 increment — the `run()` orchestrator that wires the
three building blocks (trf / postprocess / io) into an end-to-end catalog build.
**Sequential** for now; the per-contig rayon fan-out and the clap CLI are the
remaining follow-ups.

---

## What landed (in [`catalog/mod.rs`](../../../src/ssr/catalog/mod.rs))

- **`CatalogParams` unified** — the duplicate `io::CatalogParams` and
  `postprocess::PostProcessParams` (identical 4-field structs) collapsed into one
  `CatalogParams` in `mod.rs`, used by both the header and the pipeline, with a
  `Default` impl pinning the Stage-0 defaults: `min_purity = 0.8`,
  `min_score = 0` (trf-mod's own `-s 30` does the base filtering),
  `flank_bp = 50`, `bundle_threshold = 50` (`= flank_bp`, satisfying the
  `bundle_threshold ≥ flank_bp` clean-flank invariant).
- **`CatalogConfig` + `run(cfg)`** — `run()`:
  1. locates trf-mod + reads its version;
  2. streams the reference with `noodles_fasta::io::Reader` (one contig resident
     at a time), folding each contig's upper-cased bases into a single
     whole-reference md5 (the spec's upper-cased-content convention,
     recomputable by Stages 1-2 for an integrity check; `format_md5_hex` reused
     from `pop_var_caller::common`);
  3. per contig: `trf::run_on_contig` → `postprocess::build_loci`, accumulating
     loci in reference (= coordinate) order;
  4. writes the `CatalogHeader` + rows via `CatalogWriter`, then `sync_all`.
- New `CatalogError::Fasta` variant for reference read failures.

## End-to-end test

`run_builds_a_catalog_for_a_synthetic_reference` writes a `(CAG)*40` reference
with non-repetitive flanks, runs the orchestrator (driving the **installed
trf-mod**), reads the catalog back, and asserts a perfect period-3 CAG-family
locus (phase-robust: motif ∈ {CAG,AGC,GCA}, ~120 bp, pure tiling) plus a
32-hex `reference_md5` and a `Version`-bearing `trf_mod_version`. Skips
gracefully if trf-mod is absent. 22 catalog tests total, all green.

## DESIGN FINDING — homopolymers bundle-drop adjacent SSRs — RESOLVED 2026-06-16

**Resolution (PM):** drop period-1 homopolymers, matching GangSTR/HipSTR — the
target is di/tri/tetra-nucleotide microsatellites, and homopolymers are
error-prone for STR genotyping. Applied in `postprocess.rs`: `MIN_PERIOD = 2`,
so the scope filter (step 1) keeps period `2..=6` and period-1 records never
reach `drop_bundles`, sparing homopolymer-adjacent SSRs. The e2e test now uses
poly-T flanks and asserts the CAG tract survives at `[100,220)` with no period-1
survivor; a `postprocess` unit test pins the same. The original analysis below
is retained for context.

---



The first version of the e2e test used `T*100` flanks; the build dropped the CAG
locus. **Root cause is real, not a test bug:** trf-mod reports a long
homopolymer run as a period-1 repeat, and because architecture §4 **keeps**
period-1 (it lists only `period ≤ 6`, with a copy-floor of 10 for period 1), that
homopolymer record enters `drop_bundles` and, being adjacent to the SSR (within
`bundle_threshold`), drops the whole cluster — including the real SSR.

GangSTR avoids this by **dropping period-1 homopolymers before bundling**
(`2_trim.sh` step 2: `awk '$4 != 1'`). Our port deliberately did **not** copy
that step, per arch §4. The consequence: in a real genome (where poly-A/poly-T
runs are common and often abut microsatellites), keeping period-1 could
bundle-drop a meaningful fraction of legitimate SSRs.

**This is a documented-architecture decision, so it was left unchanged here** —
surfaced as an open question for arch §4: either (a) drop period-1 before
bundling (match GangSTR; loses pure-homopolymer STRs from the catalog), or
(b) keep period-1 but exclude period-1 records from *forming/forcing* bundles,
or (c) accept the interaction. Tracked in `ssr_stage1_remaining.md`.

## Deferred (remaining Stage-0)

- **Per-contig rayon fan-out** (architecture §8) — the ordered collector over a
  bounded pool sized by `--num-chroms-in-parallel`; the current build is
  single-threaded (trivially deterministic; the contig order is already sorted).
- **`ssr-catalog` clap subcommand** — `SsrCatalogArgs` + `run_ssr_catalog` wired
  into `PopVarCallerCommand` (cli.rs) + the `main.rs` dispatch + a CLI error type;
  maps the args onto `CatalogConfig` and stamps the build `date`.
- **`write_index` (CSI)** for the catalog reader's `--regions` query path.

## Validation (dev container)

- `cargo fmt --check` → clean.
- `cargo clippy --all-targets --all-features -- -D warnings` → clean.
- `cargo test --lib` → all green.
- `cargo doc --no-deps` → clean.
