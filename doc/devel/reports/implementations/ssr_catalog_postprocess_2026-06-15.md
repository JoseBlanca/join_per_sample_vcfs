# SSR Stage 0 (`ssr-catalog`) — post-processing pipeline (`catalog/postprocess.rs`)

**Date:** 2026-06-15
**Branch:** `ssr-architecture`
**Plan:** [ssr_catalog.md](../../../doc/devel/implementation_plans/ssr_catalog.md) §4;
architecture [ssr_catalog.md](../../../doc/devel/architecture/ssr_catalog.md) §4.
**Scope:** the second Stage-0 increment — the per-contig pipeline that turns
parsed TRF-mod records into catalog `Locus`es. Pure logic, no external process.
Ported faithfully from the vendored GangSTR `create_ref_panel` scripts.

---

## What landed

- [`trf.rs`](../../../src/ssr/catalog/trf.rs) — the parsed-record type
  [`TrfRecord`](../../../src/ssr/catalog/trf.rs) (`start`/`end` 0-based half-open,
  `period`, `frac_match`, `score`, `pattern`) that `postprocess` consumes, plus a
  `#[cfg(test)]` constructor. The locate/spawn/BED-parse functions are the next
  increment (trf-mod is now in the container).
- [`postprocess.rs`](../../../src/ssr/catalog/postprocess.rs) —
  [`build_loci(recs, chrom, contig_seq, params)`](../../../src/ssr/catalog/postprocess.rs)
  and the step functions, in architecture §4 order:
  1. scope + score gate (`period ∈ 1..=6`, `score ≥ min_score`, in-bounds);
  2. `is_compound` drop (`ATAT = (AT)²`) — GangSTR `is_compound` + `count_motif`;
  3. `drop_bundles` (whole-cluster drop within `bundle_threshold`) — GangSTR
     `remove_bundles.py` streaming clustering;
  4. `minimal_trim` (clean whole-motif boundaries) + per-period copy-number
     floor `{1:10, 2:5, 3:4, 4:3, 5:3, 6:3}` — GangSTR `minimal_trim`;
  5. `recompute_purity` (phase-0 perfect-tiling match fraction, spec §3.2) +
     `min_purity` floor; then embed `ref_seq` (tract + `flank_bp` each side,
     clamped at contig ends, upper-cased).
- Wired `pub mod {postprocess, trf};` into [catalog/mod.rs](../../../src/ssr/catalog/mod.rs).

## Faithful port + the two deliberate adaptations

Source: [GangSTR/scripts/create_ref_panel/scripts/](../../../GangSTR/scripts/create_ref_panel/scripts/)
(`minimal_trim.py`, `remove_bundles.py`). Ported verbatim: `count_motif`,
`is_compound` (0.8 threshold), `minimal_trim` (bound-checks-first-then-match,
returns `None` like GangSTR's `-1,-1`), `is_close`/`remove_bundles` clustering,
the copy-number floors, and the copy count computed from the **original** TRF
span (integer division) as a post-trim accept-gate.

Two adaptations from the architecture (§3/§4), both commented at the call site:
1. **Sequence source.** GangSTR reads the repeat string from a TRF column; we
   slice it from the resident contig (`contig_seq[start..end]`), upper-cased.
2. **Motif + keep-imperfect.** The catalog motif is the first `period` bases of
   the reference tract (verbatim/phase-faithful, types §5), not TRF's `pattern`.
   We keep imperfect single-motif loci above a purity floor instead of GangSTR's
   perfect-only `remove_messy` (`ref == query`).

Also: per arch §4 we keep period-1 homopolymers (with the copy-floor of 10) —
GangSTR's separate `awk '$4 != 1'` homopolymer drop is **not** ported, matching
our arch §4 which lists only `period ≤ 6`. Flagged in the module doc.

## Tests (11, all green)

Unit: `is_compound` (ATAT/ATATAT vs AT/ATC/ATCG); `recompute_purity`
(perfect = 1.0, one interruption = 7/8); `minimal_trim` (clean trim, `None`
without two clean copies); `drop_bundles` (whole-cluster drop vs isolated
survivor, all-kept when none close). End-to-end `build_loci`: clean locus with
flanks; drops period > 6 + compound; drops below copy floor; clamps flanks at a
contig end.

## Deferred (next increments)

- `trf.rs` locate/spawn/`parse_bed_line` (trf-mod available now) + a golden BED
  test from real trf-mod output.
- `run()` orchestrator (per-contig fan-out/collect) + `ssr-catalog` CLI; the
  `write_index` (CSI) for the catalog reader's region query.
- Param defaults (`flank_bp`, `min_score`, `bundle_threshold ≥ flank_bp`) pinned
  when the orchestrator/CLI lands.

## Validation (dev container)

- `cargo fmt --check` → clean.
- `cargo clippy --all-targets --all-features -- -D warnings` → clean.
- `cargo test --lib` → all green (+11 over the prior count).
- `cargo doc --no-deps` → clean.
