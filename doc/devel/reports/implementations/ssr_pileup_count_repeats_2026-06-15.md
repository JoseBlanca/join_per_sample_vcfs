# SSR Stage-1 — `count_repeats` fast-path motif counter

**Date:** 2026-06-15 · **Skill:** rust-feature-implementation · **Branch:** `ssr-architecture`

A reordered increment: the container generalization (build-order task 3,
architecture §10) is a large multi-step refactor of the production `.psp`
writer/reader and gates only the writer/round-trip, so by agreement we paused it
and built a stage module that needs **only** the already-landed `types` +
`norm_seqs` — the fast-path counter
([plan §4](../../../doc/devel/implementation_plans/ssr_pileup.md)).

## 1. Plan

New module tree `src/ssr/pileup/` with `count_repeats.rs` holding the
dependency-free core of the fast path: confirm a read's tract is a pure integer
tiling of the locus motif and, if so, return the repeat count + a base-quality
weight. The triage-typed wrapper (taking a classified read + locus and anchoring
the tract) is deferred to `triage.rs`, so this core takes raw bytes and is
testable in isolation.

## 2. Assumptions (silent choices surfaced)

- **Weight semantics.** The plan/§11 says a fast-path read emits a base-quality
  `weight`, summed across reads into `hist_weight`, "a confidence aggregate Stage 2
  uses to down-weight a length whose confident support is all low-quality" — but
  the exact per-read formula is left to the spec §4.3 column contract. I chose the
  **mean per-base probability-correct** over the tract (`1 − 10^(−Q/10)` averaged),
  a bounded `[0,1]` confidence scalar. Documented on the fn; revisit if the Stage-2
  contract wants a different aggregate (e.g. min, or a soft-count sum).
- **Empty tract → `(0, 1.0)`.** A zero-unit tract (full repeat deletion) is treated
  as a vacuous clean tiling with full confidence, rather than `None`.
- **`u16` overflow → `None`.** Counts are `u16` end-to-end (shared-types §2); a
  tract long enough to overflow is rejected, not truncated (not a real
  spanning-read allele anyway).
- **Scope: the tiling check + weight only.** No read I/O, no triage types, no
  locus anchoring — those land with `triage.rs`. This keeps the increment
  dependency-free and the function unit-testable on byte slices.

## 3. Changes made

- **New** [src/ssr/pileup/mod.rs](../../../src/ssr/pileup/mod.rs): the Stage-1
  module root (doc + `pub mod count_repeats;`).
- **New** [src/ssr/pileup/count_repeats.rs](../../../src/ssr/pileup/count_repeats.rs):
  `count_pure_tiling(tract, quals, &Motif) -> Option<(u16, f32)>` + a process-wide
  `PHRED_CORRECT` 256-entry lookup (the BAQ `Q2P` *pattern*, not its code).
- [src/ssr/mod.rs](../../../src/ssr/mod.rs): `pub mod pileup;` (the existing
  module-level `#![allow(dead_code)]` covers the not-yet-wired core).

## 4. Tests added (9, in-module)

Pure tilings (di-/tri-nucleotide, homopolymer); rejection of an interrupted tract
and of a partial trailing unit; empty-tract vacuous case; weight tracks quality
(uniform Q10 → 0.9) and is the mean of mixed qualities; malformed
`tract`/`quals` length mismatch → `None`.

## 5. Validation results

Run in the dev container (`./scripts/dev.sh`):
- `cargo fmt -- --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean (fixed two
  findings during the run: `manual_is_multiple_of`; and `private_interfaces` — the
  `pub fn` was demoted to `pub(crate)` since it takes the `pub(crate)` `Motif`).
- `cargo test --lib ssr::` — **33 passed, 0 failed** (24 types + 9 count_repeats).

## 6. Tradeoffs and follow-ups

- The triage-typed `count_fast(read, hit, locus)` wrapper is deferred to
  `triage.rs` (needs `SpanningRead`/`FastHit`, not yet built).
- Weight formula is a documented assumption pending the Stage-2 §4.3 contract.
- **Container generalization (task 3) remains open** — to return to after more
  stage modules, or as its own focused pass.
- **Next candidate:** `pair_hmm.rs` (the slow-path banded forward) — also
  dependency-free (byte slices + quals), the stage's main numerical risk.
