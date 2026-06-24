# SSR Stage-1 task 1 — allele representation in `types.rs`

**Date:** 2026-06-15 · **Skill:** rust-feature-implementation · **Branch:** `ssr-architecture`

First prerequisite of the Stage-1 `ssr-pileup` build order
([plan §0/§1](../../../doc/devel/implementation_plans/ssr_pileup.md)): the
designed-but-unbuilt allele model from
[shared-types §2–§5](../../../doc/devel/architecture/ssr_shared_types.md). No new
files — this grows [src/ssr/types.rs](../../../src/ssr/types.rs).

## 1. Plan

Add the two types and their pure reconstruction methods:
- `NormalizedSeq(Box<[u8]>)` — canonical off-ladder tract; owns the
  *bytes-are-canonical* invariant only (does **not** normalize — that is task 2's
  lifted `norm_seqs` kernel). `Eq`/`Hash`/`Debug` treat the bytes as opaque.
- `Allele { OnLadder { units: u16 }, OffLadder(NormalizedSeq) }` — two encodings
  of one tract sequence; `==` is sequence identity.
- `Allele::to_sequence(&Locus) -> Vec<u8>` — the tract bytes (no flanks).
- `Allele::repeat_count(&Locus) -> f64` — integer for on-ladder, fractional for
  off-ladder.

## 2. Assumptions (silent choices surfaced)

- **`OnLadder::to_sequence` is a clean tiling (`motif × units`).** Shared-types §2/§5
  mentions on-ladder sequence "carrying any fixed interruption structure of an
  imperfect locus." I implemented the clean tiling instead, with imperfect-locus
  interruptions represented via `OffLadder` candidates. Rationale: it matches the
  Stage-1 rung definition in [arch §6](../../../doc/devel/architecture/ssr_pileup.md)
  (`left_flank + motif×L + right_flank`, a clean tiling) and avoids the
  underspecified per-locus arm-decomposition. If interruption-aware rungs are later
  required, they belong in `candidate_generation.rs`, not in this pure method.
- **`to_sequence` returns the tract only, not flanks.** Shared-types §5 calls it
  "the tract sequence (the VCF REF/ALT bytes)"; the candidate-haplotype builder
  (`candidate_generation.rs`, task 4) wraps the locus flanks around it.
- **No canonicality validation in `NormalizedSeq::new`.** types.rs cannot verify
  canonical form (that is the normalizer's contract); the constructor trusts the
  caller, as the design specifies ("`types.rs` only owns the invariant").

## 3. Changes made

[src/ssr/types.rs](../../../src/ssr/types.rs):
- `NormalizedSeq(Box<[u8]>)` + `new(impl Into<Box<[u8]>>)`, `as_bytes()`, manual
  `Debug` (ASCII-rendering, mirroring `Motif`).
- `Allele` enum (`#[derive(Clone, Debug, PartialEq, Eq, Hash)]`).
- `Allele::to_sequence` and `Allele::repeat_count` (pure, read `locus.motif()` /
  `locus.period()`; never the FASTA).

The module-level `#![allow(dead_code)]` in `src/ssr/mod.rs` is unchanged — these
types stay dead until `candidate_generation.rs` consumes them (tracked Mi5
follow-up).

## 4. Tests added (8, all in the module)

- `on_ladder_to_sequence_tiles_the_motif` — 3/1/0 units, incl. empty tract.
- `on_ladder_to_sequence_round_trips_to_units` — `len/period == units` for a sweep.
- `on_ladder_repeat_count_is_the_integer_units`.
- `off_ladder_to_sequence_returns_canonical_bytes_verbatim`.
- `off_ladder_repeat_count_is_fractional` — 7 bp / period 2 = 3.5.
- `allele_equality_is_sequence_identity` — same/diff rung, same/diff bytes, and
  the cross-variant non-collision (on-ladder vs off-ladder with equal bytes).
- `allele_dedups_in_a_hashset` — the cohort-union primitive.
- `normalized_seq_exposes_its_bytes_and_renders_ascii`.

## 5. Validation results

Run in the dev container (`./scripts/dev.sh`):
- `cargo fmt -- --check` — clean (no diff).
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ssr::types` — **24 passed, 0 failed** (16 existing + 8 new);
  1039 other lib tests filtered out, none touched.

## 6. Tradeoffs and follow-ups

- Imperfect-locus interruption modeling deferred to `candidate_generation.rs` (see
  Assumptions); the off-ladder path is the representation for non-rung sequences.
- Types remain `#[allow(dead_code)]` until task 4 wires consumers.
- **Next (task 2):** lift `normalize_alleles` into the shared `src/norm_seqs/`
  module (SNP e2e tests as the regression gate), then it backs `NormalizedSeq`
  construction on the SSR off-ladder path.
