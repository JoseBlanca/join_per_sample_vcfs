# SSR Stage-1 — `candidate_generation` (Job 1: on-ladder rungs)

**Date:** 2026-06-15 · **Skill:** rust-feature-implementation · **Branch:** `ssr-architecture`

The translation layer between read length-estimates and typed alleles
([plan §6](../../../doc/devel/implementation_plans/ssr_pileup.md), arch §6). This
increment builds **Job 1** — the on-ladder rung candidates — which composes
task-1's `Allele::to_sequence` with the catalog `Locus`. **Job 2** (the
read-derived off-ladder candidate + its canonical normalization onto
`norm_seqs`) is deferred pending a design decision (below).

## 1. Plan

- `CandidateAllele { candidate_seq: Vec<u8>, allele: Allele }` — a haplotype the
  forward scores + the allele key the result is recorded under.
- `build_rungs(locus, observed_count, w, out)` — append one rung per `L ∈
  [observed_count − w, observed_count + w]`, each `left_flank + (motif × L) +
  right_flank` read off the catalog (no FASTA).
- `STUTTER_WINDOW_UNITS` const (window half-width; calibration placeholder).

## 2. Assumptions (surfaced)

- **Low end clamps at 0.** A read near the bottom of the ladder yields fewer than
  `2w + 1` rungs (a repeat count can't be negative) rather than wrapping.
- **`out` is appended, not cleared** — it's a reused per-worker scratch; the
  caller owns clearing between reads.
- **`STUTTER_WINDOW_UNITS = 3`** is a placeholder calibration value (arch §14);
  callers pass `w` explicitly, so the const is only the default.

## 3. Deferred — Job 2 (off-ladder normalization), and why

The off-ladder adapter `normalize_offladder(observed_tract, ref_tract) ->
NormalizedSeq` drives `norm_seqs::normalize_alleles`, which needs the variant
**bounds** (the indel span) — not derivable from raw tract bytes alone. Two
contracts are unpinned and shouldn't be guessed (they guard the cross-sample
identity invariant):
1. **What triage hands in** as the observed tract — raw bytes, or an alignment
   carrying the indel position?
2. **What `NormalizedSeq` canonicalizes to** operationally — the full tract with
   any indel left-aligned, or a minimal trimmed variant representation?

Flagged to the user as a checkpoint rather than resolved silently.

## 4. Changes made

- **New** [src/ssr/pileup/candidate_generation.rs](../../../src/ssr/pileup/candidate_generation.rs):
  `CandidateAllele`, `build_rungs`, `STUTTER_WINDOW_UNITS`.
- [src/ssr/pileup/mod.rs](../../../src/ssr/pileup/mod.rs): `pub mod candidate_generation;`.

## 5. Tests added (6, in-module)

Full window (`observed_count 5, w 2` → counts `[3..7]`); a rung is
`flank+tiling+flank` (`GGG+CACA+TTT`); the **reference rung reproduces
`locus.ref_bytes()`** (sample-independence); low-end clamps to `0..=4` not
wrapping; the zero-unit rung is just the flanks; append accumulates across calls.

## 6. Validation results

Dev container (`./scripts/dev.sh`):
- `cargo fmt -- --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib ssr::pileup::candidate_generation` — **6 passed, 0 failed**.

## 7. Tradeoffs and follow-ups

- **Job 2 (off-ladder normalization) deferred** pending the two contract
  decisions above — the first real SSR consumer of the `norm_seqs` kernel.
- The `count ± W` rungs now feed `pair_hmm::forward`; the `score_candidates`
  wrapper (read × candidate set → `Qᵣ` distribution) can be built next, joining
  `build_rungs` + `forward`.
- Container generalization (task 3) still open.
