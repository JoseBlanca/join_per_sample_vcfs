# SSR Stage-1 — `candidate_generation` Job 2 (off-ladder candidates)

**Date:** 2026-06-15 · **Skill:** rust-feature-implementation · **Branch:** `ssr-architecture`

Completes `candidate_generation` (arch §6): the read-derived off-ladder
candidate + its canonical key, after settling the off-ladder contract with the
user (A + B1, simplest form).

## 1. Contract decided (with the user)

- **A — triage hands in raw observed-tract bytes** (no alignment / indel
  position). The off-ladder generation gate (slow-path reason) is upstream;
  `build_offladder` takes just `(locus, observed_tract)`.
- **B1 — `NormalizedSeq` = the full canonical tract** (not a trimmed delta), so
  `Allele::to_sequence` means "the tract" uniformly for on- and off-ladder.

## 2. Key finding — the `norm_seqs` kernel is **not needed** for B1

Digging into B1, the canonicalization is **verbatim**: the observed tract bytes
*are* the canonical key, no left-alignment applied. This is correct, not a
shortcut:

- A literal alt **sequence** is invariant under left-alignment — left-alignment
  canonicalizes an indel's *position in a `(ref, alt)` description*, never the alt
  bytes. Storing the sequence (B1) stores a left-alignment-invariant object.
- Stage 0 guarantees **clean flanks** (arch §5.9), so the tract is
  deterministically delimited; two reads of the same molecule yield the same tract
  bytes → already byte-equal across samples.

So the cross-sample-identity invariant holds **by construction** under B1. The
shared `norm_seqs` kernel (lifted in task 2) would only be needed for a
*trimmed-variant* representation (the `(offset, ref, alt)` delta form shared-types
§4 rejected). **Consequence:** task 2's "SSR off-ladder = the kernel's second
consumer" does **not** materialize under B1 — the kernel stays validly shared and
used by the SNP path, ready if we ever adopt a trimmed form or relax clean flanks.
Surfaced here so the decision is visible, not buried.

## 3. Changes made

- [src/ssr/pileup/count_repeats.rs](../../../src/ssr/pileup/count_repeats.rs):
  extracted `pure_tiling_units(tract, motif) -> Option<u16>` (the tiling test
  without the quality weight); `count_pure_tiling` now calls it. Shared with the
  off-ladder degenerate-case check.
- [src/ssr/pileup/candidate_generation.rs](../../../src/ssr/pileup/candidate_generation.rs):
  `normalize_offladder(observed_tract) -> NormalizedSeq` (verbatim; named so the
  canonicalization contract has one home + change point) and
  `build_offladder(locus, observed_tract) -> Option<CandidateAllele>` (drops a
  pure tiling as a rung, else wraps `left_flank + tract + right_flank` with an
  `OffLadder` key). Module doc records the A + B1 contract and the invariance
  argument.

## 4. Tests added (5, in-module)

`normalize_offladder` verbatim; `build_offladder` wraps an impure tract with the
flanks + the right `OffLadder` key; a pure tiling (incl. empty) is rejected as a
rung; the `OffLadder` allele round-trips through `to_sequence`; two reads of the
same tract give an identical key (the cross-sample-identity property — and the
fact that there is only one spelling of a literal sequence *is* the proof
left-alignment is unnecessary).

## 5. Validation results

Dev container (`./scripts/dev.sh`):
- `cargo fmt -- --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib ssr::pileup::candidate_generation` — **11 passed** (6 + 5).
- `cargo test --lib ssr::pileup::count_repeats` — **9 passed** (refactor green).

## 6. Tradeoffs and follow-ups

- B1-verbatim is the simplest correct form; if a future need (trimmed-variant
  storage, or relaxed flanks) arises, `normalize_offladder` is the single change
  point to route through `norm_seqs`.
- `candidate_generation` is now complete (Jobs 1 + 2). Remaining Stage-1 work
  needs the read-input seam (`triage`/`fetch_reads`) or the container schema
  (`locus_record`).
