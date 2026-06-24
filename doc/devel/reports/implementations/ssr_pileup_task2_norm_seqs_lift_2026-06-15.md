# SSR Stage-1 task 2 — lift `normalize_alleles` to `src/norm_seqs/`

**Date:** 2026-06-15 · **Skill:** rust-feature-implementation · **Branch:** `ssr-architecture`

Second prerequisite of the Stage-1 build order
([plan §0/§2](../../../doc/devel/implementation_plans/ssr_pileup.md)): lift the
SNP caller's indel-normalization shift core into a shared, representation-neutral
module so the SSR off-ladder path (task 4) can reuse the *actual code*, not a
parallel copy (shared-types §4). A behaviour-preserving refactor with the SNP
tests as the regression gate.

## 1. Plan

- New flat module [src/norm_seqs.rs](../../../src/norm_seqs.rs) holding the
  CIGAR-agnostic kernel: `normalize_alleles(seqs, bounds, max_shift, trim)`, its
  parameter type `IndexRange`, and the three private `*_is_same` helpers.
- [src/pileup/walker/indel_norm.rs](../../../src/pileup/walker/indel_norm.rs)
  imports `{IndexRange, normalize_alleles}` and keeps everything else (the
  right-to-left CIGAR walk `left_align_cigar`, `build_cigar`, the wrappers).
- Regression gate: the existing SNP `indel_norm` test suite passes unchanged.

## 2. Assumptions / deviations from the plan

- **The kernel was less standalone than the plan assumed.** Plan §2 (echoing
  shared-types §4) described `normalize_alleles` as an "abstract `(seqs, bounds)`"
  fn that just needed to be made public. In fact it is a *shift calculator*
  returning `(start_shift, end_shift)` and mutating index ranges; its signature
  depends on the private `Range` type, and three helpers are its internals. So the
  lift had to carry `Range` + the helpers with it (the signature forces it) — still
  behaviour-preserving, just a slightly wider move than "make one fn `pub`."
- **Renamed `Range` → `IndexRange`** in the shared module. As a now-shared
  `pub(crate)` type, `Range` shadowed `std::ops::Range`; `IndexRange` is the name
  the original module doc already used ("Mirrors GATK's IndexRange") and matches
  the GATK source. Touch was 2 construction sites in `indel_norm.rs`.
- **Building the SSR adapter is *not* in this task.** Task 2 is the lift only; the
  `(off-ladder tract, ref tract) → NormalizedSeq` adapter lands in
  `candidate_generation.rs` (task 4), where it has a consumer.

## 3. Changes made

- **New** [src/norm_seqs.rs](../../../src/norm_seqs.rs): module doc (states the two
  users + the one-kernel rationale), `IndexRange` (`pub(crate)` struct + fields;
  externally-used methods `size`/`shift_left`/`shift_start_left`/`shift_end_left`
  `pub(crate)`, `shift`/`shift_start` private), the three private helpers, and
  `pub(crate) fn normalize_alleles`. Logic copied verbatim — no behavioural change.
- [src/lib.rs](../../../src/lib.rs): `pub(crate) mod norm_seqs;`.
- [src/pileup/walker/indel_norm.rs](../../../src/pileup/walker/indel_norm.rs):
  removed the moved items; added `use crate::norm_seqs::{IndexRange, normalize_alleles};`;
  renamed the 2 `Range { .. }` constructions to `IndexRange`; updated the module
  doc to say the shift core now lives in `norm_seqs` and this module is the
  CIGAR-walk wrapper.

## 4. Tests added (3 direct kernel tests in `norm_seqs`)

The 22 existing `indel_norm` tests already exercise `normalize_alleles` through
`left_align_cigar` (homopolymer / SSR / insertion / deletion / trim / convergence)
— they remain the primary gate. Added direct coverage of the now-public kernel:
- `trims_shared_prefix_and_suffix_to_the_minimal_variant` — parsimony trim of
  `TCAG`/`TGAG` down to the single-base variant at `[1,2)`, shifts `(-1, 2)`.
- `shifts_a_span_left_through_a_homopolymer_run` — a span slides left through an
  A-run, stopping at `max_shift`.
- `no_trim_no_shift_leaves_bounds_untouched` — the `trim=false, max_shift=0` no-op.

## 5. Validation results

Run in the dev container (`./scripts/dev.sh`):
- `cargo fmt -- --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **1065 passed, 0 failed, 1 ignored** (the SNP `indel_norm`
  regression gate green; +3 new `norm_seqs` kernel tests).

## 6. Tradeoffs and follow-ups

- The lift widened to carry `IndexRange` + helpers (signature-forced); documented
  above so it isn't mistaken for scope creep.
- **Next (task 3):** the container SSR schema (`registry_ssr` + `SsrLocusRecord`),
  if the container refactor hasn't landed — independent of the SSR math.
- **Task 4** wires the SSR off-ladder adapter onto this kernel and is the first
  real consumer beyond the SNP path.
