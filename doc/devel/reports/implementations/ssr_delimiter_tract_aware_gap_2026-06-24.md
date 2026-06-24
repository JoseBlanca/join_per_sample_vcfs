# Implementation: ssr-pileup delimiter — tract-aware gap penalty

**Date:** 2026-06-24
**Branch:** `ssr-cohort`
**Plan:** [ssr_delimiter_tract_aware_gap.md](../../implementation_plans/ssr_delimiter_tract_aware_gap.md)
**Investigation:** [ssr_delimiter_gap_penalty_2026-06-24.md](../research/ssr_delimiter_gap_penalty_2026-06-24.md)

## 1. Plan

Make the Stage-1 delimiter's pair-HMM gap penalty **tract-aware**: a stiff Dindel gap in the
flanks (anchors the junctions) and a soft gap inside the repeat tract (so a read's length
difference is absorbed there instead of collapsing the tract to the reference). This matches
HipSTR's flank/tract split. The delimiter only *extracts* the observed tract (the
empirical-candidate contract), so the tract gap is a flat, content-agnostic per-base value;
the in-frame-vs-out-of-frame unit distinction stays a Stage-2 scoring concern.

## 2. Domain intent

`delimit_read` extracts each spanning read's observed repeat tract by a Viterbi alignment to
the locus reference haplotype. **Contract:** the extracted bytes must be the read's actual
tract, for any allele length the read can span — verbatim, including impurities. The previous
uniform gap (`2.9e-5`, HipSTR's *flank* value) collapsed any allele ≥ ref+2 units to the
reference length, so longer alleles never reached Stage 2. The tract-aware gap restores the
contract while keeping the flank anchoring (and thus the byte-identity-across-thread-count
guarantee) intact.

## 3. Assumptions / decisions

- **Per-base, content-agnostic tract gap** (not unit-aware). For *extraction*, a flat cheap
  tract gap is enough to avoid the collapse and correctly reads out impure / out-of-frame
  tracts verbatim. Distinguishing in-frame (unit) from out-of-frame (non-unit) stutter is a
  *scoring* property and belongs to Stage-2 (`ssr-call`) — recorded as the follow-up.
- **`GAP_OPEN_PROB_TRACT = 1e-2`** — a provisional calibration constant. Chosen so the
  collapse cliff is gone across the realistic allele range (verified `ref−4 … ref+12` with no
  collapse); to be reconciled with the Stage-2 stutter slippage rate on real data.
- **Tract column predicate `left_len < j ≤ n − right_len`** keys the gap on the reference
  column a gap touches (insertion beside, or deletion of, a tract base). Junction-exact
  placement of the cheap insertions does not change the *extracted span* (the junctions are
  fixed by the flank-anchoring Match transitions), so the predicate boundary is not
  load-bearing.
- **Gap-extend / -close left uniform.** The gap-open is the dominant fixed cost; the extend is
  a documented tuning knob if very long alleles later need it cheaper.

## 4. Changes made

- [`src/ssr/pileup/alignment.rs`](../../../../src/ssr/pileup/alignment.rs):
  - New `const GAP_OPEN_PROB_TRACT: f64 = 1e-2`; `GAP_OPEN_PROB` doc clarified as the *flank*
    value.
  - `HmmModel` gains `ln_gap_open_tract`; `Default` sets it from the new const.
  - `delimit_read`: a `gap_open(j)` closure selects the tract vs flank gap-open per reference
    column; applied at the three gap-open transition sites (row-0 leading deletion, per-cell
    insertion, per-cell deletion). Col-0 insertion (leading read junk) stays flank-priced.
    Match / extend / close transitions unchanged.
- [`src/ssr/pileup/driver.rs`](../../../../src/ssr/pileup/driver.rs): re-enabled the two
  long-allele behavioural tests parked in `feaa9ef` (now that the delimiter recovers long
  alleles), adapting the drop test to the actual safe-drop reasons.
- [`src/ssr/end_to_end_tests.rs`](../../../../src/ssr/end_to_end_tests.rs): a new full-chain
  `CA×10` recovery test (the original collapse symptom).

## 5. Tests added / updated

- `alignment::tests::tract_aware_gap_recovers_long_alleles_without_collapse` — the regression:
  clean full-flank `CA×k` for `k = 4..=20` all extract `k` units (no collapse, both shorter
  and longer than the `CA×8` reference).
- `alignment::tests::tract_aware_gap_keeps_the_flanks_anchored` — a right-flank substitution is
  tolerated as a base mismatch (stiff flank gap) and the tract still reads out as the clean
  `CA×9`; the soft tract gap does not leak into the flank.
- `driver::tests::classify_read_recovers_a_window_truncated_long_allele` — an all-Match `CA×10`
  read (window-truncated) now recovers the full `CA×10` via the widened re-delimit
  (`WidenedSequence`).
- `driver::tests::classify_read_drops_an_allele_too_long_to_span_the_window` — a `CA×20`
  all-Match read is dropped (`BorderOffEnd`/`WindowTruncated`), never silently collapsed.
- `end_to_end_tests::bam_to_vcf_recovers_a_long_allele_the_uniform_gap_collapsed` — BAM→VCF:
  3× hom-`CA×8` + 3× hom-`CA×10` → a PASS variant carrying the `CA×10` ALT with the long
  samples called REPCN 10/10. The original symptom, fixed through the whole chain — and a
  bonus confirmation that Stage-2 genotypes the recovered long allele correctly.
- The pre-existing thread-invariance / byte-identity tests still pass: the per-column tract
  test is pure locus geometry, so determinism is preserved.

## 6. Validation

- `cargo fmt --check` → 0
- `cargo clippy --all-targets --all-features -- -D warnings` → 0
- `cargo test --lib --all-features` → 0, **1305 passed**; 2 ignored
- `cargo test --all-targets --all-features` → only the pre-existing `psp_writer_perf` bench
  panic (baseline; unrelated). All lib + integration tests pass.
- Performance: the DP gains one branch (or const-folded comparison) per cell — negligible
  against the existing O(m·n) Viterbi; inner-loop shape unchanged. No bench covers the SSR
  path, so no criterion comparison.

## 7. Tradeoffs and follow-ups

- **`GAP_OPEN_PROB_TRACT` is provisional** — a real-data calibration constant, to be reconciled
  with the Stage-2 stutter rate.
- **Follow-up (separate task):** confirm/strengthen Stage-2 (`ssr-call`) stutter scoring of
  **in-frame (unit) vs out-of-frame (non-unit)** length differences (HipSTR's `log_stutter_pmf`
  split), and reconcile the tract gap rate with that slippage rate. The delimiter extracts;
  Stage-2 scores.
- **`WindowTruncated` is a narrow corner** in `classify_read`: with the delimiter fixed,
  widening recovers whenever any flank is present, and a tract that fills the window yields
  `BorderOffEnd` on the first pass — so `WindowTruncated` is a defensive catch, rarely hit.
  Kept (it is still a correct, non-tallied drop) and documented.
