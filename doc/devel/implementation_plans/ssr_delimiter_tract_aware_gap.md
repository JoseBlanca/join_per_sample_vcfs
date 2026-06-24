# ssr-pileup delimiter — tract-aware gap penalty

**Status:** planned (branch `ssr-cohort`)
**Owner area:** Stage 1 (`ssr-pileup`) — `delimit_read` pair-HMM
**Investigation:** [ssr_delimiter_gap_penalty_2026-06-24.md](../reports/research/ssr_delimiter_gap_penalty_2026-06-24.md)
**Relates to:** the committed window-recovery infra (`feaa9ef`), which this unblocks.

## 1. Problem (recap)

`delimit_read` ([alignment.rs](../../../src/ssr/pileup/alignment.rs)) extracts a read's repeat
tract by a Viterbi pair-HMM alignment to the locus reference haplotype, using **one uniform
affine gap** (`GAP_OPEN_PROB = 2.9e-5`). That value is HipSTR's *flank* gap-open
(`dindel_probs[0]`); HipSTR scores tract length changes with a separate stutter model, not the
affine gap. Applied uniformly, the stiff gap collapses any allele ≥ ref+2 units to the
reference length (even with full flanks), because a multi-unit insertion loses to absorbing
the extra units as flank mismatches.

## 2. Goal

Make the gap penalty **tract-aware**, matching HipSTR's flank/tract separation: a **cheap**
gap inside the reference repeat tract (so a read's length difference is absorbed there) and
the **stiff** Dindel gap in the flanks (so they keep anchoring the junctions). The delimiter's
job is *extraction* (the empirical-candidate model — observed bytes, impurities included), so
the cheap tract gap is **per-base and content-agnostic**: we do not distinguish in-frame
(unit) from out-of-frame (non-unit) indels here. That distinction is a *scoring* concern and
belongs to Stage-2 (`ssr-call`)'s stutter model — a separate follow-up.

Non-goals: a whole-unit-jump / stutter-PMF aligner in the delimiter (HipSTR's
`StutterAlignerClass`); changing the empirical-candidate architecture; Stage-2 scoring
changes. Real-data calibration of the exact tract gap value (a starting value here; tune
later).

## 3. Design

### 3.1 Two gap regimes by reference column

The reference haplotype `hap = ref_bytes` is `left_flank (left_len)` + `ref_tract` +
`right_flank (right_len)`. Reference column `j` (1-based, `hap[j-1]`) is a **tract column**
iff `left_len ≤ j-1 < n - right_len`. In the Viterbi, the gap cost at a transition is selected
by the column the gap touches:

- **Insertion** (read base consumed, no hap base) entering/extending at column `j` → use the
  tract gap iff `j` is a tract column.
- **Deletion** (hap base `j-1` consumed) → use the tract gap iff `j` is a tract column.
- **Match** transitions are unchanged.

A read carrying a long allele then absorbs the extra bases as cheap tract insertions; a short
allele as cheap tract deletions; the flanks still pay the stiff Dindel gap and anchor.

### 3.2 The gap values

`HmmModel` gains a second gap regime (tract): `ln_gap_open_tract` (cheap) and
`ln_gap_extend_tract`. The flank regime keeps the current `GAP_OPEN_PROB = 2.9e-5` /
`GAP_EXTEND ≈ e⁻¹`.

- **`GAP_OPEN_PROB_TRACT`** — a new named constant, the per-base tract slippage rate. Starting
  value chosen so the collapse cliff is pushed past the plausible allele range (the sweep:
  `1e-2` recovers ref+5; a cheaper value recovers more). Documented as a **provisional
  calibration constant** (like the other `dev_default`s), to be tuned against real data /
  reconciled with the Stage-2 stutter rate.
- **`GAP_EXTEND`** — start by reusing the current extend for the tract too (the dominant fixed
  cost is the open); flag `ln_gap_extend_tract` as a knob if long alleles still need it
  cheaper.

### 3.3 Determinism / hot path

The per-column tract test is a pure function of the (fixed) locus geometry; the traceback
tie-break is unchanged. So the byte-identity-across-thread-count contract holds. The DP gains
one branch (or a precomputed per-column gap-open lookup) per cell — negligible; the inner loop
shape is unchanged.

## 4. Touch points

1. **`src/ssr/pileup/alignment.rs`**
   - `HmmModel`: add `ln_gap_open_tract` (+ `ln_gap_extend_tract`); new const
     `GAP_OPEN_PROB_TRACT`.
   - `delimit_read`: thread the tract column bounds (`left_len`, `n - right_len`) into the DP;
     select the gap regime per column at the insertion / deletion (and row-0 deletion / col-0
     insertion) transitions. A small `gap_open_at(j)` / `gap_extend_at(j)` helper keeps the
     inner loop readable.
2. **Tests (`alignment.rs`)**
   - `delimit_recovers_long_alleles_within_the_tract`: clean full-flank `CA×k` for
     `k = ref-4 … ref+N` all extract `k` (the cliff is gone across the plausible range) — the
     direct regression for §1.
   - `delimit_keeps_flanks_anchored`: a read with a flank **substitution / short flank indel**
     still delimits the tract correctly (the stiff flank gap still anchors; the cheap tract gap
     did not leak into the flank).
   - `delimit_extracts_an_impure_tract_verbatim`: an interrupted tract (`CACAACA`) is still
     extracted verbatim (cheap tract gap handles out-of-frame content too — the
     empirical-candidate contract).
   - Determinism: extend an existing thread-invariance check, or assert `delimit_read` is
     identical at the two ends of a `CA×k` sweep regardless of run.
3. **Re-enable the deferred window-recovery behavioural tests** (`driver.rs`): now that the
   delimiter recovers long alleles, `classify_read` genuinely recovers a window-truncated long
   allele (`WidenedSequence`) and drops an unrecoverable one (`WindowTruncated`). Restore the
   two tests parked in `feaa9ef` and the e2e `CA×10` recovery guard.
4. **Update** the window-recovery plan's gap-open note + PROJECT_STATUS; impl report.

## 5. Validation

- The delimiter sweep test (§4.2) is the acceptance signal: no collapse across `ref ± N` for a
  realistic `N`.
- Re-run the e2e BAM→VCF test with a `CA×10` cohort → the long allele now reaches the VCF.
- fmt / clippy `-D warnings` / full `cargo test` green; byte-identity across thread counts
  preserved.

## 6. Order of work

1. `HmmModel` tract gap fields + `GAP_OPEN_PROB_TRACT` const.
2. `delimit_read` per-column gap selection + the `gap_*_at(j)` helper.
3. Delimiter sweep / flank-anchor / impurity tests (red→green).
4. Re-enable the window-recovery behavioural + e2e tests.
5. Gate + docs + impl report.

## 7. Follow-up (separate task, your idea)

Stage-2 (`ssr-call`) scoring: confirm the cohort EM's stutter model scores **in-frame
(unit) vs out-of-frame (non-unit)** length differences distinctly (HipSTR's
`log_stutter_pmf`), and reconcile the delimiter's provisional `GAP_OPEN_PROB_TRACT` with that
slippage rate. The delimiter extracts; Stage-2 scores — this item lives on the scoring side.
