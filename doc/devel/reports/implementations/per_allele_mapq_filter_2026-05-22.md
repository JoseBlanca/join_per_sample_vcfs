# Per-allele MAPQ filter — implementation report (2026-05-22)

Implements [`doc/devel/implementation_plans/per_allele_mapq_tracking.md`](../../implementation_plans/per_allele_mapq_tracking.md).
Follow-up on [the caller-correctness-against-GATK report (2026-05-20)](caller_correctness_against_gatk_2026-05-20.md)
which established the GATK-comparison baseline F1 = 0.285.

## What landed

Three sequential commits on `main`:

| Commit | Phase | What |
|---|---|---|
| `a87e9a6` | A | `AlleleSupportStats.mapq_sum: u32` + `mapq_sum_sq: u64`; threaded from BAM walker through `PreparedRead` → `ReadContribution` → per-group merger compound rebuilds; two new PSP columns (`0x15 AlleleMapqSum`, `0x16 AlleleMapqSumSq`); round-trip test `mapq_scalars_round_trip`. |
| `5903aa0` | B | `INFO/MQRef` (Number=1), `INFO/MQAlt`/`MQDiff`/`MQDiffT` (Number=A); cohort-pooled mean + Welch's t computed in `record_encode.rs`; integration test `mapq_info_fields_reflect_cohort_pooled_stats`. |
| `e8f3e90` | C | `--no-mapq-diff-filter` and `--min-mapq-diff-t` flags; hard-drop at variant emission; new `records_dropped_low_mapq_diff_t` counter; integration test `mapq_diff_t_filter_decision_matches_thresholds`. |

Default `--min-mapq-diff-t = -3.0` was picked from the
2026-05-22 cross-validation against GATK `MQRankSum` on the synthetic
10-duplicate cohort (catches 25 % of GATK-extreme sites at a 2.0 %
false-flag rate on GATK-clean sites — see the implementation plan
for the full table).

## Validation against GATK — apples-to-apples on the same PSPs

PSPs regenerated under the new format
([`tmp/cohort_synth_multichrom/*.psp`](../../../../tmp/cohort_synth_multichrom/);
BAMs were first converted to CRAMs at
[`tmp/cohort_synth_multichrom_crams/`](../../../../tmp/cohort_synth_multichrom_crams/)
because `pileup` accepts CRAM only). With the new PSPs **and the
filter disabled**, every count matches the pre-regen baseline to the
integer — round-tripping the call set through BAM → CRAM → new-PSP
introduced no drift.

| metric | filter OFF (baseline) | filter ON (default `-3.0`) | Δ |
|---|---:|---:|---:|
| raw records emitted | 77,781 | **58,331** | −25.0 % |
| ours total (S0000, in regions, post-norm) | 105,251 | 72,301 | −31.3 % |
| TP | 17,759 | 14,533 | −18.2 % |
| FP | 87,492 | 57,768 | **−34.0 %** |
| FN | 1,683 | 4,909 | +191.7 % |
| precision | 16.87 % | **20.10 %** | **+3.23 pp** |
| recall | 91.34 % | 74.75 % | −16.59 pp |
| **F1** | **0.285** | **0.317** | **+0.032 (+11 %)** |

Marginal-region FP enrichment: of the 32,950 records dropped between
filter-off and filter-on, **29,724 (90.2 %) were FPs and 3,226 (9.8 %)
were TPs**. The cohort baseline was 83.1 % FP / 16.9 % TP — the
filter discards regions ~7 pp more FP-enriched than the unfiltered
cohort, which is the discriminator behaviour the plan targeted.

## Drop-counter calibration

The new `records_dropped_low_mapq_diff_t={N}` counter on the
`var-calling: …` summary line reported **19,450** records dropped.
The plan's back-of-envelope estimate was ~3,650 + ~790 ≈ 4,400.
The 4× gap comes from:

1. The counter is **raw cohort records dropped before the writer**;
   the comparison numbers above are **post-norm split-multi-allelic
   entries** for S0000 only. One raw multi-allelic record may
   expand into several S0000 lines.
2. Multi-allelic records trip the drop if **any one ALT** fails the
   Welch's t threshold; the analysis was per-(REF, ALT) pair.
3. The analysis sampled 3,000 FPs / 3,000 TPs and extrapolated;
   the actual cohort has 87 K + 17 K post-norm S0000 entries plus
   the many records that don't reach S0000 in the regions cut.

The shape (filter targets FP-enriched regions; precision up, recall
down, F1 up) matches the prediction.

## Caveats — what this report does NOT claim

* The benchmark is still the GATK 10-sample raw VCF on a 10-duplicate
  synthetic cohort. See
  [`feedback_benchmark_validity.md`](../../../.claude/projects/-home-jose-devel-join-per-sample-vcfs/memory/feedback_benchmark_validity.md)
  for why F1 deltas in low-mappability regions are partially
  benchmark-circular (GATK emits multi-mapper calls too, just with
  bad INFO annotations).
* The synthetic 10-duplicate cohort makes paralog hotspots look like
  consensus 10/10-alt calls. A real independent-sample cohort would
  shift some of these numbers — the curated dataset that's being
  prepared is the authoritative test.
* DUST stays at default; the MAPQ-diff filter is additive on top.

## Artifacts

* VCFs:
  * [`tmp/cohort_mc_mapqfilter.vcf`](../../../../tmp/cohort_mc_mapqfilter.vcf) — with default filter
  * [`tmp/cohort_mc_dust_nomapq.vcf`](../../../../tmp/cohort_mc_dust_nomapq.vcf) — filter disabled (post-regen baseline)
* Comparisons: [`tmp/gatk_compare_mapqfilter/`](../../../../tmp/gatk_compare_mapqfilter/), [`tmp/gatk_compare_nomapq/`](../../../../tmp/gatk_compare_nomapq/)
* Driver: [`tmp/gatk_compare/run_compare_dust.sh`](../../../../tmp/gatk_compare/run_compare_dust.sh)

## What's next

* Re-evaluate on the curated independent-sample cohort when ready.
  In particular, watch whether the 18 % recall loss in the synthetic
  cohort is mostly multi-mapper-hotspot positions both we and GATK
  emit (acceptable — those are the calls the filter targets) or
  whether the filter is also clipping real-sequence heterozygotes
  (would need looser threshold).
* Optionally sweep the threshold (`-2.0` more aggressive, `-4.0`
  more conservative) on the curated cohort once it's available.
  Not worth doing on the synthetic data — the benchmark validity
  is already the binding constraint.
