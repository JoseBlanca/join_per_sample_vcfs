# Per-allele mean MAPQ — tracking + INFO emission + LowMQDiffT drop

Proposal date: 2026-05-22.

## Domain intent

The current pipeline folds per-read MAPQ into the allele's `q_sum` as
`max(ln(P_err_BQ_BAQ), ln(P_err_MQ))` at the read level
([src/per_sample_pileup/pileup/open_record.rs:766](../../src/per_sample_pileup/pileup/open_record.rs#L766)).
This loses per-allele MAPQ separately — by the time the merger sees
the allele, we no longer know whether the alt-supporting reads had
systematically lower MAPQ than the ref-supporting reads.

That asymmetry is the classic multi-mapper fingerprint: a paralog
alignment carrying a divergent base reaches a position with a lower
MAPQ than the locally-unique ref-supporting reads.
[Yesterday's drill-down](../reports/implementations/caller_correctness_against_gatk_2026-05-20.md)
found this pattern at the `SL4.0ch00:556000-557000` multi-mapper
hotspot; today's MAPQ-diff analysis
([[project-mapq-diff-signal-validated-2026-05-22]]) confirmed the
signal exists in pop_var_caller's output and correlates with GATK
`MQRankSum` at Pearson `r = +0.48` (Welch's t variant).

This plan adds:

1. Per-allele per-sample **`mapq_sum: u32` + `mapq_sum_sq: u64`** in
   `AlleleSupportStats`, propagated through PSP I/O and the
   pipeline. These two summaries plus the existing `num_obs` are
   sufficient to compute mean MAPQ and sample variance for each
   allele × sample × site without storing per-read data.
2. Cohort-level **INFO/MQRef**, **INFO/MQAlt**, **INFO/MQDiff**,
   **INFO/MQDiffT** on each VCF record, the last being a
   Welch's-t-test statistic comparing alt-supporting vs
   ref-supporting MAPQ distributions across the cohort.
3. A **hard-drop filter** triggered by a configurable
   `--min-mapq-diff-t` threshold (default −3.0). Filter applies to
   the variant emission stage alongside the existing
   `records_dropped_low_qual` / `records_dropped_low_alt_obs`
   counters.

## Why now

* The signal is validated against an independent reference (GATK
  `MQRankSum`); no more measurement work is needed before
  implementation.
* Per-allele storage extension is small — two scalars per allele ×
  sample. All downstream consumers of `AlleleSupportStats` get the
  new fields automatically; the work is mostly plumbing.
* Hard-drop default matches the existing
  `--min-alt-obs-per-sample` and `--min-qual` pattern: naive user
  gets a clean VCF without learning about FILTER tags or
  `bcftools view -f PASS`.

## Non-goals

* **Exact Mann-Whitney rank-sum / true GATK `MQRankSum` parity.**
  Would require storing a per-allele MAPQ histogram (~28 bytes/
  allele×sample vs 12 for Welch's). Welch's t already gives a
  12.5× extreme-vs-clean ratio at the recommended default; the
  histogram path stays available if a future real-cohort
  evaluation finds Welch's gives too many false flags on bimodal
  MAPQ distributions.
* **Per-sample FORMAT fields.** Cohort-level INFO matches GATK
  convention and is what filtering decisions consume. Per-sample
  MAPQ would inflate the VCF with little extra signal.
* **`FILTER=LowMQDiffT` tag.** Per [[feedback-no-silent-intermediates]]
  and the user's "make life easy for a naive user" call:
  pathological records are dropped, not labeled. Power users
  disable the drop with `--no-mapq-diff-filter`.

---

## Storage — `AlleleSupportStats`

[src/per_sample_pileup/pileup/mod.rs:346–361](../../src/per_sample_pileup/pileup/mod.rs#L346)

Add two fields:

```rust
pub struct AlleleSupportStats {
    pub num_obs: u32,
    pub q_sum: f64,
    pub fwd: u32,
    pub placed_left: u32,
    pub placed_start: u32,
    pub mapq_sum: u32,      // NEW: Σ mapq over reads supporting this allele
    pub mapq_sum_sq: u64,   // NEW: Σ mapq² over reads supporting this allele
}
```

**Type choice rationale.** MAPQ is `u8` (0..255, capped at 60 in
practice by BWA-MEM). For one sample at one site:

* `mapq_sum`: max plausible value 60 × `num_obs`. With
  `num_obs: u32` (max ~4 G), `mapq_sum` headroom is ~71 M reads
  at MAPQ 60 — comfortably u32.
* `mapq_sum_sq`: max plausible value 60² × `num_obs` = 3600 ×
  `num_obs`. At extreme depths (~10⁶) the sum can exceed u32 max
  (~4 G); u64 is the safe choice.

**Mean + variance from these:**

```text
mean = mapq_sum / num_obs
var  = (mapq_sum_sq − num_obs · mean²) / (num_obs − 1)
     = (mapq_sum_sq − mapq_sum² / num_obs) / (num_obs − 1)
```

**Cohort pooling for INFO.** Sums and counts are additive across
samples, so the cohort-level Welch's t is computable from cohort-
pooled sums without revisiting per-sample data:

```text
n_a    = Σ_samples num_obs[a][s]      for allele a
sum_a  = Σ_samples mapq_sum[a][s]
ssq_a  = Σ_samples mapq_sum_sq[a][s]
mean_a = sum_a / n_a
var_a  = (ssq_a − sum_a² / n_a) / (n_a − 1)

t = (mean_alt − mean_ref) / sqrt(var_alt / n_alt + var_ref / n_ref)
```

`t` is dimensionless; negative `t` means alt MAPQs are lower than
ref MAPQs (the suspect direction).

---

## Read aggregation — where MAPQ enters

Single update site:
[src/per_sample_pileup/pileup/open_record.rs:766–773 + 840–846](../../src/per_sample_pileup/pileup/open_record.rs#L766).

The current code already has `read.mapq` available
([baq/engine.rs:401](../../src/per_sample_pileup/baq/engine.rs#L401))
and folds it into `q_sum`. Pass the same `mapq: u8` into
`add_contribution()` and add:

```rust
self.mapq_sum    += mapq as u32;
self.mapq_sum_sq += (mapq as u64) * (mapq as u64);
```

No new BAM access; no new tracking outside `AlleleSupportStats`.

---

## PSP format — extend V1_0_COLUMNS in place (no version bump)

[src/per_sample_pileup/psp/registry.rs:309–446](../../src/per_sample_pileup/psp/registry.rs#L309).

Add two column tags directly to `V1_0_COLUMNS`:

```text
0x15  AlleleMapqSum     (u32 per allele)
0x16  AlleleMapqSumSq   (u64 per allele)
```

**No version bump.** Pop_var_caller is pre-alpha and no real
analyses depend on PSP format stability yet — the existing tmp/
fixtures will need regeneration after this change, but a version
bump would just add ceremony around the same regen step.

**Fixture regeneration:** in-tree synthetic PSPs in
`src/per_sample_pileup/psp/test_fixtures.rs` are generated on the
fly. Out-of-tree PSPs at `tmp/cohort_synth_multichrom/*.psp` and
`tmp/SRR7279725.small.psp` must be regenerated as part of Phase D.

---

## VCF INFO fields

[src/var_calling/vcf_writer/record_encode.rs:265+](../../src/var_calling/vcf_writer/record_encode.rs#L265),
header at
[src/var_calling/vcf_writer/header.rs](../../src/var_calling/vcf_writer/header.rs).

Four new INFO fields. Always emitted on records that survive the
drop filter; missing (`.`) when the relevant `n` is zero or the
test is undefined.

```text
##INFO=<ID=MQRef,   Number=1, Type=Float,
  Description="Mean mapping quality of reads supporting REF (cohort-pooled)">
##INFO=<ID=MQAlt,   Number=A, Type=Float,
  Description="Mean mapping quality of reads supporting each ALT (cohort-pooled)">
##INFO=<ID=MQDiff,  Number=A, Type=Float,
  Description="MQAlt - MQRef per ALT. Negative => alt reads lower MAPQ.">
##INFO=<ID=MQDiffT, Number=A, Type=Float,
  Description="Welch's t-statistic comparing ALT vs REF MAPQ distributions per ALT. Negative => alt reads lower MAPQ. Computed across the cohort.">
```

Emit `.` for any field whose denominator is zero (no reads of the
relevant class anywhere in the cohort) or whose Welch denominator
collapses to zero variance with nonzero mean shift.

---

## Drop filter

### Mechanics

* Apply post-EM at variant emission, alongside the existing
  hard-drop reasons
  ([cohort_driver.rs:229–250](../../src/pop_var_caller/cohort_driver.rs#L229)).
* Compute cohort `MQDiffT` per ALT (same value emitted to INFO).
* Drop the record if **any** ALT has `MQDiffT < threshold` **AND**
  both REF and that ALT have enough cohort-pooled reads
  (`cohort_n_obs_ref ≥ 3` and `cohort_n_obs_alt ≥ 3`). Records
  where the test is undefined (insufficient reads on either side)
  are kept — we do not penalize sites for which the metric cannot
  speak.
* Record-count summary log gains a new counter:
  `records_dropped_low_mapq_diff_t`.

### CLI knobs

```text
--no-mapq-diff-filter           # skip the drop entirely; INFO still emits
--min-mapq-diff-t <FLOAT>       # default -3.0; threshold for the drop
```

`--no-mapq-diff-filter` is the equivalent of `--no-complexity-filter`
for this filter — a single boolean override that keeps the CLI
surface aligned. `--min-mapq-diff-t` tunes the threshold;
`-inf` is permitted (no records ever drop on this criterion).

Internal sanity constant (no CLI flag):

```rust
const MAPQ_FILTER_MIN_READS_PER_SIDE: u32 = 3;
```

### Default — rationale

From today's cross-validation, Welch's t restricted to GATK shared
sites, partitioned by GATK INFO buckets:

| `t` threshold | GATK extreme caught | GATK suspect caught | GATK clean false-flagged | extreme:clean ratio |
|--------------:|--------------------:|--------------------:|-------------------------:|--------------------:|
| ≤ −2          |              36.1 % |              23.0 % |                   3.4 %  |  10.6× |
| **≤ −3**      |          **25.0 %** |          **14.4 %** |                **2.0 %** | **12.5×** |
| ≤ −4          |              13.9 % |               9.6 % |                   0.7 %  |  19.9× |
| ≤ −5          |              11.1 % |               6.5 % |                   0.7 %  |  15.9× |

Default `−3.0` is the standard statistical threshold (p < 0.001
under one-sided normal approximation) with a 12.5× ratio of
extreme-caught to clean-false-flagged. A stricter user sets
`-4.0`; a looser user sets `-2.0`. The 2.0 % false-flag on
GATK-clean sites is acceptable because the bench is not the
ground truth in low-mappability regions ([[benchmark-validity-against-gatk]]).

### Caveat

* On the 10-duplicate synthetic cohort, hom-alt sites contribute
  `num_obs_ref = 0` across all samples, so the test is undefined
  there and the filter never triggers. Expected from the analysis
  — the signal lives at het sites with both ref and alt reads.
* DUST stays on (policy choice, [[benchmark-validity-against-gatk]]).
  Where DUST masks a site, MAPQ-diff never sees it.

### Estimated impact on the existing benchmark

From the analysis at `--min-mapq-diff-t = -3.0`:

* ≈ 4.3 % of S0000 FPs with the test defined (97 % of FPs) → ~3,650
  records dropped from the ours-only set.
* ≈ 8.7 % of S0000 TPs with the test defined (51 % of TPs; the
  other half are hom-alt with no ref reads) → ~790 records dropped
  from the shared set.

The benchmark F1 will dip slightly (more TPs lost than the absolute
drop in FPs is worth, because the benchmark itself contains
multi-mapper hotspots GATK also emits). This is intentional;
[[benchmark-validity-against-gatk]] applies — the filter targets
suspect records GATK *also* emits without filtering.

---

## Test plan

1. **AlleleSupportStats round-trip in PSP** — write a record with
   known `mapq_sum` and `mapq_sum_sq`, read it back, assert
   equality.
2. **End-to-end MAPQ tracking** — synthetic BAM with reads at
   known MAPQ values; pile up; assert that aggregated
   `mapq_sum`/`mapq_sum_sq` match expected cohort sums.
3. **VCF INFO emission** — extend
   `tests/cohort_vcf_writer_integration.rs` with a record whose
   per-allele MAPQ stats are pre-computed; assert
   `MQRef`/`MQAlt`/`MQDiff`/`MQDiffT` string emission and missing-
   field handling.
4. **Filter triggers correctly** — construct a record whose
   computed `MQDiffT = −5.0`; assert the record is dropped at the
   default threshold and emitted when `--no-mapq-diff-filter` is
   set or `--min-mapq-diff-t -inf` is passed.
5. **Format version mismatch** — write a v1.0 PSP and attempt to
   read with the v1.1 reader; assert a clean error mentions the
   format version and the regeneration instruction.
6. **Min-reads guard** — record with `n_alt = 2` and `MQDiffT =
   −10`; assert kept (test undefined).

---

## Implementation phases — status (2026-05-22)

**Phase A — landed (commit `a87e9a6`):**
`AlleleSupportStats.mapq_sum: u32` and `mapq_sum_sq: u64` flow from
the BAM walker through `PreparedRead.mapq`, `ReadContribution.mapq`,
`add_contribution`, the per-group-merger's compound-allele rebuilds,
and the PSP I/O (two new columns 0x15, 0x16). New round-trip test
`mapq_scalars_round_trip`.

**Phase B — landed (commit `5903aa0`):**
`INFO/MQRef` (Number=1), `INFO/MQAlt`/`MQDiff`/`MQDiffT` (Number=A)
emitted on every record from cohort-pooled `mapq_sum` /
`mapq_sum_sq`. New integration test
`mapq_info_fields_reflect_cohort_pooled_stats` pins both the header
declarations and an end-to-end expected `MQDiffT ≈ −5.66` on a
hand-built suspect record.

**Phase C — landed (commit `e8f3e90`):**
`--no-mapq-diff-filter` and `--min-mapq-diff-t` flags;
`record_fails_mapq_diff_t` helper at the variant-emission stage;
new `records_dropped_low_mapq_diff_t` counter on the
`var-calling: …` summary line. Tests
`mapq_diff_t_filter_decision_matches_thresholds` cover the drop /
disable / undefined-stay branches.

**Phase D — in progress:** regenerate the synthetic 10-duplicate
cohort PSPs (the new MAPQ columns are required by the reader),
re-run cohort var-calling, compare against the GATK 10-sample
reference, and validate that the actual `MQDiffT` distribution and
drop counter match the analysis prediction
(~3,650 FPs + ~790 TPs dropped at the default).
