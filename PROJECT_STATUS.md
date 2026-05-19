# Project status

<!--
ABOUT-PARAGRAPH-START — do not edit this paragraph.
Skills and agents are instructed to leave it untouched.
-->
> **About this project.** Multi-sample SNP caller:
> per-sample pileup → `.psp` artefact → DUST filter → variant grouping →
> per-group merger → posterior engine. The authoritative design document is
> [doc/devel/specs/calling_pipeline_architecture.md](doc/devel/specs/calling_pipeline_architecture.md);
> read it before anything else. Companion design context is in
> [doc/devel/specs/design_principles.md](doc/devel/specs/design_principles.md).
> All work below is graded against that spec. For the AI assistant's
> per-skill instructions on reading and
> updating this file, see [doc/devel/ia/skills/](doc/devel/ia/skills/) —
> every skill defines a "Project status protocol" section.
<!-- ABOUT-PARAGRAPH-END -->

> **Current focus.** _Maintained by skills (last-completed) and the human
> project manager (next-task)._
>
> - **Last completed task:** Cohort CLI follow-up **Wave 5**
>   (Test infrastructure + missing coverage) fixes-applied
>   2026-05-19 —
>   [cohort_cli_2026-05-19_applied_wave5.md](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave5.md).
>   Final wave of the deferred follow-up; **all 4 Wave-5
>   findings Applied**: **Mi20** (`tests/common/mod.rs`
>   consolidates fixture helpers across the two integration-test
>   binaries), **Mi23** (load-bearing
>   `estimate_contamination → var_calling` chain integration
>   test + 4 missing-coverage unit / integration tests),
>   **M1/M2 follow-up** (end-to-end CRAM-fixture test for the
>   walker-error path — `max_active_reads = 1` trips
>   `Walker(_)`, output VCF cleaned up), and the
>   only-behaviour-change finding **M5 follow-up** (real
>   FASTA → `.psp` per-contig MD5 enforcement via
>   `verify_fasta_matches_psp_chromosomes` + the typed
>   `FastaContigMd5Mismatch` / `FastaContigFetchFailed` error
>   variants on both subcommands; the v1 basename-only
>   contract is now obsolete). 890 lib tests pass (+2 from
>   the new `to_estimates_*` artefact-builder tests) + 9
>   cohort CLI integration tests (+5 net from the new
>   integration tests added this wave); fmt + clippy clean;
>   `cargo doc --no-deps` exits 101 on the **same 5
>   pre-existing errors** as every prior wave; zero introduced
>   by Wave 5.
>   **The 16 originally-Deferred findings from the 2026-05-19
>   review are now closed in full.**
>   Wave reports:
>   [Wave 1](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave1.md)
>   (commit `c7ee0c3`),
>   [Wave 2](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave2.md)
>   (commit `f44c086`),
>   [Wave 3](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave3.md)
>   (commit `248521a`),
>   [Wave 4](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave4.md)
>   (commit `d84ee8e`).
>   The reviewed slice landed earlier on 2026-05-19 in commits
>   `1523049` through `147e435` (impl report:
>   [pop_var_caller_cohort_cli_2026-05-19.md](doc/devel/reports/implementations/pop_var_caller_cohort_cli_2026-05-19.md);
>   review: [cohort_cli_2026-05-19.md](doc/devel/reports/reviews/cohort_cli_2026-05-19.md)),
>   delivering three new subcommands (`var-calling`,
>   `estimate-contamination`, `var-calling-from-bam`). Plans:
>   [pop_var_caller_cohort_cli.md](doc/devel/implementation_plans/pop_var_caller_cohort_cli.md)
>   (original slice),
>   [pop_var_caller_cohort_cli_followup.md](doc/devel/implementation_plans/pop_var_caller_cohort_cli_followup.md)
>   (five-wave deferred follow-up).
> - **Next task:** _set by human PM._ Standing candidates: the
>   manual `bcftools view` / `bcftools stats` smoke against real
>   cohort data (the synthetic fixture in the integration tests is
>   too tiny to be meaningful); write the still-pending Stage 5
>   implementation report; re-bench the full Wave-1 set on a
>   quieter host with a clean pre-perf-review checkout baseline;
>   apply the remaining Hot-path findings from the perf review —
>   H5 (`DEFAULT_BATCH_SIZE` sweep), H6 (Stage 4 bench
>   fixture-rebuild fix), H7 (cohort-size sweep); the
>   parallelisation-tuning pass deferred until after the cohort CLI
>   (rayon-over-records, `--per-group-batch-size` exposure); pick
>   up the standing items below (BED-region skip, phase-chain
>   integration tests).

---

## Pipeline stages

Stage descriptions are one-line reminders; the spec is authoritative.

### Stage 1 — per-sample pileup (BAM → `.psp`)

Stage 1 reads each BAM/CRAM once per sample and writes one `.psp` artefact.

#### CRAM input
- **Status:** shipped
- **Plan:** [per_sample_caller_cram_input.md](doc/devel/implementation_plans/per_sample_caller_cram_input.md)
- **Impl report:** [per_sample_caller_cram_input_2026-04-29.md](doc/devel/reports/implementations/per_sample_caller_cram_input_2026-04-29.md)
- **Latest review / fixes:** [per_sample_caller_cram_input_2026-04-29.md](doc/devel/reports/reviews/per_sample_caller_cram_input_2026-04-29.md), [fixes_applied_2026-05-01.md](doc/devel/reports/implementations/fixes_applied_2026-05-01.md)
- **Open:** none

#### Pileup walker
- **Status:** shipped
- **Code:** [src/per_sample_pileup/pileup/](src/per_sample_pileup/pileup/)
- **Subfeature plans & impl reports:**
  - lazy CIGAR — [plan](doc/devel/implementation_plans/pileup_lazy_cigar.md), [impl](doc/devel/reports/implementations/pileup_lazy_cigar_2026-05-07.md)
  - fold cache — [impl](doc/devel/reports/implementations/pileup_fold_cache_2026-05-07.md)
  - freebayes-style bench — [plan](doc/devel/implementation_plans/pileup_freebayes_style_benchmark.md), [impl](doc/devel/reports/implementations/pileup_freebayes_bench_c_2026-05-08.md)
  - pull-iterator walker — [plan](doc/devel/implementation_plans/pileup_pull_iterator.md), [impl](doc/devel/reports/implementations/pileup_pull_iterator_2026-05-14.md)
  - unique chain ids — [plan](doc/devel/implementation_plans/unique_chain_ids.md), [impl](doc/devel/reports/implementations/unique_chain_ids_2026-05-14.md)
- **Latest reviews:** [pileup_2026-05-11.md](doc/devel/reports/reviews/pileup_2026-05-11.md), [perf_pileup_2026-05-12.md](doc/devel/reports/reviews/perf_pileup_2026-05-12.md)
- **Latest fixes-applied:** [fixes_applied_2026-05-13.md](doc/devel/reports/implementations/fixes_applied_2026-05-13.md)
- **Open:** none

#### BAQ
- **Status:** shipped
- **Plan:** [baq.md](doc/devel/implementation_plans/baq.md)
- **Impl report:** [baq_2026-05-12.md](doc/devel/reports/implementations/baq_2026-05-12.md)
- **Latest review / perf review / fixes:** [baq_2026-05-12.md](doc/devel/reports/reviews/baq_2026-05-12.md), [perf_baq_2026-05-12.md](doc/devel/reports/reviews/perf_baq_2026-05-12.md), [fixes_applied_2026-05-12.md](doc/devel/reports/implementations/fixes_applied_2026-05-12.md)
- **Open:** none

#### Pileup → psp seam
- **Status:** shipped
- **Code:** [src/per_sample_pileup/pileup_to_psp.rs](src/per_sample_pileup/pileup_to_psp.rs)
- **Impl report:** [pileup_to_psp_seam_2026-05-14.md](doc/devel/reports/implementations/pileup_to_psp_seam_2026-05-14.md)
- **Open:** none

#### `pop_var_caller` CLI
- **Status:** Stage 1 CLI shipped (subcommands `pileup`, `psp-to-pileup`);
  cohort CLI **shipped** — every Deferred finding from the
  2026-05-19 review is Applied (Waves 1 – 5).
- **Plans:**
  - Stage 1 CLI (`pileup`, `psp-to-pileup`):
    [pop_var_caller_pileup_cli.md](doc/devel/implementation_plans/pop_var_caller_pileup_cli.md)
  - Cohort CLI (`var-calling`, `estimate-contamination`,
    `var-calling-from-bam`):
    [pop_var_caller_cohort_cli.md](doc/devel/implementation_plans/pop_var_caller_cohort_cli.md)
- **Impl report (cohort slice):**
  [pop_var_caller_cohort_cli_2026-05-19.md](doc/devel/reports/implementations/pop_var_caller_cohort_cli_2026-05-19.md)
- **Latest review (cohort slice):**
  [cohort_cli_2026-05-19.md](doc/devel/reports/reviews/cohort_cli_2026-05-19.md) —
  Request-changes: 0 Blockers, 14 Major (M1–M14), 23 Minor + grouped Nits.
- **Latest fixes-applied (cohort slice):**
  Wave 5 of the deferred follow-up,
  [cohort_cli_2026-05-19_applied_wave5.md](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave5.md) —
  **Mi20**, **Mi23**, **M1/M2 follow-up**, **M5 follow-up**
  (FASTA → `.psp` MD5 enforcement is the only behaviour change
  across the 5-wave plan). 890 lib + 9 cohort integration
  tests pass. Earlier passes:
  [Wave 4](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave4.md)
  (M13);
  [Wave 3](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave3.md)
  (Mi8 / Mi19 / M9-followup / Mi18 / M10 / M11);
  [Wave 2](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave2.md)
  (M4 / Mi2 / Mi14 / Mi21);
  [Wave 1](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave1.md)
  (M8 / Mi5 / Mi6 / Mi13);
  [original pass](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied.md)
  (20 Applied + Mi7 Disputed-with-test).
  **Status: shipped** — every originally-Deferred finding is
  Applied; the 5-wave follow-up plan is complete.
- **Code:** [src/pop_var_caller/](src/pop_var_caller/)
- **Integration tests:**
  [tests/pileup_cli_integration.rs](tests/pileup_cli_integration.rs)
  (Stage 1) and
  [tests/cohort_cli_integration.rs](tests/cohort_cli_integration.rs)
  (cohort subcommands).
- **Open (from cohort-slice review):**
  *None — the 16 originally-Deferred findings are all Applied.*
  - **Closed in Wave 5 (2026-05-19):** **Mi20**, **Mi23**,
    **M1+M2-followup**, **M5-followup**.
  - **Closed in Wave 4 (2026-05-19):** **M13**.
  - **Closed in Wave 3 (2026-05-19):** **Mi8**, **Mi19**,
    **M9-followup**, **Mi18**, **M10**, **M11**.
  - **Closed in Wave 2 (2026-05-19):** **M4**, **Mi2**, **Mi14**,
    **Mi21**.
  - **Closed in Wave 1 (2026-05-19):** **M8**, **Mi5**, **Mi6**,
    **Mi13**.
  - Selected deferred Nits — drop unused `#[from]` variants, vestigial
    `let _ = cfg;`, `Stage1RunSummary` renames, `DEFAULT_BATCH_ID`
    value-mismatch, etc.
  - `bcftools view` / `bcftools stats` manual smoke against real
    cohort data (pre-existing).

---

### Stage 2 — per-sample file (`.psp`) contract

Stage 2 is the on-disk artefact format that Stage 1 writes and Stages 3–6
consume. Not a runtime step — an interface.

#### `.psp` writer + reader
- **Status:** shipped
- **Spec:** [per_sample_pileup_format.md](doc/devel/specs/per_sample_pileup_format.md)
- **Plans:**
  - reader: [psp_reader.md](doc/devel/implementation_plans/psp_reader.md)
  - writer/reader bootstrap: [per_sample_pileup_writer_reader.md](doc/devel/implementation_plans/per_sample_pileup_writer_reader.md)
  - psp → pileup roundtrip: [psp_to_pileup.md](doc/devel/implementation_plans/psp_to_pileup.md)
- **Impl reports:** [psp_reader_2026-05-13.md](doc/devel/reports/implementations/psp_reader_2026-05-13.md), [psp_to_pileup_2026-05-15.md](doc/devel/reports/implementations/psp_to_pileup_2026-05-15.md)
- **Latest reviews:** [psp_2026-05-13.md](doc/devel/reports/reviews/psp_2026-05-13.md), [psp_reader_2026-05-13.md](doc/devel/reports/reviews/psp_reader_2026-05-13.md), [perf_psp_writer_2026-05-13.md](doc/devel/reports/reviews/perf_psp_writer_2026-05-13.md), [perf_psp_reader_2026-05-13.md](doc/devel/reports/reviews/perf_psp_reader_2026-05-13.md)
- **Latest fixes-applied:** [fixes_applied_psp_reader_2026-05-13.md](doc/devel/reports/implementations/fixes_applied_psp_reader_2026-05-13.md), [perf_psp_reader_2026-05-13_applied.md](doc/devel/reports/reviews/perf_psp_reader_2026-05-13_applied.md), [perf_psp_writer_2026-05-13_applied.md](doc/devel/reports/reviews/perf_psp_writer_2026-05-13_applied.md)
- **Open:** none

---

### Stage 3 — low-complexity (sdust) filter

Streaming per-position filter that computes sdust scores from the
reference and silently drops low-complexity records. No intermediate
mask file. Algorithm ported from `lh3/sdust` (vendored at `sdust/`,
gitignored).

- **Status:** fixes-applied
- **Spec section:** `## Stage 3 — low-complexity filter` in [calling_pipeline_architecture.md](doc/devel/specs/calling_pipeline_architecture.md)
- **Plan:** [dust_filter.md](doc/devel/implementation_plans/dust_filter.md)
- **Code:** [src/var_calling/dust_filter.rs](src/var_calling/dust_filter.rs)
- **Tests:** 38 tests in the module (algorithmic core, config
  validation, iterator plumbing, golden vector against committed
  `lh3/sdust` outputs, seeded-random invariant sweep over ~40 K
  input combinations, half-open boundary, threshold-strictness
  boundary, `pos == 0` latch, exhaustion latch). `lh3/sdust` is
  not a build- or test-time dependency.
- **Impl report:** [dust_filter_2026-05-17.md](doc/devel/reports/implementations/dust_filter_2026-05-17.md)
- **Latest review:** [dust_filter_2026-05-17.md](doc/devel/reports/reviews/dust_filter_2026-05-17.md) — Request-changes: 0 Blockers, 8 Major, 24 Minor + 9 Nits.
- **Latest fixes-applied:** [dust_filter_2026-05-17_applied.md](doc/devel/reports/reviews/dust_filter_2026-05-17_applied.md) — 36 of 41 findings Applied (7 of 8 Majors + 22 of 24 Minors + 7 Nits); 1 Applied-with-adaptation (Mi21 partially blocked on parallel `posterior_engine` WIP); 2 Won't-fix per project preference (M4 out-of-order check is the merger's responsibility; Mi6 informational tracing conflicts with the "no logs" project preference); 2 Nits won't-fix.
- **Open:**
  - **Mi21** — full `cargo doc --no-deps --all-features` with `-D warnings` to be re-run once the parallel `posterior_engine` refactor compiles; intra-doc links in `dust_filter.rs` itself have been visually verified.
  - Criterion bench — deferred until the cohort CLI runs against
    real cohort data and we can measure DUST's share of real
    wall time.
  - **Closed 2026-05-19** (cohort CLI slice): CLI parser bindings
    (`--complexity-window`, `--complexity-threshold`,
    `--no-complexity-filter`) shipped with `var-calling`;
    end-to-end integration test exists in
    [tests/cohort_cli_integration.rs](tests/cohort_cli_integration.rs).

---

### Stage 4 — grouping

Walks the filtered per-position stream and bundles overlapping / adjacent
candidate positions into `OverlappingVarGroup`s. Sequential, but emits
independent groups that downstream stages can process in parallel.

#### Multi-way per-position iterator
The k-way merge that produces the per-position stream consumed by the
grouper.
- **Status:** shipped
- **Plan:** [multi_way_per_position_iterator.md](doc/devel/implementation_plans/multi_way_per_position_iterator.md)
- **Impl report:** [multi_way_per_position_iterator_2026-05-15.md](doc/devel/reports/implementations/multi_way_per_position_iterator_2026-05-15.md)
- **Code:** [src/var_calling/per_position_merger.rs](src/var_calling/per_position_merger.rs)
- **Bench:** `var_calling_merger/*` in [benches/var_calling_perf.rs](benches/var_calling_perf.rs)
- **Latest reviews:** [per-position-merger_2026-05-15.md](doc/devel/reports/reviews/per-position-merger_2026-05-15.md), [perf_var_calling_2026-05-16.md](doc/devel/reports/reviews/perf_var_calling_2026-05-16.md)
- **Latest fixes-applied:** bundled into [fixes_applied_2026-05-16.md](doc/devel/reports/implementations/fixes_applied_2026-05-16.md) (cohort run) + [fixes_applied_2026-05-16_v2.md](doc/devel/reports/implementations/fixes_applied_2026-05-16_v2.md) (Mi14/Mi15 OutOfOrder+ChromosomeMismatch land here).
- **Open:**
  - **H6** (perf review): the Stage 4 merger bench is fixture-dominated
    (67 % inclusive in `iter_batched` setup vs 21 % in
    `PerPositionMerger::next`); rebuild fixtures once outside `b.iter`
    before any Stage 4 code-level change can be defended.
  - **L8** (perf review, depends on H6): `vec![None; self.n_samples()]`
    per emission at [per_position_merger.rs:306](src/var_calling/per_position_merger.rs#L306)
    is the candidate Stage 4 allocation site to revisit once the bench
    is unblocked.

#### Variant grouper (Stage 4 bundler)
- **Status:** fixes-applied (2026-05-16); no separate implementation report
- **Plan:** [cohort_variant_grouping.md](doc/devel/implementation_plans/cohort_variant_grouping.md)
- **Code:** [src/var_calling/variant_grouping.rs](src/var_calling/variant_grouping.rs)
- **Bench:** `var_calling_grouper/*` in [benches/var_calling_perf.rs](benches/var_calling_perf.rs)
- **Latest reviews:** [cohort_2026-05-16.md](doc/devel/reports/reviews/cohort_2026-05-16.md), [perf_var_calling_2026-05-16.md](doc/devel/reports/reviews/perf_var_calling_2026-05-16.md)
- **Latest fixes-applied:** [fixes_applied_2026-05-16.md](doc/devel/reports/implementations/fixes_applied_2026-05-16.md) + [fixes_applied_2026-05-16_v2.md](doc/devel/reports/implementations/fixes_applied_2026-05-16_v2.md) (v2 closes the entire v1 deferred list)
- **Open:**
  - Implementation report not yet written.
  - **H6** (perf review): the Stage 4 grouper bench is also
    fixture-dominated (77 % inclusive in `iter_batched` setup vs 13 %
    in `VariantGrouper::next`); rebuild fixtures once outside `b.iter`
    before any Stage 4 code-level change can be defended. Shared
    bench-shape fix with the merger above.

---

### Stage 5 — per-group processing

Allele unification + per-sample likelihood reconstruction + phase-chain
consistency on one `OverlappingVarGroup` at a time. Parallel across groups
via rayon.

#### Per-group merger (allele unification + likelihood)
- **Status:** fixes-applied (2026-05-16)
- **Plan:** [cohort_per_group_merger.md](doc/devel/implementation_plans/cohort_per_group_merger.md)
- **Code:** [src/var_calling/per_group_merger.rs](src/var_calling/per_group_merger.rs)
- **Bench:** `var_calling_per_group_merger/*` in [benches/var_calling_perf.rs](benches/var_calling_perf.rs).
  Slowest stage per element on the criterion baseline (~22 K groups/s
  biallelic, ~12 K groups/s compound) and the primary focus of the
  2026-05-16 perf review.
- **Impl report:** not yet saved.
- **Latest reviews:** [cohort_2026-05-16.md](doc/devel/reports/reviews/cohort_2026-05-16.md), [perf_var_calling_2026-05-16.md](doc/devel/reports/reviews/perf_var_calling_2026-05-16.md)
- **Latest fixes-applied:** [fixes_applied_2026-05-16.md](doc/devel/reports/implementations/fixes_applied_2026-05-16.md) + [fixes_applied_2026-05-16_v2.md](doc/devel/reports/implementations/fixes_applied_2026-05-16_v2.md) (v2 closes the v1 deferred list) + [perf_var_calling_2026-05-16_applied.md](doc/devel/reports/reviews/perf_var_calling_2026-05-16_applied.md) (Wave 1 of the perf review: L2 + H1 + H4 + H3 + H2 applied — H2 phase 2b is the headline ~14 % biallelic win, p = 0.00, confirmed by two back-to-back runs; all 571 tests pass)
- **Open:**
  - Implementation report for the Stage 5 merger has not been saved; the
    next `rust-feature-implementation` run for this feature should produce
    one and link it here.
  - Phase-chain integration tests for the likelihood calculation (see
    *Standing items* below).
  - Re-bench Wave 1 on a quieter host with a clean pre-perf-review
    checkout baseline to put an absolute number on the cumulative
    L2 + H1 + H4 + H3 + H2 effect. The in-session ~12–14 %
    biallelic_64 figure is reliable for H2 (back-to-back confirms);
    the absolute pre-vs-post number drifted with host load. See
    [perf_var_calling_2026-05-16_applied.md](doc/devel/reports/reviews/perf_var_calling_2026-05-16_applied.md).
  - Remaining Hot-path findings to apply (see
    [perf_var_calling_2026-05-16.md](doc/devel/reports/reviews/perf_var_calling_2026-05-16.md)
    §5):
    - **H5** — sweep `DEFAULT_BATCH_SIZE` (currently a self-declared
      placeholder at [per_group_merger.rs:55](src/var_calling/per_group_merger.rs#L55))
      and reset to the measured optimum.
    - **H6** — Stage 4 bench fixture-rebuild fix (the per_position_merger
      / variant_grouper benches' iter_batched setup is fixture-dominated
      per the original perf review). Unblocks Stage 4 code-level findings.
    - **H7** — cohort-size sweep at N=10/64/256/1024 samples.
  - Wave 2 / 3 / Likely / Speculative findings tracked in the report.

---

### Stage 6 — posterior engine

EM over merged records → final multi-sample VCF.

#### Posterior engine
- **Status:** fixes-applied
- **Spec sections:** `## Stage 6 — posterior engine` in [calling_pipeline_architecture.md](doc/devel/specs/calling_pipeline_architecture.md); background in [freebayes_posterior_gt_probs.md](doc/devel/specs/freebayes_posterior_gt_probs.md) and [gatk_em_calculation.md](doc/devel/specs/gatk_em_calculation.md)
- **Plan:** [posterior_engine.md](doc/devel/implementation_plans/posterior_engine.md)
- **Code:** [src/var_calling/posterior_engine.rs](src/var_calling/posterior_engine.rs)
- **Tests:** unit tests in the module + [tests/posterior_engine_integration.rs](tests/posterior_engine_integration.rs)
- **Impl reports:** [posterior_engine_2026-05-16.md](doc/devel/reports/implementations/posterior_engine_2026-05-16.md); perf history: [perf 2026-05-17](doc/devel/reports/implementations/posterior_engine_perf_2026-05-17.md), [SIMD analysis 2026-05-18](doc/devel/reports/implementations/posterior_engine_simd_analysis_2026-05-18.md), [perf wave-1 2026-05-18](doc/devel/reports/implementations/posterior_engine_perf_2026-05-18.md), [post-H4 profile](doc/devel/reports/implementations/posterior_engine_post_h4_profile_2026-05-18.md), [perf wave-2 (H4 + RecordScratch) 2026-05-18](doc/devel/reports/implementations/posterior_engine_perf_2026-05-18_v2.md)
- **Latest reviews:** [posterior_engine_2026-05-16.md](doc/devel/reports/reviews/posterior_engine_2026-05-16.md) — Request-changes: 2 Blocker test-gap findings, 11 Major, 13 Minor; [perf_posterior_engine_2026-05-18.md](doc/devel/reports/reviews/perf_posterior_engine_2026-05-18.md) — Run-experiments: 4 Hot-path (H4 homogeneous-fixation hoist subsumes H1 SIMD `log_sum_exp_2_x4 -INF` short-circuit; H2 `log_indep` cross-batch reuse; H3 GATK-style natural-log `softplus_neg` table for `log_sum_exp_2_x4`), 11 Likely (per-engine `RecordScratch` lift for 13 per-record allocs, `log_likelihoods` batch-of-4 transpose, autovec-blocking branch in `accumulate_expected_counts`, `#[inline(always)]` / `#[cold]` audits, in-bench self-validation), 11 Speculative (SQUAREM, warm-start across records, adaptive genotype pruning à la Octopus, cross-record tensor-style batching at same-shape, `f64x8` AVX-512 backend, PGO, `x86-64-v4` floor, …). Built on the May 17/18 bench evidence; the deep web + local-source (`gatk/`, `freebayes/`, `bcftools/`) research that fed it surfaced that none of the three reference callers vectorise ln/exp (GATK's `JacobianLogTable` is the closest precedent for H3) — `wide::f64x4`-native ln/exp already leads the field.
- **Latest fixes-applied:** [posterior_engine_2026-05-16_applied.md](doc/devel/reports/reviews/posterior_engine_2026-05-16_applied.md) — both Blockers fixed (B1 `NonFinitePosterior` test, B2 trivial-record row-sum bug + tests); 10 of 11 Majors applied (M3 config validation deferred — needs engine-vs-CLI boundary decision); 8 of 13 Minors applied (Mi5 perf bench, Mi6 golden test, Mi12 `with_config -> Result`, Mi13 fixture consolidation deferred as standalone follow-ups). 622 lib + 7 integration tests pass, clippy-clean on the in-scope files.
- **Open:**
  - **Mi5** — `benches/posterior_engine_perf.rs` regression-threshold criterion bench.
  - **Mi6** — golden `tests/golden/posterior_engine/*` fixtures locking emitted `f64` values across the representative-record matrix.
  - **Mi13** — fixture consolidation across the unit / integration crate boundary.
  - **Closed 2026-05-19** (cohort CLI slice): **M3** (engine-side
    config validation for `F`, pseudocounts) and **Mi12**
    (`Config::new -> Result<Self, _>`) — both shipped as
    `PosteriorEngineConfig::new` with full range validation. The
    CLI surface mirrors the same ranges in
    [src/pop_var_caller/cli/parsers.rs](src/pop_var_caller/cli/parsers.rs).
  - **Perf — applied 2026-05-18 (closed):** H1 (subsumed by H4), H2 (`log_indep` cross-batch reuse), H4 (homogeneous-fixation hoist) — see commit `013b49f`. L5–L7 (per-engine `RecordScratch` lift) — see commit `9594533`. Full session report at [posterior_engine_perf_2026-05-18_v2.md](doc/devel/reports/implementations/posterior_engine_perf_2026-05-18_v2.md).
  - **Perf — refuted by hardware-counter profile (demoted to Note):** L8/L9 (`log_likelihoods` / `posteriors` batch-of-4 transpose). Post-RecordScratch IPC = 2.97, L1 hit 99.66 %, branch-miss 0.11 % — the engine is execution-bound, not memory-bound. Cache layout has no remaining room to help.
  - **Perf — demoted to Speculative:** H3 (GATK-style `softplus_neg` table). H4 removed most `log_sum_exp_2_x4` calls; `__ieee754_log_fma` is now only ~3 % of cycles.
  - **Perf — promoted to Likely:** S5 (AVX-512 `f64x8` backend). With IPC=2.97 single-thread saturated, wider lanes are one of the few remaining single-thread levers (~2× on Skylake-X+/Zen 4-5 hosts).
  - **Perf next big lever — rayon-over-records.** Order-of-magnitude lever (16–128× on a server). RecordScratch lift is structured for per-thread ownership: each worker gets its own scratch. Separate plan; deferred from the wave-2 session per user's call.
  - **Perf remaining small-but-cheap:** L3 (`#[inline(always)]` on `log_sum_exp_*_x4`, gated by `cargo asm` check), L10 (in-bench `debug_assert!` on `diagnostics.iterations`).
  - **Closed 2026-05-19** (cohort CLI slice): end-to-end
    PspReader → … → PosteriorEngine integration test lives in
    [tests/cohort_cli_integration.rs](tests/cohort_cli_integration.rs).

#### Contamination-estimation side-pass
- **Status:** fixes-applied
- **Spec:** [contamination_estimation.md](doc/devel/specs/contamination_estimation.md)
- **Plan:** [contamination_estimation.md](doc/devel/implementation_plans/contamination_estimation.md)
- **Code:** [src/var_calling/contamination_estimation.rs](src/var_calling/contamination_estimation.rs);
  Stage 6 consumer-side wiring in
  [src/var_calling/posterior_engine.rs](src/var_calling/posterior_engine.rs)
- **Tests:** unit + 3 proptests in the module +
  [tests/contamination_estimation_integration.rs](tests/contamination_estimation_integration.rs)
- **Impl report:** [contamination_estimation_2026-05-17.md](doc/devel/reports/implementations/contamination_estimation_2026-05-17.md)
- **Latest review:** [contamination_estimation_2026-05-17.md](doc/devel/reports/reviews/contamination_estimation_2026-05-17.md) — Request-changes: 2 Blockers, 16 Major, 23 Minor + grouped Nits.
- **Latest fixes-applied:** [fixes_applied_2026-05-17.md](doc/devel/reports/reviews/fixes_applied_2026-05-17.md) — 33 Applied (both Blockers + 14 of 16 Majors + 17 of 23 Minors), 2 Applied-with-adaptation (M15 cohort_alleles per-position alloc instead of cross-position scratch; Mi11 INF-delta test inverted), 1 Already-fixed (Mi17 covered by M9), 6 Deferred (M13 / M14 / Mi9 / Mi19 / Mi21 / Mi22). 668 lib + 109 integration tests pass; clippy + fmt clean. Perf comparison: 6 of 8 benches unchanged or improved, 2 regressed (var_calling_merger/dense +1.99 %, var_calling_grouper/overlap_extension +11.31 %) — both kept per skill policy (all fixes are correctness/test-coverage grade), most likely code-layout / inlining drift from the added test fixtures in adjacent modules.
- **Open:**
  - **M13** — `..Default::default()` in 14 new test sites; mechanical fix as a focused PR before any new `ContaminationEstimationConfig` field lands. See fixes_applied §M13.
  - **M14** — `OnlineEmState` god-struct split (3 sub-structs: `EmRunningStats`, `EmFrozenParameters`, `StabilityTracker`); pair with Mi21 (`estimate_contamination` long flat state machine). See fixes_applied §M14.
  - **Mi9** — Test fixture duplication between in-module and integration tests; needs cross-crate shared test module design.
  - **Mi19** — `DEFAULT_*` constants lack spec citations; pure-doc churn across 5 sites.
  - **Mi21 / Mi22** — paired with M14 and M13 respectively.
  - Partition-parallel scans across chromosomes (spec §"Cost and
    parallelism") — add only if real cohort runs show wall time
    matters.
  - `--external-allele-frequencies` for substituting a reference-
    panel `q_b` (v2 follow-up).
  - `--contamination-estimates`-only short path that estimates
    `q_b` while freezing user-supplied `c_s`
    (`ContaminationEstimateSource::Mixed` slot exists but is unused
    in v1).
  - **Closed 2026-05-19** (cohort CLI slice): CLI parser bindings
    shipped with `pop_var_caller estimate-contamination`.
  - Compound-allele `q_b = 0` choice (Assumption 2 in the impl
    report) — calibration impact in real cohorts is open.
  - Two perf regressions surfaced by the fix-application run
    (`var_calling_merger/dense`, `var_calling_grouper/overlap_extension`)
    — likely code-layout effects from new code in adjacent modules;
    revisit if real-cohort wall time matters.

#### Cohort VCF writer (Stage 6 sink)
- **Status:** fixes-applied
- **Plan:** [cohort_vcf_writer.md](doc/devel/implementation_plans/cohort_vcf_writer.md)
- **Code:** [src/var_calling/vcf_writer/](src/var_calling/vcf_writer/)
- **Tests:** 41 unit tests in the module + 4 integration tests in
  [tests/cohort_vcf_writer_integration.rs](tests/cohort_vcf_writer_integration.rs).
- **Impl report:** [cohort_vcf_writer_2026-05-18.md](doc/devel/reports/implementations/cohort_vcf_writer_2026-05-18.md)
- **Latest review:** [cohort_vcf_writer_2026-05-18.md](doc/devel/reports/reviews/cohort_vcf_writer_2026-05-18.md)
- **Latest fixes-applied:** [cohort_vcf_writer_2026-05-18_applied.md](doc/devel/reports/reviews/cohort_vcf_writer_2026-05-18_applied.md) — 31 Applied (1 Blocker + 14 of 15 Majors + 16 of 19 Minors + Nits) + 1 Applied-with-adaptation (Mi16) + 4 Deferred.
- **Open:**
  - `bcftools view` / `bcftools stats` manual smoke against real
    cohort data — synthetic integration-test fixture is too tiny to
    be meaningful.
  - **Closed 2026-05-19** (cohort CLI slice): end-to-end exercise
    via `pop_var_caller var-calling` lives in
    [tests/cohort_cli_integration.rs](tests/cohort_cli_integration.rs).
  - Tabix `.tbi` index alongside `.vcf.gz` — out of v1 scope; add
    when random-access matters.
  - `PL` (phred-scaled likelihoods) FORMAT field — needs Stage 5 →
    `PosteriorRecord` forwarding of `log_likelihoods`.
  - Per-sample contamination fraction in INFO — wires in once the
    cohort CLI threads `ContaminationEstimates` into the writer.
  - **M15** — `benches/vcf_writer_perf.rs` (criterion bench at
    1000 samples × biallelic SNP / `emit_gp = true` /
    encode-only); deferred from the review fix pass.
  - **Mi11 + Mi12** — `WriterConfig` → `CohortVcfWriterConfig`
    + `tool_string` → `source_label`/`tool_name` public-API
    renames; pair into a coordinated naming pass before
    publishing the crate.
  - **Mi14** — three-way test-fixture deduplication via a
    `pub(crate) test_fixtures` module + `test-support` feature
    flag; deferred.
  - **M13 follow-up** — fault-injection test for `new()` tmp
    cleanup; lands when the sink-injection seam exists.
  - **Mi16 follow-up** — full `BGZF_EOF` deduplication
    (currently shared in-crate but the integration test keeps
    its own copy).

#### Posterior engine — approximate-LUT inner loop
- **Status:** not yet implemented (config flag wired only)
- **Plan section:** [posterior_engine.md §"Approximation via precomputed lookup tables"](doc/devel/implementation_plans/posterior_engine.md)
- **Open:** evaluation methodology pinned in the plan; needs the
  exact-math engine bench numbers before deciding which candidates
  to land.

---

## Standing project-wide items

Not bound to a specific stage block; pick up whenever the active feature
list is clear.

- **BED-region skip.** Add a Stage 1 CLI flag that takes a BED file of
  regions to skip (avoid pathological regions). Source:
  [doc/devel/TODO.txt](doc/devel/TODO.txt).
- **Phase-chain integration tests.** Add integration tests asserting that
  phase chains are correctly carried through the Stage 5 likelihood
  calculation. Source: [doc/devel/TODO.txt](doc/devel/TODO.txt).
- **Parallel-optimization integration perf benches.** Build end-to-end
  criterion (or equivalent) benches that measure wall time for the two
  pipeline arms whose parallelism still needs tuning: CRAM → `.psp`
  (Stage 1, `pileup` / `var-calling-from-bam`) and `.psp` → cohort VCF
  (Stages 3–6, `var-calling`). These integration benches are the
  ground truth for the deferred parallelisation-tuning pass
  (rayon-over-records, `--per-group-batch-size`, per-group batch
  sizing) — micro-benches alone can't catch end-to-end scaling
  artefacts.
  - **`.psp` → cohort VCF arm — shipped 2026-05-19:**
    [benches/cohort_e2e_perf.rs](benches/cohort_e2e_perf.rs). Two
    bench-group families: `cohort_e2e_core/*` times
    [`drive_cohort_pipeline`](src/pop_var_caller/cohort_driver.rs) in
    isolation (PSP open + FASTA verify outside the timed region) with
    sub-groups `scaling_samples` (N ∈ {10, 64, 256}), `scaling_region`
    (L ∈ {1 000, 5 000, 20 000}), and `scaling_threads` (T ∈ {1, 2,
    4, max-cores} via per-bench local `rayon::ThreadPool` +
    `pool.install(...)`); `cohort_e2e_full/*` times
    [`run_var_calling`](src/pop_var_caller/var_calling.rs) end-to-end
    over `scaling_samples` + `scaling_region`. The full group cannot
    sweep thread count within one `cargo bench` invocation —
    `configure_rayon_pool` / `ThreadPoolBuilder::build_global` is
    once-per-process; run separate invocations under
    `RAYON_NUM_THREADS=N`. Bench-surface lift:
    [`drive_cohort_pipeline`](src/pop_var_caller/cohort_driver.rs) +
    `CohortPipelineParams` promoted from `pub(crate)` to
    `#[doc(hidden)] pub`. 15 bench variants pass `cargo bench --bench
    cohort_e2e_perf -- --test`.
  - **CRAM → `.psp` arm:** still pending.

