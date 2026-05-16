# Project status

<!--
ABOUT-PARAGRAPH-START — do not edit this paragraph.
Skills and agents are instructed to leave it untouched.
-->
> **About this project.** Multi-sample SNP caller, mid-pivot from a
> gVCF-merger to a six-stage Bayesian pipeline:
> per-sample caller → `.psp` artefact → DUST filter → variant grouping →
> per-group merger → posterior engine. The authoritative design document is
> [doc/devel/specs/calling_pipeline_architecture.md](doc/devel/specs/calling_pipeline_architecture.md);
> read it before anything else. Companion design context is in
> [doc/devel/specs/design_principles.md](doc/devel/specs/design_principles.md).
> All work below is graded against that spec, not against the older
> gVCF-merger code still present in `src/` (see *Legacy* section at the
> bottom). For the AI assistant's per-skill instructions on reading and
> updating this file, see [doc/devel/ia/skills/](doc/devel/ia/skills/) —
> every skill defines a "Project status protocol" section.
<!-- ABOUT-PARAGRAPH-END -->

> **Current focus.** _Maintained by skills (last-completed) and the human
> project manager (next-task)._
>
> - **Last completed task:** Cohort criterion bench scaffold landed on
>   2026-05-16 ([benches/cohort_perf.rs](benches/cohort_perf.rs)); the
>   perf review itself was *deferred* because the session's seccomp
>   filter blocks `perf_event_open` and no sampling profile could be
>   captured. Baseline criterion numbers saved at
>   [tmp/perf_review_2026-05-16_cohort/baseline_criterion.txt](tmp/perf_review_2026-05-16_cohort/baseline_criterion.txt).
> - **Next task:** _set by human PM._ Standing candidates: re-run the
>   perf review once a sampling profiler is available (run benches
>   outside the sandbox; pair with samply / cargo flamegraph); write the
>   Stage 5 implementation report; plan Stage 3 (DUST filter); plan
>   Stage 6 (posterior engine); pick up the standing items below
>   (BED-region skip, phase-chain integration tests).

---

## Pipeline stages

Stage descriptions are one-line reminders; the spec is authoritative.

### Stage 1 — per-sample caller (BAM → `.psp`)

Stage 1 reads each BAM/CRAM once per sample and writes one `.psp` artefact.

#### CRAM input
- **Status:** shipped
- **Plan:** [per_sample_caller_cram_input.md](doc/devel/implementation_plans/per_sample_caller_cram_input.md)
- **Impl report:** [per_sample_caller_cram_input_2026-04-29.md](doc/devel/reports/implementations/per_sample_caller_cram_input_2026-04-29.md)
- **Latest review / fixes:** [reviews/per_sample_caller_cram_input_2026-04-29.md](doc/devel/reviews/per_sample_caller_cram_input_2026-04-29.md), [fixes_applied_2026-05-01.md](doc/devel/reviews/fixes_applied_2026-05-01.md)
- **Open:** none

#### Pileup walker
- **Status:** shipped
- **Code:** [src/per_sample_caller/pileup/](src/per_sample_caller/pileup/)
- **Subfeature plans & impl reports:**
  - lazy CIGAR — [plan](doc/devel/implementation_plans/pileup_lazy_cigar.md), [impl](doc/devel/reports/implementations/pileup_lazy_cigar_2026-05-07.md)
  - fold cache — [impl](doc/devel/reports/implementations/pileup_fold_cache_2026-05-07.md)
  - freebayes-style bench — [plan](doc/devel/implementation_plans/pileup_freebayes_style_benchmark.md), [impl](doc/devel/reports/implementations/pileup_freebayes_bench_c_2026-05-08.md)
  - pull-iterator walker — [plan](doc/devel/implementation_plans/pileup_pull_iterator.md), [impl](doc/devel/reports/implementations/pileup_pull_iterator_2026-05-14.md)
  - unique chain ids — [plan](doc/devel/implementation_plans/unique_chain_ids.md), [impl](doc/devel/reports/implementations/unique_chain_ids_2026-05-14.md)
- **Latest reviews:** [pileup_2026-05-11.md](doc/devel/reviews/pileup_2026-05-11.md), [perf_pileup_2026-05-12.md](doc/devel/reviews/perf_pileup_2026-05-12.md)
- **Latest fixes-applied:** [fixes_applied_2026-05-13.md](doc/devel/reviews/fixes_applied_2026-05-13.md)
- **Open:** none

#### BAQ
- **Status:** shipped
- **Plan:** [baq.md](doc/devel/implementation_plans/baq.md)
- **Impl report:** [baq_2026-05-12.md](doc/devel/reports/implementations/baq_2026-05-12.md)
- **Latest review / perf review / fixes:** [baq_2026-05-12.md](doc/devel/reviews/baq_2026-05-12.md), [perf_baq_2026-05-12.md](doc/devel/reviews/perf_baq_2026-05-12.md), [fixes_applied_2026-05-12.md](doc/devel/reviews/fixes_applied_2026-05-12.md)
- **Open:** none

#### Pileup → psp seam
- **Status:** shipped
- **Code:** [src/per_sample_caller/pileup_to_psp.rs](src/per_sample_caller/pileup_to_psp.rs)
- **Impl report:** [pileup_to_psp_seam_2026-05-14.md](doc/devel/reports/implementations/pileup_to_psp_seam_2026-05-14.md)
- **Open:** none

#### `pop_var_caller` CLI
- **Status:** shipped (subcommands `pileup`, `psp-to-pileup`)
- **Plan:** [pop_var_caller_pileup_cli.md](doc/devel/implementation_plans/pop_var_caller_pileup_cli.md)
- **Code:** [src/pop_var_caller/](src/pop_var_caller/)
- **Open:** none

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
- **Latest reviews:** [psp_2026-05-13.md](doc/devel/reviews/psp_2026-05-13.md), [psp_reader_2026-05-13.md](doc/devel/reviews/psp_reader_2026-05-13.md), [perf_psp_writer_2026-05-13.md](doc/devel/reviews/perf_psp_writer_2026-05-13.md), [perf_psp_reader_2026-05-13.md](doc/devel/reviews/perf_psp_reader_2026-05-13.md)
- **Latest fixes-applied:** [fixes_applied_psp_reader_2026-05-13.md](doc/devel/reviews/fixes_applied_psp_reader_2026-05-13.md), [perf_psp_reader_2026-05-13_applied.md](doc/devel/reviews/perf_psp_reader_2026-05-13_applied.md), [perf_psp_writer_2026-05-13_applied.md](doc/devel/reviews/perf_psp_writer_2026-05-13_applied.md)
- **Open:** none

---

### Stage 3 — low-complexity (DUST) filter

Streaming per-position filter that computes DUST from the reference and
silently drops low-complexity records. No intermediate mask file.

- **Status:** not yet planned
- **Spec section:** `## Stage 3 — low-complexity filter` in [calling_pipeline_architecture.md](doc/devel/specs/calling_pipeline_architecture.md)
- **Plan:** none yet
- **Code:** none yet
- **Open:** plan needs to be written before implementation can start.

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
- **Code:** [src/cohort/per_position_merger.rs](src/cohort/per_position_merger.rs)
- **Bench:** `cohort_merger/*` in [benches/cohort_perf.rs](benches/cohort_perf.rs)
- **Latest review:** [reviews/per-position-merger_2026-05-15.md](reviews/per-position-merger_2026-05-15.md)
- **Latest fixes-applied:** bundled into [reviews/fixes_applied_2026-05-16.md](reviews/fixes_applied_2026-05-16.md) (cohort run).
- **Open:** perf review deferred — sampling profile blocked by the
  Claude Code session sandbox (see baseline numbers in
  [tmp/perf_review_2026-05-16_cohort/baseline_criterion.txt](tmp/perf_review_2026-05-16_cohort/baseline_criterion.txt)).

#### Variant grouper (Stage 4 bundler)
- **Status:** fixes-applied (2026-05-16); no separate implementation report
- **Plan:** [cohort_variant_grouping.md](doc/devel/implementation_plans/cohort_variant_grouping.md)
- **Code:** [src/cohort/variant_grouping.rs](src/cohort/variant_grouping.rs)
- **Bench:** `cohort_grouper/*` in [benches/cohort_perf.rs](benches/cohort_perf.rs)
- **Latest review:** [reviews/cohort_2026-05-16.md](reviews/cohort_2026-05-16.md)
- **Latest fixes-applied:** [reviews/fixes_applied_2026-05-16.md](reviews/fixes_applied_2026-05-16.md)
- **Open:** implementation report not yet written; perf review
  deferred (sampling profile blocked — see Stage 3 entry above).

---

### Stage 5 — per-group processing

Allele unification + per-sample likelihood reconstruction + phase-chain
consistency on one `OverlappingVarGroup` at a time. Parallel across groups
via rayon.

#### Per-group merger (allele unification + likelihood)
- **Status:** fixes-applied (2026-05-16)
- **Plan:** [cohort_per_group_merger.md](doc/devel/implementation_plans/cohort_per_group_merger.md)
- **Code:** [src/cohort/per_group_merger.rs](src/cohort/per_group_merger.rs)
- **Bench:** `cohort_per_group_merger/*` in [benches/cohort_perf.rs](benches/cohort_perf.rs).
  Slowest stage per element on the criterion baseline; the natural
  starting point for the deferred perf review.
- **Impl report:** not yet saved.
- **Latest review:** [reviews/cohort_2026-05-16.md](reviews/cohort_2026-05-16.md)
- **Latest fixes-applied:** [reviews/fixes_applied_2026-05-16.md](reviews/fixes_applied_2026-05-16.md)
- **Open:**
  - Implementation report for the Stage 5 merger has not been saved; the
    next `rust-feature-implementation` run for this feature should produce
    one and link it here.
  - Phase-chain integration tests for the likelihood calculation (see
    *Standing items* below).
  - Perf review deferred — sampling profile blocked by the Claude Code
    session sandbox. Baseline numbers + candidate sites recorded in
    [tmp/perf_review_2026-05-16_cohort/baseline_criterion.txt](tmp/perf_review_2026-05-16_cohort/baseline_criterion.txt).
    Re-run the perf-review skill once a sampling profiler is available
    (outside the sandbox, with `samply` / `cargo flamegraph`).

---

### Stage 6 — posterior engine

EM over merged records → final multi-sample VCF.

- **Status:** not yet planned
- **Spec sections:** `## Stage 6 — posterior engine` in [calling_pipeline_architecture.md](doc/devel/specs/calling_pipeline_architecture.md); background in [freebayes_posterior_gt_probs.md](doc/devel/specs/freebayes_posterior_gt_probs.md) and [gatk_em_calculation.md](doc/devel/specs/gatk_em_calculation.md)
- **Plan:** none yet
- **Code:** none yet
- **Open:** plan needs to be written before implementation can start.

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

---

## Legacy / superseded code

Pre-pivot gVCF-merger components, kept until the new pipeline reaches
feature parity and the old CLI paths are retired. **Not in scope** of new
work unless explicitly requested.

- [src/gvcf_parser.rs](src/gvcf_parser.rs)
- [src/genotype_merging.rs](src/genotype_merging.rs)
- [src/genotype_posteriors.rs](src/genotype_posteriors.rs)
- [src/vcf_writer.rs](src/vcf_writer.rs)
