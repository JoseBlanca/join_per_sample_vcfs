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
> - **Last completed task (2026-07-06):** **SSR interrupted-repeat recall — 11-category multi-agent code review (branch `ssr-interruptions`, `550b691`…`4b86518`)**
>   Reviewed the full feature diff (Phase 1 sequence-keyed SSR alleles + Phase 2 purity→level). **Approve-with-changes:** 0 Blocker, **1 Major** (M1 — `nearest_called_by_sequence` runs `align_subst` unconditionally on the length-nearest allele: a per-slip-read `banded_forward` DP added to the whole-cohort θ_locus refit, and for composition-asymmetric *equidistant different-length* reads it can flip the slip `Δ` sign vs the old `nearest_parent`; one fix — short-circuit single-nearest + scope the composition tie-break to same-length ties — covers both), 11 Minor (mostly doc/test-coverage/tuning-knob formalism), Nits. The two hard contracts hold: cross-thread byte-identity is fully traced clean (`fit_purity_level` sorts before its float reduce; all tie-breaks lowest-index) and panic-freedom is "unusually disciplined". Diff matches the spec; the 7 "already conformant" functions verified untouched. Report: [ssr_interrupted_repeat_recall_2026-07-06.md](doc/devel/reports/reviews/ssr_interrupted_repeat_recall_2026-07-06.md). **Next (suggested follow-up: apply fixes)** — M1 + Mi1/Mi2 + the review's missing-tests set.
> - **Last completed task (2026-07-05):** **SFS genotype prior — Milestone 4: empirical-Bayes leave-one-out cohort prior, validated + accepted (branch `sfs-genotype-prior`, tip `ae6aa03`)**
>   The `p̂`-independent DM prior (M1–2) dropped cohort **strength-borrowing**; M4 restores it via the **leave-one-out** frequency update `α'_s = α_species + (E[cohort allele copies] − E[own_s])`, marginalised by the existing DM primitive and reusing the EM's `expected_counts`. **Single-sample is byte-identical to G5 by construction** (with one sample the LOO term is 0 → `α_species`), so `dispatch_e_step` routes `n_samples==1` through the unchanged §9 path — GIAB unchanged (5× SNP concordance **94.9 %**, precision/recall/FP identical). A **flat first EM step** for cohorts escapes the hom-ref fixed-point trap; the convergence guard forbids converging on it. Works at any `(ploidy, n_alleles)`; the ignored cohort test flipped on + multiallelic/polyploid cohort-sharpening tests added. Reviewed by the 4-category orchestration (no Blockers; 2 Major test-coverage gaps + minors all applied — [sfs_eb_cohort_prior_2026-07-05.md](doc/devel/reports/reviews/sfs_eb_cohort_prior_2026-07-05.md)). **tomato2 validation (no truth):** M4 raises cohort het calls **3.15×** raw / **1.6×** at GQ≥20 vs G5. Investigation (report [sfs_prior_m4_cohort_validation_2026-07-05.md](doc/devel/reports/sfs_prior_m4_cohort_validation_2026-07-05.md)): **not** frequency over-estimation (called genotypes carry *less* alt than the reads: AC/AN 0.075 vs read-AF 0.094 at DP≤5); the low-GQ, M4-only tail is **correlated-mismapping** artifact (MQDiff −10 vs −3.4; a *MAPQ* class the coverage-paralog filter misses — `PARALOG_POST` 0.037) that GQ-filters out; the **confident** (GQ≥20) excess survives both GQ **and** MAPQ-diff filters → reads as strength-borrowing **recall**, not artifact. **Owner accepted M4.** Three refinements — F-aware first step, multi-start+evidence, read-based frequency — were **tried and reverted** (all solved the mis-diagnosed over-estimation). **Open:** confident-cohort recall is *plausible* but unproven (no cohort truth) — optional freebayes/GATK cross-check + a MAPQ-diff *filter* (annotated only today). **Next:** merge to `main`; G6 synthetic-shape end-to-end tests are covered by the M4 cohort tests.
> - **Last completed task (2026-07-05):** **SFS genotype prior — generalized to all `(ploidy, n_alleles)` via the Dirichlet-multinomial + old prior paths deleted (Milestones 1–2, GIAB-validated)**
>   (branch `sfs-genotype-prior`). The biallelic-only grid prior is replaced by the general **Dirichlet-multinomial** closed form, now the caller's **sole** genotype prior. **Milestone 1** (`src/genetics.rs`, types-first, each step multi-agent-reviewed): **G1** `lgamma` via a new `libm` dep; **G2** `dirichlet_multinomial_log_priors` (any shape, cross-checked vs an independent Pochhammer oracle); **G3** `alpha_from_diversity`. **Design decision (owner, "C + choice-2"):** `α_ref = 1`, `α_alt = θ̂/(n_alleles−1)` — the clean population-genetics marginals (`P(het)=θ`, `P(hom-alt)=θ/2`, monomorphic `1−3θ/2`, θ-independent 2:1); deliberately **not** the old grid's inflated hom-ref 0.878 (an un-normalised-sum artifact ~80× too high). **Milestone 2:** **G4** wired the DM into the engine for all shapes (posterior matches the pure primitive; F via the existing Wright mixture); **G4e** flipped the driver to feed θ̂ + per-sample F. **GIAB gate PASSED** — 5× SNP concordance **94.9 %** (grid 94.6 %, freebayes 93.4 %), matches/beats the grid at every depth, precision/recall/FP unchanged ([sfs_prior_dm_giab_validation_2026-07-04.md](doc/devel/reports/sfs_prior_dm_giab_validation_2026-07-04.md)). **G5** deleted the grid (`src/var_calling/sfs_prior.rs`, `e_step_sfs_biallelic`, `sfs_prior_tables`) **and** the HWE(p̂) plug-in (owner decision A — full delete; `nucleotide_diversity` now required, DM unconditional); `m_step_p̂` + pseudocounts kept for INFO/AF + QUAL. 1563 lib + all integration green. **Deferred (owner-approved, now Milestone 4):** the DM is `p̂`-independent so it drops the plug-in's **large-cohort strength-borrowing** (N weak-but-consistent hets → conservative hom-ref; site still emitted). The fix is the **empirical-Bayes** frequency-posterior update `α' = α_species + E[cohort counts]` (spec §4a); parked as an `#[ignore]`d integration test, gated on the tomato cohort. **Next:** Milestone 4 (empirical-Bayes large-cohort sharpening) + G6 remainder (tomato inbred sanity pass, synthetic multiallelic/polyploid end-to-end tests). Plan: [sfs_genotype_prior_generalization.md](doc/devel/implementation_plans/sfs_genotype_prior_generalization.md).
> - **Last completed task (2026-07-02):** **Hidden-paralog filter — pileup-window-coverage arch, Milestone M8: tomato2 validation of M7 (byte-identical, 3.97× faster)**
>   (branch `tomato2-paralog-filter`). Regenerated all 59 tomato2 `.psp` with the current binary (the on-disk ones predated M5/M6 —
>   no windowed columns, v2 summary) and ran the M7 vs M6 comparison in the container (8 cpu / 16 GB, ref bind-mounted). **Byte-identity
>   PASS:** M7-filter-on ≡ M6-filter-on except the `##commandline` binary-path line — all 262,539 records identical, provenance identical
>   (π=0.100521, lr_cut=3.8500, converged). **Drop-profile unchanged:** 20,809 dropped = 7.34 %. **Wall: 3.97× end-to-end** (M6 178.4 s →
>   M7 44.9 s; filter marginal cost +150.3 s → +16.8 s, ~8.9×; baseline caller 28.1 s). **RSS flat** (M7 350 MB vs M6 399 MB). **σ₀
>   recorded** (owed from M4): 59/59 fit, median 0.270 (0.208–0.329) via new `examples/tomato2_sigma0.rs`. Container-capped — prod 32-core
>   re-measure belongs on the Linux dev box, but the M6:M7 ratio is the trustworthy figure. The pileup-window-coverage arch (M1–M8) is
>   complete. Report: [paralog_m7_tomato2_validation_2026-07-02.md](doc/devel/reports/implementations/paralog_m7_tomato2_validation_2026-07-02.md).
> - **Last completed task (2026-07-02):** **Hidden-paralog filter — pileup-window-coverage arch, Milestone M7: score the LR once in the caller worker (single write pass)**
>   (branch `tomato2-paralog-filter`). The main perf payoff for the +11.6× two-pass regression. Placement A (the arch doc's later
>   refinement, now trivial because M5 made the worker gather window coverage): a cohort-constant `ParalogScoringContext` (pre-pass
>   coverage models + cohort `F` slice + σ₀ + precompute) is built once and held read-only by `VariantCaller`; `call_chunk` scores each
>   locus's LR inline (from the window coverage it gathers + the record's AD) onto `CalledChunk.paralog_lr` (NaN=unscored→kept). The sink
>   folds each stored LR into the calibration histogram inline (`SinkOutput::Spill(histogram)`), so calibration is now pure
>   (`calibrate_from_histogram`, no spill read) and the write pass reads the stored LR and applies the cut — **one spill read, one scoring,
>   on the parallel workers**. The old spill-reading `calibrate`/`score_spill_record` became `#[cfg(test)]`. Also replaced the per-locus
>   binary-search window gather with a per-sample forward merge-join (monotone cursor + debug_assert). VCF byte-identity holds by
>   construction (same pure scorer, same inputs, bit-exact spill round-trip; order-independent integer-bin histogram) — verified by the
>   review's reliability pass. 1519 lib + all integration + clippy/fmt green. 7-category parallel review: 2 Major (missing worker-scoring
>   test + stale module docs) + 3 Minor, all applied. Review: [paralog_m7_2026-07-02.md](doc/devel/reports/reviews/paralog_m7_2026-07-02.md).
>   Next in plan: M8 (validate on tomato2 — drop-profile/wall/RSS — + implementation report + real tomato2 σ₀).
> - **Last completed task (2026-07-02):** **Hidden-paralog filter — pileup-window-coverage arch, Milestone M6: retire the window-spill write side + `n_tiles`→`n_positions` rename**
>   (branch `tomato2-paralog-filter`). Per-sample centred-window coverage now flows `.psp` columns → `CalledChunk::window_coverage` →
>   the record spill (landed M5), so the old sibling "window spill" (Approach A / S6c: a tile-keyed coverage file joined per locus) is
>   fully dead. M6 deleted its whole write side — `WindowSpillBuilder`/`WindowSpillWriter`/`WindowSpillRecord`/`encode_window`,
>   `WindowMeanDepthAccumulator`, `ReferenceWindowGc`, the fold-loop wiring in `pipeline.rs`, the now-dead `ChunkPlan` accessors,
>   `TwoPhaseSegment::depth_at`, `ParalogPrePass::window_bp()` (3 modules removed, +105/−1302). `reference_base_matches` (the write-pass
>   coordinate guard) relocated into `write_pass.rs`. Renamed the coverage-by-GC histogram field `n_tiles`→`n_positions` (sliding-window
>   model = one sample per covered position) and **bumped `SAMPLE_SUMMARY_VERSION` 2→3** (review finding — the on-disk `.psp` key changed).
>   1517 lib tests + paralog integration green, clippy/fmt clean, VCF byte-identity preserved (record-spill codec + calibrate/write passes
>   untouched). 7-category parallel review: 1 Major + 1 Minor + 1 missing test, all applied.
>   Review: [paralog_m6_2026-07-02.md](doc/devel/reports/reviews/paralog_m6_2026-07-02.md).
>   Next in plan: M7 (score the paralog LR once in the caller worker → single write pass) then M8 (tomato2 validation + report).
> - **Last completed task (2026-07-02):** **Hidden-paralog filter — single-individual reformulation IMPLEMENTED + validated (Steps A/B/C)**
>   (branch `tomato2-paralog-filter`). Built the reformulation: **Step A** dropped `min_samples` (the
>   LR self-gates — under-powered → LR≈0 → kept); **Step B** made `F` the caller's single cohort
>   `--inbreeding-coefficient` and **deleted `Hexp` + the main-pass accumulator + `obs_het` +
>   callable-reference** (a correct per-individual `F` isn't AF-free-identifiable, and the per-variant
>   proxy is divergence-contaminated per Premise 3; `F` is a weak knob). Only global quantity left = π
>   (over loci). **Step C re-validation on tomato2 (no truth → profile coherence):** π = 0.092 (R1
>   ≈0.09, converged), **6.64 % dropped**, drop profile **preserved** (coverage 1.29× + het-excess
>   0.116 vs 0.046), **Jaccard 0.87** vs the old drop set, filter-off byte-identical. **n=1 works:**
>   one tomato sample drops 531 loci (7.1 %) that are het at **2.3× coverage** — graceful degradation
>   confirmed, coverage carries it. `F=0` default (outbred) still works because coverage is
>   load-bearing; selfing cohorts can raise `--inbreeding-coefficient`. Full suite green, clippy clean.
>   Report: [paralog_single_sample_validation_2026-07-02.md](doc/devel/reports/implementations/paralog_single_sample_validation_2026-07-02.md).
> - **Last completed task (2026-07-02):** **Hidden-paralog filter — single-individual reformulation DESIGNED (spec + plan; supersedes Premise 3)**
>   (branch `tomato2-paralog-filter`). Settled with the owner across a long design thread: make the
>   score degrade gracefully from **one individual** up, with no cohort-wide allele-frequency
>   quantity and no `min_samples` gate. Key resolution — the evidence separates by what it needs:
>   **per-individual** coverage model (load-bearing, n=1-safe); **per-locus het expectation inferred**
>   by H1's SFS-prior marginalisation (immune to the paralog self-masking because H1's coverage term
>   is `Normal(1)` regardless of frequency, so a 2× locus is unexplainable — the coverage is the
>   discriminator); **cohort-level, weak** `F` = the caller's single `--inbreeding-coefficient` (NOT a
>   per-individual `Hexp`-derived value — that reduces to the divergence-contaminated `n_het/n_variant`
>   ratio Premise 3 already measured as inverted); **across-loci** π unchanged. Three code changes
>   (drop `min_samples` → LR self-gates; `F` = cohort coeff; **remove `Hexp` + accumulator + obs_het**
>   → two-pass flow simplifies) + profile-based T1 re-validation (no truth set). Optional SFS-prior
>   calibration + low-n soft-flag deferred. Docs:
>   [hidden_paralog_single_sample_scoring.md](doc/devel/architecture/hidden_paralog_single_sample_scoring.md)
>   + [paralog_single_sample_scoring.md](doc/devel/implementation_plans/paralog_single_sample_scoring.md).
>   **Next:** implement Step A (drop gate) → Step B (F/Hexp) → Step C (re-validate).
> - **Last completed task (2026-07-02):** **Hidden-paralog filter — hard-error on missing `.psp` summary + fix perf-work example breakage**
>   (branch `tomato2-paralog-filter`). Closed the on-by-default footgun: `require_paralog_summaries`
>   now **fails loud** (`PipelineError::ParalogSummaryMissing`, names the first sample, points at
>   `--no-paralog-filter`) when the filter is on but a `.psp` carries no parseable coverage summary —
>   no more silent unfiltered callset on pre-summary `.psp`. A summary that parses but fails the
>   coverage fit (degenerate sample) is still carried absent, not fatal. Unit tests for the guard;
>   `thread_budget_integration` now passes `--no-paralog-filter` (synthetic summary-less `.psp`). Also
>   fixed the L1-memoization commit's missed propagation to `examples/paralog_score_parity.rs`
>   (`score_locus_for_paralogy` now takes `&ParalogScorePrecompute`) so `--all-targets` builds. Full
>   suite green, clippy clean. **Open (owner):** the T2 wall cost (perf review ranked H1 parallelise
>   /L1 done/S1 spill-LR) and the **single-individual graceful-degradation reformulation** (drop
>   `min_samples` + the pop-gen `Hexp`/`F` term → per-individual obs-het excess coupled to coverage;
>   π stays as the only global, across-loci quantity).
> - **Last completed task (2026-07-02):** **Hidden-paralog filter — freebayes comparison (GIAB safety + tomato2 FP-class)**
>   (branch `tomato2-paralog-filter`). **GIAB:** the cohort filter is a structural no-op on the
>   single-sample GIAB benchmark (0 loci scored → 0 dropped; on==off byte-identical across HG002/3/4
>   @50×) — no harm to the clean human callset; GIAB also excludes seg-dups so the paralog class is
>   absent by construction. **tomato2 vs freebayes (v1.3.10, same 59 samples / 160×200 kb): freebayes
>   emits 79 % of the loci our filter drops** (99 % ALT match; **55 % of all our drops pass freebayes'
>   default QUAL≥30**; fb QUAL mean 750 vs 585 kept — the depth-inflated confidently-wrong signature)
>   → the filter removes an FP class a standard caller leaves in. Figure
>   `benchmarks/tomato2/results/paralog_vs_freebayes.png` (local/regenerable). Report:
>   [paralog_vs_freebayes_2026-07-02.md](doc/devel/reports/implementations/paralog_vs_freebayes_2026-07-02.md).
>   (Aside: pre-existing single-pass stack overflow on HG003 @default stack — untouched code, needs `RUST_MIN_STACK`.)
> - **Last completed task (2026-07-02):** **Hidden-paralog filter — performance review of the T2 wall regression (+11.6×) → root-caused + solutions ranked**
>   (branch `tomato2-paralog-filter`). Ran the `rust-performance-review` skill (5 categories in parallel) over the two
>   scoring passes + the pure kernel. Root cause: the entire +249 s is `score_locus_for_paralogy` run
>   **single-threaded, twice** (calibrate + write) after the 16-way main pass finished — RSS flat, I/O only ~15 s.
>   Verdict *Apply the listed wins*: **H1** parallelise both passes (pure scorer, order-independent histogram fold,
>   write-pass re-serialise-before-emit → byte-identity holds); **L1** memoise the two locus-invariant log-prior
>   tables (H1 Wright 200×58 = ~26 B redundant `ln`/run, H2 carrier 40×58) as a bit-identical hoist; **L2**/**L3**
>   config + scratch hoists; **S1** score-once-spill-the-LR (owner sign-off). Gate: add a paralog criterion bench +
>   one flamegraph. Build profile already optimal. Report: [perf_paralog-filter_2026-07-02.md](doc/devel/reports/reviews/perf_paralog-filter_2026-07-02.md).
> - **Last completed task (2026-07-02):** **Hidden-paralog filter — Milestone T (T1 behaviour + T2 cost gate) on tomato2 → T1 PASS, T2 memory-flat but wall +11.6× (owner decision pending)**
>   (branch `tomato2-paralog-filter`). Ran the production two-pass filter on the 59-sample tomato2
>   cohort (regenerated `.psp` with the stored summary section — the committed fixtures predate it,
>   which silently no-op'd the filter: **footgun** — recommend a loud warn/error when few samples
>   fit). **T1 PASS:** 58/59 fit, **π = 0.099** (R1 ≈0.09, EM converged), 6.98 % dropped; dropped set =
>   **coverage excess (1.27×) + het excess (0.124 vs 0.045)** = the R1 paralog signature; AF is not the
>   discriminator; introgression-like (normal-coverage) loci survive. **T2:** peak RSS **flat**
>   (336→363 MB, +8 %, independent of 376 k loci — the memory obligation holds) but **wall +11.6×**
>   (23→272 s) from the per-locus H1/H2 scorer running over 376 k loci × 58 samples **twice**,
>   **single-threaded** (calibrate + write passes); spill 3.5 GB on disk (record spill, verbatim
>   `PosteriorRecord`). **Owner decision (T2 gate):** parallelise the two scoring passes / spill the
>   LR once / or ship off-by-default — memory is fine, cost is an optimisation problem not a design
>   flaw. Report: [paralog_t1_t2_2026-07-02.md](doc/devel/reports/implementations/paralog_t1_t2_2026-07-02.md).
> - **Last completed task (2026-07-02):** **Hidden-paralog filter — Milestone S COMPLETE: var-calling wiring (S6c part 2) → the two-pass filter runs end-to-end and is ON by default (S6 checkpoint)**
>   (branch `tomato2-paralog-filter`). Approach A landed: the producer writes a second ephemeral
>   **window spill** ((chrom_id, tile, gc, [per-sample depths]) per variant window) alongside the
>   record spill; the calibrate + write passes join the two by tile key in coordinate lockstep
>   (`WindowJoin`), so the byte-identity-critical caller/chunk types stay untouched. Steps: ChunkPlan
>   + `TwoPhaseSegment::depth_at` accessors; `WindowSpillBuilder` fed inline in both fold loops with
>   a per-contig `ReferenceWindowGc` (separate streaming fetcher); `score_spilled_locus` takes the
>   joined (gc, depths); write-pass reference-base consistency guard (fail loud on coordinate drift);
>   `PARALOG_WINDOWS_WIRED` deleted → filter **on by default** at `--paralog-fdr 0.01`
>   (`--no-paralog-filter` restores the single-pass path). Adversarial review fix: `flush_ready`
>   flushes by the cohort-fed frontier (force-closing lagging samples' open windows) so `pending`
>   stays bounded near the frontier, not O(contig tiles × samples). New end-to-end two-pass
>   integration test (provenance header + spill cleanup). **1512 lib + all integration tests green,
>   clippy clean, byte-identity preserved (existing var-calling tests pass `--no-paralog-filter`).**
>   Owed: **T1** (real tomato2 drop-profile vs R1) + **T2** (cost/RSS gate — the on-by-default cost);
>   the E2E test asserts orchestration/cleanup but defers the *synthetic drop* assertion to T1
>   (`min_samples`=20 + coverage-fit make a deterministic synthetic drop brittle; drop mechanics are
>   unit-tested in `write_pass`). Plan: [paralog_varcalling_wiring.md](doc/devel/implementation_plans/paralog_varcalling_wiring.md).
> - **Last completed task (2026-07-01):** **Hidden-paralog filter — R1 data-first validation on tomato2 (checkpoint) → the pure maths is trustworthy on real data**
>   (branch `tomato2-paralog-filter`). New `examples/paralog_score_parity.rs` (Rust harness: reads
>   59 `.psp` + cohort VCF, runs Q2–Q5 end-to-end, dumps per-locus LR/post/qval + per-sample σ₀/F)
>   + `benchmarks/tomato2/src/paralog_score_parity.py` (uv companion: the checks Rust can't — LR
>   correlation vs the prototype parquet, F Spearman, flagged-set profile). **All data-first checks
>   met:** σ₀ 0.282 (≈0.28); mode/median guard passes 58/59 real samples; **LR Pearson 0.9931** vs
>   the prototype (≥0.98 target) on 268 537 shared loci; **π 8.88 %** (≈9 %); **F Spearman 0.851**
>   (≈0.86); flagged set = paralog profile (obs_het 0.115 + 1.31× coverage vs kept 0.043 + 1.04×).
>   The data caught one real bug (the `Hexp` per-position scale — must be Σ2pq/callable not
>   mean-2pq-over-variant-sites; the S1 wiring must match). Report:
>   [paralog_r1_data_validation_2026-07-01.md](doc/devel/reports/implementations/paralog_r1_data_validation_2026-07-01.md).
> - **Last completed task (2026-07-01):** **Hidden-paralog filter — Milestone Q COMPLETE (Q4 inbreeding + Q5 empirical-Bayes prior/FDR) → the pure statistics core is done (checkpoint)**
>   (branch `tomato2-paralog-filter`). Q4: `obs_het` (het rate) + `inbreeding_coefficient`. Q5:
>   bounded-RAM `ParalogLrHistogram` + `ParalogPrior::estimate` (EM for π over histogram bins) +
>   `ParalogFdrCurve` (tail-FDR `q_of_lr` + `lr_threshold_for_fdr`). Numerically verified faithful
>   to `build_paralog_eb.py` (histogram EM ≡ full-vector mean; tail-FDR-from-top ≡ rank-by-posterior
>   running mean). Fixed constructor-contract inconsistency (→ `Option`), silent EM non-convergence
>   (→ `converged` flag), duplicated binning (→ shared `UniformBinning`). **65 paralog tests**, 1463
>   lib tests, fmt/clippy clean. **Checkpoint reached: the three statistics pieces (coverage model,
>   scorer, prior/FDR) are built and unit-tested in isolation.** Next: R1 (data-first validation on
>   tomato2), then Milestone S (var-calling wiring). Report: [paralog_q5_prior_2026-07-01.md](doc/devel/reports/reviews/paralog_q5_prior_2026-07-01.md).
> - **Last completed task (2026-07-01):** **Hidden-paralog filter — Milestone Q, step Q3 (`score_locus_for_paralogy` — the H1-vs-H2 marginal likelihood ratio)**
>   (branch `tomato2-paralog-filter`). Pure per-locus scorer: H1 (real variant) marginalises allele
>   freq under the folded-SFS prior with Wright HWE genotype priors + genotype-independent Normal
>   coverage; H2 (hidden paralog) marginalises carrier config × carrier-freq, coverage `Normal(T/2,
>   σ₀√(T/2))`, `vaf=m/T<1` (hom-alt veto falls out). Result = `logL(H2)−logL(H1)`, a pure LR.
>   Two parallel review sub-agents **verified numerical faithfulness** to `build_paralog_lr.py`
>   line-by-line, then caught + fixed 3 release-mode safety holes (u32 underflow, zero-σ₀ NaN,
>   length-mismatch truncation). **44 paralog tests** incl. an absolute-value hand-computed parity
>   anchor; fmt/clippy clean. Report: [paralog_q3_locus_score_2026-07-01.md](doc/devel/reports/reviews/paralog_q3_locus_score_2026-07-01.md).
> - **Last completed task (2026-07-01):** **Hidden-paralog filter — Milestone Q, steps Q1–Q2 (pure statistics module: `ParalogModelParams` + `SingleCopyCoverageModel::fit`)**
>   (branch `tomato2-paralog-filter`). New `src/paralog/` peer module (arch Premise 0). **Q1:**
>   `ParalogModelParams` + `SfsPriorSpec`/`GridSpec` with prototype-faithful defaults. **Q2:**
>   `SingleCopyCoverageModel::fit` from the binned coverage-by-GC histogram — mode-anchored
>   `single_copy_scale` (parabolic sub-bin refine + mode/median sanity guard), per-GC weighted-median
>   `gc_bias_curve` (gap-filled + smoothed), σ₀ = 1.4826·MAD floored at histogram resolution;
>   `relative_copy_number(gc,depth)`. Two parallel review sub-agents caught + fixed a σ₀
>   band-collapse over-confidence bug, a dead error variant (→ reachable `DepthModeAtBottomBin`),
>   and a tuple→`ModeMedianRatioBounds` newtype; **24 paralog tests**, 1422 lib tests, fmt/clippy
>   clean. Reports: [paralog_q2_coverage_model_2026-07-01.md](doc/devel/reports/reviews/paralog_q2_coverage_model_2026-07-01.md).

---

## Pipeline stages

Stage descriptions are one-line reminders; the spec is authoritative.

### Stage 1 — per-sample pileup (BAM → `.psp`)

Stage 1 reads each BAM/CRAM once per sample and writes one `.psp` artefact.

#### Analysis regions from a BED file (`--regions`)
- **Status:** fixes-applied (2026-06-16, on branch `bed-regions-review`; re-review + fixes after the segment-reader retrofit); feature merged to `main` (`6ecd22c`), fix branch not yet merged
- **Plan:** [bed_regions.md](doc/devel/implementation_plans/bed_regions.md)
- **Impl report:** [bed_regions_2026-06-03.md](doc/devel/reports/implementations/bed_regions_2026-06-03.md)
- **Perf + byte-identity investigation:** [bed_regions_perf_and_byte_identity_2026-06-03.md](doc/devel/reports/bed_regions_perf_and_byte_identity_2026-06-03.md)
- **Cross-cutting:** Stage 1 (pileup) + var-calling. The pipeline always operates on a `RegionSet` — the `--regions` BED, or one full-length span per contig (whole genome). The whole-file pileup streaming path is retired; pileup now always seeks via the index.
- **Code:** [src/regions.rs](src/regions.rs) (`RegionSet` + BED parser); region-driven `run_pileup` ([src/pop_var_caller/cli.rs](src/pop_var_caller/cli.rs)) now opens pooled per-file readers once and k-way-merges per region via `SegmentMergedReads` ([src/bam/segment_merge.rs](src/bam/segment_merge.rs) + [src/bam/segment_reader.rs](src/bam/segment_reader.rs)) — the old per-region `AlignmentMergedReader::query` is **gone** (the segment-read-fetcher retrofit landed; cf. cli.rs:374-375); var-calling restriction in [src/var_calling/pipeline.rs](src/var_calling/pipeline.rs) (`restrict_intervals_to_regions`); command-line + regions provenance in [src/psp/header.rs](src/psp/header.rs).
- **Latest code review:** [regions_2026-06-16.md](doc/devel/reports/reviews/regions_2026-06-16.md) — **Approve-with-changes** (re-review of the current segment-reader path; 10 categories): 0 Blockers, 6 Major, 12 Minor + Nits. Top: M4 worker `JoinHandle`s dropped (a worker panic could silently truncate the read stream — only `thread::scope`'s re-raise saves it); M1 silent empty output on a non-selecting BED; M3 + M2 = prior M1/M3 still unfixed. Both this branch's cleanup commits verified behavior-identical. Audit trail `tmp/review_2026-06-16_regions/`.
- **Latest fixes-applied:** [fixes_applied_2026-06-16_v2.md](doc/devel/reports/reviews/fixes_applied_2026-06-16_v2.md) — 17 Applied + 2 Applied-with-adaptation, 4 Deferred, 1 Disputed (Mi11). All gates green (1165 lib + integration tests, clippy `-D warnings`, doc). Notably resolves the two carried-over 2026-06-03 Majors (counter-fold field-safety, `BedError::Io` context). Q1=empty-BED→hard error (`BedError::NoRegions`); Q2=`MappedReadSource` kept as a documented DI seam.
- **Prior review:** [bed_regions_2026-06-03.md](doc/devel/reports/reviews/bed_regions_2026-06-03.md) — Approve-with-changes on the now-superseded `query` impl; 0 Blockers, 6 Major, 9 Minor. HG002-bottle accuracy identical to `main` (SNP F1 0.9037), peak RSS −71%.
- **Latest perf review:** [perf_bed-regions_2026-06-16.md](doc/devel/reports/reviews/perf_bed-regions_2026-06-16.md) — Profile-first → **measured** on tomato1 (§8): a fragmentation sweep (80→8000 regions over identical 8 Mb = **+1.9% wall**) + a sampling profile (`probaln_glocal` ~70%; whole `--regions` read/merge path ~1–2%) **refute the per-region/read-path findings as a perf lever.** L1 (qname double-clone) + L3 (dead `Stage1Outputs` clones) applied as **perf-neutral cleanups**; rest closed/deferred. Real lever is BAQ (out of scope). Audit trail `tmp/perf_review_2026-06-16_bed-regions/`.
- **Open (from the 2026-06-03 review):**
  - **M1** — the three `merge()` counter-folds (`RunSummary`/`FilterCounts`/`BaqSkipCounts`) aren't field-addition-safe (exhaustive-destructure + 3 tests).
  - **M2** — `ContigInterval { start: 0, .. }` is constructible and aborts `query_interval`'s `Position::new(..).expect()` (checked constructor / drop `pub` fields).
  - **M3** — `BedError::Io` drops path/line and flattens its `io::Error` into `Display` (no `#[source]`).
  - **M4** — FASTA-repository construction duplicated ×3; cross-file validation copied between `new` and `load_pileup_inputs`.
  - **M5** — `--build-map-file-index` (index management) folded into this PR (split or document the coupling).
  - **M6** — whole-contig range collapse applied silently with no runtime trace.
  - **Minors:** per-region repeated work (FASTA repo rebuild + clones, undercuts "sparse BED is cheap"); pileup-vs-var-calling provenance asymmetry; empty-BED → silent empty output; `build_map_file_index` third-name / `stashed_upstream` noun; `#[from] AlignmentIndexError` collapses two origins; `#[allow(too_many_arguments)]` without justification. Plus the §8 missing tests (merge folds, `command-line` serde-default round-trip, adversarial BED coordinates, empty-BED pin).
  - **Efficiency follow-up:** indexed `.fai`-seek fetcher to replace the per-region streaming `MultiChromStreamingRefFetcher` rebuild (many sub-contig regions on a large contig).
- **Code-review fixes (2026-06-16, see [fixes_applied_2026-06-16_v2.md](doc/devel/reports/reviews/fixes_applied_2026-06-16_v2.md)):**
  - **Fixed:** M1 (empty BED → `BedError::NoRegions`), M2 + M3 (the carried-over 2026-06-03 Majors: `BedError::Io` context + counter-fold field-safety), M4 (explicit worker join + doc), M6 (cross-thread determinism test), and Mi1–Mi10 + nits. Mi11 disputed (kept as a documented DI seam).
  - **Still open — M5** (deferred): staged/inline producer interval-walk duplication (`pipeline.rs:431-595`) — regression-prone; extract a shared `drive_schedule` helper as its own byte-identity-gated effort.
  - **Still open — Mi12** (deferred): wire `cargo audit` + `deny.toml` (cargo-audit not installed; no new deps this branch).
  - **Still open — partial Mi4:** the per-line overhang-*clamp* warning (the mode-announce half landed); minor follow-ups Nit-bedlines / Nit-cache.
- **Resolved by the 2026-06-16 perf review (measured — see report §8):**
  - **Applied (perf-neutral cleanups, gates green):** L1 qname double-clone → move ([segment_merge.rs](src/bam/segment_merge.rs)); L3 dropped the unread `Stage1Outputs` `sample_name`/`contigs` fields ([stage1_pipeline.rs](src/pop_var_caller/stage1_pipeline.rs)). Neither moves wall (~1% paths) — kept as cleanups, not speedups.
  - **Closed — no gain:** L5/L6 (merge path ~1–2% of runtime); L7 region-parallelism (serial loop ~2% even at 8000 regions; cores already saturated).
  - **Deferred — untestable on pre-sliced CRAMs:** L2 `.crai` head binary-search, L4 container cache — need a full un-sliced CRAM + fragmented-BED workload (tomato1 CRAMs are pre-sliced; tiny `.crai`). Still mirrored by the segment-read-fetcher block's deferred items.
  - **Real lever (out of this feature's scope):** BAQ (`probaln_glocal` ~70% of pileup self-time) — see [perf_baq_2026-05-12.md](ia/reviews/perf_baq_2026-05-12.md).

#### Alignment-file input (CRAM + BAM)
- **Status:** fixes-applied (BAM slice); shipped (CRAM)
- **Plans:**
  - CRAM slice: [per_sample_caller_cram_input.md](doc/devel/implementation_plans/per_sample_caller_cram_input.md)
  - BAM slice: [bam_input_support.md](doc/devel/implementation_plans/bam_input_support.md)
- **Impl reports:**
  - CRAM slice: [per_sample_caller_cram_input_2026-04-29.md](doc/devel/reports/implementations/per_sample_caller_cram_input_2026-04-29.md)
  - BAM slice (2026-05-24): [bam_input_support_2026-05-24.md](doc/devel/reports/implementations/bam_input_support_2026-05-24.md)
- **Code:** [src/bam/](src/bam/) — `alignment_input.rs` (merge + filter + header validators, format-agnostic), `cram_input.rs` + `bam_input.rs` (per-format owned record-stream decoders), `index_preflight.rs` (CRAI / CSI / BAI detection + build), `errors.rs`.
- **Latest reviews:**
  - BAM slice (2026-05-24): [bam_input_support_2026-05-24.md](doc/devel/reports/reviews/bam_input_support_2026-05-24.md) — Approve-with-changes: 0 Blockers, 19 Major (M1–M19), 22 Minor (Mi1–Mi22), grouped Nits. Architecture sound (no `unsafe`, no shared mutable state, `'static + Send` iterators correct by construction); Major findings cluster around (a) rename-debt — CRAM vocabulary still in error displays + variant name `PileupCliError::CramInput` mislabels every BAM error; (b) missing `#[non_exhaustive]` on `AlignmentInputError` / `PileupCliError` / `BamIndex`; (c) `load_per_input_headers` reimplements opener work AND skips the CRAM-version gate that the per-format helper has — CRAM 4.x silently passes header-load; (d) `.csi`/`.bai` policy + CSI depth defaults live in function bodies, not named constants; (e) reliability test gaps for the new `AlignmentIndexFormatMismatch`, `UnsupportedExtension`, indexed-BAM chunk-walk, and `load_alignment_index` branches. Six missing-test specs in §8 of the report.
- **Latest fixes-applied (BAM slice):** [fixes_applied_2026-05-24.md](doc/devel/reports/reviews/fixes_applied_2026-05-24.md) — **Completed**: every Major (M1–M19) Applied; 13 of 22 Minors Applied; 7 Minors Deferred (Mi5 / Mi7 / Mi12 / Mi13 / Mi15 / Mi20 / Mi21) + M17's redesign half. 10 fix commits `9fc1df0` → `7a8569b`. **Total +10 new tests** (1 cram_input regression + 4 index_preflight `load_alignment_index` triplet + 1 alignment_input AlignmentIndexFormatMismatch + 3 bam_input [CSI arm, multi-chunk, truncated-BAM err] + 1 pileup_cli UnsupportedExtension + 1 cohort_cli MixedFormat lock-step). 924 lib + 45 integration tests pass; cargo fmt / clippy / test all clean.
- **Latest reviews / fixes (CRAM):** [per_sample_caller_cram_input_2026-04-29.md](doc/devel/reports/reviews/per_sample_caller_cram_input_2026-04-29.md), [fixes_applied_2026-05-01.md](doc/devel/reports/implementations/fixes_applied_2026-05-01.md)
- **Open (BAM slice — after the 2026-05-24 fix run):**
  - Wall-time validation on real multi-chrom BAMs (analogue of cohort H1's 3.85× at T=13 on tomato CRAMs). Picked up alongside the `examples/profile_from_bam_e2e.rs` + `benches/from_bam_e2e_perf.rs` infrastructure already deferred from the prior task.
  - **Deferred Minors** (from the fix run): **Mi5** (`process_one_chromosome_from_bam` rename — subcommand-name coupling), **Mi7** (`input_crams` PSP field rename — reach beyond `src/bam/`), **Mi12 + M17 redesign half** (lift `MixedAlignmentFileFormats` into a shared sub-enum), **Mi13** (`From<AlignmentIndexError>` 4-arm passthrough — design call), **Mi15** (`crate::bam` → `crate::alignment` module rename), **Mi20** (`OwnedIndexedBamRecords::next` phase extraction — cosmetic), **Mi21** (per-record `RecordBuf::default()` allocation — folds into parallelisation-tuning workstream), **Mi9 partial** (`tests/common::build_csi` copy stays until a test-support feature flag is introduced).
  - Lift the no-mixing restriction (CRAM + BAM in one invocation) if a real workload appears. Pre-flight gate is the only place the restriction lives; merge core is already format-agnostic.

#### Segment read fetcher (shared indexed-segment read source)
- **Status:** fixes-applied (2026-06-16; increments #1–#3; standalone primitive, no consumer wired yet)
- **Plan:** [segment_read_fetcher.md](doc/devel/implementation_plans/segment_read_fetcher.md)
- **Impl report:** [segment_read_fetcher_2026-06-16.md](doc/devel/reports/implementations/segment_read_fetcher_2026-06-16.md)
- **Code:** [src/bam/segment_reader.rs](src/bam/segment_reader.rs) — `AlignmentFile` (`BamFile`/`CramFile`), `from_input`, `get_reads_from_segment`, pooled re-seekable readers (`Mutex<Vec<Handle>>`), `MappedReadsInSegment` per-call iterators, narrow `SegmentReadFilter` config, shared `classify_segment_record` per-record cascade. Reuses `classify_pre_decode` / `record_buf_to_mapped_read` / `query_interval` / `ContigInterval::overlaps_record` (visibility lifted to `pub(super)`); two new typed errors (`MissingCramReference`, `InvalidSegment`).
- **Tests:** 21 in-module (overlap + 1-based-inclusive edges, read spanning two segments yielded by both, pool open-once / resting size, poisoned-pool recovery, `par_iter` determinism on a shared `&AlignmentFile`, empty segment, CRAM container early-stop, target-contig + MAPQ filters, both format/index mismatch arms, error paths; BAM + CRAM, `.crai` via `cram::fs::index`). All gates green (fmt / clippy `-D warnings` / doc; 1059 lib + integration tests).
- **Latest review:** [segment_reader_2026-06-16.md](doc/devel/reports/reviews/segment_reader_2026-06-16.md) — **Approve-with-changes**: 0 Blockers, 6 Major, 9 Minor + Nits. Concurrency design verified sound; happy path well tested both formats. Audit trail: `tmp/review_2026-06-16_segment-reader/`.
- **Latest fixes-applied:** [fixes_applied_2026-06-16.md](doc/devel/reports/reviews/fixes_applied_2026-06-16.md) — all 6 Majors **Applied** (M1 Mutex-poison recovery + test; M2 explicit format/index mismatch arms + test; M3 shared `classify_segment_record` chokepoint; M4 CRAM container early-stop + test; M5 narrow `SegmentReadFilter`; M6 report-location convention fixed in the skills on `main` `8c32b00` + report moved) plus Mi1 (handle-invariant docs) and Mi2 (`#[cfg_attr(not(test), allow(dead_code))]`). 1059 lib + integration tests pass; the only `--all-targets` failure is the pre-existing `psp_writer_perf` bench panic.
- **Open (deferred from the 2026-06-16 fix run):**
  - **Minors:** `get_reads_from_segment` `get_`-prefix / diverges from `query` (Mi3, naming — decide at #5); pool-never-shrinks doc note (Mi4); residual inline `Io` sites (Mi5, reader-scoped borrow constraint); cast-spelling + boundary test (Mi6); CRAM container-skip `span==0` test (Mi7, blocked on a multi-container fixture knob); `from_input` path↔header validation (Mi8, when #4 lands); test-only back-ref into `pileup::per_sample::cram_files` (Mi9). Plus the §8 missing tests not yet added (handle-leak-on-index-error, BAM error-fuse, min-length-filters-everything ×2, CRAM error paths, downstream-container completeness, `resolve_segment` boundary) — several blocked on the multi-container fixture limitation. Nits grouped for a cosmetic pass.
  - **#4** — wire the SSR Stage-1 `fetch_locus_reads` onto it (deferred: SSR module not on this branch).
  - **#5** — retrofit the SNP `--regions` path (replace per-BED-interval re-`query`), measuring the `--regions` tax improvement; converge the duplicated chunk/container-walk with `OwnedIndexed{Bam,Cram}Records` then. Remove the module-level dead-code allow when a consumer lands.
  - Deferred (plan §8): CRAM container cache; binary-search the crai head (measure-first, when #4 provides a real workload).

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

#### Indel normalization (left-alignment)
- **Status:** implemented + reviewed + restructured (branch `indel-normalization`, not yet merged to `main`)
- **Spec:** [calling_pipeline_architecture.md §"Indel normalization (left-alignment)"](doc/devel/specs/calling_pipeline_architecture.md) (commit `8d7dd17`)
- **Plan:** [indel_normalization.md](doc/devel/implementation_plans/indel_normalization.md)
- **Final approach (revised — option 2):** the code review surfaced that the
  BAM/CRAM input cascade **already** left-aligned every read's indels (the
  always-on "F3" pass `left_align_indels` in
  [src/bam/alignment_input.rs](src/bam/alignment_input.rs), pre-existing on
  `main`). Rather than add a second normalization in the BAQ prep stage, the
  feature now **replaces F3's single-forward-pass shifter with the GATK
  `leftAlignIndels` port** at F3's call site: F3's `left_align_indels` /
  `try_apply_indel_shift` (+ helpers, + its `f3_*` tests) are deleted, and
  the reader calls [`indel_norm::left_align_indels`](src/pileup/walker/indel_norm.rs)
  (GATK port: right-to-left, trim-first, dual ref+read check, **collision-merge**
  — the refinement F3 explicitly punted on). The reader's existing reference
  fetch (for the F1 mismatch filter) is reused; nothing else in the merged
  reader changes. Because the reader is upstream of everything, normalization
  is structurally mandatory — including under `--no-baq` (no BaqEngine
  wiring; the earlier prep-stage detour was reverted).
- **Code:** [src/pileup/walker/indel_norm.rs](src/pileup/walker/indel_norm.rs) (the GATK port + tests); call site at [src/bam/alignment_input.rs](src/bam/alignment_input.rs) ("F3").
- **Latest review:** [indel_normalization_f3_replacement_2026-05-29.md](doc/devel/reports/reviews/indel_normalization_f3_replacement_2026-05-29.md) — Approve-with-changes, of the F3-replacement diff. Confirmed the call-site integration is correct (`ref_seq[0]` is the read's first aligned base → `read_start=0`; reference fetch reused from the F1 filter; F3 fully deleted, port a faithful drop-in). Fixed during the review: a stale F3 comment falsely claiming "adjacent indels invariant" (the port merges them), a broken `[left_align_prepared]` doc link (a `cargo doc` regression this PR introduced), the missing `// UNREACHABLE:` comment, and the untested `left_align_indels` wrapper/debug-invariant (added `wrapper_*` tests). The [first review](doc/devel/reports/reviews/indel_normalization_2026-05-29.md) (BaqEngine detour) is superseded.
- **Open (from the F3-replacement review):**
  - **M2** — add a reader-level integration test (`AlignmentMergedReader` → `fetch_ref_for_read` → `left_align_indels`); the deleted `f3_*` unit tests have no integration equivalent, so a slice-offset/`read_start` wiring regression would be invisible.
  - **M4** — hot-path allocation regression: F3 shifted in place (zero-alloc); the port allocates 4–5 `Vec`s per indel read on the single-threaded reader loop. Add scratch reuse + a normalization criterion bench (no bench covers this path).
  - **M5** — move `indel_norm` to `src/bam/` (beside its sole consumer; removes a `bam → pileup` edge) — or a neutral peer if a second consumer appears.
  - **M6** (needs verification) — the port can merge colliding indels into an adjacent `D`/`I` pair that G2 `cigar_is_bad` (run earlier) doesn't re-check; confirm reachability + walker handling with a test.
  - **Minors:** read-consume guard emits no telemetry (Mi1); demote `left_align_cigar`/`LeftAlignResult` to private (Mi2); `indel_norm/tests.rs` dir convention (Mi3); plus carryovers `Range`→`IndexRange`, predicate-helper factoring, `build_cigar` splice→fold.
  - Deferred measurement: HG002/tomato indel-recall delta (the original acceptance signal).

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
- **Impl reports:**
  - Cohort slice:
    [pop_var_caller_cohort_cli_2026-05-19.md](doc/devel/reports/implementations/pop_var_caller_cohort_cli_2026-05-19.md)
  - `var-calling-from-bam` per-chromosome parallelism (2026-05-24):
    [var_calling_from_bam_per_chromosome_2026-05-24.md](doc/devel/reports/implementations/var_calling_from_bam_per_chromosome_2026-05-24.md)
- **Latest reviews (cohort slice):**
  [ref_fetcher_2026-05-23.md](doc/devel/reports/reviews/ref_fetcher_2026-05-23.md)
  (code review of `src/per_sample_pileup/ref_fetcher.rs` — status
  *Request-changes*; 3 Blockers, 23 Major, 19 Minor. Production
  code paths are correct; gaps are in the public-trait contract
  surface. Highlights: B1 `.fai` `line_bases=0` panic on
  attacker-influenced input; B2 `legacy_io_error` drops
  `contig_name` from every legacy-walker diagnostic; B3 the
  `OutOfPattern` contract has no production-buffer-size regression
  test; convergent Majors: stale doc link at line 282 — 6
  categories, missing `#[non_exhaustive]` on `ChromRefFetchError`,
  `Source::Memory` cfg-gating, vestigial legacy surface on
  `StreamingChromRefFetcher` after the Step-2 migration, ~80 LOC
  constructor duplication across three sites.);
  [perf_ref_fetcher_2026-05-23.md](doc/devel/reports/reviews/perf_ref_fetcher_2026-05-23.md)
  (FASTA fetcher post-Step-2-migration — verdict: *Apply the listed
  wins*, gated on H4; 4 Hot-path, 10 Likely, 4 Speculative;
  headline: `ChromRefBaseIter::next` 10.96% wall with 4.02% on
  `MutexGuard<StreamState>` drop alone — H1 drops the `Sync`
  requirement, swaps `Mutex<StreamState>` for `RefCell<StreamState>`);
  [perf_psp_to_vcf_2026-05-20.md](doc/devel/reports/reviews/perf_psp_to_vcf_2026-05-20.md)
  (end-to-end perf review against real tomato data — verdict: *Apply
  the listed wins*; 7 Hot-path findings, 12 Likely, 5 Speculative;
  headline: pipeline is essentially single-threaded, DUST 33 %,
  allocations 21 %; H1 per-chromosome parallelism is the
  order-of-magnitude lever);
  [cohort_cli_2026-05-19.md](doc/devel/reports/reviews/cohort_cli_2026-05-19.md)
  (correctness — Request-changes: 0 Blockers, 14 Major (M1–M14), 23
  Minor + grouped Nits).
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
- **Open (from 2026-05-23 ref-fetcher perf review):**
  - **H4 (gate)** — Repair the four broken benches before any other
    finding can be measured: `cohort_e2e_perf`, `var_calling_perf`,
    `baq_perf`, `pileup_walker_scaling`. API drift from the Step-2
    `ChromRefFetcher` migration (stale `SyncRefFetcher` import,
    `CohortPipelineParams.fetcher`/`.chromosomes` field drift,
    `drive_cohort_pipeline` arg count). Mechanical port; no new deps.
  - **H1** — `src/per_sample_pileup/ref_fetcher.rs:117` —
    `Mutex<StreamState>` → `RefCell<StreamState>` and relax
    `SharedRefFetcher = Arc<dyn ChromRefFetcher + Send + Sync>` to
    `Arc<dyn ChromRefFetcher + Send>`. Per-base mutex unlock is
    4.02 % of cohort wall time; the `Mutex` is uncontended by
    construction (per-worker ownership documented at L113-117).
  - **H2** — `iter_bases` returns `Box<dyn Iterator>`; GAT-monomorphic
    `type BaseIter<'a>` (or a `ChromRefFetcherTyped` subtrait that
    keeps the dyn-safe surface intact) removes per-byte vtable
    dispatch on the DUST mask scan.
  - **H3** — `StreamingChromRefFetcher::fetch` returns owned `Vec<u8>`
    where every PerGroupMerger consumer only borrows. Switch to
    `&[u8]` — entangled with H1's `RefCell` (lands together).
  - **L1 / L10** — Wrap `Source::File(File)` in `BufReader<File>` in
    both `StreamingChromRefFetcher` and `ManualEvictChromRefFetcher`.
  - **L2** — Autovec uppercase pass in `refill` /
    `read_uppercased_bases` — strip-newlines via `filter+take+extend`
    then `make_ascii_uppercase` on the destination slice.
  - **L3** — `src/var_calling/dust_filter.rs:728` —
    `io::Error::other(format!("{e}"))` → `io::Error::other(e)`
    (drops the `String` alloc, preserves `source()`).
  - **L5** — `ChromRefBaseIter::Drop` re-locks the mutex; folds
    into H1 for free.
  - **L6** — Add `[profile.release-with-debug]` to `Cargo.toml`
    that `inherits = "release"` + `debug = true`, so future
    profile captures reproduce the inlined-frame call-graph from
    a clean checkout.
  - **L7** — Wire `mimalloc::MiMalloc` as `#[global_allocator]`
    behind `alloc-mimalloc` in `src/main.rs`. Glibc allocator
    symbols sum to ~14.5 % cpu_atom / ~23.2 % cpu_core; gate
    merge on a paired A/B against the 13-thread server.
  - **L8** — Add `benches/ref_fetcher_perf.rs` with four
    `criterion` functions (`streaming_iter_bases_full_contig`,
    `streaming_fetch_per_group_window`, `manual_evict_fetch_then_evict`,
    `refill_warm_cache`) so H1/H2/H3 effects can be localised below
    the cohort_e2e noise floor.
  - **L4 / L9 / S1–S4** — see the full report.
- **Open (from 2026-05-20 perf review):**
  - **H1 — closed 2026-05-20** (per-chromosome parallelism shipped).
    Five-commit PR `309a5be` → `0b1e958`; impl report
    [cohort_per_chromosome_parallel_2026-05-20.md](doc/devel/reports/implementations/cohort_per_chromosome_parallel_2026-05-20.md).
    Realised **3.85× wall-time speedup at T=13** on the multi-chrom
    real-data fixture (106.6 s → 27.7 s). Workload imbalance
    (ch00 carries 1.4 M of 2.6 M total reads as the unplaced/decoy
    contig) gates the ceiling below the plan's predicted 6–10×
    range; L5 contention absorbs another ~23 %.
  - **L1 — closed 2026-05-20** (per-group merger inner `par_iter`
    removed, commit `309a5be`). Mandatory prep for H1 — nested
    rayon under per-chrom outer is wasteful and would have
    polluted the H1 measurement.
  - **H2 / H3** — DUST `find_perfect` inner-loop: replace
    `Vec::insert(j, …)` (O(n) memmove) and `VecDeque` indexing
    (per-element wrap + bounds check). Gated by **H7**.
  - **H4** — PSP reader: replace per-block `SeekFrom::Start` (which
    discards the 64 KiB BufReader) with `seek_relative` for the
    post-header rewind and "skip if already at offset" for the
    pre-block seek (or pull-style `fill_buf`/`consume`).
  - **H5** — M5 verify: stream the FASTA in 64 KiB windows
    feeding `Md5::update` instead of materialising whole contigs (up
    to 91 MB) and uppercasing byte-by-byte.
  - **H6** — `PerPositionMerger::next`'s `vec![None; n_samples]` per
    emit (~19M allocations / run at N=10): lending-iterator pivot,
    paired with the same fix in `DustFilter`.
  - **H7** — add an isolated `var_calling_dust_filter` criterion bench
    (gates H2 / H3 / L1 / L2 / L12 measurements).
  - **L2–L12** — see the full report; highlights: **L5 is now the
    next ceiling under H1** (`SyncRefFetcher` `RwLock<HashMap>::read()`
    on every fetch — pre-warm into `Vec<Arc<Vec<u8>>>` indexed by
    `chrom_id`, drop the noodles `Repository` runtime dep);
    L2 `AlleleObservation` `Vec` → `SmallVec`,
    L3 PSP CSR decoder, L4 `fetch_from_repository`
    `make_ascii_uppercase`, L8 `[profile.release] debug = "line-tables-only"`
    → `debug = true`, L9 `alloc-mimalloc` A/B against real-data
    workload, L10 missing drained-count assertions in benches.
  - **Per-chrom parallel follow-ups** (from the impl report): streaming
    concat (v2 — append finished fragments while slower chroms run);
    block-level bgzf concat (v2 — skip decompress+re-encode);
    sub-chromosome decomposition (push below
    `max(per-chrom-time)` on imbalanced workloads like tomato's ch00);
    `RLIMIT_NOFILE` bump (defaults to 1024; N=256 × n_chrom > 1024 fds);
    `var-calling-from-bam` + `estimate-contamination` parallelisation;
    `profile_cohort_e2e --em-convergence-threshold` knob; posterior
    `DidNotConverge` emit-with-flag long-term fix; multi-chrom
    integration-test fixture in `tests/common/mod.rs`.
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
- **Latest reviews:** [perf_psp_reader_2026-05-23.md](doc/devel/reports/reviews/perf_psp_reader_2026-05-23.md) (cohort-hot-path re-review against the post-H3 profile — verdict: *Apply the listed wins*; 3 Hot-path / 10 Likely / 4 Speculative; H1 CSR ragged-column collapse subsumes 87% of project-side allocator pressure, H2 = `SeekFrom::Current` for the 2026-05-20 H4 that was never applied, H3 = per-allele bounds-check hoist), [psp_2026-05-13.md](doc/devel/reports/reviews/psp_2026-05-13.md), [psp_reader_2026-05-13.md](doc/devel/reports/reviews/psp_reader_2026-05-13.md), [perf_psp_writer_2026-05-13.md](doc/devel/reports/reviews/perf_psp_writer_2026-05-13.md), [perf_psp_reader_2026-05-13.md](doc/devel/reports/reviews/perf_psp_reader_2026-05-13.md)
- **Latest fixes-applied:** [fixes_applied_psp_reader_2026-05-13.md](doc/devel/reports/implementations/fixes_applied_psp_reader_2026-05-13.md), [perf_psp_reader_2026-05-13_applied.md](doc/devel/reports/reviews/perf_psp_reader_2026-05-13_applied.md), [perf_psp_writer_2026-05-13_applied.md](doc/devel/reports/reviews/perf_psp_writer_2026-05-13_applied.md)
- **Open (from 2026-05-23 PSP reader perf review):**
  - **H1** — CSR collapse of `DecodedBlock`'s
    `allele_seqs: Vec<Vec<u8>>` and
    `allele_chain_ids: Vec<Vec<ChainId>>` into
    `(data: Vec<u8>, offsets: Vec<u32>)` + per-`RecordsIter`
    scratch reuse. Symmetric with the writer's existing
    `encode_list_column_csr`. Stage A keeps `AlleleObservation.seq`
    as `Vec<u8>` (one alloc per emit); Stage B (cross-cutting)
    switches it to `Box<[u8]>` or arena handle. Gate Stage B
    behind dhat data via the missing `examples/dhat_psp_reader.rs`
    (L7).
  - **H2** — `SeekFrom::Current(delta)` helper for the per-block
    seek at `reader.rs:587` and the post-block-header rewind at
    `reader.rs:807`. `SeekFrom::Start` discards `BufReader`'s
    64 KiB user-space buffer on every block transition even
    though blocks are written contiguously on disk; the 2026-05-20
    review's H4. This was never applied.
  - **H3 (gated by H1)** — leading `assert!(allele_end <=
    block.<col>.len())` for the 9 ragged columns at the top of
    `materialise_next_record`. Collapses the per-allele inner
    loop's bounds checks to one each outside the loop.
  - **L1** — `decode_list_column_pod<T: Pod>` LE-slab cast for
    `ChainId = u64` (mirror of the writer's
    `encode_list_column_csr`'s Pod fast path).
  - **L2** — drop `stream_position()` from `read_block_header` by
    threading `entry.block_offset` from the caller (saves one
    `lseek(2)` per block).
  - **L3** — `region_records` first-block seek (subsumed by H2).
  - **L4 (gated by H1)** — pack the 7 fixed-width per-allele
    columns into a single `#[repr(C)]` row for the
    `materialise_next_record` gather.
  - **L5** — `AlleleObservation` SmallVec (cross-cuts the
    2026-05-20 L2; coordinate with that finding).
  - **L6** — `benches/psp_reader_perf.rs` is `Cursor`-backed;
    add a file-backed `BufReader<File>` group so H2 / L1 / L2 /
    L3 have a microbench signal.
  - **L7** — `examples/dhat_psp_reader.rs` (deferred since
    2026-05-13 L9; the post-H3 profile justifies it as a
    prerequisite for H1's Stage B decision).
  - **L8** — bench bodies should `assert_eq!(count, expected)`
    so a silent-truncate refactor can't pass.
  - **S1** — sweep `DEFAULT_BUFFERED_IO_CAPACITY` (64 KiB
    today) after H2 lands.
  - Background-confirmed: L5 (varint fast/cold), L6 (LE-slab cast),
    L8 (BufReader doc) from 2026-05-13 are all still in place.

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

#### Hidden-paralog per-sample summaries (`.psp` metadata section)
- **Status:** shipped on branch `tomato2-paralog-filter` (A1–D3 done; the `pileup` CLI emits the per-sample summary section end-to-end; D2 coverage parity vs the tomato2 prototype confirmed — exact window counts, mean depth within histogram-bin resolution; D3 cost in the noise). Owed: an exact pileup wall/RSS before-after on host CRAM; the downstream filter model that *consumes* the summaries (curve/scale fit, H1/H2 LR, EB/FDR) is a separate follow-on plan.
- **Spec:** [hidden_paralog_filter.md](doc/devel/specs/hidden_paralog_filter.md)
- **Architecture:** [hidden_paralog_psp_integration.md](doc/devel/architecture/hidden_paralog_psp_integration.md)
- **Plan:** [paralog_psp_summaries.md](doc/devel/implementation_plans/paralog_psp_summaries.md)
- **A1 (`.psp` metadata-section container):** done — optional zstd-framed metadata section between the block index and the trailer (trailer byte-identical; located by arithmetic; zip-bomb + trailing-byte guards). Code: [src/psp/metadata.rs](src/psp/metadata.rs), `PspWriter::attach_metadata` / `PspReader::metadata`.
- **B1 (SNP summary payload + TOML serde):** done — `SampleSummary { coverage_by_gc, heterozygosity }` data model serialised as the section's TOML document, with `validate()` invariants + version policy. Code: [src/sample_summary/mod.rs](src/sample_summary/mod.rs).
- **B2 (coverage-by-GC accumulator):** done — single-pass tiled fold producing the histogram. Code: [src/sample_summary/coverage.rs](src/sample_summary/coverage.rs).
- **B3 (het accumulator):** done — binomial het-vs-hom LR three-way classification (confident het / hom-alt / ambiguous) over `SiteCounts`; amended spec §3 + arch Premise 1b/2 + B1 `HetCounts` (two→four counts). Code: [src/sample_summary/het.rs](src/sample_summary/het.rs).
- **C1+C2 (Stage-1 wiring):** done — `--gc-window-bp` flag + `SampleSummaryAccumulators` bundle wired into the [pileup_to_psp.rs](src/pileup/per_sample/pileup_to_psp.rs) seam (record→counts shaping); the `pileup` CLI emits the `.psp` summary section end-to-end (asserted in `pileup_cli_integration`). D1 reader-side accessor is satisfied by composition (`PspReader::metadata()` + `SampleSummary::from_toml_bytes`); no new code (keeps the container schema-agnostic).
- **D2/D3 (parity + cost):** done — coverage histogram re-derived from existing tomato2 `.psp` bodies matches the prototype's `window_cov.w500` (exact window counts, mean depth within 0.07%); re-derive cost ~1.25 s/sample. Tool: [examples/dump_sample_summary.rs](examples/dump_sample_summary.rs). Report: [paralog_summary_parity_2026-06-29.md](doc/devel/reports/implementations/paralog_summary_parity_2026-06-29.md).
- **Latest reviews:** [stage1_summary_wiring_2026-06-29.md](doc/devel/reports/reviews/stage1_summary_wiring_2026-06-29.md), [het_accumulator_2026-06-29.md](doc/devel/reports/reviews/het_accumulator_2026-06-29.md), [coverage_accumulator_2026-06-29.md](doc/devel/reports/reviews/coverage_accumulator_2026-06-29.md), [sample_summary_2026-06-29.md](doc/devel/reports/reviews/sample_summary_2026-06-29.md), [psp_metadata_section_2026-06-29.md](doc/devel/reports/reviews/psp_metadata_section_2026-06-29.md)
- **Latest fixes-applied:** [fixes_applied_2026-06-29_stage1_summary_wiring.md](doc/devel/reports/reviews/fixes_applied_2026-06-29_stage1_summary_wiring.md), [fixes_applied_2026-06-29_het_accumulator.md](doc/devel/reports/reviews/fixes_applied_2026-06-29_het_accumulator.md), [fixes_applied_2026-06-29_coverage_accumulator.md](doc/devel/reports/reviews/fixes_applied_2026-06-29_coverage_accumulator.md), [fixes_applied_2026-06-29_sample_summary.md](doc/devel/reports/reviews/fixes_applied_2026-06-29_sample_summary.md), [fixes_applied_2026-06-29_psp_metadata_section.md](doc/devel/reports/reviews/fixes_applied_2026-06-29_psp_metadata_section.md)
- **Open:** D2 empirical parity vs the tomato2 Python prototype + D3 cost check (both need real tomato2 CRAM/`.psp` data); A1-Mi5 (cap not inspectable — deferred); B1-Mi4 (`bad` closure extract once a 3rd use appears); B2-Mi3 (scheme/histogram field dup — deferred); het/bin params hardcoded (promote to CLI flags when calibration needs it); out-of-scope crate-wide `write_io_err` helper follow-up.

#### Hidden-paralog filter — consumer model + var-calling wiring
The follow-on that *consumes* the per-sample summaries above: the pure statistics
module `src/paralog/` (coverage-model fit, H1-vs-H2 LR, inbreeding scalar,
empirical-Bayes prior + FDR) and the var-calling wiring that scores every locus and
emits a `FILTER`/INFO verdict.
- **Status:** in progress on branch `tomato2-paralog-filter`. **P1 + Q1–Q5 + R1 done; Milestone S wiring in flight (S1–S4 done, S4 checkpoint reached — given a spill, the full calibration is produced RAM-flat). Next: S5 write pass, S6 orchestration+CLI.**
- **Spec:** [hidden_paralog_filter.md](doc/devel/specs/hidden_paralog_filter.md) (§4–§7)
- **Architecture:** [hidden_paralog_locus_statistic.md](doc/devel/architecture/hidden_paralog_locus_statistic.md) (Premises 0–6, statistics core) + [hidden_paralog_varcalling_wiring.md](doc/devel/architecture/hidden_paralog_varcalling_wiring.md) (Milestone S data flow — settled: hard removal, on-by-default, inline per-window coverage fold, ephemeral length-framed-binary spill)
- **Plan:** [paralog_filter_model.md](doc/devel/implementation_plans/paralog_filter_model.md) (P1; Q1–Q5; R1 — **done**) + the detailed wiring plan [paralog_varcalling_wiring.md](doc/devel/implementation_plans/paralog_varcalling_wiring.md) (S1–S6; T1–T2 — next)
- **P1 (callable-position total):** done — `CoverageByGcHistogram.callable_positions` (het-rate denominator) accumulated in Stage-1, round-trips through the `.psp`, `SAMPLE_SUMMARY_VERSION` 1→2. Code: [src/sample_summary/coverage.rs](src/sample_summary/coverage.rs), [src/sample_summary/mod.rs](src/sample_summary/mod.rs). Review: [paralog_p1_callable_positions_2026-07-01.md](doc/devel/reports/reviews/paralog_p1_callable_positions_2026-07-01.md).
- **Q1 (module skeleton + `ParalogModelParams`):** done — `src/paralog/` peer module; fixed model constants/grids (ε, carrier copy numbers, folded-SFS + carrier-freq grids, hom-alt veto) with prototype-faithful defaults. Code: [src/paralog/mod.rs](src/paralog/mod.rs).
- **Q2 (`SingleCopyCoverageModel::fit`):** done — mode-anchored single-copy scale (parabolic refine + mode/median guard), per-GC weighted-median bias curve, σ₀ = 1.4826·MAD floored at histogram resolution; `relative_copy_number`. Code: [src/paralog/coverage_model.rs](src/paralog/coverage_model.rs). Review: [paralog_q2_coverage_model_2026-07-01.md](doc/devel/reports/reviews/paralog_q2_coverage_model_2026-07-01.md).
- **Q3 (`score_locus_for_paralogy`):** done — pure H1-vs-H2 marginal LR (folded-SFS `p` marginal + Wright HWE for H1; carrier config × freq marginal for H2; hom-alt veto falls out). Numerically verified faithful to the prototype. Code: [src/paralog/locus_score.rs](src/paralog/locus_score.rs). Review: [paralog_q3_locus_score_2026-07-01.md](doc/devel/reports/reviews/paralog_q3_locus_score_2026-07-01.md).
- **Q4 (inbreeding scalar `F`):** done — pure `obs_het = n_het/callable` (het *rate*, not the inverting variant-site ratio) + `inbreeding_coefficient = clip(1 − obs_het/hexp, 0, 0.99)`; `Hexp` accumulated by the wiring (S1). Zero-callable → `F` inert (no per-locus data). Code: [src/paralog/inbreeding.rs](src/paralog/inbreeding.rs).
- **Q5 (empirical-Bayes prior + FDR):** done — bounded-RAM `ParalogLrHistogram` + `ParalogPrior::estimate` (EM for π over bins, `converged` flag) + `ParalogFdrCurve` (`q_of_lr` tail FDR + `lr_threshold_for_fdr`). Numerically verified faithful to the prototype. Code: [src/paralog/prior.rs](src/paralog/prior.rs). Review: [paralog_q5_prior_2026-07-01.md](doc/devel/reports/reviews/paralog_q5_prior_2026-07-01.md).
- **Milestone Q complete** — the three pure statistics pieces (coverage model, scorer, prior/FDR) built + unit-tested in isolation.
- **R1 (data-first validation):** done — `examples/paralog_score_parity.rs` + `benchmarks/tomato2/src/paralog_score_parity.py` validate Q2–Q5 on real tomato2: σ₀ 0.282, LR Pearson 0.9931 vs prototype, π 8.88 %, F Spearman 0.851, flagged=paralog profile. Caught + fixed the `Hexp` per-position scale bug (Σ2pq/callable). Report: [paralog_r1_data_validation_2026-07-01.md](doc/devel/reports/implementations/paralog_r1_data_validation_2026-07-01.md). **The maths is trustworthy before wiring (Milestone S).**
- **S1 (`ParalogPrePass` + `HexpAccumulator`):** done — the up-front per-sample state (coverage model + `obs_het` + callable count fit from each `.psp` `SampleSummary`; a missing/rejected sample carried absent, never fatal; `callable_reference()` = cohort-median callable) and the running `Hexp` accumulator (`Σ(1−Σpᵢ²)`, `finish` divides by the callable reference — per-callable-position scale, the R1 correction). Memory-flat (accumulator = two scalars). Code: [src/var_calling/paralog_filter/prepass.rs](src/var_calling/paralog_filter/prepass.rs). Reviewed (approve, doc fixes applied).
- **S2 (ephemeral spill):** done — `ParalogSpill`/`ParalogSpillWriter`/`ParalogSpillReader`: a self-contained length-framed binary store (chosen over the `.psp` container — the payload is per-sample-vector shaped), carrying the full `PosteriorRecord` verbatim + per-sample window so S5 reconstructs a bit-identical record and byte-identity holds by construction; RAII temp-file deleted on drop (success/error/panic); one record resident at a time. Code: [src/var_calling/paralog_filter/spill.rs](src/var_calling/paralog_filter/spill.rs). Reviewed (approve-with-changes; `put_support` exhaustive-destructure + bounded read applied).
- **S3 (per-window coverage pieces):** done — `ReferenceWindowGc` (shared per-window GC from the reference, streamed via the existing memory-flat `StreamingChromRefFetcher`; N-excluded) + `WindowMeanDepthAccumulator` (per-sample window mean depth, free from the `.psp` light columns; tiling byte-identical to `CoverageByGcAccumulator`) + `reference_base_matches` (coordinate-consistency guard). **Design revised with the owner:** GC is a shared *reference* property (not per-sample; FASTA not a per-position sync but a standalone window walk), depth is per-sample; per-sample depth deliberately does *not* N-exclude (documented approximation, T1-validated). Spill record reshaped to `window_gc: f32` + `window_mean_depth: Vec<Option<f32>>`. Code: [src/var_calling/paralog_filter/window_gc.rs](src/var_calling/paralog_filter/window_gc.rs), [window_coverage.rs](src/var_calling/paralog_filter/window_coverage.rs). Reviewed (approve). The producer *wiring* of these moves to S6 (the fold re-folds + is parallel, so a standalone "run unused in the fold" milestone isn't meaningful).
- **S4 (calibrate pass):** done — `score_spilled_locus` (shared, pure, deterministic per-locus scorer reused by S5 for bit-identical recompute) + `calibrate` (streams the spill → `ParalogLrHistogram` → `ParalogPrior::estimate` EM with fixed-prior fallback → `ParalogFdrCurve` + `lr_threshold`). RAM-flat (histogram + one record + O(samples) scratch). `ParalogCalibration::flags` = the FDR drop decision. Code: [src/var_calling/paralog_filter/calibrate.rs](src/var_calling/paralog_filter/calibrate.rs). Reviewed (approve-with-changes; non-finite-LR screen + fallible access applied). **S4 checkpoint: given a spill, π + FDR cut produced RAM-flat.**
- **Carried forward:** downstream `Hobs` consumers (S1) must guard the zero-callable case (`validate()` admits `callable_positions == 0` for an all-`N` sample); Q4 handles it (zero rate → inert `F`). Next: **S5** (write pass — spill → recompute LR → apply cut → VCF + header provenance + `records_dropped_paralog`), then **S6** (orchestration: main-pass spilling, producer window wiring, CLI `--paralog-fdr` / `--no-paralog-filter`, byte-identity-off). S6 also owns: the fallback-warning granularity (empty-histogram vs EM-max-iter), and the operator warning on non-convergence.
- **Latest perf review:** [perf_paralog-filter_2026-07-02.md](doc/devel/reports/reviews/perf_paralog-filter_2026-07-02.md) — the T2 wall follow-up (filter on = +11.6× wall, RSS flat). Verdict *Apply the listed wins*; 5 categories dispatched. The whole +249 s is the pure per-locus scorer run **single-threaded, twice** (calibrate + write) after the 16-way main pass finished — I/O is only ~15 s. Levers, ranked: **H1** parallelise both passes (serial read+join → rayon score-batch → serial fold/emit; histogram fold is order-independent, write pass re-serialises before the writer — byte-identity holds); **L1** memoise the two locus-invariant log-prior tables (H1 Wright 200×58 rebuilt every locus = ~26 B redundant `ln`; H2 carrier 40×58) — pure bit-identical hoist; **L2** hoist `enumerate_carrier_configs` out of the per-locus call; **L3** reuse the six per-locus scratch `Vec`s (per-worker after H1). **S1** score-once + spill the LR to disk (halves work, complements H1, needs owner sign-off — overturns recompute-not-cache). Gating deliverable: add `benches/paralog_scoring_perf.rs` (kernel + streaming) + one `cpu-clock` flamegraph to order the kernel micro-opts. Build profile already optimal (no change). Audit trail `tmp/perf_review_2026-07-02_paralog-filter/`.
- **Open (perf, from the 2026-07-02 review):** **H1** parallelise the calibrate + write passes (the +11.6× fix; RSS-vs-batch-size sweep is the gate — the batch of N owned records is a new RAM cost vs the RAM-flat one-record-at-a-time design). **L1** memoise the Wright/carrier log-prior tables (exact-float-equality gate). **L2** hoist `enumerate_carrier_configs`. **L3** `ScorerScratch` per-worker. **Bench gap:** no criterion bench or sampling profile for the paralog path yet — add both before landing the kernel micro-opts (one hypothesis per measurement).

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

#### Within-chromosome chunk-parallel rewrite (now flat under `src/var_calling/`; formerly `from_psp/`, originally `cohort_block/`)
- **Status:** **superseded (2026-06-05)** — the columnar
  `driver`/`worker`/`loader`/`columns`/`partition`/`two_pass`/`kernels`
  chain this block describes was deleted and replaced by the
  record-streaming pipeline; see the "Re-architected record-streaming
  pipeline" block below. (Historical status follows.) fixes-applied (2026-06-01 review — see "Applied" note below); reviewed (2026-06-01, full-but-prioritized subtree review); parallel block-consume shipped
  (`0d49cf8`, `51b5c63`); **streaming-columnar produce rewrite —
  Stages 1, 2 & 3 implemented**, plus the **DUST worker pool**
  (parallel DUST-ahead). Memory fix landed (Stage 2: N=26 peak RSS
  3963→550 MB); serial DUST moved off the critical path (Stage 3) then
  parallelised across the independent covered intervals (DUST pool). The
  branch now **beats `main` on both wall and RSS** at N=8 (4.1 s / 79 MB
  vs 6.0 s / 142 MB) and N=26 (7.8 s / 291 MB vs 8.9 s / 403 MB), 8
  threads. Optional follow-up: sub-span DUST jobs if a single giant
  interval ever dominates the floor.
- **Plans:**
  - Streaming-columnar produce (current): [cohort_produce_streaming_columnar.md](doc/devel/implementation_plans/cohort_produce_streaming_columnar.md)
  - DUST worker pool (implemented; parallel DUST-ahead — the measured wall floor after Stage 3): [cohort_dust_worker_pool.md](doc/devel/implementation_plans/cohort_dust_worker_pool.md)
  - Master: [cohort_within_chromosome_parallel.md](doc/devel/implementation_plans/cohort_within_chromosome_parallel.md)
  - Phase A.2 column-native EM: [cohort_within_chromosome_parallel_phase_a2_em.md](doc/devel/implementation_plans/cohort_within_chromosome_parallel_phase_a2_em.md)
  - Phase B prereq (variant-bounded chunks): [cohort_within_chromosome_parallel_phase_b1_variant_bounded_chunks.md](doc/devel/implementation_plans/cohort_within_chromosome_parallel_phase_b1_variant_bounded_chunks.md)
  - Phase B prereq (estimate-contamination migration): [cohort_within_chromosome_parallel_phase_b2_estimate_contamination.md](doc/devel/implementation_plans/cohort_within_chromosome_parallel_phase_b2_estimate_contamination.md)
  - Phase B parallel windows: [cohort_within_chromosome_parallel_phase_b_parallel_windows.md](doc/devel/implementation_plans/cohort_within_chromosome_parallel_phase_b_parallel_windows.md)
- **Impl reports:**
  - Phase A (2026-05-28): [cohort_within_chromosome_parallel_phase_a_2026-05-28.md](doc/devel/reports/implementations/cohort_within_chromosome_parallel_phase_a_2026-05-28.md)
  - Streaming produce Stage 1 — span-addressable columnar PSP reader (2026-05-31): [cohort_produce_streaming_columnar_stage1_2026-05-31.md](ia/reports/implementations/cohort_produce_streaming_columnar_stage1_2026-05-31.md)
  - Streaming produce Stage 2 — streaming fold+compact producer / memory fix (2026-05-31): [cohort_produce_streaming_columnar_stage2_2026-05-31.md](ia/reports/implementations/cohort_produce_streaming_columnar_stage2_2026-05-31.md)
  - Streaming produce Stage 3 — DUST-ahead queue / serial DUST off the critical path (2026-05-31): [cohort_produce_streaming_columnar_stage3_2026-05-31.md](ia/reports/implementations/cohort_produce_streaming_columnar_stage3_2026-05-31.md)
  - DUST worker pool — parallel DUST-ahead (2026-05-31): [cohort_dust_worker_pool_2026-05-31.md](ia/reports/implementations/cohort_dust_worker_pool_2026-05-31.md)
- **Code:** [src/var_calling/](src/var_calling/) — the chunk driver now lives flat in the `var_calling` module root: `driver.rs`, `loader.rs`, `worker.rs`, `columns.rs`, `column_span_reader.rs`, `partition.rs`, `test_helpers.rs`, plus `kernels/{mod, unify_alleles, project_scalars, compute_log_likelihoods}.rs`. (Originally `cohort_block/`, then the `from_psp/` submodule; the `from_psp/` layer was merged up on 2026-06-01 once the sibling `from_bam/` path was gone. The streaming-columnar rewrite earlier folded `pre_pass.rs` into the loader/partition path.)
- **Tests:** 88 unit tests in the module (per `cargo test --lib var_calling::cohort_block` at commit `36989d6`); 3 integration tests in [tests/cohort_cli_integration.rs](tests/cohort_cli_integration.rs) (`var_calling_emits_deterministic_vcf_across_runs`, `var_calling_byte_identical_across_worker_windows_per_chunk`, `var_calling_byte_identical_across_target_variants_per_chunk`).
- **Latest fixes-applied:** [cohort_block_2026-05-29_applied.md](doc/devel/reports/reviews/cohort_block_2026-05-29_applied.md) — **Wave 1**: all 5 Blockers Applied (B1 / B2 / B3 / B4 / B5-deferred-to-Wave-2 per Q2) + M5 (bundled with B1) + M11 (bench/example unblocker) + M14 (carryover snapshot helper) + M17 (drop trailing `..`) + M18 (delete dead `chain_id_scratch`) + M19 (`debug_assert!` on sorted `masked_intervals`). 1 026/1 026 lib pass (+3); 21/21 cohort_cli integration pass; fmt clean; criterion baseline saved. 47 findings deferred to Waves 2–3 per Q4 ("apply all structural refactors now"). 2 Won't fix per Q1 (Mi8 / Mi13). Out-of-scope edits flagged in §12.
- **Latest review:** [var_calling_2026-06-01.md](doc/devel/reports/reviews/var_calling_2026-06-01.md) — **Request-changes** (full-but-prioritized review of the whole `src/var_calling/` subtree, 27 489 LoC, commit `3c9ebf2`; prioritized the post-reorg `from_psp/` integration seams, spot-checked the mid-May-reviewed stage internals). 0 Blockers, 13 Major, 10 Minor + grouped Nits. No correctness Blockers — `unsafe_concurrency` verified the parallel-worker soundness holds after the reorg (the `Send + !Sync` REF fetcher never crosses a thread boundary; `seq_idx`-ordered emit + additive stats back byte-identity). Drivers of the verdict: the `cargo doc` CI gate is **red** (M1 — two redundant explicit intra-doc links, `column_span_reader.rs:3` + `posterior_engine.rs:925`); the rewrite's hard byte-identity contract is **unguarded by any test** (M2 — serial-vs-parallel; M9 — no committed golden VCF); and the `from_bam` removal left dead code + dead doc links (M3 dead `ChunkDriverError::PerGroupMerger` variant+`From`; M4 dead `into_shared_ref_fetcher`; M7 three dangling `from_bam` intra-doc links that render as silent dead URLs). Other Major: M5 inert/misleading `chunk_genomic_span` knob (+ phantom `MAX_CHUNK_SPAN_GROWTH` doc ref); M6 startup log omits every per-stage default; M8 `from_psp` public surface should be `pub(crate)` + its facade is bypassed; M10 no perf bench/threshold on the chunk-driver hot path; M11 `BlockQueue` hand-rolled concurrency untested; M12 `DustAheadPool::shutdown` swallows worker panics; M13 unbounded `mpsc` result/recycle channels. Four open questions gate M5/M8/M2/M9 (chunk_genomic_span fate; from_psp public-vs-internal; byte-identity oracle now that the streaming driver is deleted; keep both serial+parallel drive paths?). Per-category audit trail at `tmp/review_2026-06-01_var_calling/`.
- **Latest perf review:** [perf_var_calling_2026-06-01.md](doc/devel/reports/reviews/perf_var_calling_2026-06-01.md) — **Apply the listed wins.** Closes review-finding **M10**: two committed criterion benches added to `benches/var_calling_perf.rs` (`var_calling_chunk_loader` = produce, `var_calling_run_window` = consume kernels + EM), swept N ∈ {8,64,200}. Evidence: bench (loader ~linear in N, **`run_window` ~N^1.7**), CPU sampling profile (macOS `sample`/`samply`), DHAT (real cohort N=18: 15.5 GB / 51.4 M blocks), `cargo-show-asm`. **3 Hot-path:** **H1** the QUAL allele-count convolution `compute_qual_via_exact_af` (posterior_engine.rs:3037-3047) does a *scalar pairwise* `log_sum_exp_2` in an O(n_samples²·ploidy²) loop = ~51 % of `run_window` self-time and the mechanism behind the N^1.7 scaling (the hot `log_sum_exp_2` is the QUAL convolution, **not** the EM, which takes the homogeneous-fixation fast path under default config); **H2/H3** the two `unify_alleles` allocation sites (`build_chain_proposals_columnar` BTreeMap value-Vecs ~15.7 M blocks — M3's scratch-hoist is defeated by `BTreeMap::clear` dropping the Vecs; `project_per_position_into_scratch` byte-key clones ~7.3 M blocks). 11 Likely (mimalloc on the production binary; per-group genotype-table cache; thread oversubscription producer-rayon vs worker-threads; DUST `to_vec` copies; `for_contig` re-parsing `.fai` per interval; AoS→SoA on the LH-hot scalars; bench hardening; …), 4 Speculative, Notes. Build config already tuned (`lto=fat`, `codegen-units=1`, `panic=abort`, pinned `target-cpu`); production allocator is the one build lever. Per-category audit trail + profiles at `tmp/perf_review_2026-06-01_var_calling/`.
- **Latest fixes-applied (2026-06-01, branch `var-calling-review-fixes`):** PM-directed subset of the review, in two phases. **Phase 1** (commits `d2a1397` review-fixes + `31ddf5a` doc cleanup): **M1** applied (+ all 8 in-scope broken/redundant intra-doc links, and a bonus crate-wide pass fixing 14 *pre-existing* broken links outside `var_calling` — `cargo doc --no-deps` is now green for the whole crate; the original review under-reported this, having tailed the doc output); **M3** + **M4** applied (dead `PerGroupMerger` variant/`From` + dead `into_shared_ref_fetcher` deleted); **M5** applied (inert `ChunkSizingParams.chunk_genomic_span` field deleted — `DEFAULT_CHUNK_GENOMIC_SPAN` kept; the contamination chunked-stream genuinely uses it); **M7** applied (dangling `from_bam` doc links repointed/de-linked). Plus the PM's **always-parallel** decision: `drive_blocks_serial` removed, `drive_cohort_chunked` always runs `drive_blocks_parallel` with `n_workers = current_num_threads().max(1)` — which **moots M2** (no two paths to keep byte-identical). **Phase 2** (commit `c068af4`): the `from_psp/` → `var_calling/` merge above. **Dispositions:** **M2** moot (serial path removed); **M9** won't-fix as written — byte-identity is verified out-of-tree (human-in-the-loop benchmark + branch-vs-`main` diff during development), not via a committed in-tree golden VCF. Verification (in container): `fmt --check`, `clippy --all-targets --all-features -D warnings`, `doc --no-deps` all clean; `cargo test --lib` 1059 pass; cohort integration 17 pass. **Phase 3** (commits `79a6856` + `78a357b` + `cac838f`): **M6** (startup log now dumps every per-stage config), **M11** (5 `BlockQueue` concurrency tests), **M12** (`DustAheadPool::shutdown` logs a panicking worker), **M13** (mpsc in-flight bound documented), **Mi1** (Display source double-render), **Mi7** (`fill_block` doc), **Mi8** (Relaxed-atomic ordering comment) applied. **M8** attempted then reverted — `pub(crate)` un-masks a dead / test-only-code layer (~15 items) that needs its own removal pass first; the rationale is recorded at the module surface. Final verification: `fmt --check` / `clippy --all-targets --all-features -D warnings` / `doc --no-deps` clean; `cargo test --lib` 1064 pass; cohort integration 17 pass. **Still open:** M8 (dead-code cleanup → then internalise), M10 hot-path benches **added** (2026-06-01 perf review — `var_calling_chunk_loader`/`var_calling_run_window`), so M10 is satisfied and the perf wins it surfaced (H1–H3) are now the follow-up, Minors Mi2/Mi3/Mi4/Mi5/Mi6/Mi9/Mi10, and the out-of-scope `benches/psp_writer_perf.rs:385` panic.
- **Prior review:** [cohort_block_2026-05-29.md](doc/devel/reports/reviews/cohort_block_2026-05-29.md) — **Request-changes**: 5 Blockers (B1 writer-tmp leak on driver-error path; B2 full-chrom `Vec<u8>` materialisation in `compute_dust_mask_for_chrom` defeats per-chunk memory contract; B3 `NoSafeGap` retry is a no-op when `target_variants_per_chunk > 0`; B4 `NAllelesExceedsBitmask` silently rewritten as `DegenerateLikelihood { usize::MAX, … }`; B5 missing cross-driver byte-identity oracle test + missing unit tests for `drive_cohort_chunked` / `drive_one_chrom_generic` / `load_and_run_chunk_with_retry` / `emit_or_drop` / `compute_dust_mask_for_chrom`), 32 Major, 26 Minor, grouped Nits. `unsafe_concurrency` returned `No findings.` — the parallel-section soundness is statically enforced by the `Send + !Sync` typedef on `SharedRefFetcher`. Four open questions for the author (stable-API intent on the new pub data structs; streaming `drive_cohort_pipeline` oracle's long-term fate; sentinel-vs-`NonZero` policy for `target_variants_per_chunk` / `target_window_count`; filter-order equivalence vs streaming pipeline) gate several Major findings. Per-category audit trail at `tmp/review_2026-05-29_cohort_block/`.
- **Open (from the 2026-06-01 review — see report §6/§8 for full text + fixes):**
  - **M1** — Drop the explicit `(path)` target on the two redundant intra-doc links (`column_span_reader.rs:3`, `posterior_engine.rs:925`) to turn the `cargo doc` gate green.
  - **M2** — Add a serial-vs-parallel byte-identity test (drive the same fixture through `drive_blocks_serial` at 1 thread and `drive_blocks_parallel` at N threads; assert VCF bytes + every `ChunkDriverStats` field equal).
  - **M3 / M4** — Delete the dead `ChunkDriverError::PerGroupMerger` variant + `From` impl + import (`driver.rs:296-323`, `:52`); delete the dead `into_shared_ref_fetcher` (`worker.rs:702-715`, `mod.rs:41`).
  - **M5 / M7** — Resolve the inert `chunk_genomic_span` knob (delete vs. mark reserved; drop the phantom `MAX_CHUNK_SPAN_GROWTH` doc ref) per Q1; repoint or de-link the three dangling `from_bam` intra-doc links (`driver.rs:23,150`, `per_group_merger.rs:25`).
  - **M6** — Extend the startup config log to dump every effective per-stage default (DUST / grouper / per-group / `PosteriorEngineConfig`).
  - **M8** — `pub(crate)` sweep over the `from_psp` submodules + re-exports (gated on Q2); settle `DEFAULT_CHUNK_GENOMIC_SPAN` on one import path.
  - **M9 / M10** — Commit a golden VCF + byte-identity test (gated on Q3); add a `from_psp` chunk-driver + columnar-kernel bench with `// REGRESSION THRESHOLD: N%`.
  - **M11 / M12 / M13** — Add `BlockQueue` concurrency tests (mirror the `DustAheadPool` suite); surface/log `DustAheadPool::shutdown` worker-panic joins; document or `sync_channel`-bound the unbounded `mpsc` result/recycle channels.
  - **Minors** — `Display` source double-render (Mi1); `DownstreamFilterParams`/`ChunkSizingParams` default-binding docs (Mi2); co-dependent mapq-filter fields (Mi3); duplicated `DRIVER_PSP_BUFFER_BYTES`/`MAPQ_FILTER_MIN_READS_PER_SIDE` consts (Mi4); reconcile `chunk` vs `block` and stale `window`/`chunk.windows` vocabulary (Mi5/Mi6); `fill_block` doc says `bool` but returns `u32` (Mi7); `Relaxed`-atomic ordering comment (Mi8); direct `find_block_cut` tests (Mi9); malformed-PSP decode tests (Mi10). Plus the out-of-scope `benches/psp_writer_perf.rs:385` panic (breaks `cargo test --all-targets`) and the `genotype_order` placement question.
- **Open (from the 2026-05-29 review):**
  - **B1** — Add `CohortVcfWriter::abort()` that takes the tmp path it actually used and removes it; call from the error branch. Add an integration test injecting a mid-loop error and asserting no leftover tmp on disk.
  - **B2** — Stream `fetcher.iter_bases()?` directly into `sdust_mask_streaming` (the helper already accepts `Iterator<Item = io::Result<u8>>`). Add a regression test using a synthetic fetcher that panics on `.collect::<Vec<_>>()`.
  - **B3** — Thread `attempt_span` into the loader's `max_span` so retries can actually load more data, or push the retry inside the loader as "extend until safe gap or chrom cap". Add a test for a chunk that satisfies `target_variants_per_chunk` on the first attempt but has no safe gap.
  - **B4** — Add a dedicated `PerGroupMergerError::NAllelesExceedsBitmask { n_alleles, chrom_id, group_start, group_end }` variant; drop the `let _ = n_alleles;`. Add a test triggering the path with `cfg.max_alleles > MAX_BITMASK_ALLELES`.
  - **B5** — Add one integration test in [tests/cohort_cli_integration.rs](tests/cohort_cli_integration.rs) that runs both `drive_cohort_chunked` and the streaming `drive_cohort_pipeline` on a multi-position fixture (≥3 samples, ≥1 MNP, ≥1 LH-cap site, ≥1 hom-REF group, ≥1 below-`min_alt_obs` site, ≥1 below-`qual_phred` site, ≥1 above-`mapq_diff_t` site); assert VCF bodies byte-equal **and** field-by-field equal counter sets. If the streaming driver is scheduled for removal, capture a checked-in golden VCF first.
  - **M1** — Replace `#[derive(Default)]` on `SampleColumns` with a hand-written `impl Default for SampleColumns { fn default() -> Self { Self::empty() } }`.
  - **M2** — Rewrite `prefetch_window_ref_bytes` so the outer `Vec<Vec<u8>>` resizes-without-dropping and each inner `Vec<u8>` is cleared-in-place rather than freshly allocated.
  - **M3** — Move `detect_compound_candidates_columnar`'s two `BTreeMap`s into `UnifyAllelesScratch` for per-group reuse; add a permutation-invariance proptest.
  - **M4 / M5 / M21** — Reshape `ChunkDriverError`: drop `#[from]` on `Io` / `PspRead`; rename variants by operation (`OpenPsp`, `WriteVcf { chrom_id, start, end, … }`, `FetchRefBases`, …); drop `: {0}` interpolation; add chrom/range context to `WriteVcf`. Replace `let _ = std::fs::remove_file(...)` with a structured `tracing::warn!` event.
  - **M6** — Surface `u32_from_usize` as a typed `ChunkLoadError::CsrOffsetOverflow` error or use `try_into().expect(...)` to panic loudly with a named invariant rather than wrap silently.
  - **M7 / M9 / Mi20 / Mi21** — Either propagate `ZeroTargetWindowCount` from the pre-pass (dropping the `.max(1)`) or name the default with `pub const DEFAULT_WORKER_WINDOWS_PER_CHUNK` and emit `tracing::debug!` when applied. Same shape for `target_variants_per_chunk == 0` (name `TARGET_VARIANTS_DISABLED` or lift to `Option<NonZeroU32>`). Add startup `tracing::info!` listing every effective `ChunkDriverParams` value. Add `chunks_with_fewer_windows_than_requested: u64` counter to `ChunkDriverStats`.
  - **M8** — Either reject `max_span < initial_span` with a new `ChunkLoadError::MaxSpanBelowInitial` variant or delete the `effective_initial_span` no-op chain and document the invariant.
  - **M10 + Mi26** — Fix the 2 in-scope unresolved-link errors and 3 redundant-link warnings (`worker.rs:17`, `worker.rs:236`, `driver.rs:105`, `driver.rs:108`, `worker.rs:9`).
  - **M11** — Add the two new `VarCallingArgs` fields to `benches/cohort_e2e_perf.rs:286`, `examples/profile_cohort_e2e.rs:152`, `examples/dhat_var_calling.rs:121` with the legacy defaults (`target_variants_per_chunk: 0`, `worker_windows_per_chunk: 1`). Optional follow-up: add a `VarCallingArgs::for_profiling(...)` constructor.
  - **M12 / M13 / M14 / M15 / M16 / M17 / M18 / M19 / M20 / M22 / M23 / M24 / M25 / M26 / M27 / M28 / M29 / M30 / M31 / M32** — see the report's §6 Findings section and §8 Missing tests. Highlights: split `load_and_run_chunk_with_retry`'s 19-param body into three phase helpers; group `load_chunk_from_iters`'s span/variant knobs into `ChunkLoadExtent`; add `SampleColumns::clone_from_columns` and use it for both carryover snapshot / restore loops; split per-window counters out of `ColumnarPipelineScratch`; drop the trailing `..` from the `AlleleSupportStats` destructure at `columns.rs:116`; delete or wire `chain_id_scratch` (`#[allow(dead_code)]`); validate `masked_intervals` sorting in `partition_window`; pin filter order in `emit_or_drop` with five per-category unit tests; pin `enforce_max_alleles_columnar` tie-break against the row-shape kernel; clamp `safe_end` to `chrom_one_past_end` on the last chunk; rename `*_cfg` vs `*_config` to a single form crate-wide; rename `shared_ref_fetcher` to `into_shared_ref_fetcher`; add tests for `SampleCountMismatch` / `CarryoverLengthMismatch` (both loader and pre-pass); add a `par_iter_mut` vs sequential equivalence test.
  - **Mi-class** (~26 minors): `#[non_exhaustive]` on the new pub data structs (gated on Open Question 1); `pub mod` → `pub(crate) mod` for every submodule that has no out-of-crate consumer (gated on Open Question 1); rename `MaterialisedChunk::clear_data` → `clear`; rename `WorkerSlot.output_buf` / `WorkerSlot.scratch` to carry domain nouns; collapse `Arc::new(StreamingChromRefFetcher)` to a borrow; drop `chunk.windows.clone()`; take `&PosteriorEngineConfig` in `run_window`; remove double-clone in `push_allele_into_scratch`; demote `pub` items with no caller; merge `Ok(idx) | Err(idx)` arms; split `unify_alleles.rs`/`worker.rs`/`loader.rs` along their existing internal sub-step boundaries; move `build_overlapping_variant_group` out of `worker.rs` and into `test_helpers.rs`; convert the three `Vec<Vec<_>>` jagged arrays to CSR; consider `OneBasedPos` / `OneBasedRange` / `ChromId` newtypes; group `ChunkDriverParams` along stage boundaries; add `// REGRESSION THRESHOLD: N%` to `benches/cohort_e2e_perf.rs`; add `--ignored` should-panic regression for `u32_from_usize` overflow.
  - **Nits** — single mechanical pass to clear the 16 in-scope clippy errors (`single_range_in_vec_init` ×8, `type_complexity` ×4, `bool_assert_comparison` ×2, `doc_lazy_continuation` ×2) plus add per-call-site justification comments to the 14 `#[allow(...)]` annotations (12 `clippy::too_many_arguments` + `clippy::arc_with_non_send_sync` + `clippy::needless_range_loop`).

#### Re-architected record-streaming pipeline (replaces the chunk-parallel rewrite)
- **Status:** fixes-applied (2026-06-08) — the 2026-06-08 review applied on
  branch `var-calling-review-fixes` (6 commits: M1/M4/M6/Mi6 deletions+move; M3
  typed errors; Mi4 rename + M7/Q1 doc sweep + Minors; M5/M8 tests; **M2**
  `pub(crate)` demotion + the dead-code cleanup it unmasks — delete fully-dead
  `is_empty`/`clear`, `#[cfg(test)]`-gate the eager-decode oracle the prior M8
  named as a prerequisite). All gates green (fmt / clippy `-D warnings` / `doc`;
  1021 lib + integration tests). Open questions resolved with the PM: P7-swap
  doc debt resolved now (Q1), dead code deleted (Q2), cheap in-tree tests over
  the out-of-tree oracle (Q3), modules demoted to `pub(crate)` (Q4). Deferred:
  the harder M8 leftovers (`StalledCut` / `compact_samples` straddler /
  `dust_mask_for_interval` sub-span paths + a multi-thread e2e byte-identity
  test). Prior: reviewed (2026-06-08) —
  re-reviewed on `main` post-merge (subtree review excluding `posterior_engine`;
  verdict Approve-with-changes, 0 Blockers / 8 Major / 13 Minor; see Latest
  review). Prior: fixes-applied
  (2026-06-05) — all 8 Major review findings +
  2 Minors applied across 3 commits (`b4e767c` review, `8210f46`
  M1/M2/M5/Mi2, `ed141ff` M3/M4/M6/M8/Mi8). Shipped to production on branch
  `re-architect` (Phase 7 swap, `1d34f85`). The columnar driver/worker/
  loader/columns/partition/two_pass/kernels chain is **deleted**; the only
  cohort `.psp` → VCF path is now the three-component record-streaming
  topology (producer `CohortChunkIntegrator` → W `VariantCaller` callers →
  `VcfWriter`), wired with two bounded crossbeam channels inside a
  `std::thread::scope`. Goal was **clarity** (runtime stages map to named
  logical sections), with byte-identical calls vs `main` as the hard
  contract and memory/wall as a measured guardrail.
- **Spec / plan:**
  [re_architecture_streaming_pipeline.md](doc/devel/implementation_plans/re_architecture_streaming_pipeline.md)
  (architecture, constraints, build strategy),
  [re_architecture_module_outline.md](doc/devel/implementation_plans/re_architecture_module_outline.md)
  (module/type map),
  [re_architecture_execution_plan.md](doc/devel/implementation_plans/re_architecture_execution_plan.md)
  (the byte-identity-gated phases),
  [re_architecture_p6_measurement.md](doc/devel/implementation_plans/re_architecture_p6_measurement.md)
  (Phase 6 vs `main` at scale).
- **Code:** [src/var_calling/](src/var_calling/) — new modules `types.rs`,
  `sample_reader.rs`, `cohort_integration.rs`, `pileup_overlaps.rs`,
  `em_posterior_calc.rs`, `vcf_writer.rs`, `pipeline.rs`, rewritten
  `mod.rs`; the numeric kernels (`per_group_merger`, `posterior_engine`,
  `variant_grouping`, `per_position_merger`, `dust_filter`,
  `contamination_estimation`) are carried verbatim. CLI rewired in
  [src/pop_var_caller/var_calling.rs](src/pop_var_caller/var_calling.rs);
  contamination estimator ported to the record-based `PerPositionMerger`
  in [src/pop_var_caller/estimate_contamination.rs](src/pop_var_caller/estimate_contamination.rs).
- **Tests:** 996 lib + all cohort integration pass (in container, commit
  `ef93b67`); `fmt --check` + `clippy --all-targets --all-features -D warnings`
  clean.
- **Latest review:** [var_calling_2026-06-08.md](doc/devel/reports/reviews/var_calling_2026-06-08.md)
  — **Approve-with-changes** (orchestrator skill, 10 categories; `main` @
  `35a6b67`; subtree review excluding the PM-deferred `posterior_engine`). 0
  Blockers, 8 Major, 13 Minor + Nits. All gates green (fmt / clippy `-D warnings`
  / `doc`; 1052 lib + integration tests pass; only `--all-targets` failure is the
  pre-existing out-of-scope `psp_writer_perf` bench panic). `unsafe_concurrency`
  clean (deadlock-free channels, order-independent parallel fold,
  `VariantCaller: Sync`). Drivers: M1 dead `PileupCohortChunk`, M2 over-broad
  `pub` surface (extends prior M8), M3 `PipelineError::Config(String)`
  flattening, M4 caller→producer coupling (`merge_compacted_samples`), M5
  no worker-count byte-identity test, M6 dead/mislabeled `SamplePspChunk` API,
  M7 stale verbatim-provenance pointing at deleted modules, M8 untested
  release-guards/byte-identity helpers. Per-category audit trail at
  `tmp/review_2026-06-08_var-calling/`.
- **Previous review:** [re_architecture_pipeline_2026-06-05.md](doc/devel/reports/reviews/re_architecture_pipeline_2026-06-05.md)
  — **Request-changes** (orchestrator skill, all 11 categories; HEAD
  `ef93b67`). 0 Blockers, 8 Major, 16 Minor, grouped Nits. No correctness
  Blockers in emitted calls; `unsafe_concurrency` clean (no `unsafe`,
  sound channel-close ordering, order-independent parallel-decode/serial-fold);
  `refactor_safety` verified `derive_is_kept`/`find_cut`/`chunk_cuts`/
  `merge_block_ranges`/`emit_or_drop`/`passes_min_alt_obs`/
  `merge_group_with_ref` are line-by-line identical to `main`. Per-category
  audit trail at `tmp/review_2026-06-05_re-architect-pipeline/`.
- **Latest perf review:** [perf_var_calling_cohort_2026-06-06.md](doc/devel/reports/reviews/perf_var_calling_cohort_2026-06-06.md)
  — **Run experiments** (orchestrator skill, 6 categories; HEAD `37e02c2`).
  First CPU sampling profile collected at this commit (macOS `sample`, T=1 30 s
  + T=8 6 s, N=50 real tomato cohort, tvpc=256). **3 Hot-path:** **H1** the T≥2
  wall gap vs `main` is *scheduler oversubscription* — the producer's rayon
  compaction pool and the crossbeam caller pool are two independent pools each
  sized to `--threads` (producer thread **64 % blocked on rayon** at T=8;
  `swtch_pri`/`__ulock_wait2`/`mutexwait` elevated vs T=1; tracks the
  on-par-T=1 / +6/+12/+15 %-T=2/4/8 matrix); **H2** the chain-id dead-weight
  (REF-allele chain_ids = 65 % of in-flight payload, ~96.6 % of all chain_ids,
  never read — `per_group_merger.rs:1304` skips allele 0; verified by 3
  categories); **H3** `ln_factorial` not `#[inline]` (out-of-line `bl` on the
  likelihood triple-loop, `cargo asm`-confirmed; cheap byte-identical apply).
  7 Likely (producer serial-on-main floor; `e_step_simd` sample-axis gather;
  `records_all` per-allele clones; per-chunk REF-span `Vec`; `Vec<Option<>>`
  AoS merge scan; E-step `is_finite` branch hoist; `binary_search` keep-mask →
  merge-walk), 4 Speculative, Notes. **No cohort end-to-end criterion bench
  exists** (deleted with from-bam) — building it is measurement-plan item 1.
  Build config already tuned (`lto=fat`/`codegen-units=1`/`panic=abort`);
  allocator A/B is the one unstruck build lever. Per-category audit trail +
  profiles at `tmp/perf_review_2026-06-06_var-calling-cohort/`.
- **Applied (2026-06-05, commits `8210f46` + `ed141ff`):** **M1** (doc gate
  green: 7 dead intra-doc links repointed/removed + 5 redundant targets
  dropped; verified under `RUSTDOCFLAGS=-D warnings`), **M2** (restored
  `current_command_line()` for the `##commandline` header), **M3**
  (zero-allele `.psp` rejected at the decode boundary via new
  `BlockHeaderInvariantKind::ZeroAlleleRecord` + `validate_n_alleles_column`
  helper + 4 unit tests), **M4** (typed `ChromRefFetchError` boxed through
  `ProducerError::Ref` — new `pub type RefFetchError`), **M5** (dropped the
  crate-wide `#![allow(dead_code)]` + deleted the 3 zero-caller items
  `into_reader`/`n_positions`/`chunk_cuts`), **M6** (`ProducerError::StalledCut`
  + `WriterError::MissingChunks` replace the two release-load-bearing
  `debug_assert!`s), **M8** (named `QUEUE_DEPTH_PER_WORKER`, startup log of
  resolved `workers`/`queue_cap`/`target_variants_per_chunk`, rewrote stale
  `--threads`/`--target-variants-per-chunk` CLI docs), **Mi2** (removed the
  inert `--target-variants-per-chunk` from `estimate-contamination` + its
  vacuous test), **Mi8** (refreshed stale "Phase 4 / `!Send`" module docs).
  Verified in container: `fmt` / `clippy --all-targets -D warnings` / `doc -D
  warnings` clean; 1000 lib + cohort integration tests pass.
- **Open (from the 2026-06-08 code review — see report §6/§8 for full text + fixes):**
  - **M1** — delete the dead `pub` `PileupCohortChunk` (types.rs); retarget its 3 doc refs at `RawCohortChunk`.
  - **M2** — demote the 5 zero-external-consumer rebuilt modules (`types`/`sample_reader`/`cohort_integration`/`pileup_overlaps`/`em_posterior_calc`) to `pub(crate)`; consider re-exporting `WriterStats` from `pipeline` to demote `vcf_writer` too.
  - **M3** — replace `PipelineError::Config(String)` flattening with operation-named variants carrying each typed cause via `#[source]`.
  - **M4** — move `merge_compacted_samples` (+ `ok_record`/`KeptRecordIter`) from the producer module into the caller.
  - **M5** — add a `var_calling_byte_identical_across_worker_counts` integration test + a direct `VcfWriter` out-of-order reorder unit test.
  - **M6** — delete dead `append_kept`/`SamplePspChunk::chrom_id()`; `#[cfg(test)]`-gate or relabel the eager `from_block`/`take_*` oracle path; fix the misdirecting module doc.
  - **M7** — pin each "copied verbatim from `two_pass`/`driver`/`worker`/`loader`" provenance to the source commit hash; rewrite the migration-future-tense `mod.rs`/`types.rs`/kernel docs to present tense (P7 swap has happened); resolve the deferred `Variant`-alias decision; fix stale `doc/devel/...` doc links.
  - **M8** — add the targeted unit tests for `MissingChunks`, `StalledCut`, `emit_or_drop` ordering/counters, `overlapping_groups`, `restrict_intervals_to_regions`/`dust_mask_for_interval`, `rebuild_fold` reduce-order independence, and the `compact_samples` straddler.
  - **Minors** — Mi1 startup-log provenance tags; Mi2 `BUFFERED_IO_CAPACITY` const drift; Mi3 collapse `DownstreamFilters` mapq fields to an enum; Mi4 rename `em_posterior_calc`→`variant_caller`; Mi5 `CohortPileupRecord`/`PileupCohortChunk` near-anagram (resolved by M1); Mi6 `RefSpan::empty()` fictional-doc; Mi7 `CallStats`→`WriterStats` exhaustive-destructure; Mi8 `// PANIC-FREE:` comments on release `.expect()`s; Mi9 caller-panic→typed-error vs fatal-by-design; Mi10 `SamplePspReader::new(r,0,1,1)` placeholder comment; Mi11 `sd`/`so`/`cd`/`co` scratch names; Mi12 hot-path bench regression threshold; Mi13 `passes_min_alt_obs` layout-coupling.
  - **Out of scope (follow-ups):** `dust_filter::is_masked` release-path passthrough on an unloaded mask (verbatim kernel correctness); the pre-existing `psp_writer_perf.rs:386` bench panic; the deferred `posterior_engine` review.
- **Open (perf — from [perf_var_calling_cohort_2026-06-06.md](doc/devel/reports/reviews/perf_var_calling_cohort_2026-06-06.md), verdict Run experiments):**
  - **Bench gap (do first)** — add `benches/cohort_var_calling_perf.rs` (criterion,
    `harness=false`) sweeping N∈{1,8,50}×T∈{1,2,8} at tvpc=256; without it no
    code-level perf finding is rankable cross-commit. Also fix
    `examples/profile_cohort_e2e.rs:161` (`target_variants_per_chunk: 0` →
    `256`) so the maintained driver profiles the production shape.
  - **H1 (wall gap, highest priority)** — collapse the rayon-compaction ⟂
    crossbeam-caller two-pool oversubscription; sweep cap-rayon-below-`--threads`
    / serialize-`compact_samples` / move-compaction-onto-callers, gated on
    byte-identical calls + a re-`sample` showing the oversubscription frames fall.
  - **H2 (memory win)** — drop the never-read REF-allele chain_ids (step 1:
    `records_all`/`records_for` REF slot → empty `Vec`; step 2: CSR stores
    non-REF alleles only). DHAT-gated, zero-VCF-diff gate.
  - **H3 (cheap apply)** — `#[inline]` `ln_factorial` + `#[cold]` tail; pure codegen, byte-identical.
- **Open (deferred — lower-priority Minors + the test additions):**
  - **M7** — Add unit tests for the writer reorder buffer (permuted +
    buffered-future-chunk; the `MissingChunks` gap path is now testable),
    `emit_or_drop` per-gate ordering/counters, a multi-thread
    `read_samples`-vs-reference test, and the new `StalledCut` guard. (PM:
    byte-identity is verified out-of-tree vs the previous version + the GIAB
    benchmark, so no in-tree A/B oracle is needed.)
  - **Minors** — dropped `chunks_loaded`/`avg_variants_per_chunk` run
    summary (Mi1, needs a chunk counter in `WriterStats`); ploidy-vs-inverted-range
    error precedence test in `merge_group_with_ref` (Mi4); missing
    `// PANIC-FREE:` on the thread joins (Mi5) and `fetch_ref_span`
    `binary_search().expect()` (Mi6); `ProducerError::Merge` not `#[source]`
    (Mi7); const-via-comment drift (Mi9); two `0 ⇒ ?` target conventions
    (Mi10); `em_posterior_calc` named after a sub-step (Mi11);
    `CohortPileupRecord`/`PileupCohortChunk` near-anagram (Mi12);
    co-dependent mapq-filter fields (Mi13); pipeline → CLI `VarCallingArgs`
    back-reference (Mi14); stale Cargo.toml crossbeam comment path (Mi15);
    unguarded `max_reach - first + 1` (Mi16). Plus the out-of-scope
    `benches/psp_writer_perf.rs:386` panic (breaks `cargo test --all-targets`;
    not the CI gate, which uses `--lib --tests`).
  - **Open questions — resolved by the PM (2026-06-05):** the `##commandline`
    change was unintentional → fixed (M2); `--target-variants-per-chunk`
    should not stay a no-op → removed from `estimate-contamination` (Mi2);
    byte-identity is verified out-of-tree (vs the previous version + the GIAB
    benchmark), so no in-tree oracle is needed; dead code should be removed →
    the `#![allow(dead_code)]` and its dead surface are gone (M5).

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

## SSR/STR caller (independent pipeline)

A second, independent caller — microsatellite/STR **length** genotyping from
aligned reads — sharing only low-level alignment I/O and a few numerical kernels
with the SNP caller above, never its records or math. Mirrors the SNP shape:
`ssr-catalog → catalog → ssr-pileup → .ssr.psp → ssr-call → VCF`. Design:
[ssr_genotyping.md](doc/devel/specs/ssr_genotyping.md) (model + why) and the
architecture docs ([overall](doc/devel/architecture/ssr_genotyping_architecture.md),
[shared types](doc/devel/architecture/ssr_shared_types.md),
[Stage 0 catalog](doc/devel/architecture/ssr_catalog.md)). Decision E and the
type model are settled; built in data-flow order (types → Stage 0 → Stage 1/2).

### Stage 0 — `ssr-catalog` (reference → catalog)

#### Shared domain types (`src/ssr/types.rs`)
- **Status:** fixes-applied
- **Design:** [ssr_shared_types.md](doc/devel/architecture/ssr_shared_types.md)
- **Code:** [src/ssr/types.rs](src/ssr/types.rs), [src/ssr/mod.rs](src/ssr/mod.rs) (commit `74b5a2d`; no separate impl report — commit message carries context)
- **Tests:** 16 unit tests in the module (6 original + 10 boundary/error tests from the fix run)
- **Latest review:** [ssr_types_2026-06-12.md](doc/devel/reports/reviews/ssr_types_2026-06-12.md) — Request-changes (1 Blocker, 1 Major, 5 Minor + missing tests)
- **Latest fixes-applied:** [fixes_applied_2026-06-12.md](doc/devel/reports/reviews/fixes_applied_2026-06-12.md) — **Completed**: all 9 findings resolved (B1 + M1 + Mi1–Mi5 Applied; Nit1 kept; Nit2 already consistent). M1 → private `Locus` fields + validated `Locus::new`/`LocusError`; Mi5 → `pub(crate)` surface (with a temporary module `#![allow(dead_code)]` until Stage 0 consumers land). `cargo doc`/ssr clippy/fmt clean; 1039 lib tests pass.
- **Open:**
  - Remove the temporary `#![allow(dead_code)]` in `src/ssr/mod.rs` once `ssr-catalog` (Stage 0) wires these types up (tracked as Mi5 follow-up in the fix report).

#### Catalog builder (`src/ssr/catalog/`) — *planned*
- **Status:** planned
- **Plan:** [ssr_catalog.md](doc/devel/implementation_plans/ssr_catalog.md) (implementation sketch: files, structs, fn signatures)
- **Notes:** detector = lh3/TRF-mod (shell-out via temp files, no FFI), genome-wide; post-process drops compound/bundled loci (GangSTR-style, no split); worker-per-contig with an ordered collector; `--num-chroms-in-parallel` is a speed⇄RAM knob.

### Stage 1 — `ssr-pileup` (per-sample evidence extraction)

#### `ssr-pileup` stage
- **Status:** **implemented — rebuilt as Mark-2** (2026-06-17, branch `ssr-pileup-mark2`). The Mark-1 detail below (the reference-anchored rung model + its perf/review history) is **superseded**, kept for chronology.
- **Mark-2 rebuild (current):** [ssr_pileup_mark2_2026-06-17.md](ia/reports/implementations/ssr_pileup_mark2_2026-06-17.md). Empirical-candidate model — observed sequences, the reference is only a coordinate frame, no on/off-ladder, no Stage-1 likelihood. Design: [ssr_ladder_model.md](doc/devel/architecture/ssr_ladder_model.md) + [ssr_pileup_mark2.md](doc/devel/architecture/ssr_pileup_mark2.md); plan [ssr_pileup_mark2.md](doc/devel/implementation_plans/ssr_pileup_mark2.md). New `src/ssr/` (`types`/`catalog`/`pileup::{footprint,fetch_reads,alignment,locus_tally,driver}`) + rewritten `psp/registry_ssr.rs`; old `src/ssr_mark1/` + its bench/example deleted. fmt/clippy `-D warnings` clean; 1116 lib + integration + doctests green (end-to-end + thread-determinism). **Open:** calibrate `MIN_REGION_Q1` (15) / `MAX_READS_PER_LOCUS` on real data; a Mark-2 bench (the Mark-1 one was deleted); spec §4.2/§4.3/§5.1 amendment; then Stage 2 (`ssr-call`).
- **Latest fixes-applied:** [fixes_applied_2026-06-17.md](doc/devel/reports/reviews/fixes_applied_2026-06-17.md) — B1 (doc gate restored), M5 (shared `decode_cram_container` for refill + caching path; folds in M2 EOF-stop + Mi3/Mi6), M3 (ContainerCache eviction/dedup unit tests), Mi1/Mi4/Mi5/Mi7/Mi8/Mi9/Mi10 + doc nits. M1/M4 reworded-as-doc per the review discussion (no `.ssr.psp` baseline → H1/H2 ULP-difference is inert; cache cap is non-output-affecting). All gates green: fmt/clippy `-D warnings`/doc clean, 1171 lib (+4) + 45 integration tests. **Open:** M3 multi-container *fetch* test (needs a multi-container CRAM fixture), Mi2 (cfg-gate bench seam), Mi11 (bench baseline guard), Mi12 (driver error context — re-assess vs `par_chunks`).
- **Long-allele delimitation fix (2026-06-24):** [ssr_delimiter_tract_aware_gap_2026-06-24.md](doc/devel/reports/implementations/ssr_delimiter_tract_aware_gap_2026-06-24.md) (impl) + [ssr_delimiter_gap_penalty_2026-06-24.md](doc/devel/reports/research/ssr_delimiter_gap_penalty_2026-06-24.md) (investigation). The delimiter applied **one uniform affine gap** (`GAP_OPEN_PROB = 2.9e-5` — HipSTR's *flank* value) across the whole read↔ref alignment, collapsing any allele ≥ ref+2 units to the reference length even with full flanks. Fix: a **tract-aware gap** (`delimit_read` picks `GAP_OPEN_PROB_TRACT = 1e-2` for repeat-tract columns, the stiff Dindel value for flanks), matching HipSTR's flank/tract split. Long alleles now extract verbatim (`CA×4…20` sweep, no collapse) and recover end-to-end (BAM→VCF `CA×10` test). Determinism/byte-identity preserved (per-column test is pure locus geometry). Companion: **long-allele window-recovery infra** (`feaa9ef`) — `flank_truncated`/`widen_region` + `WidenedSequence`/`WindowTruncated` `ReadObs` + two `.ssr.psp` QC columns (`n-widened`/`n-window-truncated`) for the mapper-misaligned all-Match case. fmt/clippy `-D warnings` clean; 1305 lib tests. **Open:** `GAP_OPEN_PROB_TRACT` is a provisional calibration constant; **Stage-2 follow-up** — confirm `ssr-call` stutter scores in-frame (unit) vs out-of-frame (non-unit) length differences (HipSTR `log_stutter_pmf`) and reconcile the tract gap rate with that slippage rate.
- **Architecture:** [ssr_pileup.md](doc/devel/architecture/ssr_pileup.md) (every structural question decided)
- **Plan:** [ssr_pileup.md](doc/devel/implementation_plans/ssr_pileup.md) (implementation sketch: build order, modules, structs, fn signatures); delimiter gap fix [ssr_delimiter_tract_aware_gap.md](doc/devel/implementation_plans/ssr_delimiter_tract_aware_gap.md) + window recovery [ssr_pileup_long_allele_window_recovery.md](doc/devel/implementation_plans/ssr_pileup_long_allele_window_recovery.md)
- **Build order:** (1) allele types in `types.rs` → (2) lift `normalize_alleles` to `src/norm_seqs/` → (3) container SSR schema → (4) stage modules (`count_repeats` → `pair_hmm` → `candidate_generation` → `triage` → `fetch_reads` → `mod`).
- **Latest code review:** [ssr_pileup_2026-06-17.md](doc/devel/reports/reviews/ssr_pileup_2026-06-17.md) (branch `ssr-pileup-review`, 11 categories) — **Request-changes**: 1 Blocker, 5 Major, 12 Minor + Nits. **B1** `cargo doc` gate broken (fetch_reads.rs:9 still links the removed `AlignmentFile` — verified failure). **M1** (the crux): the H1/H2 math rewrites (`ln_1p` fold + single-pass `ln_sum_exp3`) are ULP-different from `main` (probe: 4.74% of `ln_sum_exp2` evals differ in f64 bits) and the "byte-identical" claim is **asserted, never tested vs `main`** — the in-tree bit-identity test compares the new code against itself; run the shipped `ssr_psp_concordance` with window pinned to 10 to settle it. **M2** `CachingCramReader` `.crai` walk diverges from `refill` on `decoded==0` (continues vs stops); **M3** cache eviction/multi-container path untested (single-container gate fixture); **M4** `DEFAULT_MAX_CACHED_CONTAINERS` not in the `.ssr.psp` header (uninspectable); **M5** `.crai`-walk/decode duplicated with `refill` (the source of M2's drift). The shared-prefix DP (P1) is **provably bit-identical + well-gated**, and the `map_init` per-worker concurrency is sound. Audit trail `tmp/review_2026-06-17_ssr-pileup-review/`.
- **Latest perf review:** [perf_ssr_pileup_2026-06-16.md](doc/devel/reports/reviews/perf_ssr_pileup_2026-06-16.md) (branch `ssr-pileup-review`) — **Apply the listed wins.** First measurement of the realignment hot path: new criterion bench [benches/ssr_pileup_perf.rs](benches/ssr_pileup_perf.rs) + profiling driver [examples/profile_ssr_pileup.rs](examples/profile_ssr_pileup.rs), both driving a new `#[doc(hidden)]` [bench_harness](src/ssr/pileup/bench_harness.rs) seam. The `sample` profile pins **~74.6% of self-time on `exp`+`log`** in the pair-HMM `ln_sum_exp2/3`; cost is **perfectly linear in the `2·window+1` rung count** (13.3 ms/rung). Allocation is cold (~0.02%). Audit trail `tmp/perf_review_2026-06-16_ssr-pileup/`.
- **Impl reports:**
  - Task 1 — allele representation (`Allele`/`NormalizedSeq` + `to_sequence`/`repeat_count`) in [src/ssr/types.rs](src/ssr/types.rs): [ssr_pileup_task1_allele_types_2026-06-15.md](ia/reports/implementations/ssr_pileup_task1_allele_types_2026-06-15.md)
  - Task 2 — lift `normalize_alleles` (+ `IndexRange`) to the shared [src/norm_seqs.rs](src/norm_seqs.rs); SNP CIGAR path ([src/pileup/walker/indel_norm.rs](src/pileup/walker/indel_norm.rs)) now wraps the kernel. Behaviour-preserving; full lib suite green: [ssr_pileup_task2_norm_seqs_lift_2026-06-15.md](ia/reports/implementations/ssr_pileup_task2_norm_seqs_lift_2026-06-15.md)
  - Fast-path counter (reordered ahead of task 3) — `count_pure_tiling` core in [src/ssr/pileup/count_repeats.rs](src/ssr/pileup/count_repeats.rs); needs only `types`. Triage-typed wrapper deferred. [ssr_pileup_count_repeats_2026-06-15.md](ia/reports/implementations/ssr_pileup_count_repeats_2026-06-15.md)
  - Slow-path forward — `forward` + `HmmModel` + `PairHmmScratch` in [src/ssr/pileup/pair_hmm.rs](src/ssr/pileup/pair_hmm.rs); 3-state log-space pair-HMM, Dindel emission. Constant gap-open (homopolymer-indexed table deferred) + unbanded (banding deferred), both calibration per arch §14. [ssr_pileup_pair_hmm_2026-06-15.md](ia/reports/implementations/ssr_pileup_pair_hmm_2026-06-15.md)
  - Candidate generation Job 1 — `CandidateAllele` + `build_rungs` (on-ladder rungs `left_flank + motif×L + right_flank`) in [src/ssr/pileup/candidate_generation.rs](src/ssr/pileup/candidate_generation.rs); composes `Allele::to_sequence`. Job 2 (off-ladder normalization) deferred pending a contract decision. [ssr_pileup_candidate_generation_2026-06-15.md](ia/reports/implementations/ssr_pileup_candidate_generation_2026-06-15.md)
  - `score_candidates` in [src/ssr/pileup/pair_hmm.rs](src/ssr/pileup/pair_hmm.rs) — joins `build_rungs` + `forward` into the dense per-read `Qᵣ` (one `(allele, log-lik)` per candidate, raw scores; pruning/renorm is the aggregator's job). Small addition, context in commit; 2 added tests.
  - Candidate generation Job 2 — `build_offladder` + `normalize_offladder` (off-ladder candidates) in [src/ssr/pileup/candidate_generation.rs](src/ssr/pileup/candidate_generation.rs); contract A + B1, canonical form is the **verbatim full tract** (left-alignment is provably a no-op on a full-tract key; clean flanks). **`norm_seqs` not needed for B1** — task-2's second consumer doesn't materialize. `candidate_generation` now complete. [ssr_pileup_candidate_generation_job2_2026-06-15.md](ia/reports/implementations/ssr_pileup_candidate_generation_job2_2026-06-15.md)
- **Perf (review + applied wins, 2026-06-16):**
  - **Applied (≈−54% on realignment):** H1 fold `exp(0)` in `ln_sum_exp2`; H2 single-pass `ln_sum_exp3`; **P1 shared-prefix DP** in `score_candidates` (score the rungs' longest common prefix once, continue each tail from the saved seam). **P1 is provably bit-identical** to per-candidate `forward` (gated by `score_candidates_is_bit_identical_to_per_candidate_forward`). **H1/H2 are algebraically equivalent but not bit-identical to the prior form** (`ln_1p` is more accurate near zero; the single-pass reduction reorders the sum) — f64 scores differ within ULP, absorbed by the `f32` profile storage. Since there is **no released `.ssr.psp` baseline**, this is not an output-compat concern; the first produced catalog is the baseline (review M1, [fixes_applied_2026-06-17.md](doc/devel/reports/reviews/fixes_applied_2026-06-17.md)).
  - **Fast-path investigation (closed, 2026-06-16):** [ssr_fastpath_investigation_2026-06-16.md](doc/devel/reports/reviews/ssr_fastpath_investigation_2026-06-16.md) — measured on a real ch01 catalog + tomato CRAM (built the L9 fixture: vendored `trf-mod` + `ssr-catalog`). **Verdict: do NOT build a separate exact-match fast path** — safe predicates have ~21% recall + break byte-identity; the high-recall `observed_count` predicate would emit the wrong length on the ~11% of reads where the HMM *corrects* the pre-probe count (the realign-everything rationale, confirmed real). Instead: (a) **applied — `DEFAULT_WINDOW` lowered 10 → 6** (validated on real data: 0.48% of read-calls / 0.06% of loci change, ~40% cheaper realignment; concordance via the new `ssr_psp_concordance` ignored test); (b) on real data the **fetch/catalog-walk, not realignment, is the end-to-end wall** (confirms L5/L7) — the next lever.
  - **Open — P1 follow-up:** incremental motif-by-motif seam advance would also share the per-rung tract tail (a further ~20%), more complex; deferred.
  - **Fetch path (L7) — step 1 DONE (2026-06-17), ~17× faster, byte-identical.** Design: [ssr_pileup_read_buffer.md](doc/devel/architecture/ssr_pileup_read_buffer.md). Added a per-worker [`CachingCramReader`](src/bam/segment_reader.rs) (FIFO cache of decoded CRAM containers) + `WorkerReader` (CRAM=cache, BAM=pooled), wired into `fetch_locus_reads`/`process_locus` via `map_init`. Decodes each container once instead of once-per-locus. Gate: `caching_cram_reader_matches_per_call_path` (reader-level byte-identity) + end-to-end `.ssr.psp` diff on the real fixture = **0 record diffs**. **Measured (cap A/B + decode counter):** cap=0 (no cache) 8931 decodes / 85.3s → cap=3 341 decodes / 4.23s = **~26× fewer decodes, ~20× wall, byte-identical**. **Step 2 (slice-grouping) DROPPED — measured unnecessary:** cap=3 is already at the decode floor (cap=64 → 420, no improvement), so cross-thread boundary re-decodes are noise-level; the doc's prediction held. Decode redundancy eliminated; end-to-end factor is smaller on high-coverage runs (realignment, unchanged, dominates there).
  - **Open — L1 (experiment):** band the DP (arch §5.5 `PAIR_HMM_BAND_BP`); approximation, gate on genotype concordance. With P1 landed, the remaining `exp` (~48% self-time) is the natural target for L1 / approx-math (S1).
  - **Open — B1 (build):** add an `(aarch64, linux)` `target-cpu` floor to `.cargo/config.toml` (prod misses Neoverse SIMD; macOS bench flatters prod).
  - **Open — L9 (measurement gap):** no end-to-end bench — fetch path + driver batch loop unmeasured; build one from the in-process `stage1_fixture` (no `trf-mod`).
  - Design lever (route to SSR owners, not perf): the `window` candidate multiplier touches the realign-everything contract + Stage 2 support points.
- **Open:**
  - **Task 3 (container generalization, arch §10) deferred by choice** — a large multi-step refactor of the production `.psp` writer (58KB) + reader (136KB); gates the `.ssr.psp` writer/round-trip but not the stage's compute modules. Trait-vs-builders fork (§10.7 Q1) still to be decided at its step 1.
  - `OnLadder::to_sequence` is a clean tiling; imperfect-locus interruptions deferred to `candidate_generation.rs`. The SSR off-ladder adapter onto `norm_seqs` lands with `candidate_generation` (first consumer beyond the SNP path).
  - **Architecture revised (2026-06-15): realign-everything** — v1 realigns every spanning read (pair-HMM), dropping the CIGAR-trusting two-tier fast/slow gate; the direct-count fast path (`count_repeats`, built + parked) becomes a **measured optimization the user requires be tried**. Docs amended: arch `ssr_pileup.md` §2/§14, plan §4.
  - `triage` complete (realign-everything) — `find_longest_stretch` pre-probe + `read_footprint`/`brackets`/`extract_region`/`triage_read` in [src/ssr/pileup/triage.rs](src/ssr/pileup/triage.rs): coverage classification (footprint position + clips, no flank-byte match) → region extract → window centre. Read seam = `MappedRead`. Recovers soft-clipped long alleles via the clip-included pre-probe (test: mapper's 3 units → true 6). 22 tests. [ssr_pileup_triage_2026-06-15.md](ia/reports/implementations/ssr_pileup_triage_2026-06-15.md)
  - Per-read analysis wired — `analyze_read` + `ReadOutcome` in [src/ssr/pileup/read_analysis.rs](src/ssr/pileup/read_analysis.rs): composes `triage_read` → `build_rungs` → `score_candidates` → dense `Qᵣ` (Spanning) / Flanking / InRepeat. Reuses scratch buffers. Integration-tested incl. the realign-everything win (soft-clipped 6-unit allele the mapper called as 3 → ranks rung 6). 4 tests. Off-ladder candidate generation deferred (needs anchored observed-tract isolation; rare). Context in commit.
  - `locus_record` aggregation — `aggregate` + `SsrLocusRecord` + `QcCounts` in [src/ssr/pileup/locus_record.rs](src/ssr/pileup/locus_record.rs). **Storage model revised → all-CSR** (no histogram, no weight): every spanning read stored as one pruned + renormalized `Qᵣ` profile; `AMB_LL_DROP` pruning; renormalize provably lossless (`Z_r` cancels); derive `n_spanning`/`n_flanking`/`n_frr`, take `depth`/`n_filtered`/`mapped_reads`; drop vestigial `n_flank_indel`. In-memory (CSR flattening = deferred container's job). Diverges from spec §4.3 (amend). 6 tests. [ssr_pileup_locus_record_2026-06-15.md](ia/reports/implementations/ssr_pileup_locus_record_2026-06-15.md)
  - `fetch_reads` started — the per-locus depth cap: `Reservoir<T>` (Algorithm R) + deterministic `locus_seed(chrom, start)` (FNV-1a) + inline `SplitMix64` PRNG, in [src/ssr/pileup/fetch_reads.rs](src/ssr/pileup/fetch_reads.rs). Net-new (SNP has only a column cap); byte-identity-critical (deterministic seed + caller's fixed read order, §8.3/§8.4). `MAX_READS_PER_LOCUS=1000` placeholder. 6 tests. Context in commit.
  - **Container generalization** (the deferred §10 refactor — gates Stage-1 output): implementation sketch written, [psp_container_generalization.md](doc/devel/implementation_plans/psp_container_generalization.md) — `PspSchema` trait (or two-builders, decide at step 1), 5-step sequence (parameterize → `kind` tag → interval index → `registry_ssr`+`SsrLocusRecord` → round-trip), SNP e2e as the gate each step. Reconciles with all-CSR (SSR table = `amb_*` CSR + QC scalars, no `hist_*`). Next: step 1 (writer parameterization, SNP-only).
  - Stage modules still to build (need the container writer / real I/O): the `fetch_reads` catalog-walk + `query` driver (mirror SNP `run_pileup`: load handles once, share repo, clear per contig) + admission gate (reuse triage footprint) + bundles; the driver (`mod.rs run()`); the container schema/writer (deferred refactor — flattens `SsrLocusRecord` profiles → CSR). Single-threaded semantics first, then the fetcher-thread/pool (determinism-gated). Deferred: off-ladder wiring; the measured `count_repeats` fast-path shortcut.
  - `pair_hmm` follow-ups: homopolymer-indexed gap-open + banding (calibration/optimization, arch §14).

### Stage 2 — `ssr-call` (cohort caller: `.ssr.psp` × N → VCF)

#### Genotyping + parameter pre-pass (Phases 2/3 — fused plan; Milestone A done)
- **Status:** **Milestone F shipped — ALL PLAN MILESTONES A–F COMPLETE** (2026-06-23, branch `ssr-cohort`). The algorithmic SSR cohort caller is implemented, reviewed, and fixed end-to-end: pre-pass → genotyping → FP control → VCF, parallel + byte-identical across thread counts. Remaining = empirical real-data calibration of provisional constants + documented refinements + driver wiring (TSV dump → full VCF). Checkpoints 1+2 signed off.
- **Milestone F review + fix:** [review](doc/devel/reports/reviews/ssr_call_genotyping_milestone_f_2026-06-23.md) (Approve-with-changes: 0 Blocker/Major; concurrency review clean — no `unsafe`/locks/atomics, determinism structural) + [fixes_applied_2026-06-23_v7.md](doc/devel/reports/reviews/fixes_applied_2026-06-23_v7.md) — Mi1 determinism-layer doc Applied. 1253 lib tests.
- **Milestone F impl report:** [ssr_call_genotyping_milestone_f_2026-06-23.md](doc/devel/reports/implementations/ssr_call_genotyping_milestone_f_2026-06-23.md) — F1: `run_prepass_stats` (rayon fold/`PrepassStats::merge`, integer counts) + `run_cohort_em` per-locus EM (`par_iter().zip().collect()`, order-preserving) are byte-identical at 1 vs 4 threads. F2: calibration constants provisional (real-data follow-up); correctness gates (checkpoints 1/2, M3, byte-identity, end-to-end) tested. 1253 lib tests (+2).
- **Milestone E review + fixes:** [review](doc/devel/reports/reviews/ssr_call_genotyping_milestone_e_2026-06-23.md) (Approve-with-changes: 0 Blocker/Major) + [fixes_applied_2026-06-23_v6.md](doc/devel/reports/reviews/fixes_applied_2026-06-23_v6.md) — all 5 Applied (Mi1 `Δlevel` convergence; Mi2 `add_bin` dedup; group-index assert; **end-to-end pipeline integration test** prepass→cluster→F-loop→FP→VCF). 1251 lib tests.
- **Milestone E impl report:** [ssr_call_genotyping_milestone_e_2026-06-23.md](doc/devel/reports/implementations/ssr_call_genotyping_milestone_e_2026-06-23.md) — E1 `inbreeding.rs` (`run_cohort_em`: per-locus EM with per-sample `F` + per-group level; `reduce_f` excess-homozygosity via `FixedPointAccum`, shrink+clamp ≤0.99; `reduce_level` hard-attribution refit). E2 `vcf_out.rs` (allele-balance FP defence → no-call depth-inflated false hets; emit-iff-variable; site QUAL proxy; `F_IS` warning). `em::run_locus_em_with` + `LocusCall.posterior_hom`. 1250 lib tests (+11).
- **Milestone D (D3) review + fixes:** [review](doc/devel/reports/reviews/ssr_call_genotyping_milestone_d_d3_2026-06-23.md) (Approve-with-changes: 0 Blocker/Major) + [fixes_applied_2026-06-23_v5.md](doc/devel/reports/reviews/fixes_applied_2026-06-23_v5.md) — all 4 Applied (Mi1 reference-depth → ClusterCfg; Mi2/nit docs; clustering determinism test). 1242 lib tests.
- **Milestone D (D3) impl report:** [ssr_call_genotyping_milestone_d_d3_2026-06-23.md](doc/devel/reports/implementations/ssr_call_genotyping_milestone_d_d3_2026-06-23.md) — `sample_groups.rs`: deterministic union-find clustering of per-sample `(ε,level)` → sample groups; per-`(group,period)` shape shrunk to `θ_period` (M3: two-protocol sim recovers divergent decays, single collapses); binomial-BIC `ε`-freeze check. D1 extended with per-sample slip profiles. 1241 lib tests (+3).
- **Milestone D (D1+D2) review:** [ssr_call_genotyping_milestone_d_d1d2_2026-06-23.md](doc/devel/reports/reviews/ssr_call_genotyping_milestone_d_d1d2_2026-06-23.md) — Approve-with-changes: 0 Blocker, 0 Major, 2 Minor + Nits + 2 missing tests.
- **Milestone D (D1+D2) fixes-applied:** [fixes_applied_2026-06-23_v4.md](doc/devel/reports/reviews/fixes_applied_2026-06-23_v4.md) — all 6 Applied (Mi1/Mi2/nit docs; `Rungs` import; het-contribution + determinism tests). 1238 lib tests.
- **Milestone D (D1+D2) impl report:** [ssr_call_genotyping_milestone_d_d1d2_2026-06-23.md](doc/devel/reports/implementations/ssr_call_genotyping_milestone_d_d1d2_2026-06-23.md) — `prepass.rs`: D1 `accumulate_locus` (confident-genotype slip + ε stats off the B1 gate), D2 `estimate` (ε from base-mismatch fraction, per-period shape via direction-split + decay MLE, per-sample level line). 1236 lib tests (+4). D3 (clustering + per-(group,period) shape + ε-freeze, the M3 milestone) is post-gate.
- **Milestone C review:** [ssr_call_genotyping_milestone_c_2026-06-23.md](doc/devel/reports/reviews/ssr_call_genotyping_milestone_c_2026-06-23.md) — Approve-with-changes: 0 Blocker, 0 Major, 2 Minor + Nit + 1 missing test.
- **Milestone C fixes-applied:** [fixes_applied_2026-06-23_v3.md](doc/devel/reports/reviews/fixes_applied_2026-06-23_v3.md) — Mi1 (ploidy-panic doc, kept per no-silent-fallback), Mi2 (first-of-length doc), MT-1 (moderate-stutter EM test) applied; the perf Nit deferred to F1. 1232 lib tests.
- **Milestone C impl report:** [ssr_call_genotyping_milestone_c_2026-06-23.md](doc/devel/reports/implementations/ssr_call_genotyping_milestone_c_2026-06-23.md) — C1 candidate assembly (`19b1c61`), C2 `likelihood` (`read_likelihood` Σ_Δ·Σ_v + genotype mix/λ), C3 `allele_freq_prior` (`G₀`) + `em_init` (π⁰/θ⁰ seeds), C4 `em` (slim diploid SSR EM — **Q-G2 decided: SSR-specific, not a posterior_engine graft**) + `vcf_out` (minimal `GT:GQ:REPCN`). 1231 lib tests (+24).
- **Milestone B review:** [ssr_call_genotyping_milestone_b_2026-06-23.md](doc/devel/reports/reviews/ssr_call_genotyping_milestone_b_2026-06-23.md) — Approve-with-changes: 0 Blocker, 0 Major, 3 Minor + Nits + 1 missing test.
- **Milestone B fixes-applied:** [fixes_applied_2026-06-23_v2.md](doc/devel/reports/reviews/fixes_applied_2026-06-23_v2.md) — all 6 Applied (Mi1 overflow guard; Mi2/Mi3/nits doc; pure-contraction test). 1207 lib tests.
- **Milestone B impl report:** [ssr_call_genotyping_milestone_b_2026-06-23.md](doc/devel/reports/implementations/ssr_call_genotyping_milestone_b_2026-06-23.md) — B1 `rung_ladder` (`build_rungs` + heuristic `resolve_confident_genotype`), B2 `stutter` (`s_theta`/`reach_variants`/`refine_theta_locus`, the scoring counterpart of the sim forward model), B3 `pair_hmm` (`align_subst` banded forward, substitutions-in-tract/gaps-in-flank). 1206 lib tests (+22).
- **Latest review:** [ssr_call_genotyping_milestone_a_2026-06-23.md](doc/devel/reports/reviews/ssr_call_genotyping_milestone_a_2026-06-23.md) — **Approve-with-changes**: 0 Blocker, 0 Major, 2 Minor + Nits.
- **Latest fixes-applied:** [fixes_applied_2026-06-23.md](doc/devel/reports/reviews/fixes_applied_2026-06-23.md) — all 4 findings Applied (Mi1 `FixedPointAccum` non-finite debug-assert + magnitude doc + guard test; Mi2 separated-het + per-group-shape sim tests; `below`→`index_below`; drop alloc-for-length). 1184 lib tests.
- **Plan:** [ssr_call_genotyping_and_parameters.md](doc/devel/implementation_plans/ssr_call_genotyping_and_parameters.md) — the fused Phase-2/3 plan (A1→F2), with the per-milestone execution loop (implement → commit → review → commit → fix → commit) and the two human checkpoints (after C4, after D2).
- **Architecture:** [parameters](doc/devel/architecture/ssr_call_parameters.md) (pre-pass) + [genotyping](doc/devel/architecture/ssr_call_genotyping.md) (EM/VCF); spec [ssr_cohort_mark2.md](doc/devel/specs/ssr_cohort_mark2.md) §4.2–§4.5.
- **Impl report (Milestone A):** [ssr_call_genotyping_milestone_a_2026-06-23.md](doc/devel/reports/implementations/ssr_call_genotyping_milestone_a_2026-06-23.md).
- **Code:** [src/ssr/cohort/](src/ssr/cohort/) — `param_estimation.rs` (`ParamSet`/chemistry types/`SlipProfile`/`SampleStutterStats`/`FixedPointAccum`/`MAX_SLIP`), `candidate_set.rs` (`CandidateSet`/`Admission`), `rung_ladder.rs` (`Resolution`/`ResolvedGenotype`/`PeakAllele`/`UnresolvedReason`), `stutter.rs` (`PlacementVariant`), `sim.rs` (`#[cfg(test)]` simulator + `TruthTable`, SplitMix64). 1181 lib tests (+16).
- **Open:**
  - **Milestone A review** (loop step 3) — pending.
  - Provisional calibration constants to pin in F2: `MAX_SLIP=10`, the simulator's geometric forward-model parameterization, the substitution model.
  - B2's `align`/kernel must be written to **match `sim.rs`'s documented forward model** (else D-milestone recovery tests catch the divergence).
- **Pre-existing unrelated bug surfaced:** `cargo test --all-targets` trips `benches/psp_writer_perf.rs:386` (index OOB in the bench harness) — separate fix.

#### Driver wiring (genotyper → VCF; two-pass streaming)
- **Status:** fixes-applied (2026-06-24, branch `ssr-cohort`) — the genotyper is wired into the `ssr-call` driver and emits a real VCF (plan G→H→I→J complete); the fresh-eyes multi-agent re-review returned **Request-changes**, and across two fix runs the **Blocker + all 6 Majors + 15 of 16 Minors are now Applied** (only Mi3/Mi4, both bench-gated, remain open — no SSR bench exists yet).
- **First full BAM→VCF path:** the `ssr-call` CLI wrapper ([src/pop_var_caller/ssr_call.rs](src/pop_var_caller/ssr_call.rs)) was refreshed off the stale "TSV dump / EM not built" text and its `--queue-depth` default fixed to `0` (the driver's auto-chunk sentinel; `4` had serialized the sweep). A committed lib integration test ([src/ssr/end_to_end_tests.rs](src/ssr/end_to_end_tests.rs)) now drives the **whole chain** — synthetic reference + per-sample BAMs + a written catalog (no `trf-mod`) → real `ssr-pileup` ×N → real `ssr-call` → VCF — asserting a PASS length-polymorphism (hom-ref 8/8, het 6/8, hom-alt 6/6) + a dropped monomorphic locus, guarding the catalog-md5 / chrom-table / coordinate contracts across stages. (Note surfaced by the fixture: a `+4 bp` tract insertion, `CA×10` vs a `CA×8` reference, collapsed to `CA×8` in the pileup realignment while `CA×6` recovered cleanly — a realignment-band detail to revisit in calibration, not a chain defect.)
- **Latest fixes-applied (v2 — deferred-Minor follow-up):** [fixes_applied_2026-06-24_ssr_call_driver_v2.md](doc/devel/reports/reviews/fixes_applied_2026-06-24_ssr_call_driver_v2.md) — picked up the 8 v1-deferred Minors: **6 Applied** (Mi1, Mi2, Mi5, Mi13, Mi14, Mi16), **2 Deferred** (Mi3/Mi4 — bench-gated, user decision). **Mi1**: new `attribution::nearest_parent` primitive shared by `attribute_locus`/`accumulate_locus`/`allele_balance` (one tie-break; byte-identity preserved). **Mi2**: `LocusModel` bundle drops `compute_data_ll` to 7 args (its `#[allow]` removed). **Mi5**: `LengthBin`/`AlleleCopies` structs replace the tuple bins. **Mi13**: new `SsrMergeError::InvalidAlleleByte` rejects non-ACGTN allele bytes at the merge boundary (+ `from_utf8_lossy`→`from_utf8().expect()` in `format_vcf_record`). **Mi14**: emitted records always print numeric QUAL (`0.0`, not `.`). **Mi16**: move instead of clone `level_per_group`. Gates: fmt/clippy `-D warnings`/doc clean; **1292 lib tests** (+5), only `--all-targets` failure is the pre-existing `psp_writer_perf` bench panic.
- **Prior fixes-applied (v1):** [fixes_applied_2026-06-24_ssr_call_driver.md](doc/devel/reports/reviews/fixes_applied_2026-06-24_ssr_call_driver.md) — 15 Applied (B1; M1–M6; Mi6–Mi12, Mi15), 8 Deferred. **B1**: VCF rows dense over the cohort (`format_vcf_record` takes `n_samples`, places present calls by `locus.present`, `./.:.:.` for absent) + partial-coverage test. **M1/M2**: silent `"?"` contig + `sample_chemistry` triple `unwrap_or` → loud panics. **M3**: shared `DEFAULT_G0_FALLBACK_P`. **M4**: `partial_cmp().unwrap()` → `total_cmp`. **M5**: `InvalidVcfName` typed error (tab/newline/`,<>`). **M6**: filtered-locus emit test. Each fix its own commit (`1999d82`→`f3f8d04`).
- **Prior review:** [ssr_call_driver_2026-06-24.md](doc/devel/reports/reviews/ssr_call_driver_2026-06-24.md) — **Request-changes**: 1 Blocker, 6 Major, 16 Minor + Nits (audit trail `tmp/review_2026-06-24_ssr-call-driver/`). **B1** (High): `format_vcf_record` writes per-sample columns in **present-order, not cohort-order** — any partial-coverage locus (the normal real-cohort case; cursor yields Absent + merger sparse-omits only all-absent loci) gets a short, mis-aligned data row (wrong genotypes + invalid VCF); masked by tests whose `write_cohort` gives every sample evidence at every locus. Majors: silent `"?"` contig + `sample_chemistry` triple `unwrap_or` defaults + `period_decay` duplicate `G₀` fallback (all contradict the no-silent-default invariant the module upholds loudly elsewhere); uncommented `partial_cmp().unwrap()` NaN-argmax on the parallel emit path (×2, convergent across 4 categories); no tab/control-char guard on untrusted contig/sample names; untested filtered-locus *emit* branch. Determinism/byte-identity, decision-E hard error, ploidy panic, emit/drop policy all verified sound. fmt/clippy/doc clean; 1280 lib tests (only `--all-targets` failure is the pre-existing `psp_writer_perf` bench panic).
- **Architecture:** [ssr_call_driver.md](doc/devel/architecture/ssr_call_driver.md) — settled (decisions A–E): **two-pass streaming** (materialization rejected) — bounded pre-pass/burn-in freezes the cross-locus-pooled params (`F`, group level line, `θ_period`, ε, `G₀`), then a single streaming sweep of independent per-locus EMs; per-locus stutter refined locally with shrinkage. Fit `G₀` in the pre-pass; dedicated SSR VCF header (contigs from `.ssr.psp` headers, not the catalog); hard error on samples with no confident genotype; ploidy 2.
- **Plan:** [ssr_call_driver.md](doc/devel/implementation_plans/ssr_call_driver.md) — milestones G (G₀ fit) → H (Step 1 streaming driver + VCF = task DoD) → I (Step 2 per-locus stutter adaptation) → J (parallelism), each through the implement → review → fix loop.
- **Step G1 (G₀ decay fit):** implemented `640a1e6`, [review](doc/devel/reports/reviews/ssr_call_g0_fit_2026-06-23.md) (Approve-with-changes: 0 Blocker/Major, 3 Minor) + [fixes](doc/devel/reports/reviews/fixes_applied_2026-06-23_g0_fit.md) (Mi2/Mi3 + 2 Nits Applied; Mi1 → H1). Per-period `G₀` `p` fit from the confident germline allele spread over variable loci (closed-form geometric MLE, thin-period fallback, over-tightening clamp).
- **Step H1 (`build_param_set`):** implemented `4168331`, [review](doc/devel/reports/reviews/ssr_call_build_param_set_2026-06-23.md) (Approve-with-changes: 0 Blocker/Major, 2 Minor) + [fixes](doc/devel/reports/reviews/fixes_applied_2026-06-23_build_param_set.md) (Mi2 + 2 Nits Applied; Mi1 → H4). Pre-pass → frozen `ParamSet`; decision-E hard error (`SsrCallError::UnresolvedSamples`); resolves G1-Mi1 (backfills `G0FitCfg.fallback_p` for characterized periods without a fit).
- **Step H2 (merger accessors):** implemented `f01c14e`, [review](doc/devel/reports/reviews/ssr_call_merger_accessors_2026-06-23.md) (Approve-with-changes: 0 Blocker/Major, 1 Minor) + [fixes](doc/devel/reports/reviews/fixes_applied_2026-06-23_merger_accessors.md) (md5-assert Nit Applied; Mi1 → H3/H4). `CohortMerger::chromosomes()` (name+length+md5 from the `.ssr.psp` headers, not the catalog) + `sample_names()` (basename); two-pass re-open proven.
- **Step H3 (VCF header writer):** implemented `2ff5b1b`, [review](doc/devel/reports/reviews/ssr_call_vcf_header_2026-06-23.md) (Approve-with-changes: 0 Blocker/Major/Minor, Nits) + [fixes](doc/devel/reports/reviews/fixes_applied_2026-06-23_vcf_header.md) (input-contract + md5-omission doc notes Applied). `vcf_out::write_vcf_header` — dedicated plain-text SSR header (fileformat, ##contig w/ lengths, PERIOD INFO, GT/GQ/REPCN FORMAT, SSR FILTER from `filter_text`, `#CHROM … <samples>`); warnings as `##ssrCallWarning=`.
- **Step H4 (streaming driver → VCF) — TASK DoD MET:** implemented `d17f524`, [review](doc/devel/reports/reviews/ssr_call_streaming_driver_2026-06-23.md) (Approve-with-changes: 0 Blocker/Major, 3 Minor) + [fixes](doc/devel/reports/reviews/fixes_applied_2026-06-23_streaming_driver.md) (Mi2 cross-thread-determinism test + Mi1/Mi3 doc notes Applied; filtered-record e2e + representative subset deferred to calibration). `driver::run` is now the two-pass streaming pipeline: burn-in (bounded subset → freeze chemistry+`F`+level on a `config.threads` pool) → genotyping sweep (re-open, stream, `run_locus_em_with(frozen)`, emit policy) → VCF. Sample-name uniqueness (H2-Mi1) enforced; H1-Mi1 boundary accepted; ploidy 2 (D). Integration test drives real `run()` over an on-disk cohort → PASS variant + monomorphic drop + decision-E hard error; VCF byte-identical across thread counts. TSV `write_dump` path removed. 1271 lib tests; fmt/clippy clean.
- **Step I1 (per-locus `θ_locus` shape refit):** implemented `46df902`, [review](doc/devel/reports/reviews/ssr_call_theta_locus_2026-06-24.md) (Approve-with-changes: 0 Blocker/Major, 2 Minor) + [fixes](doc/devel/reports/reviews/fixes_applied_2026-06-24_theta_locus.md) (Mi1 shared `add_slip` + pin test Applied; Mi2 π warm-start deferred). The genotype EM now adapts the stutter **shape** per locus: `run_locus_em_with` wraps the π-EM in a `θ_locus` loop (genotype under the frozen `θ_period` seed → attribute reads to called alleles → `refine_theta_locus` shrunk toward the seed → re-genotype until settled), local so the sweep stays single-pass + byte-identical. `EmCfg` gains `theta_max_rounds`/`theta_shrink`/`theta_tol`. 1276 lib tests; fmt/clippy clean.
- **Step I2 (per-locus stutter-rate refit):** implemented `fd53330`, [review](doc/devel/reports/reviews/ssr_call_level_refit_2026-06-24.md) (Approve-with-changes: 0 Blocker/Major, 2 Minor) + [fixes](doc/devel/reports/reviews/fixes_applied_2026-06-24_level_refit.md) (Mi1 rename `refit_max_rounds` + Mi2 hierarchy doc Applied). The per-locus refit loop now adapts the stutter **rate** too: `attribute_locus` → `LocusSlipFit` (slip profile + observed/expected slips); `refit_level_multiplier` = `(slipped+strength)/(expected+strength)` shrunk toward 1, applied as a multiplier on the frozen group level (group keeps the length-dependence, locus nudges the rate). `EmCfg` gains `level_shrink`/`level_tol`; `theta_max_rounds`→`refit_max_rounds`. 1278 lib tests; fmt/clippy clean. **Milestone I (Step 2 per-locus stutter adaptation) complete.**
- **Step J (chunk-parallel sweep):** implemented `8d7e4d2`, [review](doc/devel/reports/reviews/ssr_call_parallel_sweep_2026-06-24.md) (Approve-with-changes: 0 Blocker/Major, 2 Minor) + [fixes](doc/devel/reports/reviews/fixes_applied_2026-06-24_parallel_sweep.md) (Mi1 `DEFAULT_SWEEP_CHUNK` when `queue_depth=0` + test; Mi2 arch-doc realization note Applied). The Pass-2 genotyping sweep is now chunk-parallel: bounded chunks genotyped on the `--threads` pool via order-preserving `par_iter`, written in catalog order — byte-identical across thread counts (multi-chunk + 1-locus determinism tests), no `seq`-reorder needed (chosen realization over the full channel pipeline; arch §4 J-realization note). `config.queue_depth` = chunk size. 1280 lib tests; fmt/clippy clean. **Milestone J complete → the `ssr-call` plan (G→H→I→J) is DONE.**
- **Resolved by the 2026-06-24 v1 fix run (Applied):** B1 (dense VCF rows + test); M1–M6 (loud `"?"`/`sample_chemistry`, shared `DEFAULT_G0_FALLBACK_P`, `total_cmp`, `InvalidVcfName` validation, filtered-locus emit test); Mi6 (justify the `#[allow]`s), Mi7 (`LocusSlipFit.profile` literal), Mi8 (`level_multiplier`), Mi9 (`FrozenParams.chemistry`), Mi10 (stale comment), Mi11 (`inbreeding_f` doc), Mi12 (coercion docs), Mi15 (refit non-convergence test).
- **Resolved by the 2026-06-24 v2 fix run (Applied):** Mi1 (shared `attribution::nearest_parent`), Mi2 (`LocusModel` bundle — `compute_data_ll` `#[allow]` removed), Mi5 (`LengthBin`/`AlleleCopies` structs), Mi13 (`InvalidAlleleByte` non-ACGTN rejection at the merge boundary + loud `from_utf8` in `format_vcf_record`), Mi14 (always-numeric QUAL), Mi16 (move vs clone).
- **Open (deferred — Minors):**
  - **Mi3 / Mi4** — per-round (`compute_data_ll`) / per-locus (`f_present`) allocation levers; **bench-gated and no SSR bench exists yet**, so deferred by explicit decision until one does (then measure the wall cost).
  - Plus the remaining §8 missing-test specs (all-monomorphic header-only; empty-`psp_files` through `run`; `add_slip` cap boundary; `phred_gq` cap/floor; refitting-locus cross-thread byte-identity).
  - Open questions still owed (review §4): **Q1** — does Stage-1 `ssr-pileup` emit sparse per-sample `.ssr.psp` files (sets how soon partial-coverage loci appear; B1 is fixed regardless).
  - Carry-overs (accepted, not re-litigated): H4-Mi1 representative burn-in subset (positional first-cap biases chemistry + couples to decision-E); H1-Mi1 `G₀` backfill universe; I1-Mi2 π warm-start; the soft per-read responsibility split; the fully-overlapping channel pipeline (vs chunking); all `dev_default` constant *values* (real-data calibration is the separate downstream effort).

#### Reading & merge layer (Phases 0–3 done; `ssr-call` runnable; review applied)
- **Status:** fixes-applied (2026-06-21, branch `ssr-cohort`). Phases 0 (scaffolding) + 1 (cursor) + 2 (merger) + 3 (driver — single-threaded, `ssr-call` runs end-to-end → catalog-ordered TSV dump). Two-pass / prefetch-pool + the genotyping EM/VCF to follow.
- **Latest review:** [ssr_call_reading_2026-06-21.md](doc/devel/reports/reviews/ssr_call_reading_2026-06-21.md) — **Request-changes**: 1 Blocker, 5 Major, 8 Minor + Nits (audit trail `tmp/review_2026-06-21_ssr-call-reading/`).
- **Fixes applied:** [fixes_applied_2026-06-21.md](doc/devel/reports/reviews/fixes_applied_2026-06-21.md) — B1 + all 5 Majors + 7 of 8 Minors + nits (commit `c15280c`; 1165 lib tests, +6). **Open questions resolved:** same sorted catalog for all samples is a caller contract enforced with **hard errors** (B1 `LocusNotInCatalog`, M1 `UnsortedCatalog`); `.ssr.psp` treated as untrusted → coordinate validation at the **decode boundary** (M2). **Deferred:** Mi4 (tree-wide `allow(dead_code)` — tracked temporary) + cosmetic nits.
- **Spec:** [ssr_cohort_mark2.md §4.1](doc/devel/specs/ssr_cohort_mark2.md) (reading & orchestration intent, settled 2026-06-19).
- **Architecture (settled):** [ssr_call_reading.md](doc/devel/architecture/ssr_call_reading.md) — `SampleEvidenceCursor` (`held` + `last_query` monotonic guard, `evidence_at`), catalog-driven k-way merge → one `CohortLocus` at a time, shared decode-priority pool + prefetched futures (profiling-gated), two-pass re-read. Companions (drafts): [parameters](doc/devel/architecture/ssr_call_parameters.md), [genotyping](doc/devel/architecture/ssr_call_genotyping.md).
- **Plan:** [ssr_call_reading.md](doc/devel/implementation_plans/ssr_call_reading.md) — 6 incremental phases (0 scaffolding → 1 cursor → 2 merger → 3 driver/stub → 4 two-pass re-read → 5 prefetch pool).
- **Impl report (Phases 0–1):** [ssr_call_reading_phase1_2026-06-21.md](doc/devel/reports/implementations/ssr_call_reading_phase1_2026-06-21.md).
- **Code:** [src/ssr/cohort/](src/ssr/cohort/) — `types.rs` (`LocusId`/`SsrQc`/`SampleEvidence`/sparse-SoA `CohortLocus`), `reader.rs` (`SampleEvidenceCursor`: `held`+`last_query` contract, coordinate inversion, `observed→seq_counts`), `merge.rs` (`CohortMerger`: catalog-driven merge → `(seq, CohortLocus)`, `open`/`from_parts` validation: same-catalog md5 + chrom-id reconciliation), `driver.rs` (`run`/`write_dump`/`format_locus` — single-threaded, catalog-ordered TSV dump), `test_support.rs` (shared fixtures). Enabler: [src/psp/reader.rs](src/psp/reader.rs) `OwnedRecordsIter` + `PspReader::into_records_of` (SNP path untouched). `run_ssr_call` wired [src/pop_var_caller/ssr_call.rs](src/pop_var_caller/ssr_call.rs). 35 tests; 1159 lib pass.
- **Open:**
  - **Mi4 (deferred)** — narrow the tree-wide `#![allow(dead_code)]` on `src/ssr/mod.rs` once the genotyping EM consumes the cohort surface.
  - **Phase 4** — two-pass re-read (the parameter pre-pass consumes the merge stream, then genotyping re-reads it; merge must be cheap to restart).
  - **Phase 5** — profiling-gated prefetch pool (Q-R4↔Q-R6).
  - **Genotyping EM + VCF** — the real worker + output (separate plan: arch [parameters](doc/devel/architecture/ssr_call_parameters.md) / [genotyping](doc/devel/architecture/ssr_call_genotyping.md)); the driver's TSV dump is a placeholder until then.
  - Q-R3 (queue depth default, measured), Q-R6 (genomically-aligned Stage-1 blocks — Stage-1-writer follow-up).

#### Interrupted-repeat recall — sequence-keyed alleles (Phase 1) + per-allele stutter (Phase 2)
- **Status:** reviewed (2026-07-06, branch `ssr-interruptions`). Phase 1 (recall) shipped + validated; Phase 2 (purity→level) shipped guarded + inert-on-tomato (owner-accepted). 11 commits `550b691`…`4b86518`.
- **What it does:** makes `ssr-call` candidate alleles **sequence-keyed** so two same-length tracts differing only by an interior interruption (pure `(TA)7` vs interrupted `(TA)3C(TA)3`) become distinct candidates and get called — recovering length-monomorphic-but-composition-polymorphic loci the caller silently dropped. Converts the four length-keyed holdouts (nomination + §5.2 admission bar, EM seed, θ_locus slip attribution, allele-balance FP term) to key on the candidate sequence. Phase 2 adds a measurement-gated purity→stutter-level covariate.
- **Spec / arch / plan:** [specs/ssr_interrupted_repeat_recall.md](doc/devel/specs/ssr_interrupted_repeat_recall.md), [architecture/ssr_interrupted_repeat_recall.md](doc/devel/architecture/ssr_interrupted_repeat_recall.md), [implementation_plans/ssr_interrupted_repeat_recall.md](doc/devel/implementation_plans/ssr_interrupted_repeat_recall.md); conforms to [ssr_cohort_mark2.md](doc/devel/specs/ssr_cohort_mark2.md) §5/§6/§7 (§6 allele-balance amendment).
- **Code:** [src/ssr/cohort/](src/ssr/cohort/) — `rung_ladder.rs` (`RungSeq` distinct-sample tally), `candidate_set.rs` (`cohort_alleles` + §5.2 `CandidateCfg` bar), `em_init.rs` (`candidate_for_sequence`), `attribution.rs` (`nearest_called_by_sequence` + `allele_responsibilities`), `em.rs` (`SampleCall::allele_support`, `fill_allele_support`, `candidate_level`), `vcf_out.rs` (sequence-aware `allele_balance`), `param_estimation.rs` (`PurityLevel` + `interruption_count`), `prepass.rs` (`fit_purity_level`), plus the P2.0 diagnostic in `driver.rs` + `examples/ssr_slip_dump.rs`.
- **Validation reports:** [ssr_interrupted_repeat_p1_validation_2026-07-06.md](doc/devel/reports/ssr_interrupted_repeat_p1_validation_2026-07-06.md) (Phase-1: +166 PASS recovered, 116 HipSTR seq-only interruptions, het recovery 26→54%, no regression, paralog-clean), [ssr_interrupted_repeat_p20_slip_measurement_2026-07-06.md](doc/devel/reports/ssr_interrupted_repeat_p20_slip_measurement_2026-07-06.md) (P2.0 gate), [ssr_interrupted_repeat_p2_validation_2026-07-06.md](doc/devel/reports/ssr_interrupted_repeat_p2_validation_2026-07-06.md) (P2.3: inert-on-tomato + population discrepancy).
- **Latest review:** [ssr_interrupted_repeat_recall_2026-07-06.md](doc/devel/reports/reviews/ssr_interrupted_repeat_recall_2026-07-06.md) — Approve-with-changes: 0 Blocker, 1 Major (M1: `nearest_called_by_sequence` unconditional `align_subst` — hot-path DP + equidistant slip-sign flip), 11 Minor (mostly doc/test/tuning), Nits. Audit trail `tmp/review_2026-07-06_ssr-interrupted-repeat/`.
- **Open (from review + tracked follow-ups):**
  - **M1** — short-circuit `nearest_called_by_sequence` for the single-nearest case + scope the composition tie-break to same-length ties (perf + slip-sign correctness).
  - **Mi1** — `interruption_count` over-scores a first-unit interruption; guard the phase precondition.
  - **Mi2** — `allele_balance` silently returns `1.0` when `allele_support` is empty (untested; add `debug_assert!` + test).
  - Missing tests (review §8): the equidistant-composition tie, `fit_purity_level` floor + multi-`k`, Phase-2 cross-thread determinism, `sample_eps` panic, junk-read skip.
  - Tracked in `doc/devel/TODO.txt`: same-length paralog FP guard; D1 BIC confident-genotype gate (the low-depth/high-`F_IS` het-undercall tail); §5.2 threshold sweep; Phase-2 refit from the broad caller genotypes (the P2.0 population, where the effect is clean).

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
  - **`.psp` → cohort VCF arm — perf review against real data
    2026-05-20:**
    [perf_psp_to_vcf_2026-05-20.md](doc/devel/reports/reviews/perf_psp_to_vcf_2026-05-20.md).
    Used the bench above + an `examples/profile_cohort_e2e.rs`
    one-off + `perf record` on real tomato (SL4.0)
    `SRR7279725.small.psp × N=10` to produce the headline
    diagnosis: pipeline is essentially single-threaded
    (T=1=12.7s ≈ T=16=13.6s); DUST 33 %, allocations 21 %, PSP
    decode 15 %; per-group merger + posterior together <2 %. Seven
    Hot-path + 12 Likely findings; **H1 per-chromosome parallelism**
    is the order-of-magnitude lever.
  - **`.psp` → cohort VCF arm — H1 per-chromosome parallelism shipped
    2026-05-20:**
    [cohort_per_chromosome_parallel_2026-05-20.md](doc/devel/reports/implementations/cohort_per_chromosome_parallel_2026-05-20.md).
    Realised **3.85× wall-time reduction at T=13** on the multi-chrom
    real-data fixture
    (`tmp/SRR7279727.multichrom.psp` — 2 Mbp from each of 13 SL4.0
    chroms via `samtools view --regions`, N=10 cohort: 106.6 s →
    27.7 s). Workload imbalance (ch00 unplaced/decoy reads at
    13× the median per-chrom count) gates the ceiling below the
    plan's predicted 6–10×; L5 contention is the next ceiling
    (acknowledged for follow-up). Includes L1 (per-group inner
    `par_iter` removed) and a new pure-Rust bgzf-aware concat
    module (`src/var_calling/vcf_writer/concat.rs`).
  - **CRAM → VCF arm (`var-calling-from-bam`) — REMOVED 2026-06-01.**
    The direct single-sample BAM/CRAM → VCF subcommand was deleted
    entirely (see the Current-focus "Last completed task"); the only
    route to a VCF is now `pileup` → `.psp` → `var-calling`. The
    history below is retained for the record. Per-chromosome
    parallelism shipped 2026-05-24:
    [var_calling_from_bam_per_chromosome_2026-05-24.md](doc/devel/reports/implementations/var_calling_from_bam_per_chromosome_2026-05-24.md);
    plan
    [var_calling_from_bam_per_chromosome.md](doc/devel/implementation_plans/var_calling_from_bam_per_chromosome.md).
    Four-commit PR (`29b28e8` → `e5261ba` → `050da41` → `a858832`):
    new `src/bam/index_preflight.rs` + `--build-map-file-index`
    flag (opt-in `.crai` auto-build; off by default with a
    samtools-pointing error), new
    `CramMergedReader::query` indexed per-contig variant +
    `OwnedIndexedCramRecords` iterator, new
    `process_one_chromosome_from_bam` per-chrom worker, and a
    `run_var_calling_from_bam` reshape that retires the serial
    `run_cohort_pipeline_for_single_sample` + `PerChromRecordsIter`
    helpers in favour of `rayon::par_iter` over chromosomes +
    `vcf::concat::concat_fragments` (no new file-format module —
    reuses the cohort H1's pure-Rust bgzf concat). Net diff
    +1320 / -520 across the four commits; 898 lib + every
    integration test pass; clippy + fmt clean. Wall-time
    validation on real multi-chrom tomato CRAMs (the analogue of
    cohort H1's 3.85× at T=13) is deferred: the
    `examples/profile_from_bam_e2e.rs` + `benches/from_bam_e2e_perf.rs`
    infrastructure (commit 5 of the plan) was scoped out.
  - **CRAM → `.psp` arm (`pileup`): not planned** (design call
    2026-05-24). The typical PSP workflow runs N independent samples
    and the orchestrator (Snakemake, Nextflow, GNU parallel)
    already provides the per-sample parallelism; a single `pileup`
    invocation processes exactly one sample and stays serial.
    Revisit only if a real workload appears where one sample's
    `.psp` build dominates a pipeline's wall time.

