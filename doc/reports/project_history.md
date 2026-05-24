# Project history — `join_per_sample_vcfs` → `pop_var_caller`

**Span:** 2025-12-04 → 2026-05-24 (~5.5 months, 550 commits)
**Source growth:** 10 → 96 Rust files in `src/`, ~1.7 k → ~59.8 k LOC.
**Authoring model:** single primary author working with a Claude Code
assistant; commit messages and per-stage `doc/devel/reports/` and
`doc/devel/specs/` artefacts reflect that workflow.

The history breaks cleanly into **four phases**. The pivot at the
boundary of Phase 2 and Phase 3 is the key inflection point: the
project's original goal (joining per-sample gVCFs) was achieved and
then deliberately discarded in favour of a much larger goal (a
complete cohort variant caller from BAM).

## Phase 1 — gVCF parsing foundation (2025-12 → late Feb 2026)

**Objective:** build a fast, custom parser for per-sample gVCFs that
other tooling can iterate.

The repo sits dormant from the initial commit (`7050cdd`, 2025-12-04,
just a `file_is_gzipped` stub) until **mid-February 2026**, when
serious work begins. In two weeks the foundational pieces land:
magic-bytes utilities, a hand-rolled gVCF parser with QUAL, GT, and
chrom-line handling, LRU-cached GT-format parsing, and a Criterion
bench harness.

**State at end of Feb:** ~1.7 k LOC, 10 `.rs` files. A working,
perf-tuned single-sample gVCF reader with no merging logic yet.

## Phase 2 — gVCF merger (early March → mid-April 2026)

**Objective:** the project's original purpose — read many per-sample
gVCFs in parallel, group overlapping variants across samples,
reconcile alleles into unified haplotype calls, and emit a merged VCF
with recomputed genotype posteriors.

This is where the project's original *name* (`join_per_sample_vcfs`)
makes sense. By 2026-03-25 the original objective is essentially met:
per-sample gVCFs in, merged multi-sample VCF with EM-derived genotype
posteriors out, parallel and reasonably fast.

April is then spent on **code review and hardening rather than
features.** The author introduces the project's distinctive process:
write a `code-review` skill, use it on the gVCF parser, then a
`code-review-fixes` skill to systematically apply the findings.

**State at end of April (pre-pivot):** ~9.4 k LOC, 24 `.rs` files. A
complete-but-narrow gVCF merger.

## Phase 3 — The pivot: design of a full cohort caller (mid-April → early May 2026)

**Objective:** stop being a gVCF post-processor and become a
**complete cohort variant caller from BAM** — competitive in concept
with GATK `HaplotypeCaller` + `GenotypeGVCFs` and freebayes.

The pivot is signalled by `11deac3` (2026-04-17, "first ideas for
complete Variant caller") and the GATK/freebayes source trees being
gitignored a few days later. Then the author spends about two weeks
on **specifications, not code**: the six-stage
`calling_pipeline_architecture.md`, the Stage 1 pileup-walker spec,
the `phase_chain` spec, the BufferedPeekable adapter spec, and
`design_principles.md`.

The first actual code of the new pipeline lands at the end of April:
the **CRAM input slice** (`bam/`) with header validation,
peek-and-scan merge, and a filter cascade.

## Phase 4 — The cohort caller (May 2026)

**Objective:** implement, review, perf-tune, and integrate all six
stages plus a CLI.

This is by far the most productive month: ~91.7 k inserted, ~26 k
deleted, **net +65.8 k LOC in 31 days.** The author works one stage
at a time, and each stage follows the same loop: *implementation plan
→ spec → implementation → code review → fix-application → perf review
→ fix-application*.

By 2026-05-24 the six-stage pipeline is operational end-to-end on
real tomato BAMs, with active F1-vs-GATK tuning (0.285 → 0.317 after
the MAPQ filter), per-chromosome parallelism, and a final restructure
that promotes `bam/`, `baq/`, `fasta/`, `psp/`, `pileup/`, `vcf/`,
`pileup_record`, and `iter_ext` to top-level `src/` modules.

## Milestone summary

| Date       | Objective at the time                              | Status                                                | LOC (rs) | files |
| ---------- | -------------------------------------------------- | ----------------------------------------------------- | -------- | ----- |
| 2026-02-28 | Fast gVCF parser                                   | done                                                  | ~1.7 k   | 10    |
| 2026-03-31 | Cross-sample gVCF merger w/ EM posteriors          | done; review hardening underway                       | ~6.2 k   | 18    |
| 2026-04-30 | Pivot: spec a full cohort caller from BAM          | design frozen; CRAM input slice landed                | ~9.4 k   | 24    |
| 2026-05-24 | Implement six-stage cohort caller, run on real BAM | end-to-end working; tuning F1 vs GATK                 | ~59.8 k  | 96    |

## Activity chart — commits and lines of code per active day

Each row is one active day (days with no commits are omitted).
Bars are linear-scaled at one block per ~600 lines; the added-bar
caps at 40 blocks (~24 k LOC) and the removed-bar at 20 blocks
(~12 k LOC).

```
DATE         COMMITS  ADDED  REMOVED  ADDED ━━━━━━━━━━━━━━━━━━━━━━━━  REMOVED ━━━━━━━━━━━━
2025-12-04         2     255        0
2026-02-12         6     816       88  █
2026-02-13         3     671      107  █
2026-02-14         2     323      107  █
2026-02-23         4     680      104  █
2026-02-24         6     148      397                                            █
2026-02-25         2     475       25  █
2026-03-03         1     142      147
2026-03-05         5     281       32
2026-03-06         6     504       62  █
2026-03-16         8     676      176  █
2026-03-17         8    1382      196  ██
2026-03-18        12     602      262  █
2026-03-20         7     238       59
2026-03-23         5    1061      131  ██
2026-03-24         5     815      159  █
2026-03-25        14    2011     1651  ███                                       ███
2026-03-30         2     108        0
2026-04-01         5     368      111  █
2026-04-02         1     311      123  █
2026-04-09         3      73       66
2026-04-10        39     624      563  █                                         █
2026-04-12         2     304       59  █
2026-04-13         7    2012      980  ███                                       ██
2026-04-14         5     787      208  █
2026-04-15         1     866        0  █
2026-04-17         1     347        0  █
2026-04-21         5     241      169
2026-04-22         1     222        0
2026-04-24         2    2947       84  █████
2026-04-27         3    1287        0  ██
2026-04-28         6    1767      506  ███                                       █
2026-04-29         3    3561      184  ██████
2026-05-01        14    2316      368  ████                                      █
2026-05-04         4     216      101
2026-05-06        17    6720      438  ███████████                               █
2026-05-07        23    4205      525  ███████                                   █
2026-05-08        14    7684      235  █████████████
2026-05-09        15    1889      125  ███
2026-05-10        12    1796      331  ███                                       █
2026-05-11        14    2473      499  ████                                      █
2026-05-12        30    6857      596  ███████████                               █
2026-05-13        42   21058     1821  ███████████████████████████████████       ███
2026-05-14         6    5375     2765  █████████                                 █████
2026-05-15         8    4313      567  ███████                                   █
2026-05-16        23   11015     1114  ██████████████████                        ██
2026-05-17        24   17920     1175  ██████████████████████████████            ██
2026-05-18        17    8785     1312  ███████████████                           ██
2026-05-19        35   15976     9015  ███████████████████████████               ███████████████
2026-05-20        17    3356      417  ██████                                    █
2026-05-21         1     273        0
2026-05-22         6    2771      270  █████
2026-05-23        32    9109     3103  ███████████████                           █████
2026-05-24        30    7878     2461  █████████████                             ████
```

Reading the chart: nothing of substance happens before 2026-02-12;
the gVCF-merger work in March is steady but small (≤ 14 commits/day,
≤ 2 k LOC/day); the pivot in late April shows up as the first
3000-line days (specs + CRAM input); May is structurally different
— eight days exceed 5 k LOC added, and the single biggest day
(2026-05-13, the `.psp` format) lands 42 commits and 21 k LOC.

## Feature timeline — what landed and when

Only substantive feature commits are listed (not renames, formatting,
or pure review-fix passes). The two-week refactor / review windows
between feature bursts are summarised in italics.

### Phase 1 — gVCF parser foundation

| Date       | Feature                                                                         |
| ---------- | ------------------------------------------------------------------------------- |
| 2025-12-04 | initial commit; `file_is_gzipped`                                               |
| 2026-02-12 | magic-bytes utilities; gVCF parser; variant length; VCF record QUAL; GT parsing |
| 2026-02-13 | parsing of `#CHROM` header line; genotype parsing                               |
| 2026-02-14 | perf improvement; `peek_variant` replaces `peek_items`                          |
| 2026-02-23 | invariant positions retained; faster parser; Criterion bench                    |
| 2026-02-24 | only one ploidy accepted; genotype fast path for other ploidies                 |
| 2026-02-25 | overlap variant grouping (first cross-sample primitive)                         |

### Phase 2 — gVCF merger

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-03-05 | samples info attached to grouped vars; parallel analysis of variant groups                   |
| 2026-03-06 | Python draft of merger; `OverlappingVariantGroup`; sample-origins                            |
| 2026-03-16 | merge-variant test infrastructure; simple deletion merge; overlapping deletions; het deletion |
| 2026-03-17 | dict → list state migration; phase with heterozygotes; insertion fix; Rust port of merger; CLI |
| 2026-03-18 | parallelize variant-group creation (big perf gain); skip non-variant positions; decompression on another thread |
| 2026-03-20 | missing alleles expand to whole-haplotype allele; parser skips `.` as alt                    |
| 2026-03-23 | phase calculation; VCF writer; integration test scaffold                                     |
| 2026-03-24 | integrated pipeline; double-buffered decompression                                           |
| 2026-03-25 | SNP / two-SNP / deletion / ins-del integration tests; project renamed; **EM algorithm for GT posteriors**; fixation-index support in prior; PL-less VCF support; posterior probs integrated with merging |
| 2026-04-01 | overall specification document; spec dir created                                             |
| 2026-04-02 | gVCF parser handles arbitrary genotype fields                                                |
| 2026-04-09 | `BrokenHeader` line prints offending line                                                    |

*2026-04-10:* 39-commit refactor wave — phase-vector length fix, error-type cleanups, naming pass, `RuntimeError` → `unreachable`, redundant-state removal.

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-04-12 | first draft of the `code-review` skill                                                       |
| 2026-04-13 | review skill produces a report; `code-review-fixes` skill; feature-implementation skill      |
| 2026-04-14 | critical gVCF-parser correctness fixes; regression tests; explicit ploidy on `VarIterator`/`VcfWriter` |
| 2026-04-15 | implementation plan for the joining feature                                                  |

### Phase 3 — Pivot to cohort caller from BAM

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-04-17 | first ideas for complete variant caller                                                      |
| 2026-04-21 | gitignore GATK/freebayes; **podman dev container**; IA files reorganised                     |
| 2026-04-22 | per-sample binary file spec draft                                                            |
| 2026-04-24 | **calling pipeline architecture spec** + reference-caller analysis; Stage 5 statistics, contamination, inbreeding |
| 2026-04-27 | `phase_chain` spec; Stage 6 cohort-prior extension; **`BufferedPeekable`** adapter           |
| 2026-04-28 | **Stage 1 per-sample caller spec** (CRAM → per-position records); design-principles spec; CRAM input implementation plan |
| 2026-04-29 | three core design principles (side-effects-at-edges, no-API-leak, error policy); **CRAM input slice implemented** (header validation, peek-and-scan merge, filter cascade) |
| 2026-05-01 | CRAM-input review fixes (M1–M5: u64 contig length, dedicated SM-disagreement error, malformed `@SQ` M5 error, drop unmapped reads); minor follow-ups (Mi1–Mi6) |
| 2026-05-04 | naming principles (types-name-shapes; readability-over-length); design-principle refinements |
| 2026-05-06 | **pileup walker spec**; eager-closure architecture; **Stage 1 pileup walker implemented**; rustfmt/clippy in dev container; bug fixes B1, B2, M1–M5 from pileup review |

### Phase 4 — Cohort caller implementation

**Stage 1 hardening — pileup walker** (2026-05-07 → 09)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-05-07 | samtools-comparison review; **lazy-CIGAR plan + `CigarCursor`**; binary-search cursor; slot-allocator high-water warning; secondary/supplementary always dropped; `ChromBoundaryRefFetcher` (S6); per-column depth caps (S5); samtools-matching BQ math (S7) |
| 2026-05-08 | freebayes-style pileup review; per-read mismatch-fraction filter (F1); stably left-align indels (F3); min-indel-quality choice; read-N skipped at emit; mate-pair adaptor boundary clipping (G1); ill-formed CIGAR rejection (G2); GATK pileup review tracked |
| 2026-05-09 | M1: emit `expired_chains` for orphan first mates; M2: reject chrom-id regression; Mi2-Mi9 cleanups; `rust-clone-audit` skill |
| 2026-05-10 | `rust-code-review` skill split into orchestrator + per-category checklists; `rust-performance-review` skill; **pileup performance review** + L1/L6–L11 fixes; pre-allocate folded reads and chain slots |
| 2026-05-11 | concurrency-review skill expansions; `MateRole` enum; `AlleleSupportStats` rename; lifecycle-marker semantics; `ActiveReads` rename; `RefSeqFetcher` rename; `Locus` struct |

**BAQ** (2026-05-12)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-05-12 | **BAQ feature plan + implementation**: htslib `probaln_glocal` port (HMM core), `BaqEngine` driver, **rayon-parallel `BaqStream` adapter**, BAQ bench; per-target CPU floor (`x86-64-v3`, `apple-m1`); BAQ perf-review wins L1–L8, L3 by-value `MappedRead`; round-2 pileup perf review H2/H3 applied (manual SIMD guidance) |

**`.psp` per-sample file format** (2026-05-13)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-05-13 | **`.psp` byte-format spec closed**; PSP writer steps 1–5 (column registry, varint, block index/trailer, TOML header, per-shape codecs+zstd, streaming `PspWriter`); writer perf bench + review with **zstd CCtx reuse, CSR list-column layout, ACGTN lookup table, varint fast/cold split, slab-cast list-column emit on LE, `BTreeSet`→sorted `Vec` for active slots**; project methodology rules (`rust-toolchain.toml` pin to 1.95, `panic = "abort"`, opt-in `alloc-mimalloc` feature); **PSP reader implemented** (all 7 plan slices); reader perf review + fixes |

**Stage 1 wiring & cohort skeleton** (2026-05-14 → 16)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-05-14 | **pull-iterator walker + `pileup_to_psp` seam**; phase chains switched from recycled `u16` slot ids to unique `u64` chain ids; `pending_mates` 10 k cap; **`pop_var_caller pileup` subcommand** wiring CRAM → BAQ → walker → `.psp` |
| 2026-05-15 | **`psp-to-pileup` subcommand** (text dump); slot → chain-id rename across code+specs; IA reorg to `doc/devel` + `doc/help`; **linear-scan k-way merger** (`PerPositionPileups`) over per-sample `.psp` iterators; soft-masked FASTA bytes uppercased at the seam |
| 2026-05-16 | **Stage 4 variant-grouping plan + streaming variant grouper** (overlap bundling); obsolete gvcf bench removed; **Stage 5 per-group merger plan + implementation** (allele unification + likelihood reconstruction); review-fixes wave (M4 typed errors, M7 zero-config constructors deleted, M14 renames, M-batches A–D); **`cohort/` → `var_calling/`, `per_sample_caller/` → `per_sample_pileup/` rename**; **`PROJECT_STATUS.md`** as per-feature lifecycle index; **end-to-end `cohort_perf` bench** (three-stage Criterion); posterior-engine perf review + Wave 1 fixes; Stage 6 posterior-engine plan |

**Stage 6 posterior engine + contamination** (2026-05-17 → 18)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-05-17 | **Stage 6 posterior engine v1**; contamination_estimation spec + **online-EM side-pass** + Stage 6 mixture E-step; **Stage 3 sdust filter** (pivoted from BLAST DUST) + review/fixes; **`MathBackend` trait** + `ExactMath`; safe-ln hoist; genotype-shape cache + `e_step` restructure; **`InterpUnivariateMath`** (IEEE-decomp tables, tuned to 256 bins); accuracy harness; default flipped to interp; **SIMD primitives + lane-of-4 `MathBackend`** (`InterpUnivariateSimdMath`); SIMD e_step path |
| 2026-05-18 | route SIMD through `wide`'s native ln/exp; default flipped to `InterpUnivariateSimdMath`; **cohort VCF writer** (Stage 6 sink); compute_mixture per-sample alloc hoist; SIMD mixture pre-pass for HAS_LANE_4; biallelic-diploid `m_step_p_hat` fast path; cohort_vcf_writer review fixes (1 Blocker + 14 Majors + 16 Minors); H2+H4 `log_indep` + homogeneous-fixation hoist; per-engine `RecordScratch` lift (L5–L7) |

**Cohort CLI integration** (2026-05-19)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-05-19 | **cohort CLI implementation plan (Tasks 1–9)**: engine-side Config validators; `clap` `value_parser` helpers; batch-assignment TSV reader; **contamination artefact TOML I/O**; **`estimate-contamination`, `var-calling`, `var-calling-from-bam` subcommands**; Stage 1 pipeline helper extraction; **cohort CLI integration tests**; review + 5-wave follow-up fixes (`#[non_exhaustive]`, builder pattern, version gate, British-spelling rename, shared-infrastructure refactor, idempotent rayon-pool helper, test infra + **FASTA MD5 enforcement**); **legacy gVCF-merger graph deleted**; **project renamed `merge_per_sample_vcfs` → `pop_var_caller`**; **MIT license added**; **end-to-end `cohort_e2e_perf` bench** (`.psp` → cohort-VCF) |

**Real-data correctness** (2026-05-20 → 22)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-05-20 | first runs vs GATK on tomato data; `DEFAULT_CONVERGENCE_THRESHOLD` relaxed 1e-4 → 1e-3; **walker chain_id projection bug fix at finalise**; perf review + **H1 per-chromosome parallel plan**; **pure-Rust bgzf-aware concat helper** for per-chrom VCF fragments; `process_one_chromosome` per-worker driver; intra-batch `par_iter` dropped; **per-chrom parallelism live**; `DidNotConverge` → emit-with-flag; **hom-ref-only sites dropped at writer**; `--min-qual` 30 hygiene filter; **QUAL formula rewritten** to `-10·log10 P(AC=0|data)` via exact-AF 1-D convolution under a Beta-Binomial cohort prior; `--min-alt-obs-per-sample` filter (default 2) |
| 2026-05-21 | session report — caller correctness vs GATK                                                  |
| 2026-05-22 | **per-allele MAPQ scalars** (`mapq_sum` + `mapq_sum_sq`) threaded through PSP + pipeline; **`INFO/MQRef, MQAlt, MQDiff, MQDiffT`** emitted; **`--min-mapq-diff-t` Welch's-t MAPQ filter** (default −3.0, F1 0.285 → 0.317); FASTA MD5 streamed in 64 KiB windows + per-chrom parallel; per-worker `SingleChromRefFetcher` |

**Reference-fetcher unification + module restructure** (2026-05-23 → 24)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-05-23 | **`StreamingChromRefFetcher`** (1 MB sliding buffer); unused `SingleChromRefFetcher` deleted; **`ChromRefFetcher` trait + impl**; `cohort_driver` routed through new fetcher; Stage 1 pipeline routed via `WalkerLegacyAdapter`; BAQ + walker share one adapter; `var_calling_from_bam` swaps `SyncRefFetcher` → `WalkerLegacyAdapter`; legacy `SyncRefFetcher` + `ChromBoundaryRefFetcher` deleted; **`ManualEvictChromRefFetcher`** for BAQ workers; `PerChromRecordsIter` helper; DustFilter + PerGroupMerger + `drive_cohort_pipeline` migrated; **ref_fetcher perf review** (drop `Sync`, `Mutex` → `RefCell`, `fetch_into` trait method); **PSP reader perf review** (BufReader-buffer preservation across block seeks, CSR collapse, per-allele bounds-check hoist); fix-application wave (malformed `.fai` rejection, non-ACGT-byte folding to N, trait sealing, `iter_bases` borrow scope, doc `# Errors`, missing-test fixes); **legacy `RefSeqFetcher` retired** in favour of typed-error multi-chrom fetcher; `open_contig` helper; `chrom_name`/`chrom_length` unification |
| 2026-05-24 | reports moved to canonical dir; **`scripts/precommit-check.sh`**; Rust 1.95 clippy fixes; broken benches/examples ported off retired trait; **module restructure**: `BufferedPeekable` → top-level `iter_ext`; `baq/` promoted to top-level (algorithm core + pileup glue flattened); `pileup_record` promoted to top-level; `psp/` promoted to top-level; `fasta/` (ref_fetcher + `ContigList` + `MultiChromRefFetcher` trait) promoted to top-level; **`pileup/` bundled as `walker/` + `per_sample/` submodules**; `pileup_record` flattened to single file; **`vcf/` promoted; `vcf_writer` decoupled from `PosteriorRecord` via a `VcfWritable` trait**; **`bam/` (CRAM input slice) promoted to top-level**; module-structure refactoring lessons captured in skills |

## Open work as of 2026-05-24

From the current branch state (uncommitted changes plus the most
recent commits), the active work is:

- **Per-chromosome var-calling from BAM** — the new
  `doc/devel/implementation_plans/var_calling_from_bam_per_chromosome.md`
  and uncommitted edits to `src/pop_var_caller/var_calling_from_bam.rs`
  suggest the per-chromosome parallel mode (already live for the
  `.psp`-input path) is being extended to the BAM-input subcommand.
- **BAM index preflight** — new untracked file
  `src/bam/index_preflight.rs`.
- **F1 vs GATK tuning** — currently at 0.317 on a synthetic cohort
  after the MAPQ filter; the DUST filter is kept on as a policy
  choice even though `--no-complexity-filter` scored a touch better,
  on the grounds that GATK has no DUST so agreement in low-complexity
  regions isn't truth.

## Process observations

Two things are unusual and visible in the log:

1. **Skill-driven workflow.** Reviews aren't ad-hoc — the author wrote
   reusable `code-review`, `code-review-fixes`,
   `rust-performance-review`, and `feature-implementation` skills
   early and uses them mechanically on every new module. Roughly **a
   third of all commits in May are review or fix-application
   commits**, not feature code.
2. **Spec-first for the pivot.** The transition from gVCF merger to
   cohort caller was preceded by ~10 days of spec writing with no
   implementation. That spec
   (`doc/devel/specs/calling_pipeline_architecture.md`) is referenced
   by name in commit messages months later — it's the project's
   actual source of truth, not the code.
