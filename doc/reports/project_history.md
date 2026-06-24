# Project history — `join_per_sample_vcfs` → `pop_var_caller`

**Span:** 2025-12-04 → 2026-06-24 (~6.5 months, ~1040 commits)
**Source growth:** 10 → 129 Rust files in `src/`, ~1.7 k → ~80.1 k LOC.
**Authoring model:** single primary author working with a Claude Code
assistant; commit messages and per-stage `doc/devel/reports/` and
`doc/devel/specs/` artefacts reflect that workflow. The assistant model
itself moved across the project (Claude Opus 4.6 → 4.7 → 4.8, mostly
in the 1M-context variant) and the collaboration tooling matured — see
[AI tooling evolution](#ai-tooling-evolution).

The history breaks into **six phases** across two arcs. The first arc
(Phases 1–4) is the SNP/indel cohort caller; the second (Phases 5–6)
hardens it and then pivots a second time into microsatellite (SSR/STR)
genotyping. Two inflection points stand out:

1. **End of Phase 2 → Phase 3** — the project's original goal (joining
   per-sample gVCFs) was achieved and then deliberately discarded in
   favour of a much larger goal: a complete cohort variant caller from
   BAM.
2. **Within Phase 6** — with the SNP caller working end-to-end on real
   data, the author opens an entirely new front: a locus-oriented
   microsatellite genotyper (`ssr-catalog` / `ssr-pileup` / `ssr-call`)
   sharing the `.psp` container and cohort statistics but with its own
   realignment core.

*Phases 1–4 are summarised below; the detailed per-day record for them
lives in the [feature timeline](#feature-timeline--what-landed-and-when).
This revision (2026-06-24) extends the document through Phases 5–6.*

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

## Phase 5 — Re-architecture, perf, and memory (late May → mid June 2026)

**Objective:** the caller works, but on real cohorts it is too slow
and too memory-hungry. This month is spent making it *competitive* —
chunk-parallel within a chromosome, a record-streaming pipeline that
trades RAM for sample-count scaling, and a long sequence of measured
perf/memory levers — while hardening correctness against GATK and
freebayes on real tomato and human data.

The defining move is the **chunk-parallel rewrite** and then the
**`re-architect` record-streaming pipeline.** Rather than adapt the
old streaming iterators, the author treats them as a spec and builds a
new worker natively over columnar `.psp` data: column-native allele
unification and a column-native EM, per-chunk DUST masking (replacing
a serial whole-chromosome pre-pass), and a clean split between
*building* the data partitions and *running* the math on them. The
`re-architect` branch (Phase 0 → Phase 7) reaches a **byte-identity
milestone** against `main`, gets a crossbeam producer→caller→writer
topology with parallel per-sample decode, and is then swapped into
production (2026-06-04). The direct `var-calling-from-bam` path is
deleted — only `pileup → .psp → var-calling` remains.

Memory and throughput are then driven down through a documented series
of levers, most of them landing byte-identical: two-phase
column-selective decode, dropping REF-allele chain-ids, a `--low-memory`
two-pass producer, lowering the default block window 20 kb → 5 kb (the
single biggest RSS knob), and a single-pool thread budget that honours
`--threads N` with a producer/caller split. Correctness work in
parallel: indel left-alignment (a GATK-port normaliser made
mandatory), QUAL calibration that deflates scores at
systematic-artifact sites (the depth-inflation false-positive
problem), per-allele MAPQ filtering, and a BED `--regions` mode.

This phase is also where the **per-step review loop becomes
mechanical** — the giant `cohort_block` review (5 Blockers, 32 Major,
26 Minor) and its ~40 fix commits in a single day, then conventional
`feat → review → fix` triplets thereafter.

## Phase 6 — Second pivot: SSR/STR genotyping (June 2026)

**Objective:** build a microsatellite (short-tandem-repeat) genotyper
as a new, locus-oriented module — reusing the `.psp` container and the
cohort statistics machinery, but with its own realignment core,
because read mappers handle repeats badly.

The pivot is spec-first again, and heavily research-driven: a research
report digests HipSTR / GangSTR / ConSTRain, an architecture doc
settles placement (same crate, `src/ssr/`, CLI mirroring the SNP
`pileup`/`var-calling` split), and the `.psp` format is generalised
into a generic core with SNP and SSR specialisations (SNP output stays
byte-identical). The implementation then walks three stages:

- **Stage 0 — `ssr-catalog`:** build the locus catalog from a
  `trf-mod` tandem-repeat scan, dropping period-1 homopolymers and
  compound loci (HipSTR/GangSTR-style).
- **Stage 1 — `ssr-pileup`:** per-locus read genotyping. The decision
  is "realign everything" — every spanning read goes through a
  pair-HMM rather than trusting the mapper's CIGAR — with
  `count_repeats` kept only as a *measured* fast-path candidate (built,
  measured, verdict "don't build it"). A **Mark-2 rebuild** then
  replaces the first cut with an empirical-candidate allele model and a
  Viterbi delimiter, and the realignment is tuned ~54 % faster
  byte-identically.
- **Stage 2 — `ssr-call`:** the cohort statistics layer. A Mark-2 spec
  is settled end-to-end and then **hardened against a 13-finding
  adversarial review.** The reading/merge layer lands first (typed
  per-sample readers + a catalog-driven k-way merger), then the
  genotyping + parameter pre-pass: a cohort simulator, shared locus
  primitives (stutter model, alignment), a sum-over-slips likelihood,
  per-locus EM, sample-group clustering, a frozen-parameter pre-pass
  (ε, stutter shape θ per period refined per-locus, per-group stutter
  level, per-individual F, G₀ decay), and finally a chunk-parallel
  genotyping sweep with byte-identity across threads.

Shared infrastructure also advanced underneath both modules: a
`segment_reader` indexed read source replaced the old per-file query
path, and the FASTA reference reading was unified onto the pileup path.

As of 2026-06-24 the SSR genotyper runs end-to-end (`ssr-catalog` →
`ssr-pileup` → `ssr-call` → VCF) with the parameter pre-pass and
per-locus refits in place; the SNP caller continues to be tuned in
parallel.

## Milestone summary

| Date       | Objective at the time                              | Status                                                | LOC (rs) | files |
| ---------- | -------------------------------------------------- | ----------------------------------------------------- | -------- | ----- |
| 2026-02-28 | Fast gVCF parser                                   | done                                                  | ~1.7 k   | 10    |
| 2026-03-31 | Cross-sample gVCF merger w/ EM posteriors          | done; review hardening underway                       | ~6.2 k   | 18    |
| 2026-04-30 | Pivot: spec a full cohort caller from BAM          | design frozen; CRAM input slice landed                | ~9.4 k   | 24    |
| 2026-05-24 | Implement six-stage cohort caller, run on real BAM | end-to-end working; tuning F1 vs GATK                 | ~59.8 k  | 96    |
| 2026-06-07 | Re-architect for scale (record-streaming, chunk-parallel) | byte-identical swap into production; perf/memory levers landing | ~70 k    | ~110  |
| 2026-06-17 | SSR Stage 0+1 (`ssr-catalog`, `ssr-pileup`) on real data | Mark-2 rebuild done; per-locus reads genotyped       | ~75 k    | ~120  |
| 2026-06-24 | SSR Stage 2 (`ssr-call`) cohort statistics + pre-pass | end-to-end SSR VCF; per-locus param refits in EM     | ~80.1 k  | 129   |

## Commits by category

569 commits were classified by subject-line keywords (category-code
fixes like `(M3)` or `Fix Mi7`, module prefixes like `psp::reader:`,
explicit verbs like `Implement`, `Add`, `Fix`, `Apply`, `Review`,
etc.). Each commit gets one primary category — review-driven
correctness fixes count as review even when the diff is large,
perf-review fixes (the `L`/`H` category codes) count as perf, etc.

| Category                         | Commits | Share |
| -------------------------------- | ------- | ----- |
| **Implementation (new features)** | **192** | **33.7%** |
| **Code review (correctness fixes)** | **141** | **24.8%** |
| **Performance**                  | **94**  | **16.5%** |
| Refactor / rename / cleanup      | 82      | 14.4% |
| Docs / specs / plans / reports   | 47      | 8.3%  |
| Tooling / container / skills     | 10      | 1.8%  |
| Merge / revert / other           | 3       | 0.5%  |

Headline: the three categories the user asked about account for
**~75% of all commits**, split roughly 4 : 3 : 2 (implementation :
review-fix : perf). Review and perf together (**41%**) are almost as
large as raw implementation — the project spends as much effort
reviewing and tuning code as writing it new.

### Category mix by month

The mix shifts dramatically over the project's life: the gVCF-parser
months were almost pure implementation; the post-merger month
(2026-04) was the first concentrated review pass; May is where every
category fires at once.

| Month   | impl | review | perf | refactor | docs | tooling | total |
| ------- | ---: | -----: | ---: | -------: | ---: | ------: | ----: |
| 2025-12 |    2 |      0 |    0 |        0 |    0 |       0 |     2 |
| 2026-02 |   15 |      0 |    4 |        4 |    0 |       0 |    23 |
| 2026-03 |   50 |      0 |    3 |       14 |    2 |       3 |    73 |
| 2026-04 |   18 |     19 |    0 |       34 |   12 |       1 |    84 |
| 2026-05 |  107 |    122 |   87 |       30 |   33 |       6 |   387 |

What the table shows:

- **Feb–Mar 2026 (parser + merger):** zero review commits, almost no
  perf. The author is in pure-feature mode — implementations and
  small refactors only.
- **April 2026:** the first 19 review commits land. This is the
  introduction of the `code-review` and `code-review-fixes` skills
  (12 doc commits = the skill files themselves) and their first
  application against the gVCF parser. Implementation drops from 50
  → 18 — features stop, hardening starts.
- **May 2026:** review (122) actually exceeds implementation (107),
  and perf surfaces as its own column with 87 commits. This is the
  six-stage pipeline being built and *immediately* reviewed and
  tuned, stage by stage.

### Caveats on classification

- Commits often do more than one thing; each is bucketed by its
  apparent primary intent. A commit like "Apply H1 + H4 from the
  round-2 pileup perf review" counts as perf, not review.
- "Refactor" includes renames, formatting, module promotions, and
  small cleanups. Many of these are review-driven housekeeping; the
  share of "review-influenced" work is therefore higher than the
  review column alone suggests.
- The 39-commit 2026-04-10 burst (gVCF parser review) is bulk-tagged
  as review since every commit in that day was a category fix from
  the same review.

### Second period (2026-05-25 → 2026-06-24)

The classification above was a snapshot at the document's creation
(2026-05-24). The following month added **456 commits**, re-classified
with the same intent-based method (the move to conventional-commit
prefixes mid-period actually makes this *easier* — `feat`/`fix`/`perf`
scopes map almost directly):

| Category                            | Commits | Share |
| ----------------------------------- | ------- | ----- |
| **Implementation (new features)**   | **135** | **29.6%** |
| **Code review (correctness fixes)** | **91**  | **20.0%** |
| Docs / specs / plans / reports      | 81      | 17.8% |
| Tooling / bench / container         | 38      | 8.3%  |
| **Performance**                     | **36**  | **7.9%** |
| Refactor / rename / cleanup         | 33      | 7.2%  |
| Merge / revert / other              | 28      | 6.1%  |
| Tests                               | 9       | 2.0%  |

The proportions are remarkably stable versus the first period —
implementation ~30 %, review ~20 % — but **docs nearly doubles its
share** (8 % → 18 %). That is the spec-first, per-step methodology
compounding: every SSR stage gets an architecture sketch, an
implementation plan, and a per-step review report before and around
the code. Review and docs together (~38 %) now *exceed* raw
implementation.

The monthly mix for the window:

| Month        | impl | review | perf | refactor | docs | tooling | test | merge | total |
| ------------ | ---: | -----: | ---: | -------: | ---: | ------: | ---: | ----: | ----: |
| 2026-05 (25–31) |   50 |     33 |    3 |       10 |   12 |      17 |    3 |     2 |   134 |
| 2026-06         |   85 |     58 |   33 |       23 |   69 |      21 |    6 |    26 |   322 |

- **Late May** is the chunk-parallel rewrite and the `cohort_block`
  review — note the 50 impl / 33 review split inside one week, plus a
  heavy tooling column (the benchmark dashboards and dev-container
  caller comparisons).
- **June** is structurally a two-front month: the re-architecture and
  perf levers (33 perf commits — the highest of any month) on the SNP
  side, and the SSR pivot driving the 69 docs commits (research,
  architecture, specs, per-step reports) and the 26 merges (the many
  short-lived feature branches: `re-architect`, `low-memory-mode`,
  `bed-regions`, `qual-analysis`, `segment-read-fetcher`,
  `ssr-architecture`, `ssr-pileup-mark2`, `ssr-cohort`).

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
2026-05-25         2     389       36  █
2026-05-26         5    2532        6  ████
2026-05-27         2     822       93  █
2026-05-28        19   11053      426  ██████████████████                        █
2026-05-29        80   14026     4043  ███████████████████████                   ███████
2026-05-30        10    1328     1395  ██                                        ██
2026-05-31        16    4344     1261  ███████                                   ██
2026-06-01        32    3380     6238  ██████                                    ██████████
2026-06-02         4     928      230  ██
2026-06-03        22    4035      239  ███████
2026-06-04        16   21735    31091  ████████████████████████████████████      ████████████████████
2026-06-05         4     593      297  █
2026-06-06         5    1752      162  ███
2026-06-07        12    1303      122  ██
2026-06-08        20    2296     1128  ████                                      ██
2026-06-09        17    1957      181  ███
2026-06-10         7    1822     1105  ███                                       ██
2026-06-11        10    5002      604  ████████                                  █
2026-06-12        12    2802      303  █████                                     █
2026-06-14         2    1184       67  ██
2026-06-15        28    9508      950  ████████████████                          ██
2026-06-16        42    8396     4721  ██████████████                            ████████
2026-06-17        20    7959     6550  █████████████                             ███████████
2026-06-19         3    1592      418  ███                                       █
2026-06-21        13    3246      250  █████
2026-06-22         2    3465      256  ██████
2026-06-23        42   11332      459  ███████████████████                       █
2026-06-24         9    1245      149  ██
```

Reading the chart: nothing of substance happens before 2026-02-12;
the gVCF-merger work in March is steady but small (≤ 14 commits/day,
≤ 2 k LOC/day); the pivot in late April shows up as the first
3000-line days (specs + CRAM input); May is structurally different
— eight days exceed 5 k LOC added, and the single biggest day
(2026-05-13, the `.psp` format) lands 42 commits and 21 k LOC.

The Phase 5–6 tail (from 2026-05-25) reads differently: the spikes are
no longer pure growth. **2026-06-04** is the largest single day of the
whole project by *churn* (+21.7 k / −31.1 k) but barely grows the tree
— it is the `re-architect` swap, where a whole new pipeline replaces
the old one in place. **2026-05-29** (80 commits) is the `cohort_block`
review-and-fix avalanche, and **2026-06-23** (42 commits, +11 k) is the
SSR `ssr-call` genotyping milestones A–J landing in one sitting. The
high-deletion days (06-01, 06-04, 06-17) mark deliberate removals — the
direct BAM path, the old pipeline, the Mark-1 SSR module — rather than
churn from indecision.

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

### Phase 5 — Re-architecture, perf, and memory

**Chunk-parallel within-chromosome rewrite + `cohort_block` review** (2026-05-25 → 31)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-05-25 | deps update (`wide` 0.7→1.4, `lru` 0.16→0.18, drop `flate2`); `--max-alleles-lh-calc` fixes a release-build OOB on >64 alleles |
| 2026-05-26 | **dev container gains samtools + freebayes + GATK** for caller comparison; Apple `container` as podman fallback; tomato caller-comparison benchmarks promoted from `tmp/` |
| 2026-05-27 | `.psp` block target exposed on the pileup CLI; default 16 MiB → 512 KiB → 1 MiB |
| 2026-05-28 | **within-chrom chunk-parallel rewrite begins**: `drive_cohort_chunked` driver; chunk loader + pre-pass + variant-group partitioner; worker adapter reusing `per_group_merger`/`posterior_engine`; **Phase A.1 column-native allele unification** (per-position projection, compound detection, max-alleles cap, column-native log-likelihoods) |
| 2026-05-29 | **80-commit day**: Phase A.2 column-native EM (split EM along row/columnar boundary); **`cohort_block` code review (5 Blockers, 32 Major, 26 Minor)** + ~40 fix commits across Waves 1–7; **indel normalization** — pure left-alignment CIGAR rewrite (GATK-port F3 replacement), wired into the BAQ read-prep stage, made mandatory; ALT-allele pruning ported from main; accuracy/QUAL-sweep dashboards |
| 2026-05-30 | per-chunk DUST mask (replaces whole-chromosome pre-pass); **`BlockIterator` extracted — decouple block production from math**; chunk produce (data-shaping) split from consume (math); parallelize loader fold + per-sample chunk loading; drop sub-chunk windows |
| 2026-05-31 | span-addressable columnar PSP reader; streaming block loader (memory fix); DUST-ahead queue off the critical path; parallel ordered block consumption; `cohort_block` → `from_psp`, direct path → `from_bam/` renames |

**Direct-path removal, `re-architect` record-streaming pipeline, perf/memory levers** (2026-06-01 → 09)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-06-01 | **direct `var-calling-from-bam` path deleted** (only `pileup → .psp → var-calling` remains); `var_calling/` review + perf review (H1 single `log_sum_exp_slice` in QUAL convolution, H2 compound-detection churn, H3 inline SmallVec keys); `from_psp/` merged into parent module; chunk hot-path benches |
| 2026-06-02 | **`.psp` blocks cut on a fixed genomic grid** for cross-sample alignment; `run_ours` benches go through `pileup → psp → var-calling` |
| 2026-06-03 | **BED `--regions` feature** (RegionSet + parser, sub-contig query, region-driven pileup, command line + regions recorded in `.psp` header); **`--low-memory` two-pass producer** + load-time min-alt-obs filter pushdown; **linear-domain QUAL allele-count convolution** |
| 2026-06-04 | **`re-architect` record-streaming pipeline** (16-commit day, +21.7 k/−31.1 k): Phase 0 scaffold + byte-identity oracle → per-sample segment reader → cohort keep/cut math over light columns → streaming cohort producer → record-based `VariantCaller` → **BYTE-IDENTITY MILESTONE** → crossbeam producer→caller→writer topology → parallel per-sample `.psp` decode → `--regions` support → **Phase 7 swap into production**; contamination ported to the record-based merger |
| 2026-06-05 | re-architecture code review + fix waves (M1–M8, Mi2/Mi8) |
| 2026-06-06 | **two-phase column-selective decode** wired into the cohort producer (−40 % live heap); cohort perf review; thread-oversubscription diagnostic playbook |
| 2026-06-07 | **`re-architect` → `main` merge**; column-selective decode foundation; **drop REF-allele chain-ids end-to-end (H2)**; parallelize the producer fold; inline `ln_factorial` fast path (H3); independent producer rayon pool size (H1); coz causal profiler added to the container; threads × samples memory-scaling dashboard panels |
| 2026-06-08 | **default block window 20 kb → 5 kb** (the single biggest RSS knob, −58 % at N=50); `var_calling` code review + fixes (typed `PipelineError`, `pub(crate)` demotion, dead-code cleanup); split producer into fold/plan + compact + read stages; `[profile.profiling]` for symbolised dhat/perf stacks; default target-variants-per-chunk 256 → 128 |
| 2026-06-09 | **thread-budget single-pool**: honour `--threads N` with a producer/caller split (3N+c → N+c); thread-budget integration guard; first SSR research/design docs land (Stage-1 two-tier pair-HMM, file-format interfaces) |

**QUAL calibration, FASTA unification, pileup de-barrier** (2026-06-10 → 12)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-06-10 | **SSR/STR genotyping specification** + validation/Stage-2 detail; **FASTA reference readers built once per run** (fix `--regions` perf regression); walker reads reference from the shared `Repository` |
| 2026-06-11 | **QUAL deflation at systematic-artifact sites** (depth-inflation FP fix) + depth-sweep / QUAL-vs-depth dashboards (`qual-analysis` merge); QUAL-cutoff precision/recall control; bound the QUAL-refine binomial tail to avoid an O(depth) hang |
| 2026-06-12 | **SSR shared types** (`Motif`, `Locus`) + review; **Stage-1 pileup de-barriered** into a source/worker/reorder pipeline; pileup `--threads` default → 4; exact incomplete-beta QUAL tail above the n cap |

### Phase 6 — SSR/STR genotyping

**Stage 0 catalog + `.psp` container generalization** (2026-06-14 → 16)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-06-14 | Stage-1 `ssr-pileup` architecture settled + adversarial review |
| 2026-06-15 | **`.psp` container generalized** to a generic core + SNP/SSR specialisations via a `PspKind` trait (steps 1a→5, SNP byte-identical, SSR round-trips); **Stage 0 `ssr-catalog`** (format I/O, `trf-mod` spawn + BED parse, post-processing); **Stage-1 `ssr-pileup` core**: realign-everything triage, `count_repeats` fast-path (measured), pair-HMM forward scorer, on-ladder + off-ladder candidate generation, dense per-read scoring, all-CSR locus storage; `trf-mod` installed in the container |
| 2026-06-16 | **`ssr-catalog` + `ssr-pileup` CLI subcommands runnable end-to-end**; drop period-1 homopolymers; **`segment_reader` shared indexed read source** + `SegmentMergedReads` multi-file coordinate merge (replaces the old indexed-query path); per-locus read fetcher (reach gate + reservoir); Stage-1 driver parallelized (batched rayon); SNP `--regions` retrofit onto the pooled reader; bed-regions `--regions` perf + review |

**Stage 1 Mark-2 rebuild** (2026-06-16 → 17)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-06-16 | `ssr-pileup` perf: per-worker `CachingCramReader` (decode each container once), `par_chunks` warm-cache scaling, **~54 % faster pair-HMM** (byte-identical), default pair-HMM window 10 → 6; fast-path investigation verdict "don't build it" |
| 2026-06-17 | **Mark-2 rebuild**: rename `ssr/` → `ssr_mark1/`, scaffold Mark-2 (types + Stage 0 catalog), read path (`fetch_reads` + footprint), Viterbi+traceback delimiter, per-locus tally (observed sequences + counts), storage + driver + cutover (delete Mark-1); **empirical-candidate allele model**; Mark-2 code review (3 Blockers + 9 Majors) + fixes |

**Stage 2 `ssr-call` — cohort statistics** (2026-06-19 → 24)

| Date       | Feature                                                                                      |
| ---------- | -------------------------------------------------------------------------------------------- |
| 2026-06-19 | **Mark-2 cohort spec settled end-to-end** + **hardened against a 13-finding adversarial review** (empirical candidates, locus-admission motif filter, HipSTR-informed sum-over-slips likelihood, per-locus EM) |
| 2026-06-21 | **`ssr-call` reading layer** (Phases 0–3): owning typed record iterator (`PspReader::into_records_of`), `SampleEvidenceCursor` per-sample reader, catalog-driven k-way merger, single-threaded end-to-end driver; reading-layer code review + fixes |
| 2026-06-22 | `ssr-call` adversarial re-verification + fixes + confident-genotype seed; global consistency pass |
| 2026-06-23 | **`ssr-call` genotyping + parameter pre-pass** (42-commit day, Milestones A–J): A1 core types + A2 cohort simulator → B shared locus primitives (rungs, stutter, align) → C1–C4 candidate assembly + likelihood + prior/seeds + EM (**checkpoint 1**) → D1–D3 parameter pre-pass + sample-group clustering (**checkpoint 2**) → E prior-side F + stutter-level outer loop + FP control + VCF output → **F parallelism + byte-identity across threads** → G1 fit G₀ decay → H `build_param_set` (pre-pass → frozen `ParamSet`) + merger accessors + VCF header writer + two-pass streaming driver → VCF |
| 2026-06-24 | **per-locus EM refits**: I1 θ_locus shape refit, I2 per-locus stutter-rate refit; **J chunk-parallel genotyping sweep** (each with its own per-step review + fix) |

## Open work as of 2026-06-24

From the most recent commits and branch state (`ssr-cohort`), the
active work is:

- **SSR `ssr-call` genotyping** — the Stage-2 cohort caller is the
  live front: the EM now refits θ_locus shape and per-locus stutter
  rate (Steps I1/I2), and the genotyping sweep was just made
  chunk-parallel (Step J). The roadmap (`doc/devel/implementation_plans/ssr_call_roadmap.md`,
  16 steps A1→F2) still has milestones beyond J pending; the next
  checkpoints are the param-recovery and first-real-VCF gates.
- **SNP caller tuning continues in parallel** — QUAL calibration
  against the depth-inflation false-positive problem, and the
  memory-for-scaling thesis (peak-RSS levers) remain open lines on the
  `main` side.
- **Pending review debt** — the `.psp` §10 container generalization
  was flagged for a fresh-conversation code review; several
  reading-layer and spec docs are owed (reading-layer module doc,
  Mark-1 amendments, extend-vs-fork engine decision).

## AI tooling evolution

The user explicitly tracks how the *collaboration* changed alongside
the code. Three shifts are visible in the log:

**1. The assistant model moved up twice.** The `Co-Authored-By`
trailers record the model in use: **Claude Opus 4.6** (first seen
2026-04-14) → **Opus 4.7** (2026-04-21 → 2026-05-29) → **Opus 4.8**
(from 2026-05-29 onward), predominantly in the **1M-context** variant.
The model upgrades line up with the project's most ambitious bursts —
4.7 carried the six-stage build, and 4.8 carried the re-architecture
and the entire SSR module.

**2. Commit conventions tightened.** Phases 1–4 used freeform,
module-prefixed subjects (`cohort_block: …`, `psp::reader: …`,
`review fixes [M5]`). From early June the project adopts
**conventional-commit prefixes with scopes** — `feat(ssr):`,
`perf(var_calling):`, `docs(review):`, `fix(ssr): apply … review` —
which makes the *intent* of each commit machine-readable (and made the
category re-classification above almost mechanical).

**3. The review loop became a fixed, visible cadence.** Where Phase 2
introduced the `code-review` / `code-review-fixes` skills as ad-hoc
tools, by Phase 6 every increment is a rigid triplet in the log:
`feat(ssr): <step>` → `docs(ssr): code review of <step>` →
`fix(ssr): apply <step> review`. More than twenty such triplets appear
in June alone. Two heavier-weight practices also matured:

- **Adversarial spec review.** Major specs are now stress-tested
  before implementation — the Mark-2 cohort spec was hardened against
  a **13-finding adversarial review**, and the `ssr-call` design got a
  second adversarial *re-verification* pass. This is spec-first
  (Phase 3's lesson) plus an explicit attempt to break the design.
- **Research-driven design.** The SSR pivot opened with literature
  research reports (HipSTR / GangSTR / ConSTRain digested into
  `doc/devel/reports/research/`) feeding the architecture — a heavier
  research front than anything in the SNP arc.

Supporting tooling also grew: a profiling stack (coz causal profiler
in the dev container, a `[profile.profiling]` for symbolised
dhat/perf stacks), `scripts/precommit-check.sh`, the dev container
gaining samtools/freebayes/GATK/`trf-mod` for caller comparison, and a
persistent project-knowledge memory index that carries decisions
across conversations.

## Process observations

Several patterns are unusual and visible across the whole log:

1. **Skill-driven workflow.** Reviews aren't ad-hoc — the author wrote
   reusable `code-review`, `code-review-fixes`,
   `rust-performance-review`, and `feature-implementation` skills
   early and uses them mechanically on every new module. Roughly **a
   third of all commits in May are review or fix-application
   commits**, not feature code; in June the review + docs share
   actually *exceeds* implementation.
2. **Spec-first for every pivot.** The transition from gVCF merger to
   cohort caller was preceded by ~10 days of spec writing with no
   implementation. That spec
   (`doc/devel/specs/calling_pipeline_architecture.md`) is referenced
   by name in commit messages months later — it's the project's
   actual source of truth, not the code. The SSR pivot repeats the
   pattern at higher intensity: research report → architecture doc →
   per-stage spec → adversarial review, *then* code.
3. **Rewrites are byte-identity-gated, not faith-based.** The two
   biggest restructurings — the chunk-parallel rewrite and the
   `re-architect` record-streaming pipeline — were both built behind a
   byte-identity oracle that proves the new path produces the *exact*
   same VCF as the old one before it is allowed to replace it. Perf
   and memory levers are then landed one at a time, each labelled
   byte-identical (or explicitly flagged when it changes output, e.g.
   QUAL). Measurement precedes the optimisation and the verdict is
   recorded even when it is "don't build it."
4. **The project is willing to throw work away.** It has pivoted away
   from its own name twice (gVCF merger → cohort caller; then a second
   front into SSR), deleted the direct BAM path, the old pipeline, and
   the entire Mark-1 SSR module once better designs existed. The
   high-deletion days on the activity chart are deliberate, not
   thrash.
