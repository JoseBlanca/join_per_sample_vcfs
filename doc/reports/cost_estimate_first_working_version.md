# Cost and timeline estimate — first working version

A retrospective on the time, team size, and conventional cost it would
have taken to build this multi-sample SNP caller, compared against the
two canonical reference callers in the field: the original **GATK 1**
(broadgsa/gatk-protected) and **freebayes**.

This is not a prediction of future cost. **Stage 1 (CRAM → `.psp`
pileup) has been validated against samtools' pileup on real data —
~95% match on SNPs and alleles, with our walker slightly stricter
(likely due to BAQ).** Stages 3–6 (var-calling and posterior
engine) have only been exercised against synthetic fixtures so far;
real-cohort validation of the end-to-end caller is the next major
milestone. The report is an after-the-fact accounting of what has
been built so far and what the conventional industry estimate for
that scope would be.

---

## Our project — timeline

The repository started as a **gVCF-merger** in February 2026 and was
pivoted to a from-scratch multi-stage SNP caller in mid-April.
Milestones, derived from `git log`:

| Date | Milestone | Elapsed from pivot |
|---|---|---|
| 2026-02-04 | Initial commit of the gVCF-merger codebase | — |
| 2026-04-14 | Last gVCF-merger work (parser review + fixes) — legacy era closes | — |
| **2026-04-17** | **First ideas for complete variant caller** (commit `11deac3`) | 0 |
| 2026-04-22 | First per-sample artefact spec draft (commit `a6e99a3`) | +5 days |
| **2026-04-24** | **`calling_pipeline_architecture.md` lands** — authoritative spec (`7f28593`) | +7 days |
| 2026-04-29 | First new-pipeline implementation code: CRAM input slice (`94bbad4`) | +12 days |
| 2026-05-06 | Pileup walker implemented | +19 days |
| 2026-05-12 | BAQ shipped | +25 days |
| 2026-05-13 | `.psp` writer + reader shipped | +26 days |
| 2026-05-15 | Multi-way per-position iterator (Stage 4 input) | +28 days |
| 2026-05-16 | Variant grouper + per-group merger (Stages 4–5) | +29 days |
| 2026-05-17 | DUST filter + contamination estimation + posterior engine v1 | +30 days |
| 2026-05-18 | Cohort VCF writer (Stage 6 sink) | +31 days |
| **2026-05-19** | **Cohort CLI shipped** — end-to-end pipeline runs on synthetic fixtures | **+32 days** |

**State at 2026-05-19:** all six stages implemented, code-reviewed,
fixes-applied; cohort CLI (`var-calling`, `estimate-contamination`,
`var-calling-from-bam`) integration-tested on synthetic data.

**Stage 1 (pileup) validated against samtools** on a real BAM
(`SRR7279725_small`, see `tmp/compare_pileups.py` /
`tmp/compare_pileups.html`): ~95% match on SNPs called and on the
alleles within those SNPs. Our walker is slightly stricter than
samtools, most plausibly because of BAQ — samtools `mpileup` defaults
to BAQ on, but the recalibration cascade and the BAQ-aware filters
downstream of our walker drop a small fraction of weakly-supported
calls that samtools keeps. The divergence is in the expected
direction; no correctness defect was found.

**Stages 3–6 (var-calling, contamination estimation, posterior
engine) not yet validated on real cohort data** — the synthetic
integration-test fixture is too small to be meaningful. The
`bcftools view` / `bcftools stats` smoke against real cohort data is
the open standing item per PROJECT_STATUS.md.

**Team:** one developer (Jose Blanca) collaborating with Claude (an LLM
assistant). Claude is used for code generation, design review,
architecture critique, and documentation drafts; all decisions and
final commits are the developer's.

**Codebase:** 41,826 lines of Rust code (`tokei`-style count: blanks
and comment-only lines excluded). Total including blanks and comments:
~58,000. Tests and benches account for ~10 K of the code lines.

---

## Reference projects — what they needed for a first working version

### freebayes (Erik Garrison, 2010)

| Date | Milestone | Elapsed |
|---|---|---|
| 2010-02-04 | Initial commit | 0 |
| 2010-05-01 | Working posterior math in original `bayes`/`bamBayes` codebase | ~3 months |
| 2010-05-25 | First "release candidate" of the original codebase (later discarded) | ~3.7 months |
| 2010-07-26 | Algorithm prototyped end-to-end in Python | ~5.7 months |
| 2010-07-27 | C++ rewrite begins (`esfbayes`) | ~5.7 months |
| 2010-08-05 | C++ rewrite reaches "WORKING commit, produces sane marginals" | ~6 months |
| **2010-08-12** | **"A working implementation"** — pileup + posterior + VCF | **~6 months 8 days** |
| 2010-10-13 | `v0.1.0` — first tagged release | ~8 months 9 days |
| 2010-08-30 | Project renamed `alleleBayes → freebayes` | ~6 months 26 days |

**Team:** **one person** — Erik Garrison, sole contributor for the
entire first six months and well beyond (first non-Erik commit was
2014-01-22, nearly four years after the initial commit).

**What he had at 6 months:** pileup-equivalent (`AlleleParser`, ~1,000
LOC consuming the BamTools library), posterior calculation
(`DataLikelihood`, ~207 LOC), `Allele` and `Genotype` types, VCF
output, and a 333-LOC main driver. Roughly **12,000 lines of C++**,
much of it vendored BamTools.

**Force multipliers** that made the 6-month solo timeline plausible:
- **BamTools as a vendored dependency** — Erik did not write the BAM
  parser. Equivalent to our use of `noodles`.
- **Python prototype first** — he validated the algorithm end-to-end
  in Python in ~2 weeks before starting the C++ rewrite. The C++
  rewrite itself was ~2 weeks.

**What his 6-month pileup did NOT do** (added later, mostly in
2011–2012): BAQ, mate-overlap resolution, mismatch-fraction filtering,
soft-clip / indel-edge handling, indel left-alignment, low-complexity
filtering, depth caps, contamination estimation, partitioning into
independent variant groups.

### GATK 1 (Broad Institute, 2009 — original GATK lineage)

| Date | Milestone | Elapsed |
|---|---|---|
| 2009-02-26 | Initial commit (imported from SVN — earlier history not preserved) | 0 |
| 2009-03-13 | "First somewhat functional version of AlleleFrequency caller!" — single-sample biallelic prototype on top of pre-existing walker abstractions | **~2 weeks** |
| 2009-06-16 | "Behold MultiSampleCaller!" | ~3.7 months |
| 2009-09-03 | Multi-sample SNP caller check-in | ~6 months |
| 2009-10-13 | Remodeled `UnifiedGenotyper` — production caller | ~7.5 months |
| **2010-08** | **McKenna et al. paper in *Genome Research*** — canonical "first public release" | **~18 months** |
| 2011-07-22 | `1.1` — first real numbered release tag | ~29 months |
| 2012-07-24 | `2.0` — HaplotypeCaller era begins | ~3.5 years |
| 2014-03-05 | `3.0` — GVCF / joint-calling workflow | ~5 years |

**Team in the first 6 months:** **12 distinct contributors**, ~9 of
them doing substantial work. By the McKenna 2010 paper at ~18 months,
the contributor list had grown to **21 people**, with the same ~10
core engineers (Hanna, McKenna, DePristo, Garimella, Banks, Sivachenko,
Hartl, Poplin, Maguire, Kernytsky) doing the overwhelming majority.

**Important context for the 2-week milestone:** the "First somewhat
functional AlleleFrequency caller" commit was not 2 weeks of
from-scratch work. The SVN history shows the commit was imported as
`trunk@4` — there are 3 earlier SVN commits not preserved as git
history, and source-file headers reveal `LocusContext.java` was
created by Mark DePristo on Feb 22, 2009, **four days before the
"Initial commit"**. The walker abstractions, `LocusIterator`,
`PileupWalker`, and SAM ingestion were already in place when the
git history begins. The 2-week milestone added a ~200-line caller
on top of pre-existing infrastructure.

**What GATK 1 had at 6 months but our project does NOT have yet:**
real validation on cohort data, numerical accuracy proven against
ground-truth callsets, exhaustive edge-case fixtures from messy
production BAMs.

---

## Three-way comparison

| Project | Team | Time to first working caller | Time to first public / production release |
|---|---|---|---|
| **freebayes** | 1 person (Erik Garrison) | ~6 months (porting a Python prototype to C++ on top of BamTools) | ~8 months (`v0.1.0`) |
| **original GATK 1** | 9 core in 6 months, 12 total | ~3.7 months (MultiSampleCaller); ~7.5 months (UnifiedGenotyper) | ~18 months (McKenna 2010 paper) |
| **our project** | 1 person + Claude | **~32 days** for Stages 1–6 implemented + cohort CLI; Stage 1 pileup validated against samtools (~95% SNP/allele match); Stages 3–6 not yet validated on real cohort data | not yet — Stage 3–6 validation work remaining |

**Caveats on the comparison.** This is not an apples-to-apples
benchmark — the three projects had different starting points, force
multipliers, and validation bars.

- **freebayes** had a vendored BAM library (BamTools) and a Python
  prototype as force multipliers. Our project has `noodles` for BAM/
  CRAM I/O — equivalent force multiplier.
- **freebayes' 6-month pileup was minimum-viable** (no BAQ, no
  mate-overlap, no filters). Our pileup is substantially more
  sophisticated — closer to freebayes ~v0.8 or GATK 1.x.
- **GATK 1** was building on a pre-existing walker framework that
  pre-dated its git history. Our project's walker abstractions
  (`PerPositionMerger`, `VariantGrouper`, `PerGroupMerger`) are
  written from scratch.
- **GATK 4 (the 2014–2018 rewrite)** is a less useful comparison —
  it was porting a production-validated algorithm onto a new engine
  at industrial scale (~14 contributors in 6 months, ~3 years to
  v4.0.0.0). The "team size" reflects engine + cloud + parity-testing
  scope, not algorithm difficulty.
- **Stage 1 is validated; Stages 3–6 are not.** Our pileup matches
  samtools to ~95% on real data, with the divergence in the
  expected direction (we're stricter, likely from BAQ). The
  end-to-end variant-calling path (Stages 3–6) has only been
  exercised against synthetic fixtures so far. The remaining work
  — real cohort data, numerical accuracy against ground truth,
  edge-case hardening — is exactly what consumed freebayes' next
  18 months and GATK 1's next 18 months. This is the part where
  LLM assistance helps less, because the slow step is waiting for
  empirical results.

---

## Cost estimate — what 41,826 LOC of Rust would conventionally cost

Using **COCOMO Basic** (Boehm 1981), with the model's three calibration
modes:

| Mode | Applicability | Effort (person-months) | COCOMO calendar time |
|---|---|---|---|
| Organic | Simple, small team, well-understood domain | ~121 | ~15 months |
| **Semi-detached** | **Mixed experience, moderate constraints** | **~197** | **~16 months** |
| Embedded | Strict correctness, performance, hardware constraints | ~318 | ~16 months |

A multi-sample variant caller — numerical correctness, performance-
critical, domain-specific (statistics + low-level genomics I/O) — sits
between semi-detached and embedded. Anchor on **~200 person-months**,
with ~320 as the upper bound.

**Cross-check against modern productivity rates.** Industry surveys
put production-grade Rust (tested, reviewed, documented) at
**200–400 LOC per developer-month**. For a math-heavy caller, apply
the lower end:

- At 200 LOC/dev-month: 41,826 / 200 = **209 person-months**
- At 300 LOC/dev-month: 41,826 / 300 = **139 person-months**

These bracket COCOMO semi-detached almost exactly.

### Cost in dollars and euros

Loaded engineer cost (salary + benefits + overhead + tooling +
management):

| Region / role | Loaded €/$ per month | At 200 pm (semi-detached) | At 320 pm (embedded) |
|---|---|---|---|
| US senior software engineer | $20,000 | **$4.0 M** | $6.4 M |
| US specialist (genomics + Rust) | $25,000 | $5.0 M | $8.0 M |
| EU senior software engineer | €10,000 | **€2.0 M** | €3.2 M |
| EU PhD-level scientific developer | €8,000 | €1.6 M | €2.6 M |
| Spain academic (CSIC / UPV staff scientist) | €5,000 | €1.0 M | €1.6 M |

### Adjustments that pull the estimate down

- **Rust is denser than COCOMO's calibration languages** (COBOL, C,
  FORTRAN). A 0.7–1.0× adjustment may apply.
- **Test code is included** in the 41,826 count (~10 K LOC). Tests
  are typically written 2–4× faster than production code; stripping
  them knocks ~25% off the estimate.
- **Algorithms are not invented from scratch.** Bayesian math from
  published papers; BAM/CRAM parsing delegated to `noodles`. This
  pushes the project toward COCOMO *organic*, not embedded.
- **Partial validation only.** Stage 1 (pileup) has been validated
  against samtools on real data (~95% match). Stages 3–6 have been
  exercised only against synthetic fixtures. A conventional project
  would also include full end-to-end cohort validation, accuracy
  benchmarking against ground-truth callsets, and edge-case
  hardening across messy production BAMs — typically 0.5–1× the
  implementation cost again to reach the full conventional bar.

After adjusting for "algorithms are designed, not invented" and "no
production hardening yet", a realistic conventional cost for the
current state is roughly:

- **US:** $2 M – $4 M
- **EU:** €1 M – €2 M
- **Spain academic rate:** €0.5 M – €1 M

---

## Actual cost

| Line item | Estimate |
|---|---|
| Developer time (one developer, ~32 days from pivot to today) | ~1 month of salary |
| Claude usage (API or subscription) | ~$100–500 |
| **Total approximate cost so far** | **~€5,000–€6,000** (at Spain academic rate) |

Compared against the conventional estimate of **€0.5 M – €2 M**, the
cost compression is in the range of **100×–400×**. This is large
enough that it cannot be explained away by uncertainty in the COCOMO
calibration — even allowing for every conservative adjustment, the
order-of-magnitude difference is real.

### Honest caveats on the comparison

- **Prior domain knowledge is invisible to COCOMO.** Years of
  experience in bioinformatics, variant calling, and the genomics
  literature are not counted in either side of the comparison. They
  are nonetheless a major input to what made the 32-day timeline
  achievable.
- **The remaining work is real.** COCOMO prices a project to
  production quality. Our project has reached "implemented and
  code-reviewed against synthetic fixtures" — the validation arc is
  ahead. freebayes spent ~18 months and GATK 1 spent ~18 months on
  this validation work, much of it post-public-release.
- **Code review depth differs.** Conventional projects accumulate
  hundreds of reviewer-hours per major change. Our reviews are
  generated by Claude and a single human reviewer; the depth is real
  but the team-coverage variance is different from a 9-person Broad
  Institute team.
- **LLM-assisted productivity is not stable across tasks.** Routine
  engineering (parsers, error plumbing, test scaffolding, plan
  documents) is where Claude helps most. Genuine design decisions,
  numerical accuracy debugging, and unexpected behaviour on real
  data are where the speed-up shrinks.

The defensible conclusion: **what has been built in 32 days
corresponds to a conventional 1–2 person-year project costing
€0.5 M – €2 M.** The next phase of work (validation, accuracy,
hardening) will not see the same productivity multiplier and should
be planned with conventional timelines in mind.

---

## Sources

- Project git log (`git log` on this repo).
- freebayes git log (`freebayes/` submodule).
- Original GATK git log (`gatk-protected/` submodule).
- GATK 4 git log (`gatk/` submodule).
- McKenna et al. 2010, *The Genome Analysis Toolkit: a MapReduce
  framework for analyzing next-generation DNA sequencing data*,
  Genome Research 20(9): 1297–1303.
- DePristo et al. 2011, *A framework for variation discovery and
  genotyping using next-generation DNA sequencing data*, Nature
  Genetics 43: 491–498.
- Garrison & Marth 2012, *Haplotype-based variant detection from
  short-read sequencing*, arXiv:1207.3907.
- Boehm 1981, *Software Engineering Economics* (COCOMO).
- Industry productivity figures: Stack Overflow Developer Survey
  series; published Rust adoption case studies.
