# NG proposal — a step-decomposed, benchmark-driven algorithm lab

*Status: proposal / discussion draft (2026-07-10). Prompted by the GIAB HG002
single-sample STR experiment (`benchmarks/ssr_hg002/`), where **freebayes — a
general caller with no STR model at all — matched or beat both HipSTR and our
caller on single-sample detection**, and where we could only understand *why* by
taking the callers apart step by step. §1's step catalogue is grounded in a
source read of freebayes, GATK HaplotypeCaller (+DRAGstr), GangSTR, HipSTR, and
our own code (2026-07-10).*

## What the freebayes experiment taught us (the motivation)

Three lessons, one per idea below:

1. A **general, single-phase, unified** caller was competitive with two STR
   specialists on single-sample yield — so a simpler architecture is worth taking
   seriously (→ §3, the *ng* lab).
2. We only learned *where* each caller wins by **decomposing the algorithm into
   steps**: detection is not genotype accuracy; HipSTR was the most *reliable per
   call* (0.96 correct-given-detected) precisely because of the STR-specific steps
   freebayes lacks (→ §1).
3. We nearly drew the wrong conclusion because our **truth comparison** had a
   representation bug (STR indels are anchored one base before the tract), which
   mis-scored ~1000 correct calls as false positives — the comparison layer is
   itself a component that must be built once and correctly (→ §2).

The plan has three main ideas.

---

## Naming — "STR" for readers, `ssr` in the code

We use **STR** (short tandem repeat) as the canonical term in all prose,
documentation, CLI help, VCF fields, and papers: it is the field's convention and
matches every tool we compare against (HipSTR, GangSTR, GATK, TRTools). On first
mention, bridge both worlds — *short tandem repeats (STRs; a.k.a. SSRs /
microsatellites)*. Code identifiers keep **`ssr`** — module and type names,
variables, functions, the `ssr-catalog` / `ssr-pileup` / `ssr-call` subcommands,
and the `.ssr.psp` / `.ssr.catalog` artifacts — because `str` is a built-in Rust
type, so a `str` module or `Str` type would shadow it. In short: **STR is what you
read, `ssr` is what you type in Rust.**

---

## 1. The algorithm is not an algorithm

A caller is a *pipeline of steps*, each implementable several ways. Naming them is
not taxonomy for its own sake: it is that **each step gets a clean interface, so we
can swap one implementation and measure the effect end-to-end**, and **mix and
match steps across callers** (freebayes' candidate generation + HipSTR's stutter
likelihood + our cohort prior) to find the best implementation of each. That is
what turns "a list of steps" into a method.

### What we were missing (answer to "did we miss a step?")

Our earlier 10-step list described only the **per-locus genotyping core**. Reading
the five callers showed that every serious one wraps that core in **setup** and
**output** stages the list omitted — and those wrappers are where much of the
accuracy lives. The additions:

- **Read admission & QC** (before alignment) — MAPQ/BQ gates, duplicate removal,
  overlap-pair quality handling, downsampling, read pooling. All five have it as a
  distinct stage; we folded it invisibly into "read prep".
- **Candidate-locus / region determination** — *where do we even look?* For the
  generic path this is **active-region detection** (GATK's activity profile);
  for STR it is the **catalog/panel** (our old step 2, now seen as one branch of a
  general "which loci" step).
- **Parameter pre-pass (nuisance-parameter estimation)** — **the biggest gap, and it
  is really a *second caller*.** Base-quality recalibration (GATK BQSR), the
  error/stutter model (HipSTR & GATK-DRAGstr & ours), insert-size/coverage/GC
  profiling (GangSTR), inbreeding `F`, sample-group/chemistry clustering, and the
  contamination fraction are all *estimated in a pass of their own* and *frozen* as
  inputs to the likelihood and prior. And they are not estimated from nothing: a
  **rough first-pass caller** produces the confident genotypes the estimates are built
  on — so we are really building **two callers**, a rough one to learn the parameters
  and the real one that uses them (our SNP rough-het caller and STR confident-genotype
  gate both do this, for both paths). The 10-step core silently assumed these
  constants already existed. And the estimation is itself **two levels** — per-sample
  (per-individual params like `F`, frozen into the `.psp`) and then a distinct
  **cohort-gather** step that reads every sample's summary to estimate the panel-level
  params (sample-groups, per-group chemistry, the SFS `θ`). That cohort-gather step was
  the most-missed of all — it is a real pipeline stage between the per-sample artefacts
  and the cohort caller.
- **Phasing** — physical / read-backed (GATK), or SNP-based coupling of the STR
  allele to the read's haplotype (HipSTR). Neither our list nor our caller has it.
- **Uncertainty & emission mode** — confidence intervals and expansion
  probabilities (GangSTR), rich bias annotations for downstream filtering
  (freebayes, GATK), and variant-only vs **gVCF** reference-confidence output.

### The step map (setup → per-locus core → output)

**Setup (per-sample / per-cohort, before the locus loop):**

1. **Read admission & QC.** MAPQ/BQ filters, duplicate/secondary/QC-fail removal,
   overlapping-mate quality reconciliation, coverage downsampling, read pooling.
2. **Alignment / read preparation.** Trust the mapper, or realign — up to full
   **local haplotype assembly**; indel left-alignment / normalisation; BAQ. Many
   indels sit inside STR tracts, so STR-aware realignment must respect the
   motif/border division a generic small-indel realignment ignores.
3. **Candidate-locus / region determination.** Active-region detection (generic)
   or the STR **catalog/panel router** (motif / period / borders). This is the
   *router*: it decides whether a locus goes down the generic small-indel path or
   the STR-specialist path. **A general caller with no router runs every locus
   generic** — freebayes' whole world.
4. **Rough genotyping → parameter pre-pass, at *two levels* ("the two callers").** A
   cheap *rough* caller runs first; its **confident genotypes** are the substrate for
   the frozen parameters. The estimation happens at two levels, mirroring production:
   - **4a. Per sample.** From that sample's confident calls, estimate the params
     computable from *one individual* — chiefly its inbreeding `F` — plus the
     sufficient statistics the cohort step needs, and freeze them into the sample's
     artefact (in our `.psp`, after the genomic blocks).
   - **4b. Cohort-gather** *(the step the earlier list missed)*. Once, before calling:
     read *every* sample's summary and estimate the params that need the **whole
     panel** — sample-groups, per-group chemistry (error/stutter), the SFS `θ`,
     contamination — assembling the frozen parameter set the real caller consumes.

   We build two callers, not one.

**Per-locus core:**

5. **Read-class / spanning determination** *(STR path)*. Which reads inform the
   length, and how: span-vs-partial, or the richer enclosing / spanning / flanking
   / FRR classification with an insert-size fragment model that reaches lengths
   **beyond the read**.
6. **Candidate allele generation.** The allele set — assembly haplotypes, a
   repeat-length rung ladder, or observed sequences — plus iterative discovery and
   admission gates.
7. **Read likelihoods `P(read | allele)`.** The likelihood model (PairHMM,
   multinomial sampling, class-conditional, or stutter-aware HMM), the **PCR
   stutter** model on the STR branch, and **contamination as a mixture term**.
8. **Genotype priors.** From flat, through fixed-heterozygosity Dirichlet, to full
   SFS/Ewens integration, to a marginalised leave-one-out cohort prior — mixed
   with an inbreeding-`F` term.
9. **Posterior / inference.** ML grid search, per-sample genotype likelihoods then
   an EM allele-frequency calculation, or a full joint Bayesian marginalisation
   over the cohort's genotype space.

**Output:**

10. **Phasing.** Physical / read-backed, or SNP-based.
11. **Site filtering — two questions.** *(11a)* Is this a technical **artifact**?
    (hidden-duplication / paralog collapse, strand bias.) *(11b)* Is this real site
    **polymorphic enough to emit**? (the emission / QUAL decision.) Kept separate:
    they have different knobs and behaved independently in our STR work.
12. **Allele representation / normalisation.** Repeat-unit vs anchor-indel framing;
    left-alignment / trimming. The STR-path counterpart to step 2's generic
    normalisation — both canonicalise the change, routed by step 3.
13. **Quality, calibration & uncertainty.** GQ/PL/Q; site QUAL; calibration (BQSR
    upstream, or "the model posterior is the quality"); confidence intervals /
    expansion probabilities; bias annotations; emission mode (variant-only / gVCF).

**The STR-aware dependency chain.** Steps **3 → 5, 7 (stutter), 12** form a chain:
the motif from the router is what makes read-class logic, the stutter model, and
repeat-unit representation possible at all. That chain *is* the definition of
"STR-aware", and it is exactly what a general caller cannot retrofit without a
router — the striking exception being GATK's **DRAGstr**, which grafts the chain
(STR-context indel prior + STR-modulated PairHMM gap penalties) onto an otherwise
general assembly caller. Note the symmetry: **normalisation appears twice** (step 2
generic, step 12 STR) — the same job on each branch the router selects.

**Analyse each step, but measure end-to-end.** The steps are coupled (priors,
likelihood, and posterior are entangled; a missing candidate is unrecoverable), so
a step that looks optimal alone may not compose. Swap one implementation, hold the
rest fixed, judge by the whole pipeline against the standards (§2). And a step is
only "done" when its contract is narrow enough that a competing implementation —
or one lifted from freebayes/GATK/HipSTR — drops in. That is what makes both the
cross-caller bake-off and the eventual port-back (§3) real.

---

## Step details — the approaches, and who does what

For each step, the design axis and where each caller sits. `N/A` = the step does
not exist for that tool. Source reads: `freebayes/src/`, `gatk/.../haplotypecaller/`,
`GangSTR/src/`, `HipSTR/src/`, our `src/var_calling/` (SNP) + `src/ssr/` (STR).

### 1. Read admission & QC
*Axis:* how aggressively to filter reads, and whether to dedup / pool / downsample.
- **freebayes:** MAPQ (`--min-mapping-quality` 1, `standard-filters` → 30), BQ,
  mismatch-fraction/limit gates; duplicate drop; `--limit-coverage`/`--skip-coverage`
  random thinning.
- **GATK:** standard HC read filters (MAPQ ≥ 20, not-dup, not-secondary, good-cigar);
  `PositionalDownsampler` to 50 reads/start; adaptor + low-qual-end clipping;
  overlapping-mate quality capping.
- **GangSTR:** SSW score gate (≥ 0.75·2·readlen), `find_longest_stretch` copy floor;
  bound-read cap at 3×median; mate rescue by coordinate.
- **HipSTR:** N-base / low-sum-qual reject, end-qual trim, ±40 bp trim; **PCR-dup
  removal per library** (keep highest summed BQ); **read pooling** of identical
  sequences (HMM once per unique read).
- **ours:** MAPQ ≥ 20, min-len 30, drop dup/secondary/QC-fail; adaptor clip;
  samtools-style mate-overlap reconciliation. No pooling/downsampling.

### 2. Alignment / read preparation
*Axis:* **trust the mapper** vs **realign** vs **reassemble** — the biggest single divide.
- **freebayes:** trusts the mapper; decomposes CIGAR into allele primitives; only
  normalisation is per-read **indel left-alignment** (`stablyLeftAlign`). No assembly.
- **GATK:** distrusts the alignment inside active regions and does **local de-Bruijn
  read-threading assembly** (multi-kmer {10,25}, dangling-end recovery, adaptive
  pruning, ≤128 haplotypes), then realigns reads to the best haplotype. The signature.
- **GangSTR:** **realigns every read** — CIGAR fast-path if clean, else
  expansion-aware **Striped Smith-Waterman** against a synthetic repeat reference.
- **HipSTR:** left-aligns/realigns, trims to ±40 bp, then scores with a **read↔haplotype
  HMM** (seed-base split, flanks aligned in/out, log-sum over seed positions).
- **ours — SNP:** trusts the mapper + **left-align (GATK port)** + **BAQ** (htslib
  `probaln_glocal`). **STR:** per-read **Viterbi pair-HMM delimitation** against
  `flank+tract+flank` with a **tract-aware affine gap** (cheap gaps inside the tract).

### 3. Candidate-locus / region determination (the router)
*Axis:* discover loci from the data (active regions) vs a fixed STR panel; and
whether STR-typing routes into specialised steps.
- **freebayes:** generic detection windows only; repeat detection merely **widens the
  haplotype window** — **no STR router**, no period/motif output.
- **GATK:** **active-region / activity-profile** detection (per-base P(active),
  band-pass smoothed, threshold 0.002) → assembly regions. **DRAGstr** *does* type STR
  context (period 1–8, repeat-length 1–20 table) and feed steps 7–8 — a router graft.
- **GangSTR / HipSTR:** **panel/BED** defines every locus (motif/period); no de-novo
  discovery; reads gathered per locus by indexed query.
- **ours:** **`ssr-catalog`** (trf-mod + a GangSTR-style post-filter: period 2–6,
  purity ≥ 0.8, copy-number floors, de-bundling) is the STR router; the SNP path has
  no router (every position is a candidate, gated later).

**Data-driven STR discovery vs the catalog — an open question.**
The catalog is reference-based: `ssr-catalog` finds tandem structure in the
*reference* (reliable, TRF-style) and freezes the locus set before Stage 1. Its
cost is coverage. Measured against GIAB HG002 truth, our genome-wide catalog covers
**~65% of period-2..6 STR length-variants (1,412 / 2,175); ~35% are missed** —
largely loci dropped by our construction filters (copy-number floors, purity ≥ 0.8)
or where the reference tract is short/absent but the sample expanded. *(Caveats:
we should check with measuments upto which period length its better to go through the STR route; some misses —
e.g. a +716 bp `AC` expansion — are unrecoverable from 148 bp reads regardless. The
earlier freebayes benchmark did **not** measure this — it scored all callers on the
same catalog, so its recall gaps were genotyping, not coverage; this ~35% is a
separate direct measurement.)*

The tempting fix is **data-driven discovery** — route a locus into the STR path from
the *reads*, not just the reference catalog, recovering loci the catalog missed. Two
objections make the bar for adopting it high:

1. **Detecting STR-ness reliably from reads might be hard — and it is the *discovery* half,
   not the *typing* half.** Typing an already-known locus is reliable because it reads
   the tandem structure off the reference (what the catalog does). Discovering a locus
   the catalog missed means recognising a real repeat-unit length change in noisy reads
   and telling it apart from PCR **stutter** — exactly the signal that is weakest in a
   single sample.
2. **Good discovery signal needs the cohort, which fights the `.psp` architecture.**
   Real-variant-vs-stutter is a *cohort* signal (a real allele segregates; stutter is a
   consistent shoulder). But locus determination happens *before* Stage 1, so the
   per-sample pileup knows what to genotype. Cohort-dependent discovery is
   chicken-and-egg: you need the cohort to choose the loci, but Stage 1 is per-sample.
   Resolving it costs the clean per-sample → `.psp` → cohort flow — either a **cohort
   discovery pre-pass** over all BAMs that augments the catalog before pileup (an extra
   cohort scan), or **per-sample over-collection + union** (which breaks the "every
   locus genotyped in every sample" property and drags in GVCF-style re-visiting).

**Decision rule.** Adopt data-driven cohort discovery **only if it proves to beat the
catalog** on the standards (§2) — and only after the cheap, architecture-neutral option
is exhausted first: **better reference-based catalog construction** (relax the
copy-number floors and purity gate, include shorter reference tracts, the
interrupted-repeat work). Much of the measured ~35% gap is likely at loci where the
reference *does* carry a short/impure tract we filtered out — recoverable purely
reference-side, no cohort, no `.psp` change. Cohort discovery is only needed for the
residual where the reference has *no* tract: a smaller, intrinsically harder set. So
the sequencing is: (i) split the ~35% gap into *relax-construction-recoverable* vs
*needs-discovery* vs *unrecoverable-from-short-reads*; (ii) close the first bucket
cheaply; (iii) reach for cohort discovery only if the residual is large enough to
justify breaking the per-sample architecture **and** it measurably wins.

### 4. Rough genotyping → parameter pre-pass (the "two callers")
*Axis:* which nuisance parameters are learned in a dedicated pass — and the point we
had flattened: estimating them **needs genotypes**, so a *rough first-pass caller*
runs first and its **confident calls** are the substrate. **We are not building one
caller but two:** a cheap rough caller to learn the parameters, then the real caller
that uses them. This breaks the chicken-and-egg (you need parameters to genotype, but
genotypes to estimate parameters).
- **freebayes:** no learned parameters, so no rough pass (priors are analytic;
  contamination supplied, not estimated).
- **GATK:** rough passes feed calibration — **BQSR** (empirical base-quality
  recalibration by covariates) and **DRAGstr autocalibration** (`CalibrateDragstrModel`:
  per-genome STR gap-penalty/API table by ML over sampled STR loci), both upstream.
- **GangSTR:** a **sample-profile** pass — empirical insert-size PDF/CDF, coverage, GC
  bins, read length — swept from off-target regions before genotyping.
- **HipSTR:** a **length-based EM stutter pre-pass** (a rough length-genotyper) whose
  posteriors fit the per-locus stutter model before the sequence genotyper runs.
- **ours — the clearest two-caller split, and it works for BOTH paths:**
  - *SNP:* a **rough per-site het caller** in Stage-1 pileup (three binomial models
    hom-ref / het / hom-alt with a crude base-quality ε̂) produces per-sample `obs_het`,
    from which the per-individual inbreeding **`F`** is derived — all inside the
    per-sample `.psp`.
  - *STR:* a **confident-genotype gate** (a BIC 1-vs-2-allele test) over a 20k-locus
    cohort burn-in yields each sample's confident genotypes; from them the pre-pass
    freezes per-group ε, stutter shape/level, the `G₀` spread, sample-group clusters,
    and per-individual `F`.

**Considerations — the two-caller pattern.**
- The rough caller must be robust *without* the parameters it is helping estimate
  (ours uses a crude base-quality ε̂, not the fitted ε) — otherwise the bootstrap is
  circular.
- **Only confident calls may teach the parameters.** The confident-subset gate (our
  het margin / the STR 1-vs-2 BIC test) is essential; feeding shaky genotypes
  contaminates ε / stutter / `F`. That is why the confident-genotype gate is a step,
  not a knob.
- It maps onto our two-phase architecture asymmetrically: the SNP rough caller lives
  in **Stage-1 (per-sample → `.psp`)**, the STR one in a **Stage-2 cohort burn-in**.
  `F` is per-individual (fits the per-sample model); the chemistry is pooled (needs
  the cohort).
- **The pre-pass is two levels, and the cohort-gather is a distinct step** (the one we
  had missed). Per-sample estimation (per-individual `F` + sufficient stats) writes into
  the `.psp`; then a **cohort-gather** reads every sample's summary and estimates the
  panel-level params (sample-groups, per-group chemistry, SFS `θ`). The interfaces model
  this as `SampleSummarizer` (per sample) + `CohortEstimator` (the gather) —
  [`../arch/ng_step_interfaces.md`](../arch/ng_step_interfaces.md) §4.
- **For ng (§3): "single-phase" drops the `.psp` per-sample/cohort split, NOT the
  rough→real bootstrap.** The bootstrap is a within-run two-pass, orthogonal to the
  `.psp` architecture — ng keeps it.

### 5. Read-class / spanning determination (STR path)
*Axis:* binary span-vs-partial, or a multi-class scheme that reaches beyond read length.
- **GangSTR:** four classes — **enclosing** (exact count), **spanning** (insert-size),
  **flanking** (lower bound), **FRR** (Poisson count for expansions > read) — the
  defining idea; off-target FRR harvesting for large expansions.
- **HipSTR:** spanning not required by default (all overlapping reads scored);
  candidate-*generating* reads must span with ≥5 bp clean flanks.
- **ours:** reach gate (≥5 bp flank), require both flank junctions inside the read;
  **widen** the window to recover mapper-collapsed long alleles, else drop.
- **freebayes / GATK:** `N/A` (generic partial-observation support instead).

### 6. Candidate allele generation
*Axis:* assembly haplotypes vs repeat-length ladder vs observed sequences; iterative discovery.
- **freebayes:** haplotype window grown to a fixpoint; alleles admitted at
  min-count 2 / min-fraction 0.05 in some sample.
- **GATK:** events read off assembled haplotypes' CIGARs (`EventMap`); allele-based
  pre-genotyping pruning.
- **GangSTR:** integer repeat-length **grid** (bounds inferred per locus), exhaustive
  or NLopt-continuous.
- **HipSTR:** repeat sequences observed in spanning reads + **stutter-hidden alleles**
  + **flank de-Bruijn assembly** alleles, discovered iteratively.
- **ours — SNP:** overlap-bundled positions, byte-dedup, ≤6 alleles + `<OTHER>`.
  **STR:** length-keyed **rung ladder** + same-length interruption alleles behind a
  3-filter gate; admission gates (LowDepth / NotPeriodic / TooManyAlleles).

### 7. Read likelihoods `P(read | allele)` — incl. stutter & contamination
*Axis:* the emission model, whether PCR stutter is modelled, and how contamination enters.
- **freebayes:** **multinomial sampling** model over allele frequencies (no stutter);
  contamination = **per-read-group mixture**.
- **GATK:** **PairHMM** (match/insert/delete, gap-open/continue penalties);
  **DRAGstr modulates gap penalties by STR period/repeat-length** so stutter indels
  aren't over-penalised; contamination = **read downsampling**, not a mixture.
- **GangSTR:** **class-conditional** likelihoods (enclosing stutter geometric,
  spanning/FRR insert-size, flanking lower-bound) + FRR Poisson count; no contamination.
- **HipSTR:** read↔haplotype HMM with a **two-component geometric stutter model**
  (in-frame in units, out-of-frame in bp), summed over artifact sizes; no contamination.
- **ours — STR:** **Model A (HipSTR-style)** — stutter PMF (two geometrics) × substitution
  align, uniform outlier floor. **SNP:** base-quality multinomial + optional
  **contamination mixture**.

### 8. Genotype priors *(the user's example step)*
*Axis:* flat → fixed-het Dirichlet → SFS/Ewens integration → marginalised LOO cohort prior; ± inbreeding `F`.
- **GangSTR / HipSTR-seq:** **flat / uniform** over genotypes (pure ML). HipSTR's
  length pre-pass can use population allele frequencies.
- **GATK:** **fixed heterozygosity** (SNP 1e-3, indel 1.25e-4) → Dirichlet
  pseudocounts (a Watterson-θ-implicit prior); DRAGstr swaps the indel pseudocount
  for the STR-context value. **No inbreeding `F`.**
- **freebayes:** **SFS integration via the Ewens Sampling Formula** (analytic,
  `--theta`) + HWE multinomial prior. **No `F`.**
- **ours:** **marginalised leave-one-out empirical-Bayes Dirichlet-multinomial SFS**
  mixed with a **Wright inbreeding-`F`** term — the cohort borrows strength across
  samples; STR default is plug-in HWE-Wright with an opt-in marginalised prior.
  Ours is the only one with `F`; freebayes and ours are the two with real SFS priors.

### 9. Posterior / inference
*Axis:* ML point estimate vs per-sample-GL-then-EM-frequency vs full joint marginalisation.
- **GangSTR:** **ML grid / COBYLA** over the two allele lengths (flat prior).
- **HipSTR:** per-sample diploid **posterior over haplotype pairs**, coupling each
  allele to the read's **SNP phase** (`log_p1/log_p2`); MAP call.
- **GATK:** per-sample genotype likelihoods → **EM fixed-point allele-frequency**
  calculation (flat first iteration to avoid a variant-biased minimum) → P(AC=0).
- **freebayes:** **full Bayesian marginalisation over the joint cohort genotype
  space** (convergent banded combo search); QUAL = Phred P(not-variant).
- **ours:** per-record / per-locus **cohort EM** with the LOO Wright-`F` prior,
  converging on pseudocount-free expected counts; byte-identical across threads.
  (GATK's EM-AF and ours are close cousins; freebayes marginalises the whole joint.)

### 10. Phasing
*Axis:* none, physical/read-backed, or SNP-based coupling.
- **GATK:** **read-backed physical phasing** (`constructPhaseGroups` → PGT/PID).
- **HipSTR:** **SNP-based physical phasing** feeds the genotype posterior directly.
- **freebayes / GangSTR / ours:** `ABSENT` (a genuine gap for us).

### 11. Site filtering — 11a artifact, 11b emission
*Axis:* inline hard-drop vs annotations-for-downstream; and how the emit decision is made.
- **11a artifact:** **freebayes / GATK / GangSTR / HipSTR** emit **bias annotations**
  (strand, position, MQ) and defer filtering to downstream tools (freebayes →
  `vcffilter`; GATK → **VQSR / CNN**; GangSTR → `dumpSTR`). **ours** is the outlier:
  an **inline hidden-paralog likelihood-ratio filter** (+ DUST) that **hard-drops**.
- **11b emission:** everyone thresholds a site score — freebayes `pVar ≥ --pvar`,
  GATK `QUAL ≥ 30`, ours SNP min-QUAL 30. **ours STR** has a dedicated, swappable
  **emission model** (`heuristic | bic | freebayes`) — the bake-off subject.

### 12. Allele representation / normalisation
*Axis:* repeat-unit vs anchor-indel framing; where left-align/trim happens.
- **freebayes / GATK / ours-SNP:** anchor-indel VCF, **left-aligned & trimmed**
  (inline for GATK/freebayes; upstream in the merger for us).
- **GangSTR:** **repeat-copy-number** REF/ALT + `REPCN`.
- **HipSTR:** sequence-resolved REF/ALT + **`GB`** (bp diff per allele) + `BPDIFFS`.
- **ours-STR:** repeat-unit tract sequences; falls back to anchor-indel only for a
  full-tract deletion. *(This is the framing whose mishandling cost ~1000 mis-scored
  calls in the benchmark — see §2.)*

### 13. Quality, calibration & uncertainty
*Axis:* is quality the raw model posterior, or separately calibrated? what uncertainty is reported?
- **freebayes / HipSTR / GangSTR:** quality **is** the model posterior (GQ/Q); **no
  post-hoc recalibration**. GangSTR adds **bootstrap confidence intervals** (`REPCI`,
  `STDERR`) and an **expansion probability** (`QEXP`).
- **GATK:** PL/GQ from the model, but **calibration lives upstream in BQSR**;
  optional gVCF reference-confidence GQ blocks.
- **ours — SNP:** GQ + site QUAL from the **exact allele-count posterior**, then a
  **`final_qual` refinement** subtracting allele-balance/strand/position bias
  penalties. **STR:** posterior GQ + emission-model QUAL. No CIs / expansion prob.

### Cross-cutting reading (the levers ng should probe)
- **Trust vs reassemble (step 2)** is the deepest axis and the leading hypothesis for
  freebayes' detection edge; GATK's assembly and per-read STR realignment are the two
  ends we should A/B on synthetic data.
- **The router (step 3) with a DRAGstr-style graft** is the cheapest way to add
  STR-awareness to a general core — a concrete design for ng.
- **Priors (step 8)** span a clean ladder — flat → fixed-het → SFS → marginalised-LOO
  — and are individually swappable; this is the textbook per-step bake-off.
- **The parameter pre-pass (step 4)** is where our caller is richest and the general
  callers thinnest; whether that pays off is directly testable.
- **Filtering philosophy (step 11a)** splits the field: annotate-and-defer (everyone
  else) vs inline-hard-drop (us). Worth a deliberate choice, not an accident.

---

## 2. Benchmarks

To optimise the algorithms we need gold and silver standards, and we can check each
algorithmic step against all of them to learn which decisions are optimal in each
case and which are compromises.

- **Gold — GIAB (HG002).** Assembly-based truth, orthogonal to short reads. But it is
  **single-sample and high-depth**: it cannot exercise the cohort machinery or the
  low-depth lone-carrier-het regime (see `benchmarks/ssr_hg002/README.txt`).
- **Silver — tomato.** A real cohort, but **not truth** (read-grounded, partly
  circular; HipSTR-biased).
- **Synthetic — the bridge.** Each real standard covers only half the space, so
  synthetic data is not just for fast regression — it is the **only** way to (a) test
  a step **in isolation** (you injected the candidates / stutter / genotypes, so you
  *know* the answer — no truth-matching ambiguity) and (b) get a **gold *cohort*** —
  truth for the multi-sample, low-depth regime that no real dataset gives us and that
  the whole STR optimisation is stuck on.

Two concrete build items fall out of the freebayes experiment:

- **A shared, representation-normalising truth-comparison layer**, built once. Every
  ad-hoc benchmark script otherwise re-derives the anchor / normalisation bugs (that
  is how we mis-scored ~1000 calls). This is vcfeval / hap.py-style matching —
  normalise both sides, compare by haplotype and by length — worth adopting or porting
  rather than hand-rolling per script.
- **A synthetic variant generator** with known genotypes and realistic stutter/error,
  usable both for per-step isolation and as a fast regression gate.

---

## 3. ng

The current architecture is powerful but complex — the per-sample pileup → `.psp` →
cohort split, the producer/merger, the columnar streaming, the memory-scaling work —
and changing the *algorithm* inside it is laborious. So to create this
new-generation (ng) algorithm we build a **simpler caller from scratch: one thread,
a single phase** (not the per-individual / cohort split). This makes us nimble.

This does not mean writing everything from scratch. We build it as a **new module
inside `pop_var_caller`** and reuse anything useful. And — the load-bearing
constraint — **ng is scoped explicitly as an algorithm research vehicle,
single-phase *because it is a prototype*, with the understanding that the port-back
step is where scaling is re-added.**

- **The two-phase split is the scaling story, not incidental complexity.** `.psp`
  means each sample is piled up **once**, so samples can be added cheaply; a
  freebayes-style all-in-memory cohort model **does not scale** to hundreds or
  thousands of samples (freebayes/GATK OOM on large cohorts — that memory efficiency
  is one of this project's headline advantages). ng deliberately sets scaling aside.
- **The port-back step is where scaling is re-added.** The final objective, if the ng
  exercise proves to work, is to plug the winning algorithm/steps back into the
  `pop_var_caller` architecture — each behind the same interface (§1) — and take
  advantage of that work. ng feeds the production engine; it does not replace it.

**Success gate (so ng cannot become an endless side-quest).** On the standards (§2),
does ng's per-step quality match or beat the current caller, telling us *which
implementation to port back*? A spike that cannot answer that is a distraction.

**Where to start.**

- **STR first** — the smallest surface, the freshest insight (the freebayes surprise),
  and the sharpest single-sample-vs-cohort tension.
- **First probe: candidate generation + realignment (step 2/6).** Local reassembly /
  haplotype candidate generation (freebayes / GATK-style) versus our catalog-anchored
  per-locus evidence — the leading hypothesis for freebayes' detection edge, cleanly
  testable in isolation on synthetic data.

---

## Summary

Decompose the caller into steps with clean interfaces (§1) — a **setup / core /
output** pipeline with an explicit **STR-locus router** (and a DRAGstr-style graft as
the cheap way to add STR-awareness to a general core), the previously-missing
**parameter pre-pass**, and normalisation on both branches the router selects; measure
every step against gold, silver, and — crucially — **synthetic** standards through a
representation-correct comparison layer built once (§2); and iterate on the algorithm
in a nimble, single-phase **lab** with a hard success gate, then **port the winning
steps into the scalable architecture** (§3). The current architecture is not the thing
to escape — it is the thing the lab eventually feeds.
