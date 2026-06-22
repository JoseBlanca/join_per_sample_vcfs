# SSR cohort calling — Mark 2 (`ssr-call`, empirical candidates)

**Status:** draft, 2026-06-19, branch `ssr-cohort`. Built **from scratch, one
agreed premise at a time** (the way [ssr_ladder_model.md](../architecture/ssr_ladder_model.md)
and the other Mark-2 docs were grown). Sections firm up only once agreed; §9 is the
live agenda of open premises. This is the **statistics / model-intent** spec for
Stage 2 — *what* the cohort caller computes and *why*; the module/struct/signature
shape (the *how*) is deferred to a later architecture doc. Reading/orchestration
*intent* is settled in §4.1 (signature-level shape still deferred); §4.2 frames the
whole model the open sections detail, and §4.4 settles how the shared parameters are
estimated (a pre-pass freezes `ε` and seeds a per-loci-group `θ` prior refined per locus;
`F` a prior-side loop). All firmed 2026-06-19.

This is the **Mark-2** cohort spec. It is the Stage-2 companion of the Mark-2 model
([ssr_ladder_model.md](../architecture/ssr_ladder_model.md)) and Stage-1 pileup
([ssr_pileup_mark2.md](../architecture/ssr_pileup_mark2.md), as-built in
[src/ssr/pileup/](../../../src/ssr/pileup/)). It **supersedes §5 (and the §4.2/§4.3
`Qᵣ` it depended on) of the Mark-1 spec** [ssr_genotyping.md](ssr_genotyping.md) —
but it **reuses most of that §5's machinery verbatim** (§3 below). Where the two
disagree, this doc wins for Mark-2; `ssr_genotyping.md` remains the reference for the
math each piece implements, and is owed an amendment once this settles.

---

## 1. What changed from Mark-1, in one statement

> **Candidate alleles are *observed sequences*, assembled in Stage 2 from the pooled
> per-sample sequence distributions — not reference rungs, and not anything Stage 1
> scored.** Stage 1 now hands Stage 2 *raw evidence* (distinct repeat-region
> sequences + counts, uniform-quality), and **all likelihood computation moves into
> Stage 2.**

Two consequences set the whole shape of this spec:

1. **The read likelihood `Qᵣ` is computed *here*, not in Stage 1.** Mark-1's §4.2
   two-tier per-Q `Qᵣ` (on-ladder + off-ladder, computed per read during pileup) is
   gone. Mark-2 Stage 1 quality-gated the reads and stored bare sequences+counts
   (uniform-quality, [ssr_ladder_model.md](../architecture/ssr_ladder_model.md) §4),
   so Stage 2 scores each *distinct observed sequence* against each *candidate* with a
   **HipSTR-informed, flat-quality** likelihood that sums over the candidate's PCR slips
   (§6 / S3). Identical observed sequences collapse across the whole cohort, and —
   because the per-base error `ε` is **estimated once and frozen** before genotyping
   (§4.4) — the `align(obs | cand ⊕ Δ)` scores are **invariant across EM iterations** (only
   `S_θ(Δ)` re-weights them, `θ` refined per locus — a handful of numbers, not a recompute).
   That makes them **cacheable in principle**; whether to *build* a cache is a separate,
   deferred performance question (§6: post-C1 the in-tract score is a cheap substitution
   closed-form, so v1 **recomputes on demand** and we measure before caching).

2. **The candidate set is assembled from sequences, and "on/off-ladder" is gone as a
   label.** Mark-1 assembled `A_ℓ` from the cohort-aggregate integer-length profile
   (local-maximum rungs + ±1-adjacent) ∪ off-ladder-by-key. Mark-2 assembles it from
   the pooled distribution over *observed sequences* (§5 / S1); the on/off-ladder
   distinction relocates to a **relationship between candidate sequences** — "is B a
   whole-unit stutter product of A?" — used only by the stutter kernel (§4 of the
   model doc, §7 / S2 here).

Nothing else about the population model changes. Precision ≫ recall, population-first
pooling, spanning-reads-only, single uniform ploidy, per-sample `F` — all inherited
from `ssr_genotyping.md` §1 unchanged.

---

## 2. Input contract — what Stage 1 Mark-2 writes

Per sample, one `.ssr.psp`. Per locus (identity = `(chrom, start, end)`, the
TRF-flagged region; matched across the cohort by coordinate, spec §2):

- **The observed distribution:** distinct repeat-region **sequences** + their
  integer **counts** (`Vec<(Box<[u8]>, u32)>`, sorted by sequence bytes for
  byte-identity). This is the per-sample empirical ladder — no per-read rows, **no
  base qualities** (the Stage-1 first-quartile gate already made survivors
  uniform-quality).
- **QC scalars** (the lean set): `depth`, `mapped_reads`, `n_filtered`,
  `n_low_quality`, `n_border_off_end`. `n_usable` = Σ observed counts (derived). These
  counts are **duplicate-free**: PCR/optical duplicates are marked **upstream** (post-mapping
  `markdup`, BAM flag `0x400`) and **Stage 1 skips flagged reads** when tallying — so
  `ssr-call` never sees duplicates and does **no** dedup of its own (verify-fix #5).
- **From the catalog** (Stage 0, reference-only): the locus `motif` and the
  **reference tract sequence + flanks** — the coordinate frame, not an allele claim.

(As-built: [src/ssr/pileup/locus_tally.rs](../../../src/ssr/pileup/locus_tally.rs),
schema in [src/psp/registry_ssr.rs](../../../src/psp/registry_ssr.rs).)

What Stage 2 does **not** receive (vs Mark-1): no `Qᵣ` / log-likelihood profiles, no
on/off-ladder columns, no per-read length histogram with stutter pre-classified —
just sequences and counts.

---

## 3. Reuse map — what carries over from `ssr_genotyping.md` §5 unchanged

These pieces are **candidate-representation-agnostic** — they only need "a candidate
allele set, a per-(sample,locus) per-allele read likelihood, and a whole-unit step
between alleles," all of which Mark-2 still supplies. They are **mostly inherited**;
this doc restates the Mark-2-specific *adapter*, not the math — the one redesigned
piece is the π M-step's base measure (the SNP class-pseudocounts → the SSR per-candidate
`G₀`, §5.4/§5.5 rows below).

| from `ssr_genotyping.md` | Mark-2 status |
|---|---|
| **§5.3 genotype prior** — inbreeding-adjusted IBD-mixture (`F·π_i + (1−F)·π_i^ploidy`), reuse `src/var_calling/posterior_engine.rs` | **reused** — the IBD-mixture prior transfers as-is (it is allele-agnostic); feed it repeat-allele *sequences* instead of bases (adapter), and `F` as the §4.4 outer-loop input |
| **§5.4 EM topology** — E-step responsibilities, M-step for `π` (Dirichlet-smoothed); confident-homozygote seed; identifiability (shared kernel + population recurrence); convergence (non-decreasing penalised log-lik) | **reused for the per-locus loop, with pieces replaced/added.** The **E-step** (prior × likelihood → posterior) and the **IBD prior** transfer as-is; the alignments behind `Qᵣ` are constant (`ε` frozen), re-weighted each iteration by `S_θ(Δ)`. **Replaced:** the π M-step's **base measure** — the engine picks pseudocounts by SNP/indel allele *class* (REF 10 / SNP-alt 0.01 / indel-alt 0.00125, [`classify_allele`](../../../src/var_calling/posterior_engine.rs)), SNP-shaped; SSR injects the per-candidate geometric `G₀` (§5.5 row). **Added:** a **per-locus `θ_locus` M-step** (regularized stutter update, shrunk to the sample's `(group, period)` shape prior — M3 — new code, not in the engine). `ε` is frozen (§4.4 pre-pass, not an M-step); `F` is the §4.4 prior-side outer loop. SNP-only machinery (compound alleles, chain anchors, contamination) is **bypassed**. The cross-locus schedule is **new code wrapping** the engine (§4.4) |
| **§5.2 stutter kernel `S_θ`** — 3-param geometric `(u,d,ρ)`, covariate-parameterised, pooling across (sample,locus) | **kernel form unchanged**; **re-parameterised in §4.4** as `shape(group, period) × level(group, length)` (shape per `(group, period)` shrunk to a cohort-per-period parent + per-locus refinement — M3; level per sample group, linear in length) — *not* one `θ(period,length,…)` loci group; the *whole-unit step `δ`* it operates on is redefined sequence-wise (§7 / S2) |
| **§5.5 Dirichlet base measure `G₀`** — unimodal, reference-centred in the signed unit offset `Δ`; geometric v1 / Gaussian upgrade | **adapted + injected.** `Δ` is still computable (ref tract length is the frame); now **purity-agnostic** (Mark-1 `OFFLADDER_PRIOR_FACTOR` **dropped** — §4.3/§7). The engine's class-based SNP pseudocounts are **replaced** by the per-candidate geometric `G₀(Δ)`, **centred on the per-locus cohort modal allele** (reference = the Δ frame only, *not* the prior's peak — §4.3): we **generalize the engine to take a per-candidate pseudocount vector** (the psp-container precedent — extract the reusable core, write the thin SSR part), rather than contort SSR onto the SNP class scheme |
| **§5.6 small-N / single-sample**, **§5.7 ploidy & `F`**, **§5.8 FP-aversion gates** (posterior / depth / support / stutter-fraction thresholds) | **§5.6/§5.7 unchanged; §5.8 gains an allele-balance term** (binomial-tail on the deconvolved per-allele responsibilities, the SSR analog of SNP `qual_refine.rs`) feeding GQ + site QUAL — the depth-driven-het-inflation defence (§6) |
| **§5.9 output VCF** — GangSTR/TRtools-compatible, REF/ALT = actual tract sequences, derived `REPCN`/`BPDIFFS` | **adapted (§4.5)** — the GangSTR/TRtools *shell* reuses (REF/ALT sequences, REPCN/BPDIFFS, GT, allele pruning, header detection); three semantics are **redefined for SSR**: emit iff **variable** (not "non-ref"), **site QUAL = Phred(locus variable)** (not "P(all hom-REF)"), and SSR **FILTER** reasons (`notPeriodic`/`tooManyAlleles`, drop `segdup`) |

The generative model (§5.1: `π` → genotype → read via allele-pick × stutter ×
sequencing-error) is **identical**; only the three leaf definitions below change.

---

## 4. The Mark-2 Stage-2 pipeline (overview)

```
 N × sample.ssr.psp ─┐                          per locus ℓ:
 (seq, count)        │   coordinate-merge   1. POOL    per-sample dists → cohort-aggregate seq dist
 .ssr.catalog ───────┼─►  scan (spec §2)    2. ASSEMBLE candidate set A_ℓ from pooled seqs   ← S1 (§5)
 (motif, ref frame)  │                       3. SCORE   Qᵣ(obs_seq | cand) flat-error, once  ← S3 (§6)
                     │                       4. EM      π+θ_locus over A_ℓ (ε frozen §4.4; F outer); reach. δ ← S2 (§7)
                     └─►                      5. CALL    genotype + posterior + VCF (§4.5)        ─► cohort.vcf
```

(A cheap **pre-pass** — §4.4 — runs once before step 4 to estimate and **freeze**
`ε` and `θ` from the cohort's confident **genotypes** (homozygotes ∪ well-separated hets —
CG-seed §4.3/§4.4); `F` is then re-estimated in a
prior-side outer loop *around* step 4. So step 4 itself only moves `π`.)

Step 1 and the EM core of step 4 build on §3's reuse; the reading/orchestration that
*feeds* step 1 is settled in §4.1, §4.2 frames the whole model the open sections detail,
§4.4 settles how the shared parameters are estimated (`ε` frozen; `θ` a per-loci-group prior
refined per locus; `F` a prior-side loop), and §4.5 settles the **output** (step 5: emit
variable loci, site QUAL = Phred(variable), SSR FILTER reasons). Steps 2, 3, and the `δ`
inside 4 are the new work this spec must settle (§5-§7, agenda §9).

### 4.1 Reading & orchestration — assembling the cohort evidence  *(settled 2026-06-19)*

*Execution-model **intent** (topology + contracts), not module/struct/signature
shape — that stays deferred to the architecture doc (header). This settles how
step 1 (POOL) of the pipeline is fed.*

**Same-catalog precondition.** Every input `.ssr.psp` must declare the **same SSR
catalog** (header md5); a mismatch is a **hard error**, checked once at open. The
catalog is the authoritative, ordered master list of loci — every sample's file is
a *subset* of it (the loci that sample had coverage at), in coordinate order.

**Per-sample reader = a coordinate cursor with a one-block decode cache.** Each
reader holds its sample's block index in memory and **at most one decompressed
block**. The merger calls `fetch(locus)` with loci in **ascending catalog order**
(monotonic — never backwards); the reader:
1. consults its block index — is there a block whose `[first_pos, last_pos]`
   covers this locus? If none → return **absent**;
2. if the covering block isn't the one in hand, **drop it and decompress the
   covering block** (the CPU-heavy step, done lazily — "a locus from a different
   block was asked for");
3. advance a **forward within-block cursor** to the locus: a record exactly there
   → return it; the cursor already past it (the sample spans the region but had no
   reads *at this locus*) → **absent**.

So **absent = "no data for this sample here"** in all three cases. Block
decompression is the reader's private responsibility; the merger never sees blocks
or zstd.

**The k-way merger.** Walks loci in catalog order, asks **every** reader for the
current locus, and gathers the present responses into one **`CohortLocus`**.
**Sparse-omit:** a locus where *every* reader is absent is **dropped, never
emitted**; a `CohortLocus` is built only when ≥1 sample has data.

**Lockstep ⇒ bounded memory.** Because all readers are asked for the same locus in
order, each pins only the one block covering the front → resident working set ≈
**N samples × one decompressed block** (the cohort-scaling property). Today's
byte/count-gridded blocks make the lockstep approximate (~1–2 blocks/sample);
moving the Stage-1 writer to **genomically-aligned blocks** (a future task,
mirroring the SNP caller's large win) makes it exact and lets all readers refill at
the *same* boundary.

**The analysis unit — `CohortLocus`.** The catalog **frame** (motif, reference
tract + flanks) + the per-sample observed `(seq, count)` distributions (absent =
missing) + the per-sample QC scalars. Tiny per locus (a handful of distinct
sequences each), which is why **decompression**, not assembly, is the cost — and
why the decompress work is the thing worth a worker pool.

**Execution topology (mirrors the SNP cohort pipeline).**
```
 per-sample readers ─► merger (producer) ─►  bounded queue  ─► EM worker pool ─► writer
 (lazy block decode)   k-way by locus,        (locus-batches,   (one EM per        (reorder batches
                       builds CohortLocus      back-pressure)    locus)             → catalog order → VCF)
```
- **One worker pool** runs the EM (the "single pool"); the merger is one producer
  thread, the writer another (`main producer + W workers + 1 writer`, the SNP
  shape). Dedicated queue-draining workers — *not* a pool shared with decode — is
  what sidesteps the SNP path's decode-starvation.
- **Work-item = a batch of *K* `CohortLocus`** (the SSR analog of
  `target-variants-per-chunk`), to amortize channel traffic; the bounded queue
  depth caps peak resident batches.
- **Output ordering is trivial:** loci are independent and catalog-ordered, so the
  writer reorders finished batches by catalog order (a `BTreeMap` keyed on batch
  sequence, as the SNP writer does).
- **Decompression placement:** on the **producer thread** to start — a block serves
  hundreds of loci, so refills are infrequent and overlap the EM through the queue;
  fan the per-sample refills across the pool only *if* profiling shows the producer
  is decode-bound (decide empirically — arch §7).

**What SSR deliberately drops from the SNP reader, and why.** The SNP unit of
analysis (a variant group) is *discovered from the data*; the SSR unit (a locus)
is *given by the catalog*. So SSR omits three whole layers:
- **no light/heavy two-phase decode** — there is no "is this row variable?" fold;
  every covered locus is genotyped, so decode *all* columns (the SSR schema is
  already flat);
- **no straddler / safe-gap / watermark** — a locus is an atomic interval that
  lives **wholly inside one block** (an **invariant Stage-1's writer must honour**:
  never split a locus across a block boundary);
- **the merge is a plain coordinate k-way merge**, not position streaming.

### 4.2 The statistical model — the EM and its ingredients  *(framing; details in §5–§7)*

**The whole approach in plain words (for the non-statistician).** We never observe a
sample's true genotype — only a tally of sequences, blurred by PCR stutter and
sequencing error. So we **guess and check**. The **hypothesis** is a set of
**candidate alleles** (§5) — the sequences we believe are real at this locus. The
**unknowns** we guess are how common each candidate is across the population (the
**allele frequencies `π`**) and how PCR slips (the **stutter `θ`**), both *seeded*
from the confident samples (§4.3). Then we iterate two steps:
- **score each guess against the data.** For every sample, weigh each possible
  genotype by combining (a) how plausible that genotype is *a priori* — from the
  frequencies and the inbreeding `F` (the **prior**, §5.3) — with (b) how well it
  reproduces the sample's observed sequences once stutter + error are applied (the
  **likelihood**, from the HMM alignment, §6). Normalized, the product is a
  probability over genotypes for that sample (its **posterior**).
- **update the guesses** (`π`, `θ`, `F`) to match what those posteriors imply.

Each round fits the data at least as well as the last; we stop when it stops
improving (**convergence**). This guess→score→update loop is the **EM** (§3 reuse),
and the **final per-sample posteriors *are* the genotype calls**. The rest of this
section names the pieces precisely; §4.3 and §5–§7 fill them in.

**How one read is born (the generative chain).** The prior and the likelihood are
two consecutive stretches of one story:
1. the population has **allele frequencies `π`** — candidate alleles (each a repeat
   sequence) and how common each is;
2. an individual **draws a genotype** (for a diploid, an unordered allele pair)
   from `π`; the **fixation index `F`** tilts the draw toward homozygosity (excess
   over random mating);
3. for each read, **one of the individual's alleles is picked** at random;
4. **PCR slippage (stutter `S_θ`)** may miscopy it by **whole repeat units** (±1
   motif copy most often), changing the length;
5. **sequencing + alignment** add **per-base** noise — a **miscalled base**
   (substitution) and small **flank** indels / a slightly misplaced flank. A sub-motif
   indel *inside the tract* is **not** a separate process here: at the repeat the only
   length change is the whole-unit slip of step 4, and a base change that breaks the
   motif is an interrupted (impure) allele, not in-tract noise (the in-tract no-gap rule,
   §6);
6. we record the sequence; over reads, the per-sample `(seq, count)` tally.

Steps 1–3 are the **prior** (`π`, `F`, ploidy); steps 4–5 are the **likelihood**
(`S_θ` ⊗ the flat per-base model `Qᵣ`). Step 6 is what Stage 1 stored.

**The EM spine — frozen `ε`, then a per-locus loop over `π` and `θ_locus` (schedule in §4.4).**
Genotypes (step 2) are the **hidden variable** — never observed, so a genotype is
*not* a parameter. The per-base error `ε` is **estimated once and frozen** before
genotyping (§4.4); the genotyping pass then runs an **independent per-locus EM** that
estimates the local allele frequencies `π` **and** the local stutter `θ_locus` (shrunk
toward the sample's `(group, period)` shape prior, §4.4 — M3). It builds on the `posterior_engine.rs` §5.4 per-locus
loop — its **E-step and IBD prior reused as-is**, its **base measure swapped** for the
SSR `G₀` (§3, §5.5), with a **new per-locus `θ` M-step** added:
- **E-step** — given current `π, θ_locus` (and the frozen `ε` and current `F`), compute
  each individual's posterior weight over genotypes (`posterior ∝ prior × likelihood`,
  normalized per sample): the responsibilities;
- **M-step** — given those weights, re-estimate the local **allele frequencies `π`**
  (Dirichlet-smoothed by `G₀`) and the local **stutter `θ_locus`** (shrunk toward the
  sample's `(group, period)` shape, M3); `θ_locus` only **re-weights** the (iteration-invariant) alignments — `ε`
  fixed, so the `align` values don't change when `θ_locus` moves (no recompute/rebuild).

`F` is the **only** parameter estimated *across* loci, in a **prior-side outer loop**
around the per-locus genotyping (§4.4): it sits in the genotype prior, not in `align`,
so re-estimating it is a cheap deterministic global reduce that never rebuilds the
alignment cache. The per-sample genotype calls we report are a read-off of the
**final** E-step.

**The ingredients, and their status:**

| ingredient | role | status |
|---|---|---|
| **candidate set `A_ℓ`** | the support of `π` and the genotype space | **assembled upstream of EM** (S1, §5, OPEN) — gates everything: too many ⇒ stutter blobs become "alleles"; too few ⇒ missed variation |
| **allele frequencies `π`** | population prior over alleles | estimated (M-step) |
| **fixation index `F`** | excess-homozygosity knob; IBD-mixture prior `F·π_i + (1−F)·π_i^ploidy` | **per-individual `F_i`, estimated by the §4.4 prior-side outer loop** (EM-responsibility reduce over variable loci, shrunk toward the cohort mean; in the prior, **not** in the `align` cache → never rebuilds it); supplied default as the burn-in start. Hard ceiling `F_CEILING = 0.99` (no `F=1` absorbing trap) + optional lower CLI cap. Output is **apparent `F_IS`** — absorbs cohort-wide hom-excess, chiefly **Wahlund/population structure** (real, *documented warning*, not corrected; M5); null alleles largely N/A for primer-free sequencing (§9) |
| **ploidy** | sets the genotype-space size | fixed, uniform |
| **stutter `S_θ`** | whole-unit slippage kernel | **shape prior from the §4.4 pre-pass, refined per-locus in the EM**, × a **per-sample-group level** — `S_θ(Δ) = level × shape(Δ)`. **Shape:** **per `(group, period)`, shrunk toward a cohort-per-period parent `θ_period`** (M3 amend.; length is *not* a shape covariate — it drives the level), off the confident **genotypes** (homozygotes ∪ well-separated hets — CG-seed amend. §4.3/§4.4; data-rich groups that differ get their own profile; thin `(group, period)` cells borrow the parent; literature default for empty); each locus refines a `θ_locus` shrunk toward the sample's **group**-period shape (depth-driven partial pooling — Mark-1's relief valve; per-motif outliers — AT-di / homopolymers — ride the per-locus refinement). **Level:** **per sample group, linear in length** (`level_baseline + level_slope·length`, both coefficients per group), §4.4 amend. Enters `align` only via `S_θ(Δ)`, so it is **refined in the prior-side outer loop** (alongside `F`, C2 amend. 2026-06-22) — a cheap **re-weight** (no rebuild), pre-pass-seeded, **not frozen**. Motif & purity covariates deferred. The one place loci aren't independent, and what makes `θ` learnable |
| **per-base error `ε`** | flat (no per-base quals) noise inside `align` (§6) | **estimated once and frozen by the §4.4 pre-pass** — a **per-sample-group** value (data-driven soft clusters, §4.4 amend.) — not one global scalar — off the confident **genotypes** (homozygotes ∪ well-separated hets — CG-seed amend.), identifiable alongside `θ` (`θ` changes **length** in whole units, `ε` changes **composition** via substitutions — `align` admits no sub-motif indel inside the tract, §6, so this holds at period 1 too; mononucleotides flagged); seeded from the Stage-1 quality gate as the burn-in start. **Sample-group** (not per-sample) granularity *because* `ε` lives in `align`: frozen ⇒ `align` is a pure, iteration-invariant function → genotyping is **deterministic across/within thread counts** *whether recomputed or cached* (the cache is **deferred**, §6; per-sample `ε` = upgrade). A wide *within-group* spread flags the grouping is too coarse (period covariate on `ε` the documented upgrade) |
| **base measure `G₀`** | regularizing prior **on `π` itself** | reused §5.5; **cohort-mode-centred** (reference = frame only, §4.3), smooth — a **FP control** + small-N regularizer |
| **seed + convergence** | EM start + stop | reused §5.4 — confident-homozygote seed; iterate to non-decreasing penalised log-lik |

**Why the deconvolution is well-posed (identifiability).** With `π, F, θ` free over
a candidate set, a stutter peak could be explained *either* as slippage of a common
allele *or* as a real rare allele. The tie is broken by the **population**: a true
allele recurs across unrelated individuals, whereas a stutter shadow's abundance
tracks its parent's at the *learned* slippage rate. So cohort-wide pooling (and the
shared `θ`) is what makes the problem identifiable, not merely more powerful — hence
the small-N caveat (reused §5.6).

*(HipSTR's own strongest het-resolution lever — **phasing the STR against nearby
heterozygous SNPs** to assign each read to one of the two alleles — is **deliberately
excluded** (Mark-1 §1.1): SSRs are too sparse for reads to co-span an STR and an
informative SNP, and it would couple the two callers. We resolve the same merged-het /
het-vs-stutter ambiguity instead through **cohort recurrence + the per-locus `θ`** above —
a population lever rather than a physical-phasing one.)*

**Three distinct noise sources — keep them separate.** (1) **stutter** —
whole-repeat-unit length change (PCR, dominant, length-dependent); (2) **per-base
error** — substitutions / small indels / flank misplacement (the flat `Qᵣ`); (3)
**mis-assignment** — a read from a paralog/duplication or cross-sample
contamination. Stage-0/Stage-1 gates remove most of it; the residual is handled in the
genotype likelihood (§6) by **two complementary terms, because the residual has two
shapes**:
- a **uniform outlier `λ`** absorbs *randomly scattered* junk — without it the EM
  explains every oddball read as evidence for some allele, the cheapest explanation
  being an extra allele → false hets;
- an **allele-balance / overdispersion term** stops a *systematic, concentrated* minor
  signal (a paralog/contaminant/mismap at a consistent length, or under-modelled
  stutter) from inflating a false heterozygote whose confidence *grows with depth* —
  the i.i.d.-likelihood failure mode already seen on the SNP path. It mirrors the SNP
  caller's allele-balance QUAL refinement and feeds both GQ and site QUAL (§6).

Cross-sample **contamination** (the heavier dedicated mechanism) is deferred (reuse
the existing `estimate-contamination` machinery later).

### 4.3 Seeding the EM — π⁰, θ⁰, and the prior on π  *(settled 2026-06-19; `G₀` shape = geometric for v1, Gaussian upgrade open — §9)*

EM lands in whatever basin its seed points at, and this is a deconvolution (allele
vs stutter), so the seed matters. Two things are seeded — the starting frequencies
**π⁰** and the starting stutter kernel **θ⁰** — and a **prior on π** (the base
measure `G₀`, §5.5 reuse) both floors π⁰ and regularizes every later M-step.

**Putative genotype per sample → the π⁰ tally.** Give *every* sample a putative
genotype off its rungs (the level-4 clear-maximum rule, §5), then tally the alleles
(one genotype contributes `ploidy` allele-copies, regardless of depth):
- **≥ ploidy clear maxima** → genotype = its top-`ploidy` clear peaks;
- **one clear maximum** → **homozygous** for the **most-abundant sequence** on that
  peak (the allele; the rest of the peak is its per-base-error halo);
- **zero clear maxima** (too thin to peak) → no contribution.

This deliberately **mislabels merged hets** (two adjacent alleles fused into one
peak) as homozygotes — accepted, because excluding the homozygotes (the cleanest,
commonest samples) would be far worse, and the error is **self-correcting**: once
θ⁰ is seeded (below), the first E-step sees a "homozygote" whose peak is too tall or
whose −1 shoulder is too fat to be pure stutter and moves posterior weight onto the
adjacent het — **provided that masked allele was not seeded at exactly 0**, which the
pseudocount prior (below) guarantees. The seed only has to point the EM at the right
basin.

**θ⁰ — read off the confident homozygotes' skirts.** The **stutter kernel** is the
slip-size table `S_θ(Δ)` ("if the allele has N repeats, the chance PCR returned
N+Δ"; form §5.2, used in §6), shaped by `θ = (u, d, ρ)` — the **up-slip rate**, the
**down-slip rate** (`d > u` for SSRs: slippage shortens), and the geometric **decay**
of multi-unit slips. A confidently-homozygous sample seeds it *for free*: its
genotype is known `(A, A)`, so **every read that is not exactly `A` is a stutter
product of `A`**, and the skirt's shape *is* the table — the faithful fraction gives
`1 − u − d`, the down-side (`−1, −2, …`) vs up-side (`+1, +2, …`) split gives `d` vs
`u`, and the `−2`-to-`−1` ratio gives `ρ`. One skirt is noisy, so **pool every
confident homozygote of the same period** (whole cohort
+ genome) before fitting the shape `(u, d, ρ)`. The **same pass** that builds the π⁰ tally
produces these skirts, so π⁰ and θ⁰ come out of one sweep. *Fallback:* a period with
no confident homozygotes starts its shape from a generic literature stutter profile; a
sparse or high-variance period is **shrunk** toward a pooled-period
parent (§4.4). The pre-pass `θ` is only a **per-loci-group prior** — each locus's EM refines a
`θ_locus` shrunk toward it (§4.4) — so seed contamination (a masquerading het below) is
**low-stakes**: the EM can overrule a biased prior. The confident-homozygote selection
therefore only needs to be *reasonably* clean (depth gate; flag a one-sided skirt excess
— the masquerader signature); it is prior-hygiene, not the final word (Issue 5 / §9).
*(Refinement, settled 2026-06-21: make the gate **model-based** rather than heuristic.
The right test is **one-allele-vs-two-allele**, not a generic goodness-of-fit — the
question is precisely "one allele or two," since a hidden het is what poisons the
estimates. Score the reads under the best one-allele model (called length + predicted
skirt, from the burn-in's current params) vs. the best two-allele model (reusing the EM
scoring + candidate machinery); the second allele must **earn its place** via a
complexity penalty (BIC-style, no hand-picked p-value). Tune the penalty on the simulator
**for purity** — let through almost no hets even at the cost of discarding real
homozygotes (a sneaked-in het corrupts the estimates; a discarded one just costs a little
data). Plus a minimum-depth skip. Arch doc Q-P7.)*
*(Settled 2026-06-21; shape amended 2026-06-22 (M3): the kernel is
`shape(group, period) × level(group, length)`; the skirt-pooling here fits a
**per-`(group, period)` shape shrunk to a cohort-per-period parent**, while each sample
group's stutter **level** (linear in length) and `ε` are fit **per sample group** — the
per-sample estimates are only the intermediate that clusters samples into groups — §4.4 amendment. The
confident-homozygote selection **co-evolves with the parameters in the burn-in loop** —
the gate's model sharpens as the parameters settle, §4.4.)*

> **Amendment (CG-seed, 2026-06-22) — the skirt source is the confident *genotype*, not only
> the homozygote.** "Every read that is not exactly `A` is a stutter product of `A`" holds for
> a homozygote, but a **confident, well-separated heterozygote** (`A,B`, alleles **≥ 2 repeat
> units apart**, two clear peaks with a real valley) supplies the *same* labelled skirts: reads
> **below `A`** are `A`'s down-stutter, reads **above `B`** are `B`'s up-stutter, and the inner
> valley is split by the soft EM (C2). So `θ⁰` (and `ε`/level) pool over **confident genotypes
> = homozygotes ∪ separated hets**, which makes the estimate **data-rich in hyper-heterozygous
> cohorts** where homozygotes are scarce (closes m2(a)) and adds **two length anchors per het**
> for the level-vs-length slope. The one-allele-vs-two-allele gate above generalizes to a
> **1..ploidy-peak resolution** test (guards: peaks ≥ 2 apart, dosage-consistent heights, each
> allele cohort-recurrent — rejects a hom+heavy-stutter masquerading as a het); a merged
> (1-apart) het simply **fails resolution and is not seeded** (no longer mislabelled as a hom).
> Polyploid resolution is harder → lean on coded priors when peaks won't resolve. See parameters
> §2 for the full statement.

**The prior on π — pseudocounts (`G₀`).** A candidate's frequency `π_i` enters the
genotype prior **multiplicatively** (§5.3), so `π_i = exactly 0` is an **absorbing
trap**: every genotype using allele `i` then has prior 0, its posterior is 0 in every
E-step *however strongly the reads support it*, it earns no responsibility, and the
M-step leaves it at 0 — forever. And a candidate *can* land at exactly 0 in the seed,
because the **candidate set (recall-generous, S1) is broader than the seed tally
(confident genotypes only)**. The chief case is the **masked het**: an allele added
by the ±1 rescue — so it *is* a candidate — that in every carrier sits under a taller
peak and is seeded as part of a homozygote, so it is **never counted** (other minor
sources: a candidate seen only in too-thin samples, a peak beyond ploidy, a reference
allele no confident sample carried). The fix is **pseudocounts**: pretend a small
fractional count of each candidate was already seen, so every candidate starts
non-zero and the EM can grow it if the data warrant, while the pseudocount vanishes
under real evidence. This is also exactly what makes the seed's deliberate merged-het
mislabelling (above) **recoverable** — the masked allele sits at the small floor
instead of 0, so the self-correcting E-step *can* put weight on its het genotype; at
a hard 0 that correction would silently fail. So the putative-genotype seed and the
pseudocount prior are a **pair**: the seed is *allowed* to be optimistic only because
the prior keeps the alleles it skips alive. Made **unequal**, the pseudocounts
further encode prior plausibility. The shape is a **geometric decay** in the unit
offset (the geometric's heavier-than-Gaussian tail keeps real large-step alleles
callable), **centred on the per-locus cohort modal allele** — *not* on the reference —
and **symmetric** in v1, with its **decay parameter fit from pooled empirical data**
per loci group — a well-defined *parametric* prior, deliberately **not** a
data-hungry non-parametric one (pinning a shape needs many observations, so we fix the
shape as geometric and estimate its parameter, exactly as for the stutter kernel).
**The reference is only the coordinate frame** (it defines the `Δ` unit grid), not the
prior's peak: centring on the reference would re-import the very reference bias Mark-2
removed (an off-reference population's true modal allele would get a `p^|Δ|` floor so
small the EM could never grow it — defeating the masked-het recovery above), whereas
centring on the empirical mode is SMM-consistent (a population's allele-size
distribution is unimodal around *its own* mode) and degrades gracefully — with no data
the only candidate is the seeded reference, so the mode *is* the reference. It
is purely a **length-offset** prior — **agnostic to purity**: an impure allele of
length L gets the same prior as a pure allele of length L (impure alleles are *not*
specially penalized — §7). This is the §5.5 Dirichlet base measure; it is a
**false-positive control** (alongside the HMM likelihood and population recurrence),
and the **same** pseudocounts re-enter every M-step, so it regularizes throughout, not
only at the start.

**The pseudocount needs an explicit numeric floor (verify-fix #4, 2026-06-22).**
Mode-centring kills the *systematic* reference-bias floor (above), but `p^|Δ|` still
underflows to **exactly `0.0` in f64** once `|Δ|·(−ln p) ≳ 744` — and a candidate exists
only because a read supports it, so a far candidate (large `|Δ|`) is exactly the masked
allele the pseudocount must keep alive. Candidates are **read-length-capped**, so short
reads stay safe (`|Δ| ≤ read_len/period`: e.g. 150 bp / di-repeat ⇒ `|Δ| ≤ 75`, and
`75·(−ln 10⁻³) ≈ 518 ≪ 744`), but a **long read** spanning a large tract with a **steep
decay** can cross the bound (`|Δ| ~ 500` ⇒ `e⁻³⁴⁵⁰ = 0`) and re-create the `π = 0`
absorbing trap *numerically*. So the pseudocount is **floored explicitly**:
`G₀[i] = max(p^|Δ|, FLOOR)` for a tiny positive `FLOOR` (any representable `> 0` re-floors
it above hard zero — the Dirichlet M-step grows `π_i` from real evidence regardless of how
small the pseudocount is; the floor only has to keep a *zero-evidence* candidate alive for
future iterations), or — cleaner — `G₀` is carried in the **log domain** into the M-step
(the engine sums `expected_counts[a] + pseudocount[a]`, so a log-domain pseudocount is the
principled form, [`posterior_engine.rs`](../../../src/var_calling/posterior_engine.rs)
`m_step_p_hat`). This is the `π` analogue of the `F_CEILING = 0.99` clamp (§4.4) — both
keep a multiplicatively-absorbing parameter off its trap value; the `F` ceiling is already
a hard clamp (airtight), this gives `G₀` the matching hard floor.

**Combine:** `π⁰ = (putative-genotype counts + G₀ pseudocounts), normalized.` With
**no confident samples** at a locus (all under-resolved / very thin) it falls back to
the **prior alone** — π⁰ is the normalized `G₀` pseudocounts (centred on the modal
candidate, which with no data *is* the seeded reference allele) — which is also the
correct conservative behaviour at small N (§5.6). **`F⁰`** is a separate seed (a supplied / global default that is also
the burn-in start of the §4.4 prior-side `F` loop that estimates `F`). *(v1 shape =
**geometric**, settled (heavier tail keeps real large jumps callable); **open — §9:** whether a
**Gaussian** form is ever worth it — a documented upgrade, not a v1 gap, `ssr_genotyping.md` §5.5.)*

### 4.4 Estimating the shared parameters — pre-pass (freeze `ε`, prior for `θ`) + a prior-side `F` loop  *(settled 2026-06-19)*

The chemistry parameters are **characteristic of the library *protocol*** (chemistry +
DNA preservation), not of any one locus — and, when samples come from different
protocols, **not of the whole dataset either** (the 2026-06-20 amendment below corrects
this; the unamended text treats them as cohort-global). They are estimated from many
loci in a cheap **estimation pre-pass** before genotyping. But they are **not all
treated the same way** — they differ in how locus-specific they really are and in what
they cost to vary:

- **`ε` (per-base error) is frozen.** It is a global property of the chemistry — no
  biological reason it varies locus to locus — so a single value pooled over the whole
  pre-pass (**per sample group** — amendment below) beats anything one locus could
  estimate from its handful of mismatches; and
  because `ε` lives *inside* `align` (changing it **rebuilds** the alignment cache),
  freezing it is also what keeps the alignments computable once. This is the standard
  **two-stage / plug-in** pattern (and what **HipSTR** does with its stutter model).
- **`θ` (stutter) is a prior, not a frozen value.** Stutter genuinely varies (by period,
  motif, purity — and, M3, by sample group), so the pre-pass produces a **per-`(group,
  period)` *shape* prior, shrunk to a cohort-per-period parent `θ_period`** (a
  per-sample-group *level* rides on top — amendment below) and each locus's EM **refines a
  `θ_locus` shrunk toward the sample's group-period shape** (depth-driven partial
  pooling: a deep locus moves off the group-period value, a shallow one stays at it). This is
  cheap — `θ` enters only through `S_θ(Δ)`, which **re-weights** the iteration-invariant
  alignments rather than triggering a recompute (`ε` fixed) — and it is what restores Mark-1's per-locus *relief
  valve* (impure / high-depth loci that disagree with their loci-group prior), the full
  allele-vs-stutter identifiability (cohort recurrence at the locus breaks the tie),
  and, critically, makes the pre-pass `θ` **just a prior the EM can overrule** — so seed
  contamination (a masquerading het, §4.3) is **low-stakes**, not permanent (Issue 5).

`F` (and, since the C2 amendment, the per-group stutter **level**) is refined in a
prior-side outer loop rather than frozen — both sit outside `align`, so refining them
rebuilds no cache (the level is seeded by the pre-pass, `F` by a default). `F` is
handled as a cheap outer loop (below).

**Why this split is still deterministic.** With `ε` frozen, `align(obs | cand ⊕ Δ)` is a
**pure function of `(obs, cand⊕Δ, ε)`** — so it is the same on every iteration and every
worker, *whether cached or recomputed* (v1 recomputes; the cache is a deferred optimization,
§6). Determinism does **not** depend on a cache: everything that varies is **per-locus, not
per-thread** — `θ_locus` is a pure function of that locus's data + the `(group, period)`
shape prior, computed identically on whichever worker owns the locus. (This is the crucial difference
from the per-thread `ε` we rejected — *per-locus* is atomic and reproducible; *per-thread*
was order-dependent.) So Stage 2 is byte-identical across and within thread counts, carrying
the Stage-1 determinism invariant forward.

**Amendment (2026-06-20; refined & settled 2026-06-21) — chemistry varies by sample
group, inferred from the data.** *(Vocabulary home: the glossary in
[ssr_call_parameters.md §0](../architecture/ssr_call_parameters.md) — `cohort` /
`sample` / `sample group` / `chemistry` / `stutter shape` `θ_period` / `stutter level` =
`level_baseline + level_slope·length` / **loci group** = the set of loci a pooled
parameter is fit over (loci-side analogue of a sample group), **= `period` in v1** /
**geometric pseudocounts** `G₀`.)* The
premise above ("characteristic of the whole dataset") is too strong. **PCR stutter *level* and per-base error are properties of the
library protocol and DNA preservation, not of the sequencing depth:** a PCR-free library off high-MW DNA,
a high-cycle library off scant fragmented historic DNA, and everything between have
**very different** stutter (and, via deamination damage, different `ε`). Pooling them
into one cohort-global `ε`/`θ` fits no sample and **mis-calls hets in proportion to the
mismatch** — apply low-stutter (PCR-free) parameters to a high-stutter (PCR) sample and
the unexplained stutter skirt becomes a phantom second allele → **het overcall**; apply
high-stutter parameters to a PCR-free sample and a genuine minor allele is absorbed as
stutter → **het undercall**. (It also protects `F_i`: mis-modelled per-sample stutter
otherwise leaks into apparent inbreeding, §4.4 `F` caveat.) So `ε` and the stutter
**level** carry a **group axis** — the two-treatment split above (`ε` frozen, stutter a
prior) holds **per sample group**, not across the whole dataset. The stutter **shape**
**also carries a group axis** (M3 amend. 2026-06-22, below) — it is no longer assumed
cohort-shared; instead a per-`(group, period)` shape is **shrunk toward** a cohort-per-period
parent, so the data decide whether groups' profiles differ.

*The decomposition (this is what keeps it affordable).* The field reports that
**protocol/preservation move the stutter *rate* far more than its *profile*** — PCR-free
vs PCR on the same samples drops the rate ~4.6× while the distribution of error sizes
stays *about* the same (Gymrek lab 2014); the rate is **linear in repeat length** (Willems
et al. NAR 2019; forensic LUS-linear) and shrinks with **period** (mono ≫ di ≫ tetra;
Chakraborty PNAS 1997). That evidence is human PCR chemistry, so for non-model / degraded
DNA / unusual motifs we treat "profile is invariant" as a **prior expectation realized by
shrinkage, not a hard assumption** (M3). So factor stutter into a (hierarchical, group-aware)
shape and a per-group level:

```
θ(group, period, length, locus)
   = shape(group, period)         — per (GROUP, period), SHRUNK toward a cohort-per-period
                                     parent θ_period (data-rich groups that differ get their
                                     own profile; thin groups borrow the shared one — M3)
   × level(group, length)         — per SAMPLE GROUP, linear in length:
                                     level_baseline + level_slope·length (both coeffs per group)
   × θ_locus refinement(locus)    — per-locus EM toward the sample's GROUP-period shape
                                     (frozen ε; level refined per-group in the genotyping
                                     OUTER LOOP, alongside F — C2 amend. below)
ε(group)                           — per sample group (in align → frozen; per-sample = upgrade)
```

The **granularity default is the sample group for every chemistry parameter, shape
included** — a parameter is pooled up to the cohort *only* when the data show it doesn't
vary by group. Shape is **no longer the exception**: it gets the same per-group-with-
shrinkage treatment as `level` (M3 removes its former cohort-shared exemption). If the
groups' shapes genuinely agree (the PCR-free expectation), the shrinkage collapses them to
the shared per-period parent at no cost; if they diverge, the model captures it. A sample
group carries **two numbers** for its level (`level_baseline`, `level_slope`) plus its
per-period shape `(u,d,ρ)` deviation from the parent. **`ε` is frozen** at sample-group
granularity *because it lives inside `align`* — but the cache itself is **deferred** (M2:
v1 recomputes; `align` is a pure function of frozen `ε`), so per-sample `ε` is a documented
modeling upgrade, not a cache question. **The stutter
level is *not* frozen** (C2 amend. below): it enters only via the `S_θ` re-weight, *not*
the `align` cache, so it is re-estimated **per sample group in the genotyping outer loop**
(alongside `F`) at no cache cost. The per-sample (`ε`, level) estimates *are* computed in
the pre-pass, but only as the **intermediate that clusters samples into groups** (below)
plus the level's outer-loop seed (`level⁰`); the grouping is fixed once, while the
per-group level *value* adapts. The per-locus EM refines only the locus *shape* `θ_locus`
(now shrunk toward the sample's **`(group, period)`** shape, M3), so **the chemistry axis
never sees per-locus data starvation** — the group-period shape pools across loci, the
locus refinement adds the motif/purity deviation on top. Consequence: the `θ_locus`
refinement pools **level-adjusted** reads (divide out the group's level before estimating
the locus shape), so within each outer sweep the per-group `level` is **fixed** (and `ε`
fixed for the whole run); the level is revised only **at the outer-loop barrier**, exactly
as `F` is — sequencing (level before shape) holds within a sweep, adaptation happens between
sweeps.

**Amendment (2026-06-22, M3) — the stutter *shape* carries a sample-group axis.** The
earlier design held the shape **cohort-shared, keyed by period**, on the strength of human
PCR-free evidence ("protocol moves the rate, not the profile"). That evidence does not
obviously transfer to non-model organisms, degraded/ancient DNA, or unusual motifs — and
the old architecture could not express a per-group shape difference even if one existed
(the only refinement axis, `θ_locus`, is per-*locus* and pools all groups at a locus). Fix:
give the shape the same per-group-with-shrinkage treatment `level` already has —
- **estimate a per-`(group, period)` shape** (full `(u,d,ρ)`, off the confident-**genotype**
  skirts — homs ∪ resolved-het outer skirts, CG-seed §4.3 — keyed by `(group, period)`; groups are already fixed by the pre-pass clustering,
  so this is just a finer accumulator key) **shrunk toward a cohort-per-period parent
  `θ_period`.** Data-rich groups that genuinely differ get their own profile; thin groups
  borrow the shared one. Invariance becomes a **data-driven outcome** (agreement ⇒
  shrinkage collapses to the parent at no cost; divergence ⇒ captured), and the
  per-`(group, period)` spread vs the parent is a ready-made **diagnostic** of whether the
  PCR-free invariance actually holds on *this* cohort.
- **`θ_locus` refines toward the sample's `(group, period)` shape**, not the cohort one —
  the per-locus piece (motif/purity, a HipSTR-style per-locus refinement) composes on top
  of the per-group profile.

This is purely a **finer key + an extra shrinkage level** on a quantity that already enters
only through the `S_θ` re-weight (never `align`/the deferred cache, never the per-iteration
HMM), so it is **cheap, deterministic** (the skirt accumulator stays an integer order-free
reduce, now keyed by `(group, period)`), and it removes the one chemistry parameter that
was exempt from the "sample group unless the data justify pooling" default. What we keep
from HipSTR is the per-locus free-shape refinement (`θ_locus`); the **group axis is the
part beyond HipSTR** (which pools all samples per locus and so cannot separate groups).

*Sample groups — data-driven, distance-based, precision-weighted (settled 2026-06-21;
arch Q-P6).* The user does **not** label protocols (SRA provenance is unreliable, and a
*wrong* label is worse than none), so groups are **inferred from the per-sample (`ε`,
stutter-level) estimates** — and the **main reason to group is data sufficiency**: a
low-coverage sample can't support reliable per-sample chemistry, so it **borrows strength**
from similar samples. (Grouping is a pure **estimation + reporting** device with **no cache
role** — the `align` cache is **deferred**, §6 / genotyping Q-G3, so there are no cache
"hard lines" to draw either way; arch Q-P6.) Mechanism:
- **Distance-based, not k-means.** Group **close neighbours** in (`ε`, level) space:
  similar parameters ⇒ similar EM behaviour ⇒ safe to pool. The number of groups **falls
  out of a closeness threshold** (the data decides; a single-protocol cohort collapses to
  one), not a preset K. Deterministic — no random initialisation, unlike k-means — **and
  deterministic under ties (verify-fix #1): strict `<` thresholds with a defined
  equal-handling rule and a total tie-break on the sample's catalog index, so a borderline
  `distance ≈ threshold` (or two equidistant merge candidates) resolves identically across
  thread counts.** The per-sample `(ε, level)` it groups on are themselves order-independent
  fixed-point reduces of each sample's labelled observations.
- **Distances scaled by each sample's uncertainty.** A thin sample's (`ε`, level) is
  wobbly, so measure distance in units of its error bar: an uncertain sample is "close to
  many things," never confidently flung into a noise-chosen group. This is how grouping
  protects the very samples it exists to help.
- **Threshold by penalized-likelihood comparison + enough-data (M4 amend.).** Merge two
  candidate groups only when the **split does not earn its extra parameters**:
  `ℓ_pen(split) − ℓ_pen(merged)` below the BIC penalty for the additional per-group
  chemistry parameters (and each resulting group has enough data to estimate well). This
  replaces the old "merge when it doesn't move calls" — it asks whether the data *support*
  two regimes, not whether the *output* happens to flip.
- **Per-sample shrunk toward its neighbours** — partial pooling: a data-rich sample uses
  its own value; a thin one falls back to its group, *not* to a cohort average dominated by
  the clean majority. (Same machinery as `θ_locus → θ_period` and `F_i →` cohort mean.)
- **A group of size one is fine — the hierarchy handles it (m2 edge case).** A **singleton
  group** (a sample whose chemistry the BIC split test, M4, found genuinely distinct) just
  uses its **own shrunk estimate**; if it is too thin to support even that, the same
  partial-pooling **shrinks it toward its nearest neighbour / the cohort**, so it is never
  left estimating from nothing. The BIC test won't *carve out* a singleton regime unless the
  data support it, so an accidental singleton folds back in — no special-casing needed.
- **Coverage is an *uncertainty*, not a coordinate.** Depth (duplicate-free upstream, verify-fix #5) correlates with the
  chemistry parameters two ways: as **estimation variance** (∝ 1/depth, *unbiased*) and as
  a **confound through preservation** (low-input degraded DNA is *simultaneously*
  high-stutter, high-`ε`, and low-coverage). The second is **real protocol signal we must
  keep**, so coverage is **neither a clustering feature nor regressed out** — it enters
  **only as the precision of each sample's measurement** (distances are scaled by each
  sample's error bar ∝ 1/√(depth)). **Duplicate removal is upstream, not our concern
  (verify-fix #5).** PCR/optical duplicates are marked **after mapping** by the standard
  tools (`samtools markdup` / Picard set the BAM duplicate flag `0x400`); **Stage 1 simply
  skips flagged reads** when it builds the `(seq, count)` tallies. So the depth this caller
  sees is **already duplicate-free by construction** — there is no "assumed deduplicated"
  unchecked precondition and no within-caller dedup to do. Any residual PCR jackpotting is
  whatever the upstream marker missed (imperfect without UMIs) — an **input-quality matter,
  out of scope here**, not modelled by `ssr-call`. `ε` and stutter co-vary through the same regime, so
  they are clustered **jointly** into one group label per sample — a group is one
  **chemistry + preservation regime**, not two independent knobs.

*Selection and parameters co-evolve (the burn-in loop).* The pre-pass runs off confident
homozygotes, but the homozygote-confidence gate *itself* assumes a stutter level — so a
*fixed* gate would let a high-stutter sample's homozygotes look messy, fewer pass, and its
estimate bias **cleaner** than the truth (a degraded sample could then mis-assign toward a
low-stutter group). That is **selection bias, not just variance** — precision-weighting
cannot fix it. So the gate is **not fixed**: it co-evolves with the parameters in the
**burn-in loop** (below) — start from dev-computed defaults, admit homozygotes with the
current parameters, recompute, feed back, repeat until settle. "Re-select under each
sample's own fitted stutter" is just what each cycle does; the loop is deterministic
(seeded batches, frozen-params-per-batch map, barrier update — below).

**Amendment (2026-06-22) — de-bias the estimator; un-freeze the stutter level (C2
resolution).** The co-evolving *hard gate* above still has a failure mode the burn-in
cannot fix on its own: a **masquerading het** (two adjacent alleles merged into one peak)
that passes the gate is admitted as a homozygote, and its minor allele — a *different
length* — is counted as **stutter**, inflating the level; and because the gate's expected
skirt uses the *current* level, an over-estimated level admits *more* such hets, a
**positive feedback** loop. Two changes close it:

- **The hard gate becomes a *seed*, not the estimator.** The confident-**genotype**
  resolution test (§4.3, Q-P7 — homozygotes ∪ well-separated hets, CG-seed below) still **initializes** the burn-in with clean labels,
  but the **final** `ε` and level come from a **soft full-cohort EM responsibility reduce**
  over the pre-pass subset: every read contributes *fractionally* to allele-vs-slip via the
  posterior, so **population recurrence pulls a masquerading het's minor reads onto its
  *het* genotype** (the same identifiability lever genotyping uses) instead of into the
  skirt. Selection bias is gone (soft weighting never *excludes* contradicting data) and
  the feedback loop dies (no co-evolving hard threshold). Mark-1's dropped dominance /
  bimodality flag is kept as cheap **seed** insurance, not the final word.
- **The stutter level is un-frozen** (see the granularity block above). After the §6
  in-tract no-gap rule (C1), `ε` is *substitutions only*, so a masquerading het contaminates
  the **length/level**, **not** `ε` — and the level is exactly the parameter that is *cheap
  to move* (an `S_θ` re-weight, **not** in the `align` cache). So the level joins the
  prior-side **outer loop** alongside `F`: re-estimated **per sample group** from the
  cohort's soft per-allele responsibilities each outer round, accumulated by an
  **order-independent fixed-point integer reduce** (so thread count / locus-completion
  order can't change the bits — the same determinism mechanism `F` uses), rebuilding
  **no** cache. `ε` stays **frozen** (it is in `align`: rebuild cost + determinism, and no
  longer the contamination victim). Net: the contaminated, load-bearing quantity becomes
  **overrulable by the whole cohort** at ~zero cost, while the align-once + cross-thread
  determinism guarantees the freeze protected are untouched.

> **Amendment (CG-seed, 2026-06-22) — broaden the seed to confident *genotypes*, closing
> m2(a).** The seed gate is no longer homozygote-only: it admits any **confidently-resolvable
> genotype** — a homozygote (one clear peak) **or** a well-separated het (two peaks **≥ 2
> repeat units apart**, dosage-consistent heights, each allele **cohort-recurrent**), and for
> ploidy *p* up to *p* clear peaks. A well-separated het's outer skirts (below the low allele,
> above the high) are **labelled stutter** just like a homozygote's; its inner valley goes to
> the soft EM (C2). **Why this matters:** (i) it **closes m2(a)** — a hyper-heterozygous
> outbred cohort with few homozygotes *anywhere* is precisely where separated hets are
> abundant, so the chemistry seed is data-rich instead of starved; (ii) it **further kills the
> masquerade** C2 patched — a separated het is now *resolved as a het* (minor allele = allele,
> not stutter), so it is no longer mislabelled; a merged (1-apart) het simply **fails
> resolution and is not seeded**; (iii) a het gives **two length anchors** for the
> level-vs-length slope. **Guards** against a hom+heavy-stutter masquerading as a het: dosage-
> consistent peak heights + each allele cohort-recurrent (a stutter peak does not recur as an
> independent allele). **Polyploids:** the principle holds (1..*p* peaks) but dosage resolution
> is harder, so v1 seeds only the cleanly-resolvable cases and otherwise **leans on the app's
> coded priors**. **Residual m2(a)** (no confident genotype of *either* kind — all 1-apart hets
> / too thin): fall back to coded priors + soft-EM-via-recurrence, and emit a **loud "chemistry
> not estimated from data — running on literature defaults" warning** (an output obligation,
> like apparent-`F_IS`), with the confident-genotype count surfaced. Full statement: parameters §2.

*Reporting — the audit trail that replaces the missing labels.* Because protocols are
inferred, the pre-pass **reports back** so the user can audit and catch a mis-grouped
sample (which, given SRA provenance, they are best placed to spot): **per sample** — `ε`,
stutter level, depth, **#confident genotypes used** (homs + resolved hets, CG-seed — the
count m2(a) keys on when it is cohort-wide zero), the estimate's
uncertainty, and **soft group memberships**; **per group** — centroid parameters, sample
count, total information; **overall** — the per-sample parameter cloud (so
discrete-vs-continuum is visible at a glance) and, optionally, agreement with any protocol
labels the user *believes* hold (a **diagnostic cross-check only**, never a clustering
input).

*Determinism is preserved.* Per-sample/per-group estimation is a pure function of the
catalog-ordered, deterministic confident-homozygote reduce, and the clustering is a fixed
deterministic fit over those per-sample values — nothing is per-thread.

**The pre-pass: burn-in → measure → cluster.** *(Its three parts go **by name**, not by
number: the pipeline already numbers **Phase** 1/2/3 (reading / this pre-pass /
genotyping) and the algorithm numbers its **steps** S1–S5, so a third 1/2/3 scheme would
collide; all three parts live inside Phase 2.)* Both `ε` and the
`θ_period` shape priors come off the **confident genotypes** (homozygotes ∪ well-separated
hets — CG-seed §4.3) — a known-`(A,A)` homozygote makes every non-`A` read a labelled
stutter/error observation, and a resolved `(A,B)` het with alleles **≥ 2 units apart** does
the same per allele (outer skirts hard-labelled, inner valley → soft EM) (§4.3): the skirt(s)
give `θ = (u,d,ρ)`, the within-tract base mismatches give `ε`. They are jointly identifiable
because they move in different coordinates — **`θ` changes length (in whole motif units),
`ε` changes composition (substitutions)** — *provided* `align` admits **no sub-motif indel
inside the tract** (§6): without that rule a single-base in-tract indel is a second
length-changing mechanism that collapses onto stutter at period 1 (mononucleotides).

- **Burn-in (settle).** An **adaptive loop** that bootstraps the selection model from
  **dev-computed default parameters** — this is where selection and estimation co-evolve
  (it absorbs the "iterative confident-homozygote selection" above). Each cycle: **draw a
  batch** of loci in a **seeded random order** (reproducible *and* a representative length
  spread); a **worker-pool maps** each locus — the §4.3 seed gate + **soft per-allele
  responsibilities** (post-C2, the gate only *initializes*) — **using the batch's frozen
  parameters** (pure `(locus, params)`, no shared mutation); a **barrier → commutative
  reduce** in catalog order; then **update** the parameters.
  - **Stop point = the penalized marginal log-likelihood plateau (M4 amend. 2026-06-22),
    *not* "params/calls stopped moving."** Because post-C2 the pre-pass is a **soft
    full-cohort EM**, it ascends a single well-defined objective: `ℓ_pen` = the marginal
    (observed-data) log-likelihood of the subset reads under the current chemistry params
    **plus the log-priors** (the shape-shrinkage + `G₀` terms). EM/MM makes it **monotone
    non-decreasing**, and it is the **normalizer of the E-step we already compute**, so it
    is ~free. Stop when its **relative increment `Δℓ_pen/|ℓ_pen| < tol`** for the last step,
    with a **max-iteration cap**. This is principled (it measures *data fit*, not the
    model's agreement with itself — the circular "calls don't move" is retired), monotone,
    and deterministic — **`ℓ_pen` is summed by the same fixed-point integer accumulation as
    `F_i`/level (verify-fix #1), so the plateau test compares byte-identical values and the
    stop iteration is thread-count-independent**. A
    **non-monotone or never-plateauing `ℓ_pen`** is a sharp diagnostic — the loop isn't a
    proper ascent (a bug, or a model mismatch the covariates miss), *which assumption
    broke*, not *the method is wrong*.
  - **`ℓ_pen` plateau means *converged*, not *correct*** — a self-consistent wrong solution
    can sit at a high local optimum. So correctness is checked **outside** the loop, by two
    anchors: **multi-start** (run from several seeds/inits, compare final `ℓ_pen`, take the
    best; **divergent basins are a flag**, not silently resolved), and the simulator
    ground-truth recovery + known-protocol positive control (§9). `ℓ_pen` is the
    *convergence* check; those are the *correctness* gates.
  - Deterministic across thread counts because the seed fixes batch membership, the map is
    pure, parameters update **at the barrier** (never a per-locus running overwrite), and
    every reduce that feeds a decision — sufficient stats (integer), `F_i`/level *and now*
    `ℓ_pen` + the BIC log-likelihoods (fixed-point integer), the clustering (deterministic
    tie-breaks) — is order-independent (verify-fix #1). **Batch size** trades parallelism
    vs. adaptation granularity.
- **Measure.** From the settled value, estimate each parameter **per locus**
  over a **representative, stratified** sample (spanning the length range, not just the
  first/cleanest), and build the **distribution of those per-locus values**; the frozen
  value is the **average (mean)** of that distribution. The distribution's *shape* is a
  **development-time diagnostic** — tight & unimodal validates pooling; wide or multimodal
  flags unmodelled structure (covariate too coarse, mis-grouped sample, bug) — inspected
  during implementation, not at runtime.
  - **`ε`** is **frozen at the per-group mean**, licensed by a **penalized-likelihood
    model comparison** (M4 amend. — *not* the old "±1 SD flips no call"): compare
    `ℓ_pen(frozen ε)` against `ℓ_pen(ε + a covariate / per-sample ε)`; **freezing is
    justified iff the richer model does not beat the BIC complexity penalty** (i.e. the
    extra `ε` parameters don't earn their keep on the data). This references the data fit,
    not whether the *output* flips, so it is non-circular. If the richer model *does* win,
    that is the cue to give `ε` a covariate (a documented upgrade).
  - **`θ` shape** becomes the **per-`(group, period)` prior** for the per-locus refinement
    (length is a *level* covariate, not a grouping axis). The hierarchy partial-pools at
    every level: a sparse/heterogeneous `(group, period)` cell is **shrunk** toward its
    **cohort-per-period parent `θ_period`**, that period toward a **pooled-across-periods
    grandparent**, and an empty cell takes the literature default (global → period → group;
    M3). Distinguish **sampling noise** (few observations, consistent estimates → shrink)
    from **real heterogeneity** (ample data, still-wide estimates → a covariate is too
    coarse — the per-locus refinement, or a genuine per-group difference, then carries the
    load). The per-`(group, period)` spread vs the parent is the **invariance diagnostic**
    (M3): if groups agree it collapses to the shared shape for free.
  - **stutter level** → fit per sample as a line `level_baseline + level_slope·length`
    (both coefficients), from the **soft** per-allele responsibilities (C2 amend., not the
    hard gate); these per-sample lines are the input to the cluster part. The per-group
    level is **not frozen** — it is the **seed** (`level⁰`) for the genotyping outer loop,
    which refines it per group (C2 amend.).
- **Cluster into sample groups.** From the per-sample (`ε`, stutter-level)
  estimates, fit the data-driven soft sample groups; **freeze `ε` per group** and take the
  per-group level as the outer-loop **seed** (the full treatment — distance-based grouping
  of close neighbours, coverage as precision-only, reporting — is in the amendment above).

Because the pre-pass `θ` is only a prior, its quality matters for *convergence speed*
in low-depth loci (where the prior dominates) but not for *correctness* — so the
robust-selection hygiene of §4.3 (depth gate; flag a one-sided skirt excess, the
masquerading-het signature) is **prior-cleaning, no longer load-bearing**.

Phases 1–2 are **parallel and deterministic**: workers accumulate into per-thread
partials reduced in a fixed (catalog) order — a **commutative reduce**, order-independent.
(A locked *running overwrite* would re-introduce order-dependence; only the
*accumulation* is shared.) The keyed `θ` structure is respected throughout — one shape
prior per `(group, period)` (shrunk to a cohort-per-period parent — M3) + one stutter level
per sample group, never one global value.

> **Amendment (verify-fix #1, 2026-06-22) — the pre-pass *decision* floats reduce
> order-independently too, not just the sufficient statistics.** The integer
> `SlipProfile`/`SampleStutterStats` counts above are order-free, and the `F_i`/level
> reduces use fixed-point integer accumulation (M1). But the M4 amendment introduced three
> *new* floating-point quantities that **drive discrete decisions**, and their
> order-independence must be *specified*, not assumed (each is a float sum over loci that
> complete out of order, so naïve per-thread partials would be thread-count-dependent —
> the same non-associativity M1 fixed):
> 1. **`ℓ_pen` (the burn-in plateau stop).** The per-locus E-step log-normalizers are
>    summed by the **same fixed-point integer accumulation** as `F_i`/level (scale → round
>    → sum into an `i128`), so `ℓ_pen` — and hence the stop test `Δℓ_pen/|ℓ_pen| < tol` — is
>    **byte-identical across thread counts**, and the stop *iteration* cannot drift. (No
>    overflow: per-locus normalizers are bounded by `reads × O(few nats)`; even `10¹²` loci
>    × `10⁴` nats × `2⁴⁰ ≈ 10²⁸ ≪ i128::MAX ≈ 1.7e38`. Absolute precision `2⁻⁴⁰ ≈ 9e-13`
>    nats is far finer than any meaningful `tol`, so the comparison is exact.)
> 2. **The BIC model comparisons (ε-freeze; group split/merge; the number of sample
>    groups).** Each is a difference of two such fixed-point `ℓ_pen` values against a
>    *deterministic* complexity penalty — a deterministic scalar comparison, so a borderline
>    `Δℓ_pen ≈ penalty` resolves to the **same** discrete outcome on every thread count.
> 3. **The clustering distances + grouping.** Per-sample `(ε, level)` are deterministic
>    per-sample fits (fixed-point reduces of that sample's own labelled observations); the
>    `(ε, level)`-space distances are then deterministic. Grouping is **no random init**
>    (already specified) **plus deterministic tie-breaks**: strict `<` thresholds with a
>    defined equal-handling rule, and a **total tie-break on the sample's catalog index**, so
>    a borderline `distance ≈ threshold` (or two equidistant merge candidates) resolves
>    identically regardless of completion order.
>
> With this, "byte-identical across and within thread counts" covers the pre-pass's
> *decision layer* (stop iteration, freeze/merge, group count, group membership), not only
> its sufficient statistics. Applied to parameters §3/§4/§7; the schedule block below.

**`F` (and, since the C2 amend., the stutter level) — the prior-side outer loop
(per-individual `F_i`).** `F` is estimated from the cohort-wide heterozygote deficit, which
needs genotypes, so it cannot precede genotyping. But `F` lives in the **genotype prior**,
not in `align`, so re-estimating it is cheap and **never rebuilds the cache**: run a
genotyping sweep with the current `F` (from a supplied/default `F⁰`), reduce, re-sweep. This
nesting is itself an EM (`F` is a mixing weight, the per-locus `π` loops are its inner E/M)
— monotone coordinate ascent; stop on **`|ΔF| < tol`** with a **max-rounds cap**. The
**per-group stutter level rides this same loop** (C2 amend.): it too is outside `align`
(an `S_θ` re-weight), so it is re-estimated per group from the soft responsibilities in the
same reduce, rebuilding no cache. **Determinism of both reduces (M1).** They sum
per-(sample, locus) responsibilities, and loci complete out of order across threads, so a
per-thread float partial would be thread-count-dependent (floating-point addition is not
associative). Both therefore use **fixed-point integer accumulation** — scale each
contribution, round to an integer, sum into a per-individual / per-group accumulator;
integer addition is order-independent, so the result is **byte-identical across thread
counts** (the same `u64`-counts trick the pre-pass reduce already uses, §4.4 pre-pass). The
engine already takes a **per-sample `F` vector** as input, so the `F` half needs no
prior-side change, and the level half is the `S_θ` re-weight the EM already applies.

*The estimator (the M-step).* `F` is the **mixing weight of the IBD-mixture prior**
(§5.3): each genotype draw is autozygous with probability `F`, else an HWE draw. Its EM
update is therefore the **mean posterior responsibility of the autozygous branch** over
all (sample, locus) draws — a **heterozygote** contributes 0 (it can only be outbred); a
**homozygote** for allele `i` contributes `F·π_i / (F·π_i + (1−F)·π_i^ploidy)`. This is
the "het-deficit" reduce made exact (≈ `1 − Hₒ/Hₑ` in the clean limit). **Only
variable loci enter** (≥2 alleles — het opportunity) — a monomorphic locus is
self-consistent at any `F` and carries no information.

*Granularity — per-individual `F_i`, shrunk.* A single global `F` blends a structured
cohort (wild tomato: mixed mating systems) into a value that fits no one. We estimate a
**per-individual `F_i`** instead — an SSR catalog gives thousands of loci per individual,
so each `F_i` is well-determined — and **shrink it toward the cohort-mean `F`**
(hierarchical: a global mean + per-individual deviations, sparsely-typed individuals
shrinking harder). Per-individual is also *more robust* to an **idiosyncratic** bad
sample — a contaminated or mismapping-heavy individual inflates mainly its own `F_i`, an
outlier the shrinkage and a robust cohort mean contain. (It does **not** remove a
**cohort-wide** signal such as Wahlund structure, which genuinely elevates *everyone* — that
is reported as apparent `F_IS`, not localized away; see the caveat below.)

*The `F = 1` absorbing trap → a hard ceiling.* `F_i = 1` makes **every** heterozygous
genotype a-priori impossible for that individual — its hets can never be called however
strong the reads, and the M-step pins it there (the exact mirror of the `π = 0` trap the
§4.3 pseudocounts guard against). So `F_i` carries an **always-on hard ceiling
`F_CEILING = 0.99`** — a 100× prior down-weight on hets (overcome-able by strong
evidence), comfortably above genuine near-complete selfers (≈ 99.5 % selfing at
equilibrium, `s = 2F/(1+F)`). An **optional CLI cap** can lower it further for a more
conservative het-calling floor; the order is: raw `F_i` → shrink to the cohort mean →
clamp to the user cap if given → clamp to `F_CEILING`.

*Caveat — it is apparent `F_IS`, not pure inbreeding (a documented user warning, M5 amend.
2026-06-22).* The estimate **absorbs any cohort-wide source of homozygote excess** that
mimics inbreeding — chiefly **population structure (the Wahlund effect)**, plus paralog/CNV
mismapping and reference/mapping bias — so the output is **apparent `F_IS`**, and a naive
user must not read it as pure inbreeding. **Wahlund in particular is real population
structure, not an artifact**: in a structured cohort (wild tomato) every individual's `F_i`
is genuinely elevated, and there is **no clean within-caller correction** — so v1 does **not**
try to remove it; it **reports apparent `F_IS` and documents the caveat** so the user knows
structure is folded in (decompose it downstream with explicit structure/admixture analysis if
needed).

*On null alleles specifically.* Classical microsatellite work treats null alleles
(primer-binding-site mutations → allelic dropout → apparent hom-excess) as a major `F_IS`
confounder — but that is a **PCR-fragment-analysis (capillary) phenomenon**: it needs
locus-specific primers. This caller is **primer-free sequencing** (align reads, no per-locus
PCR), so the dominant null mechanism **largely does not apply**; the residual is the much
smaller **reference/mapping bias** (a divergent allele maps slightly worse), and
**capture/RAD libraries** can reintroduce probe-/restriction-site dropout. So per-locus
null-allele down-weighting stays **deferred** and that deferral is **appropriate for WGS**
(revisit only if validation on real capture/RAD data shows it bites — §9). The `F = 1`
trap guard (ceiling) and the variable-only reduce remain as sensible guards; the old
"self-reinforcing null loop" alarm is downgraded accordingly.

**The whole Stage-2 schedule, then:**
```
 pre-pass (confident GENOTYPES = homs ∪ separated hets, CG-seed; seeded batches, frozen-params-per-batch map/reduce, deterministic):
   BURN-IN (adaptive loop, from dev-computed default params; MULTI-START → best ℓ_pen): the
            1..ploidy-peak resolution gate SEEDS clean labels (CG-seed/C2); estimate = SOFT full-cohort EM reduce:
     repeat { draw seeded batch of loci → parallel MAP each locus (resolve genotype: hom OR ≥2-apart het, recurrent;
              het outer skirts hard-labelled, inner valley + leftovers → SOFT per-allele responsibilities, THIS batch's frozen params)
              → barrier → commutative reduce (catalog order) → UPDATE params }
       until Δℓ_pen/|ℓ_pen| < tol  (penalized marginal log-lik plateau — M4; NOT calls-don't-move)
       (ℓ_pen summed by FIXED-POINT integer accum like F_i/level → byte-identical stop iteration; verify-fix #1)
   MEASURE (seeded stratified subset): per-locus distributions → diagnostics; then:
       (cohort-per-PERIOD stutter SHAPE parent θ_period (groups not yet known), refined per-locus;
        per-sample ε + stutter level = level_baseline + level_slope·length, from SOFT responsibilities;
        shrink sparse periods & thin samples; literature default for empty loci groups;
        distribution shape = dev-time diagnostic)
   CLUSTER samples into SAMPLE GROUPS (distance-based on (ε,level), NOT k-means;
       distances scaled by 1/√(depth) [dup-free: markdup upstream, Stage 1 skips flag 0x400 — verify-fix #5];
       merge iff split fails BIC + enough-data (M4);
       coverage NOT a feature, NOT regressed out;
       BIC = Δ of fixed-point ℓ_pen vs deterministic penalty; ties broken on sample catalog index — verify-fix #1)
     → FREEZE per-GROUP ε (align cache DEFERRED — pure fn of frozen ε, recompute; M2);
       per-GROUP level⁰ = outer-loop SEED, NOT frozen (it's an S_θ re-weight — C2 amend.)
     → FIT per-(GROUP, PERIOD) SHAPE shrunk to the θ_period parent (groups now known; M3);
       per-(group,period) spread vs parent = invariance diagnostic
     → REPORT per-sample / per-group values + soft memberships
 genotyping outer loop, F_i⁰ + per-group level⁰ = supplied/pre-pass seed:
   repeat { per-locus EM on LEVEL-adjusted reads: refine π AND θ_locus shape
              (frozen per-group ε; CURRENT per-group level; θ_locus shrunk to (group,period) shape;
               align recomputed (cache deferred), θ_locus & level only re-weight it) with current F_i
            → per-individual F_i reduce (variable loci; mean IBD-responsibility;
              shrink to cohort mean; clamp ≤ F_CEILING=0.99)
            → per-GROUP stutter-level reduce (soft responsibilities; fixed-point integer
              accum = order-independent → byte-identical; no cache rebuild) }  until |ΔF|,|Δlevel| < tol / max rounds
 → genotype calls = final per-locus E-step
```
The alignments don't change as the EM iterates (because `ε` is frozen — only `S_θ`
re-weights them), and post-C1 each in-tract score is a cheap substitution closed-form, so
v1 **recomputes on demand** (the cache is deferred, §6); what iterates is cheap — `θ_locus`
and the per-group level re-weighting inside each locus's EM, and the prior-side `F_i` +
level reduces across loci. Identifiability holds: allele-vs-stutter is decided *within a
locus, across samples* (the per-locus EM pools the cohort's reads at that locus, and
`θ_locus` is anchored by its loci-group prior), so the per-locus `θ` refinement cannot run off
to absorb a recurrent real allele as stutter.

### 4.5 Output — the cohort VCF (the CALL step)  *(settled 2026-06-19)*

Mark-1 §5.9 deferred the output as "reuse." The **GangSTR/TRtools-compatible shell** does
reuse cleanly, but three semantics are **SNP-shaped and must be redefined** for SSR — so
§5.9 is **adapted, not pure reuse**.

**Reused as-is.** REF/ALT = **actual tract sequences** (Mark-2 stores them); per-ALT
**REPCN** (`= len/motif`, fractional for impure alleles — TRtools tolerates) and
**BPDIFFS** (`= len − ref_len`); per-sample **GT** indexing the cohort ALT list;
**allele pruning** (drop ALTs no called genotype uses — reuse `prune_unsupported_alleles`);
the **TRtools/GangSTR header detection** (the two-token rule). Diploid is fully
harmonizable; polyploid is valid but TRtools is diploid-centric.

**Emit criterion — variable, not "non-ref."** The reused writer drops a record that is
hom-REF in every sample (`is_variant_call`, keyed on genotype 0 = hom-REF). For SSR the
reference tract is **just a frame**, so that rule is wrong; replace it with **drop iff
monomorphic**: emit a locus **iff ≥2 alleles segregate in the cohort** (with support),
*regardless of frequency* — a **variable-but-rare** locus (major-allele freq ≥ 0.95,
MAF < 0.05) is real and **emitted**; only a truly **monomorphic** locus (one allele
cohort-wide) is dropped. **"Polymorphic"** (MAF ≥ 0.05) is a reported **AF annotation**,
not an emit or QUAL criterion.

**Site QUAL — Phred-confidence the locus is variable.** Not the SNP "P(all hom-REF)"
(meaningless when REF is a frame, and saturating at every polymorphic STR). Define
**`QUAL = −10·log₁₀ P(locus monomorphic)`** — the confidence that **≥2 *real* alleles
segregate**.

*The estimator (m4 amend. 2026-06-22) — sum the exact-AF convolution over candidate
alleles, not anchored on REF.* "Monomorphic" = the cohort is **fixed for exactly one
allele** — *any* allele, not REF — so decompose into the **mutually exclusive** events
"fixed for allele `i`," one per candidate, and sum (no inclusion–exclusion needed):
**`P(monomorphic) = Σ_{i ∈ A_ℓ} P(cohort fixed for allele i)`**. Each term reuses the SNP
engine's **exact-AF convolution kernel** (`convolve_ac_linear` + its Beta-Binomial-`K`
prior), with allele `i` **collapsed as the "reference"** and everything else as "alt":
`P(fixed-for-i) = P(total non-`i` allele copies = 0 | data)` = the engine's `P(K=0)`,
computed per-allele.

**This is a *generalization* of the engine's `compute_qual_via_exact_af`, not a call to it
unchanged (verify-fix #7b).** That function hard-codes the collapse on **allele-0** (its
per-sample Step-1 bucket is `ploidy − count[0]`, and it reads `α_ref = pseudocounts[0]`,
`α_alt = Σ pseudocounts[1..]`). So SSR must (1) **parametrize the collapse allele**
(`count[i]`, not `count[0]`), (2) per term set the Beta shape from `G₀` — `α_ref = G₀[i]`,
`α_alt = Σ_{j≠i} G₀[j]` — and (3) **bypass the REF-anchored `P(K=0)` outer wrapper**. Calling
the function *as-is* yields only `P(fixed-for-the-modal/REF allele)`, **not** the sum — the
kernel (`convolve_ac_linear`) and the Beta-Binomial/`K=0` math are reusable, the allele-0
indexing is not.

**Normalization — honest scope (verify-fix #7b).** Each `P(fixed-for-i)` is a posterior
under its **own** binary `(i vs rest)` collapse with its **own** partition `Z_i`, so the sum
`Σ_i P(fixed-for-i)` is **exact for a biallelic locus** and an **approximation for `≥3`
alleles** — the per-collapse `Z_i` differ, and lumping the non-`i` alleles into one Beta is
the *same* collapse the engine's own multiallelic QUAL already makes (exact only at
biallelic). This is fine for a **QUAL score**: it is accurate at the extremes (clearly fixed
→ the dominant term ≈ 1 and `Z_i ≈ Z`; clearly variable → every term ≈ 0), loose only in the
borderline middle, where a calibrated cutoff is what matters anyway. The **exact form is the
documented upgrade**: extract each term's *unnormalized* `K=0` numerator
(`log_p_ac_next[0]` before the `log_z` divide) and divide the summed numerators by a single
joint multiallelic marginal `Z` (a genuinely multi-dimensional AF convolution — a bigger
change than the kernel reuse). It is **better than the naive `Π_s` product** and **not
reference-anchored** (REF is just term `i = REF`; an off-reference population's true modal
allele dominates the sum and scores correctly), and the **§6 FP defenses feed it directly**:
a stutter/artifact false second allele has its `π` driven → 0 (recurrence + `G₀` +
allele-balance), so mass concentrates on `fixed-for-the-true-allele` → `P(monomorphic)`
rises → **low QUAL**. *(Cheaper proxy if `|A_ℓ|` convolutions ever bite on
a hypervariable locus: the naive `P(mono) = Σ_i Π_s P(GT_s = hom_i)` — a perf fallback only,
not the default.)*

It is the confidence in the emit decision. Deliberately **not** about polymorphism: a
variable-but-rare *real* locus scores **high**. Edges fall out correctly: a candidate set
of REF alone → `P(fixed-for-REF) ≈ 1` → QUAL ≈ 0 → dropped as monomorphic; an
all-thin/no-call locus → prior-dominated, conservatively high `P(monomorphic)` → low QUAL.
(**HipSTR writes `QUAL = "."`** and does everything per-sample, `seq_stutter_genotyper.cpp`
`write_vcf_record`; we go further because we are **cohort-joint** and already make the
variable/monomorphic decision — so a meaningful locus-level quality is free, and gives a
precision-first popgen pipeline a locus filter handle HipSTR lacks.)

**Site FILTER — locus-level reasons** (`PASS` + the SSR no-call reasons; a filtered locus
is **emitted with its reason, never silently dropped**, §5): **`notPeriodic`** (the §5
locus-admission filter — rung spacing isn't the motif), **`tooManyAlleles`**
(`|A_ℓ| > MAX_CANDIDATE_ALLELES`), **`lowDepth`** (cohort-wide insufficient coverage).
The SNP **`segdup`** filter is dropped (SSR delegates mappability to per-read MAPQ —
Mark-1 §3.1).

**Per-sample — GT, GQ, no-call.** GT indexes the cohort ALT list; **GQ / `Q`** = the
genotype posterior (TRtools' quality field, = HipSTR's per-sample `Q`), carrying the §6
allele-balance penalty. A sample is **`./.`** when **absent** (no data, §4.1), when its
reads are **outlier-dominated** (the `λ` term carries most of the mass), or when its best
genotype fails the §5.8 gates (posterior / depth / support). Per-sample DP/REPCN as in
§5.9. This mirrors HipSTR's per-sample-`Q` + per-sample-filter design, but routes
**locus-level** failures to site FILTER and **sample-level** ones to `GT = ./.`.

**`F_IS` reporting carries a user warning (M5).** Wherever the per-individual `F_i` (or a
cohort `F`) is surfaced — a VCF header comment, an `INFO`/sample field, or the side report —
it must be **labelled apparent `F_IS`** with the caveat that it **includes population
structure (the Wahlund effect)** and residual mapping artifacts, and is **not** pure
inbreeding: v1 does not decompose them (no clean within-caller correction), so a naive user
is explicitly told structure is folded in (decompose downstream with a structure/admixture
analysis if needed). This is the §4.4 caveat made an **output obligation**, not just design
rationale.

---

## 5. S1 — candidate assembly from observed sequences  *(shape settled 2026-06-19; thresholds open — §9)*

Evidence is assembled in three levels — **pool → rungs → candidates** — then a
coarse **locus-admission filter** decides whether the locus behaves enough like an
SSR to analyze at all. Throughout, **recall is the goal and precision is the EM's
job** (§4.2): a sequence dropped here is gone, whereas a wrongly-kept candidate is
driven to `π ≈ 0` by the population recurrence + `G₀` (§5.5). Sequences that do
*not* become candidates are **not discarded** — they re-enter the likelihood as
slip/error products of the candidates (the §6 invariant).

**Spanning-only is now a *structural* cap, not just a read filter.** Because candidates
are *observed* sequences with **no synthesis** (we never invent an allele the reads
didn't show — the precision-first no-assembly limit, §1), the callable allele range is
**hard-capped at the read length** by construction: an allele longer than a read leaves
no spanning observation, so it can never enter `A_ℓ`. For short plant SSR markers this is
the intended trade (expansions are an explicit non-goal, Mark-1 §1.3/§1.4). GangSTR's
beyond-read-length read classes (FRR Poisson-count, spanning-pair insert size, flanking
lower bounds) are the only way past it and stay **out of scope**; the one cheap
in-paradigm lever — **read-pair merging** to extend the spanning ceiling (Mark-1 §1.4) —
is **not in Mark-2** and is the first place to look if the read-length cap bites.

**1. `ObservedSeqs` — the cohort pool.** Sum the per-sample observed `(seq, count)`
distributions into one cohort-aggregate distribution (the per-sample distributions
are kept for level 3). This is the universe of sequences the locus shows.

**2. Rungs — the ladder scaffold + shared coordinate system.** A **rung** is an
occupied step of the "stutter + real allele" ladder. A length-position is admitted
as a rung by *either* signal, because they catch different things:
- **recurrence** — the length is observed in ≥ *k* samples; this admits
  the **stutter bands** (systematically present but never local maxima) and a
  **stutter-masked allele** (always a minor shoulder in its carriers), and it is
  what separates real structure from per-base **sequencing error**, which is
  sporadic and does not recur at a consistent position;
- **height** — it is a local maximum in *some* sample's distribution; a local
  maximum is strong evidence of a real allele even from a single sample (the cue a
  human uses reading a length trace), so this admits **rare private** alleles
  recurrence would miss.

A rung is **length-keyed and holds a *set* of sequences** — every distinct
same-length sequence whose **cohort frequency** clears a threshold (often one
sequence, sometimes several: interruption / substitution variants). Each such
sequence is an **independent allele**, not stutter-linked to the others (stutter
changes length, so it cannot connect same-length sequences). The **below-threshold
same-length cloud is per-base sequencing error** — *not* promoted to candidates, but
still fed to the likelihood as error products (§6); the frequency threshold is
exactly what separates a real same-length variant from that error cloud. The rung
**lengths** are the cross-sample coordinate system and the scaffold S2 (§7) reads.

**3. Locus-admission filter — is this locus a typical SSR?** Form the distribution
of **adjacent-rung length differences**. Its **mode must be the catalog motif
length** (the rung spacing is motif-dominated); if it is not — no coherent motif
period, competing periods, junk — the locus does **not** behave as the statistical
model assumes and is **no-called** (emitted as a *filtered* record with a reason,
never silently dropped). This is the data-driven analog of Stage-0 / GangSTR-style
catalog curation, but it additionally catches a locus that looks clean in the
*reference* yet shows non-periodic structure in the *population* — which
reference-only curation structurally cannot see. *(Open — §9: robustness to **empty
rungs** [allow integer *multiples* of the motif, not only the motif]; how dominant
"mostly" must be; the minimum rung count below which the test is skipped and the
locus admitted by default.)*

**4. `CandidateAlleles` — per-sample nomination, unioned.** For each sample, count
its **clear local maxima** over the rungs — a maximum is *clear* when it stands
**> 3 reads above each adjacent rung** (a prominence floor; default to confirm).
Then, conditioned on how well the sample resolves its own genotype:
- **≥ ploidy clear maxima** → the genotype is resolved: the **top-ploidy** clear
  maxima are the sample's candidate alleles, and every other band is stutter. **No
  ±1 rescue** — there is no hidden allele to recover.
- **< ploidy clear maxima** → the sample is under-resolved (a heterozygote whose two
  alleles are one unit apart merges into a single peak when the taller allele's −1
  stutter fills the valley). Add the **±1 neighbors** so the EM can resolve the
  possible hidden allele.

A rung holds a *set* of sequences (level 2), so a candidate is a **sequence**. The
±1 rescue **adds the sequences sitting on the adjacent rungs** (the observed rungs
one motif step away) — they are already in the cohort scaffold, so we *include*
them, we do **not** synthesize anything by adding/subtracting the motif. If an
adjacent length was never observed as a rung, there is nothing to add (no cohort
evidence of a hidden allele there). Slip arithmetic (the §7 `⊕`) belongs to the
likelihood, not to nomination — so no impure-placement question arises at this step.

Union across samples. The **reference allele is seeded unconditionally** (the VCF
needs a REF, `G₀` needs the reference length, hom-ref must be callable). Cap the
union at `MAX_CANDIDATE_ALLELES`; exceeding it → no-call (hypervariable / noisy).
**Off-lattice (impure) major peaks are kept**, not discarded — within an admitted
locus they are first-class candidates handled **per-allele** by the HipSTR machinery
(§6 / §7), and treated **identically to pure alleles** (no impurity penalty — §7).
The level-3 locus filter remains the only hard discard; individual impurity is
modelled, not penalized.

## 6. S3 — the read likelihood `Qᵣ(obs_seq | candidate)` — HipSTR-informed  *(shape settled 2026-06-19; details open — §9)*

We adopt HipSTR's **sum-over-slips generative form** — the law of total probability
over the latent slip size — adapted to Mark-2's collapsed-sequence evidence (distinct
observed sequences + counts, uniform quality — no per-read rows, no base qualities). The
probability of an observed sequence given a candidate allele **sums over the PCR slips**
the allele could have undergone:

> `Qᵣ(obs_seq | candidate) = Σ_Δ  S_θ(Δ) × align(obs_seq | candidate ⊕ Δ)`

- **`Δ`** is the **slip size** — how many whole motif units were added/removed,
  weighted by `S_θ(Δ)`. For an *impure* candidate, *where* those units land in the
  tract is ambiguous (slipping `(CA)₅TA(CA)₃` up gives `(CA)₆TA(CA)₃` *or*
  `(CA)₅TA(CA)₄`), so `candidate ⊕ Δ` is a **set of placement variants**. Following
  HipSTR, the **placement is marginalized by an explicit sum over those variants**
  (the faithful borrow — verify-fix #3, 2026-06-22):

  > `align(obs | cand ⊕ Δ) = Σ_{v ∈ placements(cand, Δ)}  Pr(v) · align_subst(obs | v)`

  where `placements(cand, Δ)` is the set of **placement-distinct** sequences obtained
  by adding/removing the `|Δ|` whole units in each repeat **run** of the tract,
  `Pr(v)` is a **uniform position prior** over placements (HipSTR's
  `1/(block_len + Δ + 1)`), and `align_subst(obs | v)` is the in-tract
  substitution score for that *one fixed* variant (below). The key bound: the variants
  are distinct **only across runs separated by interruptions** — adding/removing a unit
  *within* one homogeneous run yields the **same string**, so the sum has
  **≈ (#interruptions + 1) terms, not (#copies) terms** (HipSTR collapses the equal-LL
  placements within a run via its `upstream_match_lengths_` run-length trick;
  `StutterAlignerClass.cpp:84-87/133-136`). So a singly-interrupted allele costs **~2–3
  evaluations**, not a DP sweep. For a **pure** candidate the placement is degenerate
  (every run is one homogeneous block ⇒ one variant), the sum collapses to a **single**
  `cand ⊕ Δ`, and the substitution **closed-form / exact-match fast path** below applies
  directly — so pure alleles (the majority) pay nothing for this.
- **`S_θ(Δ)`** is the **stutter kernel** (§5.2 reuse): `θ_locus`, refined per locus
  under its per-loci-group prior (§4.4), the **same for pure and impure alleles** at a locus
  (slippage is a property of the tract, not the allele — as in HipSTR).
- **`align_subst(obs_seq | v)`** (one fixed placement variant `v`) is a **banded
  pair-HMM forward** with a **flat (uniform-quality) emission** — a *probabilistic*
  alignment that **sums over all the ways** the observed sequence could arise from `v`
  under a per-base error rate `ε`, yielding a genuine probability `P(obs | v)` (not a
  best-path score). It reuses the **Stage-1 SSR pair-HMM** machinery (banded,
  scratch-buffered; pattern from BAQ `probaln`), with flat emission because Mark-2
  dropped base qualities (the Stage-1 gate made survivors uniform-quality). **The slip
  placement is marginalized by the outer `Σ_v` (above), not inside this term** — each
  `align_subst` scores `obs` against *one* fixed variant; the variant set + uniform
  prior do the marginalization (verify-fix #3). An **exact-match fast path**
  (`obs == v` byte-for-byte — the clean post-gate majority) returns `(1−ε)^len`
  without running the HMM, so the HMM fires only on the error-bearing / impure
  minority. *(Affine-gap best-path was considered and rejected: it gives a score, not
  a probability, and commits to one alignment instead of marginalizing —
  inconsistent with the `Σ_Δ` sum and the placement marginalization; its only edge,
  speed, is moot since `align` is cacheable.)*

**In-tract errors are substitutions only — length changes are stutter's job (HipSTR's
repeat-block rule).** `align` admits **no sub-motif (per-base) indel inside the repeat
tract**: the *only* length-changing operations on the tract are (a) the whole-motif slip
size carried by `Σ_Δ` and (b) the whole-motif *placement* carried by `Σ_v` (above) —
both whole-unit, never sub-motif. **Within a single fixed placement variant `v`**
`align_subst` scores **substitutions only** (`(1−ε)^match · (ε/3)^mismatch` over the
aligned bases), with per-base gaps confined to the **flanks** (boundary slop / misplaced
flank). This is exactly how HipSTR
walls its stutter block off from sequencing indels — in-block insert/delete states are
`IMPOSSIBLE` and a stutter block must be followed by a match
(`HipSTR/src/SeqAlignment/{HapAligner,StutterAlignerClass}.cpp`). It is what keeps `ε`
and stutter **identifiable**: they live in genuinely different coordinates — **`ε`
changes composition, stutter changes length** — at *every* period, period 1 included (a
substitution in a homopolymer makes an interruption, not a length change; the length
change is the slip). A sub-motif indel that *does* occur in the tract is not modelled as
in-tract noise: if it breaks the motif it is an **interrupted (impure) allele** (a
first-class candidate, §5/§7); a single-base slip in a homopolymer simply *is* a
whole-unit slip and belongs to `Σ_Δ`.

**Period 1 (mononucleotides) — flagged, not specially modelled.** Here a repeat unit *is*
one base, so even with the no-gap rule the stutter term legitimately subsumes **both** PCR
slippage **and** the indistinguishable single-base sequencing indel — which is the
correct, and empirically the only separable, decomposition. Period-1 loci are therefore
**flagged**: their `ε` carries substitutions only and is not separately interpretable as a
length-error rate, and their stutter level absorbs all single-base length change.
Mononucleotides are a low-priority class for this caller, so the flag *documents* the
reduced interpretability rather than adding a model; the homopolymer-length-as-gap-process
treatment (HipSTR's Dindel form, a single in-tract indel-length model whose geometric
extension *is* the stutter decay `ρ`) is the **deferred upgrade** if mono SSRs ever matter.

**Provenance — HipSTR-*informed*, not HipSTR-identical.** Genuinely from HipSTR: the
**sum-over-slips** form, the **slip-placement marginalization** (our explicit `Σ_v`
run-collapsed placement sum, verify-fix #3 — HipSTR does it via its in-block position
sum), and the
**sequence-keyed first-class impure alleles** (§7). **Our deliberate simplifications**
(so don't read "HipSTR" as "HipSTR-exact"): (1) the stutter kernel is the **3-param
up/down/decay *in-frame*** form (GangSTR's shape) — HipSTR additionally carries an
**out-of-frame, bp-step** term, which we **fold into per-base error / off-ladder
alleles** rather than model as stutter (the deliberate simplification inherited from
Mark-1 §5.2/§12; add it back only if validation shows genuine out-of-frame slippage
being misattributed); (2) `align` is a **flat-emission forward** — HipSTR's emission is
**per-base-quality** and its flank alignment is **Viterbi (max-path)**, more accurate
per read but not cacheable; flat emission is the trade that buys our distinct-sequence
cache (see "Why HipSTR cannot cache it cheaply" below).

**The outlier component (`λ`) — for *random* junk only.** `Qᵣ` scores *one* observed
sequence against *one* candidate; the **genotype** likelihood combines it over the
genotype `G`'s alleles **plus a uniform outlier term**: `P(read | G) = (1−λ)·[mix over
G's alleles] + λ·junk`. The **junk term is uniform over the locus's distinct observed
sequences** — `junk = 1/D` where `D` is the number of distinct sequences seen at the
locus (cohort-wide). Two consequences of that choice: it is **normalized and
data-defined** (no astronomically-small `4^L` floor that `λ` could never outweigh), and
its *flatness* is exactly why `λ` catches **random** junk but not **systematic** signal —
a read matching no candidate is explained as well by `1/D` as by anything, so it goes to
junk; but a *concentrated* artifact piling up at **one** sequence is explained far better
by a real allele (≈ `(1−ε)^len`) than by `1/D`, so `λ` lets it through (→ the
allele-balance term below, not `λ`, must catch it). `(λ, junk)` is degenerate up to the
product, so we fix `junk = 1/D` and calibrate `λ`. Without the term at all, every oddball
read — and these repetitive, hard-to-assemble tracts will throw some — must be explained
as evidence for *some* allele, the cheapest explanation being an **extra allele**,
inflating **false heterozygotes**. `λ` is a conservative **fixed small value in v1**
(simulator-calibrated); estimable later **but capped** — too high dumps real minor
alleles (het *deflation*), too low lets random het inflation return.

**Depth-driven het inflation, and the allele-balance term that bounds it.** The
genotype likelihood is an **i.i.d. product over reads**, so evidence accumulates
*linearly with depth* — which means a **persistent minor fraction** at a consistent
length (the systematic junk `λ` can't catch, or stutter the kernel under-models) is
scored as a real allele with **confidence that grows as coverage rises**: a false
heterozygote whose GQ *inflates with depth*. This is the **same failure mode already
seen on the SNP path** (FP QUAL climbing with depth at a stable ~20 % VAF). The i.i.d.
multinomial is **under-dispersed**: it has no parameter for "systematic minor signal,"
so it must call one a real allele.

The fix mirrors the SNP caller's QUAL refinement
([`src/vcf/qual_refine.rs`](../../../src/vcf/qual_refine.rs), the allele-balance
binomial-tail penalty — owed to SNP **GQ** too): an **allele-balance term** that asks
whether each candidate genotype's *expected* per-allele read split (≈ 1/ploidy for a
het) is consistent with the *observed* split at this depth, and **down-weights** a
genotype whose split is implausibly skewed (an 80/20 "het"). It feeds **both the
per-sample GQ and the site QUAL**, so a depth-inflated false het is driven to **low
GQ → no-call** (precision-first), not emitted confidently. Equivalently it can be
written as **overdispersion** — a Dirichlet-multinomial on the read-to-allele
assignments with a dataset-wide concentration parameter (pre-pass, like `ε`) — so
evidence **saturates** instead of growing linearly with depth; the allele-balance
penalty is the lightweight v1 form, overdispersion the documented upgrade.

**The SSR wrinkle — balance on *deconvolved* fractions, not raw lengths.** Stutter
legitimately unbalances the *observed* length tally (a true `a`/`a+1` het has its
length-`a` pile inflated by `a+1`'s −1 stutter), so the balance test must run on the
**post-stutter per-read responsibilities** the E-step already produces (which allele
each read is attributed to *after* `S_θ`), where a real het sits near 1/ploidy — not on
the raw `(seq, count)` tally.

Together: `λ` absorbs **random** junk at the per-read level, the recurrence filter
removes junk at the **cohort** level (§5), and the **allele-balance / overdispersion**
term stops **systematic concentrated** junk from inflating false hets at depth. All
three are v1 structural FP protection, not optimizations.

**Why the per-iteration cost is small.** `Qᵣ` depends on two things, treated
differently (§4.4): the per-base error `ε` (inside `align`) is **frozen**, so the
**`align(obs | seq)` values are invariant across EM iterations** (only `S_θ(Δ)` re-weights
them as `θ` refines per locus — a handful of numbers). So `P(read | a) = Σ_Δ S_θ(Δ)·align(obs | a ⊕ Δ)`
is a cheap re-weighted sum over **iteration-invariant** `align` terms, and there is **no
per-thread state and no `δ`-gated rebuild** — every worker computes `align` against the same
frozen `ε`, so it is a pure function (deterministic whether recomputed or cached). **Whether
to materialize a cache of those invariant `align` values is a deferred, unmeasured
optimization** (*resolved 2026-06-22*): with `ε` frozen they *are* cacheable (textbook EM
practice — cache the `ε`-conditionals that don't move while `θ`/`π` do), but post-C1 the
in-tract score *within a fixed placement* is a substitution closed-form
(`(1−ε)^match·(ε/3)^mismatch`, exact-match fast path) — a **pure** allele is one such
evaluation, an **impure** allele a **small run-collapsed placement sum** (~2–3 terms,
verify-fix #3), still far cheaper than a DP sweep — so v1 **recomputes on demand** and we
**measure before building any cache** — the
first step then being a per-locus memo, only later a persistent table (§9). The factorization
itself is the law of total probability over the latent slip size — HipSTR's own generative
model.

**`ε` is estimated once and frozen (§4.4), not an in-loop EM parameter.** We
**estimate** `ε` from the data rather than hard-coding it (a wrong constant biases
every alignment), but we do it in the §4.4 pre-pass — alongside `θ`, on the confident
homozygotes — and then hold it constant through genotyping. It is identifiable
alongside `θ` because the two live in **different coordinates**: `θ` changes **length**
(whole motif units, `Δ`), while `ε`/`align` changes **composition only** (substitutions —
the in-tract no-gap rule above means `ε` cannot move length), so neither can absorb the
other's signal. (At period 1 a unit is a base, so the single-base *length* change is
stutter's, not `ε`'s — mononucleotides are flagged, above.) `ε⁰` is seeded from the Stage-1 quality gate
(`ε = 10^(−Q/10)` at the gate's quality floor; `probaln` defaults for the gap params)
as the pre-pass **burn-in start**, then settled; the simulator calibrates the seed.
Because `ε` is constant during genotyping (and `θ_locus` enters only as a cheap
re-weight), `align(obs | cand⊕Δ)` is a **pure, iteration-invariant function**, so the
calls are **deterministic across and within thread counts** *whether `align` is recomputed
or cached* — the Stage-1 byte-identity invariant carries into Stage 2 (no per-thread state,
no `δ`-gated rebuild; the genotyping pass copies the frozen `ε` into each worker, and
`θ_locus` is a deterministic per-locus quantity). `ε` is a **per-sample-group** value (one
per chemistry cluster, §4.4 amend. — *not* one cohort-global scalar) with **no per-locus
shrinkage target**: a wide *within-group* per-locus spread is not something shrinkage can
fix — it signals that the flat-`ε` assumption is strained (the cue to give `ε` a covariate,
deferred). *(This supersedes the earlier per-thread-`ε` + `δ`-rebuild
design, which traded determinism for a synchronization the pre-pass makes unnecessary;
`θ`, by contrast, is **not** frozen — it refines per locus, §4.4.)*

**Cost / risk (measure, don't assume).** What to benchmark **before** adding any
optimization: (1) the per-iteration `align` **recompute** cost — post-C1 the in-tract part
is a substitution closed-form with an exact-match fast path (pure allele) or a small
run-collapsed placement sum (impure allele, verify-fix #3), so it is plausibly cheap, but
*measure it* on a hypervariable / impurity-rich locus before deciding a cache is worth its memory (a cache
is **deferred**, not in v1 — §6 / §9); (2) the §4.4 pre-pass cost (how many loci the
burn-in needs to settle `ε` and the per-loci-group `θ` priors, and the measure-part window)
against the whole-cohort genotyping pass it amortizes over; plus the per-locus `θ` M-step's
added inner-loop iterations. The stutter / base-error **factorization** is an approximation
(the same one HipSTR makes), made safe by the composition-vs-length separation above (C1).

**Why our `align` is *cacheable* (a property we hold in reserve, not a v1 build).** HipSTR
keeps **per-base qualities** and works on **raw reads**, so its alignment emission is
*read-specific* (a low-quality base scores differently) — there is no flat `align(seq | seq)`
to reuse, and no distinct-sequence collapse to amortize one over. Our two Mark-2 decisions —
**drop base qualities** (uniform quality, §2) and **collapse reads to distinct sequences**
in Stage 1 — turn `align` into a flat function whose value is **invariant across EM
iterations** (with `ε` frozen, only `S_θ` re-weights it). That invariance is what would make
a cache *correct and cheap if we built one* — but **v1 does not build it**: it recomputes on
demand and we measure first (§6). We **buy this property with fidelity**: HipSTR's
quality-weighted per-read HMM is more accurate per read; our flat-error collapsed model is
coarser but its `align` is a pure, iteration-invariant function — the same fidelity-for-scale
trade the rest of the pipeline makes, and the reason caching (when measured to be worth it)
is a drop-in later.

*(Resolved: `align_subst` = flat-emission banded pair-HMM forward (one fixed placement);
slip placement marginalized by the **explicit `Σ_v` run-collapsed sum** with a uniform
position prior (verify-fix #3), degenerate to a single term for pure alleles; `ε`
**frozen** by the §4.4 pre-pass, `θ` a **per-loci-group prior
refined per locus** (above). Open — §9-S3: the gate-derived `ε⁰` burn-in value and simulator calibration;
the pre-pass `ℓ_pen`-plateau settle tolerance (M4) and measurement window; banding / length-pruning of the
candidate × obs × slip product.)*

## 7. S2 — stutter reachability between sequences — per-allele  *(shape settled 2026-06-19; details open — §9)*

Stutter is a **per-allele** operation, not membership of a global length lattice.
The reachability relation is:

> **`B = A ⊕ k`** — `B` is `A` after `k` whole motif units are added/removed
> **within `A`'s own repeat context, the interruption structure preserved.** For an
> impure `A` this is a **set** of placement variants (the units may land in different
> runs); §6's `Σ_v` marginalizes over them. For a pure `A` it is a single sequence.

Each candidate — pure or impure — therefore anchors **its own stutter ladder**. A
**pure** allele's ladder coincides with the global motif lattice (that lattice is
the *emergent special case*, not the architecture); an **impure** allele's ladder
runs parallel, carrying its interruption along. This is what makes impure alleles
first-class (the HipSTR property): there is no "off-ladder" branch, only candidates
that happen to carry an interruption.

- The kernel's whole-unit step **`δ` is the motif-count delta along the allele's own
  ladder** (not a global length offset), so it is well-defined for impure alleles.
- The **catalog motif** is the unit of the slip; the §5 locus-admission filter has
  already established that the motif is the dominant period before we reach here.
- Impure alleles are treated **identically to pure ones — no impurity penalty.** The
  empirical candidate set (an impure allele is in it only because it recurred/peaked),
  the HMM likelihood (which scores non-candidate sequences as errors), and population
  recurrence (which drives non-recurrent candidates' π → 0, purity-agnostic) already
  filter spurious alleles; an extra prior penalty would be redundant *and* would bias
  against genuine interrupted alleles (e.g. *FMR1* AGG). If measurement later shows
  recurrent *systematic* impure artifacts leaking through as false calls, a penalty
  can be added then, with a measured basis — not by default.

**Building S2 per-allele from the start is the load-bearing choice:** strict
behaviour (admit only pure alleles) would be a *policy* over the same machinery, and
admitting impure alleles costs only the penalty above — not a re-architecture of
S2 / S3. **Slip-site placement is *not* committed by S2:** S2 only declares that `B`
is `k` units from `A`; *which copy* the unit lands on for an impure allele is
**marginalized in §6 by an explicit sum over placement-distinct variants** — the runs
between interruptions, ≈ (#interruptions + 1) terms, uniform position prior (HipSTR's
run-collapsed `Σ_v`, verify-fix #3) — so the data picks it, and a pure allele's sum is a
single term. *(Open — §9-S2: whether reachability needs the reference frame as context;
the impure-allele penalty's form.)*

---

## 8. What Mark-1 machinery is retired

The Mark-1 `Allele::OnLadder/OffLadder` enum, `build_rungs`, the reference-ladder
window, the Stage-1 two-tier `Qᵣ` (Mark-1 §4.2) and its `(length, loglik)` /
off-ladder columns. `ssr_genotyping.md` §4.2/§4.3/§5.1 are owed the amendment once
§5–§7 settle.

**Note — impure alleles are not retired, only the Mark-1 *representation* is.** The
off-ladder *enum + normalized key* is gone; impure alleles return as first-class
candidates via the **per-allele HipSTR reachability + likelihood** (§7 / §6) —
"off-ladder done right," anchored on each allele's own stutter ladder rather than a
bolt-on key.

---

## 9. Open premises (live agenda)

**Settled since the last revision (2026-06-19):** the reading/orchestration layer
(§4.1), the model framing / ingredient list (§4.2), and the **candidate-assembly +
HipSTR likelihood shape** (§5–§7: three-level pool → rungs → candidates, the
adjacent-rung-difference locus-admission filter, per-allele stutter reachability, and
the sum-over-slips likelihood), and the **EM seed + π prior** (§4.3: every-sample
putative genotype → π⁰ tally, θ⁰ from confident **genotypes** (homs ∪ separated hets — CG-seed), `G₀` pseudocount floor),
and the **shared-parameter estimation architecture** (§4.4: a three-part pre-pass
freezes `ε` and seeds a per-loci-group `θ` prior refined per locus; `F` is a prior-side outer
loop; the whole pass is deterministic), and the **output** (§4.5: emit variable loci,
site QUAL = Phred(locus variable), SSR FILTER reasons).
**Amended 2026-06-20, refined & settled 2026-06-21:** the §4.4 chemistry premise — the
stutter **level** and `ε` are **per sample group** (not cohort-global), inferred as
**data-driven soft sample groups** (coverage as **precision-only**, not a
feature/regressor; an **iterated** confident-homozygote pre-pass). Stutter is factored
`shape(group, period) × level(group, length)`: the **shape is per `(group, period)`,
shrunk to a cohort-per-period parent** (M3 amend. 2026-06-22 — *not* assumed cohort-shared;
the PCR-free "rate not profile" expectation is realized via shrinkage, so invariance is
data-driven), refined per locus toward the sample's group-period shape;
the **level is per sample group, linear in repeat length** (a slope, never bins) and —
since the **C2 amendment (2026-06-22)** — **refined in the prior-side outer loop, not
frozen** (it is an `S_θ` re-weight outside `align`); **`ε` is frozen per sample group** (it
lives in `align`). The pre-pass estimator is the **soft full-cohort EM responsibility
reduce**, with the confident-**genotype** gate (homozygotes ∪ well-separated hets,
1..ploidy-peak resolution — CG-seed 2026-06-22) demoted to a **seed** (C2). See the §4.4
amendments + the "Sample groups" item below.

The items below are the **details those settled shapes leave to flesh out.**

**S1 — candidate assembly (§5).**
- Thresholds: the rung recurrence count *k*, the per-sample clear-maximum prominence
  floor (proposed: **> 3 reads above each neighbor**), the same-length cohort
  frequency threshold, and `MAX_CANDIDATE_ALLELES`.
- The locus-admission filter: robustness to **empty rungs** (test integer
  *multiples* of the motif, not only the motif), how dominant the motif mode must be,
  and the minimum rung count below which the test is skipped (admit by default).
- Cap accounting: do the now-kept impure candidates (§7) blow `MAX_CANDIDATE_ALLELES`
  faster? (*Resolved:* no *separate* stricter bar for impure candidates — they go
  through the identical path; §7.)

**S2 — stutter reachability + `δ` (§7).**
- *Resolved:* slip-site placement for an impure sequence is **marginalized in §6 by an
  explicit `Σ_v` sum over placement-distinct variants** (the runs between interruptions,
  uniform position prior, run-collapsed — HipSTR, verify-fix #3), not committed by S2,
  not collapsed to a single sequence, and not pre-enumerated as candidates.
- *Resolved:* **no impurity penalty** — impure alleles treated identically to pure;
  the empirical candidate set + HMM + recurrence already filter spurious alleles. Add
  a penalty only if recurrent systematic impure artifacts are *measured* to leak.
- Whether reachability needs the **reference frame** as context, or is a pure
  per-allele sequence operation.

**S3 — HipSTR likelihood (§6).**
- *Resolved:* `align(·)` = **flat-emission banded pair-HMM forward** (reusing the
  Stage-1 SSR pair-HMM), with an exact-match fast path; affine-gap rejected.
- *Resolved:* `ε` is **estimated once and frozen** by the §4.4 pre-pass
  (confident-**genotype** phase — homs ∪ separated hets, CG-seed — + sensitivity gate), held constant through genotyping — so
  `align` is a **pure, iteration-invariant function** and the calls are **deterministic
  across and within thread counts** *whether recomputed or cached* (the cache is **deferred**
  — §6). `θ` is **not** frozen: the pre-pass gives a per-loci-group
  prior and each locus refines a `θ_locus` shrunk toward it (§4.4) — cheap because `θ`
  only re-weights the cache, deterministic because it is per-locus (not per-thread), and
  it makes seed contamination low-stakes (Issue 5). Identifiable because `θ` changes
  length (whole units) and `ε` changes composition (substitutions); `align` admits **no
  sub-motif indel inside the tract**, so the separation holds at period 1 too
  (mononucleotides flagged — §6). (Supersedes the earlier per-thread-`ε` + `δ`-rebuild
  design.) Open: the gate-derived `ε⁰` burn-in value, the `ℓ_pen`-plateau settle tolerance
  + multi-start count (M4), the phase-2 measurement window, and the
  `θ_locus`←`(group,period)`←`θ_period` shrinkage strengths (M3).
- Cost control (measure): the per-iteration `align` **recompute** cost (cache **deferred**
  — §6); banding / length-pruning of the `candidate × obs × slip` product; and, *only if
  the recompute is measured to be the wall*, a cache (per-locus memo first, then a
  persistent table bounded by per-group-`ε` quantization).

**Cross-cutting.**
- **Engine reuse boundary** (§3 / §4.2 / §5.5): *Resolved* — reuse
  `posterior_engine.rs` for the per-locus **E-step** + **IBD-with-`F` prior**;
  **replace** its SNP class-based pseudocounts with the per-candidate geometric `G₀` by
  **generalizing the engine to accept a pseudocount vector** (the psp-container
  precedent: extract the reusable core, write the thin SSR-specific part — don't contort
  SSR onto the SNP scheme); **bypass** the SNP-only machinery (compound alleles, chain
  anchors, contamination). Pre-alpha, no backwards-compat: the regression gate is the
  SNP caller's end-to-end tests, not byte-identity. *(Open: extend the shared engine
  in place vs a slim SSR fork of the π-loop — a struct-shape call for the architecture
  doc; the math boundary above is fixed either way.)*
- Base measure `G₀` (§5.5): *Resolved (§4.3)* — a **geometric** length-offset decay
  **centred on the per-locus cohort modal allele** (reference = the `Δ` frame only, not
  the prior's peak — removes the reference bias Mark-2 exists to kill; degrades to the
  reference when it is the only candidate), symmetric v1, **purity-agnostic** (no
  impurity down-weight, §7), with its **decay parameter fit from pooled data** per
  loci group (parametric, not non-parametric); Gaussian is the documented upgrade.
  It seeds π⁰ as well as regularizes. Open: mode vs mean as the centre (mode chosen —
  more robust); multimodal cohorts left to the data, not the prior.
- Working-set / memory: Stage 2 holds, per locus, observed sequences + per-locus EM
  scratch (the `obs × candidate × slip` scratch of whatever is being scored, reused via
  scratch buffers) — **no persistent `align` table in v1** (caching is **deferred**, §6:
  recompute on demand, measure first). The cohort working set is bounded by §4.1's lockstep
  rule (≈ N × one block). **Open / measure-first:** is the per-iteration `align` recompute
  actually a bottleneck (post-C1 it's a substitution closed-form)? — only if so does a
  cache earn its memory, first as a **per-locus memo** (no cross-sample `ε`-key → no
  blow-up), later a persistent table bounded by per-group-`ε` quantization. The whole §6
  speed picture is **unmeasured** — that's *why* the cache is deferred rather than designed
  in.
- **Stutter shape covariate** (§4.2 / §4.4): *Resolved, amended 2026-06-22 (M3)* — the
  stutter **shape** is **per `(group, period)`, shrunk toward a cohort-per-period parent
  `θ_period`** (length is *not* a shape covariate — it drives the level), built by the §4.4
  pre-pass (confident-**genotype** skirts — homs ∪ separated-het outer skirts, CG-seed — accumulate per `(group, period)` via a deterministic
  commutative reduce, *after* clustering fixes the groups), then **each locus refines a
  `θ_locus` shrunk toward the sample's group-period shape** in the per-locus EM (depth-driven
  partial pooling — Mark-1's relief valve). `θ` enters `align` only via `S_θ(Δ)`, so the
  refinement is a cheap **re-weight** (no rebuild), and it is **deterministic** (per-locus
  /per-group, not per-thread). **Shrinkage forms the hierarchy** (`(group, period)` → period
  → global; empty cells take a literature default), so a group's shape is its own only when
  the data earn it — invariance is **data-driven, not assumed**, and the per-`(group,
  period)` spread is the invariance diagnostic. So **per-group + per-locus `θ` is in v1**
  (the partial-pooling refinement), not deferred — open parts: the `θ_locus`←group-period
  and `(group, period)`←period shrinkage strengths. **Motif** and **purity** covariates deferred —
  motif identity (AT-di / homopolymers) is the documented per-locus shape outlier the
  per-locus refinement already absorbs; purity (interruptions stabilize a tract → less
  stutter) is the Mark-2-specific one, added once impure alleles are common enough to estimate.
  **Amended 2026-06-22 (C2, M3):** the shape prior is **per `(group, period)`, shrunk to a
  cohort-per-period parent** (M3 — *not* cohort-shared; groups can differ when the data say
  so); a **per-sample-group stutter level** (a `level_baseline` + `level_slope`, i.e.
  `level_baseline + level_slope·length`, both per group, pre-pass-**seeded** and then
  **refined in the prior-side outer loop** alongside `F` — *not* frozen; it's an `S_θ`
  re-weight outside `align`, so refining it rebuilds no cache and lets the cohort overrule a
  contaminated level) rides on top. Protocol/preservation move the stutter *rate* strongly
  and the *profile* weakly (PCR-free invariance) — but, since that evidence is human PCR, the
  profile's group-invariance is left to **shrinkage to decide** (M3), not assumed; the rate
  grows linearly with length. Per-locus `θ_locus` is then refined on **level-adjusted**
  reads toward the sample's group-period shape.
- **Per-base error rate** (S3 / §4.2 / §4.4): *Resolved* — `ε` is **estimated** in the
  §4.4 pre-pass (it does **not** trade off against `θ`: `θ` changes length in whole
  units, `ε` changes composition via substitutions — `align` carries no sub-motif indel
  inside the tract, §6, so this holds at period 1 too; mononucleotides flagged) and
  **frozen** for genotyping. **Per sample group** (data-driven
  soft clusters) — *not* one global value (PCR vs PCR-free / fresh vs
  historic DNA move `ε` via deamination damage). At **sample-group** granularity (not per-sample)
  *because* `ε` lives in the `align` cache → it is built once **per (locus, sample group)**;
  per-sample `ε` is the documented upgrade. A period covariate is the further
  upgrade if within-group spread stays wide.
- **Sample groups** (§4.4 amend., settled 2026-06-21; arch Q-P6): *Resolved* — chemistry
  (`ε`, stutter level) is grouped **mainly for data sufficiency** (a low-coverage sample
  borrows strength from similar samples), inferred from the per-sample (`ε`, stutter-level)
  estimates by **distance-based grouping of close neighbours — not k-means** (similar
  params ⇒ similar EM behaviour ⇒ safe to pool; deterministic, no random init), with
  **distances scaled by each sample's uncertainty** (a wobbly thin sample is "close to
  many things," never confidently misplaced). The number of groups **falls out of a
  penalized-likelihood comparison** — split two regimes only when `ℓ_pen(split) −
  ℓ_pen(merged)` clears the BIC penalty for the extra parameters, and each group has enough
  data (M4 amend.; *not* the old "calls don't move"; not a preset K; a single-protocol
  cohort collapses to one). **No cache budget applies** — the `align` cache
  is **deferred** (§6 / genotyping Q-G3), so grouping draws no cache lines at all. (The
  earlier "cache can be per-sample because shrinkage keeps distinct `ε` small" reasoning was
  wrong — shrinkage collapses only *thin* samples, so a deep continuum kept ≈`N` distinct
  `ε`; deferring the cache moots it, and *if* a cache is ever built, per-group-`ε`
  quantization bounds it.)
  **Depth enters only as measurement precision** (distance scaling), never a
  feature/regressor; it is **duplicate-free by construction** — duplicates are marked
  upstream (post-mapping markdup) and **Stage 1 skips flag-`0x400` reads**, so `ssr-call` does
  no dedup and assumes none (verify-fix #5, §4.4); the confident-**genotype** selection **co-evolves
  with the params in the burn-in loop**; all per-sample/per-group values are **reported**
  for audit. *Remaining (the correctness gates, M4):* the **discrete-vs-continuum experiment**
  (plot the real (`ε`, level) cloud — clusters or gradient?) + a known-protocol **positive
  control** (groups recover protocol, **not** coverage) + **multi-start `ℓ_pen` agreement**
  (divergent basins flagged) — these external anchors, *not* "calls don't move," are what
  validate the grouping; whether `level_slope` later proves poolable to the cohort; whether `ε` and stutter
  ever want *different* groupings (jointly clustered in v1). **Duplicate-free input by
  construction** — markdup runs upstream (post-mapping) and Stage 1 skips flag-`0x400` reads;
  `ssr-call` does no dedup (verify-fix #5, §2/§4.4).
- **`F` granularity** (§4.2 / §4.4): *Resolved* — **per-individual `F_i`**, estimated
  by the §4.4 **prior-side outer loop** (mean-IBD-responsibility reduce over variable
  loci, shrunk toward the cohort mean; sits in the prior, never touches the `align`
  cache). Captures structured / mixed-mating cohorts without population labels and
  localizes a bad sample's effect; the engine already takes a per-sample `F` vector, so
  no prior-side change. Convergence `|ΔF| < tol` + max-rounds cap. **Hard ceiling
  `F_CEILING = 0.99`** (no `F=1` absorbing trap — the `F`-analog of the `π` pseudocount
  floor) + **optional lower CLI cap**. **Per-locus `F` rejected** (confounded with π).
  Caveat: the estimate is **apparent `F_IS`** — it absorbs cohort-wide hom-excess, chiefly
  **Wahlund/population structure** (real, documented as a user warning, *not* corrected —
  M5); null-allele contribution is largely N/A for primer-free sequencing (below).
- **Null alleles & `F_IS` (M5 amend. 2026-06-22).** Null alleles (allelic dropout → apparent
  hom-excess → inflated `F_IS`) are a **classical-microsatellite (PCR-fragment / capillary)
  confounder** — they arise mainly from **primer-binding-site mutations**, which need
  locus-specific primers. This caller is **primer-free sequencing** (align reads, no per-locus
  PCR), so the dominant mechanism **largely does not apply**; the residual is the much smaller
  **reference/mapping bias** (a divergent allele maps slightly worse), and **capture/RAD**
  libraries can reintroduce probe-/restriction-site dropout. So a dedicated per-locus
  null-allele mechanism (Micro-Checker-style hom-excess inconsistent with the cohort `F`, or
  an explicit null-allele frequency) stays **deferred — and that deferral is appropriate for
  WGS**; revisit only if validation on real **capture/RAD** data shows it bites. The honest,
  largely-unfixable residual in `F_IS` is **Wahlund/population structure** (above), handled by
  a **documentation warning to users**, not a within-caller correction. (The empirical
  null-allele literature — Brookfield 1996; van Oosterhout et al. 2004 / MICRO-CHECKER;
  Chapuis & Estoup 2007 — is all *capillary microsatellite*; there is no strong evidence nulls
  are a practical `F_IS` problem in *sequencing-based* STR genotyping.)
- **Mis-assignment / outliers** (§4.2/§6): *Resolved* — **two complementary terms** in
  the v1 genotype likelihood: a **uniform outlier `λ`** for *random* junk (junk =
  uniform over the locus's `D` distinct observed sequences, `1/D`; `λ` a conservative
  fixed small value, simulator-calibrated, estimable-with-a-cap later), and an
  **allele-balance / overdispersion term** for *systematic concentrated* junk — the
  defence against the i.i.d.-likelihood **depth-driven het inflation** (false-het GQ
  climbing with coverage, the known SNP-path failure mode). v1 uses the lightweight
  allele-balance penalty (binomial-tail on the **deconvolved** per-allele responsibilities,
  mirroring SNP [`qual_refine.rs`](../../../src/vcf/qual_refine.rs); owed to SNP GQ too),
  feeding both GQ and site QUAL; **Dirichlet-multinomial overdispersion** (dataset-wide
  concentration, pre-pass-estimated like `ε`) is the documented upgrade. Open: the
  allele-balance penalty's exact form / strength and its simulator calibration.
  **Contamination** deferred (reuse `estimate-contamination`). Close paralogs remain an
  upstream Stage-0 curation problem.
- **Output VCF** (§4.5): *Resolved* — GangSTR/TRtools shell reused; **emit iff variable**
  (≥2 alleles; drop monomorphic; "polymorphic" MAF ≥ 0.05 is an annotation, not a gate);
  **site QUAL = Phred(locus variable)** = `−10·log₁₀ P(monomorphic)` (more than HipSTR's
  `QUAL="."`, justified by cohort-joint + the emit decision); the **`P(monomorphic)`
  estimator is resolved** (m4 amend. 2026-06-22: `Σ_{i∈A_ℓ} P(fixed-for-i)`, each via the
  engine's exact-AF convolution **kernel** with allele `i` collapsed as "reference," summed —
  *not* REF-anchored; a **generalization** of `compute_qual_via_exact_af` (parametrize the
  allele-0 collapse + per-`i` `G₀` Beta), **exact for biallelic / a per-collapse-normalizer
  approximation for ≥3 alleles**, the single-joint-`Z` form a documented upgrade —
  verify-fix #7b, §4.5); FILTER reasons
  `notPeriodic`/`tooManyAlleles`/`lowDepth` (drop `segdup`); per-sample `GT`/`GQ` with
  `./.` for absent / outlier-dominated / sub-gate samples. **Still open:** the `lowDepth`
  cohort-wide FILTER threshold and whether to also emit a per-locus polymorphism flag
  (annotation) — calibration, not the QUAL math.
- **Reading layer → its own doc?** §4.1 is intent-level here; the signature-level
  module shape (reader cursor, merger, queue, batch size `K`) is owed an
  architecture/implementation note when Stage 2 is built (arch §8 roadmap item 3).
