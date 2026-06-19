# SSR cohort calling — Mark 2 (`ssr-call`, empirical candidates)

**Status:** draft, 2026-06-19, branch `ssr-cohort`. Built **from scratch, one
agreed premise at a time** (the way [ssr_ladder_model.md](../architecture/ssr_ladder_model.md)
and the other Mark-2 docs were grown). Sections firm up only once agreed; §9 is the
live agenda of open premises. This is the **statistics / model-intent** spec for
Stage 2 — *what* the cohort caller computes and *why*; the module/struct/signature
shape (the *how*) is deferred to a later architecture doc. Reading/orchestration
*intent* is settled in §4.1 (signature-level shape still deferred); §4.2 frames the
whole model the open sections detail, and §4.4 settles how the shared parameters are
estimated (a pre-pass freezes `ε` and seeds a per-cell `θ` prior refined per locus;
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
   (§4.4) — the `|distinct seqs| × |candidates| × |slips|` alignment scores are
   **computed once for the whole run**, never recomputed; each EM iteration only
   **re-weights** them by the current stutter kernel `S_θ(Δ)` (`θ` is refined per locus,
   §4.4 — a handful of numbers, not a recompute).

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
  `n_low_quality`, `n_border_off_end`. `n_usable` = Σ observed counts (derived).
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
| **§5.4 EM topology** — E-step responsibilities, M-step for `π` (Dirichlet-smoothed); confident-homozygote seed; identifiability (shared kernel + population recurrence); convergence (non-decreasing penalised log-lik) | **reused for the per-locus loop, with pieces replaced/added.** The **E-step** (prior × likelihood → posterior) and the **IBD prior** transfer as-is; the alignments behind `Qᵣ` are constant (`ε` frozen), re-weighted each iteration by `S_θ(Δ)`. **Replaced:** the π M-step's **base measure** — the engine picks pseudocounts by SNP/indel allele *class* (REF 10 / SNP-alt 0.01 / indel-alt 0.00125, [`classify_allele`](../../../src/var_calling/posterior_engine.rs)), SNP-shaped; SSR injects the per-candidate geometric `G₀` (§5.5 row). **Added:** a **per-locus `θ_locus` M-step** (regularized stutter update, shrunk to the `θ_cell` prior — new code, not in the engine). `ε` is frozen (§4.4 pre-pass, not an M-step); `F` is the §4.4 prior-side outer loop. SNP-only machinery (compound alleles, chain anchors, contamination) is **bypassed**. The cross-locus schedule is **new code wrapping** the engine (§4.4) |
| **§5.2 stutter kernel `S_θ`** — 3-param geometric `(u,d,ρ)`, covariate-parameterised `θ(period, length, motif, purity)`, Option-1 discretized cells + shrinkage, pooling across (sample,locus) in a cell | **kernel form unchanged**; the *whole-unit step `δ`* it operates on is redefined sequence-wise (§7 / S2) |
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
`ε` and `θ` from the cohort's confident homozygotes; `F` is then re-estimated in a
prior-side outer loop *around* step 4. So step 4 itself only moves `π`.)

Step 1 and the EM core of step 4 build on §3's reuse; the reading/orchestration that
*feeds* step 1 is settled in §4.1, §4.2 frames the whole model the open sections detail,
§4.4 settles how the shared parameters are estimated (`ε` frozen; `θ` a per-cell prior
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
5. **sequencing + alignment** add **per-base** noise (miscalled base, small indel,
   slightly misplaced flank);
6. we record the sequence; over reads, the per-sample `(seq, count)` tally.

Steps 1–3 are the **prior** (`π`, `F`, ploidy); steps 4–5 are the **likelihood**
(`S_θ` ⊗ the flat per-base model `Qᵣ`). Step 6 is what Stage 1 stored.

**The EM spine — frozen `ε`, then a per-locus loop over `π` and `θ_locus` (schedule in §4.4).**
Genotypes (step 2) are the **hidden variable** — never observed, so a genotype is
*not* a parameter. The per-base error `ε` is **estimated once and frozen** before
genotyping (§4.4); the genotyping pass then runs an **independent per-locus EM** that
estimates the local allele frequencies `π` **and** the local stutter `θ_locus` (shrunk
toward its `θ_cell` prior, §4.4). It builds on the `posterior_engine.rs` §5.4 per-locus
loop — its **E-step and IBD prior reused as-is**, its **base measure swapped** for the
SSR `G₀` (§3, §5.5), with a **new per-locus `θ` M-step** added:
- **E-step** — given current `π, θ_locus` (and the frozen `ε` and current `F`), compute
  each individual's posterior weight over genotypes (`posterior ∝ prior × likelihood`,
  normalized per sample): the responsibilities;
- **M-step** — given those weights, re-estimate the local **allele frequencies `π`**
  (Dirichlet-smoothed by `G₀`) and the local **stutter `θ_locus`** (shrunk toward
  `θ_cell`); `θ_locus` only **re-weights** the cached alignments (`ε` fixed → no rebuild).

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
| **fixation index `F`** | excess-homozygosity knob; IBD-mixture prior `F·π_i + (1−F)·π_i^ploidy` | **per-individual `F_i`, estimated by the §4.4 prior-side outer loop** (EM-responsibility reduce over variable loci, shrunk toward the cohort mean; in the prior, **not** in the `align` cache → never rebuilds it); supplied default as the burn-in start. Hard ceiling `F_CEILING = 0.99` (no `F=1` absorbing trap) + optional lower CLI cap. Output is apparent **`F_IS`** — absorbs null-allele homozygote excess (§9) |
| **ploidy** | sets the genotype-space size | fixed, uniform |
| **stutter `S_θ`** | whole-unit slippage kernel | **per-cell prior from the §4.4 pre-pass, refined per-locus in the EM** — **one `θ_cell` prior per covariate cell** (**period × length** in v1) off the confident homozygotes (shrunk toward the period parent for sparse cells; literature default for empty), and each locus refines a `θ_locus` shrunk toward `θ_cell` (depth-driven partial pooling — Mark-1's relief valve). Enters `align` only via `S_θ(Δ)`, so refinement is a cheap **re-weight** (no rebuild). Motif & purity covariates deferred. The one place loci aren't independent, and what makes `θ` learnable |
| **per-base error `ε`** | flat (no per-base quals) noise inside `align` (§6) | **estimated once and frozen by the §4.4 pre-pass** — a single global value off the confident homozygotes, identifiable alongside `θ` (`θ` moves whole units, `ε` moves bases); seeded from the Stage-1 quality gate as the burn-in start. Frozen, so the `align` cache is built once and never rebuilt → genotyping is **deterministic across/within thread counts**. No shrinkage target (one global scalar); a wide pre-pass spread flags flat-`ε` is strained (→ covariate, deferred) |
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

### 4.3 Seeding the EM — π⁰, θ⁰, and the prior on π  *(settled 2026-06-19; `G₀` shape open — §9)*

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
confident homozygote in the same covariate cell** (period × length bin, whole cohort
+ genome) before fitting `(u, d, ρ)`. The **same pass** that builds the π⁰ tally
produces these skirts, so π⁰ and θ⁰ come out of one sweep. *Fallback:* a cell with
no confident homozygotes starts θ⁰ from a generic literature stutter rate for its
period / length; a sparse or high-variance cell is **shrunk** toward its period-level
parent (§4.4). The pre-pass `θ` is only a **per-cell prior** — each locus's EM refines a
`θ_locus` shrunk toward it (§4.4) — so seed contamination (a masquerading het below) is
**low-stakes**: the EM can overrule a biased prior. The confident-homozygote selection
therefore only needs to be *reasonably* clean (depth gate; flag a one-sided skirt excess
— the masquerader signature); it is prior-hygiene, not the final word (Issue 5 / §9).

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
per covariate cell — a well-defined *parametric* prior, deliberately **not** a
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

**Combine:** `π⁰ = (putative-genotype counts + G₀ pseudocounts), normalized.` With
**no confident samples** at a locus (all under-resolved / very thin) it falls back to
the **prior alone** — π⁰ is the normalized `G₀` pseudocounts (centred on the modal
candidate, which with no data *is* the seeded reference allele) — which is also the
correct conservative behaviour at small N (§5.6). **`F⁰`** is a separate seed (a supplied / global default that is also
the burn-in start of the §4.4 prior-side `F` loop that estimates `F`). *(Open — §9:
the `G₀` **shape** — how fast the pseudo-mass decays from the reference, geometric (v1)
vs Gaussian, §5.5.)*

### 4.4 Estimating the shared parameters — pre-pass (freeze `ε`, prior for `θ`) + a prior-side `F` loop  *(settled 2026-06-19)*

The chemistry parameters are **characteristic of the whole dataset**, not of any one
locus, so they are estimated from many loci in a cheap **estimation pre-pass** before
genotyping. But they are **not all treated the same way** — they differ in how
locus-specific they really are and in what they cost to vary:

- **`ε` (per-base error) is frozen.** It is a global property of the chemistry — no
  biological reason it varies locus to locus — so a single value pooled over the whole
  pre-pass beats anything one locus could estimate from its handful of mismatches; and
  because `ε` lives *inside* `align` (changing it **rebuilds** the alignment cache),
  freezing it is also what keeps the alignments computable once. This is the standard
  **two-stage / plug-in** pattern (and what **HipSTR** does with its stutter model).
- **`θ` (stutter) is a prior, not a frozen value.** Stutter genuinely varies (period ×
  length, purity), so the pre-pass produces a **per-cell prior `θ_cell`** and each
  locus's EM **refines a `θ_locus` shrunk toward that prior** (depth-driven partial
  pooling: a deep locus moves off the cell value, a shallow one stays at it). This is
  cheap — `θ` enters only through `S_θ(Δ)`, which **re-weights** the cached alignments
  rather than rebuilding them — and it is what restores Mark-1's per-locus *relief
  valve* (impure / high-depth loci that disagree with their cell), the full
  allele-vs-stutter identifiability (cohort recurrence at the locus breaks the tie),
  and, critically, makes the pre-pass `θ` **just a prior the EM can overrule** — so seed
  contamination (a masquerading het, §4.3) is **low-stakes**, not permanent (Issue 5).

`F` is the one parameter that cannot be fixed up front (it needs genotypes) and is
handled as a cheap outer loop (below).

**Why this split is still deterministic.** With `ε` frozen the expensive alignments are
computed **once for the whole cohort** (§6) and never rebuilt; `θ_locus` only re-weights
them (a handful of numbers per EM iteration). Determinism holds because everything that
varies is **per-locus, not per-thread**: `θ_locus` is a pure function of that locus's
data + the `θ_cell` prior, computed identically on whichever worker owns the locus.
(This is the crucial difference from the per-thread `ε` we rejected — *per-locus* is
atomic and reproducible; *per-thread* was order-dependent.) So Stage 2 is byte-identical
across and within thread counts, carrying the Stage-1 determinism invariant forward.

**The pre-pass, in two phases.** Both `ε` and the `θ_cell` priors come off the
**confident homozygotes** — a known-`(A,A)` sample makes every non-`A` read a labelled
stutter/error observation (§4.3): the skirt gives `θ = (u,d,ρ)`, the within-tract base
mismatches give `ε`. They are jointly identifiable because they move in different
coordinates (whole units vs bases).

- **Phase 1 — burn-in (settle).** Stream confident-homozygote loci, accumulating
  sufficient statistics into a running estimate until it meets a (generous) settle
  criterion; usually quick. If a parameter **never settles** that is a diagnostic, not a
  single verdict: too little data, or genuine structure the covariates miss (e.g.
  purity), or a bug — read it as *which assumption broke*, not *the method is wrong*.
- **Phase 2 — measure.** From the settled value, estimate the parameter **per locus**
  over a **representative, stratified** sample (spanning the length range, not just the
  first/cleanest), and build the distribution of those per-locus values.
  - **`ε`** is **frozen at the mean**, licensed by a **sensitivity** check (not a raw
    variance threshold): if perturbing `ε` by ±1 SD flips no genotype call, freezing is
    provably safe — the same fact that licenses caching against a constant. A wide spread
    means flat-`ε` is strained (`ε` wants a covariate — a documented upgrade; there is
    **no shrinkage target** for a single global scalar).
  - **`θ`** becomes the **per-cell prior** for the per-locus refinement. A sparse or
    genuinely-heterogeneous cell is **shrunk** toward its period-level parent (so the
    prior itself partial-pools: period → cell), and an empty cell takes the literature
    default. Distinguish **sampling noise** (few observations, consistent estimates →
    shrink) from **real heterogeneity** (ample data, still-wide estimates → the
    covariate is too coarse — the per-locus refinement then carries the load, which is
    exactly its job).

Because the pre-pass `θ` is only a prior, its quality matters for *convergence speed*
in low-depth loci (where the prior dominates) but not for *correctness* — so the
robust-selection hygiene of §4.3 (depth gate; flag a one-sided skirt excess, the
masquerading-het signature) is **prior-cleaning, no longer load-bearing**.

Phases 1–2 are **parallel and deterministic**: workers accumulate into per-thread
partials reduced in a fixed (catalog) order — a **commutative reduce**, order-independent.
(A locked *running overwrite* would re-introduce order-dependence; only the
*accumulation* is shared.) The per-cell `θ` structure is respected throughout — one
prior per (period × length) cell, not one global value.

**`F` — the one prior-side outer loop (per-individual `F_i`).** `F` is estimated from
the cohort-wide heterozygote deficit, which needs genotypes, so it cannot precede
genotyping. But `F` lives in the **genotype prior**, not in `align`, so re-estimating it
is cheap and **never rebuilds the cache**: run a genotyping sweep with the current `F`
(from a supplied/default `F⁰`), reduce, re-sweep. This nesting is itself an EM (`F` is a
mixing weight, the per-locus `π` loops are its inner E/M) — monotone coordinate ascent;
stop on **`|ΔF| < tol`** with a **max-rounds cap**. It is the *only* parameter that
interleaves with genotyping, and the cheap one — the engine already takes a **per-sample
`F` vector** as input, so this needs no prior-side change.

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
shrinking harder). Per-individual is also *more robust*: a contaminated or null-heavy
sample inflates only its own `F_i`, not the whole cohort's.

*The `F = 1` absorbing trap → a hard ceiling.* `F_i = 1` makes **every** heterozygous
genotype a-priori impossible for that individual — its hets can never be called however
strong the reads, and the M-step pins it there (the exact mirror of the `π = 0` trap the
§4.3 pseudocounts guard against). So `F_i` carries an **always-on hard ceiling
`F_CEILING = 0.99`** — a 100× prior down-weight on hets (overcome-able by strong
evidence), comfortably above genuine near-complete selfers (≈ 99.5 % selfing at
equilibrium, `s = 2F/(1+F)`). An **optional CLI cap** can lower it further for a more
conservative het-calling floor; the order is: raw `F_i` → shrink to the cohort mean →
clamp to the user cap if given → clamp to `F_CEILING`.

*Caveat — it is `F_IS`, not pure inbreeding.* The estimate **absorbs the homozygote
excess from null alleles** (dropout / non-amplification) and residual Wahlund structure,
which mimic inbreeding — so the output is **apparent `F_IS`**. Worse, the loop is
**self-reinforcing**: null-driven apparent homozygotes raise `F_i`, the prior then
suppresses hets, raising `F_i` further. The variable-only reduce and the `0.99`
ceiling are the v1 brakes; per-locus null-allele down-weighting is deferred (§9).

**The whole Stage-2 schedule, then:**
```
 pre-pass (confident homozygotes; parallel commutative reduce, deterministic)
   → settle (burn-in) → measure → FREEZE ε  +  build per-cell θ prior
       (shrink sparse θ-cells; flag ε-spread; literature default for empty cells)
 genotyping outer loop, F_i⁰ = supplied/default for every sample:
   repeat { per-locus EM: refine π AND θ_locus (frozen ε; θ_locus shrunk to θ_cell;
              align computed once, θ_locus only re-weights it) with current F_i
            → per-individual F_i reduce (variable loci; mean IBD-responsibility;
              shrink to cohort mean; clamp ≤ F_CEILING=0.99) }  until |ΔF| < tol / max rounds
 → genotype calls = final per-locus E-step
```
Everything heavy (the alignments) is computed **once** (because `ε` is frozen); what
iterates is cheap — `θ_locus` re-weighting inside each locus's EM, and the prior-side
`F_i` reduce across loci. Identifiability holds: allele-vs-stutter is decided *within a
locus, across samples* (the per-locus EM pools the cohort's reads at that locus, and
`θ_locus` is anchored by its cell prior), so the per-locus `θ` refinement cannot run off
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
segregate**. It is the confidence in the emit decision, and the FP defenses feed it
directly: a stutter- or artifact-induced false second allele raises `P(monomorphic)` →
**low QUAL** (the §6 allele-balance term and the stutter kernel both contribute).
Deliberately **not** about polymorphism: a variable-but-rare *real* locus scores **high**.
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
  HipSTR, the **placement is marginalized inside `align`** — the aligner may
  insert/delete whole motif units anywhere in the tract and **sums over the paths**,
  so the observed sequence itself selects the consistent placement. For a *pure*
  candidate the placement is degenerate (every placement yields the same sequence)
  and this collapses to a single `candidate ⊕ Δ`.
- **`S_θ(Δ)`** is the **stutter kernel** (§5.2 reuse): `θ_locus`, refined per locus
  under its per-cell prior (§4.4), the **same for pure and impure alleles** at a locus
  (slippage is a property of the tract, not the allele — as in HipSTR).
- **`align(obs_seq | seq)`** is a **banded pair-HMM forward** with a **flat
  (uniform-quality) emission** — a *probabilistic* alignment that **sums over all the
  ways** the observed sequence could arise from `seq` under a per-base error rate
  `ε`, yielding a genuine probability `P(obs | seq)` (not a best-path score). It
  reuses the **Stage-1 SSR pair-HMM** machinery (banded, scratch-buffered; pattern
  from BAQ `probaln`), with flat emission because Mark-2 dropped base qualities (the
  Stage-1 gate made survivors uniform-quality). The forward sum is also what
  **marginalizes the slip placement** above. An **exact-match fast path**
  (`obs == seq` byte-for-byte — the clean post-gate majority) returns `(1−ε)^len`
  without running the HMM, so the HMM fires only on the error-bearing / impure
  minority. *(Affine-gap best-path was considered and rejected: it gives a score, not
  a probability, and commits to one alignment instead of marginalizing —
  inconsistent with the `Σ_Δ` sum and the placement decision; its only edge, speed,
  is moot since `align` is cached.)*

**Provenance — HipSTR-*informed*, not HipSTR-identical.** Genuinely from HipSTR: the
**sum-over-slips** form, the **slip-placement marginalization inside `align`**, and the
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
**`align(obs | seq)` values are computed once** per `(distinct seq, candidate ⊕ Δ)` —
distinct sequences collapsing across the whole cohort — and **never recomputed**; the
stutter kernel `θ` (only through `S_θ(Δ)`) is refined per locus, so each EM iteration
merely **re-weights** those cached alignments by the updated `S_θ(Δ)` (a handful of
numbers). So `P(read | a) = Σ_Δ S_θ(Δ)·align(obs | a ⊕ Δ)` is a cheap re-weighted sum,
not a recompute, and there is **no per-thread cache and no `δ`-gated rebuild** — every
worker aligns against the same frozen `ε`. The factorization itself is the law of total
probability over the latent slip size — HipSTR's own generative model — and caching the
`ε`-dependent conditionals that don't change while `θ`/`π` move is textbook EM practice.

**`ε` is estimated once and frozen (§4.4), not an in-loop EM parameter.** We
**estimate** `ε` from the data rather than hard-coding it (a wrong constant biases
every alignment), but we do it in the §4.4 pre-pass — alongside `θ`, on the confident
homozygotes — and then hold it constant through genotyping. It is identifiable
alongside `θ` because the two live in **different coordinates**: `θ` only moves things
by **whole motif units** (`Δ`), `ε`/`align` only by **sub-motif, base-level** changes,
so neither can absorb the other's signal. `ε⁰` is seeded from the Stage-1 quality gate
(`ε = 10^(−Q/10)` at the gate's quality floor; `probaln` defaults for the gap params)
as the pre-pass **burn-in start**, then settled; the simulator calibrates the seed.
Because `ε` is constant during genotyping (and `θ_locus` enters only as a cheap
re-weight, not a rebuild), the `align` cache is built **once for the whole run** and
**never rebuilt**, and the calls are **deterministic across and within thread counts** —
the Stage-1 byte-identity invariant carries into Stage 2 (no per-thread state, no
`δ`-gated rebuild; the genotyping pass copies the frozen `ε` into each worker, and
`θ_locus` is a deterministic per-locus quantity). `ε` is a single global scalar with
**no shrinkage target**: a wide per-locus spread in the pre-pass is not something
shrinkage can fix — it signals that the flat-`ε` assumption is strained (the cue to give
`ε` a covariate, deferred). *(This supersedes the earlier per-thread-`ε` + `δ`-rebuild
design, which traded determinism for a synchronization the pre-pass makes unnecessary;
`θ`, by contrast, is **not** frozen — it refines per locus, §4.4.)*

**Cost / risk (measure, don't assume).** Two things to benchmark: (1) the cached
`align` table is 3-D per locus (`obs × candidate × slip`) and can be sizable on a
hypervariable locus — its compute-for-memory trade against the cohort-scaling thesis;
(2) the §4.4 pre-pass cost (how many loci the burn-in needs to settle `ε` and the
`θ`-cell priors, and the phase-2 measurement window) against the whole-cohort genotyping
pass it amortizes over; plus the per-locus `θ` M-step's added inner-loop iterations.
The stutter / base-error **factorization** is an approximation (the same one HipSTR
makes), made safe by the unit-of-change separation above.

**Why HipSTR cannot cache it cheaply — and why we can.** HipSTR keeps **per-base
qualities** and works on **raw reads**, so its alignment emission is *read-specific*
(a low-quality base scores differently) — there is no flat `align(seq | seq)` to
cache, and no distinct-sequence collapse to amortize one over; it also folds slip +
placement + base-error into a single per-read HMM forward pass, efficient for *its*
cost. Our two Mark-2 decisions — **drop base qualities** (uniform quality, §2) and
**collapse reads to distinct sequences** in Stage 1 — are exactly what turn `align`
into a flat, reusable function **computed once and shared across the whole cohort**;
with `ε` frozen (§4.4) the `align` values are constant through genotyping (only `S_θ`
re-weights them, as `θ_locus` refines — cheap), so the expensive part is never
recomputed per iteration. We **buy the cacheability with fidelity**: HipSTR's
quality-weighted per-read HMM is more accurate per read; our flat-error collapsed
model is coarser but precomputable once — the same fidelity-for-scale trade the rest
of the pipeline makes.

*(Resolved: `align` = flat-emission banded pair-HMM forward, slip placement
marginalized inside it; `ε` **frozen** by the §4.4 pre-pass, `θ` a **per-cell prior
refined per locus** (above). Open — §9-S3: the gate-derived `ε⁰` burn-in value and simulator calibration;
the pre-pass settle criterion and measurement window; banding / length-pruning of the
candidate × obs × slip product.)*

## 7. S2 — stutter reachability between sequences — per-allele  *(shape settled 2026-06-19; details open — §9)*

Stutter is a **per-allele** operation, not membership of a global length lattice.
The reachability relation is:

> **`B = A ⊕ k`** — `B` is `A` after `k` whole motif units are added/removed
> **within `A`'s own repeat context, the interruption structure preserved.**

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
**marginalized inside the §6 alignment** (HipSTR), so the data picks it. *(Open —
§9-S2: whether reachability needs the reference frame as context; the impure-allele
penalty's form.)*

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
putative genotype → π⁰ tally, θ⁰ from confident homozygotes, `G₀` pseudocount floor),
and the **shared-parameter estimation architecture** (§4.4: a two-phase pre-pass
freezes `ε` and seeds a per-cell `θ` prior refined per locus; `F` is a prior-side outer
loop; the whole pass is deterministic), and the **output** (§4.5: emit variable loci,
site QUAL = Phred(locus variable), SSR FILTER reasons).
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
- *Resolved:* slip-site placement for an impure sequence is **marginalized inside the
  §6 alignment** (HipSTR), not committed by S2 or pre-enumerated as candidates.
- *Resolved:* **no impurity penalty** — impure alleles treated identically to pure;
  the empirical candidate set + HMM + recurrence already filter spurious alleles. Add
  a penalty only if recurrent systematic impure artifacts are *measured* to leak.
- Whether reachability needs the **reference frame** as context, or is a pure
  per-allele sequence operation.

**S3 — HipSTR likelihood (§6).**
- *Resolved:* `align(·)` = **flat-emission banded pair-HMM forward** (reusing the
  Stage-1 SSR pair-HMM), with an exact-match fast path; affine-gap rejected.
- *Resolved:* `ε` is **estimated once and frozen** by the §4.4 pre-pass
  (confident-homozygote phase + sensitivity gate), held constant through genotyping — so
  the `align` cache is built once, **never rebuilt**, and the calls are **deterministic
  across and within thread counts**. `θ` is **not** frozen: the pre-pass gives a per-cell
  prior and each locus refines a `θ_locus` shrunk toward it (§4.4) — cheap because `θ`
  only re-weights the cache, deterministic because it is per-locus (not per-thread), and
  it makes seed contamination low-stakes (Issue 5). Identifiable because `θ` moves whole
  units and `ε` moves bases. (Supersedes the earlier per-thread-`ε` + `δ`-rebuild
  design.) Open: the gate-derived `ε⁰` burn-in value, the settle criterion, the phase-2
  measurement window, and the `θ_locus`←`θ_cell` shrinkage strength.
- Cost control (measure): the 3-D `align` cache's memory, and how often `ε`
  re-estimation rebuilds it; banding / length-pruning of the `candidate × obs × slip`
  product.

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
  covariate cell (parametric, not non-parametric); Gaussian is the documented upgrade.
  It seeds π⁰ as well as regularizes. Open: mode vs mean as the centre (mode chosen —
  more robust); multimodal cohorts left to the data, not the prior.
- Working-set / memory: Stage 2 holds, per locus, observed sequences + the **cached
  3-D `align` table** (`obs × candidate × slip`, §6); the cohort working set is
  bounded by §4.1's lockstep rule (≈ N × one block) — **benchmark** the cache's
  compute-for-memory trade (esp. hypervariable loci) against the cohort-scaling thesis
  ([ssr_ladder_model.md](../architecture/ssr_ladder_model.md) §5); the §6 speed win is
  itself **unmeasured**.
- **Stutter covariate cells** (§4.2 / §4.4): *Resolved* — **one `θ_cell` prior per
  period × length cell** in v1, built by the §4.4 pre-pass (confident-homozygote counts
  accumulate into each cell via a deterministic commutative reduce), then **each locus
  refines a `θ_locus` shrunk toward its `θ_cell` prior** in the per-locus EM (depth-driven
  partial pooling — Mark-1's relief valve). `θ` enters `align` only via `S_θ(Δ)`, so the
  per-locus refinement is a cheap **re-weight** (no rebuild), and it is **deterministic**
  (per-locus, not per-thread). **Shrinkage** also forms the prior itself (sparse / hetero-
  geneous cells shrink toward the period-level parent; empty cells take a literature
  default). So **per-locus `θ` is in v1** (the partial-pooling refinement), not deferred —
  the open part is the `θ_locus`←`θ_cell` shrinkage strength. **Motif** and **purity**
  covariates deferred — purity (interruptions stabilize a tract → less stutter) is the
  Mark-2-specific one, added once impure alleles are common enough per cell to estimate it.
- **Per-base error rate** (S3 / §4.2 / §4.4): *Resolved* — `ε` is **estimated** in the
  §4.4 pre-pass (it does **not** trade off against `θ`: `θ` moves whole units, `ε`
  moves bases) and **frozen** for genotyping. A single global value in v1; a coarse
  covariate is the documented upgrade if the pre-pass spread shows flat-`ε` is strained.
- **`F` granularity** (§4.2 / §4.4): *Resolved* — **per-individual `F_i`**, estimated
  by the §4.4 **prior-side outer loop** (mean-IBD-responsibility reduce over variable
  loci, shrunk toward the cohort mean; sits in the prior, never touches the `align`
  cache). Captures structured / mixed-mating cohorts without population labels and
  localizes a bad sample's effect; the engine already takes a per-sample `F` vector, so
  no prior-side change. Convergence `|ΔF| < tol` + max-rounds cap. **Hard ceiling
  `F_CEILING = 0.99`** (no `F=1` absorbing trap — the `F`-analog of the `π` pseudocount
  floor) + **optional lower CLI cap**. **Per-locus `F` rejected** (confounded with π).
  Caveat: the estimate is apparent **`F_IS`** — it absorbs null-allele homozygote excess
  + Wahlund in a self-reinforcing loop (v1 brakes = variable-only reduce + the
  ceiling), so it is "inbreeding + null-allele bias" until null alleles get their own
  mechanism (below).
- **Null alleles** (alleles that fail to amplify / map) cause *apparent* per-locus
  homozygote excess that **inflates the estimated `F_IS`** (§4.4) — a known SSR
  confounder. v1 does **not** separate them, so the reported `F` is apparent `F_IS`, not
  pure inbreeding (the variable-only reduce + the `0.99` ceiling are the only brakes
  on the self-reinforcing loop). A dedicated mechanism — per-locus null-allele detection
  (Micro-Checker-style hom-excess inconsistent with the cohort `F`) feeding a down-weight
  in the `F` reduce, or an explicit null-allele frequency — is **deferred** until
  measured to be necessary.
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
  `QUAL="."`, justified by cohort-joint + the emit decision); FILTER reasons
  `notPeriodic`/`tooManyAlleles`/`lowDepth` (drop `segdup`); per-sample `GT`/`GQ` with
  `./.` for absent / outlier-dominated / sub-gate samples. Open: the exact
  `P(monomorphic)` estimator (from the cohort `π` posterior), the `lowDepth` cohort-wide
  threshold, and whether to also emit a per-locus polymorphism flag.
- **Reading layer → its own doc?** §4.1 is intent-level here; the signature-level
  module shape (reader cursor, merger, queue, batch size `K`) is owed an
  architecture/implementation note when Stage 2 is built (arch §8 roadmap item 3).
