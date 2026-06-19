# SSR cohort calling — Mark 2 (`ssr-call`, empirical candidates)

**Status:** draft, 2026-06-19, branch `ssr-cohort`. Built **from scratch, one
agreed premise at a time** (the way [ssr_ladder_model.md](../architecture/ssr_ladder_model.md)
and the other Mark-2 docs were grown). Sections firm up only once agreed; §9 is the
live agenda of open premises. This is the **statistics / model-intent** spec for
Stage 2 — *what* the cohort caller computes and *why*; the module/struct/signature
shape (the *how*) is deferred to a later architecture doc. Reading/orchestration
*intent* is settled in §4.1 (signature-level shape still deferred); §4.2 frames the
whole model the open sections detail. Both firmed 2026-06-19.

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
   **HipSTR-style, flat-quality** likelihood that sums over the candidate's PCR slips
   (§6 / S3). Identical observed sequences collapse across the whole cohort, so the
   `|distinct seqs| × |candidates| × |slips|` alignment scores are computed once
   (θ-independent), reused across EM iterations.

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
between alleles," all of which Mark-2 still supplies. They are **inherited, not
redesigned**; this doc restates only the Mark-2-specific *adapter*, not the math.

| from `ssr_genotyping.md` | Mark-2 status |
|---|---|
| **§5.3 genotype prior** — inbreeding-adjusted IBD-mixture (`F·π_i + (1−F)·π_i^ploidy`), reuse `src/var_calling/posterior_engine.rs` verbatim | **unchanged** — feed it repeat-allele *sequences* instead of bases; adapter only |
| **§5.4 EM topology** — E-step responsibilities, M-step for `π` (Dirichlet-smoothed), M-step for the shared stutter kernel, optional M-step for `F`; confident-homozygote seed; identifiability (shared kernel + population recurrence); convergence (non-decreasing penalised log-lik) | **unchanged in topology** — the per-read term it consumes is the Mark-2 `Qᵣ⊗S_θ` (§6–§7) |
| **§5.2 stutter kernel `S_θ`** — 3-param geometric `(u,d,ρ)`, covariate-parameterised `θ(period, length, motif, purity)`, Option-1 discretized cells + shrinkage, pooling across (sample,locus) in a cell | **kernel form unchanged**; the *whole-unit step `δ`* it operates on is redefined sequence-wise (§7 / S2) |
| **§5.5 Dirichlet base measure `G₀`** — unimodal, reference-centred in the signed unit offset `Δ`; geometric v1 / Gaussian upgrade | **adapted** — `Δ` is still computable (ref tract length is the frame); now **purity-agnostic** (the Mark-1 off-ladder `OFFLADDER_PRIOR_FACTOR` is **dropped** — §4.3/§7) |
| **§5.6 small-N / single-sample**, **§5.7 ploidy & `F`**, **§5.8 FP-aversion gates** (posterior / depth / support / stutter-fraction thresholds) | **unchanged** |
| **§5.9 output VCF** — GangSTR/TRtools-compatible, REF/ALT = actual tract sequences, derived `REPCN`/`BPDIFFS` | **unchanged** — Mark-2 *already* stores sequences, so REF/ALT fall out directly |

The generative model (§5.1: `π` → genotype → read via allele-pick × stutter ×
sequencing-error) is **identical**; only the three leaf definitions below change.

---

## 4. The Mark-2 Stage-2 pipeline (overview)

```
 N × sample.ssr.psp ─┐                          per locus ℓ:
 (seq, count)        │   coordinate-merge   1. POOL    per-sample dists → cohort-aggregate seq dist
 .ssr.catalog ───────┼─►  scan (spec §2)    2. ASSEMBLE candidate set A_ℓ from pooled seqs   ← S1 (§5)
 (motif, ref frame)  │                       3. SCORE   Qᵣ(obs_seq | cand) flat-error, once  ← S3 (§6)
                     │                       4. EM      π, θ over A_ℓ; reachability δ seq-wise ← S2 (§7)
                     └─►                      5. CALL    genotype + posterior + VCF (§5.9 reuse)  ─► cohort.vcf
```

Steps 1, 4-call are inherited (§3); the reading/orchestration that *feeds* step 1
is settled in §4.1, and §4.2 frames the whole model the open sections detail.
Steps 2, 3, and the `δ` inside 4 are the new work this spec must settle (§5-§7,
agenda §9).

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

**The EM spine (reused §5.4 topology).** Genotypes (step 2) are the **hidden
variable** — never observed, so a genotype is *not* a parameter. EM alternates:
- **E-step** — given current `π, F, θ`, compute each individual's posterior weight
  over genotypes (`posterior ∝ prior × likelihood`, normalized per sample): the
  responsibilities;
- **M-step** — given those weights, re-estimate the **population** parameters `π`,
  `F`, and the stutter `θ`.

The per-sample genotype calls we report are a read-off of the **final** E-step.

**The ingredients, and their status:**

| ingredient | role | status |
|---|---|---|
| **candidate set `A_ℓ`** | the support of `π` and the genotype space | **assembled upstream of EM** (S1, §5, OPEN) — gates everything: too many ⇒ stutter blobs become "alleles"; too few ⇒ missed variation |
| **allele frequencies `π`** | population prior over alleles | estimated (M-step) |
| **fixation index `F`** | excess-homozygosity knob; IBD-mixture prior `F·π_i + (1−F)·π_i^ploidy` | **global, estimated in the EM** (v1; settles like `ε`, but **not** in the `align` cache → no `δ`-rebuild, and a cheap *deterministic* global reduce); supplied default as override; **per-individual `F_i`** the flagged extension (wild-tomato cohorts — `F` varies by population) |
| **ploidy** | sets the genotype-space size | fixed, uniform |
| **stutter `S_θ`** | whole-unit slippage kernel | estimated (M-step); **one `θ` per covariate cell** (**period × length** in v1), pooled across the cell's loci+samples like `ε`/`F` (cheap per-cell reduce, **not** in the `align` cache); bootstrapped by `θ⁰` (§4.3). **No per-locus `θ` / no shrinkage in v1**; motif & purity covariates deferred. The one place loci aren't independent, and what makes `θ` learnable |
| **per-base error `ε`** | flat (no per-base quals) noise inside `align` (§6) | a **per-thread EM parameter** (each worker carries a running `ε`, warm-started locus-to-locus), estimated alongside `θ` (identifiable: `θ` moves whole units, `ε` moves bases); seeded from the Stage-1 gate; value keeps updating but the `align`-cache rebuild is gated by a tolerance `δ` |
| **base measure `G₀`** | regularizing prior **on `π` itself** | reused §5.5; reference-centred, smooth — where most **FP control** + off-ladder down-weight lives |
| **seed + convergence** | EM start + stop | reused §5.4 — confident-homozygote seed; iterate to non-decreasing penalised log-lik |

**Why the deconvolution is well-posed (identifiability).** With `π, F, θ` free over
a candidate set, a stutter peak could be explained *either* as slippage of a common
allele *or* as a real rare allele. The tie is broken by the **population**: a true
allele recurs across unrelated individuals, whereas a stutter shadow's abundance
tracks its parent's at the *learned* slippage rate. So cohort-wide pooling (and the
shared `θ`) is what makes the problem identifiable, not merely more powerful — hence
the small-N caveat (reused §5.6).

**Three distinct noise sources — keep them separate.** (1) **stutter** —
whole-repeat-unit length change (PCR, dominant, length-dependent); (2) **per-base
error** — substitutions / small indels / flank misplacement (the flat `Qᵣ`); (3)
**mis-assignment** — a read from a paralog/duplication or cross-sample
contamination. Stage-0/Stage-1 gates remove most of it, but a **uniform outlier
component `λ`** is included in the genotype likelihood (§6) to absorb the residual:
without it the EM must explain every oddball read as evidence for *some* allele, and
the cheapest explanation is an extra allele → **inflated false heterozygotes**.
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
period / length; the M-step refines `θ` regardless of the start.

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
offset from the reference (more mass near the reference; the geometric's
heavier-than-Gaussian tail keeps real large-step alleles callable), **reference-centred**
and **symmetric** in v1, with its **decay parameter fit from pooled empirical data**
per covariate cell — a well-defined *parametric* prior, deliberately **not** a
data-hungry non-parametric one (pinning a shape needs many observations, so we fix the
shape as geometric and estimate its parameter, exactly as for the stutter kernel). It
is purely a **length-offset** prior — **agnostic to purity**: an impure allele of
length L gets the same prior as a pure allele of length L (impure alleles are *not*
specially penalized — §7). This is the §5.5 Dirichlet base measure; it is a
**false-positive control** (alongside the HMM likelihood and population recurrence),
and the **same** pseudocounts re-enter every M-step, so it regularizes throughout, not
only at the start.

**Combine:** `π⁰ = (putative-genotype counts + G₀ pseudocounts), normalized.` With
**no confident samples** at a locus (all under-resolved / very thin) it falls back to
the **prior alone** — π⁰ is the normalized `G₀` pseudocounts (reference-centred, the
reference allele as anchor) — which is also the correct conservative behaviour at
small N (§5.6). **`F⁰`** is a separate seed (a supplied / global default, estimated
later in the optional `F` M-step, §5.3). *(Open — §9: the `G₀` **shape** — how fast
the pseudo-mass decays from the reference, geometric (v1) vs Gaussian, §5.5.)*

---

## 5. S1 — candidate assembly from observed sequences  *(shape settled 2026-06-19; thresholds open — §9)*

Evidence is assembled in three levels — **pool → rungs → candidates** — then a
coarse **locus-admission filter** decides whether the locus behaves enough like an
SSR to analyze at all. Throughout, **recall is the goal and precision is the EM's
job** (§4.2): a sequence dropped here is gone, whereas a wrongly-kept candidate is
driven to `π ≈ 0` by the population recurrence + `G₀` (§5.5). Sequences that do
*not* become candidates are **not discarded** — they re-enter the likelihood as
slip/error products of the candidates (the §6 invariant).

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

## 6. S3 — the read likelihood `Qᵣ(obs_seq | candidate)` — HipSTR-style  *(shape settled 2026-06-19; details open — §9)*

We adopt the **HipSTR generative form**, adapted to Mark-2's collapsed-sequence
evidence (distinct observed sequences + counts, uniform quality — no per-read rows,
no base qualities). The probability of an observed sequence given a candidate allele
**sums over the PCR slips** the allele could have undergone:

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
- **`S_θ(Δ)`** is the **shared stutter kernel** (§5.2 reuse): one `θ` per covariate
  cell, the **same for pure and impure alleles** at a locus (slippage is a property
  of the tract, not the allele — as in HipSTR).
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

**The outlier component (`λ`).** `Qᵣ` scores *one* observed sequence against *one*
candidate; the **genotype** likelihood combines it over the genotype `G`'s alleles
**plus a uniform outlier term**: `P(read | G) = (1−λ)·[mix over G's alleles] +
λ·junk`. Without it, every oddball read — and these repetitive, hard-to-assemble
tracts will throw some — must be explained as evidence for *some* allele, and the
cheapest explanation is to posit an **extra allele**, inflating **false
heterozygotes**. The `λ` term is the escape valve: junk goes there instead of minting
an allele, so the genotype is decided by the well-explained reads. `λ` is a
conservative **fixed small value in v1** (simulator-calibrated); estimable later **but
capped** — too high dumps real minor alleles (het *deflation*), too low lets het
inflation return. It **complements** the candidate set's recurrence filter (which
removes junk at the *cohort* level) by absorbing residual junk at the *per-read*
level. This is v1, not deferred — it is structural FP protection, not an optimization.

**Why the per-iteration cost is small.** `θ` enters `Qᵣ` only through `S_θ(Δ)`, so
each inner EM iteration merely **re-weights** the cached `|distinct seq| ×
|candidates| × |slips|` alignment table by the updated `S_θ(Δ)` (a handful of
numbers) — the expensive alignments are not recomputed; distinct sequences collapse
across the whole cohort. `ε` (the only other thing `align` depends on) is carried
**per worker thread** — each thread refines a running `ε` across the loci it
processes (it pools thousands of loci, so it is well-determined and near-identical
across threads), **not** re-estimated per inner iteration. A change in a thread's `ε`
rebuilds that thread's cached `align` tables (it is inside `align`), so rebuilds are
gated by a tolerance `δ` (below) and — since `ε` settles fast — should be rare. The
factorization itself is just the law of total probability over the latent slip size —
HipSTR's own generative model — and caching the conditional likelihoods that don't
change per iteration is textbook EM practice.

**`ε` is a full EM parameter, with a lazy cache-rebuild rule (`δ`).** We **estimate**
`ε` in the model rather than hard-coding it (a wrong constant biases every alignment).
It is identifiable alongside `θ` because the two live in **different coordinates**:
`θ` only moves things by **whole motif units** (`Δ`), `ε`/`align` only by **sub-motif,
base-level** changes — neither can absorb the other's signal. `ε`'s *value* keeps
updating each round (cheap), but the expensive **`align`-cache rebuild is gated by a
tolerance `δ`**: a thread rebuilds its cache only when `ε` has drifted **≥ `δ`** since
the cache was last built; below `δ` it treats `ε` as **settled** and reuses the cache
(alignments change negligibly for a sub-`δ` move, so the staleness is bounded by
`δ`). This is a **deliberate, acknowledged premature optimization** — recorded now
*because it is easy to forget and likely to matter*; `δ` is a calibration knob (too
tight → needless rebuilds, too loose → stale alignments). `ε⁰` is seeded from the
Stage-1 quality gate (`ε = 10^(−Q/10)`
at the gate's quality floor; `probaln` defaults for the gap params), then estimated;
the simulator calibrates it. The gate value is **only a bootstrap** for round 0 —
there is **no fixed prior** on `ε`. Each **worker thread** carries its own running
estimate, refined **locus-by-locus** as it drains its queue (the previous locus's
`ε` warm-starts the next — the running-prior idea, well-defined *within* a thread's
sequential order); the bootstrap's influence vanishes once base counts accumulate.
Each thread pools thousands of loci, so the per-thread estimates are well-determined
and converge to nearly the same value; keeping `ε` per-thread **avoids any global
synchronization**. The price: `ε` (and the rare borderline call) is **not
byte-identical across — or, with the dynamic work-queue, within — thread counts**,
accepted under the Stage-2 determinism stance (fixed iteration order + pinned
convergence, not cross-thread byte-identity); a deterministic locus→thread partition
would restore fixed-thread reproducibility if measured to matter.

**Cost / risk (measure, don't assume).** Two things to benchmark: (1) the cached
`align` table is 3-D per locus (`obs × candidate × slip`) and can be sizable on a
hypervariable locus — its compute-for-memory trade against the cohort-scaling thesis;
(2) how often `ε` crosses the rebuild tolerance `δ`, and the right `δ` (tight enough
that stale alignments don't bias calls, loose enough that rebuilds don't dominate).
The stutter / base-error **factorization** is an approximation (the same one HipSTR
makes), made safe by the unit-of-change separation above.

**Why HipSTR cannot cache it cheaply — and why we can.** HipSTR keeps **per-base
qualities** and works on **raw reads**, so its alignment emission is *read-specific*
(a low-quality base scores differently) — there is no flat `align(seq | seq)` to
cache, and no distinct-sequence collapse to amortize one over; it also folds slip +
placement + base-error into a single per-read HMM forward pass, efficient for *its*
cost. Our two Mark-2 decisions — **drop base qualities** (uniform quality, §2) and
**collapse reads to distinct sequences** in Stage 1 — are exactly what turn `align`
into a flat, reusable function that is **re-weighted (not recomputed) across `θ`
iterations** and shared across the whole cohort. We **buy the cacheability with
fidelity**: HipSTR's
quality-weighted per-read HMM is more accurate per read; our flat-error collapsed
model is coarser but precomputable once — the same fidelity-for-scale trade the rest
of the pipeline makes.

*(Resolved: `align` = flat-emission banded pair-HMM forward, slip placement
marginalized inside it, `ε` a global EM parameter (above). Open — §9-S3: the
gate-derived `ε⁰` value and simulator calibration; whether `ε` converges (freeze /
lazy-rebuild only if rebuilds are measured to be a real cost); banding /
length-pruning of the candidate × obs × slip product.)*

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
putative genotype → π⁰ tally, θ⁰ from confident homozygotes, `G₀` pseudocount floor).
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
- *Resolved:* `ε` is a **per-thread EM parameter** (each worker carries a running
  `ε`, warm-started locus-to-locus), estimated alongside `θ`, seeded from the Stage-1
  gate. Its *value* keeps updating; the expensive **`align`-cache rebuild is gated by
  a tolerance `δ`** (rebuild only when `ε` drifts ≥ `δ`) — a deliberate optimization
  recorded now, `δ` to calibrate. Caveat: not reproducible run-to-run under the
  dynamic queue; accepted, fixable via a deterministic locus→thread partition if
  needed.
- Cost control (measure): the 3-D `align` cache's memory, and how often `ε`
  re-estimation rebuilds it; banding / length-pruning of the `candidate × obs × slip`
  product.

**Cross-cutting.**
- Base measure `G₀` (§5.5): *Resolved (§4.3)* — a **geometric** length-offset decay
  (reference-centred, symmetric v1, **purity-agnostic** — no impurity down-weight,
  §7), with its **decay parameter fit from pooled data** per covariate cell
  (parametric, not non-parametric); Gaussian is the documented upgrade. It seeds π⁰
  as well as regularizes.
- Working-set / memory: Stage 2 holds, per locus, observed sequences + the **cached
  3-D `align` table** (`obs × candidate × slip`, §6); the cohort working set is
  bounded by §4.1's lockstep rule (≈ N × one block) — **benchmark** the cache's
  compute-for-memory trade (esp. hypervariable loci) against the cohort-scaling thesis
  ([ssr_ladder_model.md](../architecture/ssr_ladder_model.md) §5); the §6 speed win is
  itself **unmeasured**.
- **Stutter covariate cells** (§4.2): *Resolved* — **one `θ` per period × length cell**
  in v1, pooled per cell like `ε`/`F` (each locus's expected stutter counts accumulate
  into its cell; `θ` only re-weights the `align` cache so the update is cheap),
  bootstrapped by the `θ⁰` confident-homozygote seed (§4.3). **Per-locus `θ` /
  shrinkage** deferred (SSR stutter is well-predicted by period+length; add only if
  within-cell heterogeneity is measured). **Motif** and **purity** covariates deferred
  — purity (interruptions stabilize a tract → less stutter) is the Mark-2-specific one,
  added once impure alleles are common enough per cell to estimate it.
- **Per-base error rate** (S3 / §4.2): a fixed constant vs a coarse covariate — and
  kept **un-estimated** so it does not trade off against the stutter `θ`. Source of
  the constant is the open S3 question above.
- **`F` granularity** (§4.2): *Resolved* — **global, estimated** in the EM for v1
  (cheap deterministic global reduce, not in the `align` cache), supplied default as
  override. **Per-individual `F_i`** is the flagged extension (wild-tomato cohorts:
  `F` varies by population — `F_i` captures it without explicit population
  assignments; needs engine extension beyond the global-`F` reuse). **Per-locus `F`
  rejected** (confounded with π; its driver, null alleles, handled separately).
- **Null alleles** (alleles that fail to amplify / map) cause *apparent* per-locus
  homozygote excess that mimics high `F` — a known SSR confounder. **Not** absorbed
  into `F`; flagged for a separate explicit mechanism if measured to be necessary.
- **Mis-assignment / outliers** (§4.2/§6): *Resolved* — a **uniform outlier component
  `λ`** is in the v1 genotype likelihood (absorbs reads no allele explains, preventing
  false-het inflation); `λ` a conservative fixed small value to start
  (simulator-calibrated), estimable-with-a-cap later. **Contamination** deferred (reuse
  `estimate-contamination`). Close paralogs remain an upstream Stage-0 curation problem.
- **Reading layer → its own doc?** §4.1 is intent-level here; the signature-level
  module shape (reader cursor, merger, queue, batch size `K`) is owed an
  architecture/implementation note when Stage 2 is built (arch §8 roadmap item 3).
