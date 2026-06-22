# SSR Stage 2 — `ssr-call` parameter & prior estimation (architecture sketch)

**Status:** draft, 2026-06-19 (chemistry model & vocabulary **settled 2026-06-21**),
branch `ssr-cohort`. The second of three Stage-2 (`ssr-call`) sketches. This one is
**Phase 2: the estimation pre-pass that produces the *fixed parameters and the
priors* the genotyping EM needs** — it freezes the per-base error `ε`, fits the
**stutter shape** and per-sample-group **stutter level**, fits the `G₀` pseudocount
prior on `π`, and seeds `π⁰`/`θ⁰`/`F⁰`. Companions:

- [ssr_call_reading.md](ssr_call_reading.md) — Phase 1: the reader + k-way merge
  that feeds `CohortLocus` work-items (this phase is a *consumer* of that merge).
- [ssr_call_genotyping.md](ssr_call_genotyping.md) — Phase 3: candidate assembly,
  likelihood, the per-locus EM, the `F` loop, and the VCF (the *consumer* of
  everything this phase produces).

The *how-the-code-is-wired* companion to spec
[ssr_cohort_mark2.md §4.3 + §4.4](../specs/ssr_cohort_mark2.md) (seeding + shared-
parameter estimation) and §5.5 (the `G₀` base measure). Where they disagree on
intent, the spec wins; on code layout, this doc wins.

---

## 0. Vocabulary (settled 2026-06-21)

The chemistry model has a lot of moving terms; these are the canonical names. Read
the rest of the doc against this table.

| term | meaning |
|---|---|
| **cohort** | all the samples |
| **sample** | one individual |
| **sample group** | a **data-driven cluster of samples** sharing chemistry + provenance (PCR protocol, DNA preservation: fresh vs. degraded/ancient, …). **Unobserved** — there are no reliable user/SRA labels — so it is *inferred* from each sample's `ε` + stutter level. Replaces the earlier "library protocol" / "sample protocol" / "protocol group" / "soft group." It's a *group*, not a *protocol*, because samples cluster by provenance as much as by prep. |
| **chemistry** | the parameter bundle that characterizes a sample group: the per-base error `ε` and the stutter level. (Stutter *shape* is **not** chemistry — see below — it's shared.) |
| **`ε`** (per-base error) | within-tract **substitution** rate (C1: `align` admits no sub-motif indel inside the tract, so `ε` moves composition, not length). **Frozen per sample group** (lives in the `align` cache). |
| **stutter shape** | the *profile* of the stutter distribution — contraction bias, geometric decay `ρ`, up/down asymmetry, multi-step mass. **Per `(group, period)`, shrunk toward a cohort-per-period parent `θ_period`** (M3 amend.; *not* assumed cohort-shared), refined per locus into `θ_locus` (shrunk toward the sample's group-period shape). Protocol/preservation are expected to move the *rate* far more than the *profile* (human PCR evidence) — but for non-model/degraded DNA that invariance is left to **shrinkage to decide**, not assumed (§1). |
| **stutter level** | the overall stutter *rate* — the probability a spanning read slips by ≥1 unit. It rises with repeat length, so it is modeled **per sample group** as a line in repeat length: `level(length) = level_baseline + level_slope · length`. (Two numbers per sample group — the two rows below. Not a single scalar, because the rate depends on length.) **Pre-pass-seeded, then refined in the genotyping outer loop — *not* frozen** (it is an `S_θ` re-weight outside the `align` cache, so refining it rebuilds nothing; C2 amend.). |
| **level baseline** | the *intercept* of the stutter-level line — the stutter rate at the short-length end of the modeled range. **Per sample group.** |
| **level slope** | the *slope* of the stutter-level line — the increase in stutter rate per additional repeat unit (the rate's "length response"). **Per sample group** (point 2 below: per-group is the default; we have no evidence the length response is cohort-universal). |
| **period** | a locus's repeat period = its motif length in bp (1 = mono, 2 = di, …); from the catalog. |
| **repeat unit** | one copy of the locus's motif (`period` bp). Allele lengths, slips, and offsets are counted in **repeat units**, not in bp. |
| **unit offset** (`Δ`) | a signed count of **repeat units** between two allele lengths (`+` = longer, `−` = shorter). Used in two places: a candidate's `Δ` from the per-locus **modal allele** (the frame for `G₀`, below), and a read's slip `Δ` from a candidate allele (the frame for the stutter kernel `S_θ(Δ)`). E.g. a candidate two repeats longer than the mode has `Δ = +2`. |
| **loci group** | a set of loci grouped by shared behaviour (their repeat **covariates**), over which a *pooled* parameter is estimated — the loci-side analogue of a **sample group**. **In v1 the only grouping covariate is `period`, so a loci group = a period.** Repeat **length** is *not* a grouping axis — it enters as a continuous *level* covariate (`level_slope`), never as a bin; **motif** & **purity** are deferred covariates. The **`G₀`** decay is fit per loci group; the **stutter shape** is fit per **`(sample group, period)`** (shrunk to a cohort-per-period parent — M3), i.e. it crosses a loci group (period) with a sample group. |
| **geometric pseudocounts** (`G₀`) | the prior on candidate allele frequencies `π` (§5): a small prior "count" of mass placed on **every** candidate allele, **decaying geometrically** with the candidate's **unit offset** `Δ` (repeat units from the per-locus modal allele — see above). A *pseudocount* is prior count added to the observed counts (Dirichlet smoothing) so no candidate gets `π = 0`; *geometric* names the decay shape. It keeps every `π_i > 0` (no `π = 0` absorbing trap — **floored at a tiny `> 0`** so `p^|Δ|` can't underflow to exactly 0, verify-fix #4 / §5), so a candidate the seed missed (e.g. a masked het) stays recoverable, and re-enters **every** EM M-step as a small-N regularizer + false-positive control. Decay parameter fit **per loci group**. (= the spec §5.5 Dirichlet base measure.) |
| **confident genotype** (CG-seed, §2) | a (sample, locus) whose read-length distribution **resolves into clear, well-separated peaks** at full confidence — **one** peak (homozygote) or up to **ploidy** peaks (heterozygote), each peak's allele **cohort-recurrent**, peaks **≥ 2 repeat units apart** with dosage-consistent heights. These are the pre-pass's *labelled* observations: each read's slip attributes to a known parent allele, so the skirt is stutter and the within-tract mismatches are `ε`. **Replaces "confident homozygote"** as the chemistry seed (a homozygote is the 1-peak special case); a well-separated het contributes its two outer skirts as hard labels (inner valley → soft EM). A merged (1-apart) het is **not** a confident genotype (stutter and allele confounded). |

The kernel the likelihood uses is `S_θ(Δ) = level × shape(Δ)`: the **shape** says
*where* a slip lands (per `(group, period)`, shrunk to a cohort-per-period parent), the
**level** (`level_baseline + level_slope · length`) says *how often* a slip happens (per
group, growing with length).

**Granularity default — when in doubt, per sample group.** Every chemistry parameter
defaults to **sample-group** granularity; a parameter is pooled up to the **cohort**
*only* when the data show it does not vary by group. **All** of `ε`, level, **and shape**
follow this rule (M3 amend.): the **stutter shape** is fit **per `(sample group, period)`,
shrunk toward a cohort-per-period parent** — so if the PCR-free "rate not profile"
expectation holds the shrinkage collapses the groups to the shared shape for free, and if
it fails (non-model / degraded DNA / odd motifs) a group keeps its own profile. The
**stutter level** (both `level_baseline` *and* `level_slope`) likewise stays **per sample
group**. The cohort is the exception that must be *earned by the data* (via shrinkage), not
assumed up front — which is the M3 correction: shape was previously cohort-shared by
assumption on human-PCR evidence.

---

## 1. What this phase computes, and why it is a *separate pass*

The chemistry parameters are **properties of the dataset, not of any one locus**
(spec §4.4), so they are estimated once, from many loci, *before* genotyping — a
cheap pre-pass. The things genotyping (Phase 3) cannot start without:

| output | kind | spec | how it enters Phase 3 |
|---|---|---|---|
| **`ε`** (per-base error) | **frozen — per sample group** | §4.4, §6 | lives *inside* `align`; frozen ⇒ `align(obs|cand⊕Δ)` is a **pure, iteration-invariant function** ⇒ Stage 2 deterministic *whether recomputed or cached* (the cache is **deferred** — genotyping §4/Q-G3). Group granularity *because* it's in `align` (per-sample `ε` = upgrade) |
| **stutter shape** (per `(group, period)`, shrunk to cohort-per-period parent — M3) | **prior**, refined per locus | §4.4, §5.2 | the per-group slip profile; each locus's EM refines a `θ_locus` *shrunk toward* the sample's group-period shape (catches motif-identity outliers); enters `align` only via `S_θ(Δ)` → a cheap **re-weight**, not a rebuild |
| **stutter level** (per sample group: `level_baseline + level_slope·length`) | **seed → refined in outer loop** (C2) | §4.4 | the rate vs. length; the pre-pass **seeds** it (`level⁰`), the genotyping outer loop refines it per group from soft responsibilities; enters only the `S_θ` re-weight (no `align` rebuild). Samples are clustered into sample groups from their pre-pass `(ε, level)`; **both** `level_baseline` and `level_slope` are fit per group |
| **`G₀`** (geometric pseudocounts on `π`) | **prior**, per loci group | §4.3, §5.5 | floors `π⁰` (no `π=0` absorbing trap) and regularizes **every** M-step |
| **`π⁰`, `θ⁰`** | **seeds** | §4.3 | the EM's starting point (basin selection) |
| **`F⁰`** | **seed** | §4.3, §4.4 | burn-in start of the Phase-3 prior-side `F_i` loop (estimated *there*, not here) |

**Why both shape and level carry a group axis (level strongly, shape with shrinkage —
M3 amend. 2026-06-22).** The field reports that **protocol and preservation move the
stutter *rate* far more than its *profile*.** PCR-free vs. PCR on the same samples drops the
stutter rate ~4.6× while *"the distribution of stutter error sizes stays about the
same"* (Gymrek lab, 2014); low-input / ancient DNA is the same story — elevated rate,
roughly the same ladder. Meanwhile the **rate grows linearly with repeat length** (forensic
LUS-linear models; Willems et al. in-vitro stutter, NAR 2019) and shrinks with motif
**period** (mono ≫ di ≫ tetra; Chakraborty et al., PNAS 1997). **But that is all human PCR
chemistry**, so for non-model / degraded DNA / unusual motifs we do *not* hard-code the
profile as cohort-invariant — we let the data decide via shrinkage:

- **stutter shape** is fit **per `(group, period)`, shrunk toward a cohort-per-period
  parent `θ_period`.** If the profile really is group-invariant (the PCR-free expectation),
  the shrinkage collapses the groups onto the shared shape at no cost; if a group's
  chemistry/preservation genuinely bends the profile, it keeps its own. The per-locus EM
  then refines `θ_locus` toward the sample's group-period shape (per-*motif* outliers —
  AT-dinucleotides, homopolymers — ride that per-locus refinement, a locus property). The
  per-`(group, period)` spread vs the parent is the **invariance diagnostic**. (Previously
  shape was cohort-shared by assumption — M3 was that the architecture then *couldn't*
  express a per-group difference even if one existed.)
- **stutter level** is the part that *does* vary by group (protocol/preservation) and
  by length → estimate it **per sample group** as a line in length,
  `level_baseline + level_slope·length` (both coefficients per group). The pre-pass
  **seeds** it; the genotyping outer loop **refines** it per group (C2 amend.) — it is an
  `S_θ` re-weight outside `align`, so refining it costs no cache rebuild and lets the cohort
  overrule a level the pre-pass got wrong.
- **`ε`** is a group property, **frozen** (no biological reason to vary it
  locus-to-locus; freezing makes `align` a pure, iteration-invariant function — *cacheable*
  if a cache is ever built, but v1 recomputes (cache deferred, genotyping §4); and, post-C1,
  `ε` is substitutions-only, so it is *not* the parameter a masquerading het contaminates).

This is a cohort hierarchy of our own — per-`(group, period)` *shape* (shrunk to a
cohort-per-period parent) + per-group *level* + per-locus refinement — **distinct from**
HipSTR's stutter model (fit **per locus, pooled over all samples**, full shape re-estimated
each locus, **no per-group axis**; verified in `em_stutter_genotyper.cpp`). We borrow
HipSTR's slip-placement-marginalized `align` — an **explicit `Σ_v` sum over the
placement-distinct variants** of `cand ⊕ Δ` (runs between interruptions, uniform position
prior, equal-LL placements collapsed; verify-fix #3, genotyping §3/§4) — its in-tract no-gap
rule (C1), and its per-locus free-shape refinement idea (`θ_locus`); the **sample-group axis on the shape**
(M3) is the part beyond HipSTR — which pools all samples per locus and so *cannot* separate
groups — and the frozen-`ε` determinism is ours, justified by cohort recurrence.

```
            ┌──────────────── ESTIMATION PRE-PASS (this phase) ────────────────┐
 CohortLocus│  resolve confident GENOTYPES  (hom ∪ separated het, off rungs)  │
  stream  ─►│      │  (SEED only — final estimate = soft full-cohort EM, C2)    │ ─► FROZEN ε        (per sample group)
 (Phase 1,  │      ▼  a resolved (A,A) OR (A,B≥A+2) makes each read a LABELLED   │     stutter shape   (per (group,period), shrunk to period parent)
  subset)   │  skirt = stutter obs (→ shape+level)  +  base mismatch (→ ε)      │     stutter level⁰  (SEED → refined in outer loop)
            │      │  (het: outer skirts hard-labelled, inner valley → soft EM) │       (CG-seed amend.)
            │      │  accumulate (after grouping): (group,period)→shape         │     G₀ pseudocounts (per loci group)
            │      │                              + sample→(level,ε)            │
            │  burn-in (settle) → measure (stratified loci)                    │     π⁰, θ⁰, F⁰  seeds
            │                    → cluster samples → sample groups              │       │
            └──────────────────────────────────────────────────────────────────┘       ▼  to Phase 3 EM
```

---

## 2. The engine — confident *genotypes* (homozygotes ∪ well-separated hets)

> **Amendment (CG-seed, 2026-06-22) — the chemistry seed is the confident *genotype*, not
> only the confident *homozygote*.** Estimating `ε`/stutter needs reads whose true parent
> allele is known, so a slip or a base mismatch is a *labelled* observation. A confident
> **homozygote** gives that (every non-`A` read is a labelled slip of `A`). But so does a
> confident, **cleanly-resolvable heterozygote**: when a diploid's two alleles are **≥ 2
> repeat units apart**, its length distribution is **two clear maxima separated by a real
> valley**, and each read attributes to the nearer allele — the *outer* skirts (below the low
> allele, above the high) are cleanly labelled; the *inner* valley is a mixture handed to the
> soft EM (C2). Two consequences: (i) it **closes m2(a)** — a hyper-heterozygous outbred
> cohort, where homozygotes are scarce *everywhere*, is exactly where separated hets are
> **abundant**, so the seed is data-rich where it used to starve; (ii) a het gives **two
> length anchors** (its two allele lengths), sharpening the level-vs-length slope (§1).
> Generalizes to ploidy *p*: a confident genotype resolves into **1..p clear peaks** (*p*-ploid
> resolution is harder — below). This amendment rewrites this section and propagates to
> §1/§3/§4/§6/§9, spec §4.3/§4.4, genotyping §2/§5.

Both `ε` and the stutter parameters come off **confident genotypes** (spec §4.3, §4.4): a
genotype whose alleles are individually resolved makes **every read a labelled observation**
of *which allele it came from and how far it slipped*. The length skirt (`−2,−1,+1,…` from a
known parent) gives the stutter (shape *and* level), the within-tract base mismatches give
`ε`. They are **jointly identifiable** because they move in different coordinates: **stutter
changes length (whole motif units), `ε` changes composition (substitutions)** — `align`
admits **no sub-motif indel inside the tract** (genotyping §4 / spec §6), so the only
length-changing operation is the slip and `ε` cannot absorb it. (At **period 1** a unit *is*
a base, so the single-base length change is stutter's; mononucleotide loci are flagged,
substitutions-only `ε` — spec §6.)

**Two kinds of confident genotype, both fully labelled:**
- **Confident homozygote** (`A,A` — one clear peak): every non-`A` read is a labelled slip of
  `A`; its whole skirt is "outer" (no second allele to confuse it). The original case.
- **Confident well-separated het** (`A,B`, `B = A+k`, `k ≥ 2` — two clear peaks): reads
  **below `A`** are `A`'s down-stutter, reads **above `B`** are `B`'s up-stutter — both
  cleanly labelled. Reads **between `A` and `B`** are a mixture of `A`-up and `B`-down stutter
  and are **not hard-labelled** — the soft full-cohort EM (C2) splits them by responsibility.
  Each faithful peak is contaminated only by the *other* allele's rare `≥k`-step stutter, so
  the base-mismatch `ε` count taken off each peak is nearly clean.

**Why `k ≥ 2` (and cleaner above it).** At **`k = 1`** the peaks are adjacent: `A`'s dominant
`+1` stutter lands on `B`'s faithful peak and `B`'s dominant `−1` on `A`'s — stutter and
allele are **fully confounded** (the *masked het*), so a 1-apart het is **useless and
excluded**. At **`k = 2`** the faithful peaks straddle a valley dominated by the alleles' `±1`
stutter and are contaminated only by the rare `±2`, so they are admissible; at **`k ≥ 3`**
almost the entire signal is cleanly attributable. So `k ≥ 2` is the floor — inner-region reads
are weighted by how cleanly they attribute (outer skirts hard, valley to the soft EM), or a
higher `k` is required for the *seed* (Q-P7 calibration).

**Identifying confident genotypes** (spec §4.3 "putative genotype off rungs"): for each sample
at each locus, find the **clear local maxima** over the rungs and resolve the genotype into
**1..ploidy peaks**, admitting it as a seed only when all hold:
1. **separation** — peaks pairwise **≥ 2 units apart** (no merged/masked pair);
2. **dosage-consistent heights** — the peaks' relative heights match an integer allele dosage
   (two comparable peaks for a balanced diploid het; one dominant peak for a homozygote). This
   rejects a **homozygote-plus-heavy-stutter masquerading as a het** (one tall peak + a small
   `−2` satellite is *not* two comparable peaks);
3. **cohort recurrence** — **each** resolved allele recurs as a clear peak in **≥ k samples**.
   A genuine allele recurs across the cohort; a stutter artifact does not recur as an
   *independent* allele. This is the strongest guard against the hom-plus-stutter masquerade,
   and it is **strongest precisely in the het-heavy cohort** this serves (many samples to
   corroborate against);
4. **depth** — a minimum depth to resolve the peaks at all (depth is duplicate-free upstream,
   verify-fix #5).

A genotype that fails resolution (peaks `< 2` apart, dosage-inconsistent, non-recurrent, or
too thin) is **not a seed** — it is left to the soft EM, never *forced* into a label. (This
**supersedes** the old "one clear maximum → homozygous, deliberately mislabelling merged hets"
heuristic: we no longer *mislabel* a merged het as a homozygote — a merged 1-apart het simply
fails resolution and is not seeded. That removes, for the separated case, the very
contamination C2 was patching.)

**Model-based resolution test (Q-P7, generalized).** Make the heuristic *statistical*: score
the reads under the best **1-peak** model vs the best **2-peak** (… up to *p*-peak) model — the
called peak(s) plus their predicted stutter skirts — and admit the genotype at the peak count
that **earns its parameters** (a BIC-style complexity penalty; no hand-picked p-value). This is
the **same one-vs-two-allele machinery** the prior design already specified, now read in *both*
directions: it rejects a masked (1-apart) het that can't support two peaks **and** confirms a
genuine two-peak het that does. The expected skirts come from the burn-in's current params
(co-evolves, below). *Open (Q-P7):* the peak-count test + the separation threshold for the
*seed* (calibrated for clean labels, not maximal yield).

**Selection and parameters co-evolve — one loop, not two passes.** The gate's model itself
assumes a stutter level, so a *fixed* gate would favour clean genotypes and bias a high-stutter
sample's estimate *cleaner* than the truth (selection bias, not just variance —
precision-weighting can't fix it). So the gate is **not fixed**: it is bootstrapped by the
**burn-in loop (§4)** — start from the dev-computed defaults, admit confident genotypes with
the current parameters, recompute the parameters from what was admitted, feed them back into the
resolution test, repeat until they settle. Selection and estimation are the *same* adaptive
loop. (This replaces the earlier "iterate ×~2" framing — the loop runs to its settle criterion,
§4.)

> **Amendment (2026-06-22, C2; refined by CG-seed) — the gate is a *seed*, the estimator is the
> soft EM reduce.** A *hard* homozygote-only gate admitted **masquerading hets** whose minor
> allele was counted as **stutter**, inflating the level (and an over-estimate admitted *more*
> such hets — positive feedback). Two changes neutralize it: (i) **CG-seed** now **resolves**
> well-separated hets *as hets* — their minor allele is scored as an **allele, not stutter** —
> so the dominant masquerade (a separated het mislabelled as a hom) is **gone by construction**;
> the residual masquerade (a hom+heavy-stutter mis-read as a het) is caught by the
> dosage-height + recurrence guards above. (ii) The **final** `ε` and level still come from a
> **soft full-cohort EM responsibility reduce** over the subset — every read contributes
> *fractionally*, so population recurrence pulls any remaining mislabel onto the right genotype
> rather than the skirt. The broader, cleaner seed just starts that EM in a far better basin,
> especially in het-heavy cohorts. The contaminated quantity (the **level**, not `ε` post-C1) is
> additionally **un-frozen** into the genotyping outer loop (§1, §4) so the whole cohort can
> overrule it — at no cache cost (it's an `S_θ` re-weight).

> **Polyploids — the principle holds, the resolution is harder (CG-seed).** For ploidy *p* a
> confident genotype is **1..p clear peaks** with integer dosages summing to *p* (tetraploid
> `AAAB` = 2 peaks at 3:1, `AABB` = 2:2, `ABCD` = 1:1:1:1). The same separation + recurrence
> guards apply, but **inferring the integer dosage from peak heights is genuinely harder**
> (overlapping skirts, depth noise), so v1 seeds only from the **cleanly-resolvable** polyploid
> cases; when the peaks can't be resolved into a confident dosage it **does not seed from that
> sample/locus** — it **leans on the app's coded priors** (literature stutter/`ε` defaults) plus
> the soft EM + cohort recurrence. Polyploid seed resolution is a **documented approximation**,
> the first place to look if a polyploid cohort's chemistry estimates are poor.
> **Owner:** the coded literature `ε`/stutter defaults are the **same dev-computed defaults the
> burn-in already starts from** (§2/§4); `rungs.rs` returns "unresolved" for the sample/locus
> and `prepass.rs` simply emits **no labelled skirt** for it (the soft EM + recurrence carry it),
> so no new storage is introduced — the defaults already live as the burn-in seed.

> **m2(a) — now a far narrower corner (CG-seed).** The cohort-wide "no confident *homozygote*
> anywhere" case (a hyper-het outbred cohort) is **largely closed** — those cohorts are rich in
> confident *separated hets*, which now seed the chemistry. The genuinely degenerate residual —
> **no confident genotype of *either* kind** (essentially every sample is a 1-apart merged het
> or too thin to resolve) — falls back to the **app's coded priors** (literature `ε`/stutter),
> lets the soft EM extract what cohort recurrence allows, and emits a **loud "chemistry not
> estimated from data — running on literature defaults" warning** (an output obligation, like
> the apparent-`F_IS` warning), with the confident-genotype count surfaced as a diagnostic. This
> is the only path on which the "estimate from data" thesis falls back to constants, and it is
> now **announced, not silent**. **Owner:** `prepass.rs` detects the cohort-wide-zero
> confident-genotype condition (the count it already tracks per §3) and raises a flag on the
> Phase-2 output; the **driver / `vcf_out.rs`** emits the warning into the **VCF header + stderr**
> (mirroring where the apparent-`F_IS` label is emitted, genotyping §6), and the per-sample
> confident-genotype count rides the §4.4 reporting block.

> **The shared-primitive boundary (Q-P1, resolved).** "Build rungs from the pooled cohort
> distribution" and "find clear local maxima" are needed **both here** (confident-genotype
> resolution) **and in Phase 3** (S1 candidate assembly, spec §5) — the same `pool → rungs`
> machinery. **Resolved:** it lives in a shared `rungs.rs` both phases call (§9 Q-P1), so the
> peak definition is identical across phases — a confident genotype's peaks are peaks candidate
> assembly would also see.

---

## 3. Loci groups & the commutative reduce

The pre-pass accumulates confident-**genotype** sufficient statistics (homozygotes ∪
well-separated hets — CG-seed, §2) into **two order-independent (commutative) keyed
accumulators**, reflecting the shape/level split. **A het contributes *per resolved allele*:**
each of its two peaks is a labelled parent, so it adds skirt/level/`ε` observations for
*both* allele lengths (the two outer skirts hard-labelled; the inner-valley reads enter via
the soft-EM responsibilities, not as integer counts — §2). So one separated het feeds the
accumulators at **two** length bins — the bonus that sharpens the level-vs-length slope.

1. **`(sample group, period) → SlipProfile`** — the **per-`(group, period)` stutter shape**,
   shrunk to a cohort-per-period parent (M3). Every confident genotype contributes its
   (conditional) slip offsets — a homozygote one labelled skirt, a separated het one per
   resolved allele — to its **`(group, period)`** profile (and the parent pools over groups).
   Because the key includes the sample group, this accumulation runs **after clustering fixes
   the groups** (§4); before that, the parent per-period profile is accumulated cohort-wide.
   Length is **not** a shape covariate — it drives the *level* (§1).
2. **`sample → SampleStutterStats`** — the per-sample **level** (faithful-vs-slipped
   counts binned by repeat length, to fit `level_baseline + level_slope·length`; a het adds a
   bin at **each** allele length) and **`ε`** (base-match/mismatch off each faithful peak).
   After accumulation, samples are clustered into **sample groups** (§4) and the level slope +
   `ε` are fit **per group**.

`G₀`'s own per-loci-group accumulator (the geometric decay on `π`, §5) is separate from
either of these.

```rust
// shape sketch — not final
struct ShapeKey { group: u8, period: u8 }  // M3: shape keyed by (sample group, period);
                                           //   period alone keys the shrinkage PARENT

struct SlipProfile {                // per (group, period); the conditional profile of a slip
    down: [u64; MAX_SLIP],          // −1, −2, …  -> contraction bias + geometric decay ρ
    up:   [u64; MAX_SLIP],          // +1, +2, …  -> up rate / asymmetry
}

// per-sample: feeds (a) the per-group LEVEL slope, (b) the sample-group clustering,
//             and (c) — once groups are known — the (group,period) SHAPE key
struct SampleStutterStats {
    by_length: Vec<(u16 /*repeat length*/, u64 /*faithful*/, u64 /*slipped*/)>, // -> level_baseline + level_slope·length
    base_match: u64, base_mismatch: u64,   // within-tract -> ε (substitutions, C1)
    depth: u64,                            // PRECISION weight for clustering only, never a feature;
                                           //   duplicate-free by construction (markdup upstream, Stage 1 skips flag 0x400 — verify-fix #5)
}
// reduce (all commutative integer sums, merged in fixed catalog order -> determinism):
//   (0) period          -> SlipProfile  : cohort-per-period SHAPE PARENT θ_period (pre-grouping)
//   (1) (group, period) -> SlipProfile  : per-(group,period) SHAPE, shrunk to the parent (M3; post-grouping)
//   (2) Sample          -> SampleStutterStats : per-sample LEVEL + ε
//       samples clustered into SAMPLE GROUPS; level slope + ε per group; then (1) is keyed by group
```

**Parallel & deterministic (spec §4.4).** Workers accumulate into **per-thread
partials**, reduced in a **fixed (catalog) order** — a **commutative reduce**,
order-independent. (A locked *running overwrite* would re-introduce order-dependence;
only the *accumulation* is shared.) The keyed structure is respected throughout — one
profile per `(group, period)` (plus the per-period parent), one stat block per sample,
never one global value.

**The same discipline covers the *decision* floats (verify-fix #1), not just these integer
counts.** The M4 amendment added three floating-point quantities that drive discrete
pre-pass decisions and would otherwise be thread-count-dependent: (1) **`ℓ_pen`** (the
burn-in plateau stop), (2) the **BIC per-model log-likelihoods** (ε-freeze, group
split/merge, group count), and (3) the **clustering distances**. (1)+(2) are sums of
per-locus E-step log-normalizers reduced by the **same fixed-point integer accumulation** as
`F_i`/level (scale → round → `i128`), so `ℓ_pen` and every `Δℓ_pen` is byte-identical; the
BIC comparison is then `Δℓ_pen` vs a *deterministic* penalty. (3) is computed from the
per-sample `(ε, level)` (themselves fixed-point reduces) and resolved with **deterministic
tie-breaks** (strict `<` thresholds + a total order on the sample catalog index, §4). So the
pre-pass is order-independent at its **decision layer** — stop iteration, freeze/merge, group
count and membership — not only in its sufficient statistics.

---

## 4. The pre-pass: burn-in → measure → cluster (spec §4.4)

*(The pre-pass's three parts go **by name**, not by number, on purpose: the pipeline
already numbers **Phase** 1/2/3 (reading / this pre-pass / genotyping) and the algorithm
numbers its **steps** S1–S5, so a third 1/2/3 scheme would collide. All three parts
below live inside Phase 2.)*

- **Burn-in (settle).** An **adaptive loop** that bootstraps the selection model from
  the **dev-computed defaults** (the empirically-seeded starting parameters of §2) until
  the parameters stop moving. Each cycle:
  1. **draw a batch** of loci from the catalog in a **seeded random order** (the seed
     fixes the draw → reproducible *and* a representative spread across the length range,
     not just the first/cleanest loci);
  2. **parallel map** — a worker-pool maps each locus, running the §2 goodness-of-fit
     test and accumulating its sufficient statistics, **using the batch's frozen
     parameters** (a pure function of `(locus, params)` — workers never see each other's
     updates mid-batch);
  3. **barrier → reduce** — combine the per-worker partials in **fixed catalog order**
     (the commutative reduce of §3);
  4. **update** the parameters from the confident genotypes admitted so far (homozygotes ∪
     well-separated hets — CG-seed §2), and feed them back into the next batch's resolution model.

  **Stop point (M4 amend. 2026-06-22): the penalized marginal log-likelihood plateau.**
  Because post-C2 the pre-pass is a **soft full-cohort EM**, it ascends one well-defined
  objective — `ℓ_pen` = the marginal (observed-data) log-likelihood of the subset reads
  under the current params **plus the log-priors** (shape-shrinkage + `G₀`). EM/MM makes it
  **monotone**, and it is the **E-step's own normalizer** (`log Σ`), so it costs nothing
  extra. **Stop when `Δℓ_pen/|ℓ_pen| < tol`** with a max-iteration cap — *not* "params/calls
  stopped moving" (`ℓ_pen` measures data fit, so it is non-circular, unlike the retired
  "calls don't move"). A **non-monotone / never-plateauing `ℓ_pen`** is a sharp diagnostic
  (the loop isn't a proper ascent → a bug or model mismatch). The **batch size** is a
  speed knob (small = fine online steps; large = fills the pool), and the *stop decision*
  uses `ℓ_pen` on the full fixed subset. We **do not assert** batch size leaves the plateau
  unchanged — mini-batch EM paths *can* land in different optima — so **batch-size invariance
  of the recovered params is a simulator check** (vary it, confirm the `ℓ_pen` plateau and
  the recovered chemistry agree, F2), with **multi-start** (M4) as the guard that picks the
  best `ℓ_pen` regardless of path.
  **Determinism (verify-fix #1):** the seed fixes batch membership and the per-batch map is
  pure — but `ℓ_pen` is itself a **float sum of per-locus E-step log-normalizers over loci
  that complete out of order**, so it is reduced by the **same fixed-point integer
  accumulation** as `F_i`/level (scale → round → sum into an `i128`; §3, §7), making it
  byte-identical across thread counts. Hence the trajectory *and* the **stop iteration**
  (`Δℓ_pen/|ℓ_pen| < tol`) are **reproducible across thread counts** — not merely asserted:
  a borderline plateau can no longer drift by a thread-count-dependent last bit. (No `i128`
  overflow — per-locus normalizers are `reads × O(few nats)`, so even `10¹²` loci stay
  ≪ `i128::MAX`; precision `2⁻⁴⁰` nats ≪ any `tol`.)
  **`ℓ_pen` plateau = converged, not correct.** A self-consistent *wrong* fixed point can
  have high `ℓ_pen`, so correctness is anchored **outside** the loop: **multi-start** (run
  from several seeds, keep the best `ℓ_pen`; divergent basins = a flag, not silently
  resolved) + the simulator ground-truth recovery / known-protocol positive control (§9).
  `ℓ_pen` is the *convergence* check; those are the *correctness* gates. **Never settling
  is a diagnostic, not a verdict** — *which assumption broke* (too little data / unmodelled
  structure / a bug), not "the method is wrong."
- **Measure.** From the settled value, estimate each parameter over a
  **representative, stratified** sample of loci (span the length range, not just the
  cleanest), building the **distribution of its per-locus estimates**; the frozen
  value is the **average (mean)** of that distribution. The distribution's *shape* is
  itself a **development-time diagnostic**: tight & unimodal validates pooling at that
  granularity; wide or multimodal flags unmodelled structure (a covariate too coarse,
  a mis-grouped sample, a bug) — the "never settling is a diagnostic" principle one
  level up (inspected during implementation/validation, not at runtime). The three
  parameters:
  - **stutter shape → the per-`(group, period)` profile** (M3), the prior for the per-locus
    refinement. The hierarchy shrinks each `(group, period)` cell toward its **cohort-per-period
    parent**, and a sparse/empty period toward a pooled-across-periods grandparent. Two
    kinds of real heterogeneity are now *both* representable: per-**group** (protocol/
    preservation — the M3 axis, captured when a group's cell earns its own profile) and
    per-*motif* (AT-di / homopolymer outliers — the **per-locus refinement in Phase 3
    carries** these). The per-`(group, period)` spread vs the parent is the **invariance
    diagnostic** (does the profile actually vary by group on this cohort?). Fit **after**
    clustering fixes the groups (§3).
  - **stutter level → the per-sample rate vs. length,** fit as a line
    `level_baseline + level_slope·length` (continuous in length — no length bins; the
    literature makes rate linear in length, and bins only add edge discontinuities), from
    the **soft** per-allele responsibilities (C2, not the hard gate). These per-sample lines
    are the clustering input; both coefficients are then fit **per sample group** in the
    **cluster** part below — but the per-group level is the outer-loop **seed** (`level⁰`),
    **not** frozen: the genotyping outer loop refines it (C2, §1).
  - **`ε` → frozen per sample group,** licensed by a **penalized-likelihood model
    comparison** (M4 amend. — *not* "±1 SD flips no call"): freezing is justified iff
    `ℓ_pen(ε + covariate / per-sample ε)` does **not** beat `ℓ_pen(frozen ε)` by more than
    the BIC complexity penalty — i.e. the richer `ε` model doesn't earn its parameters on
    the data. This references the fit, not whether the output flips, so it is non-circular.
    If the richer model *does* win, that's the cue to give `ε` a covariate (a documented
    upgrade). `ε` is **per group** (not per-sample) because it lives in `align` (Phase 3
    §4/§7); per-sample `ε` is the upgrade.

  **Cluster into sample groups.** Group **mainly for data sufficiency** — a low-coverage
  sample can't support reliable per-sample chemistry, so it borrows strength from similar
  samples. Grouping is purely an **estimation + reporting** device with **no cache role**
  (the `align` cache is deferred — genotyping §4/Q-G3 — so there is no cache constraint to
  satisfy either way). Mechanism: **distance-based grouping of close
  neighbours** in (`ε`, level) space (similar params ⇒ similar EM behaviour ⇒ safe to
  pool) — **not k-means** (deterministic, no random init; number of groups falls out of a
  penalized-likelihood comparison (M4), not a preset K; single-protocol cohort collapses to one).
  **Deterministic under ties (verify-fix #1):** strict `<` thresholds with a defined
  equal-handling rule and a **total tie-break on the sample's catalog index**, so a borderline
  `distance ≈ threshold` (or two equidistant merge candidates) resolves identically across
  thread counts; the per-sample `(ε, level)` it groups on are order-independent fixed-point
  reduces, and the BIC split test is a `Δ` of fixed-point `ℓ_pen` values vs a deterministic
  penalty. No user labels (SRA provenance unreliable). **Distances are scaled by each sample's uncertainty**
  so a wobbly thin sample is "close to many things," never confidently misplaced — this is
  how grouping protects the samples it exists to help; **depth enters only here, as that
  precision** (∝ 1/√depth), **never as a feature or regressor** (the coverage↔chemistry
  correlation is real preservation signal). Depth is **duplicate-free by construction** —
  duplicates are marked upstream (post-mapping markdup) and Stage 1 skips flag-`0x400` reads,
  so this caller does no dedup (verify-fix #5; residual unmarked jackpotting is an
  input-quality matter, out of scope). Two regimes are split
  only when the split **earns its parameters** — `ℓ_pen(split) − ℓ_pen(merged)` clears the
  BIC penalty (M4 amend.; *not* "doesn't move calls") — *and* each group has **enough data**.
  `ε` and level are grouped **jointly** (one group = one chemistry + preservation regime).
  A **singleton group** (m2 edge case) is fine: it uses its own shrunk estimate, or — if too
  thin for that — the partial-pooling shrinks it toward its nearest neighbour / the cohort,
  so it never estimates from nothing; the BIC test won't carve one out unless the data
  support it. **Report** per-sample / per-group values + memberships for audit.
  *Remaining (validation, not gates):*
  the **discrete-vs-continuum** look at the (`ε`, level) cloud + a known-protocol positive
  control (Q-P6 / spec §9). The level carries **both** `level_baseline` and `level_slope`
  per sample group (the granularity default, §0); pooling `level_slope` to the cohort is a
  *future* move only if the data show the length response doesn't vary by group.

> **Q-P2 (resolved; M4 amend. 2026-06-22) — the `ε`-freeze check is a penalized-likelihood
> model comparison, not a "calls don't flip" test.** Licensing the freeze by "perturb ±1 SD,
> no call flips" was **circular** (it asks whether the model's own *output* moves, not
> whether the data support a richer `ε`). The principled check: fit `ε` at the candidate
> granularities — frozen-per-group vs. `ε + covariate` / per-sample — on the subset and
> compare their **penalized marginal log-likelihoods**; **freezing is justified iff the
> richer model does not beat the BIC complexity penalty.** Cheap because the EM (and hence
> `ℓ_pen`, the E-step normalizer) is shared code (Q-P1); run **per cohort** as a *checked
> invariant* (richer model wins → warn + point at the `ε`-covariate upgrade). References the
> data fit, so non-circular (M4).

---

## 5. The `G₀` pseudocount prior on `π` (spec §4.3, §5.5)

A candidate's `π_i` enters the genotype prior **multiplicatively**, so `π_i = 0` is
an **absorbing trap** — and a candidate *can* land at 0 in the seed (the candidate
set is recall-generous, broader than the confident-genotype tally; the chief case is
the **masked het**). The fix is **pseudocounts**: every candidate starts non-zero,
the EM grows it if the data warrant, the pseudocount vanishes under real evidence.
This is what makes the seed's deliberate merged-het mislabelling **recoverable**.

Shape (spec §4.3, settled; the *decay parameter* is fit here; the *functional form* is
**geometric for v1** (settled), with a **Gaussian form** the open §9 upgrade):

- **geometric decay** in the unit offset `Δ` (heavier-than-Gaussian tail keeps real
  large-step alleles callable),
- **centred on the per-locus cohort modal allele** — *not* the reference (the
  reference is only the `Δ` frame; centring on it would re-import the reference bias
  Mark-2 exists to kill, and degrades gracefully — with no data the only candidate is
  the seeded reference, so the mode *is* the reference),
- **symmetric** in v1, **purity-agnostic** (impure allele of length L gets the same
  prior as a pure one of length L — no impurity penalty, §7 of the spec),
- **decay parameter fit from pooled empirical data per loci group** (parametric,
  deliberately not data-hungry non-parametric),
- **floored at a tiny positive constant (verify-fix #4):** `G₀[i] = max(p^|Δ|, FLOOR)`.
  `p^|Δ|` underflows to **exactly `0.0` in f64** when `|Δ|·(−ln p) ≳ 744`; candidates are
  read-length-capped so short reads are safe, but a **long read** over a large tract with a
  **steep decay** can cross it and re-create the `π = 0` absorbing trap *numerically* — for
  exactly the far candidate (large `|Δ|`) the prior exists to keep alive. The floor (any
  representable `> 0`) re-floors it; equivalently carry `G₀` in the **log domain** into the
  M-step (the cleaner form, since the engine sums `expected_counts[a] + pseudocount[a]`).

It is a **false-positive control** *and* a small-N regularizer, re-entering **every**
M-step (not only the seed). This is the §5.5 Dirichlet base measure.

> **The engine boundary (settled spec §3 / §9, restated for the layout).** The reused
> `posterior_engine.rs` picks Dirichlet pseudocounts by SNP/indel **allele class**
> ([`classify_allele`](../../../src/var_calling/posterior_engine.rs); REF 10 /
> SNP-alt 0.01 / …). SSR **replaces** that with a **per-candidate pseudocount
> vector** computed here from `G₀`. Phase 3 owns generalizing the engine to *accept*
> the vector; this phase owns *producing* it. (Open: extend the engine in place vs a
> slim SSR fork — a Phase-3 struct-shape call.)

---

## 6. The seeds — `π⁰`, `θ⁰`, `F⁰` (spec §4.3)

**Where they're computed (Q-P3, resolved):** the per-locus seeds are produced in
**Phase 3**, at each locus's EM start, reusing the shared `rungs.rs` — *not* in this
pre-pass, which runs on a subset and so can't seed loci it never visits. This phase
supplies only the *global* inputs the seeds need (frozen `ε`, stutter shape, level, `G₀`
decay). What each seed is:

- **`π⁰` tally.** Give *every* sample a putative genotype off its rungs (§2; one
  genotype contributes `ploidy` allele-copies regardless of depth), tally the
  alleles. Combine with the prior: **`π⁰ = (putative-genotype counts + G₀
  pseudocounts), normalized`.** With **no confident samples** at a locus it falls
  back to the **prior alone** (normalized `G₀`, centred on the modal candidate =
  the reference when there is no data) — the correct conservative small-N behaviour.
  This tally rides on the **same `pool → rungs` pass Phase 3 already runs for candidate
  assembly** — the real one-pass fusion.
- **`θ⁰`** = the period **stutter shape** for the locus × the sample's frozen
  **stutter level** evaluated at that locus's length (`level_baseline + level_slope·length`)
  (§3–§4) — derived entirely from the frozen global params + the locus length, so it
  needs no pass of its own.
- **`F⁰`** = a supplied / global default; it is the **burn-in start** of the Phase-3
  prior-side `F_i` loop, **not estimated here** (it needs genotypes — Phase 3).

The seed only has to point the EM at the right basin; the pseudocount floor (§5) is
what makes an optimistic seed safe. **Seed + prior are a pair** — the seed is allowed
to be optimistic precisely because the prior keeps the alleles it skips alive.

> **Q-P3 (resolved) — seeds computed in Phase 3.** `π⁰`/`θ⁰` are *per-locus* and needed
> for **every** locus, but the pre-pass visits only a subset — so it can't seed loci it
> never sees. They're therefore computed **in Phase 3** at each locus's EM start (shared
> `rungs.rs`, Q-P1), with this phase owning only the *global* `ε`/shape/level/`G₀`
> outputs. `π⁰` fuses with candidate assembly's `pool → rungs` pass; `θ⁰` is just globals
> × the locus length.

---

## 7. Parallelism, memory, determinism

- **Parallelism** — a simple **locus-iterator → worker-pool map** (§4): each worker
  maps a locus to **sufficient statistics**, not calls, using the batch's frozen
  parameters; the reduce is the commutative merge of the two accumulators (§3). Burn-in
  runs this loop per batch; measure runs it once over the seeded subset — both cheaper
  than a full genotyping pass.
- **Memory** — holds per-thread `period → SlipProfile` partials **plus per-thread
  `Sample → SampleStutterStats` partials** (both tiny: a handful of periods + N
  samples × small stat structs) + whatever `CohortLocus`es are in flight from the
  shared merge. No per-locus `align` cache yet (that is Phase 3).
- **Determinism** — the whole pre-pass is deterministic: the burn-in loop is
  reproducible because its **batches are drawn in a seeded order**, each batch's map is a
  **pure function of `(locus, frozen params)`**, and parameters update only **at the
  barrier** (never a per-locus running overwrite) → the trajectory and settle point are
  identical across thread counts (§4); measure is the same per-thread accumulation +
  fixed-order commutative reduce (§3); and the sample-group clustering is a
  **deterministic fit over the per-sample results** (§4). **The decision floats reduce
  order-independently too (verify-fix #1):** `ℓ_pen` (the plateau stop) and the BIC
  per-model log-likelihoods (ε-freeze, group split/merge/count) are summed by the **same
  fixed-point integer accumulation** as `F_i`/level, so the *stop iteration* and the
  *freeze/merge/group-count* decisions are byte-identical across thread counts; the
  clustering resolves borderline ties on the **sample catalog index** (§3/§4) — so the
  pre-pass's *discrete decisions*, not only its statistics, are thread-count-invariant. Output (`ε` *per group* frozen,
  cohort stutter shape, per-group level **seed** + group assignments, `G₀` params) is a
  pure function of the inputs, identical across thread counts — freezing `ε` per group is
  what carries that determinism into Phase 3's per-(locus, group) `align` cache; the
  per-group level **seed** is then refined in Phase 3's outer loop (C2), which stays
  deterministic via the **same order-independent fixed-point integer reduce** as `F` —
  loci complete out of order, so the responsibility sums use scaled-integer accumulation
  (the `u64`-counts trick this pre-pass already uses), byte-identical across threads (M1;
  spec §4.4).

---

## 8. Proposed module layout (sketch)

```
src/ssr/cohort/
  rungs.rs        # SHARED with Phase 3 (Q-P1): pool -> rungs -> clear-maxima; confident-GENOTYPE
                  #   resolution (hom ∪ separated het, 1..ploidy peaks — CG-seed §2); returns the
                  #   resolved peak set + per-peak allele, or "unresolved" (merged/dosage/recurrence/thin)
  prepass.rs      # the iterated estimation: from each resolved genotype, LABEL the skirts —
                  #   hom: one skirt; separated het: TWO outer skirts (per allele) hard-labelled,
                  #   inner valley -> soft-EM responsibilities (CG-seed §2/§3); skirt/ε accumulation
                  #   (cohort shape + per-sample level/ε), burn-in settle + stratified measure +
                  #   re-select loop, commutative reduce; coded-default fallback + m2(a) zero-flag  (§2-4)
  sample_groups.rs # cluster per-sample (ε, level) -> soft sample groups (K≤~8-10 selected;
                  #   precision = 1/√(depth), dup-free upstream); per-sample shrink to group; reporting  (§4)
  stutter.rs      # S_θ(Δ) kernel = (group,period)-shape × per-group level(length); shape fit +
                  #   (group,period)->period->global shrinkage (M3); per-group level-slope fit;
                  #   reachability ⊕ + placement-variant SET enumeration for impure A⊕Δ (verify-fix #3)
                  #   (shared w/ Phase 3: re-weight + likelihood.rs sums the Σ_v)
  base_measure.rs # G₀ geometric pseudocounts: per-loci-group decay fit + per-candidate vector  (§5)
  seed.rs         # π⁰ tally (putative genotypes + G₀), θ⁰, F⁰  (§6; may fold into Phase 3 — Q-P3)
```

(`stutter.rs`/`base_measure.rs`/`seed.rs` keep the names from
[ssr_genotyping_architecture.md §4](ssr_genotyping_architecture.md); `prepass.rs` and
`sample_groups.rs` are the Mark-2 additions that the frozen-`ε` / shared-shape /
per-group-level decision introduced.)

---

## 9. Resolved design decisions + remaining calibration

The seven Phase-2 questions were worked through and **resolved 2026-06-21**; what
remains under each is *calibration* (numbers to fix on the simulator), not architecture.

- **Q-P1 (resolved) — shared primitives, thin per-phase drivers.** Phase 2 and Phase 3
  do the *same* per-(sample, locus) work (pool → rungs → peaks, then align reads vs.
  candidates with the stutter kernel + `ε`) and differ only in the **reduce** — Phase 2
  → sufficient statistics → params/seeds; Phase 3 → EM → calls. So they **share a
  locus-primitive library** (`rungs.rs` = pool/rungs/clear-maxima, the align/pair-HMM,
  `stutter.rs`, the candidate structures, `base_measure.rs`, `seed.rs`), with thin
  `prepass.rs` / genotyping drivers on top. **Shared code, not merged passes** — Phase 3
  still consumes Phase 2's frozen outputs (frozen `ε` ⇒ `align` is a pure, iteration-
  invariant function ⇒ the run is deterministic; recomputed in v1, cache deferred). This *forces* one peak definition across phases: a confident
  genotype's resolved peaks (homozygote or separated het — CG-seed §2) are peaks candidate assembly would also see.

- **Q-P2 (resolved; M4 amend.) — the `ε`-freeze check is a penalized-likelihood model
  comparison.** *Not* "run calls at {`ε`, `ε±SD`} and confirm no flip" (circular — it tests
  whether the model's output moves, not whether the data want a richer `ε`). Instead compare
  the **penalized marginal log-likelihood** of frozen-per-group `ε` vs. `ε + covariate` /
  per-sample `ε` on the subset; **freeze iff the richer model doesn't beat the BIC penalty.**
  Cheap because `ℓ_pen` is the E-step normalizer the shared EM already computes (Q-P1). Run
  **per cohort** as a checked invariant — richer model wins → warn + point at the
  `ε`-covariate upgrade. Non-circular (references the fit, M4).

- **Q-P3 (resolved) — seeds computed in Phase 3, not the pre-pass.** Decisive reason:
  Phase 2 runs on a *subset*, but `π⁰`/`θ⁰` are needed per-locus for **every** locus — a
  subset can't seed loci it never visits. `π⁰` is built off the locus's peaks, which
  Phase 3's candidate assembly already produces (the real "one pass" fusion is
  candidate-assembly + seed, *in Phase 3*); `θ⁰` derives trivially from the frozen global
  params × the locus length. Phase 2 owns only the **global** outputs (`ε`, stutter
  shape, level, `G₀` decay) + the `F⁰` default. (§6's wording reflects this.)

- **Q-P4 (resolved) — dedicated batched pool; batch size 32, fixed.** The pre-pass is a
  **batched map-reduce-with-barrier** (the barrier — update params between batches — is
  exactly what Phase 1's barrier-free pipeline lacks); reuse the reader/merge to *source*
  the subset, not the streaming pipeline to *process* it. **Batch size = 32, a fixed
  constant — never tied to thread count** (tying it to cores would make results
  machine-dependent, defeating determinism). Recompute-after-batch is cheap (combine
  running sums + refit a few numbers) and estimates are **cumulative**, so small batches
  are stable — they just check more often. (Caveat: 32 leaves no slack on ≥32-thread
  pools; bump to 128–256 for large machines.) `(seed, batch size)` define the reproducible
  run; the *stop decision* uses `ℓ_pen` on the full fixed subset, so batch size is a speed
  knob on the trajectory. We **do not assert** it leaves the plateau unchanged (mini-batch EM
  can land in different optima) — **batch-size invariance of the recovered params is a
  simulator check** (m3 / F2), with **multi-start** (best `ℓ_pen`, divergent basins flagged —
  M4) as the path-robustness guard. *Remaining:* in-memory-vs-reread of the subset → reading Q-R5.

- **Q-P5 (resolved as calibration) — fix the numbers on the simulator.** Starting
  defaults from representative real data + literature (`ε⁰` from the Stage-1 quality
  gate); **settle** = `Δℓ_pen/|ℓ_pen| < tol` (penalized marginal log-lik plateau, M4) +
  a hard iteration cap (cap-hit / non-monotone `ℓ_pen` = the "never settled" diagnostic),
  with **multi-start** (best `ℓ_pen` wins, divergent basins flagged); **measure window** = a stratified sample
  with enough loci per (period, sample); **per-locus shrinkage** strength set from how
  much the shape actually varies across loci (data-driven, not hand-picked); **`G₀` =
  geometric** for v1 (heavier tail keeps real large jumps callable), Gaussian a
  documented upgrade. None touch the architecture.

- **Q-P6 (resolved) — distance-based grouping; cache may be per-sample.** Group **mainly
  for data sufficiency** (thin / low-coverage samples borrow strength from similar
  samples), not cache cost. Mechanism: **distance-based grouping of close neighbours in
  (`ε`, level) space — not k-means** — with distances **scaled by each sample's
  estimation uncertainty** (so a noisy thin sample is "close to many things," not
  confidently misplaced); deterministic (no random init). Number of groups falls out of a
  **penalized-likelihood comparison** (M4): split two regimes only when `ℓ_pen(split) −
  ℓ_pen(merged)` clears the BIC penalty for the extra per-group parameters, guarded below by
  *enough-data* (not the old *calls-don't-move*); the ~8–10 cap is a secondary budget.
  **No cache constraint applies** —
  the `align` cache is **deferred** (genotyping §4/Q-G3), so grouping needn't draw hard
  lines for any table; it's purely an estimation + reporting device. (The earlier claim
  that "the cache can be per-sample because shrinkage keeps distinct `ε` small" was wrong —
  shrinkage only collapses *thin* samples, so a deep continuum cohort kept ≈`N` distinct
  `ε`; M2. Deferring the cache moots it. *If* a cache is ever built, per-group-`ε`
  quantization bounds it, not shrinkage.)
  `level_slope` stays per-group (pool to cohort only if the BIC comparison favours it — M4).
  *Caveat:* on a smooth gradient, naïve nearest-neighbour merging can chain into one blob —
  the **penalized-likelihood split test** breaks an over-merged group when the data support
  two regimes. *Remaining (validation, not gates):* the discrete-vs-continuum look at the
  real (`ε`, level) cloud + a known-protocol positive control (groups track
  protocol/provenance, not coverage). Ties to Phase 3 Q-G3. **Depth is duplicate-free by
  construction** — markdup runs upstream (post-mapping) and Stage 1 skips flag-`0x400` reads;
  `ssr-call` does no dedup and models no residual jackpotting (verify-fix #5 — an
  input-quality matter, out of scope).

- **Q-P7 (resolved; generalized by CG-seed 2026-06-22) — confident genotype =
  1..ploidy-peak resolution test.** *Not* a generic goodness-of-fit (χ²): score the reads
  under the best **1-peak** model vs the best **2-peak** (… up to *p*-peak) model — called
  peak(s) + predicted skirts — and admit at the peak count that **earns its parameters**
  (BIC-style, no hand-picked p-value). **CG-seed change:** the test is now read in *both*
  directions — it **admits a confident well-separated het** (two peaks **≥ 2 units apart**,
  dosage-consistent heights, **each allele cohort-recurrent**) as a seed *as well as* a
  homozygote, instead of only confirming homozygotes — because a separated het's outer skirts
  are labelled stutter too (§2). A merged (1-apart) het **fails resolution → not seeded**
  (no longer mislabelled as a hom — removes the C2 contamination for the separated case);
  a hom+heavy-stutter masquerading as a het is rejected by the dosage-height + recurrence
  guards. Behaves at low read counts (unlike χ²) and reuses existing code. *Remaining:* tune
  the penalty + the seed separation threshold on the **simulator** (clean labels over yield —
  a sneaked-in *unresolved* het corrupts the estimates; a discarded confident genotype just
  costs a little data, and the broadened seed gives plenty). Plus a minimum-depth skip; the
  expected skirt comes from the burn-in's current params (co-evolves). **Polyploids:** the
  *p*-peak resolution is harder → seed only cleanly-resolvable cases, else lean on coded
  priors (§2). *(C2 amend.: the gate **seeds**; the final `ε`/level come from the soft
  full-cohort EM, so any genotype that sneaks past is still pulled onto the right genotype by
  population recurrence, not counted as stutter.)*

**Settled (spec §4.3/§4.4/§5.5 + 2026-06-21, amended 2026-06-22 = C2, M3), not re-opened
here:** `ε` frozen **per sample group** (no per-thread `ε`, no `δ`-rebuild — retired;
per-sample `ε` = upgrade; post-C1 it is substitutions-only); **stutter shape = per
`(group, period)`, shrunk toward a cohort-per-period parent, refined per locus** (M3 —
*not* assumed cohort-shared; the PCR-free "rate not profile" expectation is realized by
shrinkage so invariance is data-driven, not an assumption); **stutter level = per sample
group, linear in repeat length** (`level_baseline + level_slope·length`, both coefficients
per group; rate, not profile, is what protocol/preservation/length move *most*) —
**pre-pass-seeded then refined in the
Phase-3 outer loop, not frozen** (C2: an `S_θ` re-weight outside `align`); the pre-pass
**estimator is the soft full-cohort EM reduce**, the confident-**genotype** gate (homozygotes
∪ well-separated hets, 1..ploidy-peak resolution — CG-seed) a **seed**
(C2); length is a **slope, never bins**; **granularity default = sample group, cohort only
when data justify** (§0); depth is **precision-only** in the clustering, never a feature, and
**duplicate-free by construction** (markdup upstream, Stage 1 skips flag `0x400`; no
in-caller dedup — verify-fix #5); `G₀` geometric, cohort-mode-centred, purity-agnostic; commutative
deterministic reduce + deterministic clustering; `F⁰` only a seed, `F` (and now the level)
estimated in Phase 3.
