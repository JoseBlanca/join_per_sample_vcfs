# SSR Stage 2 — `ssr-call` genotyping (architecture sketch)

**Status:** draft, 2026-06-19 (chemistry model & vocabulary **settled 2026-06-21**),
branch `ssr-cohort`. **A discussion starter, not a
settled design** — the third of three Stage-2 (`ssr-call`) sketches, and the one
that does the calling proper. It consumes the `CohortLocus` stream from Phase 1 and
the frozen `ε` / stutter shape / stutter level / `G₀` priors from Phase 2, and
produces the cohort VCF. Companions:

- [ssr_call_reading.md](ssr_call_reading.md) — Phase 1: reader + k-way merge → `CohortLocus`.
- [ssr_call_parameters.md](ssr_call_parameters.md) — Phase 2: the pre-pass (frozen `ε`,
  the stutter-shape / stutter-level / `G₀` priors, the `π⁰`/`θ⁰`/`F⁰` seeds).

> **Chemistry model (settled 2026-06-21; refines the 2026-06-20 amendment) — chemistry
> is per sample group (spec §4.4).** `ε` and the stutter **level** are per **sample
> group**, not cohort-global (PCR vs PCR-free / fresh vs historic DNA); the stutter
> **shape** is per `(group, period)`, shrunk toward a cohort-per-period parent, refined per locus (M3). See
> [ssr_call_parameters.md §0](ssr_call_parameters.md) for the canonical vocabulary.
> The consequence **for this phase**: `ε` is frozen at **sample-group** granularity, so
> with `ε` fixed `align(obs | cand⊕Δ)` is **invariant across iterations** (only `S_θ`
> re-weights it) — cacheable in principle, but **v1 does not build the align cache** (it's
> an unmeasured optimization; v1 recomputes on demand, §4, Q-G3). Determinism comes from
> `align` being a **pure function of `(obs, cand⊕Δ, ε)`**, not from any cached table. The
> per-group **stutter level** (linear in
> repeat length, `level_baseline + level_slope·length`) rides on top as a cheap `S_θ`
> re-weight (no extra HMM, no cache rebuild — like `θ_locus`) and is therefore **refined in
> the prior-side outer loop, *not* frozen** (C2 amend., spec §4.4): pre-pass-seeded, then
> re-estimated per group from the cohort's soft responsibilities each outer round, under
> the same order-independent fixed-point integer reduce as `F` (M1, §5/§7). Only **`ε`** is
> frozen (it is *in* the cache);
> per-sample `ε` is a documented upgrade (key the cache by `ε`, accept the cost). Sections
> 1/4/5/7/9 below carry this; read bare "frozen `ε`" as "per-sample-group-frozen `ε`".

The *how-the-code-is-wired* companion to spec
[ssr_cohort_mark2.md](../specs/ssr_cohort_mark2.md) §4.2 (model framing), §5 (S1
candidate assembly), §6 (S3 likelihood), §7 (S2 reachability), §4.4 (`F` loop), and
§4.5 (output) — **all shapes settled 2026-06-19; this doc proposes the module/struct
layout those left deferred**. Where they disagree on intent, the spec wins; on code
layout, this doc wins. Refines the Mark-1-era `cohort/` sketch in
[ssr_genotyping_architecture.md §4/§6](ssr_genotyping_architecture.md).

---

## 1. The per-locus pipeline

Per `CohortLocus` (each independent — the unit of work-pool parallelism):

```
 CohortLocus  ─►  S1 ASSEMBLE  ─►  S3 SCORE           ─►  per-locus EM        ─►  CALL
 (Phase 1)        candidates A_ℓ    align(obs|cand⊕Δ)     refine π + θ_locus      genotype + posterior
                  (pool→rungs→        flat-ε HMM,          (per-grp ε; θ_locus      + allele-balance QUAL
 per-grp ε,       nominate, admit)    in-tract subst.      shape × per-group        → VCF record
 shape, G₀,                          closed-form;         level re-weight)
 per-grp level⁰─►                     recompute (cache    seeded π⁰/θ⁰ (Phase 2)
                                      deferred — §4)
 (Phase 2 seed)                                           F_i + per-grp level from outer loop (§5)
```

Steps S1–S3 are spec §5/§6/§7; the EM core reuses
[`posterior_engine.rs`](../../../src/var_calling/posterior_engine.rs) (spec §3). The
`F_i` loop (§5 here / spec §4.4) wraps the whole per-locus pass.

---

## 2. S1 — candidate assembly (spec §5)

Three levels — **pool → rungs → candidates** — then a **locus-admission filter**.
Throughout: **recall is the goal, precision is the EM's job** — a sequence dropped
here is gone; a wrongly-kept candidate is driven to `π ≈ 0` by population recurrence +
`G₀`. Non-candidate sequences are **not discarded** — they re-enter the likelihood as
slip/error products (the §6 invariant).

```rust
// shape sketch — not final
struct CandidateSet {
    alleles: Vec<Box<[u8]>>,   // candidate sequences (each an independent allele)
    ref_idx: usize,            // the reference allele, seeded unconditionally
    admit: Admission,          // PASS | NotPeriodic | TooManyAlleles | LowDepth (-> site FILTER, §6 below)
}
```

1. **`ObservedSeqs` pool** — sum per-sample `(seq, count)` into one cohort-aggregate
   distribution (keep the per-sample dists for level 3).
2. **Rungs** — a length-position admitted by **recurrence** (length seen in ≥ *k*
   samples — catches stutter bands + stutter-masked alleles, separates structure from
   sporadic sequencing error) **or height** (a local maximum in *some* sample — admits
   rare private alleles). A rung is **length-keyed, holds a *set* of sequences** (each
   distinct same-length sequence clearing a cohort-frequency threshold is an
   independent allele; the below-threshold same-length cloud is per-base error, fed to
   the likelihood, not promoted).
3. **Locus-admission filter** — the distribution of **adjacent-rung length
   differences**; its **mode must be the catalog motif length**. If not (no coherent
   period, junk) → **no-called** as a *filtered* record (`notPeriodic`), never silently
   dropped. Catches a locus clean in the *reference* but non-periodic in the
   *population* (reference-only curation structurally cannot see this).
4. **`CandidateAlleles` — per-sample nomination, unioned.** Per sample, count **clear
   local maxima** (`> 3 reads above each neighbor`, a prominence floor — default to
   confirm): **≥ ploidy maxima** → top-ploidy are this sample's candidates, no rescue;
   **< ploidy** → add the **±1 neighbors' observed sequences** (already in the cohort
   scaffold — *include*, never synthesize) so the EM can resolve a merged het. Union
   across samples; **reference seeded unconditionally**; cap at `MAX_CANDIDATE_ALLELES`
   (exceed → `tooManyAlleles` no-call). **Impure (off-lattice) major peaks are kept**,
   first-class, treated identically to pure (no impurity discard — §3/§4 below).

> **Shared with Phase 2 (Q-G1 / Q-P1).** The `pool → rungs → clear-maxima` machinery
> is *exactly* what Phase 2 uses to **resolve confident genotypes** (homozygotes ∪
> well-separated hets — 1..ploidy clear peaks; CG-seed, parameters §2): the same peak finder,
> Phase 2 additionally testing peak count + separation + cohort-recurrence to admit a clean
> seed. **Proposed:** a shared `rungs.rs` both call. This doc and
> [ssr_call_parameters.md](ssr_call_parameters.md) must agree where it lives.

---

## 3. S2 — stutter reachability, per-allele (spec §7)

Stutter is a **per-allele** operation, not membership of a global lattice:

> **`B = A ⊕ k`** — `B` is `A` after `k` whole motif units added/removed **within
> `A`'s own repeat context, interruption(s) preserved** (the units are removed/added but
> the interrupting bases stay). For an impure `A` this is a **set** of placement variants
> (the unit can land in different runs); §6's `Σ_v` marginalizes over them.

Each candidate (pure or impure) anchors **its own stutter ladder**; a pure allele's
coincides with the global motif lattice (an emergent special case, not the
architecture), an impure allele's runs parallel carrying its interruption. The
kernel's whole-unit step `δ` is the **motif-count delta along the allele's own
ladder**, well-defined for impure alleles. The catalog **motif** is the slip unit
(the §2 admission filter already established it dominates). **Impure alleles get no
penalty** — the empirical candidate set + HMM + recurrence already filter spurious
ones; an extra penalty would bias against genuine interrupted alleles. **Slip-site
placement is NOT committed here** — S2 only declares `B` is `k` units from `A`;
*which copy* the unit lands on is **marginalized in §6 by an explicit sum over the
placement-distinct variants** of `A ⊕ k` (the runs separated by interruptions;
≈ (#interruptions + 1) terms, uniform position prior — HipSTR's run-collapsed `Σ_v`,
verify-fix #3), so the data picks it. For a **pure** allele every placement gives the
same sequence ⇒ a **single** variant ⇒ no extra cost.

This is the relocated "on/off-ladder" relationship (model doc §6) — done per-allele
from the start (the load-bearing choice: strict pure-only behaviour is then a *policy*
over the same machinery, not a re-architecture).

---

## 4. S3 — the read likelihood `Qᵣ` (spec §6, HipSTR-informed)

The **sum-over-slips** generative form, adapted to Mark-2's collapsed-sequence
evidence (distinct sequences + counts, uniform quality):

> `Qᵣ(obs_seq | candidate) = Σ_Δ  S_θ(Δ) × align(obs_seq | candidate ⊕ Δ)`
> where `align(obs | cand ⊕ Δ) = Σ_{v ∈ placements(cand,Δ)} Pr(v) · align_subst(obs | v)`

- **Placement marginalization (`Σ_v`, verify-fix #3).** `cand ⊕ Δ` is a **set** of
  placement-distinct sequences — the `|Δ|` units added/removed land in different repeat
  **runs** of the tract. The slip placement is marginalized by an **explicit sum** over
  that set with a **uniform position prior** `Pr(v)` (HipSTR's `1/(block_len+Δ+1)`).
  Crucially the set has **≈ (#interruptions + 1) members, not (#copies)**: placing the
  change anywhere inside one homogeneous run yields the *same* string, so equal-LL
  placements collapse (HipSTR's `upstream_match_lengths_` run-length trick,
  `StutterAlignerClass.cpp:84-87/133-136`). A **pure** allele ⇒ one run ⇒ **one** term ⇒
  the closed-form / fast path below applies directly (majority case, no cost); a
  singly-interrupted impure allele ⇒ **~2–3** terms.
- **`align_subst(obs | v)`** (one fixed placement variant `v`) = a **banded pair-HMM
  forward** with **flat (uniform-quality) emission** under the per-base error `ε` — a
  *probability* `P(obs | v)` (not a best-path score), summing over all ways `obs` arises
  from `v`. **Inside the repeat tract it admits no sub-motif indel** — the only
  length-changing operations on the tract are the whole-motif slip size (`Σ_Δ`) and the
  whole-motif placement (`Σ_v`), both whole-unit; *within a fixed `v`* the forward scores
  **substitutions only**, and per-base gaps are confined to the **flanks** (HipSTR's
  repeat-block rule — stutter and sequencing indels never compete on the same bases:
  `HipSTR/src/SeqAlignment/{HapAligner,StutterAlignerClass}.cpp`). This is what keeps
  `ε` (composition) and stutter (length) identifiable at every period, period 1
  included; **mononucleotide loci are flagged** (their `ε` is substitutions-only, the
  stutter level absorbs the indistinguishable single-base slip — spec §6). An
  **exact-match fast path** (`obs == v` byte-for-byte — the clean post-gate majority)
  returns `(1−ε)^len` without running the HMM. Affine-gap best-path was **rejected** (a
  score, not a probability; commits to one alignment — and would collapse the `Σ_v`). Reuses the **Stage-1 SSR
  pair-HMM** machinery (banded, scratch-buffered — pattern from BAQ
  [probaln.rs](../../../src/baq/probaln.rs)), with flat emission (Stage 1 dropped base
  qualities).
- **`S_θ(Δ)`** = the §5.2 kernel = the per-locus **shape** `θ_locus` (refined under its
  per-`(group, period)` shape prior — shrunk to a cohort-per-period parent, Phase 2, M3) **× the per-sample-group stutter level**
  (linear in repeat length `level_baseline + level_slope·length`; pre-pass-seeded and
  **refined in the outer loop**, *not* frozen — C2 amend., §5); the **same for pure and
  impure** alleles at a locus. The level makes `S_θ` *per-group*, so the re-weight
  (below) is per-group — but it is still just arithmetic over the shared align table, no
  HMM, whether the level is held fixed (within a sweep) or updated (at the barrier).

**The `align` cache is deferred — v1 recomputes, then we measure (Q-G3, resolved
2026-06-22).** With `ε` **frozen per sample group** (Phase 2), `align(obs | cand ⊕ Δ)` is
**invariant across EM iterations and outer rounds** (only `S_θ(Δ)` re-weights it), so it is
*cacheable in principle* — but a pre-materialized `(group/ε, obs, cand, slip)` align table
is an **unmeasured optimization that v1 does not build.** v1 **recomputes `align` on
demand**, for three reasons:
- **It's plausibly cheap, post-C1.** The in-tract part (per fixed placement variant) is a
  **substitution closed form** (`(1−ε)^match · (ε/3)^mismatch`) with an **exact-match fast
  path** for the clean post-gate majority; the banded pair-HMM fires only on the
  flank/impure minority. A **pure** allele is one such evaluation; an **impure** allele is
  the small **`Σ_v` placement sum** (~2–3 run-collapsed closed-form evals, verify-fix #3) —
  so the per-`(obs, cand⊕Δ)` recompute is a handful of multiplies, not a DP sweep —
  *whether this is actually the bottleneck is the thing to measure, not assume* (spec §6
  "the speed win is unmeasured"; the impure `Σ_v` is the new thing to watch on
  impurity-rich loci).
- **It dissolves the M2 cache-blowup question.** No table ⇒ no "key by per-sample vs
  per-group `ε`" decision, no `≤~8–10` vs distinct-`ε`≈`N` contradiction (former arch
  Q-G7 / Q-P6), no cohort-scaling-memory risk. Sample grouping reverts to a pure
  **estimation + reporting** device with no cache role at all.
- **Determinism is unaffected — actually simpler.** `align` is a **pure function of
  `(obs, cand⊕Δ, ε)`**; recomputing it yields **identical bits** regardless of thread or
  order, with **no per-thread cache state** to reason about. The `Σ_v` placement sum is
  itself a pure, order-fixed (catalog-run order) reduction, so it adds no nondeterminism.
  The EM re-weight stays
  `Qr_g[obs][cand] = Σ_slip Sθ_locus,g[slip] · align(obs | cand⊕slip; ε_g)` — same arithmetic,
  the `align(·)` term (its own `Σ_v` for impure alleles) recomputed (or memoized, below)
  rather than read from a table.

**If measurement shows `align` is the wall**, the *first* optimization is a **per-locus
memo** — materialize `align` lazily during one locus's EM and discard it when the locus is
done (bounded by a single locus's working set, **no cross-sample `ε`-keying**, so M2 never
arises). A persistent cross-`(group, ε)` cache is a *later* step only if even the per-locus
memo is insufficient — and *then* the per-group-`ε` quantization (distinct `ε` = #groups ≤
cap) bounds it. Per-sample `ε` remains the documented modeling upgrade, independent of
caching.

**The genotype likelihood — two FP-defense terms (spec §6, both v1 structural):**

- **uniform outlier `λ`** — `P(read | G) = (1−λ)·[mix over G's alleles] + λ·(1/D)`,
  `D` = distinct sequences at the locus. Absorbs *random* junk (normalized, data-
  defined — no `4^L` floor `λ` could never outweigh); its flatness is exactly why it
  catches random but not *systematic* signal. `λ` a conservative fixed small value in
  v1 (simulator-calibrated; estimable-with-a-cap later).
- **allele-balance / overdispersion** — stops *systematic concentrated* minor signal
  (paralog/contaminant/under-modelled stutter at a consistent length) from inflating a
  false het whose confidence **grows with depth** (the i.i.d.-likelihood failure mode
  already seen on the SNP path). Runs on the **deconvolved per-allele responsibilities**
  the E-step produces (post-`S_θ`, where a real het sits near 1/ploidy) — **not** the
  raw length tally (stutter legitimately unbalances raw lengths). Mirrors the SNP
  caller's [`qual_refine.rs`](../../../src/vcf/qual_refine.rs); feeds **both GQ and
  site QUAL** → a depth-inflated false het → low GQ → no-call. Overdispersion
  (Dirichlet-multinomial, dataset-wide concentration, pre-pass-estimated like `ε`) is
  the documented upgrade.

---

## 5. The EM — per-locus loop + the prior-side `F` outer loop (spec §3, §4.2, §4.4)

**The per-locus EM** builds on the `posterior_engine.rs` §5.4 loop:

| piece | source | Mark-2 status |
|---|---|---|
| E-step (prior × likelihood → posterior) | engine, **reused as-is** | responsibilities; alignments constant **per sample group** (per-group frozen `ε`), re-weighted by **per-group** `S_θ(Δ)` (per-group stutter level, ×length) each iter |
| IBD-mixture genotype prior `F·π_i + (1−F)·π_i^ploidy` | engine, **reused as-is** | allele-agnostic; feed it repeat-allele *sequences*, and `F_i` from the outer loop |
| `π` M-step **base measure** | engine, **replaced** | SNP class pseudocounts ([`classify_allele`](../../../src/var_calling/posterior_engine.rs)) → the per-candidate geometric `G₀` vector (Phase 2 produces it), **floored at a tiny `> 0`** so `p^|Δ|` can't underflow to hard zero (verify-fix #4 — the `π` analogue of `F_CEILING`) |
| **`θ_locus` M-step** | **new code** | regularized stutter-shape update, shrunk to the sample's `(group, period)` shape prior (M3); re-weights `align` (no rebuild — recomputed, cache deferred) |
| seed + convergence | engine §5.4, **reused** | `π⁰`/`θ⁰` seed (Phase 2); iterate to non-decreasing penalised log-lik |
| SNP-only machinery (compound alleles, chain anchors, contamination) | — | **bypassed** |

**The outer loop (spec §4.4) — `F_i` *and* the per-group stutter level.** Both live
**outside `align`** (the level is an `S_θ` re-weight, `F` is in the prior), so
re-estimating either is a cheap deterministic reduce that **never rebuilds the cache**.
The level joins the loop because freezing it made the masquerading-het contamination
permanent (C2 amend., spec §4.4) — refining it from the cohort's soft responsibilities
lets the whole cohort overrule a biased pre-pass level:

```
 F_i⁰ = supplied/default;  level_g⁰ = pre-pass seed (per group)
 repeat {
   per-locus EM over all loci (current F_i, current level_g)   # refine π AND θ_locus; align cached
   F_i  reduce over VARIABLE loci (≥2 alleles):                # mean autozygous-branch responsibility
     raw F_i  →  shrink to cohort mean  →  clamp ≤ user cap  →  clamp ≤ F_CEILING=0.99
   level_g reduce per SAMPLE GROUP (soft per-allele resp.):    # re-fit level_baseline+level_slope·length
     fixed-point integer accumulation (order-independent) → re-weight S_θ (NO align rebuild)
 } until |ΔF|,|Δlevel| < tol  or  max rounds
 genotype calls = the FINAL per-locus E-step
```

Both reduces sum per-(sample, locus) responsibilities, and loci complete **out of order**
across threads — so a per-thread float partial sum would be thread-count-dependent
(floating-point addition is not associative). The fix is **fixed-point integer
accumulation** (M1 resolution): each contribution is scaled (`× 2⁴⁰`, say), rounded to an
integer, and summed into per-individual (`F_i`) / per-group (`level`) `i128` accumulators;
integer addition is associative + commutative, so the sum is **identical regardless of
thread count or completion order**, with no buffering and `O(N_samples)` memory. This is
the *same* trick the pre-pass already uses — its `SlipProfile` / `SampleStutterStats`
sufficient statistics are `u64` counts precisely so their reduce is order-free. The final
divide (and any cohort-mean shrinkage over the fixed sample set) is then a deterministic
scalar step.

- **per-individual `F_i`**, shrunk toward the cohort mean (hierarchical; sparsely-typed
  individuals shrink harder) — captures structured/mixed-mating cohorts and localizes a
  bad sample's effect. The engine **already takes a per-sample `F` vector**
  (`log_f_per_sample` in [posterior_engine.rs](../../../src/var_calling/posterior_engine.rs)),
  so no prior-side change.
- **estimator** = mean posterior responsibility of the **autozygous branch** (a het
  contributes 0; a hom for `i` contributes `F·π_i/(F·π_i+(1−F)·π_i^ploidy)`); **only
  variable loci** enter.
- **hard ceiling `F_CEILING = 0.99`** (no `F=1` absorbing trap — the `F`-analog of the
  `π` pseudocount floor) + optional lower CLI cap.
- it is **apparent `F_IS`** (M5, spec §4.4) — absorbs cohort-wide hom-excess, chiefly
  **Wahlund/population structure** (real, **reported with a user warning, not corrected** —
  no clean within-caller fix); per-individual `F_i` localizes only an *idiosyncratic* bad
  sample, not a cohort-wide signal. **Null alleles** are a capillary-microsatellite
  (primer-dropout) confounder that **largely doesn't apply to primer-free sequencing**, so
  the per-locus null mechanism stays **deferred** (appropriate for WGS; capture/RAD may
  revisit — spec §9). The `0.99` ceiling + variable-only reduce remain as guards.

> **The engine-reuse boundary (Open Q-G2, spec §9 cross-cutting).** Reuse
> `posterior_engine.rs`'s E-step + IBD-`F` prior; **replace** its class pseudocounts
> with the per-candidate `G₀` vector by **generalizing the engine to accept a
> pseudocount vector** (the psp-container precedent: extract the reusable core, write
> the thin SSR part — don't contort SSR onto the SNP class scheme); **bypass** the
> SNP-only machinery. **Reuse the exact-AF convolution *kernel* (`convolve_ac_linear` +
> Beta-Binomial-`K`) for QUAL** (§6, m4) — call it per-candidate-allele (each collapsed as
> "reference") and sum. This is a **generalization of `compute_qual_via_exact_af`, not a
> call to it unchanged** (verify-fix #7b): that function hard-codes the collapse on
> **allele-0** (Step-1 bucket `ploidy − count[0]`; `α_ref = pseudocounts[0]`), so SSR must
> **parametrize the collapse allele** (`count[i]`), set the per-`i` Beta from `G₀`
> (`α_ref = G₀[i]`, `α_alt = Σ_{j≠i} G₀[j]`), and **bypass the REF-anchored `P(K=0)`
> wrapper** (REF = allele-0 is meaningless for an STR frame). The kernel and the `K=0` math
> are reusable; the allele-0 indexing is not. **Open: extend the shared engine in place vs a slim SSR fork of
> the π-loop** — a struct-shape call this doc owns. The math boundary is fixed either
> way; the regression gate is the SNP caller's end-to-end tests (pre-alpha, no
> byte-identity owed across format versions).

---

## 6. The CALL — cohort VCF (spec §4.5)

The GangSTR/TRtools-compatible **shell reuses** ([src/vcf/](../../../src/vcf/)); three
semantics are **redefined for SSR**:

- **Reused as-is:** REF/ALT = actual tract sequences; per-ALT **REPCN** (`len/motif`,
  fractional for impure) + **BPDIFFS** (`len − ref_len`); per-sample **GT**; **allele
  pruning** (drop ALTs no call uses — reuse `prune_unsupported_alleles`); TRtools/
  GangSTR header detection.
- **Emit iff variable, not "non-ref."** Drop iff **monomorphic** (one allele
  cohort-wide); emit iff **≥2 alleles segregate**, *regardless of frequency* — a
  variable-but-rare locus is real and emitted. "Polymorphic" (MAF ≥ 0.05) is an **AF
  annotation**, not an emit/QUAL gate.
- **Site QUAL = `−10·log₁₀ P(locus monomorphic)`** — confidence ≥2 *real* alleles
  segregate. **Estimator (m4 resolved 2026-06-22; scope amended verify-fix #7b):**
  `P(monomorphic) = Σ_{i ∈ A_ℓ} P(fixed-for-i)`, the **mutually-exclusive** per-allele
  "fixed" events summed; each `P(fixed-for-i)` reuses the engine's **exact-AF convolution
  kernel** with allele `i` collapsed as "reference" (so `P(fixed-for-i) = P(K=0)` over non-`i`
  copies). **Generalizes `compute_qual_via_exact_af`, not a call to it unchanged**
  (verify-fix #7b): parametrize the allele-0 collapse (`count[i]`), per-`i` Beta from `G₀`,
  bypass the REF-anchored `P(K=0)` wrapper (Q-G2). **Honest scope:** each term uses its own
  binary-collapse normalizer, so the sum is **exact for biallelic / a per-collapse-normalizer
  approximation for ≥3 alleles** (accurate at the extremes, loose only mid-range — fine for a
  QUAL score); the single-joint-`Z` exact form (unnormalized `K=0` numerators ÷ one joint
  marginal) is the documented upgrade. *Not* reference-anchored. The §4 FP defenses feed it
  directly: a stutter/artifact false second allele → its `π → 0` → mass on
  `fixed-for-the-true-allele` → high `P(monomorphic)` → **low QUAL**. (Cheaper proxy: naive
  `Σ_i Π_s P(GT_s=hom_i)` — perf fallback only.) (More than HipSTR's `QUAL="."`, justified because we are cohort-joint and
  already make the variable/monomorphic decision.)
- **Site FILTER** = `PASS` + SSR locus reasons (`notPeriodic`, `tooManyAlleles`,
  `lowDepth`); a filtered locus is **emitted with its reason, never silently dropped**.
  SNP `segdup` **dropped** (SSR delegates mappability to per-read MAPQ).
- **Per-sample** = GT + **GQ/`Q`** (the genotype posterior carrying the §4 allele-
  balance penalty) + DP/REPCN; **`./.`** when absent (§4.1), outlier-dominated (`λ`
  carries the mass), or sub-gate (posterior/depth/support, §5.8).
- **`F_IS` carries a user warning (M5, spec §4.5).** Wherever `F_i` / cohort `F` is surfaced
  (header comment, field, or side report), label it **apparent `F_IS`** — it **includes
  population structure (Wahlund)** + residual mapping artifacts, **not** pure inbreeding; v1
  does not decompose them (no clean correction), so a naive user is told structure is folded
  in. (`vcf_out.rs` owns emitting this label/warning.)

> **Q-G4** — the `P(monomorphic)` estimator is **resolved** (m4: sum-over-alleles exact-AF
> convolution *kernel*, generalized off allele-0, exact-biallelic / approx-≥3 with the
> single-joint-`Z` exact form deferred — verify-fix #7b, above). *Still open:* the `lowDepth`
> cohort-wide FILTER threshold and whether to also emit a per-locus polymorphism flag
> (annotation) — both calibration, not the QUAL math (spec §9).

---

## 7. Parallelism, memory, determinism

- **Parallelism** — embarrassingly parallel **across loci** (each `CohortLocus` is an
  independent EM); one EM per worker. The outer loop (§5) is the only cross-locus
  coupling — global **reduces between sweeps** (`F_i` *and* the per-group stutter level,
  C2), not per-locus dependencies. The pre-pass froze `ε` and seeded stutter
  shape/level/`G₀` (Phase 2); the level is then refined in the outer loop, so genotyping's
  only cross-locus state is those two between-sweep reduces.
- **Memory** — per locus, v1 holds only the observed sequences + per-locus EM scratch;
  there is **no persistent `align` table** (deferred, §4), so the per-locus hot allocation
  is just the `(obs × candidate × slip)` *scratch* of whatever is being scored, reused via
  scratch buffers. The hypervariable-locus compute-for-memory trade and any later
  `align`-cache memory are an explicit **measure-first** item (spec §9), not a v1 cost.
  Cohort working set bounded by Phase 1's lockstep rule (≈ N × one block).
- **Determinism** — per-sample-group frozen `ε` ⇒ `align(obs|cand⊕Δ)` is a pure function;
  recomputing it (no cache) yields identical bits with **no per-thread cache state** to
  reason about; the per-iteration motion (per-group `S_θ` re-weight, `π`/`θ_locus` updates)
  is **per-locus/per-group, not per-thread** — a pure function of the locus's data + the
  global priors + the frozen `ε`, computed identically on whichever worker owns it. The two **between-sweep** reduces (`F_i` and the per-group level, C2) are
  the only cross-locus float sums, and loci complete out of order — so they use
  **fixed-point integer accumulation** (scale → round → sum into `i128`; the pre-pass's
  `u64`-counts trick), which is associative + commutative ⇒ identical regardless of thread
  count or completion order (M1 resolution, §5). With that, Stage 2 is **byte-identical
  across and within thread counts**, carrying the Stage-1 invariant forward (spec §4.4);
  the loop's `|ΔF|,|Δlevel| < tol` + max-rounds cap is a fixed, order-independent stop.

---

## 8. Proposed module layout (sketch)

```
src/ssr/cohort/
  rungs.rs        # SHARED w/ Phase 2 (Q-G1): pool -> rungs -> clear-maxima
  candidate_set.rs# S1: nominate + union + locus-admission filter -> CandidateSet  (§2)
  stutter.rs      # SHARED w/ Phase 2: S_θ(Δ) kernel + reachability ⊕ (S2, §3): emits the
                  #   placement-variant SET for an impure A⊕Δ (runs between interruptions) + θ_locus M-step
  likelihood.rs   # S3: align = Σ_v placement sum over run-collapsed variants (uniform prior, ~2-3
                  #   for impure, 1 for pure — verify-fix #3); per variant: in-tract subst. closed-form
                  #   + exact-match fast path, banded HMM for flanks/impure; + Σ_Δ re-weight + λ outlier
                  #   + allele-balance term  (§4; align recomputed on demand — cache deferred)
  pair_hmm.rs     # the banded forward kernel align_subst(obs|v) for ONE fixed placement variant
                  #   (reuse Stage-1 SSR pair-HMM; flat emission; in-tract = substitutions only,
                  #   gaps in flanks only — §4)
  em.rs           # per-locus EM: drives the engine E-step + G₀ M-step + θ_locus M-step  (§5)
  prior.rs        # thin adapter onto posterior_engine IBD-F prior; per-candidate G₀ vector inject (§5, Q-G2)
  inbreeding.rs   # the outer loop: F_i reduce (variable-loci + shrink + clamp) AND the
                  #   per-group stutter-level reduce (soft resp.; fixed-point i128 accum,
                  #   order-independent → byte-identical across threads)  (§5, C2, M1)
  vcf_out.rs      # the CALL: variable-emit, P(monomorphic) QUAL, SSR FILTER/GQ  (§6)
```

(Names extend the [ssr_genotyping_architecture.md §4](ssr_genotyping_architecture.md)
`cohort/` sketch; `likelihood.rs`/`inbreeding.rs` are Mark-2 additions for the
sum-over-slips HMM and the `F` loop.)

---

## 9. Open items (Phase-3 agenda)

- **Q-G1** — shared `rungs.rs` boundary with Phase 2 (§2; ties to Q-P1).
- **Q-G2** — engine reuse: **extend in place vs slim SSR fork** of the π-loop (§5);
  the per-candidate-pseudocount-vector generalization of `posterior_engine.rs`.
- **Q-G3** (**resolved 2026-06-22 — defer the cache**) — the pre-materialized `align`
  table is an **unmeasured optimization v1 does not build**: v1 **recomputes `align` on
  demand** (in-tract = substitution closed-form + exact-match fast path, post-C1; banded
  HMM only for flanks/impure), which keeps `align` a pure function of `(obs, cand⊕Δ, ε)`
  (determinism, no per-thread cache) and **dissolves the M2 cache-keying question
  entirely** (§4). *Measure first*; if `align` is the wall, add a **per-locus memo** (lazy,
  discarded at locus end — no cross-sample `ε`-key, so still no M2), and only then a
  persistent cross-`(group, ε)` cache (bounded by per-group-`ε` quantization: distinct `ε`
  = #groups ≤ cap; a cache keys on the placement variant `v`, not just `cand⊕Δ`). The
  in-tract band carries **substitutions only** (no sub-motif indel states — C1, §4); within
  a fixed placement variant a fixed `Δ` aligns same-length sequences, so the recompute is
  cheap and the period-1 flag rides here.
- **Q-G3b** (**resolved 2026-06-22 — verify-fix #3**) — the slip-placement marginalization
  is an **explicit `Σ_v` sum** over the placement-distinct variants of `cand ⊕ Δ` (the runs
  between interruptions, uniform position prior, equal-LL placements collapsed —
  HipSTR's `upstream_match_lengths_`), **not** a single committed sequence. Pure alleles ⇒
  1 term (no cost); a singly-interrupted impure allele ⇒ ~2–3 terms. This restores the
  faithful HipSTR borrow (a single substitution-only forward against one `cand⊕Δ` could
  *not* marginalize placement for impure alleles). `stutter.rs` enumerates the variant set;
  `likelihood.rs` sums. *Open (calibration):* whether multi-interruption / compound impure
  alleles need a richer placement enumeration than "one run per interruption," and whether
  the uniform position prior should ever be informed — both deferred to the simulator.
- **Q-G7** (**superseded 2026-06-22 by Q-G3**) — this used to ask whether the `align` cache
  is keyed per-sample vs per-group-`ε` (and claimed shrinkage keeps distinct `ε` small "on
  its own" — the M2 flaw: shrinkage only collapses *thin* samples, so a deep continuum
  cohort kept ≈`N` distinct `ε` and blew the cache). **With the cache deferred (Q-G3) the
  question is moot:** no table ⇒ no `ε`-keying decision. Sample grouping is now a pure
  **estimation + reporting** device with **no cache role**, and the discrete-vs-continuum
  question is validation only. The per-group stutter level enters only the `S_θ` re-weight
  and is **refined in the outer loop** (C2); the shape is per `(group, period)` shrunk to a
  cohort-per-period parent (M3). If a cache is ever added (Q-G3), per-group-`ε` quantization
  — not shrinkage — is what bounds it.
- **Q-G4** — `P(monomorphic)` estimator **resolved** (m4: sum-over-alleles exact-AF
  convolution *kernel*, generalized off allele-0; exact-biallelic / approx-≥3, single-joint-`Z`
  exact form deferred — verify-fix #7b, §6); *still open:* `lowDepth` cohort-wide FILTER
  threshold + optional polymorphism flag (calibration — §6; spec §9).
- **Q-G5** (**resolved 2026-06-21**, arch Q-P3) — the per-locus seed (`π⁰`/`θ⁰`) is
  computed **here in Phase 3**, fused with candidate assembly's `pool → rungs` pass — *not*
  in the pre-pass, which runs on a subset and so can't seed every locus.
- **Q-G6** (spec §9 details) — rung recurrence *k*, prominence floor, same-length
  frequency threshold, `MAX_CANDIDATE_ALLELES`; admission-filter robustness (empty
  rungs / motif multiples / min rung count); `λ` value + allele-balance form/strength +
  simulator calibration; `θ_locus`←`(group,period)`←`θ_period` shrinkage strengths (M3).

**Settled upstream (spec §5/§6/§7/§4.4/§4.5), not re-opened here:** empirical
candidates (no synthesis, read-length-capped); per-allele reachability, no impurity
penalty, **placement marginalized by an explicit `Σ_v` run-collapsed sum** (verify-fix #3 —
pure ⇒ 1 term, impure ⇒ ~2–3, uniform position prior); flat-`ε` forward (per fixed variant)
+ exact-match fast path, affine-gap rejected; **in-tract substitutions-only (no sub-motif
gaps; length change is the `Σ_Δ` slip + `Σ_v` placement, both whole-unit — HipSTR's
repeat-block rule), mononucleotides flagged** (spec §6);
per-sample-group-frozen `ε` ⇒ `align` is a pure iteration-invariant function ⇒
deterministic (cache **deferred** — recompute, M2); per-group stutter level (linear in
length) as an `S_θ` re-weight **refined in the outer loop alongside `F` — not frozen** (C2
amend., spec §4.4); stutter shape **per `(group, period)` shrunk to a cohort-per-period
parent** (M3 — not cohort-shared); the pre-pass estimator is the **soft full-cohort EM
reduce** (confident-**genotype** gate = seed — homozygotes ∪ well-separated hets, CG-seed
2026-06-22; C2); `λ` + allele-balance as
the two FP terms; per-individual `F_i` with `0.99` ceiling; emit-variable + QUAL =
Phred(variable) + SSR FILTER reasons.
