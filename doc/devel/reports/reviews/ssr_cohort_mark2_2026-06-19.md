# Adversarial review — `ssr_cohort_mark2.md` (Stage-2 `ssr-call` model spec)

**Status:** ✅ **all 13 findings resolved in the spec** (2026-06-19, same day) — see the
resolution table below; the per-issue write-ups that follow are the original review, kept
as the rationale trail. Branch `ssr-cohort`. Target:
[doc/devel/specs/ssr_cohort_mark2.md](../../specs/ssr_cohort_mark2.md). Companion
docs read for grounding: [ssr_genotyping.md](../../specs/ssr_genotyping.md) (Mark-1,
the reuse targets in §5), [ssr_ladder_model.md](../../architecture/ssr_ladder_model.md),
[ssr_pileup_mark2.md](../../architecture/ssr_pileup_mark2.md) (as-built Stage 1),
[ssr_offladder.md](../../architecture/ssr_offladder.md). Code read:
[src/var_calling/posterior_engine.rs](../../../../src/var_calling/posterior_engine.rs)
(the IBD-prior reuse target). Reference implementations mined:
`HipSTR/`, `GangSTR/` (vendored under the main repo).

This is a skeptic's pass, not a sign-off. The spec is strong on the parts it treats
as *new* (empirical candidate assembly, the sum-over-slips likelihood, the
cacheability argument). The defects cluster in two places: **things waved through as
"reuse" that the reuse target does not provide**, and **the EM control structure** —
claimed "settled end-to-end" but never specified and not supplied by the reuse. The
"HipSTR-style" provenance is also overstated.

The through-line: tighten the boundary between *reuse* and *new work*. Findings 1–5
all stem from crediting `posterior_engine.rs` / Mark-1 §5.4 / §5.9 with capabilities
they don't have (cross-locus pooling, F estimation, λ, an SSR-appropriate QUAL).

Each issue below is self-contained, with the exact section, the evidence (quoting the
spec and the contradicting doc/code/tool), and a proposed fix — so we can take them
one at a time. Severity: **BLOCKER** (must resolve before "settled"), **MAJOR**
(genuine defect / large gap), **MINOR** (wording / smaller gap). Each is tagged
**defect** (objectively wrong/missing) or **reconsider** (defensible but worth a
second look).

---

## Resolution status (all closed 2026-06-19)

Discussed one-by-one with the author; each landed in the spec. Where the resolution
diverged from the review's first proposal, that's noted — the discussion improved several.

| # | sev | finding (short) | resolution | spec § |
|---|---|---|---|---|
| 1 | BLOCKER | EM schedule undefined; "reuse §5.4 unchanged" false | Two-level schedule written: pre-pass → per-locus EM → `F` outer loop. Reuse boundary made explicit: E-step + IBD prior reused; base measure **replaced**; `θ_locus` M-step **added** | §4.4; §3 (§5.4 row); §4.2 spine |
| 2 | MAJOR | per-thread `ε` breaks cache-once + determinism | `ε` **frozen** (pre-pass); per-thread/`δ`-rebuild retired; Stage-2 now byte-identical across/within threads. *(Diverged: chose freeze over the global-`ε`-with-lock first floated.)* | §4.4; §6; §4.2 (`ε` row) |
| 3 | MAJOR | `F` "estimated" but reuse target can't | **Per-individual `F_i`**, prior-side outer loop, EM-responsibility estimator, `0.99` ceiling, apparent-`F_IS` caveat. *(Diverged: per-individual, not global.)* | §4.4; §4.2 (`F` row); §9 |
| 4 | MAJOR | output/QUAL/FILTER SNP-shaped | New §4.5: **emit-variable / drop-monomorphic**, **QUAL = Phred(locus variable)**, SSR FILTER (`notPeriodic`/`tooManyAlleles`, drop `segdup`), per-sample no-call | §4.5; §3 (§5.9 row) |
| 5 | MAJOR | `θ⁰` seed dominance gate dropped (worse once `θ` frozen) | `θ` **demoted to a per-cell prior, refined per-locus** → seed contamination low-stakes; selection downgraded to prior-hygiene. *(Diverged: solved by `θ`-in-EM, not robust-selection machinery.)* | §4.4; §4.3 (`θ⁰`); §4.2 (stutter row) |
| 6 | MAJOR | "HipSTR-style" mischaracterized | Relabeled **"HipSTR-informed"** + provenance note: genuinely-HipSTR vs our simplifications (GangSTR in-frame stutter; out-of-frame folded; flat-emission forward ours) | §6; §1 |
| 7 | MAJOR | depth-driven het inflation; λ doesn't fix | **Allele-balance / overdispersion term** (mirrors SNP `qual_refine.rs`) → GQ + site QUAL; λ scoped to *random* junk only | §6; §4.2; §3 (§5.8 row); §9 |
| 8 | MAJOR | reference-centred `G₀` re-imports ref bias | `G₀` **centred on per-locus cohort modal allele**; reference = coordinate frame only | §4.3; §3 (§5.5 row); §9 |
| 9 | MINOR | λ "junk" support undefined | **junk = uniform over the locus's `D` distinct observed sequences (`1/D`)** | §6; §9 |
| 10 | MINOR | "computed once" wording | **Moot** — with `ε` frozen, "alignments computed once for the whole run" is literally true | §1; §6 |
| 11 | reconsider | spanning-only is structural, restate | Restated in §5 as a **structural callable-range cap** (no synthesis); read-pair merging flagged as the lever | §5 |
| 12 | reconsider | SNP-phasing excluded but solves merged-het | One-paragraph note: excluded (sparsity + caller coupling), merged-het resolved via cohort recurrence + per-locus `θ` instead | §4.2 (identifiability) |
| 13 | reconsider | stop attributing joint-EM to HipSTR | Folded into #6's provenance note; the pre-pass freeze is now the genuinely-HipSTR-faithful piece | §6 |

**Net:** the per-thread-`ε`/`δ`-rebuild design is gone; the EM schedule and engine-reuse
boundary are explicit; `ε`/`θ`/`F` estimation is fully specified (frozen `ε`,
per-cell-prior-refined `θ`, prior-side per-individual `F_i`); FP control is a clean
three-layer story (λ random / recurrence cohort / allele-balance systematic); the output
is defined; `G₀` is cohort-mode-centred; HipSTR provenance is honest. The remaining
*open* items (now genuinely details, not gaps) live in spec §9 — thresholds,
calibrations, shrinkage strengths, and the extend-vs-fork engine-shape call.

---

## Lead findings (most threaten correctness)

### Issue 1 — BLOCKER / defect — The hierarchical EM schedule is undefined, and "reuse §5.4 topology unchanged" is false

**Sections:** §3 (reuse table, "§5.4 EM topology … **unchanged in topology**"),
§4.2 ("The EM spine (reused §5.4 topology)"), §4.2 ingredient table.

**Problem.** The spec assigns each parameter a *shard* — π per-locus, θ per
period×length cell (pooled across loci+samples), ε per-thread (pooled across loci),
F global, λ — but never defines the **loop nesting**: the order of M-steps, how the
per-locus π E/M interleaves with the cross-locus θ/ε/F reduces, how many outer
rounds, and what convergence means across the hierarchy. This is *model intent*
(how parameters at different shards interleave), not struct shape, so "deferred to
the architecture doc" doesn't cover it — and it is absent from the §9 open agenda,
in a doc whose commit message says "model settled end-to-end."

The reuse claim makes it worse. The cited target runs **per-locus-independent EM**:

> "Each upstream record is processed by an **independent per-record EM loop**."
> — [posterior_engine.rs:13-14](../../../../src/var_calling/posterior_engine.rs#L13)

It estimates only `p̂` and `f̂_C` *within one record*, with **no cross-record
coupling of any parameter**. The Mark-2 model needs θ, ε, F shared across loci —
exactly what the reused engine does not and cannot do. So the cross-locus machinery
is **new orchestration**, not "unchanged topology."

**Fix.** Specify the two-level schedule explicitly: inner per-locus EM over
π/genotypes given *frozen* global (θ-cell, ε, F); outer global M-steps re-estimating
θ/ε/F from accumulated per-locus responsibilities; outer convergence on the global
parameters + aggregate penalised log-lik. State that this outer loop is **new code
wrapping** the per-locus engine — the engine supplies only the inner per-locus
E/M-for-π.

---

### Issue 2 — MAJOR / defect — Per-thread ε breaks both the "computed once across the cohort" cache claim and the only determinism guarantee

**Sections:** §1 (point 1), §5 (step 3), vs §6 ("`ε` is a full EM parameter, with a
lazy cache-rebuild rule (`δ`)").

**Problem — two contradictions from one decision.**

*Cache-once is false.* §1:

> "Identical observed sequences collapse **across the whole cohort**, so the
> `|distinct seqs| × |candidates| × |slips|` alignment scores are **computed once**
> (θ-independent), reused across EM iterations."

But §6 makes ε **per worker thread**, and ε lives *inside* `align`:

> "A change in a thread's `ε` rebuilds **that thread's** cached `align` tables."

So each thread holds its *own* cache keyed on its *own* ε — the same
`(obs, cand, slip)` triple is aligned up to **once per thread**, rebuilt on every
ε-drift past δ. Not "once across the whole cohort."

*Determinism is sacrificed.* Stage 1 hard-guarantees byte-identity across thread
counts ([ssr_pileup_mark2.md §7](../../architecture/ssr_pileup_mark2.md)); §6 accepts
Stage 2 is "**not byte-identical across — or, with the dynamic work-queue, within —
thread counts**." For a popgen tool whose genotype calls feed downstream
diversity/association analysis, non-reproducible *calls* (§6 admits "the rare
borderline call" flips) is a real cost the spec under-weights, and it diverges from
the SNP path's reproducibility.

**Why the asymmetry is unjustified.** θ and F are *global / per-cell deterministic
reduces* (§4.2). ε has the same structure; it is made per-thread only to "avoid any
global synchronization." But the δ-tolerance already makes rebuilds rare — a
**global ε with a δ-gated rebuild** is one synchronized reduce per outer round
(identical cost structure to θ's per-cell reduce), making the cache *genuinely*
once-per-epoch *and* deterministic. The spec itself calls the per-thread δ-rebuild "a
deliberate, acknowledged **premature optimization**."

**Fix.** Make ε a **cohort-global EM parameter** (one deterministic reduce per outer
round, δ-gated rebuild). Resolves both contradictions at once and aligns ε with θ/F.
Demote per-thread ε to a measured optimization only if a global reduce is shown to
bottleneck.

---

### Issue 3 — MAJOR / defect — F is claimed "estimated in the EM (v1)," but the reuse target does not estimate F

**Sections:** §4.2 ingredient table ("**F** … global, **estimated in the EM** (v1) …
a cheap **deterministic global reduce**"), §9 ("**F granularity**: *Resolved* —
global, estimated"), vs §3 ("reuse `posterior_engine.rs` **verbatim** … adapter
only").

**Problem.** `posterior_engine.rs` treats the fixation index as a **fixed input**:
`fixation_index_default` + per-sample `fixation_index_overrides`
([posterior_engine.rs:228-237](../../../../src/var_calling/posterior_engine.rs#L228)),
used only to *evaluate* the prior. There is **no F M-step** — the M-steps are on `p̂`
and `f̂_C`, where `f̂_C` is the **compound-allele cohort frequency, not inbreeding**
(module doc lines 21-23; `m_step_f_hat_compound` at
[posterior_engine.rs:2890](../../../../src/var_calling/posterior_engine.rs#L2890)). A
repo-wide search for any heterozygote-deficit / F estimator returns nothing. Mark-1
§5.4 itself only ever called F-estimation "optional" and never built it. So
"estimated F via a cheap global reduce, reusing the engine verbatim" is **new,
unimplemented machinery mislabelled as reuse**.

**Fix.** Either (a) drop F-estimation from v1 and use supplied/default F — exactly
what the reuse supports, the honest v1; or (b) keep it but specify it as new code:
the heterozygote-deficit estimator, where it sits in the outer loop (Issue 1), and
that it is *not* `posterior_engine` reuse.

---

### Issue 4 — MAJOR / defect — Site QUAL and no-call/filter surfacing are undefined for SSR; "§5.9 reuse, unchanged" imports SNP-shaped output that does not fit

**Sections:** §3 / §4.1 / pipeline step 5 ("CALL … VCF (§5.9 reuse)").

**Problem.** The reused QUAL is "Phred of `Π_s P(hom-ref)_s`"
([posterior_engine.rs:720-723](../../../../src/var_calling/posterior_engine.rs#L720))
— "probability the site is variant," anchored on **hom-REF**. For an STR the
reference tract "is **not an allele claim**"
([ssr_ladder_model.md §1](../../architecture/ssr_ladder_model.md)) and at a
polymorphic locus may be carried by **zero** samples, so `Π_s P(hom-ref) → 0` and
site QUAL saturates to ~∞ at essentially every callable locus. The metric is
meaningless for SSR.

Separately, Mark-2 introduces **no-call reasons the reused FILTER vocabulary lacks**:
the §5 locus-admission filter ("emitted as a *filtered* record with a reason") for
non-periodic loci, and the `MAX_CANDIDATE_ALLELES` no-call — neither maps to §5.9's
`lowGQ/lowDepth/segdup/lowSupport`. And the λ-outlier's effect on a call (a sample
whose reads are mostly junk) has no defined surfacing. So §5.9 is *not* a free reuse.

**Fix.** Define SSR site QUAL afresh (e.g. Phred that the locus is **polymorphic in
the cohort**, or per-sample GQ only with no site QUAL); add FILTER reasons
(`nonSSR` / `tooManyAlleles` / `outlierDominated`); define how per-sample-varying
allele sets reconcile into cohort GT/REPCN/BPDIFFS. Move §5.9 from "reuse" to "open."

---

### Issue 5 — MAJOR / defect — The θ⁰ seed regresses from Mark-1: the dominance gate that protected the stutter kernel from merged-het contamination is silently dropped

**Sections:** §4.3 ("**θ⁰ — read off the confident homozygotes' skirts**" … "every
read that is not exactly `A` is a stutter product of `A`"), using §5's "clear maximum
= **> 3 reads above each adjacent rung**."

**Problem.** Mark-1 §5.4 spent its longest passage on exactly this hazard: a **skewed
adjacent-unit het** (e.g. 85/15) "can pass the dominance gate and masquerade as a
homozygote, **feeding its real minor allele into the stutter estimate**," guarded by
`SEED_MIN_DOMINANCE_FRAC ≈ 0.9` plus an optional per-cell bimodality auto-detector.
Mark-2's seed uses only a **prominence floor** (">3 reads above neighbor") with **no
dominance fraction and no masquerader analysis**. A 60/40 adjacent het easily clears
">3 reads above each adjacent rung" → mislabelled homozygote → its partner allele
contaminates θ⁰, **inflating u/d**.

The spec's "self-correcting" argument addresses only the **π⁰** mislabelling
(pseudocount keeps the masked allele alive); it says nothing about **θ⁰
contamination**, which is a distinct harm and the entire reason Mark-1 built the
dominance machinery.

**Fix.** Carry forward a dominance gate for *which* samples seed θ⁰ (separate from the
π⁰-tally rule), or restate Mark-1's `SEED_MIN_DOMINANCE_FRAC` / per-cell
auto-detection. At minimum, acknowledge the θ⁰-inflation direction (toward *more*
stutter → het-deflation) and that it compounds with Issue 7.

---

## Further defects / gaps

### Issue 6 — MAJOR / defect — "HipSTR-style" is mischaracterized; the stutter model is actually GangSTR's, and HipSTR's out-of-frame term is silently dropped

**Sections:** §1, §6 (title "HipSTR-style"), §6 "Why HipSTR cannot cache it cheaply,"
§7.

**Evidence (from HipSTR / GangSTR source).**
- HipSTR's stutter is **6-param** — a separate *in-frame* repeat-unit geometric **and**
  an *out-of-frame* bp geometric, each with up/down (`HipSTR/src/stutter_model.h:15-31`).
- HipSTR's alignment uses **per-base qualities** (`base_log_correct/wrong`,
  `HipSTR/src/SeqAlignment/HapAligner.cpp:582`) and **Viterbi/`std::max` in the
  flanks** (`HapAligner.cpp:149-153`), not flat-forward.
- HipSTR's stutter is a **separate frozen length-based EM pre-training**, not joint
  with the sequence genotypes (`genotyper_bam_processor.cpp` → `EMStutterGenotyper`,
  then the frozen model handed to `SeqStutterGenotyper`).
- The spec's "3-param geometric `(u,d,ρ)`, flat emission" is in fact **GangSTR's**
  enclosing-read stutter (`stutter_up/down/p`, `GangSTR/src/enclosing_class.cpp:64-71`).

What *is* genuinely HipSTR: slip-placement marginalized inside `align` with a uniform
position prior (confirmed, `HipSTR/src/SeqAlignment/StutterAlignerClass.cpp:65`), and
sequence-keyed first-class impure alleles. So the provenance is mixed and the headline
is wrong.

Also: Mark-1 (§5.2/§12) explicitly flagged "we fold out-of-frame into error/off-ladder
— a deliberate simplification … add a separate term if validation shows
misattribution." Mark-2 inherits this omission but **never restates the caveat**.

**Fix.** Relabel: "**HipSTR-inspired** sequence-keyed candidates + slip-placement
marginalization; **GangSTR-style** 3-param in-frame stutter; flat-error *forward* (our
cacheable simplification, coarser than HipSTR's quality-weighted Viterbi)." Restate
the dropped out-of-frame term as a known simplification.

---

### Issue 7 — MAJOR / reconsider — The SSR likelihood inherits the SNP path's depth/het-inflation blind spot, and λ does not fix it

**Sections:** §4.2 / §6 (the λ outlier rationale).

**Problem.** The known SNP failure mode (project memory
`qual_fp_depth_inflation`): the i.i.d.-over-reads likelihood over-rewards persistent
low-VAF signal, so FP QUAL inflates with depth and FPs sit at ~20% VAF. The SSR
genotype likelihood
`P(reads|G)=Π_r[(1−λ)Σ_a(1/ploidy)P(r|a)+λ·junk]` is **the same i.i.d. product**. λ
is a **uniform** term — it absorbs *random* junk but **not a systematic, recurrent
low-VAF shadow at a consistent length**, which is precisely what mints a false het and
what scales with depth. The genuine defenses are the stutter kernel (only as good as
its rate — see Issue 5) and cohort recurrence (a *cohort* defense, not a
within-sample-at-high-depth one). There is **no allele-balance / overdispersion
term**. So λ *masks* rather than *fixes* the blind spot, contrary to the §4.2 framing.

**Mitigation already present (uncredited):** Stage 1 reservoir-caps reads at
`MAX_READS_PER_LOCUS` ([ssr_pileup_mark2.md §2/P6](../../architecture/ssr_pileup_mark2.md)),
bounding per-locus depth — but **site QUAL still inflates with cohort size N**, and
the cap means §5's ">3 reads" / "≥k samples" thresholds operate in *capped-count*
space (the spec never mentions this).

**Fix.** Consider a beta-binomial / overdispersion or explicit allele-balance penalty
on the per-genotype likelihood (the same fix the SNP path needs), and state that λ
handles random junk only. Acknowledge the count-cap's interaction with §5 thresholds.

---

### Issue 8 — MAJOR / reconsider — Reference-centred geometric G₀ re-imports the reference bias that Mark-2 was created to remove

**Sections:** §4.3 / §3 ("`G₀` … **reference-centred** and **symmetric**," geometric
in `Δ=(L−ref)/period`).

**Problem.** Mark-2's premise is demoting the reference to "a coordinate frame … **no
claim about which alleles are real**"
([ssr_ladder_model.md §2](../../architecture/ssr_ladder_model.md)), because a
reference-anchored model is reference-biased when "the population sits off the
reference's whole-unit ladder." Yet G₀ **centres the prior mass at the reference** and
re-enters every M-step. For a population shifted off-reference — the exact case Mark-2
targets — a reference-centred geometric **down-weights the true modal allele** at
small N / thin loci. The one place Mark-2 keeps reference-anchoring is the prior,
undercutting its own thesis.

**Fix.** Centre G₀ on the **cohort-empirical mode** (available from the
candidate-assembly pool), not the reference — consistent with Mark-2's empirical
philosophy. If reference-centring is kept, justify why the bias removed elsewhere is
acceptable here.

---

### Issue 9 — MINOR / defect — The λ outlier "junk" term is undefined (support / normalization)

**Section:** §6 ("`λ·junk`").

**Problem.** A uniform outlier needs a *support* to be uniform over (all
byte-sequences? the observed distinct sequences? the candidate set?) and a
normalization — and that constant sets the likelihood floor at which a read is deemed
junk vs a minor-allele read. Undefined, yet it directly controls the
het-inflation/deflation trade the spec leans on. (Neither HipSTR nor GangSTR has such
a term to borrow from — this is wholly new.)

**Fix.** Define the junk support and its normalization explicitly; tie λ's
calibration to it.

---

### Issue 10 — MINOR / defect — Headline "cache computed once" wording

**Sections:** §1, §5.

Even if Issue 2's global-ε fix is taken, the cache is "once **per ε-epoch**," not
literally "computed once … reused across EM iterations" (it rebuilds on ε-drift past
δ). The §6 body is honest about this; the headlines are not. Align the wording.

---

## Defensible, but worth reconsidering

### Issue 11 — reconsider — Spanning-only is now *structural*, not a filter choice; restate it in §5

Because Mark-2 candidates are *observed* sequences with **no synthesis**, the callable
allele range is hard-capped at read length **by construction**. This is stated upstream
(Mark-1 §1.3/§1.4) but the empirical-candidate design makes it stronger, and §5 should
say so. GangSTR's FRR-Poisson-count, spanning-pair insert-size, and flanking
lower-bound classes (`GangSTR/src/frr_class.cpp:169-222`,
`spanning_class.cpp:115-159`, `flanking_class.cpp`) are the only ways past it —
**correctly out of scope** for plant SSRs (expansions are a non-goal), do *not* adopt
now. But note Mark-1's cheap in-paradigm lever (read-pair merging to extend the
spanning ceiling, Mark-1 §1.4) is silently absent from Mark-2.

### Issue 12 — reconsider — HipSTR SNP-phasing is excluded, but it solves the exact problem the spec struggles with

HipSTR's biggest het-resolution lever is phasing STRs against nearby het SNPs
(`HipSTR/src/snp_phasing_quality.cpp`, soft-assigning each read to one of the two STR
alleles). The spec excludes it (Mark-1 §1.1: SSRs too sparse, couples the callers) —
defensible. But the merged-het / het-inflation problems (Issues 5, 7) are *precisely*
what phasing resolves. Worth one sentence acknowledging the trade rather than only
"excluded by design."

### Issue 13 — reconsider — Stop calling the joint cohort EM "the HipSTR approach"

HipSTR deliberately *decouples* stutter (separate frozen pre-training) from genotyping
(per-sample MAP). The spec's joint cohort EM is arguably better (cohort recurrence buys
identifiability) but heavier and carries the θ⇄genotype circularity the seed must tame.
Justify the joint choice on its own merits; don't attribute it to HipSTR (see Issue 6).

---

## Spot-checks that came back clean (no action)

- **Input contract (§2) matches the as-built Stage 1** — `Vec<(Box<[u8]>,u32)>` +
  the exact QC scalar set (`depth, mapped_reads, n_filtered, n_low_quality,
  n_border_off_end`) is what `SsrLocusObs` writes
  ([ssr_pileup_mark2.md §5](../../architecture/ssr_pileup_mark2.md)). ✓
- **`OFFLADDER_PRIOR_FACTOR` removal is clean and justified** (§3/§4.3/§7: impure
  alleles first-class, no penalty; off-ladder enum retired in §8). No dangling
  reference. ✓
- **"No impurity penalty" is coherent with the locus-admission filter** — locus-level
  periodicity discard vs per-allele impurity-without-penalty are distinct and
  consistent (§7). ✓
- **Genotype enumeration is specified** (integer partitions of ploidy over `A_ℓ`,
  reused via `genotype_order`) — not a gap. ✓
- **Small-N / single-sample / zero-coverage** — covered (§4.3 prior-alone fallback;
  §4.1 absent = no data; Mark-1 §5.6). ✓

---

## Suggested discussion order

1. Issue 1 (EM schedule) — blocks "settled"; everything else assumes a schedule.
2. Issue 2 (per-thread ε → global ε) — resolves a contradiction *and* simplifies.
3. Issue 3 (F estimation: drop, or specify as new).
4. Issue 5 (θ⁰ seed dominance gate) — couples to Issue 7.
5. Issue 7 (depth/het-inflation; allele balance / overdispersion).
6. Issue 4 (SSR output / QUAL / FILTER).
7. Issue 8 (G₀ centring).
8. Issues 6, 13 (HipSTR provenance) + 9, 10 (λ definition, wording).
9. Issues 11, 12 (scope restatement, phasing note).
