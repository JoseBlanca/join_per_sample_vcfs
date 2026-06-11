# Adversarial review — SSR/STR genotyping specification

**Date:** 2026-06-10
**Target:** [`doc/devel/specs/ssr_genotyping.md`](../../specs/ssr_genotyping.md)
**Method:** close read of the spec plus verification of every concrete claim it
makes about *other* code — the SNP caller's genotype prior
([`src/var_calling/posterior_engine.rs`](../../../../src/var_calling/posterior_engine.rs)),
the vendored TRTools harmonizer (`TRTools/trtools/utils/tr_harmonizer.py`), and
HipSTR's stutter model (`HipSTR/src/`). Those three were checked against source;
findings are quoted with file:line below.

The spec is strong overall: internally cross-referenced, honest about deferrals,
and with a clean separation of stutter (Stage 2) from sequencing error (Stage 1).
This review is deliberately adversarial — it records only where the spec appears
wrong, circular, or claims "settled" over something that is not. Findings are
grouped by severity. Three of the four high-severity items are assertions the
spec makes about external code that verification contradicts.

Status of each item is tracked in the **Resolution** line; update it as we work
through the points and amend the spec.

---

## High severity

### H1 — Integer-unit allele model contradicts "keep imperfect loci"

**Resolution:** _resolved 2026-06-10 — hybrid allele model._ Allele identity is now
the **sequence**, not the integer count (glossary `allele (here)` + `on-ladder /
off-ladder`; §1). On-ladder alleles are stored compactly as integer rungs (sequence
reconstructible from the catalog scaffold); off-ladder alleles (1 bp indels, partial
units, variable interruptions) carry their **normalized sequence** explicitly
(§4.2, §4.3 `offl_*` columns). The integer repeat count is demoted to a derived
coordinate used only by stutter (§5.2) and the prior (§5.5). Imperfect loci keep
their fractional `ref_copies` as a derived annotation; the reference identity is the
reference tract sequence, fractional remainder pinned as fixed scaffold (§3.2). VCF
emits true REF/ALT sequences + `BPDIFFS`; copy number derived & fractional-tolerant
in TRtools (§5.9). No assembly → off-ladder sensitivity limit documented (§4.2,
§12-adjacent). Off-ladder + imperfect-locus stress added to Bucket-1 tests (§7).
Chosen per discussion: HipSTR-style "keep full information" so any downstream
analysis is possible; hybrid keeps Stage 1 efficient. _Decided over: pure-integer
(GangSTR-style, collapses 12 vs 12+1 bp) and all-sequence (drops the integer ladder
entirely — purer but heavier)._

The whole length-genotyper math assumes alleles are **integer repeat counts on a
clean motif tiling**:

- haplotypes are `H_L = outer_flank + (motif × L) + inner_context` (§4.2) — only
  integer `L`;
- candidate lengths, `hist_lengths`, `amb_lengths` are `int16` **repeat units**
  (§4.3);
- `ref_copies = (end−start)/period`, `Δ = (L − ref)/period` (§3.2, §5.5);
- REPCN and `len/motif` copy numbers in the VCF (§5.9).

But Stage 0 **deliberately keeps imperfect/interrupted single-motif loci** and
genotypes them (§3.1). For an imperfect locus, `end − start` is generally **not**
a multiple of `period`, so:

- `ref_copies` is non-integer → the `Δ`-offset grid in §5.5 is misaligned, and the
  "reference-centred" base measure has no integer centre;
- the real alleles include partial-unit and interruption variants that the
  integer-`L` grid **cannot represent** — yet §4.2 routes exactly these
  (impure/interrupted) reads to the slow path and then stores `Qᵣ(L)` only over
  integer units;
- the TRTools harmonizer computes copy number as `len(allele)/len(motif)` (float
  division, `tr_harmonizer.py:740,757`) and **does not read REPCN at all**. So an
  imperfect REF tract whose length is not divisible by the motif yields a
  **fractional** TRtools copy number.

The spec never reconciles "we keep and genotype imperfect loci" with "an allele
is an integer repeat count." Either the allele model needs a
non-integer/interruption-aware representation, or imperfect loci must be
genotyped only on their integer-tiling component with the interruption pinned as
fixed context — and the spec must say which, and how `ref` is defined for them.

### H2 — The GangSTR auto-detection mechanism described does not work

**Resolution:** _resolved 2026-06-10._ Verified the full `InferVCFType` logic
(`tr_harmonizer.py:180-244`). Two corrections landed in §5.9:
- **Detection bullet rewritten.** Real rule is `command=` **and** `gangstr` both
  present in the lowercased raw header (two tokens, any lines). The honest recipe:
  our real `##command=ssr-genotype …` line + a truthful `##source=…GangSTR-compatible`
  line — no spoofed GangSTR command. Anti-collision avoidances listed
  (`hipstr`/`longtr`/`source=advntr`/`source=popstr`/symbolic `##ALT=<ID=STR\d+>`;
  also no INFO `VID`/`VARID`). **Corrected the `--vcftype gangstr` trap:** it is *not*
  a bypass — inference still runs and the flag only succeeds if the header is already
  gangstr-detectable, so the two tokens are mandatory either way.
- **Required-fields bullet split** into harmonizer-mandatory (INFO `RU`; FORMAT
  `GT`; `##FORMAT=<ID=Q,…>` header) vs GangSTR-conventional-but-unread
  (`END`/`PERIOD`/`REF`/`DP`/`REPCN`), so we stop implying TRtools enforces them.
- **Added a live-test requirement:** assert `TRRecordHarmonizer(our.vcf)` infers
  gangstr uniquely and round-trips, guarding against header regressions and TRtools
  detection changes.

Note: the REPCN-derived-copy-number half of the original H2 was already fixed under
H1 (§5.9 REF/ALT bullet). The gangstr path being minimal + fractional-tolerant
(copy number = len/motif) makes it a *good* target for the hybrid alleles, so the
gangstr-compat goal was kept, not reopened.

§5.9 says detection is "triggered honestly via a `##source` line stating
'GangSTR-compatible'." The harmonizer actually requires **both** `command=`
**and** `gangstr` in the lowercased raw header (`tr_harmonizer.py:210-213`):

```python
header = vcffile.raw_header.lower()
if 'command=' in header and 'gangstr' in header:
    possible_vcf_types.add(VcfTypes.gangstr)
```

A `##source=GangSTR-compatible` line will **not** auto-detect — there is no
`command=`. The only honest path is `--vcftype gangstr` (which works);
auto-detection requires emitting a `##command=…gangstr…` line, i.e. spoofing
GangSTR's invocation, which contradicts the "honestly" claim. Drop the
auto-detection-via-`##source` claim and standardise on `--vcftype gangstr`, or
accept the spoof and stop calling it honest.

Related: §5.9 leans on `REPCN` as the carrier of per-allele copy numbers, but the
harmonizer ignores REPCN and derives copy number from REF/ALT sequence length
(`tr_harmonizer.py:740,757`). The real hard requirement is therefore **REF/ALT
lengths must be exact integer multiples of the motif** — which loops back to H1
for imperfect loci. REPCN is only useful to *other* consumers, not TRtools.

### H3 — Confident-homozygote init is circular for the non-selfing case

**Resolution:** _resolved 2026-06-10._ Reframed from "settled clean defence" to
"honest seed + safeguards," after establishing two mitigations: (a) the
contamination bias is **precision-safe** (inflates `u/d` → suppresses hets), and
(b) it's **only a seed** — the M-step re-estimates the kernel from all loci, so it
shifts the starting basin, not the converged answer. Edits:
- **§5.4** — dropped "by construction every non-modal read is a stutter product";
  honest skewed-adjacent-het contamination caveat; named gates `SEED_MIN_DEPTH` +
  `SEED_MIN_DOMINANCE_FRAC` (~0.9, the primary control); a **per-cell measured
  cut-off diagnostic** (discriminating statistic = one-sided adjacent-neighbour
  excess, *not* raw dominance) that recommends the threshold and warns under strong
  bimodality, exposed as `--seed-dominance auto` — **shipped off by default, its
  reliability flagged as an open question to validate on real data**, promotable to
  default only if validation is convincing; a **seed fallback hierarchy**
  (cell-local → period-level → fixed default) so the kernel is always seedable.
- **§5.6** — single-sample seeds regardless of mating system (homozygosity is
  per-locus); it is the weakest identifiability regime (no population recurrence),
  precision-first no-call as safety valve; stripped the selfing-only flavour.
- **§5.8 + §7** — documented the **separate** call-time He/Ho downward bias (a
  property of the converged kernel + thresholds, not the seed), worst at small N,
  recovering with cohort size; added an Ho/He-vs-N validation check and a
  seed-cut-off-detection check.
- **§6/§8** — added the seed constants and the `--seed-dominance` flag.

User's framing accepted: the stringent dominance fraction is the key lever, and
since this only sets priors the residual risk is low; the auto-cut-off is a nice-to-
have pending empirical validation.

§5.4 seeds the stutter kernel from "confident homozygotes" (deep coverage +
single dominant length) and asserts "every non-modal read is, by construction, a
stutter product of that single allele." But those selection criteria are
**exactly** satisfiable by an **adjacent-unit heterozygote** (alleles differing
by one unit) with allele-balance skew from mapping/coverage — the very
hom-vs-adjacent-het ambiguity §5.4 admits is unresolvable within a locus. Such
loci contaminate the "labelled stutter observations" with real minor alleles,
biasing initial `u/d` upward.

The defense only holds when the cohort is overwhelmingly homozygous. The spec
quietly assumes this ("Selfing species supply abundant homozygotes," §5.6) — but
§1 targets popgen/diversity/GWAS broadly, which includes **outcrossing** cohorts,
and §5.6 makes **single-sample a first-class case**. Single-sample outcrosser is
the worst case for this seed and gets no treatment beyond "multi-restart
deferred." At minimum: state that init quality degrades with cohort
heterozygosity, and specify the fallback when confident homozygotes are scarce
(the genome-wide motif-class kernel of §5.6 should be the *primary* seed there,
not the per-cohort homozygote pass).

### H4 — The covariate-pooled kernel is an unproven bet sold as settled

**Resolution:** _resolved 2026-06-10._ Split into the two problems and fixed both.
Verified all three callers' stutter parameterization (GangSTR = fixed global
`0.05/0.05/0.9` + optional per-locus table; HipSTR = per-locus EM; popSTR2 =
covariate regression) — they bracket our design.
- **Problem A (bet sold as settled):** §5.2 now frames pooling as an explicit
  tradeoff against per-locus, names the regime where it's weakest — **impure loci**,
  because a scalar `purity` can't encode interruption structure (user's insight) —
  and adds a **per-locus relief valve + impure-locus policy** (per-locus when depth
  allows, else **no-call**; never force an impure locus onto a distrusted cell
  rate), with a **fixed-default-kernel** fallback (GangSTR's numbers).
- **Problem B (undefined functional form):** §5.2 now gives a concrete *provisional*
  model — the link-function invariant (multinomial-logit over no-slip/up/down +
  logit for ρ) plus **two structural options** with the choice deferred to data:
  Option 1 (cells + empirical-Bayes shrinkage, = v1) vs Option 2 (continuous
  hierarchical GLM, = upgrade). User deferred the modelling choice to me and asked to
  write options down and defer.
- **Honesty:** §11 now lists the kernel form as one of two genuinely-open *modelling*
  cores (with `--seed-dominance auto`); "settled" downgraded to architecture +
  generative model. §7 gains a covariate-kernel-fit study (per-cell fit + residual
  locus variation across size/GC/purity) to settle Option 1 vs 2 and the impure
  cutoff. §6 gains `DEFAULT_STUTTER_(U,D,RHO)` and `MIN_PERLOCUS_STUTTER_DEPTH`.

Decision deliberately deferred (per user): Option 1 vs 2, the impure per-locus/
no-call cutoff, and whether per-locus relief is needed — all pending real-data
behaviour across repeat size, GC, and impurity.

§5.2/§9 frame the shared covariate kernel
`θ(period, allele_length, motif, purity)`, pooled across loci, as an improvement
over HipSTR's per-locus EM. HipSTR genuinely fits a **separate stutter model per
locus** (`seq_stutter_genotyper.cpp:1466-1505`, an independent `EMStutterGenotyper`
per block). That per-locus choice is presumably compensating for stutter
variation **not captured by those four covariates** (flank context, local
secondary structure, alignment idiosyncrasies). The spec asserts pooling "pins it
down" without addressing what HipSTR's per-locus modeling was buying.

Worse, the functional form is undefined: "regression + hierarchical shrinkage"
with no link functions for the constrained parameters (`ρ ∈ (0,1)`, `u,d ≥ 0`,
`u+d ≤ 1`), no covariate transforms, no shrinkage hierarchy. This is a large,
genuinely unsettled modeling component, yet §11 says "everything else in this
spec is settled." It is not — this is the statistical core.

---

## Medium severity

### M1 — Candidate allele set is never defined

**Resolution:** _resolved 2026-06-11._ §5.1 now fully specifies `A_ℓ` via a
**threshold-free peak + ±1-adjacent rule** (converged with the user over two turns):

> `A_ℓ` (on-ladder) = evidenced rungs that are a **local maximum** of the
> cohort-aggregate length profile (support ≥ both neighbours) **or adjacent (±1) to
> one**; ∪ supported off-ladder sequences (cross-sample, normalized key).

The insight (user's): a real allele one repeat from a bigger one is **topologically
identical to a stutter satellite**, so don't try to distinguish them at prune time —
the **±1-adjacent clause keeps it unconditionally** and hands the real-vs-stutter
call to the EM (shared kernel + population recurrence), which *can* make it; the
**local-max clause** keeps ≥2-away real alleles (they peak). Everything pruned is a
non-peak ≥2 from every peak = unambiguous stutter tail/valley. Pruning removes a
genotype *candidate* but **not its reads** (they remain stutter observations for the
kept peak via the convolution), so the kernel is unaffected. Aggregate profile →
population recurrence protects real alleles (per-sample-peak admission is the
deferred escape hatch for ultra-rare alleles). §5.4 "self-prune" now means the EM
resolving the kept ±1 shoulders (`π → 0` + `MIN_ALLELE_SUPPORT_FRAC`). Cost backstop:
`MAX_CANDIDATE_ALLELES` → **no-call + log** (Q1 choice; never silent top-K), the real
defence against the hypervariable × high-ploidy blow-up (§5.7). Honest residual cost:
a real allele that *never* peaks anywhere and sits ≥2 from every peak is lost — at
the sensitivity floor, documented. Dropped the earlier `MIN_CANDIDATE_SUPPORT`
threshold (no tuning needed). New const: `MAX_CANDIDATE_ALLELES` (§6).

Stage 2's convolution `P(read|a,θ)=Σ_L Qᵣ(L)·S_θ(L|a)` and the genotype space
(integer partitions of ploidy) both require a per-locus **candidate allele set**,
but no section says how it is constructed. This is not a detail: stutter means a
true allele may appear only through its ladder (must be included even if not
directly observed as a confident length), and the set size drives both
correctness (missing a true allele) and cost (genotype enumeration is
combinatorial in |candidates| × ploidy). Needs its own subsection.

### M2 — Soft-clip recovery fights the spanning/flank requirement

**Resolution:** _resolved 2026-06-11._ Separated two notions §4.1 conflated:
- **"Spanning" is now a property of recovered *content*, not the CIGAR** — a read is
  spanning if its sequence carries ≥ `MIN_FLANK_BP` clean flank on both sides of the
  tract, whether aligned or **recovered from a soft-clip**. So soft-clipped reads
  carrying long alleles (the most valuable length evidence, previously disqualified
  by "matched flanks both sides") are now correctly counted as spanning.
- **Soft-clip recovery is explicitly a slow-path (realignment) job, not a fast-path
  count** — the fast path now requires both flanks *cleanly aligned* (no soft-clip),
  so soft-clipped reads route to the slow path by construction; the realignment both
  recovers the length and discriminates a real long-allele clip from junk
  (adapter/low-qual → no clean flank).
- **Accept/reject rule stated:** a clip yielding a clean flank → spanning long allele;
  a clip yielding no clean flank → "allele ≥ read length" → counted, not used (the
  §1.4 scope boundary). The slow-path candidate window centres on the *recovered*
  count, not the reference.

This recovers real long-allele sensitivity *within* the spanning-only scope (pure
win, no precision cost). **Knock-on to M4:** soft-clip recovery is always slow-path,
so a long-allele-rich locus does more pair-HMM work — fold into the M4 cost
discussion.

§4.1 correctly notes the aligner soft-clips the bases that carry a long allele,
and says to re-examine the whole read. But §4.1/§4.2 also require a spanning read
to have `≥ MIN_FLANK_BP` clean **matched** bases on both sides. For an allele
longer than ref, the downstream flank is often *in* the soft-clip — the reads
most diagnostic of length variation are the ones most likely to fail the "clean
matched flank both sides" test. And "re-examine the whole read" to recover the
count *is* local realignment, so the fast path's "exact motif count,
O(read length)" claim quietly depends on the slow-path machinery for exactly
these reads. The fast/slow boundary and the soft-clip handling must be worked out
jointly.

### M3 — Synchronized positional scan: brittleness + the "compresses to ~nothing" claim

**Resolution:** _resolved 2026-06-10 — dissolved by the format switch._ Stage-1
evidence moved off Parquet onto the project's `.psp` columnar-block container
(now `.ssr.psp`, sibling of the SNP `.snp.psp`; §2/§4.3). That format is **sparse
and coordinate-keyed**: no-coverage loci are simply absent (no dense empty rows, so
the "compresses to ~nothing" claim is moot), and cross-sample alignment comes from
the **genomic-window block grid** (every sample's blocks cover the same intervals)
rather than from a fragile row-index identity. The cohort scan merges by
`(chrom, start)` with no join and no empty rows. Decided in discussion: reuse `.psp`
for the same memory/perf control it already buys the SNP path, abstracting the
read/write/index/compression code into a shared container both callers ride on
(structural placement folded into decision E, §11). Cost accepted: no third-party
(Parquet) interop on an internal intermediate, one extra column registry to
maintain. Catalog↔evidence binding still via header md5s; **still worth** a defined
failure semantic when md5s mismatch (minor, fold into Stage-1 plan).

The one-row-per-catalog-locus-including-no-coverage design (§4.3) is elegant for
the N-file scan but: (a) it hard-couples every evidence file to one exact catalog
version/order — no subsetting, re-sorting, or incremental catalog edits without
rebuilding all evidence; the spec mentions `catalog_md5` binding but not the
failure semantics. (b) "no-coverage loci … compress to ~nothing" is asserted, not
shown. A genome-wide SSR catalog can be millions of loci; at hundreds of samples
that is millions × N rows of mostly-empty Parquet. Given this project's
memory/perf focus, the empty-row overhead at real cohort scale deserves a
measurement or a sparse alternative, not an assertion.

### M4 — Stage 1 cost dismissed with an irrelevant comparison

**Resolution:** _resolved 2026-06-11._ Replaced the "negligible next to the
alignment" line (§4.2) with an honest cost model + a measurement commitment: the
relevant baseline is Stage-1's own wall time (and the per-sample SNP pileup it
parallels), not the upstream alignment; cost = BAM/CRAM decode + per-read work
(fast-path counting vs slow-path pair-HMM); the **slow-path fraction is the driver**
and rises on repeat-rich/impure genomes and with soft-clip recovery (M2 knock-on).
Now framed as "benchmark on a repeat-rich genome, reporting Stage-1 wall + fast/slow
split + slow-path share," not asserted. Kept the (reasonable) intuition that it's
cheap — only the unmeasured assertion and wrong baseline were the problem.

§4.2 closes: "Performance is negligible … next to the whole-genome alignment
already paid for the BAM." Alignment was paid *once, before this tool runs* — it
is not re-paid, so it is not the relevant baseline. The relevant cost is Stage
1's own wall time: re-reading the BAM and running a banded pair-HMM over `2W+1`
haplotypes for every ambiguous read across (millions of loci × the reads at
each). For a project that measures everything, this hand-wave is out of character
and likely wrong on repeat-rich genomes. Flag it as "to be measured," not
"negligible."

### M5 — Polyploid prior formula is incomplete

**Resolution:** _resolved 2026-06-11._ §5.3 now gives the **general** form over a
genotype's allele counts `(k_a)`: fully-homozygous = `F·π_i + (1−F)·π_i^ploidy`;
≥2-distinct-allele (incl. partial polyploid hets `AABB`/`AAAB`) = the multinomial
HWE term with no `F` contribution. Shown to reduce to the diploid `P(ii)`/`P(ij)`
special case, confirmed = what the engine computes (enumerates ploidy-tuples,
`π^ploidy` term), and named as the **single-parameter polysomic** approximation
(`F` = full autozygosity only; partial IBD / double reduction deferred §5.7/§12).

§5.3 writes only the full-homozygote and (diploid) heterozygote terms:
`P(ii)=F·π_i+(1−F)·π_i²`, `P(ij)=(1−F)·2π_iπ_j`. The verified
`posterior_engine.rs` does generalize (it uses `π^ploidy` and enumerates all
non-decreasing ploidy-tuples — `genotype_order(ploidy, n_alleles)` in
`per_group_merger.rs:522`, homozygote term `(1−F)·π_i^ploidy + F·π_i` at
`posterior_engine.rs:3309-3331`), but the spec's *written* formula does not cover
**partial** polyploid genotypes (tetraploid `AAAB`, `AABB`), which need the full
polysomic IBD model, not a hom/het dichotomy. Since §5.7 claims polyploid
support, the spec's stated prior should either show the general form or
explicitly point at the engine's generalization and the single-parameter
polysomic approximation's limits.

---

## Lower severity / polish

### L1 — "Everything settled" vs the parameter table

**Resolution:** _resolved 2026-06-11 (via H4)._ §11 now says "'Settled' elsewhere
means the architecture and the generative model; several model defaults remain to be
swept (§6/§7)," and explicitly lists the two open modelling cores. No further edit.

§11 says everything but repo placement is settled, yet `AMBIGUITY_THRESHOLD`,
`alpha (α)`, base-measure decay `p`, the §5.2 covariate regression form, and the
M1 candidate-set rule all have `—`/no source. Downgrade §11's claim.

### L2 — `AMBIGUITY_THRESHOLD` is a bundle, not a threshold

**Resolution:** _resolved 2026-06-11._ Renamed to `FAST_PATH_GATE` throughout
(§4.2/§4.3/§6); §6 row now reads "predicate … a gate, not a tunable scalar / n/a".

The table calls it a threshold; §4.2 says it "bundles" three boolean gates
(flanks, zero interior indels, base-qual). It is a predicate, not a tunable
number — rename (e.g. `FAST_PATH_GATE`) or it will be mistaken for a sweepable
scalar.

### L3 — Out-of-frame stutter exists in the model you are borrowing

**Resolution:** _resolved 2026-06-11._ §5.2 already flagged it as a deliberate
simplification (added during H1); §12 now carries the matching deferred-item entry
(add an out-of-frame term if validation shows misattribution).

§5.2's "only in-frame; any non-unit-multiple change is Stage-1 sequencing error"
is a deliberate simplification — HipSTR itself carries a separate **out-of-frame**
stutter component in bp (`stutter_model.cpp:29-49`). Genuine out-of-frame
slippage will be misattributed to sequencing error under your split. That may be
acceptable, but present it as a known approximation with a cost, not as a clean
fact.

### L4 — M-step estimator wording

**Resolution:** _resolved 2026-06-11._ §5.4 stutter M-step now states the
normalisation explicitly: `u`/`d` = posterior-weighted gain/loss counts over **all**
events (no-slip included), not longer-vs-shorter among slips; `ρ` = geometric MLE.

§5.4/§5.6 describe HipSTR's updates as "`u`=fraction longer, `d`=fraction
shorter." HipSTR actually normalizes `u,d` over **all** posterior-weighted events
including no-stutter and out-of-frame (`em_stutter_genotyper.cpp:117-121`), not
just longer-vs-shorter among stutter events. Tighten the wording so the
implementer copies the right denominator.

### L5 — `reference_md5` is underspecified

**Resolution:** _resolved 2026-06-11._ Glossary `md5` entry now defines
`reference_md5` = md5 of the concatenated upper-cased per-contig sequence (CRAM `M5`
convention), not file bytes (invariant to wrapping/compression); `catalog_md5` over
the data rows.

md5 of the FASTA file, the concatenated sequence, or per-contig M5 tags? Tools
disagree; pin it so catalog/evidence/VCF binding is reproducible.

### L6 — Read-pair merging (§1.4)

**Resolution:** _resolved 2026-06-11._ §1.4 now notes merging is not free: combined
base-quality model at the overlap, conflict handling, and chimeric-merge-inside-the-
repeat risk → merged fragments get the same slow-path scrutiny as soft-clipped reads,
not blind concatenation.

Merging changes the base-quality model at the overlap (combined quals) and can
introduce chimeras precisely inside the repeat — non-trivial for the pair-HMM
emission model. "Cheap, in-paradigm" undersells it; note the quality-model
implication.

---

## Verified-correct claims (no action)

These spec claims were checked against source and **hold** — recorded so we do
not re-litigate them:

- The SNP caller's prior **is** the exact F-mixture the spec reuses, **is**
  genuinely multiallelic, and **does** generalize to ploidy ≤ 8 with a reusable
  `(π, F)` interface (`posterior_engine.rs:3309-3331`, `per_group_merger.rs:522`).
  Caveat folded into M5 (written formula is diploid-specific).
- HipSTR's stutter kernel **is** the 3-parameter geometric `(u, d, ρ)` form the
  spec reproduces, in whole repeat units for the in-frame component
  (`stutter_model.h:17-28`, `stutter_model.cpp:29-52`). Caveats folded into L3/L4.
- TRTools reads genotype quality from the scalar `Q` FORMAT field for GangSTR
  (`tr_harmonizer.py:333,1610,1733`); the spec's "TRtools reads only the scalar
  `Q`" is correct.

---

## Suggested order of work

H1 and H2 can invalidate concrete outputs — settle them before any implementation
plan. H3/H4 are methodological risks to re-frame as bets-with-fallbacks rather
than settled decisions. M1 (candidate set) is a missing core subsection. The rest
is tightening.
