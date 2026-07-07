# Heuristic / ratio-based SNP calling — inspiration for the Stage-1 rough caller

**Date:** 2026-07-07
**Scope:** the *rough* per-site SNP genotype call in Stage-1 pileup
([`src/sample_summary/het.rs`](../../../../src/sample_summary/het.rs)),
whose only downstream product is the observed-het **rate** feeding the
inbreeding coefficient `F` and (via the `.psp`) candidate-site flagging.
SSR deferred to a follow-up.

## 0. What the rough caller does today, and what it is *for*

`HetAccumulator::observe_site` scores each candidate site with three
binomial models and one margin `M`:

```
hom-ref :  alt ~ Binomial(n, ε)        LL = k·ln ε      + r·ln(1−ε)
het     :  alt ~ Binomial(n, ½)        LL = n·ln ½
hom-alt :  alt ~ Binomial(n, 1−ε)      LL = k·ln(1−ε)   + r·ln ε
```

- **Variant gate:** admit iff `max(LL_het, LL_hom_alt) − LL_hom_ref > M`.
- **Split:** `logLR = LL_het − LL_hom_alt`; het if `> +M`, hom-alt if
  `< −M`, else ambiguous.

This is already "simple arithmetic" (a few multiplies, the binomial
coefficient cancels). Its **inputs are the two raw fragment counts**
`(ref_obs, alt_obs) = num_obs` and **one global flat `ε`**. Everything
else about the reads is discarded.

**Why accuracy here matters narrowly.** The consumer is
[`inbreeding.rs`](../../../../src/paralog/inbreeding.rs):
`obs_het = n_het_sites / callable_positions`, then
`F = clip(1 − obs_het / Hexp, 0, 0.99)`. So the caller is graded almost
entirely on **`n_het_sites`**:

- a spurious het site (strand artifact, mismapping, indel-adjacent
  noise, collapsed paralog) **inflates `obs_het` → deflates `F`**;
- `n_hom_alt` and `n_ambiguous` barely matter (they are *not* in the
  numerator — this is the deliberate "het rate, not het/variant ratio"
  choice, whose inversion was already measured: Spearman −0.28 vs +0.86).

So the single highest-value thing a heuristic caller can teach us is
**cheap rejection of false heterozygous sites**. Sensitivity (recall of
true hets) matters second; the depth-aware ambiguous bucket already
handles the low-depth tail gracefully.

## 1. The load-bearing finding: the signals are already in the record

The `.psp` per-allele [`AlleleSupportStats`](../../../../src/pileup_record.rs#L44-L66)
**already carries every field the heuristic-FP-filter literature uses** —
computed once in Stage-1, freebayes-style, and today read *only* by the
Bayesian cohort caller:

| Field | What it gives us for free | Heuristic-caller analogue |
|---|---|---|
| `num_obs` | ref/alt counts (used today) | VarScan min-reads2, freebayes minAltCount |
| `q_sum` = Σ max(ln P_err(BQ,BAQ), ln P_err(MQ)) | **quality-weighted** error mass per allele | bcftools QS, LoFreq per-base quality |
| `fwd` | forward-strand count → 2×2 strand table | bcftools FS/SP, freebayes SAP/SRP, VarScan strandedness |
| `placed_left`, `placed_start` | read-placement distribution | freebayes RPP/EPP, bcftools RPBZ, VarScan readpos/dist3 |
| `mapq_sum`, `mapq_sum_sq` | per-allele MAPQ mean + variance → Welch's t | bcftools MQBZ, VarScan mapqual-diff (the paralog/mismap signal) |

**Reframing:** "improve the rough caller" is *not* a plumbing job. Two of
the three highest-value upgrades below cost only arithmetic on fields the
record already holds — the same fields, notably, that the existing
paralog filter and MAPQ-diff work already trust downstream.

## 2. How the surveyed heuristic callers actually decide (the menu)

### VarScan 2 — pure ratio + Fisher gate
Genotype straight from VAF bands after fixed evidence gates. Defaults
(`mpileup2snp`): min-coverage 8, **min-reads2 2**, min-avg-qual 15,
**min-var-freq 0.01**, **min-freq-for-hom 0.75**, p-value 0.99. Rule:
VAF < min-var-freq → ref; min-var-freq ≤ VAF < 0.75 → **het**; ≥ 0.75 →
hom-alt. The call is admitted by a **Fisher exact test** of
(ref, alt) counts vs. the sequencing-error baseline. No likelihood, no
prior — a threshold ladder plus one 2×2 test.

### VarScan `fpfilter` (via bam-readcount) — the false-positive catalogue
This is the richest cheap-heuristic set found, and it is *entirely* a
false-het rejector. Defaults:

- `min-var-count 3` (min-var-count-lc 1 at low coverage)
- `min-var-basequal 30`
- `min-strandedness 0.01` — alt must appear on **both strands**
- `min-var-readpos 0.15`, `min-var-dist3 0.15` — alt reads not clustered
  at read ends / 3′ ends
- `max-mapqual-diff 10` — |mean MAPQ(ref) − mean MAPQ(alt)| bounded
  (**mismapping/paralog signal**)
- `max-var-mmqs 100`, `max-ref-mmqs 50`, mmqs-diff — mismatch-quality-sum
  bound (**alignment-noise / paralog signal**)
- `max-rl-diff 0.05` — read-length skew between ref/alt reads

Every one of these maps to a field we already have (`fwd`, `placed_*`,
`mapq_sum{,_sq}`, `q_sum`) **except** mmqs and read-length, which we do
not store.

### LoFreq — quality-*aware* counts, no flat ε
Core idea worth stealing conceptually: instead of a single global `ε`,
each observation contributes its **own** error probability; the "alt is
all error" null becomes a Poisson-binomial tail, not a fixed-`ε`
binomial. Our `q_sum` is exactly the pre-summed ingredient for a
cheap version of this (see §3.1). LoFreq then strand-bias-filters
(`lofreq filter`) as a separate pass.

### SNVer — one-sided binomial strand test
Rejects a candidate if the alt allele is strand-skewed via a one-sided
binomial on (alt-fwd, alt-rev). A one-liner FP guard.

### GATK hard filters / allele-balance literature
The germline hard-filter thresholds are the community's distilled cheap
signals: **FS** (Fisher strand, phred), **SOR** (strand odds ratio,
end-of-target-safe alternative to FS), **MQRankSum** (MWU of MAPQ
ref-vs-alt, one-sided < −12.5), and **allele balance AB ≈ 0.5** modeled
as Binomial(D, ½) — recurrent AB deviation is a documented systematic
false-het source.

### bcftools / htslib mpileup model (vendored, exact arithmetic)
Grounded in `bcftools/bam2bcf.c` + `htslib/errmod.c`:

- **Neighbor-quality reduction** (`bam2bcf.c:430`): a base's Q is capped
  by its neighbours' Q + δ (δ=30) — cheap correlated-error damp.
- **Quality composition** (`:437-463`): `q = min(BQ, seqQ, MQ)`, MQ0
  reads floored to BQ 4, clamp `[4,63]`. `QS[allele] += q` is the
  load-bearing evidence number.
- **`fk[n]` geometric discount** (`errmod.c:77`):
  `fk[n] = (1−depcorr)ⁿ·(1−η) + η`, η=0.03. The *n*-th identical
  (base,strand) read is geometrically down-weighted — a stack of
  identical same-strand reads **cannot** accumulate unbounded evidence.
  This is the principled version of "don't trust 40 identical reads as
  40× the evidence," directly relevant to collapsed-paralog piles.
- **Strand bias**: Fisher exact on the 2×2 (`bam2bcf.c:1156`); per-sample
  `SP` computed **only when every margin ≥ 2** (`:1361`) — the guard that
  stops the test firing on tiny counts.
- **MWU-Z biases** (`calc_mwu_biasZ`, `:817`): MQBZ (one-sided, alt
  lower MAPQ), RPBZ (read position), BQBZ, MQSBZ (MAPQ×strand), NMBZ
  (mismatch count), SCBZ (soft-clip). Closed-form `Z=(U−m)/√var`.

### freebayes (pre-Bayesian gates only)
`minAltFraction 0.05` **and** `minAltCount 2` on a single sample
(`AlleleParser.cpp:3927`) — the primary false-het rejector — plus a
read-level mismatch-fraction reject (alignment-noise reads dropped
before they vote), and Hoeffding-bound phred bias annotations
(SAP/SRP strand, RPP/EPP placement, ABP allele-balance):
`hoeffdingln(succ, trials, ½) = ln½ − 2(trials·½ − succ)²/trials`.

## 3. Concrete candidate improvements, ranked by (value ÷ cost)

"Value" = expected reduction in false-het contamination of `obs_het`.
"Cost" assumes byte-identity re-baselining of the `.psp` het counts.

### 3.1 Replace the flat `ε` with the record's `q_sum` — **do first**
Today `ln_eps` is a CLI constant applied to *every* site regardless of
the actual base/mapping quality of the reads there. But `q_sum` already
holds `Σ ln P_err` per allele (BQ+BAQ+MQ combined, freebayes prodQout).
A site's effective error rate is `ε̂ = exp(q_sum_over_all_alleles /
num_obs_total)` — a **per-site, quality-aware** ε at zero new plumbing.
This directly imports LoFreq's central idea in its cheapest form and
makes the hom-ref/variant gate tighter where quality is high and more
forgiving where it is low. Likely the single biggest precision lever.
*Open question:* whether to feed `q_sum` per-allele into per-model LLs
(closer to bcftools QS) or just derive one site `ε̂` (minimal change).

### 3.2 A strand-bias veto on het admission — **cheap, high value**
Before counting a site as het, compute Fisher exact (or the cheaper
freebayes Hoeffding `SAP`) on the `(ref_fwd, ref_rev, alt_fwd, alt_rev)`
table from `fwd` — **guarded by "every margin ≥ 2"** (bcftools `:1361`)
so it never fires on thin counts. If the alt allele is strand-confined,
demote het → not-counted (or → ambiguous). Kills the classic
strand-artifact false het. `fwd` is already in the record.

### 3.3 A MAPQ-diff / MWU-Z veto — **reuses existing paralog signal**
`mapq_sum`, `mapq_sum_sq` give mean+variance per allele → Welch's t or a
one-sided MWU-Z of alt-vs-ref MAPQ (bcftools MQBZ; VarScan
`max-mapqual-diff 10`; GATK MQRankSum < −12.5). Alt reads with
systematically lower MAPQ = mismapping/paralog → demote the het. This is
the **same signal** the downstream MAPQ-diff filter and coverage-paralog
work already rely on (see memory: `mapq_diff_filter_fn`,
`sfs_prior_m4` correlated-mismapping tail) — applying it *at the rough
stage* would stop those artifacts from ever entering `obs_het`/`F`.

### 3.4 A min-alt-count + min-VAF floor on het — **trivial**
freebayes/VarScan both gate on `alt_count ≥ 2` AND `VAF ≥ ~0.05`. Our
depth-aware binomial gate already suppresses a lone alt read at high
depth, but an explicit `alt_obs ≥ 2` floor cheaply removes the singleton
class entirely and matches the field's consensus. One comparison.

### 3.5 Geometric down-weighting of stacked identical reads — **research**
`errmod.c`'s `fk[n]` discount is the principled fix for collapsed-paralog
piles inflating het/hom-alt confidence. We do **not** currently store the
per-(base,strand) multiplicity needed to apply it exactly, but the chain
IDs / `fwd` split give a coarse approximation. Higher cost; park until
3.1–3.4 are measured.

### 3.6 Not now: mmqs, read-length skew, soft-clip, NM bias
VarScan's mismatch-quality-sum and freebayes' NM/soft-clip biases are
strong FP signals but need fields the record does **not** carry
(per-read mismatch sum, read length, NM). Adding them is a Stage-1
pileup change with its own byte-identity cost — defer unless 3.1–3.4
leave a measurable false-het residue traceable to alignment noise.

## 4. Suggested next step

Prototype **3.1 (q_sum → per-site ε̂)** and **3.4 (alt-count floor)**
first — both are pure `het.rs` changes on fields already present, both
are individually measurable against the tomato/GIAB `F` and het-rate
baselines, and neither touches the `.psp` format. Then add **3.2 strand
veto** and **3.3 MAPQ-diff veto** as separate, independently-gated
increments. Each is a small, testable step; validate the `F` /
het-rate delta before stacking the next (per the incremental-steps
working rule).

## Sources

Web: VarScan manual & FAQ (`varscan.sourceforge.net`), VarScan
`fpfilter` defaults (`dkoboldt/varscan`, `ckandoth/variant-filter`),
LoFreq (Wilm et al. 2012, NAR 40:11189), SNVer manual, GATK
hard-filtering germline short variants, allele-balance-bias
(PMC6587442), Nielsen et al. 2011 genotype/SNP-calling review
(PMC3593722).

Vendored code (paths under repo root): `bcftools/bam2bcf.c`,
`htslib/errmod.c`, `htslib/realn.c`, `bcftools/mpileup.c`,
`freebayes/src/{AlleleParser,ResultData,Parameters,Utility}.cpp`.
Full extraction in the research agent's reference sheet (this session).
