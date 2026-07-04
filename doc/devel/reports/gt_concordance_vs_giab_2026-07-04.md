# Genotype-concordance gap vs freebayes (GIAB per_sample) — root cause

**Date:** 2026-07-04
**Scope:** Why our per-sample genotype (GT) concordance at true-positive
variants trails freebayes on the GIAB benchmark, in two specific ways:
low-coverage SNPs, and indels at all depths.
**Data:** `benchmarks/giab/results/per_sample/<cov>/{high-recall,freebayes}`,
GIAB v4.2.1 truth, HG002/HG003/HG004. Single-sample var-calling (n = 1).
**Analysis scripts:** `tmp/gt_confusion.py`, `tmp/indel_reffrac.py`.

## TL;DR

Both failure modes are the **same observable transition** — true hom-alt
(`1/1`) called het (`0/1`) — but they have **two distinct root causes**:

1. **SNP, low coverage (≤ 10×): prior-driven het over-call.** The
   single-sample EM's REF Dirichlet pseudocount (`α_ref = 10`) pins the
   estimated ALT frequency near `p̂_alt ≈ 0.08` regardless of the read
   data, and the Wright/HWE genotype prior built from that `p̂` then makes
   `P(het) ≫ P(hom-alt)`. With few reads the likelihood can't overcome it,
   so a genuinely-sampled hom-alt (AD `0,2`) is called het with **low GQ**.
   Resolves by ~15× once there are enough ALT reads to swamp the 10-count
   prior. **Fix lives in the genotyper/prior.**

2. **Indels, all depths: likelihood-driven het call from inflated REF
   support.** At homopolymer/short-tandem-repeat sites (86 % of the
   mismatches), ~15–20 % of reads at a true hom-alt indel are (mis)assigned
   to the REF allele — stutter/shifted alignments that our short local-indel
   representation buckets as REF. That genuine intermediate VAF makes the
   genotyper **correctly** call het (high GQ, confidently wrong). Freebayes
   represents the same site as a long haplotype window and assigns those
   reads to ALT, so it sees ~0 % REF. Depth-independent because the
   mis-assigned fraction is systematic, not sampling. **Fix lives upstream in
   indel allele assignment / an indel error (stutter) model — not the
   genotyper.**

## Characterization — confusion matrix over mismatched TP sites

Truth-GT → called-GT among TP sites (cohort-summed HG002/3/4), from
`tmp/gt_confusion.py`. Dominant transition and GQ of the mismatched calls:

### SNPs
| cov | ours mismatch | dominant | GQ(mm med) | freebayes mismatch | fb dominant |
|----:|--------------:|----------|-----------:|-------------------:|-------------|
| 5×  | 239/1455 (16.4%) | `1/1→0/1`: 214 | **5** | 80/1206 (6.6%) | `0/1→1/1`: 75 |
| 10× | 42/1881 (2.2%) | `1/1→0/1`: 33 | 5 | 28 (1.7%) | `0/1→1/1`: 20 |
| 15× | 8 (0.4%) | `1/1→0/1`: 7 | 10 | 14 (0.8%) | `1/1→0/1`: 10 |
| 30× | 5 (0.2%) | `1/1→0/1`: 5 | 5 | 11 (0.6%) | `1/1→0/1`: 11 |
| 300×| 4 (0.2%) | `1/1→0/1`: 4 | 99 | 10 (0.5%) | `1/1→0/1`: 10 |

The SNP gap is **entirely** the 5–10× `1/1→0/1` over-call at low GQ. It
collapses to ≈ freebayes by 15×. (The minority `0/1→1/1` at 5× is sampling —
at 5× a true het is often sequenced all-ALT; freebayes does this *more*, 75
vs our 25.) The handful of high-GQ `1/1→0/1` SNPs at 300× (AD like `56,103`)
are a separate, tiny collapsed-repeat/paralog residue, not the gap.

### Indels
| cov | ours mismatch | dominant | GQ(mm med) | freebayes mismatch |
|----:|--------------:|----------|-----------:|-------------------:|
| 5×  | 43/175 (24.6%) | `1/1→0/1`: 42 | 34 | 21/169 (12.4%) |
| 10× | 54/255 (21.2%) | `1/1→0/1`: 54 | 52 | 13 (5.7%) |
| 30× | 69/307 (22.5%) | `1/1→0/1`: 69 | 99 | 9 (3.0%) |
| 50× | 72/312 (23.1%) | `1/1→0/1`: 72 | 99 | 6 (2.0%) |
| 300×| 82/310 (26.5%) | `1/1→0/1`: 82 | **99** | 5 (1.6%) |

Flat ~21–27 %, **all** `1/1→0/1`, GQ rising to 99 with depth (confidently
wrong). Freebayes falls monotonically to 1.6 %.

## Evidence — SNP low-coverage mechanism (the prior)

`src/var_calling/posterior_engine.rs`: `ref_pseudocount = 10.0`,
`snp_alt_pseudocount = 0.01`, `indel_alt_pseudocount = 0.00125`,
`fixation_index_default = 0.0` (strict HWE). The M-step
(`m_step_p_hat`, line ~2828) is a Dirichlet posterior mean:
`p̂_alt = (e_alt + α_alt) / (Σe + α_ref + α_alt)`. The E-step
(`e_step`, ~2531) builds the per-genotype Wright/HWE prior from `p̂`.

Trace for a single hom-alt sample, AD `0,2` (2 ALT, 0 REF):
- Iter 1 (flat `p̂ = 0.5`): likelihood favors hom-alt ~4:1 → posterior
  `P_AA ≈ 0.67, P_RA ≈ 0.33` → `e = [0.33, 1.67]`.
- M-step: `p̂_alt = (1.67 + 0.01)/(2 + 10.01) ≈ 0.14`.
- Iter 2 (`p̂ = [0.86, 0.14]`): HWE prior `[0.74, 0.24, 0.02]`;
  posterior ∝ `[~0, 0.24·0.25, 0.02·1] = [~0, 0.06, 0.02]` →
  **`P_RA ≈ 0.75` → het**, GQ ≈ 6. Matches observed (`GQ=7, 0/1`).

**Fingerprint:** at het-called sites the AF INFO is a constant
`AF = 0.08342881` across sites with completely different AD — exactly
`(1 + 0.00125)/(12.00125)` for indels (`(1+0.01)/12.01 ≈ 0.084` for SNPs),
i.e. the value produced when `expected_counts = [1,1]` (het) and the α_ref=10
Dirichlet fixes the frequency. The population-level prior fully controls the
n = 1 allele frequency.

Because n = 1, there is no cohort to overpower the prior — the 10-count
REF pseudocount behaves like 10 phantom REF chromosomes against ~2 real
reads. This is the "strong prior + few reads pulls hets" hypothesis, and it
lands on hom-alt→het (not het→hom-ref).

## Evidence — indel mechanism (upstream allele assignment)

From `tmp/indel_reffrac.py` (ref-fraction `AD[0]/(AD[0]+AD[1])` at OURS indel
TP sites):

| cov | matched indels | `1/1→0/1` mismatches |
|----:|----------------|----------------------|
| 5×  | rf med 0.20, rep-ctx 96% | rf med 0.17, rep-ctx 83% |
| 30× | rf med 0.52, rep-ctx 93% | rf med **0.16**, rep-ctx 86% |
| 300×| rf med 0.56, rep-ctx 93% | rf med **0.15**, rep-ctx 88% |

(Matched indels sit at rf ≈ 0.5 because most matched TP indels are true
hets.) The mismatches are true **hom-alt** sites carrying a persistent
~15–20 % REF fraction, **86 %** in homopolymer/STR context — depth-invariant.

Side-by-side, same sites, 300× HG002 (ours vs freebayes AD):
```
chr1:1388565  TC-repeat del   ours GTCTC>GTC  AD 82,208 → 0/1   fb (20bp window) AD 1,207 → 1/1
chr1:4927622  polyA ins       ours TA>TAA     AD 38,208 → 0/1   fb AD 4,202 → 1/1
chr1:21684679 repeat del      ours CTTACAT>C  AD 28,238 → 0/1   fb AD 0,361 → 1/1
chr1:68898168 polyAT ins      ours C>CAT      AD 24,225 → 0/1   fb AD (window) → 1/1
```
Our raw calls at these sites are **biallelic** (single ALT) — we don't emit
the stutter alleles; the shifted/stutter reads collapse into REF. Freebayes'
window/haplotype observation captures them on ALT. Given AD `82,208`, a het
call is statistically *correct* (82 REF reads are impossible under hom-alt at
ε≈0.01) — so the defect is the AD, not the genotyper. GQ=99 follows.

## Confounders ruled out
- **Allele-balance filter:** it removes calls (→ FN), it does not flip GT;
  GT concordance is measured on survivors. The het mismatches pass the
  filter, and it is depth-gated inert < ~30× so it can't drive the 5× SNP
  gap.
- **min-qual gating:** applied before this comparison; concordance is on
  survivors only.
- **Multiallelic split (`norm -m -any`):** handled — comparison recodes GT
  phase-insensitively per split record on both sides.
- **Not cohort AF/HWE across many samples:** this *is* single-sample; the
  prior problem is specifically the n = 1 regime.

## Proposed fixes (ranked, not yet implemented)

### SNP low-coverage
1. **Make the REF Dirichlet concentration cohort-size-aware** so it can't
   dominate n = 1 (e.g. scale `α_ref` toward ~1 as n → 1, or treat the
   pseudocount as a fixed *concentration* spread over chromosomes rather than
   an absolute count). Lowest-risk, directly targets the mechanism.
2. **Decouple the per-sample genotype prior from the EM AF at low evidence**
   — use a flatter genotype prior (or a floor on `p̂_alt` for the
   genotype-prior evaluation) when the site has few reads.
3. Validate any change against **FP/precision at 5–30×** — `α_ref = 10` is
   REF-favoring and currently helps low-VAF FP suppression; the risk is
   trading the het-recall win for FP. Quick experiment: rerun 5×/10×/15× with
   a reduced/scaled `α_ref` and re-read both the GT panel and the FP panel.

### Indels (the bigger, depth-independent gap)
1. **Stutter-aware indel genotype likelihood** (GATK/DRAGEN/HipSTR style):
   model the elevated REF-looking read rate as a function of repeat length so
   ~15–20 % REF at an STR is not read as het. Contained to the indel
   likelihood; no pileup rearchitecture. The SSR module already has a
   HipSTR-style stutter model (`src/ssr/cohort/read_model/`) to borrow from.
2. **Haplotype-window allele observation for indels** (freebayes-style):
   assign shifted/stutter reads to the ALT haplotype at the source. Correct
   but a larger change to the pileup/merge path.
3. **Left-alignment is not the culprit** here (the F3 GATK `leftAlignIndels`
   port already runs) — the issue is per-read allele *bucketing* in repeats,
   not normalization of the called allele.

## Reproduce
```
uv run tmp/gt_confusion.py      # confusion matrix + GQ + examples, all covs
uv run tmp/indel_reffrac.py     # indel ref-fraction, match vs mismatch, rep-ctx
```
Panel of record: "Genotype concordance vs GIAB (TP variants)" in
`benchmarks/giab/src/freebayes_comparison_dashboard.py`.
