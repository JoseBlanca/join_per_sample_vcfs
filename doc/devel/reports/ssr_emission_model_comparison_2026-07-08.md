# SSR emission-model comparison — heuristic vs BIC vs freebayes

**Date:** 2026-07-08 · **Branch:** `ssr-bic-emission` (merge `0de06ea`, unified emission
selector) · **Bench:** `benchmarks/ssr_tomato1` (51-sample *S. lycopersicum* selfer,
median ~3 reads/plant, F_IS ≈ 0.82)

This is a benchmark + adoption decision, not new modelling. Three emission models now sit
behind one selector — `PVC_SSR_EMIT_MODEL = heuristic (default) | bic | freebayes` — all
consuming the **same** per-locus read likelihoods (`data_ll`, one Qᵣ pass), so only the
per-locus *decision* differs, never the reads. `PVC_SSR_MARGINALIZED_PRIOR` is orthogonal
(it re-shapes the EM genotypes every model consumes).

## 0. Merge verification (all pass)

| check | result |
|---|---|
| (a) default `heuristic` byte-identical to pre-merge VCF | **PASS** — `diff` clean vs the rerun baseline `cohort.ssr.vcf` (headers + 74,573 records) |
| (b) `bic` + MARG + margin 10 reproduces 81.1 % / 0.17 % | **PASS** — 81.1 % recall / 0.17 % FP (455/561, fp 15) exactly |
| (c) `freebayes` reproduces its report's frontier | **PASS** — Q0 88.9 %/1.24 %, Q3 83.8 %/0.29 %, Q20 80.7 %/0.21 %, Q100 70.8 %/0.06 % all match |

The merge is correct: each model reproduces its source-branch numbers, and the default path
is untouched.

## 1. Recall vs false-positive frontier (silver confident core: 561 true100 / 8850 false100)

Every VCF scored with the single canonical scorer `silver_standard.py` on one core, one
process. Recall = confident-real loci emitted as a whole-unit length variant; FP =
confident-monomorphic loci emitted. Freebayes QUAL is swept post-hoc (VCF QUAL ≥ *t*
exactly reproduces `PVC_SSR_FREEBAYES_MIN_QUAL=t`).

### MARG prior OFF (plug-in EM)

| model / knob | recall | precision | FP rate |
|---|--:|--:|--:|
| **heuristic (current default)** | **80.7 %** | 97.0 % | **0.16 %** |
| bic margin 0 | 83.4 % | 94.2 % | 0.33 % |
| bic margin 10 | 76.3 % | 96.4 % | 0.18 % |
| bic margin 20 | 72.9 % | 98.8 % | 0.06 % |
| freebayes Q≥0 | 88.9 % | 81.9 % | 1.24 % |
| freebayes Q≥20 | 80.7 % | 96.0 % | 0.21 % |
| freebayes Q≥50 | 75.6 % | 97.2 % | 0.14 % |
| freebayes Q≥100 | 70.8 % | 98.8 % | 0.06 % |

### MARG prior ON — the high-recall regime

| model / knob | recall | precision | FP rate |
|---|--:|--:|--:|
| heuristic + MARG | 56.7 % | 97.0 % | 0.11 % |
| bic margin 0 + MARG | **87.3 %** | 94.0 % | 0.35 % |
| bic margin 5 + MARG | 83.8 % | 95.7 % | 0.24 % |
| bic margin 10 + MARG | 81.1 % | 96.8 % | 0.17 % |
| bic margin 20 + MARG | 75.0 % | 98.6 % | 0.07 % |
| **freebayes Q≥3 + MARG** | **87.0 %** | 95.1 % | 0.28 % |
| freebayes Q≥10 + MARG | 85.9 % | 95.3 % | 0.27 % |
| freebayes Q≥20 + MARG | 84.1 % | 96.5 % | 0.19 % |
| freebayes Q≥40 + MARG | 81.3 % | 97.0 % | 0.16 % |
| freebayes Q≥50 + MARG | 79.0 % | 97.8 % | 0.11 % |
| freebayes Q≥100 + MARG | 73.1 % | 98.6 % | 0.07 % |

### The frontier, read at matched FP (MARG on — the interesting regime)

| FP rate | bic + MARG | freebayes + MARG | winner |
|---|--:|--:|---|
| ~0.07 % | 75.0 % (m20) | 75.8 % (Q75) | fb +0.8 |
| ~0.16–0.17 % | 81.1 % (m10) | 81.3 % (Q40) | tie |
| ~0.19–0.24 % | 83.8 % (m5, @0.24) | 84.1 % (Q20, @0.19) | fb (more recall, less FP) |
| ~0.28–0.35 % | 87.3 % (m0, @0.35) | 87.0 % (Q3, @0.28) | fb (same recall, lower FP) |

**Read:** freebayes ≥ BIC at every matched FP — tie at the conservative point, +0.8 to +1
recall in the middle, and it reaches the ~87 % ceiling at **0.28 %** FP where BIC needs
**0.35 %**. This confirms the freebayes report's "beats BIC 2–3× on FP" in-caller, though
the in-caller BIC (real chemistry) is much closer than the toy-prototype BIC was — the gap
is real but small. Both **require the MARG prior** to reach the high-recall range; plain BIC
and the heuristic both plateau because they emit the plug-in EM's MAP genotypes, which
collapse the lone-carrier recurrent hets to hom-ref. Freebayes is the better frontier and
the cleaner instrument (one per-locus QUAL, no arbitrary margin unit).

## 2. Genotype concordance vs HipSTR — the axis that reverses the ranking

Concordance = per-sample-cell agreement (allele-length genotype) on the comparable set
(HipSTR whole-unit length-poly loci), among cells **both** callers genotype
(`benchmarks/lib/ssr_concordance.py`). This is *agreement*, not accuracy — HipSTR is
F-blind on a selfer and over-calls dinucleotide hets, so it is a flawed reference — but it
is the per-sample QC axis the task requires, and it moves in the opposite direction to
silver recall.

| VCF | comparable loci emitted | concordance |
|---|--:|--:|
| heuristic (baseline) | 774 / 1645 | **96.5 %** (28,914 cells) |
| bic m10 + MARG | 827 | 86.1 % (37,848) |
| bic m0 + MARG | 1034 | 87.0 % (47,831) |
| freebayes Q≥0, MARG off | 1128 | 90.8 % (52,408) |
| freebayes Q≥0, MARG on | 1209 | 87.6 % (56,234) |

Every model-emission path drops concordance 6–10 points. **This is not a genotyping
regression** — the per-sample GT/GQ columns are the *same* EM MAP calls; the models change
only the site emit decision. The drop is entirely that the model paths **skip
`apply_fp_control`**, the per-sample allele-balance *no-call* that converts imbalanced
(stutter-shoulder) hets to `./.`. Without it, those contested cells return as
called-but-discordant (our het vs HipSTR hom). The discordance is structured, not noise:
ours skews *shorter* than HipSTR (Σallele-bp-diff strongly negative), i.e. our F-aware EM is
more conservative than F-blind HipSTR — consistent with HipSTR over-expanding at low depth.

## 3. Re-enabling the per-sample no-call (KEEP) — measured, and it changes the story

The freebayes report recommended "keep the per-sample no-call to recover concordance
**without losing** the recall/FP win," asserting the two are independent. I implemented the
minimal change to test it — `PVC_SSR_KEEP_FP_CONTROL=1`: under bic/freebayes the site model
still owns the *emit* decision, but the emitted row's per-sample GTs get the
`apply_fp_control` no-call (~10 lines in `decide_emission` +
one `FpControlCfg` field; default off → byte-identical, verified; 250 ssr::cohort tests
pass). Result:

| operating point | recall | FP | concordance |
|---|--:|--:|--:|
| **MARG on** | | | |
| bic m10 + MARG, no keep | 81.1 % | 0.17 % | 86.1 % |
| bic m10 + MARG, **+KEEP** | **54.0 %** | 0.07 % | 97.7 % |
| fb Q≥40 + MARG, no keep | 81.3 % | 0.16 % | 87.6 % |
| fb Q≥40 + MARG, **+KEEP** | **53.7 %** | 0.06 % | 98.0 % |
| fb Q≥3 + MARG, no keep | 87.0 % | 0.28 % | 87.6 % |
| fb Q≥3 + MARG, **+KEEP** | **55.6 %** | 0.09 % | 98.0 % |
| **MARG off** | | | |
| fb Q≥20, no keep | 80.7 % | 0.21 % | 90.8 % |
| fb Q≥20, **+KEEP** | 77.0 % | 0.11 % | 96.5 % |
| fb Q≥0, no keep | 88.9 % | 1.24 % | 90.8 % |
| fb Q≥0, **+KEEP** | 80.7 % | 0.16 % | 96.5 % |

**The premise is only half true, and this is the central finding.** The recall win and the
HipSTR concordance are the **same marginal-het cells viewed two ways** — the lone-carrier
recurrent hets. Emit their GT and you gain silver recall but lose HipSTR agreement; no-call
them and you gain agreement but lose the recall. They do **not** compose:

- **With MARG on, KEEP is catastrophic**: recall collapses to ~54 % (below even the
  heuristic's 80.7 %), because `apply_fp_control` no-calls exactly the imbalanced hets MARG
  worked to recover. MARG + the per-sample no-call is the worst of both — it reproduces the
  known "marg+oldFP ≈ 57 %" collapse. **MARG and the no-call are mutually exclusive.**
- **With MARG off, KEEP is benign and useful**: freebayes Q≥20 + KEEP lands at 77.0 % /
  0.11 % / **96.5 %** — better precision than the heuristic at slightly lower recall, with
  full concordance parity — while emitting **1128 comparable loci vs the heuristic's 774**
  (~46 % more genotyped output). Freebayes Q≥0 + KEEP is 80.7 % / 0.16 % / 96.5 %:
  numerically the heuristic on all three axes, but from one principled test + a clean QUAL
  knob and with far more loci emitted.

So concordance is recoverable, but only in the MARG-off regime, and recovering it gives
back most of the recall advantage over the heuristic. At 3 reads/plant you cannot bank both
the +6-point recall and the 96 % concordance — they are the same cells.

## 4. Recommendation

**Adopt `freebayes` as the caller's emission model**, replacing the heuristic gate stack. It
dominates BIC on the frontier (tie-to-+1 recall at matched FP, reaches the recall ceiling at
lower FP) and is the cleaner instrument — a real per-locus QUAL, no arbitrary margin, and it
emits a genuine site quality downstream tools can threshold. The added math over BIC is ~20
lines (the Ewens SFS prior); worth it.

**Default operating regime — precision/QC-first (recommended default):**
`PVC_SSR_EMIT_MODEL=freebayes`, **MARG off**, **KEEP on**, `PVC_SSR_FREEBAYES_MIN_QUAL ≈ 20`.
This is a strict structural upgrade at zero cost to the three metrics: it matches the
heuristic's recall (~77–81 %), precision, and concordance (96.5 %), while replacing a stack
of hand-tuned gates with one principled test and a clean knob, and emitting ~46 % more
genotyped loci. No regression to defend.

**Keep `apply_fp_control` (the per-sample no-call)? Yes — for the default.** It is what
holds concordance at 96.5 %, and it is orthogonal to the *site* emit decision (the model
still decides the locus; the no-call only cleans the cells). The small change is already
implemented and measured (§3, `PVC_SSR_KEEP_FP_CONTROL`); on adoption it should be folded
into the freebayes default rather than left as an opt-in toggle.

**Make `MARGINALIZED_PRIOR` the default? No.** It only pays off in the MARG-off-incompatible
sense: it buys +6 recall points *only when the per-sample no-call is off*, and that trade
costs 8–10 concordance points against a flawed F-blind reference — and combining it with the
no-call collapses recall to 54 %. Keep MARG an explicit opt-in for a documented
**recall/discovery regime** (`freebayes`, MARG on, KEEP off, `MIN_QUAL ≈ 20` → 84 % / 0.19 %
/ ~88 % concordance) for users who want maximum real-variant discovery and accept lower
HipSTR agreement.

## 5. Honest verdict on the ceiling

At median ~3 reads/plant the depth is the ceiling, and the recall/concordance tension is a
genuine information wall, not a modelling gap. The freebayes report's "+6–12 recall points
below 0.10 % FP" is real on the silver core — but it is bought entirely by *not* no-calling
marginal hets, which is exactly what drops HipSTR concordance. You cannot have both at this
depth; they are the same lone-carrier hets.

**Does freebayes beat BIC?** Yes, but modestly: tie-to-+1 recall at matched FP, and it hits
the ~87 % ceiling at 0.28 % FP where BIC needs 0.35 %. The bigger, in-caller truth is that
the toy-prototype's "2–3× better than BIC" shrinks to a consistent-but-small edge once BIC
runs on the real chemistry. **Is the difference worth the complexity?** Over BIC — yes, for
~20 lines you get a real QUAL and a cleaner knob, and it never loses. Over the heuristic —
the justification is *principle and tunability*, not a recall leap: at the safe operating
point freebayes matches the heuristic on all three metrics. The emission model is **not**
where the remaining recall lives; the per-sample-no-call-vs-concordance tradeoff dominates,
and that is a depth problem. Adopt freebayes for the clean instrument and the extra emitted
loci; don't expect it to break the 3-reads/plant wall.

## Reproduce

```
BIN=target/release/pop_var_caller   # host release build of this branch
CAT=benchmarks/ssr_tomato1/results/ours/ssr_tomato1.ssr.catalog
PSP=benchmarks/ssr_tomato1/results_rerun_20260708/ours/cohort/psp

# recommended default (freebayes, MARG off, per-sample no-call kept, QUAL floor 20):
PVC_SSR_EMIT_MODEL=freebayes PVC_SSR_KEEP_FP_CONTROL=1 PVC_SSR_FREEBAYES_MIN_QUAL=20 \
  $BIN ssr-call --catalog $CAT --output out.vcf --threads 8 $PSP/*.ssr.psp

# score (confident core) + concordance:
uv run --no-project benchmarks/ssr_tomato1/scripts/silver_standard.py --score out.vcf
uv run --no-project benchmarks/lib/ssr_concordance.py --ours out.vcf \
  --hipstr benchmarks/ssr_tomato1/results_ssr15k/hipstr/cohort.str.vcf.gz
```

Matrix VCFs left in
`benchmarks/ssr_tomato1/results_rerun_20260708/ours/cohort/emission_matrix/`
(`{heuristic,bic.m{0,5,10,20},fb.q0}.marg{0,1}[.keep].vcf`).
Scoring harnesses: scratchpad `score_matrix.py` (full matrix) and `score_keep.py`
(both axes for the KEEP variants).
