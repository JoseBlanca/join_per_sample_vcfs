# Investigation: the SSR delimiter collapses multi-unit alleles (gap-penalty root cause)

**Date:** 2026-06-24
**Branch:** `ssr-cohort`
**Trigger:** the end-to-end BAM→VCF test surfaced a `CA×10` allele collapsing to the
reference `CA×8`; first attributed to window truncation, then traced deeper.

## 1. Symptom

The Stage-1 delimiter ([`src/ssr/pileup/alignment.rs`](../../../../src/ssr/pileup/alignment.rs)
`delimit_read`) extracts a read's repeat tract by a Viterbi pair-HMM alignment to the
locus's reference haplotype, reading the tract off the two flank junctions. With **clean,
full-flank input** against a `CA×8` reference it extracts correctly from `CA×4` up to `CA×9`,
but **any allele ≥ ref+2 units collapses to exactly the reference length**:

```
clean GGGGGG + CA×k + TTTTTT  →  extracted units
  k = 4..9   →  k        (correct; −4 … +1)
  k = 10..16 →  8        (collapsed to the reference)
```

A short allele (`CA×6`, −2 units) extracts fine; the cliff is one-sided (insertions).

## 2. Cause

The collapse scales monotonically with the HMM gap-open probability (sweep on the same
`CA×8` locus, full flanks):

| `GAP_OPEN_PROB` | longest allele recovered |
|---|---|
| `2.9e-5` (current) | ref **+1** (`CA×9`) |
| `1e-3` | ref **+2** |
| `1e-2` | ref **+5** |

The delimiter uses **one uniform affine gap** (`GAP_OPEN_PROB = 2.9e-5`, `GAP_EXTEND ≈ e⁻¹`)
across the *entire* read↔reference alignment. A long allele needs a multi-unit insertion of
read bases against the reference tract; the affine gap charges a fixed `gap_open` (≈ −10.4
nats) plus per-base `gap_extend`, so for ≥ +2 units the insertion path loses to an
alternative path that absorbs the extra units as flank mismatches — collapsing the tract to
the reference length. STR length variation *is* multi-unit tract indels, so a flank-grade
affine gap is the wrong model for the tract.

## 3. What the reference tools do (`../pop_var_caller/`)

**HipSTR** (`HipSTR/src/`) — the tool our architecture doc follows — separates the two
regions explicitly:

- **Flanks + base substitutions:** a Dindel affine gap. `AlignmentModel.cpp` literally sets
  `dindel_probs[0] = 2.9e-5` — **our constant is HipSTR's *flank* gap-open**, scaled by local
  homopolymer length (`LOG_MATCH_TO_INS[homop_len]`).
- **The STR tract:** a dedicated **repeat block** (`RepeatBlock`, `RepeatStutterInfo`,
  `StutterAlignerClass`). Length changes are *not* affine-gapped; they are scored by a
  **stutter model** (`stutter_model.cpp::log_stutter_pmf`): a geometric distribution in
  **repeat units** — `P(+n units) = up · geom · (1−geom)^(n−1)` (in-frame), with a separate
  bp-level geometric for out-of-frame, capped at ±6 units (`MAX_STUTTER_REPEAT_INS`). A whole
  unit step is **one event**, not `period` base insertions — there is no per-base `gap_open`
  inside the tract.

So HipSTR borrowed the same `2.9e-5`, but **only for the flanks**; the tract gets unit-level
stutter scoring. We borrowed it and applied it to the whole alignment, including the tract.

**GangSTR** (`GangSTR/src/`) confirms the principle through a different architecture:
read-class likelihoods (enclosing / flanking / fully-repetitive) with a stutter/slippage
model — STR length is a stutter distribution, never a generic affine gap.

## 4. Why this matters for *our* delimiter specifically

Our Mark-2 delimiter's job is **extraction**, not genotyping: it must read each spanning
read's *observed* tract bytes so Stage-2 (`ssr-call`) can genotype them under its own stutter
model. If the delimiter collapses `CA×10 → CA×8` at extraction time, the `CA×10` evidence
never reaches Stage 2 — the stutter model downstream can't recover what extraction threw
away. So the extraction alignment must place the flank junctions correctly even when the read
tract is many units longer than the reference. That requires tract-internal length changes to
be cheap; the flanks must stay stiff so they still anchor.

## 5. Options

- **A — Tract-aware gap (recommended; HipSTR-aligned).** In the Viterbi, use a cheap gap for
  reference columns *inside* the tract (between the left/right flank junctions, which the
  delimiter already knows from `left_len`/`right_len`), keeping the stiff Dindel gap for flank
  columns. The read's tract then aligns fully (length difference absorbed as cheap tract
  indels) while the flanks anchor. The tract gap rate is a calibration constant tied to the
  slippage/stutter rate. A refinement makes whole-`period` steps cheapest (true stutter); for
  a *delimiter*, a per-base cheap gap in the tract is enough to avoid the collapse.
- **B — Globally soften `GAP_OPEN_PROB`.** One-line, but softens the flanks too (HipSTR keeps
  them deliberately stiff) — junction-misplacement risk in noisy flanks. A blunt instrument,
  not the real fix.
- **C — Length-dependent (full Dindel table) gap.** Doesn't help di-/tri-nucleotide motifs
  (Dindel scales by *homopolymer* length = 1 for `CACACA`), which is exactly why HipSTR needs
  the separate repeat block. Insufficient alone.

## 6. Recommendation

**Option A.** Make the per-reference-column gap penalty tract-aware: cheap inside the repeat
tract (a slippage-rate constant, ultimately calibrated with the Stage-2 stutter model), stiff
in the flanks (keep the Dindel `2.9e-5`). This fixes the collapse at the root, matches the
canonical tool's architecture, and keeps the delimiter's flank anchoring intact. The
committed window-recovery infra (`feaa9ef`) becomes its complementary defence — once the
delimiter can score long alleles, widening the window genuinely recovers the
mapper-misaligned all-Match case.

Open calibration question (defer to real data): the exact tract gap rate, and whether to
model whole-unit steps (period-aligned) explicitly vs a flat cheap per-base tract gap. Both
are far better than the current uniform stiff gap; the per-base cheap tract gap is the
minimal first step.
