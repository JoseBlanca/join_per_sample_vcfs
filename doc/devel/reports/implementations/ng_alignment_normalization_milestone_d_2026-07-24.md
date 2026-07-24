# ng alignment — normalization (plan 3 of 3), Milestone D — the screen + its result

**Date:** 2026-07-24
**Plan:** [alignment_normalization.md](../../ng/impl_plan/alignment_normalization.md) — step D1 (**plan 3 complete**).
**Spec/arch:** [spec §6/§10.3](../../ng/spec/alignment.md), [arch §Test & bench shape](../../ng/arch/alignment.md).
**Status:** implemented → reviewed → fixes-applied. **Plan 3 (and the alignment module) complete.** At **Checkpoint D**.

---

## 1. Plan

D1 is the **differ-at-all screen** — the cheapest discriminating measurement in the normalization
comparison: run all three left-aligners over the same real reads and count, per pair, how many
outputs disagree. If near zero, normalization *placement* cannot explain a difference in calling and
the avenue closes for one run (spec §6). The screen supplies the measurement; adoption is a different
plan's call.

Added at the owner's request: a **synthetic stress-read generator** that engineers reads to make the
three disagree, so the screen's discriminating power is proven, not assumed.

## 2. Changes made

- [examples/ng_normalizer_screen.rs](../../../../examples/ng_normalizer_screen.rs): reads a reference
  (`.fai` for the contig table, `WindowedRefSeq` for per-read windows) and one or more BAM/CRAMs
  through ng's real ingestion (`SampleReads::reads_in_region` → `MappedRead`), builds an `Alignment`
  per indel-bearing read, normalizes it with 1a/1b/1c (1c under `catch_unwind`, since it panics), and
  tallies per-pair disagreements + cap-hits + panics + a `moved_by_normalization` counter.
- [examples/ng_synthesize_stress_reads.rs](../../../../examples/ng_synthesize_stress_reads.rs):
  writes a reference + coordinate-sorted BAM whose reads place a single indel at the rightmost end of
  an A-homopolymer of length L (shift `L−1` deletion / `L` insertion), lengths straddling 1b's 20-pass
  cap.

## 3. The result — recorded at Checkpoint D

**GIAB HG002 10x** (`benchmarks/ssr_hg002/bam/10x` + GRCh38, run on the host, 8.3 s):

| | count |
|---|---|
| reads seen | 942,502 |
| reads with an indel | 63,757 |
| moved by normalization | 6 |
| **1a vs 1b disagree** | **0** |
| **1a vs 1c disagree** | **0** |
| **1b vs 1c disagree** | **0** |
| 1b hit its cap | 0 |
| 1c panicked | 0 |

**Zero disagreements** across all real indel reads. Caveat: the BAM is novoalign-aligned and
novoalign pre-left-aligns, so only 6/63,757 reads were non-canonical — the GIAB data weakly *exercises*
the normalizers (a run on a non-left-aligning aligner's output, e.g. the tomato bwa CRAMs, would
harden the claim; the tomato reference is not present locally, so that run is deferred).

**Synthetic stress set** (64 reads engineered to disagree):

| | count |
|---|---|
| reads with an indel | 64 (100% moved) |
| **1a vs 1b disagree** | **28** |
| **1a vs 1c disagree** | **0** |
| **1b vs 1c disagree** | **28** |
| 1b hit its cap | 36 |
| 1c panicked | 0 |

So the screen **is** discriminating; the *only* way the three differ on well-formed reads is 1b's
20-pass cap on a shift >20 bp (1a/1c reach leftmost in one pass, 1b stops short). Real mappers do not
misplace indels that far — hence the GIAB 0.

## 4. Interpretation (bounded by spec §6)

- **The normalizer choice does not change calling on real reads** — normalization placement is **not**
  the lever behind production's indel deficit. That avenue closes; the reference-favouring genotype
  prior and the filters (outside this module) remain the candidates.
- **Recommended normalizer: 1a** (`StructuredLeftAligner`, the GATK/production port) — always leftmost
  (like 1c, unlike 1b), cheapest (one pass, no loop, no panic). **1c ≡ 1a** confirmed (0 disagreements,
  0 panics on 63,821 reads — 1a is a true fixpoint). **1b** (freebayes) is the one that can be wrong
  (the cap). 1b/1c are retained as comparators / regression-guards, not production defaults (owner
  decision, 2026-07-24). The **1b complex-D/I trim** flagged at Checkpoint C is left as-is: 1b is not
  the default, and complex D/I never survives the filters on real reads (0 in GIAB).

## 5. Validation

- `cargo fmt --check` clean · `cargo clippy --all-targets --all-features -- -D warnings` clean ·
  `cargo test --example ng_normalizer_screen` 5/5 · `cargo test --example ng_synthesize_stress_reads`
  3/3 · `cargo test --lib` 2331 passed (container); both examples built + run on the host.

## 6. Review and fixes

- Review [ng_normalizer_d1](../reviews/ng_normalizer_d1_2026-07-24.md): 1 Bl / 2 Maj / 3 Min — the
  measurement logic and generator confirmed correct (numbers reproduced by hand); all findings were
  test-coverage or diagnostics. Fixes [fixes_applied_ng_normalizer_d1](../reviews/fixes_applied_ng_normalizer_d1_2026-07-24.md).

## 7. Follow-ups

- A stronger real-reads run on a non-left-aligning aligner's output (tomato bwa CRAMs) once their
  reference is reachable — would move `moved_by_normalization` well above 6 and re-confirm 0
  disagreement with the normalizers genuinely exercised.
- Adoption of 1a on the generic read-preparation path — read-preparation's plan (D1's result is its input).
