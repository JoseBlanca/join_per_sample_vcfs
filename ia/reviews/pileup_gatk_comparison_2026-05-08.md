# Algorithm Review: pileup walker vs. GATK

**Date:** 2026-05-08
**Reviewer:** Claude (algorithm-comparison study)
**Module reviewed:** `src/per_sample_caller/pileup` plus the
`cram_input` filter cascade that gates it.
**Reference codebase:** `gatk/` (in-tree copy of the
Hellbender / GATK4 Java source).
**Status:** Advisory — open backlog to triage. Companion to
[pileup_samtools_comparison_2026-05-07.md](pileup_samtools_comparison_2026-05-07.md)
and [pileup_freebayes_comparison_2026-05-08.md](pileup_freebayes_comparison_2026-05-08.md).
GATK is the most distant relative of the three: its production
caller (HaplotypeCaller) is assembly-based, and its pileup layer
serves ancillary tools rather than the main calling path. Lessons
are correspondingly more "implementation hygiene from a mature
codebase" than "core algorithmic ancestry."

---

## 1. Scope

This is **not** a defect-finding review. It is a comparative
study of GATK's pileup-equivalent layer against our Stage 1
walker, looking for techniques worth importing — and equally,
places where our model is already cleaner than GATK's so we
don't accidentally regress.

In-scope project files:

- [src/per_sample_caller/pileup/mod.rs](../../src/per_sample_caller/pileup/mod.rs)
- [src/per_sample_caller/pileup/walker.rs](../../src/per_sample_caller/pileup/walker.rs)
- [src/per_sample_caller/pileup/active_set.rs](../../src/per_sample_caller/pileup/active_set.rs)
- [src/per_sample_caller/pileup/cigar_cursor.rs](../../src/per_sample_caller/pileup/cigar_cursor.rs)
- [src/per_sample_caller/pileup/decompose.rs](../../src/per_sample_caller/pileup/decompose.rs)
- [src/per_sample_caller/pileup/open_record.rs](../../src/per_sample_caller/pileup/open_record.rs)
- [src/per_sample_caller/cram_input.rs](../../src/per_sample_caller/cram_input.rs)

GATK files studied (Java; paths under
`gatk/src/main/java/org/broadinstitute/hellbender/`):

- `utils/locusiterator/LocusIteratorByState.java` — the
  per-position walker (analogous to our `walker.rs`).
- `utils/locusiterator/AlignmentStateMachine.java` — per-read
  CIGAR cursor (analogous to our `CigarCursor`).
- `utils/locusiterator/{ReadStateManager, PerSampleReadStateManager}.java`
  — active-read containers, per-sample partitioning.
- `utils/pileup/{PileupElement, ReadPileup, PileupBasedAlleles}.java`
  — per-locus output objects.
- `utils/baq/BAQ.java` — GATK's port of Heng Li's BAQ HMM.
- `utils/downsampling/{LevelingDownsampler, AlleleBiasedDownsamplingUtils}.java`
  and 10 sibling files — the GATK-specific downsampling layer.
- `engine/filters/ReadFilterLibrary.java` plus 36 sibling
  filter classes.
- `utils/read/{ReadUtils, CigarUtils}.java` — supporting helpers
  (adaptor-boundary computation, "good CIGAR" predicate).

Out of scope:

- `tools/walkers/haplotypecaller/` — assembly-based caller, with
  local de Bruijn reassembly + PairHMM scoring. We explicitly
  rejected this design (Option C in
  [calling_pipeline_architecture.md §"Per-read likelihood
  quality"](../specs/calling_pipeline_architecture.md)). The
  active-region detector (`PileupBasedAlleles`,
  `ReferenceConfidenceModel`, `IsActive`) is bounded by the
  same rejection.
- Mutect2 — somatic-specific.
- BQSR (`BaseRecalibrator`) — per-base recalibration; upstream
  of our pipeline. We accept whatever quality scores reach
  Stage 1 and let BAQ refine them.
- VQSR / variant filtering — post-call.
- Spark variants (`LocusWalkerSpark` etc.) — implementation
  framework, not algorithm.

## 2. Summary verdict

GATK's pileup layer is well-engineered but serves a different
output contract than ours. The three useful findings are
implementation-hygiene items GATK has refined over years:
**adaptor-region per-base filtering**, a **CIGAR-wellformedness
read filter**, and a **BAQ parameter reference** for whenever
we implement BAQ. The six rejections are mostly "GATK serves a
different contract" issues — per-(read, position) `PileupElement`
allocation, skip-empty-positions emission, post-pileup
quality-mutation for mate overlap, allele-aware downsampling
for somatic, HaplotypeCaller's reassembly path itself, and a
soft-clip-ratio filter that duplicates the work MAPQ + F1
already do.

The single most reassuring finding is `A1`: GATK's
mate-overlap math
([ReadPileup.java:327-360](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/pileup/ReadPileup.java#L327-L360))
is byte-identical to samtools' and to what we adopted in `S7`
of the samtools review — sum-and-cap on agreement,
0.8-haircut on disagreement, zero-the-loser on both.
Three independent codebases converging on the same formula
is as close as this kind of work gets to validation.

## 3. Where we already match or beat GATK

These are noted up-front so we don't accidentally regress
them while acting on the findings below.

- **`A1` — Mate-overlap math.**
  GATK's `fixPairOverlappingQualities`
  ([ReadPileup.java:327-360](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/pileup/ReadPileup.java#L327-L360))
  applies the same formula we use in our walker:
  agreeing pairs sum quals capped at the SAM max
  (`MAX_SAM_QUAL_SCORE`); disagreeing pairs scale the
  higher-qual side by `SAMTOOLS_OVERLAP_LOW_CONFIDENCE = 0.8`
  ([line 32](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/pileup/ReadPileup.java#L32))
  and zero the loser. Identical to what we adopted in `S7`
  ([pileup_samtools_comparison_2026-05-07.md `S7`](pileup_samtools_comparison_2026-05-07.md#s7--adopt-samtools-bq-combining-math-for-match-only-mate-overlap)).
  Three independent implementations (samtools, GATK, ours)
  converge on the same math.

- **`A2` — Stage 1 indel left-alignment.**
  GATK does **not** left-align indels at the pileup layer.
  Their pipeline relies on either upstream tooling
  (`LeftAlignAndTrimVariants` runs *post*-call, on VCFs) or
  HaplotypeCaller's local reassembly to converge ambiguous
  homopolymer-context indel placements. For a pileup-shaped
  output like ours, that path is unavailable — our `F3`
  left-alignment in `cram_input`
  ([commit `3b9e579`](../../src/per_sample_caller/cram_input.rs))
  is the right move and isn't mirrored by GATK at the
  same layer.

- **`A3` — Pre-aggregated per-allele scalars.**
  GATK's `PileupElement`
  ([PileupElement.java:20-63](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/pileup/PileupElement.java#L20-L63))
  carries a single (read, base, qual, CIGAR-state) tuple per
  (read, locus) pair. Depth at a position is a
  `list.size()` over its `ReadPileup`
  ([ReadPileup.java:258-260](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/pileup/ReadPileup.java#L258-L260));
  per-allele counts come from `getBaseCounts()`
  ([ReadPileup.java:281-295](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/pileup/ReadPileup.java#L281-L295))
  iterating the list every time. Our walker pre-aggregates the
  five scalars per allele at fold time
  ([calling_pipeline_architecture.md §"The five per-allele
  scalars"](../specs/calling_pipeline_architecture.md)), so
  Stage 5 reconstructs likelihoods without re-iterating
  per-read state. The right choice for the `.psf` artefact
  shape; don't regress to per-read storage.

- **`A4` — Emit at every covered position natively.**
  `LocusIteratorByState` only yields an `AlignmentContext`
  when the pileup is non-empty
  ([LocusIteratorByState.java:322-324](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/locusiterator/LocusIteratorByState.java#L322-L324));
  filling in empty-coverage positions requires wrapping with
  `IntervalAlignmentContextIterator`. Our walker emits a
  `PileupRecord` for every covered position unconditionally
  ([calling_pipeline_architecture.md §"Stage 1 — per-sample
  caller"](../specs/calling_pipeline_architecture.md)), which
  is what the `.psf`'s `delta_pos` encoding expects.
  Native fit; don't regress.

## 4. Findings — backlog to act on

Each finding has a stable id (`Gi` for "GATK-inspired") so we
can reference them in commits and a one-by-one work plan.

Priority key (same as the prior reviews):
- **High** — clear user-visible benefit, low risk.
- **Medium** — usability or hygiene; nice to have.
- **Low** — only matters under conditions we don't currently
  target, or speculative.

### `G1` — Adaptor-region per-base filter

- **Priority:** Medium → High once paired-end data is the
  primary input.
- **Effort:** Small (compute adaptor boundary on
  `PreparedRead` admission; one branch in the cursor's
  Match-event emission).

**Observation.** The single per-base filter `LocusIteratorByState`
applies inside the iterator is `ReadUtils.isBaseInsideAdaptor`
([LocusIteratorByState.java:337-339](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/locusiterator/LocusIteratorByState.java#L337-L339)):

```java
private static boolean dontIncludeReadInPileup(final GATKRead rec, final long pos) {
    return ReadUtils.isBaseInsideAdaptor(rec, pos);
}
```

The boundary is computed from mate insert size
([ReadUtils.java:1066-1072](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/read/ReadUtils.java#L1066-L1072)):
when a paired read's alignment runs *past* its mate's start
(forward strand) or *before* its mate's end (reverse strand),
the bases beyond the boundary were sequenced *through the
adaptor* on the flow cell. The aligner often soft-clips
adaptor sequences but not always — when it places them as
`M`-op bases, they are noise mapped to spurious reference
positions.

We have nothing equivalent. Our `F1` mismatch-fraction filter
catches some adaptor-readthrough cases (a 30%-mismatch tail
trips the threshold), but only at whole-read granularity. A
read with a clean prefix and a 50bp adaptor tail passes `F1`
overall (1.5% mismatch fraction across 150 bp) while still
contributing 50 bp of adaptor noise.

**Proposal.**

1. Compute the adaptor boundary for paired reads at
   `cram_input` admission time, attach it to `PreparedRead`
   (one optional `u32` field).
2. In the `CigarCursor`'s Match-event emit sites
   ([cigar_cursor.rs](../../src/per_sample_caller/pileup/cigar_cursor.rs)),
   skip events whose `ref_pos` falls past the boundary. Same
   pattern as the `F5` read-`N` skip — skip at emit, not at
   fold.
3. Single-end reads and pairs without computable boundaries
   (mate unmapped, weird insert size) are unaffected.
4. Surface a counter in `RunSummary`.

**Rationale.** Cheap defence against a class of false positives
that BAQ doesn't address (BAQ adjusts confidence per base
within an alignment; it does not detect adaptor-runthrough
placement).

**Risk.** Low. The boundary computation is well-understood;
GATK's implementation is the reference. The skip-at-emit
pattern is the same we use for `F5`.

### `G2` — `GoodCigar`-style read-level rejection

- **Priority:** Medium
- **Effort:** Small (one predicate, one branch in the
  `cram_input` filter cascade)

**Observation.** GATK's
[`GoodCigarReadFilter`](../../gatk/src/main/java/org/broadinstitute/hellbender/engine/filters/ReadFilterLibrary.java#L55-L58)
rejects whole reads whose CIGAR is "ill-formed" by the rule in
[`CigarUtils.isGood`](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/read/CigarUtils.java#L96-L106):

```java
//  - has no consecutive I/D elements
//  - does not start or end with deletions (with or without preceding clips).
```

Both conditions are sentinel-grade alignment artefacts: no
biological event produces an immediate I-then-D in a single
read, and a deletion at a read boundary has no flanking
evidence on the missing side.

Our coverage of the "deletion at boundary" case is partial. Our
`decompose`/`CigarCursor` drops indel *events* at first/last
CIGAR op
([decompose.rs:106](../../src/per_sample_caller/pileup/decompose.rs#L106),
[decompose.rs:120](../../src/per_sample_caller/pileup/decompose.rs#L120)),
but the read itself stays in the active set and contributes
matches. Our `F3` left-alignment can shift an interior
deletion *to* a boundary, where the existing rule then drops
it. So the deletion event is dealt with — but we don't catch
the "consecutive I/D" case at all, and we don't reject the
whole read for either pattern.

**Proposal.**

1. Add a `cigar_is_good` predicate to `cram_input` (or its
   `PreparedRead` admission path):
   - Reject if the CIGAR contains any adjacent `I`/`D` pair
     (in either order).
   - Reject if the CIGAR's first or last op (ignoring leading
     soft/hard clips) is a deletion.
2. Surface a counter in `FilterCounts.bad_cigar`.
3. Default-on, no opt-out — same posture as `F3`. Both rules
   fire on alignments that no biological event produced.

**Rationale.** Defence-in-depth at the input boundary. The
"consecutive I/D" case is a small fraction of real data but
when it fires, the alignment is genuinely confused; keeping
the read produces noise events that BAQ and our F1 mismatch
filter don't catch.

**Risk.** Low. The rejection criteria are objective and
sentinel-grade. Tests should pin both branches (adjacent I/D,
boundary D) plus the negative case (normal CIGAR passes).

### `G3` — Soft-clip ratio read filter

**Rejected after discussion** (see §5 for the reasoning). Kept
in this section as a placeholder so the `Gi` numbering stays
stable across cross-references in commits and prior conversation;
the rejection rationale lives in §5.

### `G4` — When we implement BAQ, mirror samtools / GATK parameters

- **Priority:** Low (documentation only, until BAQ implementation lands)
- **Effort:** Tiny (one paragraph in
  [calling_pipeline_architecture.md §"Per-read likelihood
  quality"](../specs/calling_pipeline_architecture.md))

**Observation.** GATK's
[`BAQ.java`](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/baq/BAQ.java)
is a port of Heng Li's HMM-Glocal BAQ algorithm with the
following defaults
([BAQ.java:76-78](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/baq/BAQ.java#L76-L78)
plus
[BAQ.java:88-142](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/baq/BAQ.java#L88-L142)):

- `DEFAULT_GOP = 40` (gap-open penalty, Phred-scaled — i.e.
  P(gap) = 1e-4).
- `DEFAULT_BANDWIDTH = 7` (DP band half-width).
- `e = 0.1` (gap extension probability).
- `minBaseQual = 4` (any base below Q4 is raised to Q4 before
  the HMM, to prevent very-low-quality bases from dominating).

The code is explicitly "synchronized with samtools repo"
([BAQ.java:174-175](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/baq/BAQ.java#L174-L175)).
Two independent ports converge on the same parameters.

We don't implement BAQ yet. When we do, port from samtools' or
GATK's reference rather than re-derive — the parameters are
empirically calibrated and have been stable for ~15 years.

**Proposal.** Add a sentence to
[calling_pipeline_architecture.md §"Why BAQ in-process"](../specs/calling_pipeline_architecture.md)
naming the four parameter values and stating that the future
implementation should match them, citing both samtools' source
and GATK's port for cross-reference.

**Rationale.** Cheap insurance against a future engineer
re-deriving the parameters and getting it slightly wrong.

**Risk.** None — documentation only.

## 5. Findings explicitly *rejected*

For the record, so we don't relitigate them:

- **`R1` — Per-(read, position) `PileupElement` allocation.**
  Rejected.
  [`PileupElement`](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/pileup/PileupElement.java)
  is a heavy object (read reference + offset + CIGAR state +
  qual fields, ~40 bytes plus the read pointer). At 30×
  coverage on a 50 kb window that's 1.5M PileupElements
  allocated and freed during a walker pass. The lazy_cigar
  refactor
  ([pileup_lazy_cigar_2026-05-07.md](../reports/implementations/pileup_lazy_cigar_2026-05-07.md))
  measured the eager-event-vector approach as the dominant
  memory cost (~178× speedup at L=5000 after dropping it).
  Re-introducing per-position objects on the hot path would
  undo that win for unclear benefit. Our `CigarCursor` queries
  on demand and does not materialise per-position state.

- **`R2` — Skip emission at non-variant positions.**
  Rejected. `LocusIteratorByState`
  ([line 322-324](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/locusiterator/LocusIteratorByState.java#L322-L324))
  yields no `AlignmentContext` for empty positions; getting
  every covered position requires wrapping with
  `IntervalAlignmentContextIterator`. Our `.psf` requires
  every covered position natively (the `delta_pos` encoding
  in [Stage 2](../specs/calling_pipeline_architecture.md)
  represents gaps implicitly, but only by counting positions
  the walker actually emitted), so emitting unconditionally
  is the right native shape. Don't import GATK's
  empty-skip-by-default pattern.

- **`R3` — `AlleleBiasedDownsamplingUtils` for somatic calling.**
  Rejected.
  [`AlleleBiasedDownsamplingUtils`](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/downsampling/AlleleBiasedDownsamplingUtils.java)
  preferentially downsamples reads supporting specific alleles
  for Mutect2's contamination-correction workflow. It is
  *not* integrated into the default `LocusIteratorByState`
  path — GATK's default downsampler is the allele-blind
  [`LevelingDownsampler`](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/downsampling/LevelingDownsampler.java),
  which is the same pattern we adopted in
  [`S5`](pileup_samtools_comparison_2026-05-07.md#s5--per-column-depth-cap-originally-framed-as-per-allele).
  Allele-aware downsampling would bias allele-frequency
  estimates in germline cohort calling, which is what we
  rejected in `S5`'s revision history. Don't pull forward.

- **`R4` — Mate-overlap as post-pileup mutation.**
  Rejected. GATK's
  [`fixOverlaps`](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/pileup/ReadPileup.java#L313-L319)
  mutates the read's `byte[]` qual array *in place* after
  pileup construction. The math is what we adopted in `S7`,
  but the mutation pattern is destructive — a downstream
  consumer holding the same `GATKRead` reference sees the
  modified quals. Our walker resolves mate overlap at fold
  time using a per-position decision; the `PreparedRead`'s
  qual buffer is never mutated. Cleaner for a streaming
  walker; don't import the post-construction mutation pattern.

- **`R5` — HaplotypeCaller's local reassembly + PairHMM.**
  Rejected. Already explicitly rejected in
  [calling_pipeline_architecture.md §"Per-read likelihood
  quality"](../specs/calling_pipeline_architecture.md) as
  Option C. GATK builds a de Bruijn graph of haplotypes from
  reads in an active region, scores each read against each
  haplotype with PairHMM, and runs a Bayesian model on the
  result. At our 2-10× target coverage it is both
  computationally infeasible and statistically unsound (too
  few reads per haplotype for reassembly to propose reliable
  novel haplotypes). Documented for completeness; the
  active-region detection that triggers it (`PileupBasedAlleles`,
  `IsActive`, `ReferenceConfidenceModel`) is bounded by the
  same rejection.

- **`G3` — Soft-clip ratio read filter.** Rejected after
  discussion (originally proposed in §4; relocated here so the
  `Gi` numbering stays stable across earlier commits and
  references). Two reasons, both decisive:

  1. **No principled default exists.** GATK ships
     [`SoftClippedReadFilter`](../../gatk/src/main/java/org/broadinstitute/hellbender/engine/filters/SoftClippedReadFilter.java)
     without a default value. The "right" threshold depends on
     read length, library prep, sequencing technology, and
     organism (ancient-DNA libraries with very fragmented
     molecules can legitimately produce many high-soft-clip
     reads, for example). A tunable that no user can calibrate
     with confidence is a footgun: shipped default-off, most
     users won't enable it; shipped default-on, the threshold
     misfires on legitimate data classes we care about.

  2. **The signal duplicates work other filters already do.**
     Soft-clipping is the aligner's confidence call. Filtering
     post-hoc on the soft-clip ratio amounts to second-guessing
     that call, which is a mapper-trust issue: if you don't
     trust the aligner's judgement about which bases to keep,
     the principled response is to use a different aligner or
     pre-process the input — not to layer heuristics on its
     output. MAPQ already encodes the aligner's overall
     confidence, `F1` already catches the case where the kept
     M-portion is itself noisy, and `G2` already catches the
     unambiguous CIGAR pathologies. Soft-clip ratio sits in
     the middle as a fuzzy signal duplicating what those
     filters already cover.

  Revisit only if real-data analysis shows
  `MAPQ ≥ 20` + `F1` + `G2` letting through a class of false
  positives that a soft-clip-ratio filter would catch.

## 6. Why we keep our model

The samtools and freebayes reviews each have a §6 explaining
why we don't import the rival's *data structure* even where it
looks tempting. The GATK review's §6 has a different shape:
GATK's pileup layer and our walker serve **different output
contracts at different points in their respective pipelines**.

**GATK's pileup is intermediate state.** Tools like
DepthOfCoverage, PileupSpark, and CollectAllelicCounts consume
it directly. HaplotypeCaller — GATK's flagship caller —
*bypasses* it, driving its assembly engine off raw reads in an
active region instead. The pileup layer is therefore an
ancillary tool's data structure, not the artefact the calling
pipeline rests on.

**Our pileup *is* the artefact.** The `.psf` is the durable
per-sample output, written to disk and consumed by every
downstream stage without re-touching the BAM. Per-allele
scalars are pre-aggregated because Stage 5 reconstructs
likelihoods from them. Phase chain ids span positions because
Stage 5 needs cross-position linkage to call compound
haplotypes without reads.

These two contracts diverge in everything that follows:

- **Allocation:** GATK can afford per-(read, position)
  `PileupElement` objects because they live and die within a
  single tool's iteration. Our walker can't, because the
  per-position output flows through a hot path on the way to
  disk.
- **Aggregation:** GATK's lazy `getBaseCounts()` is correct for
  ad-hoc queries. Our pre-aggregation is correct for "every
  position writes five scalars per allele to disk."
- **Empty-position emission:** GATK skips them by default
  because most tool consumers want only positions of interest.
  We emit them because the merger and posterior engine need
  the complete coverage shape.
- **Mate-overlap mutation:** GATK can mutate the read's qual
  array because the pileup is a transient in-memory view;
  the read object is reconstructed each tool invocation. We
  hold the `PreparedRead` once and walk it; mutating its
  bytes would couple unrelated stages.

So we import GATK's *individual ideas* (the mate-overlap math,
the adaptor-region check, the GoodCigar predicate) without
importing the *shape* of their pileup. That's the freebayes
review's §6 conclusion in different clothes: pick what's
useful from a related pipeline without adopting design choices
that serve a contract we don't share.

The single place where GATK's design *would* be a better fit
than ours is if we ever decided to emit only per-variant
positions to disk (giving up the lossless `.psf`). That
decision rests on the cohort-recall property the architecture
spec commits to in
[constraint 5](../specs/calling_pipeline_architecture.md), so
we keep our model.

## 7. Suggested execution order

Each finding is independent and can be its own small commit.
None of them block each other; pick whichever is most useful
on the current sample backlog.

1. ✅ `G1` — adaptor-region per-base filter, landed in commit
   `864dd98`.
2. ✅ `G2` — `GoodCigar`-style read filter, landed in commit
   `a3c441d`.
3. ✅ `G4` — BAQ parameter reference, documented in
   [calling_pipeline_architecture.md §"Why BAQ in-process"](../specs/calling_pipeline_architecture.md).
   Names samtools and GATK as the two reference ports, plus
   the four parameter values, so a future BAQ implementation
   can port from one and cross-check against the other.
4. ❌ `G3` — soft-clip ratio filter, rejected. See §5.

If profiling later shows real-data hot spots, the place to
look is *not* GATK — their pileup is allocation-heavier than
ours by design, and the algorithmic ideas worth importing are
the three remaining findings above, not the data structures.
