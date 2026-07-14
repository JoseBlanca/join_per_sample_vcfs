# ng step 2 (generic path) — read preparation: left-align + BAQ

*Status: design spec, 2026-07-14. The **generic (SNP/indel) path** of step 2. Inherits the
shared contract and discipline from the preamble
[`read_preparation.md`](read_preparation.md) — read that first; this spec covers only what is
specific to the generic path. Grounded in the production `process_read` fold
([pileup/per_sample/read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs)).
**No code yet.***

---

## 1. Scope — goals and non-goals (path-specific)

**Goal:** prepare a generic (SNP/indel) read by **trusting the mapper's placement** but
canonicalising it — left-align its indels, cap its base qualities by alignment confidence — into
a still-decomposable read the pileup walker can turn into per-position evidence.

**Non-goals** (beyond the shared preamble's): this is the *trust-the-mapper* implementation. It
does **not** realign or reassemble; local reassembly is a separate, deferred `ReadPrep` sibling
(§6). It does not decompose the CIGAR into events — that is the pileup walker's job (§6).

The output type is **`PreparedReadNg`** and the consumer is the **pileup walker** (preamble §2).

---

## 2. The transform — two steps, ported from `process_read`

Production runs one per-read fold, `process_read`, whose stages are
`G2 bad-CIGAR → F3 left-align → F1 mismatch-fraction → BAQ`. ng already assigned the two
*rejects* — `G2` and `F1` — to **step 1** (they are filters #9 and #8 in `read_filtering.md`
§3). What remains for step 2 is the two *transforms*, in order:

1. **Indel left-alignment.** Shift every indel to its leftmost equivalent reference position and
   merge adjacent indels, so the same event placed at different offsets by the aligner collapses
   to one allele. This rewrites **only the CIGAR** — the read's bases and qualities are untouched.
   It is a GATK `leftAlignIndels` port wrapping a representation-neutral kernel that also serves
   the STR off-ladder path (§7).
2. **BAQ (base alignment quality).** A banded HMM re-aligns the read to the reference to estimate,
   per base, the probability it is *mis-aligned*, and caps each base quality at that confidence
   (`bq = min(base_quality, BAQ)`); bases the HMM cannot place are set to quality 0. BAQ
   down-weights bases in and near ambiguous indels **without rewriting the alignment**. It is an
   htslib `probaln_glocal` port. BAQ can decline a read outright (HMM overflow; reference window
   past the contig end; no aligned `M` op) → `None` (§4).

**Jargon, once — BAQ.** *Base Alignment Quality* is a per-base confidence that a base is aligned
to the right reference position (as opposed to the base-*call* quality, confidence in the letter).
A base sitting in an ambiguously-placed indel gets a low BAQ and is de-weighted, so a
mis-alignment cannot masquerade as a confident mismatch.

### Why the step-1/step-2 split is safe (a decision, and the alternative)

The obvious worry: production runs the mismatch-fraction filter (`F1`) *after* left-alignment
(`F3`), but ng runs mismatch filtering in step 1, *before* step-2 left-alignment. **This is
safe because left-alignment provably preserves the mismatch count** — a debug-assert in the
production `left_align_indels` guarantees the mismatch total is unchanged by the re-placement.
So ng's order gives the identical verdict to production's, and the bad-CIGAR check (`G2`/#9) sees
the raw decoded CIGAR in both. *Alternative considered:* keep the whole fold (both filters and
both transforms) in one step. Rejected — it would re-merge filtering and preparation that ng
deliberately separated (`read_filtering.md` §2.5), and the mismatch-count invariant makes the
split free.

---

## 3. The output type — `PreparedReadNg`

```rust
/// The generic prepared read — the output of the generic ReadPrep impl (`type Prepared`).
/// A still-*decomposable* read: the pileup walker turns it into per-position events. A NEW ng
/// type; the production `PreparedRead` name is not reused (ng_step_interfaces §6).
pub struct PreparedReadNg {
    pub chrom_id: ContigId,
    pub alignment_start: Position,
    pub alignment_end: Position,          // cached, so the walker never re-walks the CIGAR for span
    pub cigar: Vec<CigarOp>,              // LEFT-ALIGNED (unlike MappedRead.cigar)
    pub seq: Vec<u8>,                     // uppercase ACGTN
    pub bq_baq: Vec<u8>,                  // BAQ-capped min(BQ, BAQ) (unlike MappedRead.qual)
    pub mq_log_err: f64,                  // derived ln(P_err) from MAPQ (not on MappedRead)
    pub mapq: MapQual,                    // raw, preserved
    pub is_reverse_strand: bool,          // decoded from flag
    pub mate_role: MateRole,              // Solo | FirstOfPair | SecondOfPair (from flag bits)
    pub adaptor_boundary: Option<Position>, // carried through, applied later by the walker (§6)
}
```

What `PreparedReadNg` carries that a raw `MappedRead` does not: a **left-aligned** CIGAR,
**BAQ-capped** qualities, a precomputed alignment end, a precomputed `mq_log_err`, decoded strand
and mate role. It carries **no** per-base overlap adjustment — that is pairwise and happens in
the walker (§6), which is what keeps a `PreparedReadNg` self-contained and pairwise-independent.

The generic implementation:

```rust
pub struct TrustMapperPrep { /* config: BAQ on/off, thresholds */ }
impl ReadPrep for TrustMapperPrep {
    type Prepared = PreparedReadNg;
    fn prepare_read(&self, read: &MappedRead, window: &LocusWindow) -> Option<PreparedReadNg> { /* left-align → BAQ */ }
}
```

---

## 4. Reference dependency and the "no observation" reasons

**Reference (reuse the `RefSeq` split, `ref_seq.md`).** Production reads the reference two ways
here, and ng reuses exactly that split: **raw, case-preserving bytes** (`RawRefSeq`) for
left-alignment — matching the aligner's own view — and **canonical uppercased bytes** (`RefSeq`)
for the BAQ HMM. No new accessor.

**`None` reasons (tallied):** `Baq` — BAQ declined the read (HMM overflow, reference window past
the contig end, no aligned `M` op, an `N`/ref-skip in the CIGAR). With BAQ disabled by config,
the transform is left-align only and never returns `None`.

---

## 5. Cross-cutting concerns

- **Performance / parallelism.** Pairwise-independence (preamble §3) makes generic prep
  embarrassingly parallel with deterministic output. BAQ (the per-read HMM) is the cost; reuse a
  per-worker scratch buffer for its matrices rather than allocating per read (the project's
  scratch-buffer discipline). `PreparedReadNg` is owned per read; reusing its allocations is the
  same locus-stream-level question deferred in `read_filtering.md` §6, not a per-read decision.
- **Determinism.** No mate/column context here means the same read prepares to the same
  `PreparedReadNg` regardless of thread interleaving — a property downstream determinism relies on.

---

## 6. Deferred, with a recommended home

- **Adaptor-mask application, mate-overlap reconciliation, CIGAR decomposition → the pileup
  walker.** `PreparedReadNg` is *decomposable*, not decomposed; it *carries* the adaptor boundary
  but does not apply it (preamble §5).
- **Local haplotype reassembly → a future generic `ReadPrep` sibling.** GATK-style de-novo
  assembly + realign-to-haplotype is the "reassemble" pole of step 2's axis (`ng_proposal.md`
  §2) — a genuine bench competitor to `TrustMapperPrep`, added as a sibling `impl ReadPrep` when
  we test it, not part of v1.

**Alternatives to bench (recorded, not built now):** freebayes-style **trust + `stablyLeftAlign`
only** (no BAQ) — the same trust-the-mapper approach minus the BAQ cap; and GATK-style **local
reassembly** (above). Both are `ReadPrep` siblings producing a `PreparedReadNg`-shaped output, so
the pileup walker consumes any of them unchanged.

---

## 7. Reuse over rewrite — the map to production

The parity oracle is the production prepared read: a ported generic impl is correct when its
`PreparedReadNg` is byte-identical to the production `PreparedRead` on a fixture — same
left-aligned CIGAR, same BAQ-capped qualities.

| what | existing code | ng reuse |
|---|---|---|
| the per-read prep fold | `process_read` ([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs)) | model for `TrustMapperPrep` — **only** its F3 + BAQ stages (G2/F1 are step-1 filters) |
| indel left-alignment | `left_align_indels` / `left_align_cigar` ([pileup/walker/indel_norm.rs](../../../../src/pileup/walker/indel_norm.rs)) → `normalize_alleles` ([norm_seqs.rs](../../../../src/norm_seqs.rs)) | call directly (GATK `leftAlignIndels` port over a representation-neutral kernel) |
| BAQ | `BaqEngine::process` ([pileup/per_sample/baq_engine.rs](../../../../src/pileup/per_sample/baq_engine.rs)), htslib `probaln_glocal` | call directly; `None` on BAQ-skip |
| the prepared read | `PreparedRead` + `mapped_to_prepared` ([pileup/walker/mod.rs](../../../../src/pileup/walker/mod.rs)) | model for `PreparedReadNg` — **new ng type**, name not reused |
| reference | `RawRefSeq` (left-align) + `RefSeq` (BAQ) ([ref_seq.md](ref_seq.md)) | reuse as-is |

---

## 8. Open questions

- **Does `PreparedReadNg` equal the production `PreparedRead` exactly? — leaning yes.** The survey
  shows `PreparedRead` *is* precisely the F3+BAQ output the walker consumes. If a later generic
  impl (reassembly) needs a richer payload, `PreparedReadNg` widens then; v1 mirrors `PreparedRead`.
- **`LocusWindow` shape.** The generic impl needs only reference bases (raw + canonical) for its
  window; the exact type is co-owned with the router spec (`LocusKind`).
- **Where generic prep is *invoked* — leaning: by the pileup, per read, as it walks a non-STR
  stretch** (matching production, where `process_read` runs as the walker ingests reads). Confirm
  the call site when `pileup/` is specced; it is what makes the compose-not-subsume resolution
  (preamble §2) concrete.
