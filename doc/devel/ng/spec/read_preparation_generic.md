# ng step 2 (generic path) — read preparation: left-align + optional BAQ

*Status: design spec, 2026-07-14. The **generic (SNP/indel) path** of step 2. Inherits the
shared contract and discipline from the preamble
[`read_preparation.md`](read_preparation.md) — read that first; this spec covers only what is
specific to the generic path. Grounded in the production `process_read` fold
([pileup/per_sample/read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs)).
**No code yet.***

---

## 1. Scope — goals and non-goals (path-specific)

**Goal:** prepare a generic (SNP/indel) read by **trusting the mapper's placement** but
canonicalising it — left-align its indels and, optionally, cap its base qualities by alignment
confidence — into a still-decomposable read the pileup walker can turn into per-position evidence.

**Both modes ship in v1**, not one now and one later:

- **left-align + BAQ** — the production default;
- **left-align only** (BAQ off) — the freebayes-style preparation. This is a real run mode we
  need, not a hypothetical bench alternative.

BAQ is a **config toggle on the one preparer**, not a second implementation (§2).

**Non-goals** (beyond the shared preamble's): this is the *trust-the-mapper* implementation. It
does **not** realign or reassemble; local reassembly is a separate, deferred `ReadPreparer` sibling
(§6). It does not decompose the CIGAR into events — that is the pileup walker's job (§6).

The output type is **`PreparedRead`** and the consumer is the **pileup walker** (preamble §2).

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
2. **BAQ (base alignment quality) — optional.** A banded HMM re-aligns the read to the reference
   to estimate, per base, the probability it is *mis-aligned*, and caps each base quality at that
   confidence (`bq = min(base_quality, BAQ)`); bases the HMM cannot place are set to quality 0.
   BAQ down-weights bases in and near ambiguous indels **without rewriting the alignment**. It is
   an htslib `probaln_glocal` port. BAQ can decline a read outright (HMM overflow; reference
   window past the contig end; no aligned `M` op) → `None` (§4). **With BAQ off the transform is
   step 1 only** — the raw qualities pass through uncapped and the preparer never returns `None`.

**Jargon, once — BAQ.** *Base Alignment Quality* is a per-base confidence that a base is aligned
to the right reference position (as opposed to the base-*call* quality, confidence in the letter).
A base sitting in an ambiguously-placed indel gets a low BAQ and is de-weighted, so a
mis-alignment cannot masquerade as a confident mismatch.

### BAQ is a toggle, not a second implementation (a decision, and the alternative)

Left-align-only is a **first-class v1 mode** — the freebayes-style preparation, needed for real
runs — and it is expressed as a **config toggle** (`baq: Option<BaqConfig>`) on the one preparer,
not as a sibling `LeftAlignOnlyPreparer`. Production already models it exactly this way: a
`--no-baq` path whose `prepare_passthrough` copies the raw qualities into `bq_baq` uncapped,
sharing the same fold.

*Alternative considered:* two sibling impls, so a recipe names "freebayes-style" vs "ours-SNP" as
an implementation swap rather than a config change. **Rejected** — the two differ *only* by
whether the final cap is applied, so a second preparer would buy one bake-off row at the price of
two code paths for one algorithm, and it would depart from the production shape we are porting.
The step-1 precedent is identical: `ReadFilterConfig`'s `max_read_mismatch_fraction: None` turns
filter #8 off without a second `ReadFilter` type.

*Consequence for the name:* `LeftAlignBaqPreparer` names the full transform it is capable of;
with BAQ off it performs the prefix. Same convention as step 1 — the type is named for its
capability, the config says what is active.

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

## 3. The output type — `PreparedRead` (reused from production)

The generic output is the **existing production `PreparedRead`**
([pileup/walker/mod.rs](../../../../src/pileup/walker/mod.rs)), **reused as-is** — it already
*is* the F3+BAQ output, field for field (§8), exactly as step 1 reuses `MappedRead`. It is
reproduced here for the reader, **not redefined**:

```rust
/// The generic prepared read — the generic ReadPreparer's `type Prepared`. A still-*decomposable*
/// read: the pileup walker turns it into per-position events. REUSED from production, not a new
/// ng type (§8).
pub struct PreparedRead {
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

What `PreparedRead` carries that a raw `MappedRead` does not: a **left-aligned** CIGAR,
**BAQ-capped** qualities, a precomputed alignment end, a precomputed `mq_log_err`, decoded strand
and mate role. It carries **no** per-base overlap adjustment — that is pairwise and happens in
the walker (§6), which is what keeps a `PreparedRead` self-contained and pairwise-independent.

The generic implementation:

```rust
pub struct LeftAlignBaqPreparer<Raw: RawRefSeq, Canon: RefSeq> {
    raw_ref: Raw,        // raw, case-preserving bytes — left-alignment (the aligner's own view)
    canon_ref: Canon,    // canonical uppercased bytes — the BAQ HMM (unused when BAQ is off)
    ref_buf: Vec<u8>,    // reused fetch scratch, as step 1's ReadFilter does
    /// `None` = BAQ off — the freebayes-style left-align-only mode (§2). Both modes are v1;
    /// `Option` marks "absent", never a sentinel (the step-1 config convention).
    baq: Option<BaqConfig>,
}

impl<Raw: RawRefSeq, Canon: RefSeq> ReadPreparer for LeftAlignBaqPreparer<Raw, Canon> {
    /// The generic path needs no locus: it prepares the read against the reference around the
    /// read's own span (production's `process_read` takes no window or locus either).
    type Locus = ();
    type Prepared = PreparedRead;
    fn prepare_read(&self, read: &MappedRead, _locus: &()) -> Option<PreparedRead> {
        /* fetch around read span → left-align → BAQ if configured */
    }
}
```

The preparer **holds** both accessors rather than receiving a materialised window — the two
reference views are exactly why a single passed-in span could not work (preamble §3).

---

## 4. Reference dependency and the "no observation" reasons

**Reference (reuse the `RefSeq` split, `ref_seq.md`).** Production reads the reference two ways
here, and ng reuses exactly that split: **raw, case-preserving bytes** (`RawRefSeq`) for
left-alignment — matching the aligner's own view — and **canonical uppercased bytes** (`RefSeq`)
for the BAQ HMM. No new accessor. With **BAQ off**, only the raw accessor is touched — the
canonical fetch happens solely for the HMM, so the left-align-only mode needs no canonical
reference at all.

**`None` reasons (tallied):** `Baq` — BAQ declined the read (HMM overflow, reference window past
the contig end, no aligned `M` op, an `N`/ref-skip in the CIGAR). With BAQ disabled by config,
the transform is left-align only and never returns `None`.

---

## 5. Cross-cutting concerns

- **Performance / parallelism.** Pairwise-independence (preamble §3) makes generic prep
  embarrassingly parallel with deterministic output. BAQ (the per-read HMM) is the cost; reuse a
  per-worker scratch buffer for its matrices rather than allocating per read (the project's
  scratch-buffer discipline). `PreparedRead` is owned per read; reusing its allocations is the
  same locus-stream-level question deferred in `read_filtering.md` §6, not a per-read decision.
- **Determinism.** No mate/column context here means the same read prepares to the same
  `PreparedRead` regardless of thread interleaving — a property downstream determinism relies on.

---

## 6. Deferred, with a recommended home

- **Adaptor-mask application, mate-overlap reconciliation, CIGAR decomposition → the pileup
  walker.** `PreparedRead` is *decomposable*, not decomposed; it *carries* the adaptor boundary
  but does not apply it (preamble §5).
- **Local haplotype reassembly → a future generic `ReadPreparer` sibling.** GATK-style de-novo
  assembly + realign-to-haplotype is the "reassemble" pole of step 2's axis (`ng_proposal.md`
  §2) — a genuine bench competitor to `LeftAlignBaqPreparer`, added as a sibling `impl ReadPreparer` when
  we test it, not part of v1.

- **freebayes' extra alignment passes around indels → a future refinement.** freebayes does more
  than a single left-align pass to keep indel placement stable. The mechanism we know of is
  `stablyLeftAlign`, which iterates `leftAlign` **to convergence** rather than aligning once; our
  port currently does one GATK-style `leftAlignIndels` pass. Whether that iteration is the whole
  of it — or freebayes does further local realignment around indel-bearing reads — needs pinning
  against the vendored source (`freebayes/src/LeftAlign.cpp`, `AlleleParser::registerAlignment`)
  **before** we implement, rather than guessing from the taxonomy. Then measure whether the extra
  passes actually move indel placement on our data. (freebayes is MIT-licensed, so reading and
  porting it is fine — unlike the AGPL TRF-mod case.)

**Alternatives to bench (recorded, not built now):** GATK-style **local reassembly** (above) — a
`ReadPreparer` sibling producing a `PreparedRead`-shaped output, so the pileup walker consumes it
unchanged. *(The freebayes-style trust + left-align-only preparation is no longer listed here: it
is a v1 config mode, §2 — not a deferred alternative.)*

---

## 7. Reuse over rewrite — the map to production

The parity oracle is the production prepared read, **in both modes**: a ported generic impl is
correct when its `PreparedRead` is byte-identical to the production `PreparedRead` on a fixture —
same left-aligned CIGAR, and same qualities under **BAQ on** (vs production's default) *and*
**BAQ off** (vs production's `--no-baq` passthrough). Two parity fixtures, one per mode; the
BAQ-off one also proves the left-align stage in isolation, since nothing else touches the
qualities.

| what | existing code | ng reuse |
|---|---|---|
| the per-read prep fold | `process_read` ([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs)) | model for `LeftAlignBaqPreparer` — **only** its F3 + BAQ stages (G2/F1 are step-1 filters) |
| indel left-alignment | `left_align_indels` / `left_align_cigar` ([pileup/walker/indel_norm.rs](../../../../src/pileup/walker/indel_norm.rs)) → `normalize_alleles` ([norm_seqs.rs](../../../../src/norm_seqs.rs)) | call directly (GATK `leftAlignIndels` port over a representation-neutral kernel) |
| BAQ (on) | `BaqEngine::process` ([pileup/per_sample/baq_engine.rs](../../../../src/pileup/per_sample/baq_engine.rs)), htslib `probaln_glocal` | call directly; `None` on BAQ-skip |
| BAQ (off) — the v1 left-align-only mode | `prepare_passthrough` (the `--no-baq` path, same file) | call directly — copies the raw `qual` into `bq_baq` uncapped; the precedent for BAQ-as-a-toggle (§2) |
| the prepared read | `PreparedRead` + `mapped_to_prepared` ([pileup/walker/mod.rs](../../../../src/pileup/walker/mod.rs)) | **reuse as-is** — it *is* the F3+BAQ output (§8); may want hoisting out of `pileup/walker/`, since step 2 produces it and the pileup only consumes it |
| reference | `RawRefSeq` (left-align) + `RefSeq` (BAQ) ([ref_seq.md](ref_seq.md)) | reuse as-is |

---

## 8. Resolved decisions

*Nothing on the generic path is open. Each entry records the choice, its evidence, and the
alternative it beat — kept here (rather than deleted) so the reasoning survives the decision.*

- **The prepared read — resolved: reuse production's `PreparedRead`, don't mint a parallel type.**
  The survey shows production's `PreparedRead` *is* precisely the F3+BAQ output the walker
  consumes, field for field — the same concept, not a different one. So step 2 reuses it as-is,
  exactly as step 1 reuses `MappedRead`. (An earlier draft called it `PreparedReadNg` to dodge a
  name clash with a type that turned out to be the same thing; the `Ng` suffix named *where the
  type lived*, not what it is, and would have become a lie on port-back. Rust namespaces by
  module anyway — `ng::read::PreparedRead` and `pileup::walker::PreparedRead` could coexist, so
  the clash never justified a suffix.) **Fork trigger:** if a later generic impl (reassembly)
  needs a richer payload, define an ng-owned `ng::read::PreparedRead` *then* — widening a
  production type to serve the lab is the thing to avoid. Until that bites, reuse. Reuse may
  want `PreparedRead` hoisted out of `pileup/walker/` to a shared home (§7).
- **The reference argument — resolved: there isn't one.** The preparer **holds** its `RawRefSeq`
  + `RefSeq` accessors and fetches around each read's own span, rather than receiving a
  materialised window (`RefWindow`/`LocusWindow`, both now dropped). Evidence: production's
  `process_read(read, baq, raw_ref, cfg)` takes no window or locus and holds the fetchers itself;
  the generic transform needs *two* reference views (raw for left-align, canonical for BAQ) that
  one span cannot carry; and step 1's `ReadFilter` already sets the precedent (`reference: R` +
  a reused `ref_buf`). The generic path's `type Locus = ()` — it needs no locus at all
  (preamble §3).
- **Where generic prep is invoked — resolved: by the pileup, per read, as it walks a non-STR
  stretch.** We follow production, where `process_read` runs as the walker ingests reads. This is
  what makes the compose-not-subsume resolution (preamble §2) concrete: the pileup **calls**
  `prepare_read` on each read it ingests and consumes the resulting `PreparedRead` — it does not
  absorb the preparation. The `pileup/` spec inherits this as a given, not a question to reopen.
  *Alternative considered:* the pileup subsumes the preparation into its walk (the
  `module_layout.md` "subsume or compose?" pole). Rejected — production already draws the seam at
  `PreparedRead`, and subsuming would dissolve step 2's bake-off surface (a reassembly preparer
  could no longer be swapped in behind the same contract).
