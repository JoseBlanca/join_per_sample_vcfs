# SSR shared types — `src/ssr/types.rs`

**Status:** first draft, 2026-06-11 — a discussion starter for the SSR
caller's **shared domain model**, the first per-module pass after the overall
[architecture doc](ssr_genotyping_architecture.md) (its §8 module 0).

`types.rs` is the **spine every stage shares**: the catalog defines loci,
`ssr-pileup` emits per-read likelihoods over candidate alleles, `ssr-call`
unions them across the cohort, the VCF writes them back. Get these types right
and each stage pass is clean; get the **allele representation** wrong and the
fix touches all four stages — which is why it is settled first.

This document designs the *Rust types*. The spec
([`ssr_genotyping.md`](../specs/ssr_genotyping.md) §1, §4.2, §4.3, §5.1)
settles the *concept*; here we make it concrete. **Proposals are flagged to
argue with**, not conclusions.

---

## 1. The invariants these types must encode (spec recap)

Four load-bearing rules from the spec, restated as type obligations:

1. **An allele's identity is its tract *sequence*, never a bare count** (§1).
   Two molecules differing by a non-motif amount (12 clean repeats vs 12 + 1 bp)
   are **distinct alleles**. The "count" this rules out is the **lossy
   projection** `sequence → floor(len/period)`, which maps *both* of those onto
   `12` and so fuses two distinct alleles — that projection is never a key. It
   is **not** the on-ladder `units` field of §2, which is a *lossless* encoding
   of the sequence; the *identity vs. encoding* note there draws the line.
2. **Two encodings of that sequence** (§4.2): **on-ladder** = a clean rung
   `ref ± k` units, sequence reconstructible from `(reference tract, k)` — so it is the
   *same allele in every sample by construction*; **off-ladder** = a non-rung
   sequence carried **explicitly** as its canonical normalized sequence (§4).
3. **Cross-sample identity is by sequence-key** (§4.2, §5.1). The union that
   builds the candidate set `A_ℓ` requires: the *same* on-ladder rung, or the
   *same* normalized off-ladder sequence, in two samples → *one* allele. This is the
   single property that makes the whole cohort model work.
4. **Alleles are ordered by length** (§5.4) — no label-switching in the EM.

> **"on-/off-ladder" vs "pure/impure" — two orthogonal axes, kept distinct on
> purpose.** Both sound like "clean vs messy," so they are easy to fuse — don't.
> The **ladder** is the integer ladder of clean rungs `… ref−1, ref, ref+1 …`
> (the field's *allelic / stutter ladder*); the two axes are:
>
> - **pure / impure** (a.k.a. purity, perfect/imperfect) — is the *reference
>   tract of a locus* a clean motif tiling, or interrupted? A property of the
>   **catalog locus**, fixed once at build time (spec glossary, §3, §5.2).
> - **on-ladder / off-ladder** — is *an allele's* length a whole-motif **rung**
>   of that locus's reference tract (`ref ± k` units), or off by a non-motif
>   amount? A property of an **allele**, per sample (spec §4.2). On-ladder ⇒
>   sequence reconstructible from `(reference tract, k)`, so the same allele in
>   every sample by construction; off-ladder ⇒ carried explicitly.
>
> They are **independent**. The case needing both at once: an **on-ladder allele
> at an impure locus** — interrupted reference `(CAG)₅ CAA (CAG)₄` (impure
> *locus*) with allele `(CAG)₅ CAA (CAG)₆`, two clean units added on one side
> (on-ladder — a clean rung of the interrupted reference tract, the same `CAA`
> pinned in, spec §3.1). Saying that is impossible if both axes share a word.
>
> **Division of labour:** `types.rs` models only the **on-/off-ladder** axis
> (the `Allele` enum, §2); **pure/impure** lives on the catalog side (`Locus`'s
> `purity_fraction`, §5) and is consumed by the Stage-2 stutter model, not by
> `Allele`.

A reassuring anchor: the **SNP path already lives by rule 1** —
[`AlleleObservation.seq: Vec<u8>`](../../src/pileup_record.rs) *is* the literal
allele bytes, identity-by-sequence. SSR is the same principle plus a compact
rung encoding for the on-ladder common case.

---

## 2. The central type — `Allele`

The crux. Proposal: an enum that is **two encodings of one sequence**, not two
kinds of thing.

```rust
/// A repeat allele *at a given locus*. Identity is the tract sequence; the
/// variants are two encodings of it. Locus-relative (see §3): an `Allele`
/// is only meaningful — and only comparable — within one locus's candidate set.
#[derive(Clone, PartialEq, Eq, Hash)]
enum Allele {
    /// Clean rung: sequence = the locus's reference tract tiled to `units` motif
    /// copies (carrying any fixed interruption structure of an imperfect locus,
    /// §3). Reconstructible from (reference tract, units); identical across
    /// samples by construction. `units` is the absolute repeat count, not the
    /// signed Δ.
    OnLadder { units: u16 },

    /// Non-rung sequence, carried explicitly as its canonical normalized
    /// sequence (§4) — same bytes in two samples ⇒ same allele (`==`).
    OffLadder(NormalizedSeq),
}
```

> **Identity vs. encoding — the one subtlety to internalize.** Rule 1 says
> identity is the *sequence*, never a count — yet `OnLadder` stores `units`, a
> count. These do **not** conflict, because `units` is a *lossless encoding* of
> the sequence, not the *lossy projection* rule 1 forbids:
>
> - The forbidden count is `sequence → floor(len/period)` — **many-to-one**. It
>   maps on-ladder "12 clean" *and* off-ladder "12 + 1 bp" both onto `12`.
> - `units` is a **rung coordinate** valid only on the on-ladder set, where —
>   with the locus fixed — `units ↔ sequence` is a **bijection**. Given the
>   locus, `to_sequence` recovers the exact bytes; nothing is lost.
>
> The collapse rule 1 warns about never happens **because off-ladder alleles are
> kept out of the integer ladder entirely** — "12 + 1 bp" is
> `OffLadder(NormalizedSeq)`, not `OnLadder { 12 }`, so it can never fuse with
> the clean rung:
>
> | molecule | `Allele` value |
> |---|---|
> | 12 clean repeats | `OnLadder { units: 12 }` |
> | 12 repeats + 1 bp | `OffLadder(NormalizedSeq(…))` — a different variant, never `==` |
>
> So the enum's `==` **is** sequence-identity (the two variants partition
> sequence space disjointly; `OnLadder` is injective in `units`), just *computed*
> via a `u16` compare in the common case instead of comparing tract bytes.
> `units` is a fast path for the identity check, not a rival notion of identity.

Why this shape:

- **`OnLadder` is just an integer** — the compact, sample-independent rung. Per
  the note above, its `Eq`/`Hash` *are* sequence-identity for the on-ladder
  majority (same rung ⇔ same sequence, locus fixed) at the cost of a `u16`
  compare. The optimization over the SNP path's "store every allele as bytes":
  the on-ladder bulk — ~99% of alleles — never materializes a tract sequence
  until asked, and reconstructs losslessly when it is.
- **`OffLadder` carries its `NormalizedSeq`**, whose `Eq`/`Hash` carry rule 3
  for the minority. on-ladder vs off-ladder are never equal — an off-ladder
  allele is non-rung *by construction* (§5.2), so the variants can't collide.
- **`#[derive(Eq, Hash)]` gives the cohort union for free** — `A_ℓ` assembly
  (§5.1) becomes a `HashSet<Allele>` insert, *provided* `NormalizedSeq` holds a
  canonical form (§4). The derive is correct only if those bytes are canonical;
  that is the one thing §4 must guarantee.

> **Decided — `units` is `u16` end-to-end.** A repeat count is non-negative, so
> `u16` is the honest type *and* the storage type: the `.ssr.psp` on-ladder
> length columns are **`uint16`** (spec §4.3, updated to match) — no signed
> storage, no signed↔unsigned hop at the boundary, and unsigned counts encode as
> plain varint (tighter than the zig-zag a signed type would use). The only
> signed length quantity, **Δ = (units − ref_units)**, is a *derived* coordinate
> (§5), computed when the prior/stutter need it, never stored as identity.

---

## 3. `Allele` is locus-relative — the key scoping decision

**Proposal: `Allele` does NOT carry its locus.** `units` is meaningful only
against *this* locus's reference tract; an off-ladder `NormalizedSeq` is
canonicalized relative to *this* locus's ref tract. Rung-5 at one locus and
rung-5 at another are different sequences — but they are never compared, because
**the candidate set `A_ℓ` is per-locus** (§5.1), so the locus is *ambient
context*, supplied by the container, not stamped on every allele.

Consequences (all desirable):

- `Allele` stays tiny (a `u16` or a short `NormalizedSeq`) — it is stored millions of times
  across a cohort's `.ssr.psp` files; carrying a locus would bloat it.
- `Eq`/`Hash` are *within-locus* by contract. A `HashSet<Allele>` is always a
  per-locus set; the type system doesn't enforce the "same locus" precondition,
  so it lives in a doc-comment invariant (and debug asserts), exactly as the
  candidate-set type (§6) is the thing that owns a locus.
- Sequence reconstruction and repeat-count derivation are **methods that take
  the `Locus`** (§5), not stored fields.

> **Alternative considered:** make `Allele` self-describing (carry a
> `LocusId` or the whole `Locus`). Rejected: it bloats the hot type and
> duplicates context the candidate set already holds — doubly so now that
> `Locus` carries its reference bytes (§5). Flag if a cross-locus structure ever
> genuinely needs globally-comparable alleles.

---

## 4. `NormalizedSeq` — the canonical normalized sequence (DECIDED)

The hard part of rule 3. An off-ladder allele is a non-rung sequence (an
in-frame count + a 1 bp indel, a partial unit, a variable interruption). Two
samples that carry the *same* such allele must produce the *same* bytes, or the
union splits one biological allele into N.

**Decision (2026-06-11): the key is the canonical normalized *sequence*
itself.**

```rust
/// An off-ladder allele's identity: its tract written out, in canonical
/// left-aligned normalized form. Byte-equal ⇔ same allele (rule 3). The
/// canonicalization is the whole contract — two samples' encodings of one
/// biological allele MUST normalize to identical bytes.
struct NormalizedSeq(Box<[u8]>);   // interning is an optimization, not identity
```

Why the sequence, not a compact `(offset, ref, alt)` delta:

- It is the **most obviously-correct equality key** — a canonical byte string
  with no "did both samples normalize the tuple the same way" subtlety to get
  wrong; equality is just `[u8] == [u8]`.
- It **matches what's already true elsewhere**: the SNP path's
  identity-is-its-bytes (`AlleleObservation.seq`), and the `.ssr.psp`
  `offl_seqs` **dict** column (spec §4.3) stores these sequences directly — so
  serialization is the identity map, and the dict column de-dups the repeats
  the `Box<[u8]>` would otherwise cost.
- The size downside is moot: off-ladder columns are **empty at the vast
  majority of loci** (§4.3), so the cost is paid only where real off-ladder
  signal exists.

(A compact delta tuple was the alternative; rejected — it trades the bulletproof
byte-equality key for marginal memory we don't need, since the off-ladder set is
sparse and dict-compressed anyway.)

**Canonicalization — share the SNP normalizer, do not re-implement it
(DECIDED).** The spec mandates "**canonical left-aligned form, reusing the SNP
caller's indel-normalization discipline**" (§4.2), and the right reading of
"reuse" is *the actual code*, not a parallel copy of the rule — a hand-rolled
SSR normalizer is a divergence-bug farm (a fix in one place that never reaches
the other), and the cross-sample-equality invariant (§4) is exactly where that
bites. Grounding (from a read of
[indel_norm.rs](../../src/pileup/walker/indel_norm.rs)):

- The left-alignment **core is already a representation-neutral kernel** —
  `normalize_alleles(seqs: &[&[u8]], bounds, max_shift, trim)` — operating on
  byte sequences, *not* CIGARs (the CIGAR functions are a wrapper on top). So
  the primitive the off-ladder path needs already exists; it is just private.
- It is a **port of GATK `leftAlignIndels` + `normalizeAlleles`, cross-checked
  against freebayes** — the `bcftools norm` / `vt normalize` canonical form GIAB
  truth uses. Reference-validated, subtle code: precisely what must not be
  re-derived.
- It solves the **identical problem**: the module's own doc comment describes
  one biological indel landing at different offsets and "fragmenting support"
  across the cohort merge — swap "reads" → "samples" and that is §4's off-ladder
  union. One canonicalizer ⇒ SNP and SSR keys are produced by identical logic,
  and fixes propagate to both.

**Build-time shape** (settled when building `pileup`): lift `normalize_alleles`
(+ its `Range`/helpers) into a **shared, public, representation-neutral module**,
relocated out of `pileup/walker/` now that there are **two real users**; the SNP
CIGAR path keeps its wrapper and calls the kernel (regression gate: SNP e2e
tests), and the SSR off-ladder path is a **thin adapter** that builds
`(seqs, bounds)` from `(off-ladder candidate, ref tract)` and calls the same
kernel — not a CIGAR-faking shim, not a copy. `types.rs` only owns the
**invariant** that `NormalizedSeq` bytes are already canonical, so
`Eq`/`Hash`/serialization treat them as opaque.

---

## 5. Supporting types

```rust
/// A repeat unit, ≤ 6 bp (SSR scope). `period = bytes.len()`. Stored
/// **verbatim** — the reference-strand, phase-faithful unit as it appears at the
/// locus, *not* canonicalized: canonicalizing (rotating) would break tiling, and
/// reconstruction reads phase-correct bytes from `Locus::ref_bytes` anyway. The
/// canonical *class* (lexicographically-min rotation over the motif and its
/// reverse complement) is **derived on demand** as the stutter-model pooling key
/// (spec §5.2 "motif-or-GC class"), never stored.
struct Motif { bytes: ArrayVec<u8, 6> }   // or SmallVec; period ≤ 6 is the cap

/// One catalog locus — and it carries its own reference bases, so reconstruction
/// (and every stage past the catalog builder) is FASTA-free. (chrom, start, end)
/// are the ref tract [start, end); motif gives the period.
struct Locus {
    chrom:  ContigId,     // dict id, as in the .psp contig table
    start:  u32,          // 0-based, half-open — tract start
    end:    u32,          //                      tract end
    motif:  Motif,
    /// Fraction of the tract that is a clean motif tiling, in [0.0, 1.0]
    /// (1.0 = perfect, < 1.0 = interrupted) — a *degree*, not a flag (hence the
    /// `_fraction` suffix). A continuous covariate of the stutter model (§5.2)
    /// and a build-time filter threshold (§3.2). The perfect/imperfect *bool*
    /// is derived, never stored: `is_perfect() == (purity_fraction == 1.0)`.
    purity_fraction: f32,
    /// Embedded reference bases: the tract plus a bounded flank margin each side
    /// (clamped at contig ends). `ref_bytes_start` is the genomic coordinate of
    /// `ref_bytes[0]`, so the tract is `ref_bytes[(start - ref_bytes_start)
    /// .. (end - ref_bytes_start)]` and the flanks are what's left either side.
    /// Carrying these is what makes Stages 1–2 never open the reference FASTA;
    /// the flank margin is a catalog-build parameter sized to Stage 1's
    /// anchoring + pair-HMM band need (decided when building `pileup`).
    ref_bytes:       Box<[u8]>,
    ref_bytes_start: u32,
}
```

- **Reconstruction is a pure `Locus` method — no FASTA, no `Scaffold` type.**
  Because the locus carries its own bases, `to_sequence` is a pure function of
  `(allele, locus)`. For `OnLadder { units }` it tiles the motif `units` times,
  reading any fixed interruption structure straight off `ref_bytes` (so the
  imperfect-locus case that drove the old `Scaffold` idea is handled with no
  external lookup); for `OffLadder` it returns the `NormalizedSeq` bytes. The
  per-locus arm decomposition an impure tract needs is computed on demand from
  `ref_bytes` + `motif` — cheap, since a tract is tens of bp and the `2W+1`
  candidate haplotypes are built once per locus, not per read.

  ```rust
  impl Allele {
      /// This allele's tract sequence (the VCF REF/ALT bytes). Pure; reads
      /// `locus.ref_bytes`, never the reference FASTA.
      fn to_sequence(&self, locus: &Locus) -> Vec<u8>;
      /// Derived repeat count — integer for on-ladder, **fractional** for
      /// imperfect/off-ladder (§3.2/§5.9). Feeds VCF `REPCN` (rounded)/`BPDIFFS`.
      fn repeat_count(&self, locus: &Locus) -> f64;
  }
  ```

  This also settles the old "does `types.rs` depend on the reference reader?"
  fork: **it does not.** `Locus` is pure data that happens to *contain* bytes;
  the only thing that ever reads the FASTA is the Stage-0 catalog builder, which
  bakes `ref_bytes` in. (Trade-offs — the catalog grows by a sequence column and
  this is a deliberate denormalization of §3.2's minimal schema — live on the
  catalog/spec side, not here.)

---

## 6. `CandidateSet` (`A_ℓ`) — the per-locus allele universe

```rust
/// The candidate alleles at one locus the EM works over. Owns the locus; its
/// alleles are length-sorted (rule 4 — no label-switching). Built by the
/// cohort union (§5.1); the cap `MAX_CANDIDATE_ALLELES` is enforced here.
struct CandidateSet {
    locus:   Locus,
    alleles: Vec<Allele>,   // sorted; index is the EM's allele coordinate
}
```

This is the type that makes §3's "locus is ambient" concrete: a `CandidateSet`
*is* the locus context, and the EM addresses alleles by **index into
`alleles`** (a dense `usize` coordinate over a small set), while *identity*
across samples stays the `Allele` key. The two-level addressing — stable
`Allele` key for cross-sample union, dense index for the per-locus math — is
the same split the SNP posterior engine uses over its allele set.

---

## 7. Mapping to storage & VCF (consistency checks, not new design)

- **`.ssr.psp` (§4.3 / arch §10):** `OnLadder.units` → the `uint16`
  `hist_lengths` / `amb_lengths` columns; `NormalizedSeq` → the `offl_seqs` dict
  column (the normalized-sequence form makes this the identity map). The
  type ↔ column mapping is the `pileup` module's `locus_record.rs` job —
  `types.rs` just has to *serialize cleanly* to those shapes.
- **VCF (§5.9):** `allele.to_sequence(&locus)` → REF/ALT bytes (on- and
  off-ladder round-trip losslessly); `repeat_count` → `REPCN`/`BPDIFFS`. The
  VCF is sequence-first, so the reconstruction methods (§5) are what it needs —
  and because `Locus` carries its bases, even the VCF writer is FASTA-free.

---

## 8. Decisions (all settled)

The type model is closed. Settled:

- **`Allele` enum** — `OnLadder { units }` | `OffLadder(NormalizedSeq)`, two
  encodings of one sequence; the enum's `==` *is* sequence-identity (§2).
- **Locus-relative scoping** — `Allele` doesn't carry its locus; the candidate
  set is ambient context (§3).
- **`NormalizedSeq`** = the canonical normalized off-ladder sequence (§4).
- **`units` = `u16`** end-to-end (`.ssr.psp` columns are `uint16`; §2).
- **Reconstruction** — `Locus` carries its own reference bytes, so
  `to_sequence`/`repeat_count` are pure `Locus` methods; no `Scaffold` type, no
  FASTA past Stage 0 (§5).
- **Indel-norm** — share the SNP `normalize_alleles` kernel (lift it to a
  neutral module; SSR off-ladder is a thin adapter), do not re-implement (§4).
- **Motif** — stored **verbatim** (reference-strand, phase-faithful); the
  canonical class is **derived on demand** as the stutter pooling key, never
  stored (§5).

The two items that remain are not type questions but build-time mechanics,
already pointed where they belong: the exact `normalize_alleles` extraction (the
container/`pileup` refactor, §4) and `flank_bp` sizing (the catalog/`pileup`
build, spec §3.2).
