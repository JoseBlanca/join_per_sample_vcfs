# ng step 2 — read preparation (the shared discipline)

*Status: design spec, 2026-07-14. The **shared preamble** for step 2, read preparation. Read
preparation has two implementations that share a contract but produce different evidence to
different consumers, so each path has its own spec and this doc holds only what they have in
common. The two path specs:*

- *generic (SNP/indel): [`read_preparation_generic.md`](read_preparation_generic.md) —
  left-align + BAQ.*
- *STR: [`read_preparation_ssr.md`](read_preparation_ssr.md) — tract extraction.*

*Under [`ng_proposal.md`](ng_proposal.md) (§2) and the arch docs
[`../arch/ng_step_interfaces.md`](../arch/ng_step_interfaces.md) (the `ReadPreparer` trait sketch)
and [`../arch/module_layout.md`](../arch/module_layout.md) (the `read/` module). Follows
[`read_filtering.md`](read_filtering.md) (step 1). **No code yet.** Naming: **STR** in prose,
`ssr` in code.*

---

## 1. What read preparation *is* — goals, non-goals, and what it is not

Read preparation is the **per-read transform** that turns a filtered `MappedRead` into a
locus-ready observation. Read filtering (step 1) decided *which* reads survive; read
preparation decides *what each surviving read looks like as evidence at its locus* — realigning
or delimiting it against the reference and reconciling its per-base qualities, so the downstream
evidence-gatherer sees a canonical, comparable read.

It is **marker-aware**, and that is the one structural difference from step 1: the router (step
3) classifies each locus as generic (SNP/indel) or STR, and read preparation runs the matching
implementation. The two implementations share the discipline below but almost nothing else —
different algorithms, different output types, different consumers — which is why they are two
specs. This doc is the contract they both obey.

### Goals

- Turn each filtered read into a **per-read, pairwise-independent** prepared observation that the
  downstream gatherer can turn into locus evidence.
- **Reuse the production preparations** behind a shared trait (§3), one implementation per path.
- Keep read preparation **separate from the gatherer it feeds**, so the per-step bake-off
  (`ng_proposal.md` §2) can swap a preparation implementation and hold the rest.

### Non-goals (deliberately excluded — could be goals, are not)

- **The per-base, pairwise evidence-gathering** — applying the adaptor mask, reconciling
  overlapping mate qualities, decomposing the CIGAR into per-position events. These need the
  locus-column context and are **not per-read**; they live downstream in the gatherer (§5).
- **Read likelihood (step 7).** Preparation produces observations, never `Lr = P(read | allele)`.
  On the STR path in particular, prep emits the tract *bytes* and scoring is a separate step and
  a separate pair-HMM; fusing them would foreclose the read-likelihood bake-off.
- **Candidate generation (step 6) and the router (step 3).** Preparation is handed an
  already-routed locus; it does not decide STR-ness or enumerate alleles.
- **Local reassembly in v1.** GATK-style haplotype reassembly is a *future* generic-path
  alternative to bench, not the first implementation.
- **Base-quality recalibration (BQSR).** Production does not do it here; neither does ng.

**Read preparation is not filtering.** It never *drops* a read for a whole-read property — that
was step 1. It *may* return "no usable observation" for a read the transform itself defeats, and
those are tallied, not silent (§5).

---

## 2. Where read preparation sits — compose, not subsume

The most important shared architectural fact is one the production code already settled on
**both** paths: **read preparation is a distinct step, composed with the gatherer, not fused
into it.** In production the generic per-read fold (`process_read`) produces a `PreparedRead`
that the **pileup walker consumes**; the STR Stage-1 pipeline produces an extracted-tract
observation (`ReadObs::Sequence`) that the **tract tally + cohort likelihood consume**. Neither
gatherer *is* the preparation — each takes a prepared read as input.

This resolves the standing open question in `module_layout.md` (*"does `pileup/` subsume the
generic path's step 2, or is it built from it?"*): **built from it.** `ReadPreparer` is a real step;
the pileup (generic) and the tract tally (STR) are its consumers, one step upstream of where the
two paths converge at `LocusEvidence`.

```
generic:  MappedRead ─prepare▶ PreparedRead ─▶ pileup walker ─▶ LocusEvidence
STR:      MappedRead ─prepare▶ SsrTractObs    ─▶ tract tally    ─▶ LocusEvidence ─▶ Lr (step 7)
                                                                  └── the paths converge here ──┘
```

Because the paths only converge at `LocusEvidence`, the prepared read is **inherently
pre-convergence and path-specific** — which is exactly why forcing a single output type here
(an earlier draft's `LocusRead` enum) was wrong, and why the trait keeps the output path-owned
(§3). (Contrast step 6, candidate generation, where generic and STR legitimately share one
output type, `AlleleCandidates`, because that step sits *after* convergence.)

---

## 3. The shared contract — the `ReadPreparer` trait

One trait states the shared shape; **each path owns its output type** via an associated type. The
generic and STR outputs go to different consumers and are never held interchangeably (the router
picks the path before prep runs), so the output is never unified — only the contract is.

The names read as one idea: a **`ReadPreparer`** does **`prepare_read`** and yields a
**`Prepared`** observation. Note what that does *not* promise: `PreparedRead` is the **generic
path's** output specifically, not the universal step-2 output. The trait yields `Self::Prepared`
— `PreparedRead` on the generic path, `SsrTractObs` on the STR path — so the tidy
`ReadPreparer → PreparedRead` pairing holds for the generic path only. That asymmetry is
deliberate (§2): do not read it as an invitation to re-unify the two outputs.

```rust
pub trait ReadPreparer {
    /// What this implementation needs to know about the locus — path-owned:
    /// `()` on the generic path (it needs no locus at all — it prepares the read against the
    /// reference around the read's *own* span) or `SsrLocus` on the STR path (motif, borders,
    /// flanks — so the delimiter's gaps can be tract-aware).
    type Locus;

    /// The per-read prepared observation this implementation emits — path-owned:
    /// `PreparedRead` (generic) or `SsrTractObs` (STR). The two converge only downstream, at
    /// `LocusEvidence`, never here.
    type Prepared;

    /// Prepare one filtered read. The implementation **holds its own reference accessors**
    /// (`RefSeq` / `RawRefSeq`) as fields and fetches what it needs — there is no window
    /// argument (see below). `None` = the read produced no usable observation here (§5) — a
    /// tallied per-read outcome, not a run error.
    fn prepare_read(&self, read: &MappedRead, locus: &Self::Locus) -> Option<Self::Prepared>;
}
```

### No window argument — the preparer holds its accessors (a decision, and the alternative)

An earlier sketch passed a pre-materialised reference window (`window: &RefWindow`, then a
bundled `LocusWindow`). Production shows that is wrong on **both** paths:

- **Generic:** `process_read(read, baq, raw_ref, cfg)` takes **no window**. It holds the
  reference fetchers and reads around the read's own span. It also needs *two* views of the
  reference — raw, case-preserving bytes for left-alignment and canonical bytes for the BAQ HMM
  — which a single materialised span cannot carry.
- **STR:** `delimit_read(region_seq, region_qual, locus, model, scratch)` takes the **locus**,
  not a window — the tract structure is the whole point.

So the reference is a **field on the preparer**, fetched per read — exactly as step 1's
`ReadFilter` holds `reference: R` plus a reused `ref_buf` scratch — and the only per-call context
is the routed locus, which the generic path does not need at all (`type Locus = ()`).

*Alternative considered:* keep a window type carrying bases plus tract structure. **Rejected** —
it needs a name for a concept that isn't one ("a window of *what*?"); it duplicates production's
existing `RefSpan` (the sequence-carrying span); it would hand the generic path STR fields it
ignores — the same papering-over we rejected for the output (§2); and it still could not carry
the raw+canonical duality the generic transform requires.

Three properties every implementation upholds:

- **Per-read and pairwise-independent.** Every read is prepared in isolation — no mate, no other
  read, no locus-column context. Production relies on this to prepare reads in parallel with
  deterministic output, and it is why the genuinely pairwise work (mate-overlap reconciliation)
  is downstream, not here (§5).
- **Content-canonicalising, not content-inventing.** Preparation re-*places* what the read
  already says (left-aligning an indel that is already there, capping a quality already measured,
  extracting a tract already sequenced). It does not assemble new sequence. (Local reassembly,
  which *does* invent haplotypes, is the one deliberately deferred alternative — §5.)
- **`None` means "this read, unusable here"** — a normal, tallied per-read result, never a hidden
  drop and never a masked run error (§5).

The `read/` module holds every step-2 implementation side by side (`module_layout.md` principle
1): the generic impls (`left_align_baq.rs`, later `reassembly.rs`) and the STR impl
(`pair_hmm.rs`) are siblings, each an `impl ReadPreparer`. STR-ness is a property of the
*implementation*, not a separate subtree (principle 2).

---

## 4. Shared vocabulary and reuse

Common to both paths (path-specific reuse lives in each path spec's reuse map):

| what | existing code | ng reuse |
|---|---|---|
| the read itself | `MappedRead` ([bam/alignment_input.rs](../../../../src/bam/alignment_input.rs)) | reuse as-is — the step-2 input |
| reference access | `RefSeq` + `RawRefSeq` ([ref_seq.md](ref_seq.md)) | reuse as-is — held **as a field on the preparer** and fetched per read (§3), raw bytes where the aligner's view matters, canonical bytes for HMM alignment |

There is **no window/span type** to define: the preparer fetches the reference itself (§3), so
nothing needs to bundle bases with a locus. The only per-call context is the routed locus
(`Self::Locus`), whose STR shape is co-owned with the router spec (`LocusKind`) and is `()` on the
generic path. A per-sample running tally (the step-2 analogue of `ReadFilterCounts`) records every
`None` by reason — the "no silent caps" discipline — with the reasons enumerated per path (§5).

---

## 5. Shared error model, and deferred-with-a-home

**Error model (identical on both paths).** Two distinct outcomes:

- **A read produces no usable observation** — `prepare_read` returns `None`, the reason is
  tallied. Normal per-read result, not a run error. (The reasons differ by path and are listed
  in each path spec.)
- **A reference fetch fails** — a contig mismatch, a window past the contig end. A run-level
  configuration error, **fatal**, surfaced up front by window validation, never swallowed into a
  per-read `None` (as in `read_filtering.md` §5). So a `None` always means "this read, unusable"
  and never hides a broken reference.

**Deferred, with a home (shared):**

- **Adaptor-mask *application* and overlapping-mate quality reconciliation → the gatherer** (the
  pileup walker on the generic path; subsumed by tract extraction on the STR path). Per-base and,
  for mate overlap, pairwise — needing the locus-column context. Consistent with
  `read_filtering.md` §6, which deferred the same two operations. (The generic prepared read
  *carries* the adaptor boundary as an annotation but does not apply it — the walker does.)
- **Read likelihood `P(read | allele)` → step 7.** Preparation emits observations, not scores.
- **Local haplotype reassembly → a future generic `ReadPreparer` alternative** (see the generic spec).

Everything else — the transforms, the output types, the consumers, the path-specific reuse and
open questions — is in the two path specs.
