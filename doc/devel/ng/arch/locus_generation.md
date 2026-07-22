# ng — locus generation: types & interfaces (the shared shape)

*Architecture draft (2026-07-22), the code-facing companion to
[`../spec/locus_generation.md`](../spec/locus_generation.md) — the middle arrow of the caller
spine (reference → typed regions → **a sample's loci** → calls). It gives the types and
signatures as they appear in code; **every *why* points back to the spec**, which is the
authority and is not re-argued here. It **extends** the shared arch docs, it does not restate
them: the domain newtypes and `LocusKind` live in
[`ng_step_interfaces.md`](ng_step_interfaces.md) (§1, §2) and are referenced, not re-typed;
[`module_layout.md`](module_layout.md) places `locus_generation/` in the tree. The first generator has its
own doc, [`locus_generation_ssr.md`](locus_generation_ssr.md). Naming: **STR** in prose, `ssr`
in code. The `SampleReads` / `IngestError` citations come from the read-ingestion line, now merged
into this branch (via `main`, commit 3ec3f97), so they resolve in this worktree. Signatures are
illustrative; the **contract** is the deliverable.*

## Module home

`src/ng/locus_generation/` — a **folder**, per `module_layout.md` principle 1: a step whose implementations
compete keeps its trait and every implementation side by side, and ng expects several generators
per region kind (spec §5). Defer the tree-wide rules to the module-layout doc; this unit owns:

```
src/ng/locus_generation/
  mod.rs   – SampleLocusObservations, ObservedSequence, ReadCoverage, LocusGenerator,
             the dispatcher, NoLoci, the iterator + its config/counts/error
  ssr.rs   – the STR generator (locus_generation_ssr.md)
  pileup/  – the generic generator (deferred, spec §11)
```

## 1. The locus type

The unit every generator produces and every later step consumes. **New code shape** — it
supersedes and renames `ng_step_interfaces.md` §2's `LocusEvidence` sketch (which borrowed its
reads and carried a redundant `sample`; both wrong — spec §3). `LocusKind` / `SsrDetail` are the
one part *not* re-typed here: they are shared vocabulary in `ng_step_interfaces.md` §2, defined
concretely in spec §3, and `SsrDetail`'s fields get their code home where the STR generator mints
them (`locus_generation_ssr.md` §1).

```rust
/// One sample's locus: the stretch of genome it covers, and what that sample's reads
/// showed there. **Owned, no lifetimes** — a cohort stage merges these across samples and
/// an artifact writes them, so it must outlive every buffer it was built from (spec §3).
pub struct SampleLocusObservations {
    /// The stretch this locus covers — one base for a candidate SNP, the tract for a
    /// microsatellite.
    pub region: GenomeRegion,
    /// The REF bases over `region` — what a wider-span projection needs when samples merge.
    pub reference_bases: Box<[u8]>,
    /// The distinct sequences the reads showed, each with its support. Observations, not
    /// alleles.
    pub observed_sequences: Vec<ObservedSequence>,
    /// Reads that covered the locus but produced no observation. A scalar with no
    /// positions: "coverage that said nothing" ≠ "no coverage" (spec §3).
    pub reads_without_observation: u32,
    /// Reads a depth cap discarded. Non-zero ⇒ the support counts are a subsample, not the
    /// depth (spec §3).
    pub reads_discarded_by_cap: u32,
    /// What kind of locus this is, and the extras that kind carries. `LocusKind` +
    /// `SsrDetail`: ng_step_interfaces.md §2 / spec §3 (authoritative). Read the kind off
    /// the type; never infer it from which fields are populated.
    pub kind: LocusKind,
}

/// One distinct sequence the reads showed, with its support. The moments between `num_obs`
/// and `chain_ids` are the per-read statistics the SNP filters read; an STR model reads
/// only `num_obs` (spec §3). Ported field-for-field from production's per-allele shape
/// (reconciliation below).
pub struct ObservedSequence {
    /// The observed bases — allele content, in **read** coordinates.
    pub bases: Box<[u8]>,
    /// How much of the locus a read of this sequence spanned. **Part of the identity**: a
    /// `Complete` and a `PartialLeft` of the same bases are different evidence and stay
    /// separate entries (spec §3).
    pub read_coverage: ReadCoverage,
    /// How many reads showed this sequence — the whole support on the STR path.
    pub num_obs: u32,
    /// Forward-strand reads (strand bias).
    pub num_fwd: u32,
    /// Σ per-read log-error — the freebayes per-read error term.
    pub q_sum: f64,
    /// Σ MAPQ and Σ MAPQ² — the multi-mapper Welch's-t filter recovers mean and variance.
    pub mapq_sum: u32,
    pub mapq_sum_sq: u64,
    /// Phase-chain ids of the reads folded here — what lets a later step chain observations
    /// into a haplotype. Read-position-bias fields (`placed_left`/`placed_start`) are
    /// **not** here: anchor-relative, degenerate on a tract, so the pileup generator adds
    /// them (spec §3, §11).
    pub chain_ids: Vec<ChainId>,
}

/// How much of the locus a single read spanned — **one read's span, not depth**.
/// `Complete` = it reached both borders; a partial ran off its own end partway, carrying
/// how many of the locus's positions it reached from that border, in **locus** coordinates
/// (spec §3).
pub enum ReadCoverage {
    Complete,
    PartialLeft(u16),
    PartialRight(u16),
}

impl SampleLocusObservations {
    /// Read depth at each position of `region`, in order — **derived, not stored**. A
    /// `Complete` counts at every position, a `PartialLeft(n)` at the leftmost `n`. Length
    /// = `region.len()`. This is *observation* depth and only exact per locus; the paralog
    /// filter owns the covering-but-unobserved and overlapping-loci caveats (spec §3, §11).
    pub fn num_obs_along_locus(&self) -> Vec<u32>;

    /// The observations a likelihood may score directly — the `Complete` ones. A partial is
    /// a lower bound that mis-scores as a short allele until step 7 models censoring, so
    /// reaching the partials is a deliberate act — this method is the guard (spec §3).
    pub fn complete_observations(&self) -> impl Iterator<Item = &ObservedSequence> + '_;
}
```

## 2. The generator contract

A generator takes one typed segment plus a way to read a sample's reads and yields that sample's
loci there — zero, one, or many, one at a time (spec §4).

```rust
/// Generates a sample's loci from one segment of kind `S`, streaming **one locus at a
/// time**. `S` is a parameter on the *contract*, not an associated type inside each impl,
/// so two generators for the same kind stay swappable behind `Box<dyn LocusGenerator<S>>`
/// (§4 decision; spec §4).
pub trait LocusGenerator<S> {
    /// Start a new segment: reset progress. Gathers nothing; cannot fail.
    fn begin_segment(&mut self, region: GenomeRegion);

    /// The next locus, or `None` once the segment has no more. `&mut self` because a
    /// generator owns reusable scratch (alignment matrices, sampling buffers) that must not
    /// be reallocated per segment. `segment` is borrowed **per call** (`&S`) rather than
    /// stored, so the generator stays lifetime-free. `Ok(None)` immediately is a normal
    /// outcome, not a failure.
    fn next_locus(
        &mut self,
        segment: &S,
        reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError>;
}
```

**Contract.** Loci come out in coordinate order (spec §2), one resident at a time regardless of
how many a segment yields — the memory guarantee `next_locus` exists to hold (spec §4, §9). A
generator holds its own accessors (reference, aligner, scratch) as fields, not arguments (spec
§4). Reads are not capped at this seam; a generator that needs a per-locus depth bound owns that
decision and its determinism (spec §4).

## 3. The dispatcher, and the generator that makes nothing

The dispatcher is an exhaustive `match` over `RegionKind` (reconciliation below), handing each
branch's payload to the generator the run supplies for it. **Every branch resolves to a
generator** — no `_ => {}`, no `TODO`; an unfilled kind falls back to `NoLoci`, so this shape run
alone produces nothing (spec §5).

```rust
/// A generator that produces no loci and counts what it passed over. One impl covers every
/// kind because it ignores the segment (`impl<S> LocusGenerator<S> for NoLoci`).
pub struct NoLoci {
    pub reason: UnhandledReason,
}

/// **Why** a kind produces no loci — a boundary we chose vs. a gap we have not filled. The
/// two must not be added together: they answer different questions (spec §5, §7).
pub enum UnhandledReason {
    /// Deliberately outside the caller's scope — e.g. satellite arrays. Permanent.
    OutOfScope,
    /// No generator written yet. Temporary by construction.
    NotImplemented,
}
```

*Impl-time confirmation (not an open design item).* The generator set the dispatcher routes over — one slot per
region kind, each a `Box<dyn LocusGenerator<…>>` or `NoLoci` — has no shape the spec pins; a
`GeneratorSet` struct with a per-kind field is the natural form. Confirm when the second generator
lands, since only then does more than one slot hold real work.

## 4. The public surface

An iterator over a typed-region stream, mirroring `TypedRegionIterator` (reconciliation below) —
same lazy shape, same running counts, same in-stream fatal error.

```rust
/// Lazily turns a typed-region stream into a sample's loci. Holds **no buffer of loci**: it
/// calls `next_locus` until the generator returns `None`, then pulls the next region and
/// calls `begin_segment` — one locus resident at a time (spec §6, §9).
pub struct SampleLocusObservationsIterator<T, G> { /* … */ }

impl<T, G> Iterator for SampleLocusObservationsIterator<T, G>
where
    T: Iterator<Item = Result<TypedRegion, TypedRegionError>>,
{
    /// A fatal condition is yielded once as `Some(Err(_))` then `None`, so `?` makes it
    /// un-ignorable rather than a silent end of stream (spec §6; read_filtering.md §5). A
    /// read that yields no observation is a tallied per-read outcome, never an error.
    type Item = Result<SampleLocusObservations, LocusGenerationError>;
}

impl<T, G> SampleLocusObservationsIterator<T, G> {
    pub fn new(regions: T, reads: SampleReads, generators: G, config: LocusConfig) -> Self;
    /// The running tally — current at any point, final once the stream is exhausted.
    pub fn counts(&self) -> &LocusCounts;
}

/// Fatal, run-level. `#[non_exhaustive]`, mirroring `TypedRegionError`'s in-stream-fatal
/// shape. Every variant wraps an upstream failure — this step mints no error of its own.
#[non_exhaustive]
pub enum LocusGenerationError {
    /// The upstream typed-region walk failed.
    TypedRegion(TypedRegionError),
    /// A read query or the alignment input failed.
    Reads(IngestError),
    /// A reference fetch failed — a broken reference, or a region past a contig end.
    Reference(RefSeqError),
}
```

## 5. Config and counts

```rust
/// The dispatcher's own knobs — and only those. Each generator owns its knobs and takes
/// them at construction, so adding a generator never widens this struct (spec §7). Minimal
/// in v1; defaults as named `pub const`s, `Option` = absent, when a field arrives.
pub struct LocusConfig { /* … */ }

/// Running tally — "no silent caps": every region and every base is accounted for, so
/// "how much genome does this caller not cover, and how much of that is temporary?" is
/// answerable from the counts alone. The base counters are the other half of why `SsrBundle`
/// and `Satellite` exist as types rather than holes (spec §7). Per-generator counts (reads
/// fetched/dropped and why) belong to the generator — no meaning shared across kinds.
pub struct LocusCounts {
    pub regions_in: u64,
    pub loci_emitted: u64,
    pub unhandled_not_implemented: u64,
    pub unhandled_not_implemented_bp: u64,
    pub unhandled_out_of_scope: u64,
    pub unhandled_out_of_scope_bp: u64,
}
```

## 6. Decisions — decided (why in the spec)

- **One unified `SampleLocusObservations`; `LocusKind` tags the kind** — the evidence is uniform
  across kinds, so nothing else in the record names it; a downstream step matches on `kind`
  rather than guessing from populated fields (spec §3).
- **Owned, no lifetimes** — it is what a cohort stage merges and an artifact writes. `chain_ids`
  are kept at a known cost the coder must plan for: production's REF chain-ids were ~31% of peak
  live heap, so a generator filling this field pays that unless it skips REF (spec §3).
- **Stream one locus at a time, not a filled `Vec`** — a memory decision (a `Generic` segment is
  up to ~100k loci) that also matches how the walker closes loci as it passes; a filled `Vec`
  would be a buffer the interface imposes, not the algorithm (spec §4).
- **The segment type is a trait *parameter*, not an associated type** — the fork that decides
  whether one generator is swappable for another at all. *Rejected:* taking the whole
  `TypedRegion` (forces an impossible branch into every generator), and an associated input type (a
  trait object cannot name it, so "some generator chosen at run time" is inexpressible) (spec §4).
- **Segment borrowed per call; `Box<dyn Iterator>` rejected** — holding the borrow forces a
  lifetime that spreads to the dispatcher and iterator; a boxed iterator is opaque to the optimiser
  and allocates per segment for nothing (spec §4).
- **Every branch resolves to a generator; `NoLoci` carries a reason** — routing an unbuilt kind
  to a counted generator-with-a-reason beats a silent `_ => {}`, and keeps the two kinds of
  nothing (`OutOfScope` vs `NotImplemented`) separate so the manager can ask them apart (spec §5).
- **Only `Satellite` is resolved here** (`OutOfScope`, permanent); every real generator plugs in
  from its own spec — this shape is deliberately independent of the STR generator (spec §2, §5).
- **Item is a `Result`, error is fatal/run-level, in-stream** — same shape as `read_filtering.md`
  §5 and `TypedRegionIterator`; this step mints no error of its own, only wraps upstream ones
  (spec §6).
- **A generator sees exactly one region at a time — decided (owner, 2026-07-22).** Each genomic
  segment is resolved independently, so no generator needs cross-segment context; the gap-free
  tiling means a neighbour is generated next and the merge reconciles by overlap. Confirms §2's
  contract and closes the one assumption a later generator could have invalidated (was spec §12's
  open question).
- **Reference and `@SQ` match exactly, enforced at read ingestion — decided (owner, 2026-07-22).**
  Reads were mapped against the reference, so each file's `@SQ` list must equal the reference
  contig list — names, lengths, order — a check read ingestion runs at file open
  ([open_bam.rs:243](../../../../src/ng/read/input/open_bam.rs#L243)). A contig in the reference
  but absent from a sample's `@SQ` therefore **cannot reach this step**: there is no benign-mismatch
  case to absorb, and §4's fatal-error model is never at risk from one. Supersedes spec §8 item 4 /
  §12's "zero reads" lean.

## 7. Open items

- `OPEN:` **One locus type carrying N samples vs. a single-sample type the cohort composes** —
  leaning to compose (`CohortLocus` holds `Vec<Option<SampleLocusObservations>>`, production's
  shape). Settled by whether ng's generic path projects alleles onto a group span; confirm when
  the pileup generator is designed — does not block v1 (spec §12).
- *Impl-time confirmations, not design items:* the `GeneratorSet` shape (§3); `LocusConfig`'s
  active fields (none in v1 — generators own their knobs).

## 8. Reconciliation with existing code

Every row read at the cited line (2026-07-22). Convergence on existing types, not new types
alongside them.

| ng name | existing code | action |
|---|---|---|
| `RegionKind` / `TypedRegion` (dispatch input) | [region_typing/mod.rs:170](../../../../src/ng/region_typing/mod.rs#L170), [:146](../../../../src/ng/region_typing/mod.rs#L146) | consume as-is; the `match` is exhaustive over its four variants (`SsrSegment`/`SsrBundle`/`Generic`/`Satellite`) |
| iterator + counts shape | `TypedRegionIterator` [region_typing/mod.rs:789](../../../../src/ng/region_typing/mod.rs#L789), `TypedRegionCounts` [:269](../../../../src/ng/region_typing/mod.rs#L269) | mirror — same lazy `Item = Result<_,_>`, same `counts()` accessor |
| `GenomeRegion` / `Position` / `ContigId` | [types.rs:77](../../../../src/ng/types.rs#L77), [:32](../../../../src/ng/types.rs#L32), [:11](../../../../src/ng/types.rs#L11) | reuse as-is (1-based inclusive, `u64`); `region.len()` derives `num_obs_along_locus()`'s length |
| `ObservedSequence` | `AlleleObservation` [pileup_record.rs:138](../../../../src/pileup_record.rs#L138) + `AlleleSupportStats` [:44](../../../../src/pileup_record.rs#L44) | **model** — the moments (`num_obs`/`q_sum`/`fwd`/`mapq_sum`/`mapq_sum_sq`) map field-for-field; ng drops `placed_left`/`placed_start` to the pileup generator and adds `read_coverage` |
| `ChainId` | [pileup_record.rs:30](../../../../src/pileup_record.rs#L30) (`type ChainId = u64`) | reuse as-is; production's REF-chain-id drop (§8 there) is the memory lesson §6 records |
| `CohortLocus` composing `SampleLocusObservations` | `CohortPileupRecord.per_sample: Vec<Option<PileupRecord>>` [var_calling/types.rs:39](../../../../src/var_calling/types.rs#L39) | **model for the cohort type** — not this step; the merge groups by overlap (spec §3, §11) |
| `SampleReads` (read access) | [read/input/mod.rs:345](../../../../src/ng/read/input/mod.rs#L345), `reads_in_region` [:440](../../../../src/ng/read/input/mod.rs#L440) | reuse as-is; §8 there is this step's feedback (per-query `Arc` accessor) — read on `main`, not this worktree |
| `LocusGenerationError::Reads` | `IngestError` [read/input/mod.rs:584](../../../../src/ng/read/input/mod.rs#L584) | wrap — read on `main` |
| `LocusGenerationError::Reference` | `RefSeqError` [ref_seq.rs:39](../../../../src/ng/ref_seq.rs#L39) | wrap |
| `LocusGenerationError::TypedRegion` | `TypedRegionError` [region_typing/mod.rs:1640](../../../../src/ng/region_typing/mod.rs#L1640) | wrap |
| reference accessor a generator holds | `RefSeq`/`RawRefSeq`/`EvictableRefSeq` [ref_seq.rs:142](../../../../src/ng/ref_seq.rs#L142), [:180](../../../../src/ng/ref_seq.rs#L180), [:230](../../../../src/ng/ref_seq.rs#L230) | a generator's own field (spec §4) |
| `MappedRead` (generator's read input) | [bam/alignment_input.rs:78](../../../../src/bam/alignment_input.rs#L78) | reuse as-is |
| `LocusKind` / `SsrDetail` | [ng_step_interfaces.md §2](ng_step_interfaces.md), spec §3 | referenced, not re-typed — shared vocabulary |

## Test & bench shape

The regression anchor is the dispatch-and-accounting test (spec §13), provable independent of any
generator: the `match` is exhaustive (compiler-proven), `regions_in` reconciles against loci +
both unhandled counters, the two kinds of nothing stay in separate counters, order is preserved,
and `num_obs_along_locus()` derives correctly from mixed complete/partial observations. **No
`bench/`** at this seam — the dispatcher has no bake-off; the generators do. Tests live beside the
code in `locus_generation/mod.rs`; the STR generator's own test is in its spec (`locus_generation_ssr.md` §9).
