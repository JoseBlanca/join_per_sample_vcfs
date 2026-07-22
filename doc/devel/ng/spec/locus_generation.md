# ng — locus generation: typed regions → a sample's loci (the shared shape)

*Status: design spec, 2026-07-21. **No code yet.** The middle arrow of the caller spine
(`ng_proposal.md`): reference → typed regions → **a sample's loci** → calls. This is the shared
preamble; each segment kind's generator gets its own spec, the first being
[locus_generation_ssr.md](locus_generation_ssr.md). Stands on the typed-region generator
([typed_regions.md](typed_regions.md)) and on read ingestion ([sample_reads.md](sample_reads.md),
[alignment_file.md](alignment_file.md) — built on the `ng-read-ingestion` line; not in this worktree's
tree, so its citations were verified against that worktree).
Naming: **STR** in prose, `ssr` in code. Code-facing companion:
[../arch/locus_generation.md](../arch/locus_generation.md) (the types and interfaces).*

---

## 1. What locus generation *is* — goals, non-goals, and what it is not

The typed-region generator says *what the reference sequence is* at every position. Read ingestion
says *what reads a sample has* over any interval. Neither knows about the other. This step joins
them: **it turns each typed region into loci, and a locus is a stretch of the genome with one
sample's evidence attached.**

That stretch is often a single base — a candidate SNP is a locus one base wide. It can be a hundred,
when the locus is a microsatellite tract. And the evidence is **the same kind of thing either way**:
the distinct sequences the reads showed there, each with how much support it had. A microsatellite
differs by carrying a repeat unit, not by needing a different sort of evidence — so one type serves
every kind of locus (§3).

This spec settles that shape. It does not settle how any particular kind of region becomes loci —
that is a **locus generator**, one per segment kind, each with its own spec.

**Goals.**

1. Define the locus type: the unit every generator produces and every later step consumes.
2. Define the `LocusGenerator` contract, so a new generator can be written — or an existing one
   replaced and re-measured — without touching the machinery around it.
3. Dispatch every typed region to a generator, with **no kind left unassigned**.
4. Account for every region and every base, including the ones no generator handles yet.
5. Be the first consumer of `SampleReads`, and report back what it needs (§8).

**Non-goals — could be goals, deliberately are not.**

- **Any generator's algorithm.** This doc owns the seam, not the work. If a statistical model appears
  here, it is in the wrong file.
- **Anything cross-sample.** One sample per stream. Reads are never merged across samples; loci are
  what merge later (§3).
- **Calling.** This step emits evidence, not genotypes — no likelihood, no prior, no QUAL.

**What it is not.** It is not a step with one implementation. ng is an experimental caller: several
generators are expected per kind, and the experiment is to swap one, hold the rest of the pipeline
fixed, and measure which calls better (§4).

---

## 2. Where it sits

```
reference ─▶ TypedRegionIterator ─▶ SampleLocusObservationsIterator ─▶ (calls)
                                          │
                          dispatch on the region's kind (§5)
          ┌───────────────┬───────────────┼───────────────┐
          ▼               ▼               ▼               ▼
      SsrSegment       Generic        SsrBundle       Satellite
          │               │               │               │
     [a generator]   [a generator]   [a generator]      NoLoci
      plugged in      plugged in      plugged in     (out of scope)
```

**This spec ships one generator — `NoLoci` (§5) — and nothing else.** Every kind that does real work
is handled by a generator with its **own** spec and module, plugged into the run. So *this* document
defines the socket and the wiring, not the work: it is deliberately independent of the STR generator,
the pileup, or any other. Which generator fills a kind's slot is the run's choice (§4); an unfilled
slot falls back to `NoLoci`.

| region kind | handled by | specced in |
|---|---|---|
| `SsrSegment` | a locus generator | [locus_generation_ssr.md](locus_generation_ssr.md) (the STR generator) |
| `Generic` | a locus generator | deferred (§11) |
| `SsrBundle` | a locus generator | deferred (§11) |
| `Satellite` | `NoLoci` — never handled | here (§5) |

**`Satellite` is the one kind this spec resolves.** It is out of scope for the whole caller, so it is
wired to `NoLoci` permanently. Every other kind's behaviour lives in that kind's own generator spec;
this document neither describes nor depends on it — including the STR generator that v1 supplies for
`SsrSegment`.

**Lazy, stream to stream.** The whole step is an iterator wrapping the typed-region iterator, not a
function you hand one region — because the typed-region walk is itself lazy over a windowed
reference, and a per-region entry point would hand the caller a reference window to manage.

**Loci come out in coordinate order** — the typed-region stream is already in that order and nothing
here reorders it. That order is part of the output contract, so later parallelism must reassemble
results rather than merely produce the same set.

---

## 3. The locus type

**What it is.** A locus: the stretch of genome it covers, and what one sample's reads showed there.

**One type for every kind of locus**, not a per-kind payload. A microsatellite and a candidate SNP
record the *same* thing — a table of distinct sequences the reads showed, each with its support.
Production already stores exactly that on both sides: `Vec<AlleleObservation>` for the generic path
([src/pileup_record.rs:138](src/pileup_record.rs#L138)) and `Vec<(sequence, count)>` for STR
([src/ssr/cohort/types.rs:60](src/ssr/cohort/types.rs#L60)) — the same table, STR carrying only the
count.

```rust
/// One sample's locus. Owned throughout: a cohort stage merges these across samples and a
/// future artifact writes them, so it must outlive every buffer it was built from.
pub struct SampleLocusObservations {
    /// The stretch of genome this locus covers. One base for a candidate SNP; several for
    /// an indel; the tract for a microsatellite.
    pub region: GenomeRegion,
    /// The reference bases over `region` — the REF sequence, and what a wider-span
    /// projection needs when samples merge (below).
    pub reference_bases: Box<[u8]>,
    /// The distinct sequences the reads showed, each with its support. **Observations, not
    /// alleles** — they become alleles when something calls them.
    pub observed_sequences: Vec<ObservedSequence>,
    /// Reads that covered this locus and produced no observation at all. *No coverage* and
    /// *coverage that said nothing* are different states, and only one means "look at the
    /// mapping". The per-reason breakdown is the generator's to report; that they existed
    /// belongs to the locus.
    pub reads_without_observation: u32,
    /// Reads a depth cap discarded. **Non-zero means the support counts are a subsample,
    /// not the depth** — a model reading a capped count as true depth is being misled.
    pub reads_discarded_by_cap: u32,
    /// What kind of locus this is, and the extras that kind carries. Explicit, so a
    /// downstream step reads the kind off the type rather than inferring it from which
    /// fields are populated.
    pub kind: LocusKind,
}

/// The kind of locus, plus whatever that kind adds to the shared fields above.
/// `non_exhaustive` because a kind's payload fills in as the generator that mints it is
/// written (§5).
#[non_exhaustive]
pub enum LocusKind {
    /// A SNP/indel candidate site. Its evidence is the observed alleles; no extras.
    Generic,
    /// A microsatellite tract — carries the motif and the flanks the read model needs.
    Ssr(SsrDetail),
    /// A repeat cluster with no clean flanks, coarser than a single tract. What it carries
    /// is the bundle generator's to decide (deferred, §11).
    SsrBundle,
}

/// What an `Ssr` locus carries — grouped, so a repeat's motif and flanks are present or
/// absent together, never half.
pub struct SsrDetail {
    /// The repeat unit.
    pub motif: Motif,
    /// The reference flanks the read model aligns the tract against — left and right of
    /// `region`. Lengths are `.len()`; clamped at contig ends, so the two can differ.
    pub left_flank: Box<[u8]>,
    pub right_flank: Box<[u8]>,
}

/// One distinct sequence the reads showed, with its support. The fields between `num_obs`
/// and `chain_ids` are the per-read moments the SNP filters read; the STR model reads only
/// `num_obs`.
pub struct ObservedSequence {
    /// The observed bases — allele content, in read coordinates.
    pub bases: Box<[u8]>,
    /// How much of the locus a read of this sequence spanned — the whole thing, or only
    /// part (below). **Part of the identity**: a `Complete` `ATAT` and a `PartialLeft(4)`
    /// `ATAT` are different evidence and stay separate entries.
    pub read_coverage: ReadCoverage,
    /// How many reads showed this sequence. The whole support on the STR path, and the one
    /// field every model on both paths reduces to.
    pub num_obs: u32,
    /// Reads on the forward strand — strand bias.
    pub num_fwd: u32,
    /// Σ per-read log-error (production's `q_sum`, [pileup_record.rs:47](src/pileup_record.rs#L47)) —
    /// the freebayes per-read error term.
    pub q_sum: f64,
    /// Σ MAPQ and Σ MAPQ²: the MAPQ Welch's-t multi-mapper filter recovers the mean and
    /// variance from these.
    pub mapq_sum: u32,
    pub mapq_sum_sq: u64,
    // Read-position-bias fields (production's `placed_left` / `placed_start`) are NOT here:
    // they are anchor-relative — meaningful on the generic path, degenerate on a tract — so the
    // pileup generator adds them when it is specced (§11), rather than the shared type carrying
    // generic-only fields nothing fills in v1.
    /// Phase-chain ids of the reads folded here — what lets a later step chain observations
    /// at neighbouring loci into a haplotype.
    pub chain_ids: Vec<ChainId>,
}

/// How much of the locus a single read spanned. `Complete` = it reached **both** borders;
/// the partials = it ran off its own end partway, carrying how many of the locus's positions
/// it reached from that border. (One read's span of the locus — not depth.)
pub enum ReadCoverage {
    Complete,
    PartialLeft(u16),
    PartialRight(u16),
}

impl SampleLocusObservations {
    /// Read depth at each position of `region`, in order — **derived, not stored**. A
    /// `Complete` observation counts its `num_obs` at every position; a `PartialLeft(n)` at
    /// the leftmost `n`, a `PartialRight(n)` at the rightmost `n`. Length is `region.len()`.
    pub fn num_obs_along_locus(&self) -> Vec<u32>;

    /// The observations a likelihood may score directly — the `Complete` ones. A partial is
    /// a lower bound and mis-scores as a short allele until step 7 models censoring
    /// ([locus_generation_ssr.md](locus_generation_ssr.md)), so reaching the partials is a
    /// deliberate act.
    pub fn complete_observations(&self) -> impl Iterator<Item = &ObservedSequence> + '_;
}
```

**The `kind` names the locus as a type, not a convention.** The evidence is uniform, so nothing else
in the record says whether this is a candidate site, a microsatellite, or a repeat cluster —
`LocusKind` says it outright, and a downstream step matches on it rather than guessing from which
fields are populated. This is `ng_step_interfaces.md`'s "STR-awareness is a *type*, not a convention",
now naming every kind rather than only repeat-or-not.

**How it relates to the region's kind.** Region typing produces a `RegionKind`
([src/ng/region_typing/mod.rs:170](src/ng/region_typing/mod.rs#L170)) — `SsrSegment`, `SsrBundle`,
`Generic`, `Satellite`. `LocusKind` is its output-side counterpart: the generator for each region
kind mints loci of the matching locus kind. **Two enums, two stages, deliberately not merged**, and
they differ in three ways worth stating:

- **`Satellite` has no locus kind** — it produces no loci (§5), so it drops out.
- **Cardinality differs** — one `Generic` *region* becomes many `Generic` *loci*; the region is the
  input, the locus the output.
- **Payloads differ in nature** — a `RegionKind` is *reference-defined* (what the sequence is:
  `SsrSegment` holds coordinates, motif, purity, no bases); a `LocusKind` is *evidence-defined* (what
  the generator gathered: `Ssr(SsrDetail)` holds the motif and the flank **bytes** it fetched). The
  region says where to look; the locus says what was found.

The correspondence is fixed by which generator handles which region (§5), not by the type system — a
generator sets its output kind — but in practice it is `SsrSegment → Ssr`, `Generic → Generic`,
`SsrBundle → SsrBundle`.

**`read_coverage` carries two facts at once, and on the STR path they are the same fact.** Whether a read
reached both borders of a tract or ran off its own end partway answers two questions together.
*Censoring:* a partial observation is a **lower bound**, the sequence at least this long, and a
likelihood scoring it as complete would read a long allele as short — `complete_observations()` is the
guard (step 7 owns the censored term, [locus_generation_ssr.md](locus_generation_ssr.md)). *Depth:*
summing each observation's covered positions — a `Complete` at every position, a `PartialLeft(n)` at
the leftmost `n` — gives `num_obs_along_locus()`, with nothing stored. The `u16` is that covered extent
in **locus** coordinates, a separate axis from `bases` (the allele, in **read** coordinates); both are
needed. On STR the two facts coincide cleanly because tracts do not overlap and a partial read is one
physical event.

**`num_obs_along_locus()` is *observation* depth, and only exact per locus.** Two caveats the paralog
filter (§11) must handle, not this step:
- **It omits covering-but-unobserved reads.** A read that spans tract positions but anchors no border
  yields no observation (it lands in `reads_without_observation`, a scalar with no positions), so it
  is invisible to the derived depth. The count is observed depth, not fragment depth — and it reads
  low exactly where tracts are long or hard.
- **Loci overlap on the generic path.** A multi-base deletion locus contains the per-position loci
  inside it, so composing one genome-wide profile by summing `num_obs_along_locus()` across loci would
  double-count the interior — the exact hazard production's walker guards against
  ([pileup_to_psp.rs:96](src/pileup/per_sample/pileup_to_psp.rs#L96)). Whether a deletion read counts
  at its deleted positions (span coverage) or not (production's base coverage), and how overlapping
  loci compose, is unsettled and belongs to the paralog filter. This spec neither decides nor tests
  it — it defines only the per-locus quantity.

**One sample's observations — a `CohortLocus` will *contain* it, not be it.** Region and evidence are
what *this* sample's reads showed; merging samples yields cohort loci carrying content no per-sample
locus has — a region unioned across samples, alleles re-projected onto it. So the cohort type composes
this one as `per_sample: Vec<Option<SampleLocusObservations>>`, production's shape
([src/var_calling/types.rs:39](src/var_calling/types.rs#L39)); §12 weighs the alternative of one type
carrying N samples. This is also why per-sample regions may differ and that is fine: where some plants
carry a five-base deletion the locus spans one base in a non-carrier and six in a carrier, and the
merge groups by overlap and projects onto the union span
([variant_grouping.rs:99](src/var_calling/variant_grouping.rs#L99),
[per_group_merger.rs:262](src/var_calling/per_group_merger.rs#L262)). The only requirement it places
on a generator: carry enough to be projected onto a wider span — for SNP/indel, the locus's own
reference bases (microsatellites never vary this way).

**One sample's observations — a `CohortLocus` will *contain* it, not be it.** Region and evidence are
what *this* sample's reads showed; merging samples yields cohort loci carrying content no per-sample
locus has — a region unioned across samples, alleles re-projected onto it. So the cohort type composes
this one as `per_sample: Vec<Option<SampleLocusObservations>>`, production's shape
([src/var_calling/types.rs:39](src/var_calling/types.rs#L39)); §12 weighs the alternative of one type
carrying N samples. This is also why per-sample regions may differ and that is fine: where some plants
carry a five-base deletion the locus spans one base in a non-carrier and six in a carrier, and the
merge groups by overlap and projects onto the union span
([variant_grouping.rs:99](src/var_calling/variant_grouping.rs#L99),
[per_group_merger.rs:262](src/var_calling/per_group_merger.rs#L262)). The only requirement it places
on a generator: carry enough to be projected onto a wider span — for SNP/indel, the locus's own
reference bases (microsatellites never vary this way).

**Owned, with no lifetimes** — this is what a cohort stage merges and an artifact writes. *Owned is
necessary but not sufficient:* `chain_ids` are meaningful only beside their neighbours, since the
merger resolves them *across* records to chain compound haplotypes
([src/var_calling/per_group_merger.rs:1296](src/var_calling/per_group_merger.rs#L1296)). They are
kept here because phasing an STR allele against a neighbouring SNP is as real as phasing two SNPs —
with one cost worth recording, learned the expensive way: in the production cohort path chain-ids were
**~31% of peak live heap**, almost all of it REF observations the merger then skips. A generator that
fills this field pays that unless it skips REF, as production learned to.

---

## 4. The `LocusGenerator` contract

**What a generator is.** Something that takes one typed segment plus a way to read a sample's reads,
and yields that sample's loci there — zero, one, or many, **one at a time**.

```rust
/// Generates a sample's loci from one segment of kind `S`.
///
/// `S` is the segment payload the generator consumes — `SsrSegment` for the STR
/// generator. It is a parameter on the *contract* rather than a fixed type inside each
/// implementation, so that two generators for the same kind stay interchangeable (below).
pub trait LocusGenerator<S> {
    /// Start a new segment: reset progress. Does no gathering and cannot fail.
    fn begin_segment(&mut self, region: GenomeRegion);

    /// The next locus of the segment begun, or `None` once it has no more. Called
    /// repeatedly with the same `segment` until it returns `None`; returning `None`
    /// immediately is a normal outcome, not a failure.
    ///
    /// `&mut self` because a generator owns reusable scratch (alignment matrices,
    /// sampling buffers) that must not be reallocated per segment.
    fn next_locus(
        &mut self,
        segment: &S,
        reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError>;
}
```

**One locus at a time, not a filled `Vec` — decided, and this is a memory decision.** A `Generic`
region is whatever lies between typed repeats, so it can run to tens or hundreds of kilobases, and
the pileup emits a locus per covered position. A 100 kb segment is therefore ~100,000 loci, each
owning a vector of observations with their bases and chain-ids — tens of megabytes resident for **one
segment**, against a caller whose headline claim is memory efficiency. Handing back a filled `Vec`
would make that the default.

It also matches how the work actually happens: production's walker advances position by position over
an active read set and closes each locus as it passes, so it never needs the whole segment before it
can hand over the first locus. The `Vec` would be a buffer imposed by the interface, not by the
algorithm.

*Two consequences worth knowing, because they shaped the signature.* **The generator borrows the
segment on each call instead of storing it.** Either way it is a borrow — `&S`, a pointer, never a
copy — but a generator that *held* the borrow would need a lifetime parameter, and that parameter
then spreads to the dispatcher, the iterator, and anything owning a set of generators. Borrowing per
call costs nothing at run time and keeps the generator a plain lifetime-free type. (Taking the
segment **by value** is the option to avoid: `SsrSegment` owns a boxed contig name, so that would
allocate on every call.) And **`Box<dyn Iterator<Item = SampleLocusObservations>>` was rejected** — the
read-ingestion work already recorded why it chose an enum over a box for its read streams (a box is
opaque to the optimiser and stops the chain inlining into the consumer's loop), and here it would
also cost an allocation per segment while buying nothing these two methods do not.

**Why the segment type is a parameter on the trait — decided.** Two alternatives were live, and the
choice matters because it decides whether one generator can be swapped for another at all — which is
the experiment this module exists to serve.

- *Every generator takes the whole `TypedRegion`.* One flat contract, but a generator then receives a
  container that might hold any kind and must handle a case the dispatcher makes impossible — an
  untestable branch in every generator, forever. **Rejected** for that dead code.
- *Every generator declares its own input type* (an associated type). Precise, and it re-creates the
  exact defect one step upstream: `ReadPreparer` declares its input as an associated type, which is
  why `CallerRecipe`'s `read_preparer: Box<dyn ReadPreparer>` (`ng_step_interfaces.md:418`, sketched
  there as the intended field) does not in fact compile as written — a trait object must name every
  associated type. When each implementation keeps its input private, "some generator, chosen at run
  time" is not expressible — and that is precisely what swapping one implementation for another
  requires. **Rejected**: it buys precision and sells the swap.
- *The segment type is a parameter on the contract* — chosen. `LocusGenerator<SsrSegment>` and
  `LocusGenerator<SsrBundle>` are distinct contracts, each generator's signature names the payload it
  actually eats, and `Box<dyn LocusGenerator<SsrSegment>>` stays swappable.

The reason this is the right shape rather than a trick: **generators are only ever alternatives to
others handling the same kind.** Nobody substitutes the satellite generator for the STR one — they
take different inputs. So interchangeability is needed *within* a kind, which is exactly what this
gives, and the parameter costs nothing at run time.

**A generator holds its own accessors.** The reference accessor, the read preparer, the scratch
buffers are fields, not arguments — the convention `read_preparation.md` §3 already set. It is also
what keeps the contract free of the associated types that would break the swap.

**Reads are not capped here.** A generator that needs to bound per-locus depth owns that decision and
its determinism consequences; the STR generator does
([locus_generation_ssr.md](locus_generation_ssr.md) §4).

---

## 5. The dispatcher, and the generator that makes nothing

The dispatcher is a `match` over the region's kind, handing each branch's payload to the generator the
run supplies for it (§4). **Every branch resolves to a generator** — the match is exhaustive and no
arm is a `TODO`. A kind with no real generator supplied falls back to `NoLoci`, so *this* spec, run on
its own, produces nothing: the work arrives when a run plugs generators in.

That is what the empty generator is for:

```rust
/// A generator that produces no loci and counts what it passed over. One implementation
/// covers every kind, because it ignores the segment entirely.
pub struct NoLoci {
    pub reason: UnhandledReason,
}

/// **Why** a kind produces no loci. Not cosmetic: it separates a boundary we chose from a
/// gap we have not filled, and those answer different questions (§7).
pub enum UnhandledReason {
    /// Deliberately outside the caller's scope — e.g. satellite arrays.
    OutOfScope,
    /// No generator written yet. Temporary by construction.
    NotImplemented,
}

impl<S> LocusGenerator<S> for NoLoci { /* counts the segment and its bases, emits nothing */ }
```

**Why this beats an unhandled arm.** A region kind nobody has built for is otherwise either a silent
`_ => {}` or a `TODO`, and both are holes that stop being visible. Routing it to a generator makes
"we produce nothing here" a **configuration with a reason attached**, counted like everything else,
and makes plugging in a real generator later a one-line change at the dispatch site rather than a
new code path.

**Two kinds of nothing, and they must not be added together.** `NoLoci` carries a reason because a
kind can reach it two ways. `Satellite` is `OutOfScope` — a tandem array longer than `max_str_len` is
not what this caller is for, a boundary we chose, permanent. A kind whose real generator simply is not
supplied in this run is `NotImplemented` — a gap we have not filled, temporary. Reporting both as
"skipped" would conflate the two, and "how much genome does this caller not cover *yet*?" is a
question the manager needs answered separately from "what will it never cover?" In v1 the run supplies
the STR generator, so `SsrSegment` does real work while `Generic` and `SsrBundle` fall to
`NotImplemented` and `Satellite` to `OutOfScope`.

**Module layout.** `src/ng/locus_generation/` — a folder, per `module_layout.md` principle 1: a step whose
implementations compete keeps its trait and every implementation side by side.

```
src/ng/locus_generation/
  mod.rs      – SampleLocusObservations, ObservedSequence, LocusGenerator, the dispatcher, NoLoci
  ssr.rs      – the STR generator (locus_generation_ssr.md)
  pileup/     – the generic generator (deferred, §11)
```

*This reverses an earlier call in this design* that the step was a file, on the grounds that it had no
competing implementations. It has: ng is an experimental caller, and alternative generators per kind
are the point.

---

## 6. The public surface

An iterator over a typed-region stream, mirroring `TypedRegionIterator`
([src/ng/region_typing/mod.rs:789](src/ng/region_typing/mod.rs#L789)) — same lazy shape, same running
counts, same in-stream fatal error.

```rust
pub struct SampleLocusObservationsIterator<T, G> { /* … */ }

impl<T, G> Iterator for SampleLocusObservationsIterator<T, G>
where
    T: Iterator<Item = Result<TypedRegion, TypedRegionError>>,
{
    type Item = Result<SampleLocusObservations, LocusGenerationError>;
}

impl<T, G> SampleLocusObservationsIterator<T, G> {
    /// `regions` is the typed-region stream, `reads` the sample's reads, `generators` the
    /// set the dispatcher routes to (§5).
    pub fn new(regions: T, reads: SampleReads, generators: G, config: LocusConfig) -> Self;
    /// The running tally — current at any point, final once the stream is exhausted.
    pub fn counts(&self) -> &LocusCounts;
}
```

The iterator holds **no buffer of loci**: it calls `next_locus` until the generator returns `None`,
then pulls the next region and calls `begin_segment`. That is what §4's streaming contract buys — one
locus resident at a time no matter how many a segment yields, where a filled-`Vec` contract would
have forced the same `queue` shape `TypedRegionIterator` needs.

**Errors: the item is a `Result`**, per `read_filtering.md` §5. A fatal condition yields
`Some(Err(_))` once, then `None`, so `?` makes it un-ignorable rather than a silent end of stream.

```rust
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

The line is the established one: **a read that yields no observation is a tallied per-read outcome**,
never an error. An error means the run is broken.

---

## 7. Config and counts

**Config** mirrors `ReadFilterConfig` / `TypedRegionConfig`: one field per active knob, no dormant
levers, `Option` meaning absent rather than a sentinel, defaults as named `pub const`s. The shared
config holds only what the dispatcher itself needs; **each generator owns its own knobs** and takes
them at construction, so that adding a generator does not widen a shared struct.

**Counts** follow the "no silent caps" discipline — every region and every base is accounted for.

```rust
pub struct LocusCounts {
    pub regions_in: u64,
    pub loci_emitted: u64,
    /// Regions that produced no loci because no generator exists yet, and their bases.
    pub unhandled_not_implemented: u64,
    pub unhandled_not_implemented_bp: u64,
    /// Regions deliberately outside scope, and their bases.
    pub unhandled_out_of_scope: u64,
    pub unhandled_out_of_scope_bp: u64,
}
```

The base counters are what make "how much genome does this caller not cover, and how much of that is
temporary?" answerable from the counts alone. `SsrBundle` exists as a type precisely so its bases
cannot become a hole nobody accounts for; counting them — and satellites' — is the other half of
that. Per-generator counts (reads fetched, reads dropped, and why) belong to the generator and are
reported by it, since they have no meaning shared across kinds.

---

## 8. What this step needs from `SampleReads`

This is the first consumer of the read-ingestion interface. Its code is built and these requirements
were verified against the `ng-read-ingestion` worktree; they stand as feedback to that work whether or
not it has merged.

1. **A per-query reference accessor that costs nothing to make — the one real gap.**
   `reads_in_region` takes a factory `F: FnMut() -> R` and calls it once per file per query
   (`src/ng/read/input/mod.rs:440`). A per-locus generator issues on the order of 10⁶ queries. Neither
   accessor is cheap to construct: `WindowedRefSeq::new` holds no resident contig, so a fresh one per
   query re-reads its window at every locus and defeats the sliding buffer; and `ResidentRefSeq::new`
   takes a `ContigList` by value ([src/ng/ref_seq.rs:367](src/ng/ref_seq.rs#L367)), which is a
   `Vec<ContigEntry>` of `String` names ([src/fasta/mod.rs:60](src/fasta/mod.rs#L60)) — one `String`
   allocation per contig, per query, per file. On a scaffold-rich assembly that is the dominant cost
   of the step.
   **Suggested fix, small:** `impl RefSeq for Arc<T> where T: RefSeq` (and the same for `RawRefSeq`).
   The factory then returns `Arc::clone(&shared)` — one atomic increment. The blanket impl on `&T`
   that `ref_seq.md` ruled out is not needed; `Arc` sidesteps the lifetime problem it hit.
2. **Backward queries must stay legal.** They are today — the order guard lives in the per-query
   `OrderVerified`, not on the handle. This step sweeps forward, but a future batched sweep will
   re-query, so nothing should later "tidy" that state onto the handle.
3. **Depth capping stays out.** Confirming, not requesting: capping is a generator's decision (§4).
   Two caps would compose into a silently different kept set.
4. **Reference and `@SQ` must match exactly — resolved (owner, 2026-07-22), and already enforced
   upstream.** The reads were mapped against this reference, so each file's `@SQ` list must equal the
   reference contig list — names, lengths, order — which read ingestion checks at file open
   ([src/ng/read/input/open_bam.rs:243](src/ng/read/input/open_bam.rs#L243)). A contig in the
   reference but absent from a sample's `@SQ` is therefore rejected there and never reaches this
   step, so there is no per-locus "zero reads vs. error" case for this step to handle. *This
   supersedes the earlier "zero reads" lean*, which assumed header subsets were ordinary; the
   ingestion contract requires exact equality instead.

**Nothing here blocks this step.** Item 1 is a cost, not a correctness problem, and can land on
either branch.

---

## 9. Cross-cutting concerns

**Memory.** **One locus resident at a time**, plus whatever the active generator holds internally —
which is what §4's `next_locus` contract exists to guarantee, since a generic segment can cover
hundreds of kilobases and yield a locus per position. The generators fetch forward, which drives
`evict_before` behind the current region, so the reference window stays proportional to one region
plus the next read's forward reach — not to a contig. What a generator holds *internally* is its own
to bound: the STR generator's reservoir cap is that bound
([locus_generation_ssr.md](locus_generation_ssr.md) §4).

**Determinism.** Output is a deterministic function of (reference, config, sample files). The
dispatcher contributes no nondeterminism of its own; where a generator samples or aggregates, it owns
the guarantee and states it.

**Concurrency.** Regions are independent, so the work parallelises naturally. Two costs: results must
be **reassembled in coordinate order** before emission (§2), and a shared `WindowedRefSeq` is `Send`
but not `Sync`, so each worker needs its own accessor. Neither is paid in v1; the interface precludes
neither.

---

## 10. Reuse over rewrite

| what | existing code | ng reuse |
|---|---|---|
| region stream | [src/ng/region_typing/mod.rs:789](src/ng/region_typing/mod.rs#L789) | consume as-is; mirror its iterator + counts shape |
| read access | `SampleReads::reads_in_region` | reuse as-is (§8) |
| locus/cohort split | [src/var_calling/types.rs:39](src/var_calling/types.rs#L39) | model for `CohortLocus` composing `SampleLocusObservations` (§3) |
| cross-sample reconciliation | [src/var_calling/variant_grouping.rs:99](src/var_calling/variant_grouping.rs#L99) | **not this step** — the merge groups by overlap (§3) |

---

## 11. Deferred, with a recommended home

- **The pileup generator** for `Generic` regions — `src/ng/locus_generation/pileup/`, its own spec. It also owns
  the `chain_ids` question of §3 (carry cross-locus read identifiers, or drop compound haplotypes) and
  adds the read-position-bias fields (`placed_left`, `placed_start`) to `ObservedSequence` — generic-only,
  omitted from the shared type because a tract makes them degenerate (§3).
- **The windowed statistics** — mean depth and GC over a window centred on a locus, the paralog
  filter's inputs. Computed downstream by sliding a window over `num_obs_along_locus()` (§3), not
  stored here. **Home: the paralog filter** (step 11a). It inherits every depth caveat §3 names: that
  the derived count is *observation* depth (covering-but-unobserved reads are invisible); whether GC
  covers the window's covered positions (production) or all of them; and — on the generic path —
  span-vs-base coverage and how to compose one profile from overlapping loci without double-counting a
  deletion's interior.
- **A generator for `SsrBundle`** — it records at least depth per position, so a bundle's bases stop
  being a hole in the depth profile the windowed statistics slide over. Whether it also emits
  observed sequences is open. Routed to `NoLoci { NotImplemented }` meanwhile, with its bases counted,
  so the gap stays visible until it is filled. **Home: `src/ng/locus_generation/`**, beside the STR generator.
- **The cohort merge** — many samples' loci into cohort loci, by overlap (§3). **Home: the cohort
  spec.** Its one requirement on this step is that a locus carry enough to be projected onto a wider
  span.
- **The cohort artifact.** ng adopts production's two-level shape (per-sample stage → artifact →
  cohort gather) but **not the `.psp` file yet**: the seam is the load-bearing thing, in-memory
  first, serialization when memory forces it or when these types stop churning. This reversed
  `module_layout.md`'s "no `.psp` split, single-phase", and that doc has been revised to match
  (*Crate boundary and the port-back*).
- **The assembly check.** `SampleReads` only forwards the `@SQ` MD5 tags; someone must call
  `check_assembly` after joining the reference verification, before output is committed. **Home: the
  pipeline driver** (unspecced). Noted because nobody else will, and a wrong-assembly run otherwise
  wastes its whole work.
- **Parallelism** — deferred whole (§9).

---

## 12. Open questions

- **One locus type carrying N samples, instead of a single-sample type the cohort composes? — open,
  leaning to compose (§3).** The unified type buys one word for the concept in both stages, no
  parallel cohort hierarchy, and one serialization for both — and for a microsatellite-only caller it
  is plainly better, since tract regions never differ between samples; production's `CohortLocus`
  ([src/ssr/cohort/types.rs:79](src/ssr/cohort/types.rs#L79)) is essentially that design. Against it:
  a cohort locus carries merged content no per-sample locus has, so those fields would sit empty
  through the whole per-sample stage; where sample regions differ, one type must hold both the group
  span and each sample's own; and the sample axis is redundant in a stage whose stream has exactly
  one sample. *What would settle it:* whether ng's generic path projects alleles onto a group span,
  as production does, or keeps per-sample allele sets and reconciles at call time. **Confirm when the
  pileup generator is designed;** it does not block v1.
- **How to carry the kind-specific fields — resolved: an explicit `LocusKind` enum.** This was open
  while the only kind-specific field was `motif`, where a lone `Option<Motif>` was enough. Two things
  closed it: the read-model flanks arrived (a second and third repeat-only field, so they grouped into
  `SsrDetail`), and the unified evidence left nothing else naming the kind — a consumer could no longer
  tell `Generic` from `SsrBundle` at all. So the kind is now a first-class `LocusKind` naming every
  kind, each variant carrying its own extras (§3) — which also restores `ng_step_interfaces.md`'s
  "STR-awareness is a type, not a convention". Purity joins `SsrDetail` when the interrupted-repeat
  work needs it, with no further type change.
- **Does a generator ever need to see more than one region at a time? — resolved: no (owner,
  2026-07-22).** Each genomic segment is resolved independently, so no generator needs cross-segment
  context: the typed-region tiling is gap-free, the neighbouring region is generated next, and the
  merge reconciles by overlap. A pileup wanting to call across a segment boundary (an indel spanning
  the join between a `Generic` stretch and an adjacent tract) reconciles at the merge, not by
  widening this contract. This confirms §4's one-segment contract.
- **A contig in the reference but absent from a sample's `@SQ` header — resolved: cannot occur.**
  Reference and `@SQ` must match exactly, enforced at read ingestion — stated in full as §8 item 4.
  Noted here too so it is visible to a reader scanning the open questions.

---

## 13. Acceptance test

The shape is done when the dispatch and the accounting are provable, independent of any generator's
output. The STR generator's own test is in its spec.

1. **Every region kind is assigned.** The `match` is exhaustive, so the compiler proves it; a test
   asserts each kind reaches the generator it should, using `NoLoci` for the unbuilt ones.
2. **Nothing is unaccounted for.** Over a fixture: `regions_in` equals regions producing loci plus
   `unhandled_not_implemented` plus `unhandled_out_of_scope`, and the two base counters sum to the
   bases those regions cover.
3. **The two kinds of nothing stay separate.** A fixture with both a satellite and a generic region
   reports them in different counters — the check that §5's distinction survives contact with code.
4. **Order is preserved.** Emitted loci are in coordinate order across a multi-region, multi-kind
   fixture, including where one region yields several loci.
5. **Depth derives correctly from read-coverage.** `num_obs_along_locus()` yields `region.len()`
   values; a `Complete` observation raises every position, a `PartialLeft(n)` only the leftmost `n` —
   checked on a fixture mixing complete and partial observations, so §3's read-coverage-to-depth rule
   is exercised.
6. **A partial and a complete of the same bases stay distinct.** Two observations with identical
   `bases` but different `read_coverage` are separate entries — the check that a lower bound cannot be
   silently folded into a called allele.
