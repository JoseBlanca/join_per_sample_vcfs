# ng — common types and step interfaces

*Status: architecture draft (2026-07-10), companion to
[`../spec/ng_proposal.md`](../spec/ng_proposal.md). This doc defines the **shared
vocabulary** (domain newtypes and the data that flows between steps) and the **step
traits** that let one implementation of a step be swapped for another and measured
end-to-end. It is a starting point to iterate on, not a frozen API. Naming follows
[`ai/skills/rust-code-review/code_review/naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md):
domain nouns for types, verbs for functions, **newtypes for domain scalars**,
`str` in prose ↔ `ssr` in code.*

## Why interfaces, not just modules

The ng plan (spec §1) is to try several implementations of each step and keep the
best. That only works if a step is a **trait with a narrow contract**: swap the
impl, hold the rest fixed, measure the whole pipeline against gold/silver/synthetic
(spec §2). So this doc's two jobs are:

1. **Common types** — a shared vocabulary every step speaks, so a candidate set from
   implementation *A* of step 6 is consumable by implementations *X, Y, Z* of step 7.
   Newtypes here do double duty: they document intent at every call site, and they
   make the compiler reject unit/quantity transpositions (e.g. a base-pair delta
   used where a repeat-unit delta is expected — the exact class of bug that produced
   the ~1000 mis-scored benchmark calls).
2. **Step traits** — one per swappable step, so implementations are interchangeable
   behind the same signature.

Two traits already exist and this design builds on them, not around them:
`ReadLikelihoodModel` ([`src/ssr/cohort/read_model/`](../../../../src/ssr/cohort/read_model/))
is step 7 for STR; `GenotypeEmModel`
([`src/var_calling/posterior_engine.rs`](../../../../src/var_calling/posterior_engine.rs))
is step 9. The ng work generalises the rest to the same shape. (The crate-level
`src/genotype_em/` hoist for that trait is still pending — it lives in
`posterior_engine.rs` today.)

---

## 1. Domain newtypes — the vocabulary

All are `#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]`
(float-carrying ones drop `Eq/Ord/Hash`), expose `.get()`, and live in one module
(`ng::units`) so the whole pipeline shares them. The point of separate types for
same-primitive quantities is that **the compiler refuses to mix them**.

### Coordinates & sequence

```rust
pub struct ContigId(pub u32);       // index into the contig table
pub struct Position(pub u64);       // 1-based reference position (matches VCF/SAM/IGV and the production engine)
pub struct GenomeRegion { pub contig: ContigId, pub start: Position, pub end: Position } // 1-based inclusive [start, end]
pub struct SsrMotif(Box<[u8]>);     // the repeat unit, phase-faithful (STR)
pub struct SsrPeriod(u8);           // motif length in bp, the repeat period (STR); >= 1, validated
```

### Lengths & deltas — the unit distinction that bit us

Base pairs and repeat units are **different quantities**; conflating them is the
representation bug of spec §2. Give them different types.

```rust
pub struct Bp(pub u32);             // a length in base pairs (allele/tract/flank/read) — generic
pub struct BpDelta(pub i32);        // signed length change in base pairs — generic
pub struct SsrRepeatUnits(pub u16); // a tract length in whole repeat units (STR)
pub struct SsrUnitDelta(pub i32);   // signed length change in whole repeat units (STR)
```

`Bp` / `BpDelta` stay **unprefixed**: base pairs are the generic length currency both
the SNP/indel and STR paths speak. Only the *repeat-unit* quantities are STR-specific,
so only they carry `Ssr`. The `BpDelta` → `SsrUnitDelta` conversion is thus exactly the
generic→STR bridge, and it is **only** valid with an `SsrPeriod` in hand, so it lives on
`SsrPeriod` (`period.units_of(bp_delta) -> Option<SsrUnitDelta>`, `None` when the bp
change is not a whole-unit multiple). That one method is the choke point the benchmark
bug slipped through when it was implicit.

### Depth, counts, identity

```rust
pub struct ReadDepth(pub u32);      // reads at a locus for one sample
pub struct AlleleObsCount(pub u32); // reads supporting one allele
pub struct ReadWeight(pub f64);     // a read's likelihood weight (>= 0): pooling count or partial fraction
pub struct SampleId(pub u32);       pub struct SampleGroupId(pub u32);
pub struct LocusId(pub u64);        pub struct AlleleId(pub u16);   // index within a locus' candidate set
```

### Probability & quality — keep log-space explicit

```rust
pub struct LogProb(pub f64);        // natural-log probability; the internal currency
pub struct Phred(f32);              // -10 log10 p, for I/O (QUAL, GQ); >= 0, validated
pub struct BaseQual(pub u8);        pub struct MapQual(pub u8);
```

`LogProb` vs `Phred` are different scales with different bases; making them distinct
types stops the classic "added a Phred to a ln-prob" error. Conversions are named
functions (`Phred::from_log_prob`, `LogProb::ln`), never `as`.

### Population genetics — constrained, so validated

```rust
pub struct AlleleFreq(f64);         // in [0, 1]   — private field, checked constructor
pub struct InbreedingF(f64);        // Wright's F_IS, in [0, 1)
pub struct Theta(f64);              // Watterson/Ewens diversity, > 0

impl AlleleFreq {
    /// The only constructor — enforces the [0, 1] invariant, so every use site can
    /// trust it. `InbreedingF`, `Theta`, `SsrPeriod` (>= 1), `Phred` (>= 0) follow suit.
    pub fn try_new(p: f64) -> Result<Self, DomainError> {
        (0.0..=1.0).contains(&p).then_some(Self(p)).ok_or(DomainError::AlleleFreq(p))
    }
    pub fn get(self) -> f64 { self.0 }
}
```

**Validation on constrained newtypes.** Types with a domain range (`AlleleFreq`
[0,1], `InbreedingF` [0,1), `Theta` > 0, `SsrPeriod` ≥ 1, `Phred` ≥ 0) **hide their
field and construct through a checked constructor**, so an illegal value is
*unrepresentable* rather than merely discouraged. Policy by the value's source:

- **Untrusted input** (a parsed VCF / catalog field, a CLI arg): `try_new -> Result`
  and fail loudly — never silently coerce.
- **Internally computed** (an EM frequency, an estimator's `F`): the value is in range
  *by construction*, so a `new` that `debug_assert!`s the bound and clamps only a
  float-epsilon overrun (an EM yielding `1.0000001` → `1.0`) is fine — but a *gross*
  out-of-range value stays a loud bug, never a silent clamp. (The existing code already
  clamps `F` to a `0.99` ceiling; the type is where that intent should live.)

Unconstrained newtypes (`Bp`, `ReadDepth`, the ids, `BpDelta`) keep `pub` fields —
there is no invariant to protect, so a checked constructor would be ceremony.

> **Ergonomics.** Wrappers cost a `.get()` at arithmetic sites; for the few types with
> heavy arithmetic (`LogProb`, `BpDelta`) implement the needed operators so hot code
> reads naturally. Do **not** blanket-impl `Deref` to the primitive — that would
> re-open the transposition hole the newtype exists to close.

---

## 2. The data that flows between steps (wire types)

These are the payloads passed step→step; every step impl consumes and produces them.

```rust
/// One sample's evidence at one locus — the per-locus, per-sample unit, and the
/// currency of the **locus stream** (spec §1, *The locus stream*). A locus enters the
/// stream from one of two mints, but downstream steps consume this type identically —
/// that uniformity is what keeps SNP / indel / STR at one level:
///   • an **STR** locus is reference-defined (the catalog blesses a tract 1:1);
///   • a **generic** locus is data-defined (the `pileup/` module walks a non-STR
///     stretch, splits it into loci, and gathers each one's reads).
pub struct LocusEvidence<'a> {
    pub sample: SampleId,
    pub reads: &'a [LocusRead],   // step-2 prepared reads (STR path) / pileup output (generic)
    pub depth: ReadDepth,
}

/// A read after step 2 (realigned / delimited). The STR path stores the observed
/// tract bytes; the generic path stores the aligned span.
pub struct LocusRead {
    pub map_qual: MapQual,
    pub observed: Observation,       // enum: generic aligned allele, or an STR tract
    pub weight: ReadWeight,          // scales this read's likelihood: 1.0 normally; a
                                     // pooling count (>=1); or a partial-observation
                                     // fraction (<1) — the freebayes lever we mean to test
}

/// The router's verdict for a locus (spec step 3). This is what makes STR-awareness
/// a *type*, not a convention: only `Ssr` carries a motif/period/borders.
pub enum LocusKind {
    Generic(GenomeRegion),
    Ssr(SsrLocus),
}
pub struct SsrLocus {
    pub id: LocusId,
    pub span: GenomeRegion,           // the tract, without flanks
    pub motif: SsrMotif,
    pub period: SsrPeriod,
    pub ref_units: SsrRepeatUnits,
}

/// The candidate allele set at a locus (step 6). REF is always present at `ref_idx`.
pub struct AlleleCandidates {
    pub alleles: Vec<Allele>,        // sequence-resolved
    pub ref_idx: AlleleId,
    pub admit: Admission,            // enum { Ok, LowDepth, NotPeriodic, TooManyAlleles }
}

/// One individual's genotype: a multiset of allele-ids of size = ploidy (order-free,
/// repeats allowed). The field is PRIVATE so ploidy is not baked into the surface —
/// diploid is simply `len() == 2`. A polyploid impl changes only the constructor and
/// this type; the traits that consume `Genotype` never see the difference.
pub struct Genotype(Box<[AlleleId]>);   // sorted; len() == ploidy
impl Genotype {
    pub fn alleles(&self) -> &[AlleleId] { &self.0 }
    pub fn ploidy(&self) -> u8 { self.0.len() as u8 }
    pub fn is_homozygous(&self) -> bool { self.0.iter().all(|a| *a == self.0[0]) }
}

/// Per-sample genotype log-likelihoods at a locus (step 7 → step 9). The single
/// read-level quantity every prior/posterior/emission model consumes (spec's `Lg`).
pub struct GenotypeLikelihoods {
    pub genotypes: Vec<Genotype>,    // the enumerated genotypes over the candidates (diploid v1)
    pub log_lik: Vec<LogProb>,       // parallel to `genotypes`, per sample handled by caller
}

/// A called genotype for one sample at one locus.
pub struct GenotypeCall {
    pub genotype: Genotype,          // opaque multiset; diploid v1
    pub quality: Phred,              // GQ
}

/// The caller's model parameters — per-group error + stutter, sample-group
/// assignment, per-sample F, contamination. Estimated once by the pre-pass (spec
/// step 4, the "first caller") and held fixed: consumed as *inputs* by steps 7-9,
/// never re-fit downstream.
pub struct ModelParams {
    pub error_by_group: Vec<PerBaseError>,     // ε per sample-group
    pub stutter_by_group: Vec<StutterModel>,   // STR
    pub group_of_sample: Vec<SampleGroupId>,
    pub f_by_sample: Vec<InbreedingF>,
    pub contamination: Option<ContaminationModel>,
}
```

---

## 3. The step traits — swappable implementations

One trait per swappable step. Signatures are illustrative; the contract (what goes
in, what comes out, what invariants hold) is the deliverable. Steps 1 and 10–13 that
are less about competing algorithms are sketched briefly; the algorithmic hot-spots
(2, 3, 4, 6, 7, 8, 9, 11) get the real interfaces.

### Step 1 — read filtering

Its "real interface" lives in a dedicated companion, [`read_filtering.md`](read_filtering.md)
(spec: [`../spec/read_filtering.md`](../spec/read_filtering.md)) — not a swappable-algorithm
step but a fixed prelude, so it gets its own doc rather than a trait sketch here. In short: a
`RecordSource` fills one reused record buffer, a nine-filter cascade runs #1–#6 (flag/MAPQ)
*before* decode and #7–#9 after, and a `ReadFilter` iterator yields the surviving
`MappedRead`s while carrying a running `ReadFilterCounts`. It reuses the
[`bam/alignment_input.rs`](../../../../src/bam/alignment_input.rs) predicates and constants.

### Step 2 — read preparation / realignment
```rust
pub trait ReadPrep {
    /// Realign/delimit one read against the reference window (STR impls also take the
    /// tract so gaps can be tract-aware). Returns None if the read is unusable here.
    fn prepare_read(&self, read: &MappedRead, window: &RefWindow) -> Option<LocusRead>;
}
```
*Impls to bench:* trust-mapper+left-align (freebayes-style), local reassembly
(GATK-style), per-read pair-HMM to the tract (ours-STR / HipSTR-style).

### Step 3 — the locus router
```rust
pub trait LocusRouter {
    /// Classify a reference region: is it an STR (with motif/period/borders) or generic?
    fn route_locus(&self, window: &RefWindow) -> LocusKind;
}
pub trait LocusSource {
    /// Enumerate the candidate loci to genotype — a catalog iterator, or active regions.
    fn loci(&self) -> impl Iterator<Item = GenomeRegion>;
}
```
*Impls to bench:* fixed STR catalog (ours), active-region detection (GATK), the
catalog+discovery hybrid (spec §step-3). Keeping `LocusSource` and `LocusRouter`
separate lets a data-driven source feed the same reference-based router.

**Step 3 is the head of the locus stream** (spec §1, *The locus stream*). It segments
the genome into STR / non-STR stretches; from there **two mints feed one stream** of
`LocusKind`-carrying `LocusEvidence`:

- an **STR** stretch is blessed as a locus 1:1 — reference-defined, minted before any
  read is read;
- a **non-STR** stretch is handed to the **`pileup/`** module, which walks it, splits it
  into loci, and gathers each one's evidence — data-defined.

`pileup/` is **not** a step trait here: it is infrastructure (the reused
`pileup/walker/`), the generic-path counterpart of the STR read-class/gather, and it has
no bake-off surface of its own. What *is* still swappable inside it — the non-STR
locus/window definition (single position vs active-region vs haplotype-fixpoint) — is the
`LocusSource`/active-region axis above. See `module_layout.md` (*The locus stream*) for
the open subsume-or-compose question.

### Step 4 — the rough caller, the per-sample summary, and the cohort estimator

Two levels, mirroring the production `.psp` engine: per-sample estimation runs in
Stage 1 and lands in the `.psp`; a **cohort-gather** step then computes the panel-level
parameters before calling. (That cohort-gather step was missing from the earlier step
list — see the note below.)

```rust
pub trait Caller {
    /// Call one sample's genotype at one locus. The pre-pass drives a cheap
    /// implementation that must NOT depend on the parameters being estimated — it uses
    /// a crude base-quality ε̂, not the fitted ε. `None` where it makes no call; the call
    /// carries its own quality, so "confident" is a property of the call, not a return type.
    fn call(&self, evidence: &LocusEvidence) -> Option<GenotypeCall>;
}

/// Per sample (Stage 1). From this sample's confident rough calls, compute the
/// per-INDIVIDUAL parameters (e.g. its inbreeding `F`) and the sufficient statistics the
/// cohort step will need. In production this is written into the `.psp` after the genomic
/// blocks — the sample's summary section.
pub trait SampleSummarizer {
    fn summarize(&self, confident: &[ConfidentGenotype]) -> SampleSummary;
}

/// Once, cohort-wide, before calling — THE STEP WE HAD MISSED. Gather every sample's
/// summary and estimate the parameters that need the whole panel: sample-groups,
/// per-group chemistry (error/stutter), the SFS `Theta`, contamination. Assembles the
/// final `ModelParams`, folding in each sample's `F` from its summary.
pub trait CohortEstimator {
    fn estimate(&self, summaries: &[SampleSummary]) -> ModelParams;
}
```
*Contract.* **(1) Per sample:** a cheap concrete `Caller` — a `RoughHetCaller` (SNP) or
the STR 1-vs-2 BIC gate — makes confident rough calls, and `SampleSummarizer` turns them
into the per-individual params (`F`) + sufficient stats, frozen into the `.psp`. **(2)
Cohort-wide, once:** `CohortEstimator` gathers all summaries and estimates the panel-level
params (sample-groups, per-group chemistry, SFS `Theta`) into the final `ModelParams`.
Only confident calls teach the parameters; **"Rough" is the concrete struct's quality,
not the trait's.** Splitting at exactly these two levels is what lets the interface follow
the two-phase engine (per-sample → `.psp` → cohort gather) instead of assuming the whole
panel is in memory at once. `SampleSummary` here is the ng counterpart of the existing
`.psp` `SampleSummary` section.

### Step 6 — candidate allele generation
```rust
pub trait CandidateGenerator {
    fn generate_candidates(&self, evidence: &[LocusEvidence], locus: &LocusKind,
                           params: &ModelParams) -> AlleleCandidates;
}
```
*Impls to bench:* assembly haplotypes, the repeat-length rung ladder, observed
sequences + iterative stutter/flank discovery.

### Step 7 — read likelihood (exists: `ReadLikelihoodModel`)
```rust
pub trait ReadLikelihood {
    /// P(read | allele), in log space, given the frozen error/stutter model.
    fn read_log_lik(&self, read: &LocusRead, allele: &Allele,
                    params: &ModelParams) -> LogProb;
}
```
Contamination enters here as a mixture over the source distribution (spec's choice),
not as a separate pass. The STR stutter model is the swappable part.

### Step 8 — genotype prior
```rust
pub trait GenotypePrior {
    /// log P(genotype) given the cohort frequency estimate and the sample's F.
    fn genotype_log_prior(&self, genotype: &Genotype, freq: &[AlleleFreq],
                          f: InbreedingF) -> LogProb;
}
```
*Impls to bench (the spec's ladder):* flat, fixed-het Dirichlet, SFS/Ewens
integration, marginalised leave-one-out cohort prior.

### Step 9 — posterior / inference (exists: `GenotypeEmModel`)
```rust
pub trait GenotypeModel {
    /// Run inference over a locus' per-sample likelihoods + prior, returning the
    /// per-sample calls and the site-level frequency posterior.
    fn infer(&self, lik: &[GenotypeLikelihoods], prior: &dyn GenotypePrior,
             f: &[InbreedingF]) -> LocusInference;
}
```
*Impls to bench:* ML grid (GangSTR), per-sample-GL→EM-AF (GATK/ours),
full joint cohort marginalisation (freebayes).

### Step 11 — site filtering (two traits, two questions)
```rust
pub trait ArtifactFilter {                       // 11a
    fn artifact_score(&self, site: &SiteEvidence) -> ArtifactVerdict; // drop vs annotate
}
pub trait EmissionModel {                        // 11b
    fn decide(&self, inference: &LocusInference) -> EmitDecision;     // { Emit(Phred), Drop }
}
```
*Impls to bench:* 11a inline hard-drop (ours' hidden-paralog LR) vs annotate-and-defer
(everyone else's bias stats → VQSR/CNN/dumpSTR); 11b heuristic vs BIC vs
freebayes-marginal.

### Steps 12–13 — representation and quality
```rust
pub trait AlleleRepresentation {   // step 12 — repeat-unit vs anchor-indel; routed by LocusKind
    fn to_vcf(&self, call: &GenotypeCall, cand: &AlleleCandidates, locus: &LocusKind) -> VcfAlleles;
}
pub trait QualityModel {           // step 13
    fn genotype_quality(&self, call: &GenotypeCall, inference: &LocusInference) -> Phred;
    fn site_qual(&self, inference: &LocusInference) -> Phred;
}
```

---

## 4. How swapping works — the recipe

A run selects one implementation per step. The pipeline is generic over the traits;
a `CallerRecipe` names the chosen impls, and the bake-off (spec §2) swaps exactly one
field, holds the rest, and re-measures.

```rust
pub struct CallerRecipe {
    pub read_prep:    Box<dyn ReadPrep>,
    pub router:       Box<dyn LocusRouter>,
    pub rough_caller: Box<dyn Caller>,
    pub summarizer:   Box<dyn SampleSummarizer>,   // per-sample -> .psp
    pub cohort_estimator: Box<dyn CohortEstimator>, // cohort gather -> ModelParams
    pub candidates:   Box<dyn CandidateGenerator>,
    pub likelihood:   Box<dyn ReadLikelihood>,
    pub prior:        Box<dyn GenotypePrior>,
    pub model:        Box<dyn GenotypeModel>,
    pub artifact:     Box<dyn ArtifactFilter>,
    pub emission:     Box<dyn EmissionModel>,
    pub representation: Box<dyn AlleleRepresentation>,
    pub quality:      Box<dyn QualityModel>,
}
```

`Box<dyn _>` for the lab (clarity, per-run selection, cheap to swap); the port-back
into the scaling engine (spec §3) can monomorphise the winners to generics where the
per-call vtable cost matters. The recipe *is* the experiment: "freebayes candidate
generation + HipSTR stutter likelihood + our cohort prior" is one `CallerRecipe`.

---

## 5. Design decisions (resolved during review; add new open items here with `OPEN:`)

- **The pipeline spine is a locus stream, marker-agnostic — decided.** SNP, indel, and
  STR sit at one level: step 3 segments the genome into STR / non-STR stretches, an STR
  stretch is a locus 1:1 (reference-defined), a non-STR stretch is split into loci and
  evidenced by the `pileup/` module (data-defined), and both feed one uniform stream of
  `LocusEvidence` that the per-locus core consumes identically. The router is the *only*
  principled fork; SNP and indel are both generic loci, not separate branches. `pileup/`
  is first-class infrastructure (the reused `pileup/walker/`), not a step trait — see the
  step-3 note above and `module_layout.md`. ng is **not** "STR-first" as a design; that
  phrase is only about experiment ordering (spec §3). `OPEN:` whether `pileup/` subsumes
  the generic path's step-2/window traits or is built from them.
- **Reference access is one trait, `RefSeq`, with both access patterns — decided.**
  A single `RefSeq` trait ([`../spec/ref_seq.md`](../spec/ref_seq.md)) consolidates
  production's `ChromRefFetcher` + `MultiChromRefFetcher`, with three impls (resident,
  streaming, in-memory-synthetic) reusing the `src/fasta` fetcher machinery. ng carries
  **both** the whole-contig-resident and the streaming/sub-range patterns (memory *and*
  speed, and a clean port-back), and matches production's raw-vs-canonical byte split so
  the ported filters behave identically. It is foundational infra (`ref_seq.rs`), not a
  step. **Coordinate base — decided: ng is 1-based** (matching production, VCF/SAM/IGV,
  and cutting conversion seams): `Position` is 1-based, `GenomeRegion` 1-based inclusive
  (§1), and `RefSeq` fetches by `(chrom_id, start_1based, length)`. Only BED input
  converts, at the boundary, as production already does.
- **Cohort vs per-locus granularity — decided: thread the frequency in.**
  `GenotypeModel::infer` takes a locus' worth of per-sample likelihoods *plus* the
  cohort allele frequency as an **input** (`&[AlleleFreq]`, computed by a shared outer
  loop), rather than owning the whole cohort at the locus and estimating the frequency
  itself. Owning the cohort is cheaper in the single-phase lab, but threading the
  frequency keeps the *same* trait working under the `.psp` streaming engine (which
  can't hold every sample at once) — so the winner ports back without a redesign.
  Frequency estimation becomes its own small shared step.
- **Where the parameter pre-pass sits — decided: two explicit levels (follow production).**
  Per sample (Stage 1), the rough `Caller` + `SampleSummarizer` compute the
  per-individual params (`F`) and write them into the `.psp` after the genomic blocks;
  then a **cohort-gather** step (`CohortEstimator`) reads every sample's summary and
  estimates the panel-level params (sample-groups, per-group chemistry, SFS `Theta`)
  into the final `ModelParams`, before calling. The cohort-gather step was **missing
  from the earlier step list** — it must also be added to the spec's step map (§1) and
  its "what we missed" list.
- **`Genotype` beyond diploid — decided: ploidy-agnostic type now, diploid impl.**
  `Genotype` is an **opaque multiset** of allele-ids (private field, reached via
  `alleles()` / `ploidy()` / `is_homozygous()`), so ploidy is not baked into its
  surface — diploid is just `len() == 2`. v1 does the diploid math, but a future
  polyploid impl changes only the constructor and this type, not the traits that
  consume it. (Polyploidy is common in plants, so this preparation is worth the tiny
  cost up front.)
- **Read pooling / partial observations — decided: add the weight now.** `LocusRead`
  carries a `weight: ReadWeight` (1.0 normally), which serves both a HipSTR-style
  pooling count (>= 1) and a **freebayes-style partial-observation fraction (< 1)**.
  We add it up front because partial observations are a **modeling lever we intend to
  test**, not just the pooling perf optimization; simple likelihood impls treat weight
  as 1.0 and never notice.

---

## 6. Existing-types reconciliation

The names above are the **target** clean vocabulary. Several already exist in the code
under a different name — or, worse, under the *same* name for a *different* concept — so
when ng lands these are **consolidation points**, not new types to invent alongside the
old. Verify against the code when implementing; mappings marked `≈` are approximate and
were not freshly re-read.

| ng name | existing code | action |
|---|---|---|
| `GenomeRegion` | `Region` ([regions.rs](../../../../src/regions.rs)), `ContigInterval` ([bam/alignment_input.rs](../../../../src/bam/alignment_input.rs)) | consolidate the coordinate-span types into one |
| `RefWindow` | ≈ `RefSpan` ([var_calling/types.rs](../../../../src/var_calling/types.rs)) | the sequence-carrying span |
| `RefSeq` + `RawRefSeq` (traits) | `ChromRefFetcher` + `MultiChromRefFetcher` + `RepositoryRefFetcher` + `StreamingChromRefFetcher` + `ManualEvictChromRefFetcher` ([fasta/fetcher.rs](../../../../src/fasta/fetcher.rs)) | **consolidate** into `RefSeq` (universal canonical fetch) + the `RawRefSeq` capability + an inherent `evict_before` (no silent no-ops); reuse the fetcher impls behind them. Spec: [`../spec/ref_seq.md`](../spec/ref_seq.md) |
| `MappedRead` | `MappedRead` ([bam/alignment_input.rs](../../../../src/bam/alignment_input.rs)) | reuse as-is (the step-2 input) |
| `LocusRead` | — (new) | **distinct** from `PreparedRead` ([pileup/walker/](../../../../src/pileup/walker/)), which is a *different* concept (a decoded walker read) — do **not** reuse that name |
| `AlleleCandidates` | `CandidateSet` ([ssr/cohort/candidate_set.rs](../../../../src/ssr/cohort/candidate_set.rs)) | rename |
| `SampleSummary` | ≈ the `.psp` `SampleSummary` ([sample_summary/](../../../../src/sample_summary/)) | reuse / align |
| `ModelParams` | ≈ the SSR chemistry param set ([ssr/cohort/param_estimation.rs](../../../../src/ssr/cohort/param_estimation.rs)) + per-individual `F` | assemble from both levels |
| `Genotype` | ≈ SNP/STR genotype reprs | replace with the ploidy-agnostic multiset |
| `ReadLikelihoodModel` | ([ssr/cohort/read_model/](../../../../src/ssr/cohort/read_model/)) | **exists** — build on it (step 7) |
| `GenotypeEmModel` | ([var_calling/posterior_engine.rs](../../../../src/var_calling/posterior_engine.rs)) | **exists** — build on it (step 9; the crate-level `genotype_em/` hoist is still pending) |

The lesson from this review: the codebase already has ~5 genomic-span types and ~5
read types, so the real risk when implementing ng is creating a *sixth* of each rather
than converging. This table is the convergence map.
