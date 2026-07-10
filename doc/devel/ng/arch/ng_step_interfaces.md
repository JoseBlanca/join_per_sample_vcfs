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
pub struct Position(pub u64);       // 0-based reference position
pub struct RefInterval { pub contig: ContigId, pub start: Position, pub end: Position } // half-open
pub struct SsrMotif(Box<[u8]>);     // the repeat unit, phase-faithful (STR)
pub struct SsrPeriod(pub u8);       // motif length in bp, the repeat period (STR)
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
pub struct SampleId(pub u32);       pub struct SampleGroupId(pub u32);
pub struct LocusId(pub u64);        pub struct AlleleId(pub u16);   // index within a locus' candidate set
```

### Probability & quality — keep log-space explicit

```rust
pub struct LogProb(pub f64);        // natural-log probability; the internal currency
pub struct Phred(pub f32);          // -10 log10 p, for I/O (QUAL, GQ)
pub struct BaseQual(pub u8);        pub struct MapQual(pub u8);
```

`LogProb` vs `Phred` are different scales with different bases; making them distinct
types stops the classic "added a Phred to a ln-prob" error. Conversions are named
functions (`Phred::from_log_prob`, `LogProb::ln`), never `as`.

### Population genetics

```rust
pub struct AlleleFreq(pub f64);     // in [0,1]
pub struct InbreedingF(pub f64);    // Wright's F_IS, in [0,1)
pub struct Theta(pub f64);          // Watterson/Ewens diversity for the SFS prior
```

> **Ergonomics.** Wrappers cost a `.get()` at arithmetic sites; for the few types with
> heavy arithmetic (`LogProb`, `BpDelta`) implement the needed operators so hot code
> reads naturally. Do **not** blanket-impl `Deref` to the primitive — that would
> re-open the transposition hole the newtype exists to close.

---

## 2. The data that flows between steps (wire types)

These are the payloads passed step→step; every step impl consumes and produces them.

```rust
/// One sample's evidence at one locus — the per-locus, per-sample unit.
pub struct LocusEvidence<'a> {
    pub sample: SampleId,
    pub reads: &'a [PreparedRead],   // output of step 2
    pub depth: ReadDepth,
}

/// A read after step 2 (realigned / delimited). The STR path stores the observed
/// tract bytes; the generic path stores the aligned span.
pub struct PreparedRead {
    pub map_qual: MapQual,
    pub observed: Observation,       // enum: generic aligned allele, or an STR tract
}

/// The router's verdict for a locus (spec step 3). This is what makes STR-awareness
/// a *type*, not a convention: only `Ssr` carries a motif/period/borders.
pub enum LocusKind {
    Generic(RefInterval),
    Ssr(SsrLocus),
}
pub struct SsrLocus {
    pub id: LocusId,
    pub span: RefInterval,           // the tract, without flanks
    pub motif: SsrMotif,
    pub period: SsrPeriod,
    pub ref_units: SsrRepeatUnits,
}

/// The candidate allele set at a locus (step 6). REF is always present at `ref_idx`.
pub struct CandidateSet {
    pub alleles: Vec<Allele>,        // sequence-resolved
    pub ref_idx: AlleleId,
    pub admit: Admission,            // enum { Ok, LowDepth, NotPeriodic, TooManyAlleles }
}

/// Per-sample genotype log-likelihoods at a locus (step 7 → step 9). The single
/// read-level quantity every prior/posterior/emission model consumes (spec's `Lg`).
pub struct GenotypeLikelihoods {
    pub genotypes: Vec<Genotype>,    // the enumerated diploid genotypes over the candidates
    pub log_lik: Vec<LogProb>,       // parallel to `genotypes`, per sample handled by caller
}

/// A called genotype for one sample at one locus.
pub struct GenotypeCall {
    pub genotype: Genotype,          // allele-id pair (diploid v1)
    pub quality: Phred,              // GQ
}

/// The frozen output of the parameter pre-pass (spec step 4) — the "first caller's"
/// legacy, consumed by steps 7/8. Every field is an *input*, never re-fit downstream.
pub struct FrozenParams {
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

### Step 2 — read preparation / realignment
```rust
pub trait ReadPrep {
    /// Realign/delimit one read against the reference window (STR impls also take the
    /// tract so gaps can be tract-aware). Returns None if the read is unusable here.
    fn prepare(&self, read: &MappedRead, window: &RefWindow) -> Option<PreparedRead>;
}
```
*Impls to bench:* trust-mapper+left-align (freebayes-style), local reassembly
(GATK-style), per-read pair-HMM to the tract (ours-STR / HipSTR-style).

### Step 3 — the locus router
```rust
pub trait LocusRouter {
    /// Classify a reference region: is it an STR (with motif/period/borders) or generic?
    fn route(&self, window: &RefWindow) -> LocusKind;
}
pub trait LocusSource {
    /// Enumerate the candidate loci to genotype — a catalog iterator, or active regions.
    fn loci(&self) -> impl Iterator<Item = RefInterval>;
}
```
*Impls to bench:* fixed STR catalog (ours), active-region detection (GATK), the
catalog+discovery hybrid (spec §step-3). Keeping `LocusSource` and `LocusRouter`
separate lets a data-driven source feed the same reference-based router.

### Step 4 — the rough caller and the parameter estimator (the "two callers")
```rust
pub trait RoughCaller {
    /// A cheap first-pass genotype that must NOT depend on the parameters being
    /// estimated (spec: it uses a crude base-quality ε̂, not the fitted ε). Returns
    /// only *confident* calls — the subset allowed to teach the parameters.
    fn call_confident(&self, evidence: &LocusEvidence) -> Option<ConfidentGenotype>;
}
pub trait ParameterEstimator {
    /// Freeze the nuisance parameters from the rough caller's confident genotypes.
    fn estimate(&self, confident: &[ConfidentGenotype]) -> FrozenParams;
}
```
*Contract:* `RoughCaller` is the confident-subset gate (our het-margin / STR
1-vs-2 BIC test); `estimate` never sees a non-confident call. `F` is per-sample,
chemistry is pooled per group (spec §4 considerations).

### Step 6 — candidate allele generation
```rust
pub trait CandidateGenerator {
    fn candidates(&self, evidence: &[LocusEvidence], locus: &LocusKind,
                  params: &FrozenParams) -> CandidateSet;
}
```
*Impls to bench:* assembly haplotypes, the repeat-length rung ladder, observed
sequences + iterative stutter/flank discovery.

### Step 7 — read likelihood (exists: `ReadLikelihoodModel`)
```rust
pub trait ReadLikelihood {
    /// P(read | allele), in log space, given the frozen error/stutter model.
    fn read_log_lik(&self, read: &PreparedRead, allele: &Allele,
                    params: &FrozenParams) -> LogProb;
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
    fn to_vcf(&self, call: &GenotypeCall, cand: &CandidateSet, locus: &LocusKind) -> VcfAlleles;
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
    pub rough_caller: Box<dyn RoughCaller>,
    pub estimator:    Box<dyn ParameterEstimator>,
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

## 5. Open questions for the interfaces

- **Cohort vs per-locus granularity.** `GenotypeModel::infer` takes a locus' worth of
  per-sample likelihoods; the cohort prior needs cross-sample frequency. Decide whether
  the frequency estimate is threaded in (`&[AlleleFreq]` computed by an outer EM) or
  the model owns the whole cohort at a locus. The single-phase ng makes the latter
  cheap; the port-back to `.psp` streaming may force the former.
- **Where the parameter pre-pass sits.** `RoughCaller` + `ParameterEstimator` produce
  `FrozenParams` once per run; in the two-phase engine the rough caller is split
  (SNP per-sample, STR cohort burn-in). The interface should not bake in *where* it
  runs, only that its output is frozen before steps 7–9.
- **`Genotype` beyond diploid.** The types say diploid (allele-id pair). Ploidy is a
  documented follow-up; `Genotype` should be an opaque multiset behind an accessor so
  a polyploid impl doesn't ripple through every trait.
- **Read pooling / partial observations** (HipSTR / freebayes) change `LocusEvidence`
  from `&[PreparedRead]` to weighted/unique reads. Decide if that is a `PreparedRead`
  weight field now or a later refinement.
