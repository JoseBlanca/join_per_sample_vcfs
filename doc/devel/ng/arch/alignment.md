# ng — the alignment module: types & interfaces

*Status: architecture draft (2026-07-23), companion to the spec
[`../spec/alignment.md`](../spec/alignment.md) (the design and every "why"), and to the shared arch
docs [`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary + step traits) and
[`module_layout.md`](module_layout.md) (the `src/ng/` tree). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): domain nouns for types,
verbs for functions, newtypes for domain scalars, **STR** in prose ↔ `ssr` in code. Signatures are
illustrative; the **contract** is the deliverable. See the spec for the reasoning behind every
decision recorded here.*

## Module home

`src/ng/alignment/` — a **folder**, because it holds competing implementations that get compared
([`module_layout.md`](module_layout.md) principle 1): the trait lives in `mod.rs`, each algorithm in
its own file beside its siblings.

It is **not a pipeline step**. It is called by step 2 (read preparation) and step 7 (read
likelihood), and knows about neither — it takes sequences and returns alignments or probabilities
(spec §1). That is why it sits outside the step-per-module rule rather than breaking it.

```
src/ng/alignment/
├── mod.rs                            – the three traits + the shared types below
├── emission.rs                       – the swappable per-base scoring component (§2.3)
├── stutter.rs                        – the stutter model: parameters + distribution (§2.4)
├── left_align_structured.rs          – algorithm 1a: one structured pass; 1c wraps it to a fixpoint
├── left_align_repeated.rs            – algorithm 1b: repeated simple passes, capped
├── affine_best_path.rs               – algorithm 2: banded affine matrix
├── affine_marginal.rs                – (not in the initial set, spec §5.3) algorithm 2's matrix,
│                                       summed instead of maximised
├── ssr_best_path_flat_gap.rs         – algorithm 3: two gap regimes, one flat gap in the repeat
├── ssr_best_path_unit_slip.rs        – algorithm 4: whole-unit slips priced apart
├── ssr_marginal_sequence.rs          – algorithm 5: one sequence vs another, gaps at the ends only
└── ssr_marginal_whole_read.rs        – algorithm 6: forward over the whole read
```

**The file name says which trait the algorithm implements and whether it is repeat-aware**, because
both are load-bearing: `best_path` and `marginal` implement *different traits* (§3), and every
`ssr_` algorithm has an affine counterpart it must not be confused with. `affine_` rather than
`generic_` names what the algorithm is, not which caller path uses it — this module knows no callers
(§4).

Two files carry no `ssr_` prefix although only repeat-aware algorithms use one of them: a prefix
should *disambiguate*, and there is no non-repeat stutter model to disambiguate `stutter.rs` from.
`emission.rs` is shared by both families.

Algorithm numbering is the spec's (§10.1); file names are indicative. Only algorithm 4 reads the
one `stutter.rs` model — it is written once, not per algorithm (spec §4.2). `affine_marginal.rs`
shares its matrix with `affine_best_path.rs` through a private helper rather than reimplementing it
(spec §10.2).

## 1. The domain scalar this module seeds

`LogProb` is **already specified but not yet written**: `ng_step_interfaces.md` §1 declares it in the
shared vocabulary, and `module_layout.md` already assigns it to `src/ng/types.rs` — but `grep` over
`src/` finds nothing. So this module does not introduce it; it is its **first user**, and lands it
where the vocabulary already put it.

```rust
/// A probability held as its natural logarithm. `NEG_INFINITY` means impossible — a real value,
/// not an error. Unconstrained (any finite f64 or −∞ is a valid log-probability), so the field
/// is public with a `.get()`, per the types.rs convention.
pub struct LogProb(pub f64);
```

Its whole purpose is that the compiler refuses to mix it with an ordinary probability — the
transposition this module is most exposed to, since production's existing code returns linear
probabilities (§5) and the boundary conversion is one call site.

## 2. Types

### 2.1 What a best-path aligner returns

Two output shapes, one per marker kind — which is why `BestPathAligner::Output` is an associated
type rather than a fixed one (§3).

The affine aligner returns a placement plus the operations that got there. **It reuses production's
`CigarOp`** rather than minting a parallel operation type (§5):

```rust
/// A read placed against a reference stretch: where it starts, and the operations from there.
pub struct Alignment {
    /// Offset into the reference stretch the aligner was given — not a genome position.
    pub reference_offset: usize,
    pub ops: Vec<CigarOp>,
}
```

The repeat-aware aligners return where the read's repeat sits, and **which flanks held it** — the
part production's `Delimited` does not carry, and the part the STR preparer needs to tell a measured
repeat from a lower bound (spec §8; `../spec/read_preparation_ssr.md` §3):

```rust
/// Where a read's repeat lies, in read coordinates, and what anchored it.
/// A span is only a *measurement* when both flanks anchored; otherwise it is a lower bound.
pub enum RepeatSpan {
    /// Both flanks anchored: the span pins the repeat's length.
    Between(Range<usize>),
    /// Only the left flank anchored — the read ran off its own end inside the repeat.
    FromLeft(Range<usize>),
    /// Only the right flank anchored.
    FromRight(Range<usize>),
    /// Neither flank anchored: the read lies wholly inside the repeat. No per-read fact.
    Unanchored,
}
```

### 2.2 The per-call context — what a repeat-aware algorithm is told each time

Every repeat-aware algorithm needs to know where the repeat is and how it slips. Both change at every
locus, so they are a **per-call argument, not constructor state**: holding them would force a new
aligner per locus across millions of loci, and `Motif` is an allocation. This mirrors production,
whose read-likelihood model is a stateless value taking a per-call `ReadScoringContext` (§5).

```rust
/// What varies from one repeat to the next. `BestPathAligner::Context` and
/// `MarginalAligner::Context` for the repeat-aware algorithms; `()` for the affine ones,
/// which need nothing locus-specific.
pub struct RepeatContext<'a> {
    pub geometry: &'a RepeatGeometry,
    /// Only the algorithms that price whole-unit slips read this (algorithm 4; §2.4).
    pub stutter: &'a StutterModel,
}

/// Where the repeat sits inside the reference stretch the aligner is given, measured in bases
/// from its start. Both flanks may be short or absent when a repeat is near a contig end, so
/// these are measured, never assumed equal.
pub struct RepeatGeometry {
    pub left_flank_len: usize,
    pub right_flank_len: usize,
    /// The repeat unit. Needed only by the algorithms that price whole-unit slips.
    pub motif: Motif,
}
```

### 2.3 Emission — the swappable per-base scoring

Per-base quality or one flat error rate. A **component**, so the two are configurations of one
algorithm rather than two algorithms (spec §3). A type parameter, not a trait object: this is called
once per matrix cell.

```rust
/// Scores one read base against one reference base, in log space.
pub trait Emission {
    /// `quality` is ignored by implementations that do not use it.
    fn emit_ln(&self, read_base: u8, reference_base: u8, quality: BaseQual) -> f64;
    /// An inserted read base, which has no reference base to score against.
    /// NOTE the two implementations do not have a common source for this. Production's
    /// per-quality path scores it against a uniform base composition, ln(1/4); its flat path
    /// has no emission for it at all — the value it uses there is a *transition* cost, not an
    /// emission. The flat implementation's value is therefore a decision, not a port.
    fn insert_ln(&self) -> f64;
}
```

Two implementations: one reading a per-quality table (production's, §5), one a flat rate.
**Contract:** pure and total — no per-base quality may yield a non-finite score, including quality
zero, which production floors precisely so it cannot annihilate a path (spec §4.2).

### 2.4 The stutter model

The seven parameters and the distribution over length changes, spelled out in spec §5.2. Used by
algorithm 4 inside this module, and by the genotyping likelihood outside it — algorithms 5 and 6 do
not use it (they are pure sequence comparison).

```rust
/// How likely each length change is, for one repeat. Two regimes: whole-unit slips (in frame)
/// and everything else (out of frame), each split by direction. Built **per locus** — from the
/// locus's stutter shape and its REFERENCE allele length, never from a candidate's length (§5:
/// the per-candidate slip level is the genotyping's, not this module's). See spec §5.2 for the
/// distribution and its two conversion traps.
pub struct StutterModel {
    pub equal: f64,
    pub in_up: f64,
    pub in_down: f64,
    /// Geometric *success* probability per unit — NOT a continuation probability (spec §5.2).
    pub in_geom: f64,
    pub out_up: f64,
    pub out_down: f64,
    pub out_geom: f64,
}

impl StutterModel {
    /// `P(length change)` for a change of `bp_diff` bases on a repeat of period `period`.
    /// Zero beyond the slip cutoff, so an implausibly large slip is not explained away.
    pub fn probability(&self, bp_diff: i32, period: u8) -> f64;
}
```

**Contract:** `equal` is *defined* as one minus the four direction masses, **but it is floored**, so
when the floor binds the five values sum to slightly more than one. That is deliberate — the floor is
what stops a hostile parameter combination producing a negative probability — and it means "the five
sum to one" must **not** be written as a test. (HipSTR instead *asserts* the masses sum below one at
construction and never clamps; both disciplines are defensible, but they differ, and only the floor
matches the code being ported.) The geometric probabilities are likewise held strictly inside (0, 1)
(spec §5.2).

## 3. The interfaces

Three traits, one per shape of transformation — what goes in and what comes out, with nothing about
who calls them (spec §7).

```rust
/// Line a read up against a reference sequence the single most probable way.
pub trait BestPathAligner {
    /// Reused matrices and traceback buffer — allocated per worker, never per read.
    type Scratch: Default;
    /// `Alignment` for the affine aligner; `RepeatSpan` for the repeat-aware ones.
    type Output;
    /// `()` for the affine aligner; `RepeatContext<'_>` for the repeat-aware ones (§2.2).
    type Context;

    fn align(&self, read: &[u8], quality: &[u8], reference: &[u8],
             context: &Self::Context, scratch: &mut Self::Scratch) -> Self::Output;
}

/// The total probability the reference produced this read, summed over every line-up.
pub trait MarginalAligner {
    type Scratch: Default;
    type Context;

    fn marginal_probability(&self, read: &[u8], reference: &[u8],
                            context: &Self::Context, scratch: &mut Self::Scratch) -> LogProb;
}

/// Rewrite an existing alignment into its canonical (left-most) spelling, in place.
pub trait AlignmentNormalizer {
    /// Takes the whole `Alignment`, not just its ops: left-alignment can shift a leading
    /// deletion off the front, which moves `reference_offset`. Production's does exactly
    /// that, so a signature over the ops alone would silently foreclose it.
    fn normalize(&self, alignment: &mut Alignment, read: &[u8], reference: &[u8]);
}
```

**Contract, common to all three.** Pure per call: the output is a function of the arguments and the
value's constructor state, with no hidden mutation beyond `Scratch` — so results never depend on
call order or thread count, which the cohort's byte-identity guarantee rests on. Scratch is
caller-owned and reused; an implementation that allocates per call is a defect, not a slow path.
Every implementation is selected by generic type parameter, never behind `dyn` (§4).

**Errors: none — but the caller owes invariants.** No operation here returns a `Result`, deliberately:
a fallible alignment would push error handling onto the caller's hottest path for cases that are
answers, not failures. Degenerate *inputs* have defined results — an empty reference gives
`RepeatSpan::Unanchored`, an unreachable line-up gives `LogProb(NEG_INFINITY)`.

That is not the same as "nothing can go wrong", and the difference must be stated rather than
implied. Two conditions are the **caller's** to uphold, because nothing here can check them cheaply:

- the repeat geometry must fit the reference it is used with (`left_flank_len + right_flank_len` no
  greater than the reference length) — production computes a boundary by subtraction and would
  underflow if this were violated, and is safe only because its locus type enforces it upstream;
- the read and its quality slice must be the same length.

State them as **debug assertions plus documented preconditions**, not as error variants — and note
that a debug assertion compiles out of the release build this project actually runs, so a violated
precondition is a wrong answer rather than a crash. If either turns out to be reachable from
untrusted input, it becomes a checked constructor on the context type, not a `Result` on the hot call.

**Determinism.** The best-path aligners must break ties by a fixed rule and assign a gap on a
flank/repeat boundary to the block on its 5′ side. **Two reads of the same molecule must measure the
same repeat**; without both rules, equally-scoring line-ups let them measure
different repeats from identical input (spec §4.2).

## 4. Design decisions — decided

The spec argued each of these; the record carries the code shape and a cross-ref, not the argument.

- **A folder with the traits in `mod.rs`, one algorithm per file — decided.** Eight entries get
  compared — three normalizers and five aligners; that is the bake-off shape
  (`module_layout.md` principle 1).
- **Not a pipeline step, and it knows none — decided.** Two steps call it. Its interfaces are keyed
  on the shape of the transformation, never on the caller's purpose (spec §1, §7).
- **Three interfaces, split by what goes in and out — decided.** Normalization rewrites an
  alignment; the two aligners take sequences and differ in what they return. The two aligners share
  a noun because they are one recurrence under two reductions (spec §7).
- **The marginal returns `LogProb` — decided.** Callers multiply across reads; an ordinary
  probability underflows. This constrains the *returned value* only — an implementation may run its
  matrix in linear space and take one logarithm at the end (spec §7).
- **Static dispatch only — decided.** Per-read, cohort-wide, so a virtual call per read is not
  affordable. Implementations are generic parameters (spec §7).
- **Emission is a component, not a variant — decided.** With-quality and without-quality are two
  configurations of one algorithm, which is what makes the comparison a swap (spec §3).
- **The locus travels per call; only the emission model is constructor state — decided.** Geometry
  and stutter change at every repeat, so holding them would mean a new aligner per locus across
  millions of loci, with an allocated motif inside. The emission model stays fixed because it is the
  experiment's configuration and carries the sample group's error rate. This is production's own
  shape — a stateless model plus a per-call context (spec §7). *(An earlier draft made the locus
  constructor state, to keep the call to two sequences; that argument does not survive the
  per-locus cost.)*
- **The stutter model is written once and shared — decided.** Algorithm 4 uses the same
  parameters and distribution as the genotyping likelihood that composes this module; two copies would
  drift (spec §4.2, §5.2).
- **Per-base alignment confidence (BAQ) is not provided — decided.** It is a third output shape and
  would need a third interface; production experience says it is not worth carrying (spec §1). If
  reopened, it is the same machinery a repeat-boundary posterior would need.
- **Wavefront, bit-vector and difference-recurrence cores are closed, not deferred — decided.** Each
  computes a best path only, and none can consume the per-base qualities the affine aligner uses.
  **This closes those three cores, not vectorisation as a technique** — a vectorised *forward* is a
  different kernel and remains open (spec §4.1).
- **Normalization ships three algorithms and is graded against a definition — decided.** Neither
  production's structured pass nor freebayes' capped iteration is provably leftmost, so there is a
  real fork; and because "leftmost" is a definition, a property test grades all three against ground
  truth rather than against each other. *(An earlier draft recorded this as one algorithm with
  nothing to compare, on an inaccurate description of both incumbents — spec §6.)*
- **`RepeatSpan` carries which flank anchored — decided.** Production's `Delimited` does not, and
  the STR preparer cannot otherwise tell a measurement from a lower bound (spec §8).
- **Algorithm 3 ships unbanded; banding is a later, separately-proved change — decided.**
  Production's delimiter fills its whole matrix, so only an unbanded port can be shown byte-identical
  to it — and that parity is this module's only hard oracle. Banding then lands on its own, at the
  width spec §3 sets: a **per-read floor of |read − reference| length difference**, plus headroom.
  Not a compile-time constant, and **not** derived from the slip cutoff — that would be numerically
  too narrow for production's own long-allele test *and* would let a scoring parameter blind the
  ruler (spec §3, §10.3).

## 5. Reconciliation with existing code

Every row read and cited. New code is marked; everything else converges on an existing type or
function rather than duplicating it.

| ng name | existing code | action |
|---|---|---|
| `LogProb` | declared in [`ng_step_interfaces.md`](ng_step_interfaces.md) §1; assigned to `types.rs` by [`module_layout.md`](module_layout.md) | **specified, not yet written** — `grep` over `src/` finds none. This module is its **first user**, not its author (§1); land it where the vocabulary already put it, beside `Bp` ([types.rs:151](../../../../src/ng/types.rs#L151)) |
| `Alignment.ops`, `AlignmentNormalizer::normalize` | `left_align_indels(cigar: &mut Vec<CigarOp>, seq, ref_seq)` ([indel_norm.rs:408](../../../../src/pileup/walker/indel_norm.rs#L408)) | reuse `CigarOp`. The trait takes the whole `Alignment` rather than this function's `Vec<CigarOp>` **deliberately** — `left_align_cigar` ([:254](../../../../src/pileup/walker/indel_norm.rs#L254)) can drop a leading deletion and move the alignment's start, which an ops-only signature could not express (§3) |
| algorithm 1a (structured pass) | `left_align_indels` ([indel_norm.rs:408](../../../../src/pileup/walker/indel_norm.rs#L408)), `left_align_cigar` ([:254](../../../../src/pileup/walker/indel_norm.rs#L254)), `normalize_alleles` ([norm_seqs.rs:108](../../../../src/norm_seqs.rs#L108)) | port. **Not a naive single pass** — it merges consecutive indels and propagates one across alignment blocks; port those, they are the point |
| algorithm 1b (repeated passes) | freebayes `stablyLeftAlign` (vendored, `freebayes/src/LeftAlign.cpp:385`, cap 20 at `LeftAlign.h:118`) | port the bounded loop. Note freebayes **returns false on exhaustion and its caller ignores it** — do not copy that; 1c exists precisely to fail loudly instead |
| algorithm 1c (fixpoint) | — | **new**; a loop around 1a, not a third shift implementation |
| algorithm 3 (two gap regimes) | `delimit_read` ([alignment.rs:171](../../../../src/ssr/pileup/alignment.rs#L171)), `HmmModel` ([:80](../../../../src/ssr/pileup/alignment.rs#L80)), `ViterbiScratch` ([:111](../../../../src/ssr/pileup/alignment.rs#L111)) | port; **the parity anchor**. **Unbanded** — the inner loop runs the full row and `ViterbiScratch` allocates all `(m+1)×(n+1)` backpointer cells. Port it that way, so parity is provable; band later (§4) |
| `RepeatSpan` | `Delimited` ([alignment.rs:141](../../../../src/ssr/pileup/alignment.rs#L141)) | port and **widen** — `Delimited::BorderOffEnd` does not say which flank held |
| gap constants | `GAP_OPEN_PROB` 2.9e-5 ([:38](../../../../src/ssr/pileup/alignment.rs#L38)), `GAP_OPEN_PROB_TRACT` 1e-2 ([:50](../../../../src/ssr/pileup/alignment.rs#L50)) | import as named `pub const`s; **provisional calibration**, not findings (spec §4.2) |
| `Emission` (per-quality impl) | `EMISSION_LN` ([alignment.rs:59](../../../../src/ssr/pileup/alignment.rs#L59)) | port the table, including its quality-zero floor |
| `Emission` (flat impl) | the flat `eps` arm of `align_subst` ([pair_hmm.rs:52](../../../../src/ssr/cohort/pair_hmm.rs#L52)) | port |
| algorithm 5 (repeat marginal) | `align_subst` ([pair_hmm.rs:52](../../../../src/ssr/cohort/pair_hmm.rs#L52)) and its `banded_forward` ([:83](../../../../src/ssr/cohort/pair_hmm.rs#L83)) | port **this function only** — not `HipstrModel`, which is the genotyping model that calls it (row below). **Returns a linear probability; this interface returns its logarithm** — convert at the boundary. Its banding differs from algorithm 3's source: the *computation* is band-limited but the full matrix is still allocated — the two ports do not share a banding stance. **Interior gaps are forbidden** (`FLANK_SLOP` = 2, [:28](../../../../src/ssr/cohort/pair_hmm.rs#L28)) and that restriction is load-bearing, not an optimisation (spec §5.1) |
| **NOT this module:** the stutter half | `HipstrModel::q_r` ([hipstr.rs:85](../../../../src/ssr/cohort/read_model/hipstr.rs#L85)) | the genotyping likelihood (step 7). It resizes the candidate to the observed length, so `align_subst` **always** takes its equal-length early return and `banded_forward` is unreachable from here — a parity test built around this model never exercises the forward pass (spec §5.1) |
| `StutterModel` | `HipstrParams` + `stutter_pmf` ([hipstr.rs:168](../../../../src/ssr/cohort/read_model/hipstr.rs#L168)) | port; hoist from private struct to this module's shared type |
| `StutterModel` construction | `hipstr_params` ([hipstr.rs:131](../../../../src/ssr/cohort/read_model/hipstr.rs#L131)) from `StutterShape` ([param_estimation.rs:46](../../../../src/ssr/cohort/param_estimation.rs#L46)) | **adapt, do not port straight.** `hipstr_params` reads a *per-candidate* slip level, because stutter rises with allele length — that dimension belongs to the genotyping, which scores against candidates. This module measures, so its model is built from the locus and the **reference** allele length. `OUT_FRAME_REL` = 0.05 ([hipstr.rs:44](../../../../src/ssr/cohort/read_model/hipstr.rs#L44)) is a **placeholder, not an estimate** (spec §5.2) |
| the per-call context | `ReadScoringContext` ([read_model/mod.rs:45](../../../../src/ssr/cohort/read_model/mod.rs#L45)) | the precedent for `RepeatContext` (§2.2) — production already passes the varying facts per call and keeps the model itself a stateless value |
| slip cutoff | `MAX_SLIP` = 10 ([param_estimation.rs:21](../../../../src/ssr/cohort/param_estimation.rs#L21)) | import as a named const |
| placement sum | `reach_variants` ([stutter.rs:140](../../../../src/ssr/cohort/stutter.rs#L140)), `PlacementVariant` ([stutter.rs:21](../../../../src/ssr/cohort/stutter.rs#L21)) | port — the interrupted-repeat placement enumeration (spec §5.2) |
| scratch shape | `HmmScratch` ([pair_hmm.rs:32](../../../../src/ssr/cohort/pair_hmm.rs#L32)), `ViterbiScratch` ([alignment.rs:111](../../../../src/ssr/pileup/alignment.rs#L111)) | model for each algorithm's `Scratch`; per-algorithm, not one shared buffer |
| consumer (step 7) | `ReadLikelihoodModel` ([read_model/mod.rs:63](../../../../src/ssr/cohort/read_model/mod.rs#L63)), `ReadScoringContext` ([:45](../../../../src/ssr/cohort/read_model/mod.rs#L45)) | **not** reimplemented here — step 7's trait composes a `MarginalAligner` |
| consumer (step 2) | `MappedRead` ([alignment_input.rs:78](../../../../src/bam/alignment_input.rs#L78)) | never seen by this module — the preparer unpacks it to slices (spec §7) |
| algorithms 2, 4, 6 | — | **new**; no production counterpart |

Production is a **read-only oracle**: ng ports from it and never edits it.

## 6. Open items

- **`OPEN:` the affine aligner's `Output`.** `Alignment` above assumes a re-placement is expressible
  as an offset plus `CigarOp`s. If the re-align mode needs to carry more than the mapper's record
  does, this widens — and `../spec/read_preparation_generic.md` §8 already names that as its fork
  trigger. Blocked behind the same decision that gates algorithm 2 (spec §10.3).
- **`OPEN:` how much headroom above the band floor.** The floor is the per-read length difference
  and needs no decision; the headroom above it is a real open parameter (spec §9). Production's
  forward uses `FLANK_SLOP` = 2 ([pair_hmm.rs:28](../../../../src/ssr/cohort/pair_hmm.rs#L28)); whether
  that suits the delimiter is untested. It lands as a named `pub const` with its reasoning, never a
  literal — a too-small value loses long alleles silently. **Do not reach for
  [`MAX_SLIP`](../../../../src/ssr/cohort/param_estimation.rs#L21)**: it is a scoring cutoff, and an
  earlier draft's band derived from it was both too narrow and a layering violation (spec §3).
- *Impl-time confirmation, not a decision:* the exact `Emission` signature (whether quality arrives
  per call or as a pre-resolved row) resolves when the first two implementations exist.

*(An earlier draft filed "is `StutterModel` per locus or per read?" here as an impl-time
confirmation. It was neither — it decides the trait signature, and it is now decided: the locus
travels per call, §2.2 and §4.)*

The spec's remaining open questions (§9) are design questions settled by measurement, not interface
questions — they do not change the shapes above.

## Test & bench shape

- **Unit tests beside each algorithm.** For the repeat-aware best-path aligners: a clean repeat
  measures exactly; a longer allele does not collapse to the reference; an interrupted repeat comes
  out verbatim; each of `RepeatSpan`'s four cases is reachable. For `StutterModel`: the distribution
  reproduces the published formula term by term, and both conversion traps are pinned by a test that
  fails if the geometric is inverted (spec §5.2).
- **The regression anchor is algorithm 3 against production.** Its measured repeat must match
  `delimit_read` ([alignment.rs:171](../../../../src/ssr/pileup/alignment.rs#L171)) byte for byte on
  a shared fixture. It is the only *aligner* with a parity oracle: algorithms 2, 4, 5 and 6 are
  measured, not verified. (The three normalizers have their own, better oracle — the property test
  below.) The ported version is **unbanded**, matching production (§4) — and the later banding
  commit carries its own test: the same fixture, the same output, which is what makes banding a
  performance change rather than a behaviour change.
- **A property test for the normalizers, which is stronger than their bake-off.** For every indel in
  an output, assert it cannot shift one base left and still represent the same read against the same
  reference. This grades 1a, 1b and 1c against the *definition* of leftmost rather than against each
  other, and it is the only place in this module where ground truth is available without a rival
  implementation (spec §6). Pair it with the differ-at-all screen: run all three over the same reads
  and count disagreements before investing in anything larger.
- **`bench/` — yes, unlike read filtering.** This module *is* a bake-off, so it carries the
  comparison harness: synthetic reads with known truth, scoring **calibration as well as accuracy**,
  and running the algorithm-3-vs-4 comparison at **period 1 as well as period 2 and above** — the two
  regimes exercise different parts of the model, and period 1 is where the indel deficit lives. At
  period 1 the fixtures must separate an indel of the repeat's own base from an indel of a different
  base, because only the first is genuine stutter and averaging the two cancels the effect
  (spec §4.2, §10.3).
