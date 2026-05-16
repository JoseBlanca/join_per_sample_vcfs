# Posterior engine (Stage 6 — EM over allele frequency, HWE-with-F prior, optional contamination)

Proposal date: 2026-05-16.

> **Status: draft.** First-pass scaffolding for review. Open design
> questions are flagged inline (search for `**Open question**`); the
> "Decisions to confirm before implementation" section at the bottom
> collects them in one place. Iterate before locking the plan.

## Domain intent

Stage 6 is the last cohort-side stage. It consumes the
`MergedRecord` stream produced by Stage 5 and emits, per record, a
per-sample, per-genotype posterior table plus the site-level
quantities (`p̂`, `QUAL`) that the VCF writer needs:

```
… ─► PerGroupMerger ─► MergedRecord stream ─► PosteriorEngine ─► PosteriorRecord stream ─► VcfWriter
```

Stage 6's contract from
[calling_pipeline_architecture.md §"Stage 6 — posterior engine"](../specs/calling_pipeline_architecture.md#L1815-L1908):

1. **Run EM** over the per-record allele frequency `p`, the
   per-compound cohort frequency `f_C` (one per chain-anchored
   compound in the record), and — when contamination correction is
   enabled — the per-batch contamination source `q_b` and per-sample
   contamination fraction `c_s`.
2. **At convergence**, emit per-sample per-genotype posteriors
   under a HWE-with-`F` prior, plus the site-level `QUAL` (Phred of
   `Π_s P(hom-ref)_s`) and the per-sample `GT` (argmax) + `GQ`
   (Phred of `1 − P(best)`).
3. **Stage 6 owns the math, not the VCF.** The VCF writer is a
   separate module downstream; Stage 6 emits a structured
   `PosteriorRecord` carrying everything the writer needs.

This is the smallest of the cohort-side stages in code volume
(closed-form M-steps, no allele unification, no scalar projection)
but the most subtle in EM bookkeeping — particularly when
contamination is enabled, because `c_s` and `q_b` are estimated
**across records** while `p̂` and `f_C` are estimated **within each
record**. See §"Per-record vs cross-record M-steps" for the
v1 architecture decision and the alternatives considered.

## Why now

Stage 5 has shipped (fixes-applied 2026-05-16,
[src/var_calling/per_group_merger.rs](../../src/var_calling/per_group_merger.rs)),
so its output type — `MergedRecord` carrying merged alleles, scalar
table, chain-anchor flags, and per-sample per-genotype
log-likelihoods under `c_s = 0` — is fixed and the test fixtures
that exercise it are proven. Stage 6 is the final hop required to
emit a multi-sample VCF; once it lands, the new pipeline reaches
end-to-end functional parity with the legacy gVCF path and the
deferred gVCF cleanup ([cohort_per_group_merger.md §"Deferred
cleanup"](cohort_per_group_merger.md#L1401-L1426)) can proceed.

A perf review of Stage 5 is in flight in parallel; it does not
gate this plan — Stage 6's correctness contract is independent of
Stage 5's internal layout.

## What's already in place

- **Stage 5's `MergedRecord`** —
  [src/var_calling/per_group_merger.rs:133-163](../../src/var_calling/per_group_merger.rs#L133-L163).
  Every field Stage 6 needs is there:
  - `alleles: Vec<MergedAllele>` — REF at index 0; `is_compound`
    bit per allele.
  - `scalars[sample_idx][allele_idx]: AlleleSupportStats` — the
    five freebayes scalars per (sample, allele). Required for the
    contamination-augmented E-step in mixture mode; unused in the
    no-contamination path (Stage 5 has already precomputed the
    `c_s = 0` log-likelihoods).
  - `other_scalars[sample_idx]: AlleleSupportStats` — pooled
    scalars for alleles dropped by Stage 5's `max_alleles` cap;
    needed by the contamination E-step's error-cost recomputation.
  - `chain_anchor_flags[sample_idx][allele_idx]` — `true` when
    the sample is chain-broken at a compound allele; flips the
    sample to the cohort-prior path for `f_C` (see Stage 5
    §"Chain-broken compound fallback").
  - `log_likelihoods[sample_idx][genotype_idx]` — natural-log
    likelihood under `c_s = 0`. The starting input to the EM in
    no-contamination mode; the iteration-0 input in
    contamination mode (subsequent E-steps recompute from
    `scalars`).
  - `ploidy: u8`, `chrom_id`, `start`, `end` — metadata.
- **Genotype enumeration** —
  [src/var_calling/per_group_merger.rs:299-321](../../src/var_calling/per_group_merger.rs#L299-L321).
  `genotype_order(ploidy, n_alleles)` returns the canonical VCF/GATK
  ordering shared by all samples in a record. Stage 6 indexes
  posteriors by the same order.
- **Architecture spec** — the 9-step EM algorithm at
  [calling_pipeline_architecture.md:1837-1894](../specs/calling_pipeline_architecture.md#L1837-L1894).
  The prior shape (Dirichlet pseudocounts) and HWE-with-`F`
  formulation are pinned at
  [calling_pipeline_architecture.md:1581-1729](../specs/calling_pipeline_architecture.md#L1581-L1729).
- **Background reading.** The freebayes posterior derivation
  ([freebayes_posterior_gt_probs.md](../specs/freebayes_posterior_gt_probs.md))
  and the GATK GenotypeGVCFs algorithm
  ([gatk_em_calculation.md](../specs/gatk_em_calculation.md))
  cover the two reference algorithms. The pipeline's choice is
  Dirichlet + EM (GATK-shaped) with a per-record scope and a
  HWE-with-`F` prior that GATK does not have. The rationale is
  in [calling_pipeline_architecture.md:1615-1663](../specs/calling_pipeline_architecture.md#L1615-L1663)
  ("Why Dirichlet rather than Ewens").

## What is *not* in place (and is deliberately not being carried forward)

- **The legacy gVCF posterior module** —
  [src/genotype_posteriors.rs](../../src/genotype_posteriors.rs).
  Obsolete code from the pre-pivot gVCF path; slated for deletion
  as part of the deferred gVCF cleanup
  ([cohort_per_group_merger.md §"Deferred cleanup"](cohort_per_group_merger.md#L1401-L1426)).
  Stage 6 does **not** extend, refactor, or copy from it. The
  PriorConfig / EM / QUAL machinery there used `log10` and consumed
  gVCF-style PLs; the new engine works in natural-log on
  `MergedRecord` natural-log likelihoods. The legacy file is
  referenced here only so a reader looking for "doesn't this
  already exist?" knows the answer is "yes, and it's being deleted,
  not generalised."

## Algorithmic alternatives considered

### Algorithm 1 — per-record independent EM (chosen for v1, no-contamination mode)

**Shape.** Each `MergedRecord` runs its own EM loop, independent of
all others. The EM converges `p̂` and `f_C` for that record only;
when done, the posteriors are emitted and the record is
discarded. Memory is bounded by one record at a time.

**Strengths.**

- **Streaming-friendly.** Iterator-shaped Stage 6 fits the rest of
  the pipeline's pull-iterator pattern.
- **Embarrassingly parallel across records.** Trivial future
  optimisation: rayon over the record stream. The architecture
  spec defers parallelism here as Stage 6 is cheap relative to
  Stage 5 ([calling_pipeline_architecture.md:1765-1768](../specs/calling_pipeline_architecture.md#L1765-L1768)),
  but the per-record-independent shape leaves the door open.
- **Closed-form per-record M-step on `p̂`.** One pass over
  per-sample posteriors per iteration. Convergence in 3–5
  iterations per record (per
  [gatk_em_calculation.md §"Repeat until stable"](../specs/gatk_em_calculation.md#L182-L192)).
- **`f_C` updates only need data inside the record.** Compound
  frequency is intrinsically a within-group quantity — every
  sample contributing evidence sits in the same `MergedRecord`.

**Weaknesses.**

- **Cannot host cross-record M-steps as-is.** Contamination
  parameters `c_s` and `q_b` are estimated from many sites'
  data; they cannot converge inside a single record's EM. See
  §"Per-record vs cross-record M-steps" below for how v1 handles
  this when contamination is enabled.

### Algorithm 2 — cohort-wide EM over the whole record stream

**Shape.** Buffer the entire `MergedRecord` stream, then run a
single EM loop that touches every record per iteration.
Per-record `p̂` and `f_C` are still per-record parameters; the
loop's outer iteration covers `c_s` and `q_b` (cross-record) plus
all per-record parameters in the same outer step.

**Strengths.**

- Cleanest formulation of `c_s` / `q_b` convergence — they update
  alongside `p̂` and `f_C` in the same iteration loop.

**Weaknesses.**

- **Buffers the genome in RAM.** A cohort of 1000 samples × 10M
  variants × per-record posteriors (a few hundred bytes per
  sample per record after compression) is on the order of a TB.
  Violates the architecture's "memory bounded at every stage,
  cohort size and genome size never enter together"
  contract ([calling_pipeline_architecture.md:1770-1800](../specs/calling_pipeline_architecture.md#L1770-L1800)).
- **Wrong shape for streaming pipelines.** Breaks the iterator
  contract everywhere downstream.

**Verdict.** Rejected.

### Algorithm 3 — two-pass EM with scratch-file replay (chosen for contamination mode)

**Shape.** When `--contamination-batches` is supplied:

1. **Pass 1** (no `c_s`, no `q_b`): per-record independent EM
   exactly as Algorithm 1. As each record is processed, two
   things happen:
   - Per-batch posterior-weighted allele-count statistics for
     `q_b` and per-sample "odd reads at confidently-homozygous
     sites" statistics for `c_s` accumulate into a fixed-size
     in-RAM aggregator (memory
     `O(cohort_size × n_batches × allele_classes)` — does not
     grow with genome length).
   - The `MergedRecord` itself is **serialised to a scratch
     file** before being discarded. The scratch file is a
     length-prefixed stream of `MergedRecord`s in genomic
     order, ideally zstd-framed to match the project's `.psp`
     compression style.
2. **`q_b` and `c_s` M-step**: once the record stream is
   exhausted, perform the cross-record M-steps from the
   accumulated statistics. One closed-form M-step on each
   per-batch `q_b`; one 1D maximisation per sample for `c_s`.
3. **Pass 2** (with frozen `c_s` and `q_b`): **replay the
   scratch file** of `MergedRecord`s. Each record's E-step now
   uses the mixture likelihood derived from `scalars` + `c_s` +
   `q_{batch(s)}`. Per-record `p̂` and `f_C` re-converge under
   the new likelihoods. **This pass** emits the posteriors that
   go into the VCF.
4. **Scratch cleanup**: delete the scratch file on successful
   Pass 2 completion (or on error, unless `--keep-scratch` is
   set for debugging).

**Where does Pass 2 get its records?** Three sub-options were
considered:

- **3a. Re-read `.psp` files** — re-run the full Stage 1→Stage 5
  chain (PspReader → PerPositionMerger → VariantGrouper →
  PerGroupMerger) per pass. Lowest disk footprint (no scratch
  file), highest IO + CPU recompute cost. Wastes the Stage 4/5
  work entirely.
- **3b. Hold all `MergedRecord`s in RAM between passes** — same
  memory problem as Algorithm 2; doesn't fit at cohort scale.
- **3c. Materialise `MergedRecord`s to a scratch file during
  Pass 1, replay in Pass 2** — chosen. Pays serialisation +
  scratch disk space; avoids both the `.psp` re-read and the
  RAM blow-up. Scratch is typically **smaller than the
  combined `.psp` input** because per-variant records are
  sparser than per-covered-position records (variants are a
  small fraction of covered positions on most cohorts).

3c wins because Stage 4 and Stage 5 are the cohort side's
heaviest per-record work; replaying their output is far cheaper
than recomputing it, and the disk cost is bounded and
typically below the input footprint.

**Strengths.**

- **Bounded memory** — Pass-1 aggregators are
  `O(cohort × n_batches × allele_classes)`; the scratch file is
  on disk, not RAM. Does not grow with cohort × genome
  simultaneously.
- **Streaming-shaped** — two streaming passes over a sequential
  on-disk stream.
- **No Stage 4/5 recompute** — the expensive cohort-side work
  (allele unification, scalar projection, likelihood
  reconstruction) runs exactly once.
- **Scratch disk footprint is bounded and typically modest**
  relative to input `.psp`s.

**Weaknesses.**

- **Requires scratch disk space.** Order of `n_variants ×
  cohort_size × per-(sample,allele) bytes` after compression.
  Estimate: 10M variants × 1000 samples × 4 alleles × ~40 B per
  (sample, allele) ≈ 1.5 TB before compression, perhaps 100–300
  GB with zstd. Must be sized into the cohort-CLI documentation
  and exposed via `--scratch-dir` (defaults to project-local
  `tmp/` per the CLAUDE.md scratch-space convention; users
  pointing at a fast SSD partition is common).
- **Needs a serialisation format for `MergedRecord`** — new
  on-disk schema, even if internal to a single CLI run.
  Length-prefixed + zstd-framed `bincode`/`postcard`/similar.
  Adds a small surface area to maintain.
- **No outer iteration.** Pass 2 takes `c_s` and `q_b` as
  fixed, even though strictly a third pass would re-converge
  them against the new posteriors. In practice the architecture
  spec treats this as a fixed-point that does not need a full
  second outer iteration ("typically 3–5 rounds" of the inner
  per-record EM is enough). If real data shows the single
  outer step is insufficient, an outer iteration would re-read
  the same scratch file a third time — much cheaper than
  reopening the `.psp` chain. Add
  `--contamination-outer-iterations` if motivated — see
  §"Out-of-scope follow-ups".

### Algorithm 4 — interleaved EM (per-record + cross-record updates within one outer loop)

**Shape.** A single outer EM loop. Each outer iteration:

1. For every record, run *one* E-step under current `p̂_r`,
   `f_{C,r}`, `c_s`, `q_b`. Emit per-record soft posteriors.
2. Per-record M-step on `p̂_r` and `f_{C,r}` (one per record;
   closed-form).
3. Cross-record M-step on `q_b` and `c_s` from the cohort-wide
   accumulated statistics.

**Strengths.**

- Tightest theoretical coupling between cross-record and
  per-record parameters.

**Weaknesses.**

- **Requires keeping per-record state in RAM across outer
  iterations** — every record must remain accessible for the next
  outer iteration's E-step. Same memory problem as Algorithm 2.
- **Or rereads the scratch file per outer iteration** —
  multiplies Algorithm 3's IO by the outer-iteration count.
  Cheap per pass (scratch replay only) but earns marginal
  calibration gain.

**Verdict.** Rejected for v1; the cost/calibration trade-off does
not pencil. If real data shows Algorithm 3's single-outer-pass is
insufficient, the cheapest fix is one additional scratch-file
replay (a third pass), exposed as `--contamination-outer-iterations`.

### Algorithm 5 — VCF-replay contamination correction (deferred — full evaluation not done)

**Shape.** Treat contamination correction as a **separate
user-facing step** on top of an already-emitted, non-contamination
VCF, rather than as an internal Stage 6 mode:

1. **Standard pipeline run**: Stages 1→6 run with `c_s = 0`,
   producing a draft VCF.
2. **Contamination-recompute subcommand** (separate `pop_var_caller`
   subcommand, e.g. `pop_var_caller recompute-contamination`):
   reads the draft VCF + the original `.psp` files + the
   `--contamination-batches` definition, runs the cross-record
   M-step on `c_s` and `q_b` using VCF-derived posteriors and AD,
   then re-emits a final VCF.

**Two sub-flavours, depending on what the recompute subcommand
consumes:**

- **5a — extended VCF carries scalars.** The draft VCF's FORMAT
  fields are extended with custom entries for `q_sum`, `fwd`,
  `placed_left`, `placed_start` per (sample, allele). The
  recompute subcommand needs only the VCF — no `.psp` access.
  Drawback: invents a non-standard VCF dialect; the per-(sample,
  allele) scalar storage is text-bulky compared to the bincode
  scratch file; reinvents the scratch file as VCF.
- **5b — standard VCF + `.psp` re-read.** The draft VCF carries
  only standard FORMAT fields (GT, AD, PL/GP). The recompute
  subcommand reads it for posteriors + AD (enough to drive the
  `c_s` / `q_b` M-step), but then re-reads the `.psp` files to
  recover `q_sum` and the other scalars for the Pass 2 E-step.
  Effectively re-runs Stages 1 (already done — `.psp` exists)
  through 6 on the back of cohort-level `c_s` / `q_b` derived
  from the VCF. Drawback: Stage 4/5 work runs twice (once for
  the draft VCF, once again for the recompute).
- **5c — standard VCF + approximate likelihood reweighting.**
  The recompute subcommand uses **only** the draft VCF, treating
  the draft VCF's `PL` as `L_own` and reweighting
  multiplicatively under the mixture. This loses the
  `q_sum`-driven error-cost reweighting; the approximation
  error is small at typical `c_s < 0.05` but grows non-linearly
  at higher contamination. **Whether the approximation is
  acceptable is an open empirical question that requires real
  cohort data to answer.**

**Strengths (across sub-flavours).**

- **Decouples contamination correction from the main pipeline.**
  A user can run the no-contamination pipeline first, inspect
  the draft VCF, then opt in to contamination correction
  without re-running Stages 1–4.
- **Reusability for already-emitted cohorts.** A historical
  cohort's VCF can be re-genotyped under contamination
  correction later, without holding the original `MergedRecord`
  scratch file (which Algorithm 3 would have deleted after Pass
  2 completed).
- **Debuggability.** The draft VCF is human-readable and
  inspectable.
- **Composability.** The recompute subcommand is a self-contained
  CLI artefact that fits other workflows (e.g. running on a
  collaborator's VCF without their `.psp`s — viable only in
  flavour 5a or 5c).

**Weaknesses (across sub-flavours).**

- **No scenario where it strictly beats Algorithm 3 on
  efficiency.** Flavour 5a reinvents the scratch file with
  parsing overhead; 5b runs Stages 4/5 twice; 5c trades exact
  math for an unmeasured approximation.
- **The cross-record M-step from VCF-derived posteriors is
  doable** (GP + AD + GT + AF suffices for the spec's `c_s`/`q_b`
  shape — modulo the still-unresolved `c_s` sufficient-statistic
  question that gates Algorithm 3 as well), but **the Pass 2
  E-step is where the information loss bites** — either by
  inventing a VCF dialect (5a), re-doing Stage 5 (5b), or
  approximating (5c).

**Where Algorithm 5 genuinely shines** — as a *user-facing
product feature* rather than as a *Stage 6 internal design*:
post-hoc contamination correction on an emitted VCF is a real
tooling use case (especially flavour 5c if its approximation
budget pencils). Algorithm 3 and Algorithm 5 are **not mutually
exclusive** — Stage 6 could ship the scratch-file internal path
for fresh runs, and a separate `recompute-contamination`
subcommand for post-hoc.

**Verdict.** Recorded as a real alternative; **full evaluation
deferred** along with all other contamination-mode work (see
§"Out of scope (v1)"). The contamination follow-up plan will
revisit Algorithms 3 and 5 side by side with the benefit of
real cohort data — including the open question of whether
flavour 5c's approximation is calibration-acceptable.

### Decision

- **No-contamination mode (default):** Algorithm 1
  (per-record-independent EM, streaming). **Implemented in v1.**
- **Contamination mode (`--contamination-batches` supplied):**
  choice between Algorithm 3 (two-pass with scratch-file replay)
  and Algorithm 5 (VCF-replay, possibly as a separate
  subcommand) is **deferred** to the contamination follow-up
  plan, alongside the unresolved `c_s` sufficient-statistic
  question. Algorithm 3 is the current best candidate for the
  in-engine path; Algorithm 5 is the candidate for a separate
  post-hoc subcommand; the two are not mutually exclusive and
  may both ship. **Not implemented in v1.**

The split is invisible from the iterator's external shape — the
v1 engine emits a `PosteriorRecord` per `MergedRecord` in genomic
order. Whatever contamination machinery later lands plugs in
behind that shape.

## API shape

```rust
// src/var_calling/posterior_engine.rs
use crate::var_calling::per_group_merger::{MergedRecord, PerGroupMergerError};

pub struct PosteriorEngine<I>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
{ /* fields private */ }

impl<I> PosteriorEngine<I>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
{
    /// Construct a no-contamination posterior engine.
    pub fn new(upstream: I) -> Self;

    pub fn with_config(upstream: I, config: PosteriorEngineConfig) -> Self;

    pub fn config(&self) -> &PosteriorEngineConfig;
}

impl<I> Iterator for PosteriorEngine<I>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
{
    type Item = Result<PosteriorRecord, PosteriorEngineError>;
}

#[derive(thiserror::Error, Debug)]
#[non_exhaustive]
pub enum PosteriorEngineError {
    #[error("upstream: {0}")]
    Upstream(#[from] PerGroupMergerError),

    /// EM produced a non-finite posterior (NaN, +∞, or all -∞)
    /// for the named sample. The closed-form math is finite for
    /// all valid `MergedRecord` inputs; surfacing this signals an
    /// internal bug, not a data condition.
    #[error(
        "non-finite posterior at chrom {chrom_id} {start}-{end} \
         for sample_idx {sample_idx}: {kind:?}"
    )]
    NonFinitePosterior {
        chrom_id: u32,
        start: u32,
        end: u32,
        sample_idx: usize,
        kind: NonFiniteKind,
    },

    /// EM did not converge within `max_iterations`. Configurable
    /// max defaults to 50 (per-record); a hard cap exists so that
    /// pathological records cannot stall the engine.
    #[error(
        "EM did not converge within {max_iterations} iterations \
         at chrom {chrom_id} {start}-{end} (last delta: {last_delta})"
    )]
    DidNotConverge {
        chrom_id: u32,
        start: u32,
        end: u32,
        max_iterations: u32,
        last_delta: f64,
    },
}
```

**Errors latch.** Per the merger/grouper precedent
([per_position_merger.rs:130-135](../../src/var_calling/per_position_merger.rs#L130-L135)):
once any error surfaces, subsequent `next()` calls return `None`.

**Contamination-mode API.** When `--contamination-batches` is in
play, the engine takes the same single-pass upstream iterator
but materialises each `MergedRecord` to a scratch file during
Pass 1, then replays the scratch file in Pass 2. The
caller-facing shape:

```rust
pub struct ContaminationAwarePosteriorEngine<I>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
{ /* fields private */ }

impl<I> ContaminationAwarePosteriorEngine<I>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
{
    pub fn new(
        upstream: I,
        contamination_config: ContaminationConfig,
        scratch_dir: PathBuf,
    ) -> Self;

    /// Consumes the upstream (Pass 1), runs the cross-record
    /// M-step, then opens an iterator that replays the scratch
    /// file (Pass 2).
    pub fn run(self) -> impl Iterator<Item = Result<PosteriorRecord, PosteriorEngineError>>;
}
```

The upstream iterator is consumed exactly once — no re-open
trait required. Pass 2 reads from the scratch file written
during Pass 1.

**Scratch file format.** Length-prefixed `MergedRecord`s serialised
via `bincode` (or `postcard`), framed in zstd-compressed blocks
matching the `.psp` writer's convention. The exact format is
internal to a single Stage 6 run (it never persists across
invocations, and never crosses tool versions), so schema
evolution is not a concern — the same binary writes and reads
it within the same `pop_var_caller` process. Detailed format
specification deferred to the contamination-mode follow-up plan.

## Output shape

```rust
pub struct PosteriorRecord {
    pub chrom_id: u32,
    pub start: u32,
    pub end: u32,
    pub alleles: Vec<MergedAllele>,    // forwarded from MergedRecord
    pub ploidy: u8,
    /// Estimated per-allele frequencies at this site after EM.
    /// Length equals `alleles.len()`. Sums to 1.
    pub allele_frequencies: Vec<f64>,
    /// Estimated cohort frequency for each chain-anchored compound
    /// in the record, indexed by `alleles.iter().position(|a| a.is_compound)`.
    /// `None` for non-compound alleles.
    pub compound_frequencies: Vec<Option<f64>>,
    /// Per-sample per-genotype posterior. `posteriors[sample_idx]`
    /// is a vector of length `genotype_count(ploidy, alleles.len())`
    /// in canonical [`genotype_order`] order. Each row sums to 1.
    pub posteriors: Vec<Vec<f64>>,
    /// Per-sample best genotype (argmax of `posteriors[sample_idx]`).
    /// Indexed by genotype-enumeration position.
    pub best_genotype: Vec<usize>,
    /// Per-sample genotype quality, Phred-scaled `−10 log10 (1 − P(best))`.
    pub gq_phred: Vec<f64>,
    /// Site-level QUAL: Phred of `Π_s P(hom-ref)_s`.
    pub qual_phred: f64,
    /// EM diagnostics — useful for benchmarking and bug hunting.
    /// Kept opaque (boxed) so PosteriorRecord doesn't grow unnecessarily.
    pub diagnostics: EmDiagnostics,
}

pub struct EmDiagnostics {
    pub iterations: u32,
    pub final_max_delta_p: f64,
    pub converged: bool,
}
```

**Open question — VCF writer interface.** The VCF writer is a
separate module (out of scope for this plan). The fields above are
chosen to cover everything VCF emission needs (GT from
`best_genotype`, GQ from `gq_phred`, QUAL from `qual_phred`, AF from
`allele_frequencies`, PL from `−10 log10 posteriors`). The writer
may also want raw per-sample scalars (DP, AD) — those come from
forwarding `MergedRecord.scalars`. Decide at writer-plan time
whether `PosteriorRecord` should carry forwarded scalars or whether
the writer joins on (chrom, pos) against a parallel stream.

## Algorithm details — no-contamination mode (Algorithm 1)

Per `MergedRecord`:

### Step 1 — initialise

- `p̂[a] = 1 / n_alleles` for every allele `a` (flat).
- `f̂_C = 1 / n_alleles` for every chain-anchored compound `C` in
  the record (flat). **Open question** — flat or
  Dirichlet-pseudocount-derived? GATK uses flat for the first
  iteration to "avoid suppressing real variants before they
  accumulate evidence"
  ([gatk_em_calculation.md §"Convergence"](../specs/gatk_em_calculation.md#L80-L84)).
  Mirror that.
- Compute per-allele SNP-vs-indel classification once (allele
  byte-length and content vs the REF span). This routes the
  per-allele Dirichlet pseudocount (`α_alt = 0.01` SNP,
  `α_alt = 0.00125` indel) for the M-step.

### Step 2 — E-step

For sample `s` and each genotype `G` in `genotype_order(ploidy, n_alleles)`:

```
log_posterior_unnorm(s, G) = log_likelihood(s, G)
                           + log_prior(G | p̂, F_s, f̂_C, ca_flags[s])
```

where:

- `log_likelihood(s, G)` comes directly from
  `MergedRecord.log_likelihoods[s][G]` — Stage 5 has already
  computed it under `c_s = 0`.
- `log_prior(G | p̂, F_s, f̂_C, ca_flags[s])` is the HWE-with-`F`
  prior at allele frequencies `p̂`, with the per-sample inbreeding
  coefficient `F_s` and — for compound alleles in `G` — the
  cohort-derived compound frequency `f̂_C` substituted for the
  per-allele `p̂[a]` of that compound.
  - When `ca_flags[s][a] = true` (sample is chain-broken at
    compound allele `a`), the prior uses `f̂_C` as the
    compound's frequency parameter; the
    constituent-product likelihood path in
    `log_likelihoods[s][G]` already accounts for the
    decomposition (per Stage 5
    [§"Chain-broken-compound case"](cohort_per_group_merger.md#L696-L721)).

Normalise via log-sum-exp:

```
log_Z(s) = logsumexp_G(log_posterior_unnorm(s, G))
posterior(s, G) = exp(log_posterior_unnorm(s, G) − log_Z(s))
```

### Step 3 — M-step on `p̂`

Per-record allele-frequency update from posterior-weighted allele
counts plus Dirichlet pseudocounts:

```
E[n_k] = Σ_s Σ_G posterior(s, G) · (copies of allele k in G)
α_k = α_ref if k == 0 else (α_snp_alt if allele k is SNP else α_indel_alt)
p̂_k ← (α_k + E[n_k]) / Σ_j (α_j + E[n_j])
```

Closed-form, one pass over samples × genotypes.

### Step 4 — M-step on `f̂_C` (one per chain-anchored compound)

For each compound allele `C` in `alleles`:

```
E[n_C] = Σ_s Σ_G posterior(s, G) · (copies of C in G)
n_C_chromosomes_total = 2 · n_samples       // diploid case; generalises
f̂_C ← (α_compound + E[n_C]) / (α_compound + α_NotC + n_C_chromosomes_total)
```

with `α_compound = 0.001` (default; overridable per spec). Note that
both chain-evident and chain-broken samples contribute their
posteriors to `E[n_C]` — chain-evident with full-resolution
likelihood, chain-broken via the constituent-product fallback
already baked into `MergedRecord.log_likelihoods`.

### Step 5 — convergence check

```
max_delta = max over alleles a of |p̂_new[a] − p̂_old[a]|
if max_delta < CONVERGENCE_THRESHOLD: stop
```

Default `CONVERGENCE_THRESHOLD = 1e-4` (revisit against real data;
GATK uses 0.1 on raw allele *counts*, which is roughly equivalent
at depth ~20 samples). Hard iteration cap
`MAX_ITERATIONS = 50`; exceeding it returns
`PosteriorEngineError::DidNotConverge`.

### Step 6 — emit `PosteriorRecord`

- `posteriors[s]` = final per-sample posterior vector.
- `best_genotype[s]` = argmax of `posteriors[s]`.
- `gq_phred[s] = −10 · log10(1 − posteriors[s][best_genotype[s]])`,
  clamped at some max (e.g. 99) to avoid `+∞` from
  posterior = 1 exactly.
- `qual_phred = −10 · log10(Π_s posteriors[s][hom_ref_index])`,
  computed in log-space as
  `−10 · log10_e · Σ_s log(posteriors[s][hom_ref_index])`.

### Edge cases

- **Single-allele record (every sample hom-REF after EM)** —
  Stage 5 already drops these
  ([cohort_per_group_merger.md §"Step 6 — emit the merged record"](cohort_per_group_merger.md#L730-L738)),
  so Stage 6 never receives them. Belt-and-braces: if `alleles.len()
  < 2` slips through, emit a `PosteriorRecord` with trivial
  posteriors and QUAL = 0 (no variant).
- **Sample with zero observations across the record** — its
  `log_likelihoods[s]` is all-zeros (uniform over genotypes);
  the prior alone determines its posterior. Behaves correctly
  with no special-casing.
- **`posteriors[s][best_genotype[s]] == 1.0` exactly** — would
  divide by zero in the Phred calculation. Clamp the input to
  `1.0 − f64::EPSILON` before taking the log.
- **HWE-with-`F` at higher ploidies.** The architecture spec
  gives the diploid formula explicitly
  ([calling_pipeline_architecture.md:1664-1689](../specs/calling_pipeline_architecture.md#L1664-L1689)).
  For polyploids, the prior generalises to the
  Wright–Fisher partition: P(genotype = `a^{k_a}`) under HWE-with-`F`
  is a mixture of "all `n` copies IBD" (probability `F`,
  homozygous-for-some-allele weighted by `p̂`) and "all `n` copies
  drawn independently" (probability `1 − F`, multinomial on `p̂`).
  **Open question** — does v1 actually need polyploids, or is
  diploid-only acceptable for the first ship? The walker
  produces ploidy-agnostic scalars, and Stage 5 supports the
  configured ploidy; locking Stage 6 to diploid would be the
  simplest v1.

## Algorithm details — contamination mode (Algorithm 3)

### Pass 1 — first-pass per-record EM, accumulate cross-record statistics, materialise records

Per `MergedRecord` (same as Algorithm 1, Steps 1–5), but:

- `c_s = 0` initially → use `MergedRecord.log_likelihoods`
  directly (no mixture-likelihood recompute needed).
- After convergence on `p̂` and `f_C`, **serialise the
  `MergedRecord` to the scratch file** (zstd-framed bincode or
  postcard, length-prefixed, written in genomic order) and
  accumulate into cross-record aggregators:

  ```rust
  pub struct CrossRecordAggregator {
      /// Per-batch posterior-weighted allele counts: indexed by
      /// (batch_idx, allele_class). `allele_class` is a tiny
      /// fixed-size enum (REF, SNP_alt, INDEL_alt) — `q_b` is
      /// over allele classes, not per-record allele indices
      /// (alleles are not comparable across records).
      pub q_b_allele_counts: Vec<[f64; N_ALLELE_CLASSES]>,

      /// Per-sample "odd reads at confidently-homozygous sites"
      /// statistics for c_s estimation. Per the architecture
      /// spec (line 1883-1886), c_s is driven by reads
      /// disagreeing with the sample's soft-assigned homozygous
      /// genotype.
      pub c_s_evidence: Vec<CSEvidence>,
  }
  ```

  **Open question — what exactly accumulates for `c_s`?** The
  architecture spec says "1D maximisation per sample of the
  mixture likelihood against `q_{batch(s)}` and the sample's
  soft-assigned genotypes." Pinning the exact sufficient
  statistic needs a separate side-derivation that this draft
  plan does not do. Candidates:

  1. **Per-(sample, allele_class) posterior-weighted counts of
     non-genotype-supporting reads.** Lightest. Loses the
     mixture-likelihood shape (`c_s` would be solved on a
     simplified objective).
  2. **Per-sample list of (record_idx, scalars,
     soft_posterior)** — full mixture-likelihood recompute at
     `c_s` M-step time. Most faithful but heaviest in
     pass-1 RAM.
  3. **Per-sample "expected error fraction" running stats** — a
     pre-aggregated form somewhere between (1) and (2).

  Decide in a side spec; until then, this plan ships v1 with the
  contamination machinery scoped *out* (see §"Out of scope
  v1") and only the no-contamination path (Algorithm 1) lands
  in code.

### Cross-record M-steps (after Pass 1)

(Scoped out of v1; documented here so the plan has the full
architecture in writing.)

- `q_b ← Dirichlet-normalised allele class counts per batch`,
  with `α_q` pseudocount per class
  (`--contamination-source-pseudocount`, default = SNP alt
  pseudocount).
- `c_s ← argmax_{c ∈ [0, 1]} L_s(c | q_b, soft_genotypes)`,
  1D over a fixed grid (e.g. 0.001-step) or via golden-section
  search. Held at 0 for singleton batches and batches below
  `--min-batch-size-for-contamination` (default 5).

### Pass 2 — replay scratch file with frozen `c_s` and `q_b`

Per `MergedRecord` (deserialised from the scratch file in
genomic order):

- E-step at iteration 0 recomputes per-sample log-likelihoods
  from `MergedRecord.scalars`, `MergedRecord.other_scalars`,
  `c_s`, and `q_b`, using the **mixture likelihood**:

  ```
  L_s(read | G, c_s, q_b) = (1 − c_s) · L_own(read | G)
                          + c_s · L_contam(read | q_b)
  ```

  with `L_own` reconstructed from `scalars` exactly as Stage 5
  did internally (the `c_s = 0` precomputed
  `log_likelihoods` table is the special case; Stage 6 needs
  the general form here). **Open question** — does this require
  a shared helper that lives in Stage 5 (so both stages call the
  same scalar→likelihood code), or does Stage 6 carry its own
  copy? The shared-helper route is cleaner but couples the two
  stages tighter.
- Re-converge `p̂` and `f_C` against the new likelihoods.
- Emit `PosteriorRecord`s. **These** are the values that go to
  the VCF writer.

On successful completion of Pass 2 the scratch file is deleted.
On error the file is retained iff `--keep-scratch` was set
(useful for debugging the contamination-mode path); otherwise
it is deleted by the engine's `Drop` impl.

## Configurable parameters

```rust
pub const DEFAULT_CONVERGENCE_THRESHOLD: f64 = 1e-4;
pub const DEFAULT_MAX_ITERATIONS: u32 = 50;
pub const DEFAULT_REF_PSEUDOCOUNT: f64 = 10.0;
pub const DEFAULT_SNP_ALT_PSEUDOCOUNT: f64 = 0.01;
pub const DEFAULT_INDEL_ALT_PSEUDOCOUNT: f64 = 0.00125;
pub const DEFAULT_COMPOUND_ALT_PSEUDOCOUNT: f64 = 0.001;
pub const DEFAULT_INBREEDING_COEFFICIENT: f64 = 0.0;
pub const MAX_GQ_PHRED: f64 = 99.0;

#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct PosteriorEngineConfig {
    pub convergence_threshold: f64,
    pub max_iterations: u32,
    pub ref_pseudocount: f64,
    pub snp_alt_pseudocount: f64,
    pub indel_alt_pseudocount: f64,
    pub compound_alt_pseudocount: f64,
    /// Per-sample inbreeding coefficient vector. Length equals
    /// the cohort size; entries default to the CLI scalar
    /// `--inbreeding` (default 0.0). Per the architecture spec
    /// [§"From likelihood to posterior" item 2](../specs/calling_pipeline_architecture.md#L1690-L1704)
    /// we keep this per-sample from day one even though the
    /// user-facing flag is a scalar, so per-sample F via a
    /// supplied file is a UI change later, not a refactor.
    pub fixation_index_per_sample: Vec<f64>,
    pub max_gq_phred: f64,
    /// When `true`, use precomputed lookup tables for the EM's
    /// pure inner-loop functions (multinomial coefficients,
    /// HWE-with-`F` log-prior at quantised `(p, F)` grids, etc.).
    /// Trades a small calibration-error budget for inner-loop
    /// speed. See §"Approximation via precomputed lookup tables
    /// (evaluation)" — the flag is wired from day one but the
    /// approximating implementation lands only after the
    /// evaluation gates it.
    pub approximate_posterior_calculation: bool,
    pub contamination: Option<ContaminationConfig>,
}
```

CLI binding (when the cohort subcommand lands):
`--approximate-posterior-calculation` (default **`true`** —
flipped explicitly to opt out with
`--approximate-posterior-calculation=false`). Default-on assumes
the evaluation below confirms the speedup is real and the error
budget is tight; if evaluation falls short, flip the default to
`false` before the cohort CLI ships.

**`ContaminationConfig` is scoped out of v1** — see §"Out of scope".

## Out of scope (v1)

- **All contamination machinery and its evaluation.** Both
  Algorithm 3 (in-engine scratch-file replay) and Algorithm 5
  (VCF-replay / post-hoc subcommand) are drafted above so the
  architecture is visible end-to-end, but **neither is
  implemented in v1**, and **the choice between them is also
  deferred**. `PosteriorEngineConfig::contamination` exists as
  a typed `Option`, but `Some(_)` returns
  `PosteriorEngineError::ContaminationModeUnimplemented` (a new
  error variant). The full design — including the choice
  between Algorithms 3 and 5, the `c_s` sufficient-statistic
  question, and whether flavour 5c's approximation is
  calibration-acceptable — lands in a follow-up plan driven by
  real cohort data.
- **VCF emission.** Separate module
  (`var_calling::vcf_writer` or similar). Stage 6 emits
  `PosteriorRecord`; the writer turns it into VCF lines.
- **CLI wiring.** `--inbreeding`, `--ref-pseudocount`, etc. get
  bound when the cohort CLI subcommand lands.
- **Parallelism across records.** v1 is sequential per the
  architecture spec
  ([calling_pipeline_architecture.md:1765-1768](../specs/calling_pipeline_architecture.md#L1765-L1768)).
  Trivial rayon-over-records is a follow-up if profiling shows it
  matters (each record's EM is independent in non-contamination
  mode).
- **Cohort-only inference of compounds with no chain anchor
  anywhere.** v2 extension; the prior-side machinery is the same
  as the chain-anchored case but the allele-unification step
  (Stage 5) does not propose these candidates. See
  [calling_pipeline_architecture.md §"Future extension: cohort-only inference"](../specs/calling_pipeline_architecture.md#L1910).
- **Variational-Bayes upgrade for rare-variant calibration.**
  Mentioned in the architecture spec
  ([calling_pipeline_architecture.md:1905-1908](../specs/calling_pipeline_architecture.md#L1905-L1908))
  as a future upgrade if needed; not in v1.
- **Per-sample ploidy.** Inherits Stage 5's single-cohort-ploidy
  constraint.

## Test strategy

Stage 6 tests live in
`src/var_calling/posterior_engine.rs`'s `#[cfg(test)]` module and
use synthetic `MergedRecord` fixtures (no Stage 5 round-trip in
unit tests; integration tests cover that).

### Required fixture builders

Add to the existing var_calling test-helpers module:

```rust
fn merged_record_simple(
    chrom_id: u32,
    pos: u32,
    alleles: Vec<&[u8]>,
    ploidy: u8,
    likelihoods: Vec<Vec<f64>>,    // per-sample, per-genotype
) -> MergedRecord;

fn merged_record_with_compound(
    chrom_id: u32,
    pos: u32,
    alleles: Vec<&[u8]>,
    compound_indices: Vec<usize>,
    chain_anchor_flags: Vec<Vec<bool>>,
    likelihoods: Vec<Vec<f64>>,
    ploidy: u8,
) -> MergedRecord;
```

### Unit tests (no-contamination mode)

1. **Single-sample, single-SNP, strong evidence.** Two alleles,
   one sample with `log_likelihoods = [0.0, -50.0, -50.0]` (REF/REF
   strongly preferred). EM converges to `p̂ ≈ (1, 0)` (clamped by
   pseudocount); posterior is REF/REF; QUAL low.
2. **Single-sample, strong alt evidence.** `log_likelihoods =
   [-50.0, -50.0, 0.0]` (homozygous alt). Posterior is alt/alt;
   QUAL high.
3. **Single-sample, ambiguous.** `log_likelihoods = [-1.0, -1.0,
   -1.0]` (uniform). Posterior matches the HWE-with-F prior at
   the flat-initialised `p̂`.
4. **Two samples, opposite alt evidence.** S0 strong REF/REF, S1
   strong alt/alt. `p̂` converges to (0.5, 0.5). Posteriors
   reflect each sample's data.
5. **Pooling effect across cohort.** Twenty samples; 19 strong
   REF/REF, 1 weak het. Compare two scenarios:
   - Sample 20's likelihood vector strongly suggests het.
   - Sample 20's likelihood vector weakly suggests het.
   In both cases, the rare-allele pseudocount pulls `p̂[alt]`
   down; the weak-evidence sample's posterior shifts toward
   REF/REF as a result. The strong-evidence sample's posterior
   does not.
6. **EM convergence in 3–5 iterations** (assert on `EmDiagnostics.iterations`).
7. **Hard non-convergence** — a pathological likelihood matrix
   (e.g. all `-∞`) returns `PosteriorEngineError::DidNotConverge`.
   **Open question** — can valid `MergedRecord` inputs actually
   produce non-convergence? If not, this test exercises only the
   error path.
8. **F = 0 vs F = 1 prior effect.** Same input likelihoods; F = 0
   produces HWE-shaped priors; F = 1 produces homozygote-only
   priors (no het mass). Posteriors shift accordingly.
9. **SNP vs indel pseudocount routing.** A site with both a SNP
   alt and an indel alt. After EM, `p̂[snp_alt]` is pulled toward
   the SNP pseudocount and `p̂[indel_alt]` toward the indel
   pseudocount.
10. **Chain-anchored compound, all samples chain-evident.** EM
    converges on `f̂_C` from the compound posteriors; per-sample
    posteriors use `f̂_C` in the prior.
11. **Chain-anchored compound, mixed chain-evident +
    chain-broken samples.** Chain-broken samples (`ca_flags = true`)
    use the cohort-derived `f̂_C` as their compound prior; their
    posteriors lean on the cohort estimate.
12. **All samples chain-broken at a compound.** `f̂_C` estimated
    entirely from the constituents-independent likelihood;
    converges but with high uncertainty. Assert numerical
    stability (no NaN, no +∞).
13. **QUAL calculation.** Site with N samples all weak hom-ref
    (e.g. P(hom-ref) = 0.95 each); QUAL = `-10 · log10(0.95^N)`.
    Direct check against the closed-form.
14. **GQ calculation.** Sample with posterior = (0.99, 0.005,
    0.005); GQ = `-10 · log10(0.01) = 20`.

### Property tests (proptest)

1. **Posteriors sum to 1 per sample** (within `1e-9`) for every
   `MergedRecord` shape proptest generates.
2. **`p̂` is a valid simplex** (non-negative, sums to 1, within
   `1e-9`).
3. **Permuting samples does not change `p̂`, QUAL, or
   sample-set-invariant posteriors.** EM is symmetric in sample
   order.
4. **Pseudocount monotonicity.** Increasing `α_ref` (while
   holding other pseudocounts) cannot increase `p̂[alt]` for any
   alt allele.

### Integration tests

Add to `tests/`:

1. **End-to-end (PspReader → ... → PosteriorEngine)** on a
   small synthetic cohort (3 samples × 1 chromosome × 100 bp).
   Assert: posterior records emitted in genomic order; expected
   genotypes for the synthetic input.
2. **Cohort-pooling test.** 5 samples × 1 site, 1 sample with
   weak alt evidence, 4 with strong REF — assert the weak-evidence
   sample is called REF (cohort prior dominates). Mirror as 5
   samples with weak alt evidence — assert all 5 called het
   (pooled evidence overcomes the rare-allele prior).
3. **Reproducibility.** Running the engine twice on the same
   input produces bit-identical output (no time- or
   thread-dependent randomness).

### Reproducing legacy-test scenarios

The legacy `genotype_posteriors.rs` tests
([src/genotype_posteriors.rs:tests](../../src/genotype_posteriors.rs))
encode a handful of EM correctness scenarios that v1 should
reproduce conceptually. These do **not** port verbatim (different
input shape, natural-log vs log10, no PL fallback), but the
scenarios are:

- EM converges flat-`p̂` to data-driven `p̂` in a few iterations.
- Posteriors shift between iteration 0 and convergence in the
  expected direction (weak alt + cohort REF → REF; consistent
  cohort het → het).
- QUAL distinguishes "every sample weak REF" from "few samples
  strong alt".

Pick representative scenarios and reproduce as Stage 6 unit tests
on synthetic `MergedRecord` input. **No formal porting checklist**
(unlike Stage 5's gVCF-merger port) because the legacy
`genotype_posteriors.rs` is being deleted regardless of v1's test
coverage — its scenarios are background reference, not a
correctness oracle.

## Validation

Inside the dev container (`./scripts/dev.sh`):

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib --tests`
- `cargo build --examples --benches`

End-to-end validation: run the full pipeline (Stage 1 → ... →
Stage 6) on a few small test cohorts and inspect the resulting
posteriors. Compare against the legacy gVCF path on the same
cohort, *after* converting per-sample BAMs through both
pipelines. Disagreement is expected at:

- **Chain-anchored compounds.** The new path calls them; the old
  path does not.
- **Sub-threshold single-read observations.** The new path keeps
  these; the old path filters them.
- **Inbreeding-affected cohorts.** The old path has no F prior;
  the new path with `--inbreeding > 0` shifts het calls.

Agreement should be tight on cohorts that exercise none of the
above (outcrossing, common-SNP-only, no compounds).

A criterion bench is not in v1 of this plan; Stage 6 in
no-contamination mode is structurally cheaper than Stage 5 per
record. Add the bench when the cohort CLI lands and there is a
realistic end-to-end pipeline to drive it.

## Assumptions / silent choices

- **Natural-log throughout.** `MergedRecord.log_likelihoods` is
  natural-log; the engine stays in natural-log. Phred conversions
  (GQ, QUAL) convert to log10 at the boundary using
  `log10_e = 1.0 / 10.0_f64.ln()`. The legacy `genotype_posteriors.rs`
  used log10 internally; we deliberately do not.
- **Flat `p̂` initialisation.** Matches GATK
  ([gatk_em_calculation.md §"Convergence"](../specs/gatk_em_calculation.md#L80-L84)).
  Avoids the prior dominating the first iteration before evidence
  accumulates.
- **Per-record EM, never cross-record (in no-contamination mode).**
  Each `MergedRecord` is independent. Contamination mode breaks
  this; documented above.
- **HWE-with-`F` per-sample F vector.** Internal representation
  stored per-sample even when the CLI flag is a scalar, per the
  architecture spec.
- **Phred clamping at GQ 99.** Avoids `+∞` from posterior = 1
  exactly. Standard convention in VCF outputs.
- **Posterior sums normalised via logsumexp.** Numerical
  stability at low-likelihood genotypes.
- **`max_iterations = 50`.** Belt-and-braces cap. Real records
  converge in 3–5 per GATK's empirical experience.
- **SNP vs indel classification per merged allele.** Computed
  once at EM init, cached for the duration of the record.

## Risks

- **`c_s` sufficient-statistic question is unresolved.** Pins
  the contamination mode scope. The no-contamination path is
  unaffected.
- **HWE-with-`F` at higher ploidies.** The diploid formula is
  pinned in the spec; the polyploid generalisation is not. v1
  may need to ship diploid-only and defer higher ploidies until
  the formula is written down.
- **Numerical stability at very deep coverage.**
  `log_likelihoods` values can be very negative (e.g. `-10000`
  at depth 10000); logsumexp must handle this. The standard
  shift-by-max trick suffices.
- **Convergence on adversarial inputs.** A `MergedRecord` whose
  log-likelihood matrix has heavy symmetry could oscillate. Hard
  cap protects from infinite loops; surfacing
  `DidNotConverge` lets ops see it. **Open question** — should
  the engine fall back to emitting the last iterate as a
  best-effort answer, with a flag, instead of erroring? GATK
  emits the last iterate; freebayes does too.
- **Hidden coupling with Stage 5's likelihood internals.**
  Contamination mode re-implementing the
  scalar→likelihood path duplicates logic that Stage 5 already
  owns. The shared-helper option (see §"Pass 2") would couple
  the stages; the duplicated-code option drifts. Decide before
  the contamination plan.

## Approximation via precomputed lookup tables (evaluation)

The EM inner loop is dominated by a handful of **pure functions**
— same input always gives the same output, no I/O, no shared
mutable state. Some of them take integer-typed inputs from a
small fixed range; some take floats from a bounded range. The
classic perf move at this shape is to precompute a quantised
lookup table once (either as a `const` baked into the binary via
`build.rs` or as a `OnceLock` populated on first use) and replace
the hot-path computation with a table read + linear
interpolation. CLI flag
`--approximate-posterior-calculation` (default `true`, see
§"Configurable parameters") toggles the approximate path.

**This section is an evaluation task, not a v1 implementation
commitment.** The exact-math path ships first; the LUT path is
opt-in once benchmarking proves it pays.

### Candidate pure functions, ranked by expected payoff

1. **Multinomial coefficient `C(ploidy; n_1, …, n_k)`.** Pure
   integer function. With diploid/triploid/tetraploid and
   `max_alleles = 6`, the input space is the genotype-tuple set
   from `genotype_order(ploidy, n_alleles)` — fewer than 200
   entries total across all (ploidy, n_alleles) combinations the
   engine ever sees. **Trivially a static table** baked at
   `genotype_order` construction; no quantisation needed; no
   precision loss. Almost certainly worth it.
2. **`lgamma(n)` for integer `n` in `[0, MAX_DEPTH]`.** The
   multinomial expansion uses `lgamma(n + 1)` for non-negative
   integer `n` bounded by per-allele observation counts.
   `MAX_DEPTH` is `max_snp_column_depth = 8000` from the walker
   config. A `[f64; 8001]` table is 64 kB — fits in L2. Exact
   for the integer-input case; no quantisation. Likely worth it.
3. **HWE-with-`F` log-prior per genotype.** For diploid
   biallelic, the prior is `P(AA) = (1−p)² + F·p(1−p)`,
   `P(AB) = 2p(1−p)(1−F)`, `P(BB) = p² + F·p(1−p)` — a function
   of `(p, F)`. Quantise `p ∈ [0, 1]` to 1024 steps (10 bits)
   and `F ∈ [0, 1]` to 64 steps (6 bits): 3 genotypes × 1024 ×
   64 × 8 B = **1.5 MB**. Linear interpolation in `p` for
   accuracy. For triallelic / higher the `p` lives on an
   `n_alleles − 1` simplex and the table size blows up
   combinatorially — only the diploid-biallelic case is
   table-friendly; for ≥3 alleles fall back to direct
   computation. Worth evaluating; the diploid-biallelic site is
   the cohort's bulk by record count.
4. **Phred conversion `−10 · log10(1 − p)` for GQ.** Pure
   function of `p ∈ [0, 1 − ε]`. A `[f64; 4096]` table with
   linear interpolation is 32 kB and gives ≈3 Phred decimal-digit
   accuracy. Cheap, but called only once per (sample, record),
   not in the hot inner loop — payoff is small in absolute terms.
5. **`xlogy(n, p) = n · ln(p)` (with `0 · ln 0 = 0`).** Used in
   the multinomial. Tabulating helps only if both arguments
   quantise well; `n` is bounded integer (good) but `p` is a
   quantised `p̂` allele frequency. Same caveat as item 3 —
   only the biallelic case fits in a small table.

### Non-candidates

- **`log(x)`, `exp(x)`, `logsumexp(values)`.** Already
  hardware-fast via the FPU's `log`/`exp` intrinsics; a
  table-with-interpolation typically loses to hardware on
  modern x86_64. Leave alone.
- **Convergence delta and argmax.** Trivial scalar comparisons.

### Approximation budget

How much calibration loss can the engine tolerate? The output
quantities that matter:

- **GT (argmax genotype):** unaffected as long as the
  approximation does not flip the rank of the top-two
  posteriors. Empirically the gap between best and second-best
  posterior at a confident call is many decimal digits; a
  3rd-decimal-digit approximation error cannot flip the
  argmax.
- **GQ (Phred of `1 − P(best)`):** the integer Phred output
  rounds to whole units in the VCF, so any approximation that
  preserves the first decimal digit of `P(best)` is invisible.
- **QUAL:** integrates over the whole cohort; per-sample
  approximation errors partially cancel. A ≤1 Phred-unit budget
  cohort-wide is the bar.
- **AF (allele frequency `p̂`):** EM iterates on `p̂` itself,
  so approximation error in the HWE prior feeds back into `p̂`.
  Has the highest sensitivity. Worth a dedicated test (compare
  exact vs approximate `p̂` on the same input cohort).

### Evaluation steps

Before committing to the LUT path:

1. **Ship the exact-math engine** (this plan's v1 scope).
   Lock the API surface and the test suite.
2. **Add a criterion bench** for the EM inner loop on a
   realistic record (e.g. 1000 samples × 4 alleles × ploidy 2).
   Measure where time goes — `cargo flamegraph` or `samply`
   profile of `posterior_engine::run_em_for_record` on the
   bench input. If the hot path is not in the candidates
   above (e.g. it turns out to be in the logsumexp normalisation
   that is already hardware-optimal), the LUT idea is not
   worth pursuing; document and move on.
3. **Prototype** each promising candidate as a feature flag,
   guarded by `approximate_posterior_calculation`. Multinomial
   coefficient first (item 1; no precision loss, smallest code
   change). Then `lgamma` (item 2). Then HWE-with-`F` for
   diploid-biallelic (item 3 — the trickiest).
4. **Measure two things per candidate:**
   - **Speed.** Criterion: `posterior_engine_em/{exact,
     approximate}/{record_shape}` benches; expect a real
     win on the inner-loop-bound shapes.
   - **Calibration.** A side test that runs the exact and
     approximate engines on the same `MergedRecord` input
     and asserts: GT agreement 100 %, GQ within 1 Phred,
     QUAL within 1 Phred cohort-wide, `p̂` within
     `5e-4` per allele. Tighten the bounds if real data
     allows.
5. **Decide per candidate.** Land the ones that pay; drop the
   ones that do not. Document each decision in the
   implementation report.
6. **Set the CLI default.** If every landed candidate is
   calibration-clean, default `--approximate-posterior-calculation`
   to `true` (matching this plan's configured default). If any
   landed candidate has a non-trivial calibration cost, flip
   the default to `false` and document the cost in the CLI
   help text.

### Why this is an evaluation, not a commitment

The biggest unknown is item 3's calibration cost. The HWE
prior's `p̂` quantisation feeds into the M-step's reweighting,
which in turn updates `p̂` — closed-loop error amplification is
possible at marginal-evidence sites. The other candidates
(items 1, 2, 4) are either exact or have no closed-loop path;
they are likely safe wins. The plan keeps the door open for all
of them via the flag, and the evaluation tells us which to
actually implement.

### File touch for the LUT path (when it lands)

- New module `src/var_calling/posterior_engine_luts.rs`
  housing the precomputed tables (constants or `OnceLock`s),
  the quantise + interpolate helpers, and unit tests for
  table-vs-exact accuracy at the budget above.
- `posterior_engine.rs` gains a thin branch on
  `config.approximate_posterior_calculation` in each candidate
  call site; the exact-math implementation stays as the
  fallback and as the calibration oracle.
- `build.rs` only if a candidate's table is large enough
  (≫ 100 kB) that compile-time generation beats first-use
  initialisation. Default to `OnceLock`-on-first-use until
  there is a reason otherwise.

## Out-of-scope follow-ups

- **Contamination machinery** — full plan in a follow-up. Pins
  the `c_s` sufficient-statistic question, the
  `ContaminationConfig` shape, and the shared-likelihood-helper
  question.
- **Rayon-parallel records.** Each record's EM is independent
  in no-contamination mode; trivial rayon-iter over the upstream.
  Add when profiling motivates it.
- **Variational-Bayes upgrade for rare-variant calibration**
  ([calling_pipeline_architecture.md:1905-1908](../specs/calling_pipeline_architecture.md#L1905-L1908)).
- **Polyploid HWE-with-`F` prior** (if v1 ships diploid-only).
- **Cohort-only compound inference**
  ([calling_pipeline_architecture.md §"Future extension"](../specs/calling_pipeline_architecture.md#L1910)).
- **Outer iteration in contamination mode**
  (`--contamination-outer-iterations`) if real cohorts show
  Pass 2's frozen `c_s` is insufficient. Implemented as
  additional scratch-file replay passes — cheap per pass since
  Stage 4/5 work never re-runs.
- **CLI wiring** — `--inbreeding`, `--ref-pseudocount`,
  `--snp-alt-pseudocount`, `--indel-alt-pseudocount`,
  `--compound-alt-pseudocount`. Ships with the cohort CLI
  subcommand.
- **`PosteriorRecord` ↔ VCF writer** — separate plan.

## File touch list

**This plan's commit** (no-contamination posterior engine):

- `src/var_calling/mod.rs` — add `pub mod posterior_engine;`.
- `src/var_calling/posterior_engine.rs` — new file:
  `PosteriorEngine`, `PosteriorRecord`, `EmDiagnostics`,
  `PosteriorEngineConfig` (with `DEFAULT_*` consts),
  `PosteriorEngineError`, `NonFiniteKind`, the EM loop,
  the HWE-with-`F` prior, the per-record M-steps on `p̂` and
  `f̂_C`, full `#[cfg(test)]` module covering the no-contamination
  unit tests above.
- `tests/posterior_engine_integration.rs` (new) — end-to-end
  tests on synthetic `MergedRecord` fixtures and on a
  small-cohort PspReader→...→PosteriorEngine pipeline.
- (No changes to `src/lib.rs` — `pub mod var_calling;` already
  exported.)

**Deferred follow-up commits** (not in this plan):

- Contamination mode (separate plan).
- VCF writer (separate plan).
- gVCF cleanup including
  [src/genotype_posteriors.rs](../../src/genotype_posteriors.rs)
  deletion — already gated by the Stage 5 plan
  ([cohort_per_group_merger.md §"Deferred cleanup"](cohort_per_group_merger.md#L1401-L1426)).
  No changes here.

## Decisions to confirm before implementation

The draft above flags these inline; collected here for the review
pass:

1. **Diploid-only for v1 or full polyploid support?** Stage 5
   supports the configured ploidy; Stage 6's HWE-with-`F` prior
   is only spec'd at diploid. Diploid-only ships faster.
2. **Plan name and module path.** Plan uses
   `posterior_engine.md` and `src/var_calling/posterior_engine.rs`
   — confirm or rename.
3. **Initialise `p̂` flat or from Dirichlet pseudocount?** GATK
   flat (first iteration), then prior kicks in from iteration 2.
   Mirror.
4. **Convergence threshold and max-iterations defaults.** Plan
   defaults: `1e-4` on `max |Δp̂|`, 50 iterations. Confirm or
   revise — GATK uses `0.1` on raw counts (depth-dependent).
5. **`PosteriorRecord` shape.** Currently carries posteriors,
   GT, GQ, QUAL, `p̂`, `f̂_C`, diagnostics. Should it also
   forward `MergedRecord.scalars` for the VCF writer's DP/AD
   fields, or should the writer join on (chrom, pos)?
6. **Non-convergence behaviour.** Error out, or emit last
   iterate with a `not_converged: true` flag? GATK and
   freebayes emit; the engine here errors. Decide.
7. **Scope of v1.** No-contamination mode only, per the draft;
   confirm contamination mode (and the choice between
   Algorithms 3 and 5) is a follow-up plan rather than a v1
   must-have.
8. **LUT-based approximate posterior calculation** — keep
   `--approximate-posterior-calculation` default `true`, or
   wait on default until the evaluation actually lands? Plan
   currently defaults to `true` on the assumption the
   evaluation confirms the wins (matches the user's intent for
   "default on"), with a documented "flip the default if
   evaluation falls short" escape hatch. See §"Approximation
   via precomputed lookup tables (evaluation)" for the full
   evaluation methodology.
