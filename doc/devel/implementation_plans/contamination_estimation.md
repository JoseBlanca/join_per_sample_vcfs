# Contamination estimation side-pass — implementation plan

Proposal date: 2026-05-17.

## Domain intent

Implement the contamination side-pass specified in
[contamination_estimation.md](../specs/contamination_estimation.md):
a separate streaming pass over the per-sample `.psp` files that
estimates per-sample contamination fractions `c_s` and per-batch
contamination-source allele-class distributions `q_b` via block-
wise online EM on a read-level mixture model, and hands both to
the Stage 6 posterior engine as frozen inputs.

```
sample_A.psp ─► PspReader ─┐
sample_B.psp ─► PspReader ─┼─► PerPositionMerger ─► ContaminationEstimator ─► ContaminationEstimates
sample_C.psp ─► PspReader ─┘                                                          │
                                                                                       ▼
                                          PerGroupMerger ──► PosteriorEngine (Stage 6) [frozen c_s, q_b]
```

The side-pass and the main pipeline both consume the `.psp`
files; the side-pass runs first when `--contamination-batches`
is supplied, and its output is consumed once by Stage 6 at
engine construction. No on-disk artefact between the two stages
by default.

Pipeline-overview placement and the Stage 6 consumer-side
changes are in
[calling_pipeline_architecture.md §"Contamination estimation side-pass (optional)"](../specs/calling_pipeline_architecture.md)
and the revised
[Stage 6 algorithm](../specs/calling_pipeline_architecture.md);
the side-pass spec
([contamination_estimation.md](../specs/contamination_estimation.md))
is the authoritative source for the algorithm, parameter list,
and behavioural contract. **This plan does not re-derive the
spec.** When the plan and the spec disagree, the spec wins —
flag the divergence on this plan instead of forking the design.

## Why now

The posterior engine (Stage 6 v1) shipped on 2026-05-16 with
contamination explicitly scoped out — the contamination-mode
follow-up plan called out the side-pass design as a "leading v2
candidate" ([posterior_engine.md §"Algorithm 6"](posterior_engine.md)),
and the side-pass spec landed on 2026-05-17 with the online-EM
+ mutually-exclusive-stopping refinements. Implementing the
side-pass now is the cleanest path to a contamination-aware
end-to-end pipeline: it requires no changes to Stages 1–5, only
the minimal Stage 6 plumbing to consume frozen `c_s` / `q_b`.

The work also unblocks the cohort-CLI subcommand: that
subcommand needs to know whether contamination correction is
enabled, and the only way to make that decision concrete is to
have the side-pass behind a real flag.

## What is already in place

- **Spec.** [contamination_estimation.md](../specs/contamination_estimation.md)
  — authoritative.
- **Architecture-level wiring.** [calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md)
  shows the side-pass as a parallel branch off the `.psp` files
  feeding Stage 6; the Stage 6 section is already revised so
  `c_s` and `q_b` are frozen inputs (no in-EM M-steps on them).
- **`.psp` reader.** [src/per_sample_pileup/psp_reader.rs](../../src/per_sample_pileup/psp_reader.rs)
  (and the `psp` module surface in
  [src/per_sample_pileup/mod.rs](../../src/per_sample_pileup/mod.rs))
  provides per-sample streaming access to `.psp` records with
  the per-allele scalars and BQ aggregates the side-pass needs.
- **Multi-way per-position iterator.** [src/var_calling/per_position_merger.rs](../../src/var_calling/per_position_merger.rs)
  is the cohort-wide k-way merge that the side-pass consumes —
  same iterator Stage 3 will sit on top of. The side-pass uses
  its unfiltered output (no DUST in the path).
- **`PosteriorEngineConfig` skeleton.** [src/var_calling/posterior_engine.rs](../../src/var_calling/posterior_engine.rs)
  already carries a `contamination: Option<ContaminationConfig>`
  field (currently an empty struct that returns
  `PosteriorEngineError::ContaminationModeUnimplemented` when
  set). The wiring point for `ContaminationEstimates` is
  exactly there.
- **`MergedRecord` likelihood scalars** in
  [src/var_calling/per_group_merger.rs](../../src/var_calling/per_group_merger.rs)
  carry the per-(sample, allele) scalars Stage 6 needs to
  recompute the mixture likelihood. **Stage 5 already exposes
  the scalar projection**; what is missing is the helper that
  applies frozen `c_s, q_b` to those scalars in the Stage 6
  E-step — see §"Stage 6 consumer-side changes".

## What is *not* in place (and is being built here)

- The side-pass itself: new module
  `src/var_calling/contamination_estimation.rs`.
- The `ContaminationEstimates` hand-off type (lives in the new
  module, re-exported into `posterior_engine.rs`'s public API
  surface so config consumers see it from the same place).
- Stage 6's consumption of frozen `c_s` / `q_b` in its
  mixture-likelihood E-step. The current Stage 6 E-step calls
  the `c_s = 0` path; the contamination-aware path needs to
  apply `(1 − c_s) · L_own + c_s · L_contam` at the per-read
  level, derived from the `MergedRecord` scalars.
- CLI parser bindings for the side-pass flags. Out of scope for
  this plan (lands with the cohort subcommand); the plan ships
  library API only and configures the side-pass programmatically
  in tests.

## Module placement and file layout

New module: `src/var_calling/contamination_estimation.rs`,
exported from `src/var_calling/mod.rs` as a sibling of
`per_position_merger`, `variant_grouping`,
`per_group_merger`, `posterior_engine`.

Naming: `contamination_estimation` (matches the spec file name;
parallel to other Stage-Nth module names — verbose but
unambiguous).

Public surface re-exports the `ContaminationEstimates` type from
the new module under `var_calling::posterior_engine` as well, so
config consumers can write
`use posterior_engine::ContaminationEstimates` without crossing
module boundaries.

## API shape

```rust
// src/var_calling/contamination_estimation.rs

use crate::per_sample_pileup::psp_reader::{PspReader, PspReaderError};
use crate::var_calling::per_position_merger::{PerPositionMerger, PerPositionMergerError};

/// Allele-class enum matching the q_b indexing in the spec.
/// Mirrors the existing q_b discussion in
/// `calling_pipeline_architecture.md §Stage 6` and the spec
/// §"Step 4 — block-wise online EM".
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum AlleleClass {
    Ref = 0,
    SnpAlt = 1,
    IndelAlt = 2,
}

pub const N_ALLELE_CLASSES: usize = 3;

/// Hand-off type from the side-pass to Stage 6.
///
/// Construction sites:
/// - `estimate_contamination(...)` — the side-pass.
/// - `ContaminationEstimates::from_user_supplied(...)` —
///   `--contamination-estimates` / `--contamination-source-distributions`.
/// - `ContaminationEstimates::zero(...)` — used when
///   `--contamination-batches` is not supplied; Stage 6 takes
///   this path without entering the contamination-aware E-step.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct ContaminationEstimates {
    /// Per-cohort-sample c_s. `None` for samples in singleton
    /// or below-floor batches (Stage 6 reads this as c_s = 0).
    pub c_s_per_sample: Vec<Option<f64>>,
    /// Per-batch q_b. Indexed by batch_idx; each entry is a
    /// probability vector over allele classes.
    pub q_b_per_batch: Vec<[f64; N_ALLELE_CLASSES]>,
    /// cohort-sample-idx → batch-idx mapping (parsed from
    /// `--contamination-batches`).
    pub sample_to_batch: Vec<usize>,
    /// How this estimate was produced. Carried for provenance
    /// and for the optional TSV artefact.
    pub source: ContaminationEstimateSource,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ContaminationEstimateSource {
    /// Produced by `estimate_contamination` with the given mode.
    SidePass {
        mode: StoppingMode,
        sites_processed: u32,
    },
    /// Loaded from `--contamination-estimates` (and possibly
    /// `--contamination-source-distributions`).
    UserSupplied,
    /// A combination — user supplied c_s for some samples;
    /// side-pass filled in the rest.
    Mixed,
    /// All zeros — `--contamination-batches` was not supplied.
    Disabled,
}

/// Stopping-criterion mode for the side-pass.
///
/// Mutually exclusive at the config level; the type makes the
/// mutual exclusivity unrepresentable in code.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub enum StoppingMode {
    /// Convergence mode — run until per-block delta < tolerance
    /// for `stability_blocks` consecutive snapshots.
    Convergence {
        tolerance: f64,
        stability_blocks: u32,
    },
    /// Fixed-N mode — process exactly N informative sites.
    FixedSites { num_sites: u32 },
}

impl Default for StoppingMode {
    fn default() -> Self {
        StoppingMode::Convergence {
            tolerance: DEFAULT_STABILITY_TOLERANCE,
            stability_blocks: DEFAULT_STABILITY_BLOCKS,
        }
    }
}

#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct ContaminationEstimationConfig {
    pub stopping_mode: StoppingMode,
    pub block_size: u32,
    pub min_depth: u32,
    pub min_major_fraction: f64,
    pub min_cohort_minor_count: u32,
    pub min_cohort_minor_fraction: f64,
    pub min_batch_size_for_contamination: u32,
    pub ref_pseudocount: f64,
    pub snp_alt_pseudocount: f64,
    pub indel_alt_pseudocount: f64,
    /// Optional output-file paths for the artefact and the
    /// per-block diagnostics. `Some` only when the user passed
    /// the corresponding flag; library tests usually leave both
    /// `None`.
    pub estimates_out: Option<std::path::PathBuf>,
    pub diagnostics_out: Option<std::path::PathBuf>,
}

/// Default constants — mirror the spec's parameter-table
/// defaults. Per project convention these live as `pub const`
/// alongside the config struct.
pub const DEFAULT_STABILITY_TOLERANCE: f64 = 1e-3;
pub const DEFAULT_STABILITY_BLOCKS: u32 = 3;
pub const DEFAULT_BLOCK_SIZE: u32 = 1000;
pub const DEFAULT_MIN_DEPTH: u32 = 10;
pub const DEFAULT_MIN_MAJOR_FRACTION: f64 = 0.95;
pub const DEFAULT_MIN_COHORT_MINOR_COUNT: u32 = 2;
pub const DEFAULT_MIN_COHORT_MINOR_FRACTION: f64 = 0.005;
pub const DEFAULT_MIN_BATCH_SIZE_FOR_CONTAMINATION: u32 = 5;
// REF / SNP-alt / INDEL-alt pseudocount defaults are imported
// from the posterior engine; the side-pass uses the same values.

/// The side-pass entry point.
///
/// Library-shaped function (not an iterator): the side-pass is
/// a one-shot computation. Callers (the future cohort
/// subcommand, tests) construct readers + a batch mapping +
/// config, then await the resulting estimates.
pub fn estimate_contamination(
    psp_readers: Vec<PspReader>,
    sample_to_batch: Vec<usize>,
    n_batches: usize,
    config: ContaminationEstimationConfig,
) -> Result<ContaminationEstimates, ContaminationEstimationError>;

#[derive(thiserror::Error, Debug)]
#[non_exhaustive]
pub enum ContaminationEstimationError {
    #[error("upstream: {0}")]
    PspReader(#[from] PspReaderError),

    #[error("upstream merger: {0}")]
    Merger(#[from] PerPositionMergerError),

    /// Convergence mode, .psp exhausted without convergence.
    /// Carries the actionable recommendation the spec defines.
    #[error(
        "contamination estimate did not converge after {sites_processed} \
         informative sites; max per-sample delta in last two blocks = {max_delta} \
         (tolerance = {tolerance}). Recommended next run: {recommendation}"
    )]
    DidNotConverge {
        sites_processed: u32,
        max_delta: f64,
        tolerance: f64,
        per_sample_deltas: Vec<f64>,
        recommendation: NextRunRecommendation,
    },

    /// Fixed-N mode, .psp exhausted before N sites found.
    #[error(
        "contamination estimate found only {found} informative sites; \
         {requested} were requested. Recommendation: {recommendation}"
    )]
    InsufficientSites {
        requested: u32,
        found: u32,
        recommendation: NextRunRecommendation,
    },

    /// Sanity: e.g. batch-index out of range, or sample_to_batch
    /// length doesn't match cohort size. Signals a caller bug,
    /// not a data condition.
    #[error("bad input: {0}")]
    BadInput(String),

    /// I/O on the optional artefact / diagnostics outputs.
    #[error("I/O on output file: {0}")]
    Io(#[from] std::io::Error),
}

#[derive(Debug, Clone)]
pub struct NextRunRecommendation {
    pub message: String,
    /// Concrete suggested value when applicable
    /// (e.g. recommended `--contamination-num-sites N` in
    /// convergence-mode non-convergence; recommended lower N in
    /// fixed-N insufficient-sites).
    pub suggested_num_sites: Option<u32>,
}
```

### Errors latch

Per the merger/grouper precedent in this crate
([per_position_merger.rs:130-135](../../src/var_calling/per_position_merger.rs#L130-L135)),
once `estimate_contamination` returns an `Err` it does not
retry; the caller decides whether to surface the recommendation
to the user (CLI path) or to use it as a test fixture.

### Why a function and not an iterator

`estimate_contamination` is a one-shot computation: it consumes
the entire upstream stream up to its stopping criterion and
returns a single value. There is no downstream stage that wants
to interleave with it. The function shape matches that.

A future variant might want streaming progress callbacks for
long runs (e.g. for a CLI progress bar); we add those as a
config field (`progress: Option<Box<dyn FnMut(&ProgressEvent)>>`)
when the cohort subcommand actually needs them, not now.

## Algorithm details

### Informative-site filter (Step 1 in the spec)

For each `PerPositionPileups` item from the upstream merger:

1. **1a (cohort-wide):** sum allele counts across all samples
   present at this position. Compute observed allele fractions
   and pick the second-most-common. Reject if
   `cohort_minor_count < min_cohort_minor_count` OR
   `cohort_minor_fraction < min_cohort_minor_fraction`.
   - Implementation: a small fixed-size array of (allele,
     count) sorted descending. Allocation-free.
2. **1b (per-sample):** for each sample with a `PileupRecord`
   at this position:
   - Reject if `sum(allele_counts) < min_depth`.
   - Reject if `max_allele_count / depth < min_major_fraction`.
   - On retain: emit one `InformativeObservation { sample_idx,
     g_major: AlleleClass, allele_counts_by_class: [u32;
     N_ALLELE_CLASSES], mean_base_error: f64 }` into the
     online-EM update (Step 4 below).

The per-sample loop yields zero, one, or many observations per
qualifying position. **A position contributes to the
informative-site count only when at least one sample passed
1b** — see §"Open question — site counting" below.

### Online-EM core (Step 4 in the spec)

Two parallel state structures:

```rust
struct OnlineEmState {
    // Per-sample running sufficient statistics.
    s_per_sample: Vec<f64>,      // S_s
    n_per_sample: Vec<u32>,      // N_s

    // Per-batch running sufficient statistics.
    s_per_batch: Vec<[f64; N_ALLELE_CLASSES]>,  // S_b[a]
    n_per_batch: Vec<u32>,                       // N_b

    // Frozen parameter snapshots used for this block's γ
    // computations. Refreshed at block boundaries.
    c_s_frozen: Vec<f64>,
    q_b_frozen: Vec<[f64; N_ALLELE_CLASSES]>,

    // Stability tracking (convergence mode only).
    c_s_last_snapshot: Vec<f64>,
    consecutive_within_tolerance: u32,

    sites_processed: u32,
}
```

Per-observation update:

```rust
fn process_observation(state: &mut OnlineEmState, obs: &InformativeObservation, ...) {
    let batch = sample_to_batch[obs.sample_idx];
    for class in [AlleleClass::Ref, AlleleClass::SnpAlt, AlleleClass::IndelAlt] {
        let n_reads = obs.allele_counts_by_class[class as usize];
        if n_reads == 0 { continue; }
        // Per-read mixture likelihood:
        let p_own = own_genotype_probability(class, obs.g_major, obs.mean_base_error);
        let p_contam = state.q_b_frozen[batch][class as usize];
        let c = state.c_s_frozen[obs.sample_idx];
        // γ for each read at this (sample, site, class) is the same,
        // so multiply by n_reads.
        let γ = (c * p_contam) / ((1.0 - c) * p_own + c * p_contam);
        state.s_per_sample[obs.sample_idx] += γ * (n_reads as f64);
        state.n_per_sample[obs.sample_idx] += n_reads;
        state.s_per_batch[batch][class as usize] += γ * (n_reads as f64);
        state.n_per_batch[batch] += n_reads;
    }
}
```

`own_genotype_probability` is small and pure:
- If `class == g_major`: returns `1.0 − mean_base_error`.
- Else: returns `mean_base_error / (N_ALLELE_CLASSES − 1)` (uniform
  distribution of error reads over non-genotype classes).
  This matches the spec Step 3's `ε̄_{s,i}` use; the divisor
  reflects that base errors are not biased toward any
  particular non-genotype class.

`mean_base_error` per observation is derived from the
`.psp`-stored `Σ max(ln_BQ, ln_MQ)` aggregate per sample at
this site: `mean_base_error = exp(Σ_lnBQ / n_reads_at_site)`.

At block boundary (every `block_size` *informative* sites; see
§"Open question — site counting"):

1. Refresh parameters from running sufficient statistics:
   ```rust
   for s in 0..n_samples {
       state.c_s_frozen[s] = if state.n_per_sample[s] > 0 {
           state.s_per_sample[s] / state.n_per_sample[s] as f64
       } else {
           CS_INIT  // 0.02 — no data yet for this sample
       };
   }
   for b in 0..n_batches {
       let alpha = [ref_pseudocount, snp_alt_pseudocount, indel_alt_pseudocount];
       let alpha_sum: f64 = alpha.iter().sum();
       let s_sum: f64 = state.s_per_batch[b].iter().sum();
       let denom = alpha_sum + s_sum;
       for a in 0..N_ALLELE_CLASSES {
           state.q_b_frozen[b][a] = (alpha[a] + state.s_per_batch[b][a]) / denom;
       }
   }
   ```
2. **Convergence mode only**: compute
   `max_delta = max_s (state.c_s_frozen[s] − state.c_s_last_snapshot[s]).abs()`;
   bump or reset `consecutive_within_tolerance`; if it reaches
   `stability_blocks`, the side-pass returns success. Copy
   `c_s_frozen` to `c_s_last_snapshot`.
3. **Fixed-N mode**: if `sites_processed == num_sites`, return
   success.
4. If `--contamination-diagnostics-out` is set, append one row
   to the TSV.

### Stopping-criterion machinery (Step 2 in the spec)

The two stopping modes share the per-observation update path
and the block-boundary refresh; they diverge only in the
"should we stop?" test:

- **Convergence mode** stops on
  `consecutive_within_tolerance >= stability_blocks`. Returns
  `DidNotConverge` on stream exhaustion.
- **Fixed-N mode** stops on `sites_processed == num_sites`.
  Returns `InsufficientSites` on stream exhaustion.

The CLI mutex (the type-level `StoppingMode` enum makes it
impossible to construct a config with both) is the
implementation of the "passing both is a CLI parse-time error"
spec requirement; the actual CLI parser (out of scope for this
plan) maps the two flags into the enum and errors at parse time
if both are set.

### Finalisation (Step 5 in the spec)

Once the loop returns success:

1. For singleton batches and batches below
   `min_batch_size_for_contamination`, set `c_s_per_sample[s]
   = None` (Stage 6 reads `None` as `c_s = 0`); zero out
   `q_b_per_batch[b]`.
2. Apply `--contamination-estimates` user overrides (if any
   were supplied alongside this run — see §"Hand-off to Stage
   6" in the spec).
3. Build `ContaminationEstimates { c_s_per_sample,
   q_b_per_batch, sample_to_batch, source }` and (if
   `estimates_out` is set) write the TSV.

### Open question — site counting

The spec says block_size and num_sites are counted in
*informative sites* (positions, not per-sample tuples). The
edge case: a position passes Step 1a (cohort-wide) but no
sample passes Step 1b. Two interpretations:

- **(A)** Count it: every Step-1a-passing position ticks the
  counter, even if zero samples contributed.
- **(B)** Count only positions where at least one sample
  contributed.

(B) matches the spec's word "informative" more closely — those
positions actually advanced the estimate — and gives more
predictable convergence behaviour at edge-case cohorts where
many positions pass 1a but few pass 1b. **This plan ships with
(B);** flag if real cohort data shows a reason to switch.

## Configurable parameters

All parameters come from the spec's Parameters table; the plan
adds no new tunables. The defaults mirror the spec table's
defaults (1e-3 tolerance, 3 stability-blocks, 1000 block size,
depth 10, major-fraction 0.95, cohort-minor-count 2,
cohort-minor-fraction 0.005, batch floor 5). The `--ref-`,
`--snp-alt-`, `--indel-alt-pseudocount` flags share their
defaults with the existing Stage 6 constants
(`DEFAULT_REF_PSEUDOCOUNT = 10.0`, etc. — see
[posterior_engine.rs:59-105](../../src/var_calling/posterior_engine.rs#L59-L105)).

## Stage 6 consumer-side changes

The side-pass produces frozen `c_s` and `q_b`; Stage 6 has to
consume them in its mixture-likelihood E-step. The current
engine ships no-contamination only and returns
`PosteriorEngineError::ContaminationModeUnimplemented` when
`PosteriorEngineConfig.contamination` is set. The minimal
changes:

1. **`PosteriorEngineConfig.contamination` type change.**
   Replace the placeholder `Option<ContaminationConfig>` with
   `Option<ContaminationEstimates>`. The empty struct goes
   away.
2. **E-step mixture-likelihood path.** When `contamination` is
   `Some(estimates)`, the per-(sample, genotype) likelihood
   computation in Stage 6 evaluates
   `(1 − c_s) · L_own(reads | G) + c_s · L_contam(reads | q_b)`
   per read, where `L_own` is reconstructed from
   `MergedRecord.scalars` exactly as the `c_s = 0` path does
   today, and `L_contam` uses the per-(sample, allele-class)
   summary from `MergedRecord` evaluated against
   `q_{batch(s)}`. The shared helper for the per-read
   reconstruction is the Stage 5 scalar-projection function —
   see [per_group_merger.rs](../../src/var_calling/per_group_merger.rs).
3. **Singleton-batch / below-floor handling.**
   `ContaminationEstimates.c_s_per_sample[s] == None` is read
   as `c_s = 0` for that sample; the E-step for that sample
   collapses to the own-DNA path. No special-casing beyond
   "if None, take the existing c_s = 0 branch."

These changes are smaller in code than the side-pass itself;
they live in `posterior_engine.rs`. They are tested through
the side-pass integration tests (a sample with known c_s
flowing through both stages) and through new direct unit tests
on the E-step with hand-built `ContaminationEstimates`.

The `ContaminationModeUnimplemented` error variant is removed
in the same change.

## Out of scope (v1)

Lifted from the spec's "What the side-pass does not do" and
the broader §"Cohort-derived vs external allele frequencies":

- **Parallelism across chromosome partitions.** The spec
  describes the cancellation protocol for partitioned scans;
  v1 of this plan ships sequential. The cost analysis in the
  spec already justifies this — the side-pass is cheap
  relative to the main pipeline regardless. Add parallelism if
  real cohort runs show it matters.
- **`--external-allele-frequencies`** for swapping in a
  reference-panel `q_b`. v2 feature; see spec §"Cohort-derived
  vs external allele frequencies".
- **The `--contamination-estimates`-only short path** that
  estimates `q_b` while freezing `c_s` at user-supplied values.
  Half-built: the same machinery can do it (run the loop with
  `c_s_frozen` never refreshed and only `q_b` updating), but
  the configuration surface and the partial-finalisation logic
  add real complexity. Defer until a user asks.
- **CLI parser bindings.** The cohort CLI subcommand lands
  separately; until then the side-pass is tested through the
  library API directly. The mutual-exclusivity check is
  *type-level* in this plan (the `StoppingMode` enum); the CLI
  parser maps the two flags into the enum at parse time.
- **Per-(sample, allele-class) contamination rates.** Spec
  scope-out; would break the existing Stage 6 mixture
  likelihood shape.

## Test strategy

Unit and property tests live in
`src/var_calling/contamination_estimation.rs`'s `#[cfg(test)]`
module; integration tests in
`tests/contamination_estimation_integration.rs` (new file).

### Fixture builders

```rust
// In the same test module, alongside the existing var_calling
// test helpers in per_group_merger.rs / posterior_engine.rs.

fn pileup_record_simple(
    chrom_id: u32,
    pos: u32,
    sample_idx: usize,
    alleles: Vec<&[u8]>,
    counts_per_allele: Vec<u32>,
    mean_bq: f64,
) -> PileupRecord;

fn cohort_observations_at_position(
    chrom_id: u32,
    pos: u32,
    per_sample: Vec<(usize, Vec<u32>, f64)>,  // (sample_idx, counts, mean_bq)
) -> PerPositionPileups;

/// Build a complete synthetic per-position stream with
/// configurable cohort polymorphism and per-sample noise.
fn synth_cohort_stream(
    n_samples: usize,
    n_sites: u32,
    cohort_maf: f64,
    contamination_per_sample: Vec<f64>,
    sample_to_batch: Vec<usize>,
) -> Vec<PerPositionPileups>;
```

`synth_cohort_stream` is the main workhorse — it simulates the
read-level mixture model so tests have ground truth. The
ContaminationEstimator is tested against `synth_cohort_stream`
output with known `c_s` and recovered estimates compared to
ground truth.

### Unit tests

1. **Step-1a filter rejects monomorphic cohort.** All samples
   carry only REF at a site → site is dropped (no observations
   emitted).
2. **Step-1a filter accepts at exactly the cohort-minor
   thresholds** (boundary test).
3. **Step-1b filter rejects low-depth (sample, site) pairs.**
4. **Step-1b filter rejects ambiguous-genotype (sample, site)
   pairs** (major fraction just below 0.95).
5. **Step-1b filter accepts confidently-hom-major
   observations.**
6. **Online-EM single-batch, single-sample, ground-truth
   c_s = 0.** Side-pass run on synthetic uncontaminated stream
   returns `c_s ≈ 0` within tolerance.
7. **Online-EM single-batch, two samples, one contaminated.**
   Ground-truth c_s = 0.03 for sample 1 → side-pass recovers
   within 1e-3 tolerance after convergence.
8. **Online-EM convergence in convergence mode.** Same as #7,
   but assert on `sites_processed` reasonably bounded
   (e.g. < 5000 sites for default tolerance).
9. **Fixed-N mode hits exactly N sites.** Construct a stream
   with abundant informative sites; pass
   `StoppingMode::FixedSites { num_sites: 2000 }`; assert
   returned `ContaminationEstimateSource::SidePass { mode,
   sites_processed: 2000 }`.
10. **Type-level mutual exclusivity.** A compile-fail test
    (using `compile_fail` doc-comment) demonstrating that
    `StoppingMode` cannot be constructed in both modes
    simultaneously.
11. **Singleton-batch finalisation.** A two-sample cohort with
    each sample in its own batch → both samples get `c_s =
    None`; both batches' `q_b` zeroed.
12. **Below-floor batch finalisation.** Two samples in
    batch_small (floor = 5) → both `c_s = None`.
13. **Estimate-source provenance.** `ContaminationEstimates.
    source` is `SidePass {...}` for an unmodified side-pass run,
    `UserSupplied` when constructed from a user file, `Mixed`
    when both contribute, `Disabled` for the zero-init
    constructor.
14. **`q_b` is a probability simplex** (sum to 1, all
    non-negative) on every successful run.
15. **`c_s ∈ [0, 1]`** on every successful run.

### Property tests (proptest, 64 cases each)

1. **`c_s ∈ [0, 1]`** for every cohort proptest generates.
2. **`q_b` is a simplex per batch.**
3. **Permuting samples within a batch leaves `q_b[batch]`
   unchanged.**
4. **Block-size independence within tolerance.** For a
   converged cohort, running the side-pass at `block_size`
   ∈ {500, 1000, 2000} yields per-sample `c_s` agreeing within
   `tolerance` of each other. Validates the block-wise online
   EM design — block size should change *how often* stability
   is checked, not the estimator's limit.

### Integration tests

In `tests/contamination_estimation_integration.rs`:

1. **End-to-end, convergence mode, synthetic cohort.** Build
   `.psp` fixtures with synthetic contamination injection;
   open with real `PspReader`; run side-pass through
   `estimate_contamination`; recovered c_s within tolerance of
   ground truth.
2. **End-to-end, fixed-N mode, same cohort.** With
   `num_sites = 5000`; c_s within tolerance of mode #1 result.
3. **Non-convergence error path (convergence mode).** Tiny
   cohort (4 samples, 1 batch), low coverage everywhere →
   `DidNotConverge` returned with a non-empty `NextRunRecommendation`.
4. **Insufficient sites error path (fixed-N mode).** Same
   cohort, `num_sites = 50_000` → `InsufficientSites { requested:
   50_000, found: M, .. }`.
5. **Side-pass → Stage 6 round-trip.** Synthetic cohort with
   known c_s; run side-pass; pass output into Stage 6 via
   `PosteriorEngineConfig.contamination`; assert Stage 6's
   genotype calls match the no-contamination ground-truth more
   closely than a Stage 6 run with `ContaminationEstimates::
   zero(...)` would. This is the test that proves the
   end-to-end contamination correction works.
6. **Round-trip via `--contamination-estimates-out` and
   `--contamination-estimates`.** Run 1 writes the artefact;
   Run 2 loads it via `ContaminationEstimates::from_user_supplied`
   and skips the side-pass; final VCFs match bit-for-bit.
   (Defers the TSV format to the spec — this test only checks
   the round-trip property.)

### What integration tests *don't* cover in v1

- Real-cohort data. Synthetic injection is the v1 ground-truth
  source; a real cohort with externally-measured (e.g.
  VerifyBamID2) c_s values is a v1.5 validation activity, not
  a test-suite gate.
- The CLI parser. Tested when the cohort subcommand lands.
- Parallel partition cancellation. Sequential v1 has no
  cancellation path.

## Validation

Beyond the tests:

- **Manual run on a known-contaminated public cohort.** After
  v1 lands, run the side-pass on (e.g.) a 1000-Genomes subset
  known to have a few mis-pooled samples and compare against
  the VerifyBamID2 numbers published for that cohort. Target:
  per-sample `c_s` agreement within ~1 % absolute error. This
  goes in the implementation report, not the test suite.
- **Stability vs block-size on a real cohort.** Run at
  `block_size` ∈ {500, 1000, 5000} on the same input; confirm
  the per-sample c_s estimates agree within `tolerance`.
- **Memory profile.** Confirm steady-state memory is in the
  expected `O(n_samples + n_batches × n_classes)` range
  (10s of KB at 1000 samples × 10 batches).

## Assumptions / silent choices

- **AlleleClass discriminator on observed reads.** A read's
  class is determined by the allele's mutation type at that
  site (REF if matches reference, SNP_alt if a single-base
  substitution, INDEL_alt otherwise). The `MergedRecord` /
  `PileupRecord` types already carry the per-allele type tag;
  the side-pass reads that tag rather than re-deriving from
  reference comparison.
- **Read base-error rate is uniform across non-genotype
  classes.** When a base error happens, the spec model says
  the erroneous base is uniformly distributed over the
  non-true-genotype possibilities. This is the standard
  assumption used by Stage 5's likelihood reconstruction; the
  side-pass shares it.
- **`mean_base_error` uses the `Σ max(ln_BQ, ln_MQ)` aggregate
  divided by depth.** This is the same per-position summary
  Stage 5 uses; reasonable because the side-pass needs only a
  scalar error rate, not per-read fidelity.
- **CS_INIT = 0.02** for samples with no observations yet.
  Spec value; matches typical sequencing-batch contamination
  prior.
- **Numerical safety on γ denominator.** The mixture
  denominator `(1 − c_s) · p_own + c_s · p_contam` can be
  tiny when both `c_s ≈ 0` and the read's allele has `p_own
  = mean_base_error` (small) AND the contam allele has
  `q_b[class] = 0` (rare but possible at uniform init for
  unobserved classes). We floor the denominator at `f64::MIN_POSITIVE`
  to avoid divide-by-zero; in practice the Dirichlet
  pseudocounts on `q_b` keep it well clear of zero.
- **`PspReader` panics propagate as
  `ContaminationEstimationError::PspReader(...)`.** No special
  retry / recovery. Same pattern as the merger.

## Risks

1. **Bias from early-block γ values.** Online EM has a small
   bias from γ values computed against init parameters. The
   spec bounds this at `O(c̄² / N)`; v1 ships with no
   correction. **Mitigation:** the block-size-independence
   proptest and the side-pass → Stage 6 integration test
   together detect any bias material to calibration. If real
   cohorts show drift > tolerance after convergence, a
   one-block warmup pass that throws away the first block's
   contributions to the running sums is a five-line fix.
2. **Stage 6 mixture-likelihood implementation mistakes.**
   The E-step change touches the hot per-record EM loop in
   `posterior_engine.rs`. **Mitigation:** the side-pass →
   Stage 6 integration test has hand-built ground-truth
   `MergedRecord` fixtures that exercise the mixture path
   directly; pre-existing unit tests cover the `c_s = 0` path
   regression.
3. **Performance on the per-position scan.** The side-pass
   walks every cohort position past the cohort-MAF filter,
   which is one full `.psp` traversal across all samples. On
   a 1000-sample cohort this is non-trivial wall time. **Not
   mitigated in v1** beyond the early-stop in convergence
   mode; the cost-and-parallelism section of the spec calls
   out partition-parallel as the v2 follow-up if real runs are
   too slow.
4. **PspReader open-cost at large N.** Opening 1000+
   `PspReader`s is an N-times open-file cost the existing
   pipeline already pays in Stage 3 onwards. **No new risk;**
   the side-pass uses the same machinery.

## Implementation order

Land in small, individually-testable PRs. Each phase has
working tests; nothing is half-built across phases.

**Phase 0 — Stage 6 plumbing for frozen inputs (smallest,
unblocks everything).**
- Replace `PosteriorEngineConfig.contamination:
  Option<ContaminationConfig>` with `Option<ContaminationEstimates>`.
- Add `ContaminationEstimates::zero(n_samples, n_batches,
  sample_to_batch)` constructor.
- Stage 6 E-step learns to branch: `None` → existing c_s = 0
  path; `Some(estimates)` → new mixture path (call into
  unimplemented! for now, just to confirm the wiring compiles
  and the `c_s = 0` regression is intact).
- All existing posterior-engine tests pass; new test:
  passing `ContaminationEstimates::zero(...)` produces
  identical output to passing `None`.

**Phase 1 — Stage 6 mixture-likelihood E-step.**
- Implement the contamination-aware E-step branch.
- Unit tests on the E-step with hand-built `MergedRecord`
  fixtures and known `ContaminationEstimates`.
- Direct ground-truth: at `c_s = 0`, mixture path returns
  same result as own-only path.

**Phase 2 — `contamination_estimation` module skeleton.**
- New file with the public types (`AlleleClass`,
  `ContaminationEstimates`, `ContaminationEstimationConfig`,
  `StoppingMode`, `ContaminationEstimationError`).
- The function signature for `estimate_contamination`,
  returning `unimplemented!()`.
- Compile-fail test for the `StoppingMode` mutual exclusivity.

**Phase 3 — informative-site filter.**
- Implement `step_1a_cohort_filter` and `step_1b_sample_filter`
  as pure functions on `PerPositionPileups` and `PileupRecord`.
- Unit tests #1–#5 from the test list.

**Phase 4 — online-EM core, convergence mode.**
- Implement `OnlineEmState`, the per-observation update, the
  block-boundary refresh, the convergence stability check.
- Unit tests #6–#8 and integration test #1.

**Phase 5 — fixed-N mode and error paths.**
- Implement the `FixedSites` branch.
- Implement `DidNotConverge` and `InsufficientSites` with
  populated `NextRunRecommendation`.
- Unit tests #9 and integration tests #2–#4.

**Phase 6 — finalisation (Step 5).**
- Singleton-batch and below-floor handling.
- `--contamination-estimates` user-override merge.
- Unit tests #11–#13.

**Phase 7 — proptests.**
- Land the four proptests; tune seed counts to match the
  project's `64 cases each` precedent.

**Phase 8 — end-to-end with Stage 6.**
- Integration test #5 (side-pass → Stage 6 round-trip).
- Integration test #6 (TSV round-trip via
  `--contamination-estimates-out`).
- Implementation report write-up: side-pass + Stage 6
  consumer changes together.

Optional follow-ups (separate PRs, not gating v1):
- Diagnostics TSV output.
- Artefact TSV writer.
- CLI parser binding when the cohort subcommand lands.

## Out-of-scope follow-ups

- **Parallel partition scans** — spec §"Cost and parallelism"
  has the cancellation protocol. Add when bench shows wall
  time matters.
- **`--external-allele-frequencies`** for swapping in a
  reference-panel `q_b`.
- **`--contamination-estimates`-only mode** that estimates
  `q_b` while freezing `c_s` at user-supplied values.
- **Streaming progress callback** for long CLI runs.
- **One-block warmup pass** to remove the early-block γ bias.
  Land only if calibration testing on real cohorts shows the
  bias matters.

## File touch list

New files:
- `src/var_calling/contamination_estimation.rs`
- `tests/contamination_estimation_integration.rs`

Modified files:
- `src/var_calling/mod.rs` — `pub mod contamination_estimation;`
- `src/var_calling/posterior_engine.rs`
  - `PosteriorEngineConfig.contamination` type change
  - E-step mixture-likelihood branch
  - Remove `ContaminationModeUnimplemented` error variant
  - Re-export `ContaminationEstimates`
- `PROJECT_STATUS.md` — new entry for the side-pass under
  Stage 6's group, marking the contamination-mode work as
  shipped.

Documentation:
- Implementation report at
  `doc/devel/reports/implementations/contamination_estimation_<DATE>.md`
  after Phase 8.

## Decisions to confirm before implementation

1. **Site-counting interpretation (A) vs (B)** — see §"Open
   question — site counting". Plan ships (B); confirm before
   Phase 4.
2. **AlleleClass derivation source** — confirm
   `MergedRecord` / `PileupRecord` carries the per-allele
   `AlleleType` tag (SNP / INDEL) the side-pass needs. If
   not, derive from reference comparison in Step 1b.
3. **`mean_base_error` aggregate** — confirm the `Σ
   max(ln_BQ, ln_MQ)` summary in `.psp` is exposed by
   `PspReader` at the per-(sample, site) granularity the
   filter needs. If the aggregate is per-allele (not
   per-site), the per-site mean is a weighted sum across
   alleles.
4. **`PspReader` ownership in `estimate_contamination`.**
   Takes `Vec<PspReader>` by value (consumed) vs
   `&mut [PspReader]` (borrowed). Default to by-value: the
   side-pass is one-shot, and the caller almost certainly
   discards the readers afterward. Confirm against the future
   cohort-subcommand call site once it exists.

If any of these turns out wrong during implementation, treat
this plan as the design constraint and the discovered reality
as the bug — push back to the spec author rather than silently
adapting.
