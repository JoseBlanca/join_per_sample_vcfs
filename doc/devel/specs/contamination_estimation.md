# Contamination estimation — `.psp` side-pass

**Status:** Design spec, 2026-05-17. Defines the optional pre-pass
that estimates per-sample contamination fractions `c_s` and
per-batch contamination-source allele-class distributions `q_b`
from raw `.psp` data, before the posterior engine runs. The side-
pass exists so that Stage 6 can take `c_s` and `q_b` as **frozen
inputs** rather than estimating them inside its EM loop.

The architecture-level overview lives in
[calling_pipeline_architecture.md §"Contamination estimation
side-pass"](calling_pipeline_architecture.md); this document is
the authoritative spec for the algorithm, the parameters, the
artefact format hand-off to Stage 6, and the design rationale.

For the broader history of how the design got here, see the
"algorithmic alternatives considered" section of the posterior-
engine implementation plan
([posterior_engine.md](../implementation_plans/posterior_engine.md))
— Algorithms 3, 5, and 6 there. This side-pass supersedes all
three for the in-engine path.

## Purpose

A fraction `c_s` of sample `s`'s reads at any given site are not
sample `s`'s DNA — they are contamination, typically from index
hopping within a sequencing batch. Uncorrected, contamination
biases per-sample likelihoods (false-positive het calls in
homozygous samples) and through the cohort allele-frequency
estimate `p̂` it weakly biases the whole cohort. The Stage 6
posterior engine consumes `c_s` and `q_b` to correct this via a
per-read mixture likelihood.

This side-pass produces those two parameter sets by running a
direct maximum-likelihood estimator on raw `.psp` read piles,
bypassing Stages 3, 4, and 5 entirely. It is **optional**: when
`--contamination-batches` is not supplied, no side-pass runs and
Stage 6 uses `c_s = 0` for every sample. When the user instead
supplies `--contamination-estimates` directly, the side-pass is
skipped and the user-supplied values are used.

## Shape: a side-pass, not a stage of the main pipeline

```
.psp_i ─┐
.psp_j ─┼─► [contamination side-pass] ─► (c_s vector, q_b table) ─┐
.psp_k ─┘            ▲                                              │
                     │                                              ▼
                     └── reference FASTA (optional, for AF prior)   Stage 6 init
                                                                    │
.psp_i ─┐                                                           │
.psp_j ─┼─► Stages 3 → 4 → 5 ────────────────────────────────────► Stage 6 EM ──► VCF
.psp_k ─┘
```

The side-pass and the main pipeline both consume the same `.psp`
files. They run in sequence: side-pass first, main pipeline
second, with the side-pass's output handed to Stage 6 in memory
as a config field on `PosteriorEngineConfig`. There is no on-disk
artefact between them by default — see §"Optional artefact
persistence" for the opt-in to write the estimates out.

The side-pass does *not* go through Stages 3, 4, or 5. It reads
`.psp` files at the per-position level (same multi-way merge as
Stage 3 consumes, minus the DUST filter) and operates on raw
per-(sample, position) read counts and base-quality aggregates.
No DUST filter, no grouping, no per-group merger, no per-record
EM.

## Why a side-pass rather than inside Stage 6

Three forces push the estimation out of Stage 6:

1. **Contamination parameters are global, not per-record.** `c_s`
   is a property of a sample's library, not of any single site;
   `q_b` is a property of a sequencing batch. Estimating them
   inside Stage 6's per-record EM either (a) breaks the streaming
   shape (Algorithm 2 / 4 in [posterior_engine.md](../implementation_plans/posterior_engine.md)
   — buffers the genome), (b) needs a scratch file (Algorithm 3 —
   100–300 GB at cohort scale), or (c) estimates from an in-
   engine subsample (Algorithm 6 — still pays per-record EM and
   Stage 4/5 cost on the estimation pass).
2. **Reads, not records, are the natural sufficient statistic
   for `c_s`.** The contamination model is a per-read mixture:
   each read is either own or contaminant. The site-level
   contribution is a marginalisation over read-level latent
   variables. A direct read-level MLE on `.psp` aggregates is
   simpler, faster, and statistically equivalent to (or better
   than) any path that goes through per-record EM first.
3. **This is what existing tools do.** VerifyBamID(2) and ContEst
   are both raw-read mixture estimators that bypass the variant
   caller entirely. They are the field's reference for
   contamination estimation, and the design here is a direct
   adaptation — the only differences are that `q_b` is estimated
   per-batch from the cohort itself rather than supplied from an
   external reference panel (see §"Cohort-derived vs external
   allele frequencies" for the trade-off), and that the
   estimator is shaped to drop out and pass `c_s = 0` cleanly
   when contamination is not requested.

Stage 6's EM loop loses two M-steps as a result (the existing
architecture's `q_b` and `c_s` M-steps; see [calling_pipeline_architecture.md §Stage 6](calling_pipeline_architecture.md)).
It now runs only `p` and `f_C` M-steps over the merged-record
stream, with frozen `c_s` and `q_b` baked into the E-step's
mixture likelihood.

## Algorithm — read-level mixture MLE

### Step 1 — informative-site filter (per (sample, site) pair)

A `(sample s, position i)` pair enters the estimator when **all**
of the following hold, computed from raw `.psp` aggregates with
no per-record EM:

- **Sample depth at the site** `n_{s,i} ≥ MIN_DEPTH` (default 10).
  Low-depth pairs cannot confidently call the sample's genotype,
  and contribute too little signal to be worth the bookkeeping.
- **Sample's observed major-allele fraction**
  `(max_a n_{s,i,a}) / n_{s,i} ≥ MIN_MAJOR_FRACTION` (default 0.95).
  The sample is confidently homozygous-major at this site. Any
  minor reads are candidate contamination (or base errors —
  step 3 separates them).
- **Cohort has minor-allele variation at the site.** Summed
  across all samples, the second-most-common allele has count
  `≥ MIN_COHORT_MINOR_COUNT` (default 2) and observed cohort
  minor-allele fraction `≥ MIN_COHORT_MINOR_FRACTION` (default
  0.005). Sites that are 100 %-major across the cohort carry no
  information about `c_s` — any minor reads would be pure base
  errors regardless of contamination.

The per-(sample, site) hom-major call uses raw observed
fractions, not a per-record EM estimate. Raw fractions are good
enough: the post-EM call is the same quantity up to pseudocount
shrinkage and contamination-induced bias, neither of which
matters for "is the sample confidently homozygous here?"

### Step 2 — reservoir-sample informative sites until convergence or cap

The estimator processes `.psp` files in genomic order via the
same multi-way per-position iterator Stage 3 uses (minus the
DUST filter). At each position that passes the cohort-minor-
variation check, it iterates samples and emits one (sample,
site, allele_counts, BQ_aggregates) tuple per (sample, site) pair
passing the per-sample filters.

Tuples accumulate into in-RAM per-sample and per-batch
aggregators. The estimator stops when either of:

- **Convergence**: at every `--contamination-block-size` tuples
  (default 1000), re-estimate `c_s` from all data accumulated so
  far. The largest per-sample delta `max_s |c_s_new − c_s_prev|`
  must fall below `--contamination-stability-tolerance` (default
  1e-3) for `--contamination-stability-blocks` (default 3)
  consecutive blocks.
- **Site cap**: total tuples processed reaches
  `--contamination-max-sites` (default 10 000 sites, counted as
  positions — not per-sample tuples).

If the genome ends before convergence and before the site cap,
that means there were not enough informative sites in the cohort
to estimate contamination at the requested tolerance. The pass
returns `ContaminationEstimateDidNotConverge` carrying the
final-block per-sample deltas and a recommended next-run
configuration (typically a wider tolerance or a lower
MIN_MAJOR_FRACTION). See §"Non-convergence behaviour".

### Step 3 — separate base errors from contamination signal

For each retained tuple, the expected minor-read count *under
the c_s = 0 hypothesis* is

```
E[n_minor | c_s = 0] = Σ over reads r at (s, i) of 10^(−BQ_r / 10)
                     ≈ n_{s,i} · ε̄_{s,i}
```

where `ε̄_{s,i}` is the depth-weighted per-base error rate at the
site (recovered from the `.psp` `Σ ln_BQ` scalar). Excess minor
reads above this expectation are contamination signal; subtract
the error contribution before passing to the mixture MLE so the
estimator doesn't conflate the two.

### Step 4 — joint MLE on `c_s` (per sample) and `q_b` (per batch)

Notation:
- `n_{s,i,a}` = observed read count for sample `s` at site `i`
  carrying allele class `a ∈ {REF, SNP_alt, INDEL_alt}` (matching
  the existing `q_b` allele-class enum in
  [calling_pipeline_architecture.md §Stage 6](calling_pipeline_architecture.md)).
- `g_{s,i}` = sample s's confident homozygous-major genotype at
  site `i` (decided in Step 1).
- `ε_{s,i,r}` = per-read base-error probability from BQ.
- `c_s ∈ [0, 1]` = unknown.
- `q_b ∈ Δ^{K−1}` (probability simplex, K = number of allele
  classes) = unknown.

Per-read mixture likelihood:

```
P(read r with allele a | s, i, c_s, q_{b(s)})
    = (1 − c_s) · P(a | own genotype g_{s,i}, error ε_{s,i,r})
    +     c_s   · q_{b(s)}[class(a)]
```

The first term encodes "read came from sample s with genotype
`g_{s,i}`; mismatches happen at base-error rate"; the second
encodes "read came from a contaminant whose allele frequency is
`q_b` for the batch s belongs to".

Joint MLE: maximise `Π_{(s, i, r)} P(read r | ...)` over `{c_s}`
and `{q_b}` simultaneously. Coordinate ascent converges in a few
iterations because the parameter coupling is weak:

1. **Init**: `c_s = 0.02` for every sample (rough prior on
   typical sequencing-batch contamination), `q_b` uniform over
   allele classes.
2. **Update `q_b`** (closed form, per batch):
   ```
   q_b[a] ∝ α_a + Σ over (s ∈ batch b, i, r) of c_s · 1[class(read r) = a]
   ```
   Dirichlet pseudocounts `α_a` match the Stage 6 per-position
   pseudocounts (`--ref-pseudocount`, `--snp-alt-pseudocount`,
   `--indel-alt-pseudocount`).
3. **Update `c_s`** (1D maximisation per sample over [0, 1]):
   maximise
   ```
   L_s(c_s) = Σ over (i, r) for s of
              ln[(1 − c_s) · P(r | g_{s,i}, ε_{s,i,r}) + c_s · q_{b(s)}[class(r)]]
   ```
   via golden-section search or a fixed grid at 0.001-step
   resolution. The objective is concave in `c_s` for fixed `q_b`
   so either method works; the grid is preferred for its
   reproducibility.
4. **Iterate** until `max_s |c_s_new − c_s_prev| < 1e-5` (an
   *inner* convergence threshold on the joint MLE, distinct from
   the *outer* block-stability check in Step 2). Typically 5–10
   iterations.

### Step 5 — apply small-batch and explicit-override floors

After Step 4 converges:

- **Singleton batches**: any batch containing exactly one sample
  cannot have its `c_s` and `q_b` distinguished from each other
  (the math degenerates). For samples in singleton batches,
  force `c_s = 0` and leave `q_b` unallocated for that batch.
- **Small batches**: any batch with fewer than
  `--min-batch-size-for-contamination` samples (default 5) gets
  the same treatment — `c_s = 0` for its samples, `q_b` zeroed.
  At small batch sizes the MLE is too noisy to be useful and
  the conservative thing is to skip the correction.
- **User-supplied overrides**: if `--contamination-estimates`
  was supplied, the user-provided `c_s` values for the named
  samples replace whatever the MLE produced for those samples.
  Other samples keep their MLE values. This is the escape hatch
  for cohorts where the user has external contamination
  estimates from VerifyBamID(2) or another tool.

The side-pass returns `ContaminationEstimates { c_s_per_sample,
q_b_per_batch }` to Stage 6.

## Convergence and non-convergence behaviour

### Convergence is a property of the subsample, not of the MLE

The Step 4 inner MLE always converges (concave 1D objectives,
closed-form Dirichlet update). What can fail is the **outer**
finite-sample-stability check (Step 2): the `c_s` estimates may
keep drifting as more sites are added, indicating the subsample
is too small relative to the cohort's signal-to-noise.

A stable estimate after K = 3 blocks of M = 1000 sites each
means: across the last 3 000 informative sites added, the
largest per-sample `c_s` movement was below 1e-3. That is what
the user is buying when they take the `c_s` numbers at face
value.

### When stability fails — actionable error

`ContaminationEstimateDidNotConverge` carries:

- The final per-sample deltas across the last two blocks.
- The total number of informative sites processed.
- A recommended `--contamination-max-sites` for the next run,
  computed by extrapolating the observed delta-vs-N curve to the
  N at which the largest delta would fall under tolerance.
  Formula: if the last-two-blocks delta scaled roughly as
  N^(−1/2) (the expected MLE scaling), then
  `N_recommended ≈ N_current · (delta_observed / tolerance)^2`.
  Capped at 10× the current cap so the recommendation does not
  spiral.

The user re-runs with the recommended cap (or a wider
tolerance, or a relaxed MIN_MAJOR_FRACTION to admit more
sites). **No silent fallback to `c_s = 0`** — that would
produce a confident-looking VCF whose contamination correction
is missing. Per Design principle 3 (no silent defaults), the
user explicitly opts in if they want contamination skipped.

### Convergence diagnostics

When `--contamination-diagnostics-out FILE` is set, the side-
pass writes a small TSV with one row per block:
`block_idx`, `sites_processed_so_far`, `max_c_s_delta`,
`mean_c_s_so_far`, `min_c_s_so_far`, `max_c_s_so_far`. Useful
for tuning the tolerance against real cohort data.

## What the side-pass consumes

- **All `.psp` files in the cohort.** Same files Stages 3–5
  consume; opened independently here. (Future optimisation:
  share the open file handles via a hand-off from a single
  `.psp`-opening machinery shared with Stage 3 — not in scope
  for v1 of this design.)
- **`--contamination-batches` mapping.** One line per sample
  naming its sequencing batch. Required for `q_b` to be well-
  defined.
- **No reference FASTA.** Unlike Stage 3, the side-pass does
  not need the reference at all — it operates entirely on
  per-sample observed allele counts. (If we ever need to filter
  by reference low-complexity, that could be added — but at the
  10 000-site target the marginal effect of low-complexity
  contamination is below the tolerance and not worth the cost.)

## What the side-pass does *not* do

- **Does not run Stages 3, 4, or 5.** No DUST filter, no
  grouping, no per-group merger, no per-record EM.
- **Does not modify the `.psp` files.** Read-only.
- **Does not write an on-disk artefact by default.** The
  estimate is handed to Stage 6 in memory. Opt-in to persistence
  via `--contamination-estimates-out FILE` (see below).
- **Does not call variants.** It identifies confidently
  homozygous (sample, site) pairs from raw counts; that is not
  a variant call (no posteriors, no QUAL, no GQ).
- **Does not estimate per-(sample, allele_class) contamination
  rates.** `c_s` is a single scalar per sample under the
  assumption that the contaminant is well-modelled as a mixture
  whose composition is captured by `q_b`. Per-class
  contamination would be statistically meaningful but breaks
  the existing Stage 6 mixture-likelihood shape and is out of
  scope.

## Parameters and opt-out

| Flag | Default | Effect |
|---|---|---|
| `--contamination-batches FILE` | — | **Enables the side-pass.** TSV mapping `sample_id → batch_id`. Without this flag, no side-pass runs; Stage 6 uses `c_s = 0`. |
| `--contamination-estimates FILE` | — | Skip the side-pass entirely; use the user-supplied `c_s` values. Format: TSV `sample_id → c_s`. `q_b` is still estimated by a short side-pass restricted to `q_b`-only updates (Step 4 with `c_s` frozen) unless `--contamination-source-distributions FILE` is also supplied. |
| `--contamination-source-distributions FILE` | — | Skip the side-pass entirely; use user-supplied `q_b` values. Requires `--contamination-estimates` to also be set. |
| `--contamination-estimates-out FILE` | — | Write the side-pass output to disk for inspection / reuse. TSV with per-sample `c_s` and per-batch `q_b`. |
| `--contamination-diagnostics-out FILE` | — | Write the per-block convergence trace (see §"Convergence diagnostics"). |
| `--contamination-max-sites N` | 10 000 | Hard cap on positions processed. |
| `--contamination-block-size N` | 1000 | Sites per stability-check block. |
| `--contamination-stability-tolerance T` | 1e-3 | Max per-sample `c_s` delta between consecutive blocks for convergence. |
| `--contamination-stability-blocks K` | 3 | Number of consecutive within-tolerance blocks required to declare convergence. |
| `--contamination-min-depth D` | 10 | Per-sample minimum read depth at a site to consider the (sample, site) pair. |
| `--contamination-min-major-fraction F` | 0.95 | Per-sample minimum observed major-allele fraction to call hom-major. |
| `--contamination-min-cohort-minor-count N` | 2 | Cohort-summed minimum minor-allele read count for a site to be informative. |
| `--contamination-min-cohort-minor-fraction F` | 0.005 | Cohort-summed minimum minor-allele fraction for a site to be informative. |
| `--min-batch-size-for-contamination N` | 5 | Batches below this floor get `c_s = 0` regardless. Singleton batches always get `c_s = 0` (the math degenerates). |

Following Design principle 2, every threshold is explicit. The
defaults are taken from the VerifyBamID family where they have
empirical grounding (`MIN_MAJOR_FRACTION = 0.95`, `MIN_DEPTH =
10`) and from the existing architecture's prior choices where
applicable (`--min-batch-size-for-contamination`, the Dirichlet
pseudocounts in Step 4).

## Cost and parallelism

**Memory** is `O(cohort_size + n_batches × n_allele_classes +
sites_processed × cohort_size × n_allele_classes_observed)`.
At the default 10 000 sites and a 1000-sample cohort with ~3
allele classes per kept tuple, that is ~30 MB of aggregator —
trivial.

**CPU** is dominated by the multi-way per-position scan over
`.psp` files (the same cost as Stage 3 minus DUST, except the
side-pass also stops early once convergence fires — typically
after a few thousand sites, well before whole-genome scan).
Step 4's MLE iteration runs on the small accumulated tuples and
is negligible.

**Parallelism**. The per-position scan parallelises across
chromosome partitions the same way Stage 1 parallelises across
samples. Each partition produces partial aggregator state; merge
at the end. The block-stability check requires reasoning about
global state, so it runs sequentially over the merged
aggregator after each partition completes a full block. Cheap.

For early-stop to work cleanly with partitioned scans, partition
work is dispatched in genomic-order chunks and the side-pass
joins partitions in order; once early-stop fires, in-flight
partitions are cancelled. This gives reproducible results at the
cost of slightly less parallelism than full-genome processing
would offer — fine, because the side-pass is cheap relative to
the main pipeline regardless.

## Hand-off to Stage 6

In-memory hand-off via `PosteriorEngineConfig.contamination:
Option<ContaminationEstimates>`:

```rust
pub struct ContaminationEstimates {
    /// c_s indexed by cohort sample index. None for samples in
    /// singleton or below-floor batches (treated as c_s = 0).
    pub c_s_per_sample: Vec<Option<f64>>,
    /// q_b indexed by batch index, with the batch_id <-> batch_idx
    /// mapping carried alongside. Each entry is a probability vector
    /// over allele classes (REF, SNP_alt, INDEL_alt, …).
    pub q_b_per_batch: Vec<[f64; N_ALLELE_CLASSES]>,
    /// Sample-to-batch mapping (cohort-sample-idx → batch-idx).
    pub sample_to_batch: Vec<usize>,
    /// Estimator provenance — was this MLE-fitted by the side-pass,
    /// loaded from --contamination-estimates, or both?
    pub source: ContaminationEstimateSource,
}
```

Stage 6 reads this once at engine construction and uses it in
every per-record E-step's mixture-likelihood calculation. The
side-pass's output is immutable inside Stage 6.

## Optional artefact persistence

When `--contamination-estimates-out FILE` is set, the side-pass
writes a small TSV after Step 5 completes:

```
# pop_var_caller contamination-estimates v1
# generated 2026-05-17T09:15:22Z
# cohort_sha256 = <hash of sorted sample IDs>
sample_id   batch_id   c_s
SAMPLE_001  batch_A    0.0231
SAMPLE_002  batch_A    0.0184
...

# batch contamination-source distributions
batch_id   ref_frac    snp_alt_frac    indel_alt_frac
batch_A    0.812       0.183           0.005
batch_B    ...
```

The same file is acceptable as input to a future run via
`--contamination-estimates` (read the per-sample section) plus
`--contamination-source-distributions` (read the per-batch
section). The cohort_sha256 line allows the loader to refuse
mismatched cohorts.

This is *not* the on-disk artefact that gates the side-pass on /
off. It is purely a debugging / reuse convenience.

## Cohort-derived vs external allele frequencies

Existing tools (VerifyBamID2, ContEst) estimate `q_b` against an
external reference panel (1000 Genomes, gnomAD). This pipeline
estimates `q_b` from the cohort itself — the cohort *is* the
allele-frequency reference.

The trade-off:

- **Cohort-derived `q_b` is biased by the contamination it is
  trying to estimate.** If sample s is contaminated by sample t,
  reads from t leak into s and inflate the apparent minor-allele
  count at sites where t is het and s is hom-major. The Step 4
  joint MLE accounts for this partially (the bias enters
  `q_b`'s update through `c_s`-weighted counts) but a small
  residual bias remains. The bias is `O(c̄²)` where `c̄` is the
  cohort-average contamination — well below the tolerance at
  typical `c̄ ≤ 0.05`.
- **Cohort-derived `q_b` is more accurate when the cohort is
  unrepresented in reference panels.** Non-European
  populations, plant cohorts, breeding lines — these are
  exactly the populations this pipeline targets, and external
  panels would mis-estimate `q_b` there. The cohort-derived
  approach is the right default.
- **External panels are still useful when the cohort is tiny.**
  At cohort sizes below a few dozen samples per batch, the
  cohort-derived `q_b` is noisy. A future extension would
  expose `--external-allele-frequencies FILE` to swap in a
  pre-computed `q_b` from a reference panel. Not in scope for
  v1 of the side-pass.

## Test strategy

### Unit tests

1. **Single-sample, no contamination, all-REF cohort.** No
   informative sites (cohort minor-allele check fails everywhere).
   Side-pass returns `c_s = 0` with provenance "no informative
   sites".
2. **Two-sample cohort, batch_A contains both, one sample seeded
   with synthetic 3 % contamination from the other.** Side-pass
   recovers `c_s ≈ 0.03` for the contaminated sample within
   tolerance.
3. **Singleton batch.** Sample in batch_solo gets `c_s = 0`
   regardless of input data.
4. **Below-floor batch.** Two samples in batch_small (floor = 5)
   both get `c_s = 0`.
5. **User-supplied `--contamination-estimates` overrides MLE for
   the named samples; unnamed samples still get MLE values.**
6. **Convergence diagnostics output format.** TSV header and row
   shape match the spec.

### Property tests (proptest)

1. **`c_s ∈ [0, 1]` and `q_b` is a simplex** for every cohort
   shape proptest generates.
2. **Permuting samples within a batch does not change `q_b` for
   that batch.**
3. **Adding more sites to a converged subsample does not move
   `c_s` beyond the tolerance** (i.e. the convergence criterion
   is consistent with itself).

### Integration tests

1. **End-to-end with `--contamination-batches` enabled.** A
   synthetic cohort with known per-sample contamination
   fractions; side-pass output matches ground truth within
   tolerance; Stage 6 consumes the output and emits a VCF whose
   genotype calls match the no-contamination ground-truth more
   closely than a Stage 6 run with `c_s = 0` would.
2. **Non-convergence on a deliberately underspecified cohort.**
   A 4-sample cohort with 1 batch, low-coverage everywhere. The
   side-pass returns `ContaminationEstimateDidNotConverge` with
   a sensible recommended next-run cap.
3. **Round-trip via `--contamination-estimates-out` and
   `--contamination-estimates` + `--contamination-source-
   distributions`.** Run 1 produces the artefact; Run 2 loads
   it and skips the side-pass; final VCFs match bit-for-bit.

Integration tests use the same per-sample-pileup fixtures the
rest of the pipeline uses, with synthetic contamination
injection at the `.psp` level so the ground truth is known.

## Relationship to other documents

- [calling_pipeline_architecture.md](calling_pipeline_architecture.md)
  carries the architecture-level overview (one-paragraph
  summary and the side-pass's place in the pipeline diagram)
  and the Stage 6 revisions that follow from `c_s` / `q_b`
  being frozen inputs rather than EM parameters.
- [posterior_engine.md](../implementation_plans/posterior_engine.md)
  is the implementation plan for Stage 6. Its "algorithmic
  alternatives considered" section discusses Algorithms 3, 5,
  and 6 — all three are superseded by this side-pass for the
  in-engine path. The plan's Algorithm 6 section is preserved
  as design history and points here for the authoritative
  design.
- [per_sample_pileup_format.md](per_sample_pileup_format.md)
  is the `.psp` format spec. The side-pass reads the same
  required columns Stages 3–5 read (per-allele counts, base-
  quality aggregates); it does not add any new column
  requirements.
- [design_principles.md](design_principles.md) — the no-silent-
  defaults rule (Principle 3) drives the explicit failure on
  non-convergence rather than a silent fallback to `c_s = 0`.
