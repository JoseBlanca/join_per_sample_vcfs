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

The filter operates at **two levels** and a `(sample s, position
i)` pair enters the estimator only when **both** sets of
conditions hold. All conditions are evaluated from raw `.psp`
aggregates with no per-record EM.

The two levels point in opposite directions on the
major-allele-fraction axis, and that is deliberate — see below.

#### 1a. Cohort-wide site filter (applied once per site)

The site qualifies if the **cohort-aggregated major-allele
fraction is LOW enough** that the site is polymorphic. Concretely,
summed across all samples in the cohort:

- The second-most-common allele has count
  `≥ MIN_COHORT_MINOR_COUNT` (default 2).
- The observed cohort-minor-allele fraction is
  `≥ MIN_COHORT_MINOR_FRACTION` (default 0.005) — equivalently,
  the cohort major-allele fraction is `≤ 0.995`.

Sites that fail (i.e. the cohort is essentially monomorphic at
that position) carry no information about `c_s`: at a site where
the whole cohort shows only allele A, contamination from another
cohort sample also brings allele A reads, so contamination is
undetectable. The few reads carrying a different allele at such
sites are pure base errors regardless of contamination.

This is the "low major-allele frequency" criterion in the
broader Algorithm 6 discussion ([posterior_engine.md](../implementation_plans/posterior_engine.md)),
phrased here as a positive lower bound on the cohort *minor*
fraction.

#### 1b. Per-sample filter (applied to each sample at a qualifying site)

For each sample `s` at a site that already passed 1a, the
sample's observation qualifies if the **sample's own
major-allele fraction is HIGH enough** to call the sample
confidently homozygous-major:

- `n_{s,i} ≥ MIN_DEPTH` (default 10). Low-depth pairs cannot
  confidently call the sample's genotype.
- `(max_a n_{s,i,a}) / n_{s,i} ≥ MIN_MAJOR_FRACTION`
  (default 0.95). At least 95 % of this sample's reads at this
  site carry the same allele — the sample is unambiguously
  homozygous-major.

Samples that fail (low depth, or genuinely heterozygous, or
ambiguous) are skipped for this site, but other samples at the
same site may still contribute. The site itself is not
discarded; only the (sample, site) tuples that don't meet
1b are.

#### Why the two filters point in opposite directions

- The **cohort** filter wants polymorphism (low cohort major
  fraction) so that contamination from another cohort sample is
  likely to introduce a *different* allele than sample s's own
  genotype — otherwise contamination is invisible.
- The **per-sample** filter wants confidence in the sample's
  genotype (high per-sample major fraction) so that minor reads
  in this sample's pile are interpretable as "must be
  contamination or base error." If the sample were genuinely
  heterozygous, minor reads would just be the second allele and
  carry no contamination signal.

A useful (sample, site) tuple is therefore one where **the
cohort is polymorphic but this particular sample happens to be
homozygous-major** — exactly the regime VerifyBamID and ContEst
exploit.

#### Why raw fractions are enough

The per-sample hom-major call uses raw observed fractions, not a
per-record EM estimate. Raw fractions are good enough: the
post-EM call is the same quantity up to pseudocount shrinkage
and contamination-induced bias, neither of which matters for
"is the sample confidently homozygous here?" at the 0.95 threshold.

### Step 2 — stream informative sites; stop on convergence or fixed-N target

The estimator processes `.psp` files in genomic order via the
same multi-way per-position iterator Stage 3 uses (minus the
DUST filter). At each position that passes the cohort-wide
filter (Step 1a), it iterates samples and feeds each (sample,
site) pair passing the per-sample filter (Step 1b) through the
online-EM update in Step 4.

**No tuples are stored.** All per-(sample, site) contributions
fold into a small set of running sufficient statistics whose
size is `O(n_samples + n_batches × n_allele_classes)` — KB at
realistic cohort sizes. See Step 4 for the precise state.

The estimator stops by exactly one of two stopping criteria,
**selected at CLI parse time and mutually exclusive**:

#### Convergence mode (default)

Active when `--contamination-stability-tolerance T` is set
(default `1e-3`); silently active when neither stopping flag is
passed.

At every `--contamination-block-size` informative sites
(default 1000), snapshot the current `c_s` vector derived from
the running sufficient statistics. Convergence is declared when
the largest per-sample delta `max_s |c_s_new − c_s_prev|` stays
below `T` for `--contamination-stability-blocks` consecutive
snapshots (default 3).

If the `.psp` files are exhausted before convergence fires, the
side-pass returns
`ContaminationEstimateDidNotConverge` carrying the final-block
per-sample deltas, the number of informative sites processed,
and a recommended next-run configuration (typically a wider
tolerance, a relaxed `MIN_MAJOR_FRACTION`, or a switch to
fixed-N mode with `N` set to the observed count). See
§"Non-convergence behaviour".

#### Fixed-N mode

Active when `--contamination-num-sites N` is set.

The side-pass processes exactly `N` informative positions (i.e.
positions passing Step 1a, regardless of how many samples
contribute at each) and returns whatever estimate the running
sufficient statistics yield at that boundary. **No convergence
check, no stability tolerance.** The user has explicitly
declared "use these N sites, give me the resulting estimate" —
the mode is reproducible by construction, which makes it the
right tool for regression tests and cross-version benchmarks.

If the `.psp` files are exhausted before `N` informative sites
are found, the side-pass returns `ContaminationInsufficientSites
{ requested: N, found: M, recommendation: ... }`. See
§"Non-convergence behaviour".

#### Mutual exclusivity

Setting both `--contamination-stability-tolerance` and
`--contamination-num-sites` is a CLI parse-time error. The user
picks one stopping discipline; the spec does not silently
combine them. This is the only CLI shape that satisfies Design
principle 3 (no silent defaults) — if both were allowed and
both fired, the user would not know which one decided the run.

Convergence mode is the default for users who have no informed
opinion on `N`; fixed-N mode is the explicit opt-in for
reproducibility.

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

### Step 4 — block-wise online EM on `c_s` and `q_b`

The maximum-likelihood estimator runs as **block-wise online
EM** over the stream of reads produced by Steps 1 and 2. State
is kept entirely in compact running sufficient statistics — **no
tuples and no per-block read buffer are retained**.

#### Notation

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

#### Per-read mixture likelihood

```
P(read r with allele a | s, i, c_s, q_{b(s)})
    = (1 − c_s) · P(a | own genotype g_{s,i}, error ε_{s,i,r})
    +     c_s   · q_{b(s)}[class(a)]
```

The first term encodes "read came from sample s with genotype
`g_{s,i}`; mismatches happen at base-error rate"; the second
encodes "read came from a contaminant whose allele frequency is
`q_b` for the batch s belongs to".

#### Running sufficient statistics

Per sample `s`:
- `S_s = Σ over reads r for s of γ_r`
- `N_s = Σ over reads for s of 1`

Per batch `b`, per allele class `a`:
- `S_b[a] = Σ over reads r in b of γ_r · 1[class(r) = a]`
- `N_b = Σ over reads in b of 1`

Where `γ_r = P(read r came from contamination | currently
frozen c_s, q_b)`. Total state: `O(n_samples × 2 + n_batches ×
(K + 1))` scalars — a few KB at typical cohort sizes.

#### Init

- `c_s_0 = 0.02` for every sample (rough prior on typical
  sequencing-batch contamination).
- `q_b_0` uniform across allele classes (renormalised through
  the Dirichlet pseudocounts `α_a` from the per-position
  pseudocounts: `--ref-pseudocount`, `--snp-alt-pseudocount`,
  `--indel-alt-pseudocount`).
- All running sufficient statistics zeroed.
- The "frozen" parameters used to compute the first block's
  `γ_r` values = the init values.

#### Block update (fires every `--contamination-block-size` sites)

The block is a **heartbeat**, not a buffer: reads are processed
one at a time as they arrive, but the parameters used to compute
`γ_r` stay frozen for the duration of the block, only refreshing
at block boundaries.

For each incoming `(sample s, site i, read r)`:

1. **E-step (γ at frozen parameters):**
   ```
   γ_r = c_s_frozen · q_{b(s)_frozen}[class(r)]
       / [(1 − c_s_frozen) · P(r | g_{s,i}, ε_{s,i,r})
          + c_s_frozen · q_{b(s)_frozen}[class(r)]]
   ```
2. **Sufficient-statistic update:** add `γ_r` into the running
   sums for `s` and for `b(s)`. Increment the matching `N`
   counters.

At end of every block (every `--contamination-block-size`
informative sites):

3. **Refresh parameter estimates** from the updated running
   sufficient statistics:
   ```
   c_s = S_s / N_s
   q_b[a] = (α_a + S_b[a]) / (Σ_a α_a + Σ_a S_b[a])
   ```
   The Dirichlet pseudocounts `α_a` appear in the `q_b`
   normalisation so `q_b` is well-defined even on cold-start
   batches; `c_s` has no pseudocount (it is a fraction, not a
   distribution).
4. **Snapshot `c_s`** for the Step-2 stability check.
5. **Re-freeze** `c_s_frozen ← c_s, q_b_frozen ← q_b` for the
   next block's `γ` computations.

No buffer is allocated, cleared, or grown between blocks. The
running sufficient statistics are cumulative across the entire
side-pass.

#### Why block-wise (and not pure per-read online EM)

A strict per-read variant would refresh the parameters after
each read, so every `γ_r` is computed against parameters that
have just shifted by one read's worth of evidence. That is the
classical online-EM setting that requires a stepwise learning
rate (Cappé & Moulines 2009, "Online EM") to converge cleanly.

Block-wise refresh keeps `γ` values consistent with a single
parameter snapshot per block, removes the need for any learning-
rate schedule, and aligns naturally with the stability snapshots
Step 2 needs anyway. The `--contamination-block-size` knob is
the same heartbeat for both EM refresh and stability check.

#### Why this converges to the offline MLE

As the side-pass advances, the running sufficient statistics
converge to their population expectations under the stationary
input stream; the parameter estimates derived from them converge
to the offline (batch) MLE; and the per-block stability shrinks
below any chosen tolerance. The first block's `γ` values are
biased (computed at init parameters), but for cohorts with
`c_s ≤ 0.05` the asymptotic bias is `O(c̄² / N)` and falls well
below the default `1e-3` stability tolerance after a few
thousand sites. The same bound bounds the difference between
this online estimator and the batch MLE that the previous
version of the spec described.

### Step 5 — apply small-batch and explicit-override floors

Once the Step 2 stopping criterion fires (convergence in
convergence mode; site target in fixed-N mode), the current
parameter estimates are finalised with the following floors and
overrides applied:

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

### Convergence is a property of the running estimate, not of an inner MLE

The Step 4 online-EM update has no inner MLE — there is no
inner loop to converge. What "convergence" means in this spec
is the **outer block-stability** check from Step 2: the per-
block snapshots of `c_s` derived from the running sufficient
statistics stop moving more than the tolerance between
consecutive blocks.

A stable estimate after K = 3 blocks of M = 1000 sites means:
across the last 3 000 informative sites added, the largest
per-sample `c_s` movement was below the tolerance. That is what
the user is buying when they take the `c_s` numbers at face
value.

The block-wise online EM is guaranteed to converge to the
offline (batch) MLE in the limit (see Step 4); what can fail is
reaching the chosen tolerance before exhausting the cohort's
informative sites. Noisy cohorts, very low contamination
signal, or a very tight tolerance can leave the per-block
deltas oscillating.

### `ContaminationEstimateDidNotConverge` — convergence mode only

Raised when the `.psp` files are exhausted before convergence
fires. Carries:

- The final per-sample deltas across the last two blocks.
- The total number of informative sites processed.
- A recommended next-run configuration: extrapolating the
  observed delta from the last two blocks under the expected
  MLE scaling `delta ∝ N^(−1/2)`, the side-pass suggests
  `N_recommended ≈ N_current · (delta_observed / tolerance)²`,
  capped at 10× the current input size to avoid runaway
  recommendations. The user can act on this by relaxing the
  tolerance, by switching to fixed-N mode with
  `--contamination-num-sites N_recommended`, or by relaxing
  `MIN_MAJOR_FRACTION` (admits more sites per genome).

### `ContaminationInsufficientSites` — fixed-N mode only

Raised when `.psp` files are exhausted before
`--contamination-num-sites N` informative sites are found.
Carries:

- `requested: N` — the user-supplied target.
- `found: M` — actual count of informative sites observed.
- A recommended next-run configuration: typically "lower
  `--contamination-num-sites` to ~M", or switch to convergence
  mode with `--contamination-stability-tolerance`, or relax
  filtering thresholds (e.g. lower `MIN_MAJOR_FRACTION`) to
  admit more sites.

### No silent fallback

Neither error type silently falls back to `c_s = 0`. Per Design
principle 3 (no silent defaults), the user must explicitly opt
out of contamination correction — either by omitting
`--contamination-batches` entirely, or by supplying
`--contamination-estimates` with explicit zero values for the
affected samples. A confident-looking VCF whose contamination
correction is missing is a calibration trap the pipeline must
not produce.

### Convergence diagnostics (convergence mode only)

When `--contamination-diagnostics-out FILE` is set, the side-
pass writes a small TSV with one row per block: `block_idx`,
`sites_processed_so_far`, `max_c_s_delta`, `mean_c_s_so_far`,
`min_c_s_so_far`, `max_c_s_so_far`. Useful for tuning the
tolerance against real cohort data.

Setting `--contamination-diagnostics-out` together with
`--contamination-num-sites` is a CLI parse-time error: the
per-block delta is meaningless in fixed-N mode (no stability
check runs), and silently emitting an empty or
delta-less file would hide that fact.

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

**Side-pass enable / replace**

| Flag | Default | Effect |
|---|---|---|
| `--contamination-batches FILE` | — | **Enables the side-pass.** TSV mapping `sample_id → batch_id`. Without this flag, no side-pass runs; Stage 6 uses `c_s = 0`. |
| `--contamination-estimates FILE` | — | Skip the side-pass entirely; use the user-supplied `c_s` values. Format: TSV `sample_id → c_s`. `q_b` is still estimated by a short side-pass restricted to `q_b`-only updates (Step 4 with `c_s` frozen) unless `--contamination-source-distributions FILE` is also supplied. |
| `--contamination-source-distributions FILE` | — | Skip the side-pass entirely; use user-supplied `q_b` values. Requires `--contamination-estimates` to also be set. |

**Stopping criterion (mutually exclusive — passing both is a CLI parse-time error)**

| Flag | Default | Effect |
|---|---|---|
| `--contamination-stability-tolerance T` | `1e-3` (default mode when neither stopping flag is passed) | **Convergence mode.** Run until `max_s |c_s_new − c_s_prev| < T` for `--contamination-stability-blocks` consecutive snapshots. |
| `--contamination-num-sites N` | — | **Fixed-N mode.** Process exactly `N` informative positions and return the resulting estimate. No convergence check. Reproducible by construction. |

**Convergence-mode-only knobs**

| Flag | Default | Effect |
|---|---|---|
| `--contamination-stability-blocks K` | 3 | Consecutive within-tolerance snapshots required to declare convergence. |
| `--contamination-diagnostics-out FILE` | — | Per-block convergence trace (see §"Convergence diagnostics"). **CLI error** when combined with `--contamination-num-sites`. |

**Both modes**

| Flag | Default | Effect |
|---|---|---|
| `--contamination-block-size N` | 1000 | Sites per heartbeat — controls both the EM parameter refresh and (in convergence mode) the stability snapshot. |
| `--contamination-estimates-out FILE` | — | Write side-pass output to disk for inspection / reuse. TSV with per-sample `c_s` and per-batch `q_b`. |
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

**Memory** is `O(n_samples + n_batches × n_allele_classes)` —
the running sufficient statistics from Step 4 and a handful of
scalar parameters per sample. At a 1000-sample cohort with ~10
batches and ~3 allele classes, that is on the order of **20–30
KB** of state regardless of how many sites the side-pass
processes. There is no per-tuple or per-block buffer to size.

**CPU** is dominated by the multi-way per-position scan over
`.psp` files (the same cost as Stage 3 minus DUST). In
convergence mode the side-pass stops early once stability fires
— typically after a few thousand sites, well before whole-
genome scan. In fixed-N mode it stops at exactly `N` sites. The
online-EM per-read work (one closed-form `γ` plus four scalar
adds) is negligible against the `.psp` read cost.

**Parallelism**. The per-position scan parallelises across
chromosome partitions the same way Stage 1 parallelises across
samples. Each partition contributes incremental updates to the
running sufficient statistics; the merge at partition boundaries
is a scalar-wise sum. The block-stability check (convergence
mode) runs sequentially over the merged statistics after each
partition completes a block. Cheap.

For early-stop to work cleanly with partitioned scans (both the
convergence-mode early-stop and the fixed-N-mode exact-stop),
partition work is dispatched in genomic-order chunks and the
side-pass joins partitions in order; once the stopping
criterion fires, in-flight partitions are cancelled. This gives
reproducible results at the cost of slightly less parallelism
than full-genome processing would offer — fine, because the
side-pass is cheap relative to the main pipeline regardless.

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
   In convergence mode the side-pass exhausts the input without
   ever updating the running statistics and returns
   `ContaminationEstimateDidNotConverge` (or, when the cohort
   has no informative sites at all, a dedicated "no informative
   sites" provenance variant).
2. **Two-sample cohort, batch_A contains both, one sample seeded
   with synthetic 3 % contamination from the other.** Side-pass
   recovers `c_s ≈ 0.03` for the contaminated sample within
   tolerance, in convergence mode at the default `1e-3` tolerance.
3. **Singleton batch.** Sample in batch_solo gets `c_s = 0`
   regardless of input data.
4. **Below-floor batch.** Two samples in batch_small (floor = 5)
   both get `c_s = 0`.
5. **User-supplied `--contamination-estimates` overrides MLE for
   the named samples; unnamed samples still get MLE values.**
6. **Convergence diagnostics output format.** TSV header and row
   shape match the spec.
7. **Fixed-N mode hits exactly N sites and returns.** Synthetic
   cohort with abundant informative sites; pass
   `--contamination-num-sites 2000`; side-pass processes exactly
   2000 informative positions, runs no stability check, returns
   the resulting estimate.
8. **CLI mutual-exclusivity error.** Passing both
   `--contamination-stability-tolerance` and
   `--contamination-num-sites` fails at parse time before any
   work happens.
9. **CLI diagnostics-in-fixed-N error.** Passing
   `--contamination-diagnostics-out` together with
   `--contamination-num-sites` fails at parse time.

### Property tests (proptest)

1. **`c_s ∈ [0, 1]` and `q_b` is a simplex** for every cohort
   shape proptest generates, in both modes.
2. **Permuting samples within a batch does not change `q_b` for
   that batch.**
3. **Adding more sites to a converged subsample does not move
   `c_s` beyond the tolerance** (i.e. the convergence criterion
   is consistent with itself) — convergence mode.
4. **Online-EM block-size independence (within tolerance).**
   Running the side-pass on the same input at two different
   `--contamination-block-size` values produces `c_s` estimates
   within `tolerance` of each other once both converge. Probes
   the block-wise refresh design — the heartbeat rate should
   change only how often stability is checked, not the estimator's
   limit.

### Integration tests

1. **End-to-end with `--contamination-batches` enabled,
   convergence mode.** Synthetic cohort with known per-sample
   contamination fractions; side-pass output matches ground
   truth within tolerance; Stage 6 consumes the output and
   emits a VCF whose genotype calls match the no-contamination
   ground-truth more closely than a Stage 6 run with `c_s = 0`
   would.
2. **End-to-end, fixed-N mode.** Same cohort as test 1, with
   `--contamination-num-sites 5000` instead of the default
   tolerance. `c_s` estimates from the two modes agree within
   tolerance.
3. **Non-convergence on a deliberately underspecified cohort
   (convergence mode).** A 4-sample cohort with 1 batch,
   low-coverage everywhere. The side-pass returns
   `ContaminationEstimateDidNotConverge` with a sensible
   recommended next-run config.
4. **Insufficient sites in fixed-N mode.** Same kind of
   underspecified cohort, but invoked with
   `--contamination-num-sites 50000`. The side-pass returns
   `ContaminationInsufficientSites { requested, found,
   recommendation }`.
5. **Round-trip via `--contamination-estimates-out` and
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
