# The per-locus hidden-paralog statistic (H1-vs-H2 likelihood ratio) and the prior that turns it into a probability

**Status:** draft, 2026-06-30. Architecture working doc for the
`tomato2-paralog-filter` branch, continuing
[`hidden_paralog_psp_integration.md`](hidden_paralog_psp_integration.md).
That doc settled how the two per-sample summaries (coverage-by-GC
histogram + observed-het counts) are produced in Stage-1 and carried in
the `.psp` metadata section (shipped, A1–D3). **This doc settles the
consumer**: the per-locus likelihood-ratio statistic that combines those
summaries with the per-locus allele/coverage signals, plus the
data-estimated **prior probability** of a hidden duplication that turns
that raw likelihood ratio into a per-locus probability and a
false-discovery-controlled cut.

Model *intent* lives in the spec [`hidden_paralog_filter.md`](../specs/hidden_paralog_filter.md)
(§5 the likelihood ratio, §6 the empirical-Bayes cut, §8 the pre-pass
architecture); the validated
maths reference is the tomato2 prototype
([`build_paralog_lr.py`](../../../benchmarks/tomato2/src/build_paralog_lr.py),
[`build_paralog_eb.py`](../../../benchmarks/tomato2/src/build_paralog_eb.py),
[`build_gc_normalization.py`](../../../benchmarks/tomato2/src/build_gc_normalization.py)).
Grown one agreed premise at a time. **Types-first; no production wiring
into var-calling in this slice** (spec: emitting the FILTER verdict is a
later slice — §6 of this doc).

---

## Scope of this slice

In: the *statistic and its prior*, built as a small module whose
functions are each **pure** — given the same inputs they always return
the same output, and they touch nothing else (no files, no global state,
no randomness). That is what lets each one be unit-tested and
parity-checked against the prototype in isolation:

1. **`SingleCopyCoverageModel`** — from a sample's stored
   `CoverageByGcHistogram`, fit (a) its single-copy depth level and
   (b) its GC-bias curve, so that at any locus the observed window depth
   becomes a **relative copy number** (`1.0 = one copy`) via
   `observed_depth / (single_copy_scale · gc_bias_curve(GC))`.
2. **The scoring function** — `LocusObservations → ParalogScore` (a pure
   function of each sample's `(relative_copy_number, alt_reads,
   total_reads, inbreeding_coefficient)`). Its output is the
   **`paralog_log_likelihood_ratio`**:
   `log P(data | hidden paralog) − log P(data | real variant)` — a *pure*
   likelihood ratio, which is what lets the prior/FDR step (piece 3) treat
   it as a Bayes factor. Divergence (MQDiff) is deliberately excluded from
   the score (it is shared with introgression — spec §2, §3) and lives in
   the VCF as an INFO field.
3. **The prior probability of a hidden duplication** — turn each locus's
   raw score into an actual probability and a false-discovery rate. The score tells us which
   explanation fits better, but not how *likely* the locus is to be a
   paralog; for that we also have to know how common paralogs are across
   the genome. We estimate that genome-wide rate **from the data itself**
   (an iterative loop that repeatedly re-scores every locus and re-averages
   the rate until it stops moving), fold it into each locus's score to get
   the probability that the locus is a paralog, and rank the loci by that
   probability to assign a q-value — so a cut can be chosen that removes
   paralogs while sacrificing only a small, controlled fraction of real
   variants. (Exact formulas and symbols: Premise 4.)
4. **A check that the Rust code behaves on real data.** A small runnable
   program (a Rust `examples/` binary, like the existing
   `dump_sample_summary.rs` we already use to cross-check the coverage
   summaries) that runs pieces 1–3 on the real tomato2 test data. The
   prototype is a *reference draft*, not a bit-exact oracle (the production
   model has deliberately moved past it), so validation is against the
   **data** — π in the expected range, the flagged set paralog-like, the
   single-copy peak anchored at `relative_copy_number = 1.0` — with a **loose
   Python↔Rust LR correlation** as a porting sanity-check. (Details: Premise 5.)

Out (later slices, §6): reading every sample's summary up-front in
`var_calling::pipeline::run_var_calling`, computing each locus's
relative copy number from the `.psp` body at scoring time, and emitting a
`FILTER` verdict into the VCF.

---

## Premise 0 — where the consumer lives

A **new top-level module `src/paralog/`**, sibling to `src/sample_summary/`
(which holds the producer side). It holds only the statistics functions
and their types; it depends on `sample_summary` (for
`CoverageByGcHistogram`, `HetCounts`) and on nothing in `var_calling`.
Layout:

```
src/paralog/
    mod.rs            // re-exports, the model-params struct, top-level docs
    coverage_model.rs // SingleCopyCoverageModel::fit + relative_copy_number
    locus_score.rs    // LocusObservations, SampleObservation, ParalogScore + score_locus_for_paralogy
    prior.rs          // estimate the prior probability, then posterior + FDR q-values
examples/
    paralog_score_parity.rs // runnable tomato2 parity check (Scope item 4)
```

**Why a standalone module, not inside `var_calling`:** these functions
are pure and have no dependency on the cohort plumbing
(`CohortPileupRecord`, `MergedRecord`, the producer/caller/writer
stages). Keeping them separate (a) makes them unit-testable and
parity-checkable in isolation, (b) lets the later var-calling wiring
slice depend *inward* on a settled interface, and (c) mirrors how
`sample_summary` is a standalone module the pileup driver calls into.
The var-calling integration (§6) will construct `LocusObservations` from
`PosteriorRecord.scalars` + the fitted per-sample
`SingleCopyCoverageModel`s and call the same function.

---

## Premise 1 — `SingleCopyCoverageModel`: what single-copy depth looks like in one sample *(proposed)*

**What it is.** A per-sample model that answers "at a window of this GC
content, what read depth would *one copy* produce in this sample?" It
holds two fitted quantities:

- `single_copy_scale` — the sample's typical one-copy depth level (it
  averages 6 reads/base, say). This divides out **how deeply the sample
  was sequenced**.
- `gc_bias_curve` — a multiplier indexed by GC content (×0.8 where
  GC-poor, ×1.1 mid-GC). This divides out the **GC sequencing bias**.

Their product, `single_copy_scale · gc_bias_curve(GC)`, is the
**expected single-copy depth** at that GC. Dividing an observed window
depth by it gives the **relative copy number** — how many copies' worth
of reads are actually there, with `1.0 = one copy`:

```
relative_copy_number = observed_depth / (single_copy_scale · gc_bias_curve(GC))
```

**Why fit it from the histogram, not raw windows.** The prototype fits
from the *raw* per-window `(GC, depth)` table
([`build_gc_normalization.py`](../../../benchmarks/tomato2/src/build_gc_normalization.py));
we fit from the stored **binned 2-D histogram** `[gc_bin][depth_bin]`
(architecture Premise 2 kept the histogram raw precisely so the consumer
can re-fit). The fit, step by step:

1. **`single_copy_scale`.** The **mode** (dominant peak) of the per-tile
   depth distribution in the histogram — the densest depth bin, or a
   smoothed kernel-density / half-sample-mode peak over the bin counts. In
   the post-filter, covered-base distribution the single-copy windows form
   that peak (mapping-quality filtering strips repeat coverage, so repeats
   lose coverage or drop out rather than showing excess; collapsed
   duplications are a small high tail), so the peak marks the one-copy
   level. We use the mode, not the prototype's median: the median assumes
   single-copy windows are an outright *majority*, which can fail on a
   duplication-rich or polyploid genome, whereas the mode only needs
   single-copy to be the *most common* state. We take the **global**
   (dominant) mode, *not* the lowest/"first" peak: the low end carries small
   bumps from het-deletions and low-mappability edges (they read ~0.4–0.6),
   and a first-peak rule could latch onto one of those and mis-scale upward.
   The dominant peak is also the right reference for a true polyploid (it
   lands on the baseline ploidy, so `1.0` = one baseline dose). The one case
   the global mode could get wrong — a genome so duplicated that a >1×
   peak becomes the plurality — is not reachable with realistic data
   (mapping-quality filtering strips the repeat coverage that would build a
   rival peak) and is caught, not left silent, by a **sanity guard**: a
   healthy sample has `mode/median ≈ 0.83–0.98` (measured on tomato2), so we
   reject the fit (`CoverageModelError`) when `mode/median` falls outside
   `[0.5, 1.5]` — a mode near twice the median means it landed on a 2× peak.
   The mode's precision is bounded by `depth_bin_width`, so that width is kept
   a small fraction of the single-copy level (≈ 0.05× or finer) to pin the
   scale to a couple of percent. *(Divergence from the prototype — it
   rescales every `relative_copy_number` by a per-sample constant, so the
   per-tile correlation is preserved.)* The
   prototype also filters `breadth > 0.8`; we **do not need to**. The stored
   depth is a *covered-bases mean* (`depth_sum / covered_positions`) —
   identical to the prototype's own `mean_depth` (`sum_depth / n_cov`;
   verified equal on tomato2), so there is no depth-definition difference to
   reconcile. And the breadth filter is unnecessary: low-breadth windows read
   *low*, not single-copy (on tomato2, covered-bases mean ≈ 1.7 vs ≈ 5.7
   normal — poor mappability depresses breadth and per-covered-base depth
   together), so they self-segregate to the low end and never touch the
   single-copy peak. Excluding them (5.7% of windows) leaves the single-copy
   mode unchanged and tightens its robust SD by only ~3% — not worth the
   breadth-plumbing (the accumulator sees only covered positions, so breadth
   would need the tile's non-`N` length supplied from the reference).
2. **Per-tile relative depth.** `rel = tile_depth / single_copy_scale`
   (1.0 = single copy).
3. **`gc_bias_curve`.** For each GC bin, the multiplier is the weighted
   median of `rel` over that bin's tiles, using only tiles in the
   single-copy band `0.4 < rel < 1.6` (so paralog/deleted windows don't
   distort the bias estimate — mirrors the prototype). GC bins with too
   few tiles (default `< 50`) are filled by linear interpolation from
   their neighbours, then a 5-point median smooth flattens noise.
   `gc_bias_curve(GC)` interpolates linearly between GC-bin centres.
4. **`single_copy_depth_sd`** (σ₀ in the maths). The robust spread
   (`1.4826 · MAD`) of `relative_copy_number` across single-copy tiles —
   how much a true one-copy window's coverage scatters. It feeds the
   Normal coverage term in the LR (Premise 2). **Fit per sample from the
   (mode-anchored) data** — ≈ 0.28 on tomato2, higher than the prototype's
   hardcoded `0.26` because the mode anchoring stretches the distribution
   slightly; a constant override stays available but is not the default.

The model carries only `single_copy_scale` + `gc_bias_curve` +
`single_copy_depth_sd`; the actual `(GC, depth)` at a locus is measured
*fresh* from the `.psp` body at scoring time. The histogram is just the
training data.

```rust
/// What one copy's read depth looks like in a single sample, as a
/// function of GC content. `relative_copy_number == 1.0` means one copy.
pub struct SingleCopyCoverageModel {
    window_bp: u32,             // inherited from the histogram (= analysis window)
    single_copy_scale: f64,     // the sample's typical one-copy depth level
    gc_bins: u32,
    gc_bias_curve: Vec<f64>,    // per-GC-bin coverage multiplier (smoothed, gap-filled)
    single_copy_depth_sd: f64,  // robust SD of a one-copy window's relative depth (σ₀)
}
impl SingleCopyCoverageModel {
    pub fn fit(hist: &CoverageByGcHistogram, cfg: &CoverageFitConfig)
        -> Result<Self, CoverageModelError>;
    /// Expected depth of a single-copy window at this GC content
    /// = `single_copy_scale · gc_bias_curve(gc_fraction)`.
    pub fn expected_single_copy_depth(&self, gc_fraction: f64) -> f64;
    /// Observed window depth expressed in copies (1.0 = one copy):
    /// `observed_depth / expected_single_copy_depth(gc_fraction)`.
    /// `gc_fraction ∈ [0,1]`, `observed_depth ≥ 0`.
    pub fn relative_copy_number(&self, gc_fraction: f64, observed_depth: f64) -> f64;
    pub fn single_copy_depth_sd(&self) -> f64;
    pub fn window_bp(&self) -> u32;
}
```

`CoverageFitConfig`: `min_bin_count` (50), `single_copy_lo`/`hi` (0.4,
1.6), `smooth_window` (5), `single_copy_depth_sd_override: Option<f64>`,
`mode_median_ratio_bounds` (default `[0.5, 1.5]`). Failure modes
(`CoverageModelError`): histogram with no single-copy tiles (degenerate
sample), all-zero `single_copy_scale`, or `mode/median` outside
`mode_median_ratio_bounds` (the mode landed on the wrong copy-number peak —
raised rather than silently mis-scaling).

---

## Premise 2 — the scoring function: per-locus input/output types *(proposed)*

The scoring computation is a pure function, run on **every** locus to
decide whether it looks like a hidden paralog (the name carries no claim
that it already is one). Per sample at a locus it consumes the
`relative_copy_number` (Premise 1), the allele read counts (AD):
`alt_reads` out of `total_reads`, and the sample's
`inbreeding_coefficient` (the F of Premise 3). It consumes **no divergence
signal**: MQDiff is shared with introgression and cannot safely separate
the two (spec §2, §3), so it is kept out of the score and emitted as a VCF
INFO field instead. That keeps the score a pure likelihood ratio (below).

```rust
/// One sample's evidence at a locus — computed for EVERY locus we score,
/// not only known paralogs.
pub struct SampleObservation {
    pub relative_copy_number: f64,   // window depth in copies (1.0 = one copy)
    pub alt_reads: u32,              // reads supporting the ALT allele(s) (AD alt)
    pub total_reads: u32,            // total reads at the site            (AD sum)
    pub inbreeding_coefficient: f64, // the sample's F (selfer↔outbred), [0, 0.99]
}

/// All samples' evidence at one locus.
pub struct LocusObservations<'a> {
    pub samples: &'a [Option<SampleObservation>], // None = no usable data for that sample
}

/// The paralog verdict for one locus.
pub struct ParalogScore {
    /// log P(data | hidden paralog) − log P(data | real variant). A *pure*
    /// likelihood ratio: > 0 favours hidden paralog; < 0 favours a real
    /// single-copy variant.
    pub paralog_log_likelihood_ratio: f64,
    pub samples_used: usize,             // samples with usable data
    pub confident_homalt_carriers: usize,// confident hom-alt carriers (hom-alt veto signal)
    // diagnostics, exposed for the parity check:
    pub log_likelihood_real_variant: f64,    // log P(data | real variant)  (H1)
    pub log_likelihood_hidden_paralog: f64,  // log P(data | hidden paralog) (H2)
}

/// The fixed model grids/constants (defaults reproduce the prototype).
pub struct ParalogModelParams {
    pub pseudocount_vaf: f64,        // 0.01  VAF floor for hom-ref / non-carrier
    pub max_relative_copy_number: f64, // 4.0  winsorise relative_copy_number for the Normal tail
    pub carrier_copy_numbers: Vec<u32>, // [3,4,6,8]  hidden-paralog total copies (T ≤ 8)
    pub allele_freq_prior: SfsPriorSpec, // folded SFS 1/(p(1-p)) on [1/2N, 1-1/2N], ~200 pts
    pub carrier_freq_grid: GridSpec,    // 40 points, 0.004..0.6 (carrier freq, flat). The 0.6
                                        // ceiling can't reach fixation, but that's fine: a fully
                                        // fixed duplication is caught by absolute coverage, not the
                                        // carrier-freq term (raising the ceiling to 0.9/0.99 barely
                                        // moves π — measured, spec §7).
    pub homalt_vaf_threshold: f64,   // 0.9  veto: VAF threshold for "confident hom-alt"
    pub homalt_min_depth: u32,       // 5    veto: depth threshold
}

/// Score one locus. The per-sample `single_copy_depth_sd` (σ₀) comes from
/// each sample's `SingleCopyCoverageModel` (or the params override); it is
/// passed alongside the observations so the function depends only on its
/// arguments (no hidden state).
pub fn score_locus_for_paralogy(
    observations: &LocusObservations,
    single_copy_depth_sd: &[f64], // per-sample σ₀ (one-copy relative-depth SD)
    params: &ParalogModelParams,
) -> ParalogScore;
```

**The maths the function implements** (faithful to `build_paralog_lr.py`).
The derivation uses the conventional single-letter symbols for the
quantities being *integrated out* — `p` = allele frequency, `q` = carrier
frequency, `T`/`m` = total/mutant copies — defined inline below; the
durable field names (`alt_reads`, `total_reads`, `inbreeding_coefficient`,
…) are the interface above.

- **`log_likelihood_real_variant` (H1 — real, single-copy variant).**
  Marginalise over the allele frequency `p` under a **folded site-frequency-
  spectrum prior** (`∝ 1/(p(1−p))` on `p ∈ [1/2N, 1−1/2N]`, `N` = samples),
  not a flat grid — real variants are mostly rare, and this removes the
  arbitrary low-frequency cutoff (measured to make the paralog rate robust
  to it; spec §5.2, §7). Per sample, the genotype prior is the Wright
  inbreeding-adjusted HWE with the sample's `inbreeding_coefficient` `F`:
  `P(homref)=q²+Fpq, P(het)=2pq(1−F), P(homalt)=p²+Fpq` (here `q=1−p`);
  the allele likelihood is
  `alt_reads·ln(vaf) + (total_reads−alt_reads)·ln(1−vaf)` over
  `vaf ∈ {pseudocount, ½, 1−pseudocount}`; coverage is `~ Normal(1, σ₀)`,
  **independent of genotype**. Marginalise genotypes per sample by
  log-sum-exp, sum over samples, marginalise `p` by `LSE − ln(grid)`, add
  the coverage term.
- **`log_likelihood_hidden_paralog` (H2 — hidden paralog).** Per carrier
  total copies `T ∈ carrier_copy_numbers` (`{3,4,6,8}`), coverage is
  `~ Normal(T/2, σ₀·√(T/2))` and `vaf = m/T`, keeping only the single-PSV
  (`m=1`) and balanced (`m≈T/2`) configurations with `1 ≤ m ≤ T−2`. A
  non-carrier is `Normal(1, σ₀)` with `vaf = pseudocount`. The
  carrier-vs-non-carrier prior is the per-sample dosage HWE in the carrier
  frequency `q`; marginalise carrier status per sample by `logaddexp`, sum
  over samples, then marginalise configurations × `q` by `LSE − ln(size)`.
  *Note the top config's mean `T/2 = 4` equals `max_relative_copy_number` (4.0),
  so any ≥ 4× coverage — a T ≥ 8 collapse or a high-copy artefact — winsorises
  to 4 and reads as "looks like the top config." That is fine for detection (all
  such windows are paralog-like) but means the model does not distinguish
  4× from 8×+, and a T=8 carrier's upward scatter is clipped. If the config set
  or the cap is ever retuned, keep the winsor cap above the top carrier mean +
  a few σ (≈ `T_max/2 + 3σ₀·√(T_max/2)`) so the top config is not pinned at the
  clip boundary.*
- **The result:**
  `paralog_log_likelihood_ratio = log_likelihood_hidden_paralog −
  log_likelihood_real_variant`. Nothing is added to it: keeping it a pure
  likelihood ratio is what makes the `sigmoid(ratio + prior log-odds)`
  posterior of Premise 4 a valid Bayes-factor update. Divergence (MQDiff)
  is deliberately excluded (spec §2, §3) and emitted as a VCF INFO field.
- **Hom-alt veto** falls out of H2's capped VAF; `confident_homalt_carriers`
  (samples with `vaf > homalt_vaf_threshold` and
  `total_reads ≥ homalt_min_depth`) is reported for diagnostics.

All sums use a stable log-sum-exp helper. The function never allocates
per sample beyond small fixed-width scratch.

---

## Premise 3 — the inbreeding scalar `F` from a *correctly computed* obs het *(proposed — resolved)*

The Wright genotype prior needs a per-sample `F`. It comes from the
sample's observed heterozygosity and one cohort-wide expected
heterozygosity: `F_s = clip(1 − Hobs_s / Hexp, 0, 0.99)`.

**`Hobs_s` is the het *rate*, `n_het / callable_positions`** — het count
over *all* callable positions, a single-sample observable from the Stage-1
pileup. An earlier draft used `n_het/(n_het+n_hom_alt)` (het over *variant*
sites) and mapped `F = 1 − Hobs`; that is **wrong** and was measured so:
that ratio is dominated by the hom-alt count, i.e. by the sample's
divergence from the reference, so it *inverts* — on tomato2 it correlates
with the true cohort `F` at Spearman **−0.28**, and the most-outbred
samples come out at `F ≈ 0.83`. The het *rate* fixes this (Spearman **0.86**
with the cohort `F`; outbred samples ≈ 359 het/Mb vs selfers ≈ 47 het/Mb).
The only new stored ingredient is the **callable-position total**, which the
coverage accumulator already counts per tile — it just keeps the grand total
(a single `u64`).

**`Hexp` is one cohort scalar (`mean 2pq`), accumulated inside var-calling's
existing per-locus genotype pass** — *not* a separate pre-pass. Var-calling
already visits every locus and computes the cohort allele frequencies, so
`Σ 2p(1−p)` falls out for free; the completed `Hexp` is ready before the
paralog scoring pass (which is already global — it needs every locus's LR
for the prior/FDR). This handles every population type: an all-selfing
cohort still gets a correct `Hexp` from real allele frequencies (a
percentile-of-obs-het shortcut would underestimate it there), an F2 gets
`Hobs ≈ Hexp → F ≈ 0` (correct — F2 sibs are outbred).

**Why this is enough.** `F` is a *secondary* knob under coverage: on
tomato2, swapping the (broken) proxy for the exact cohort `F` moved the
flagged set by < 10% (Jaccard 0.92) and left π unchanged. So even a rough
`F` barely changes the flagging — but we compute it correctly anyway, so
the `inbreeding_coefficient` field is trustworthy across population types,
not just where reference divergence and inbreeding happen to align.

**Validation (Premise 5):** the `Hobs`/`Hexp` → `F` derivation is checked on
its own terms — Spearman agreement with the prototype's cohort `F` on tomato2
(0.86) — not by reproducing it exactly (the prototype uses per-locus
`2p_iq_i`; production uses one genome-wide `Hexp`). The end-to-end check is
that the flagged set on tomato2 is stable when `F` is formed this way.

*(Producer-side change owed in the implementation plan: the pileup stores
the callable-position total alongside `HetCounts`; `Hobs` is then formed
downstream, not `n_het/(n_het+n_hom_alt)`.)*

---

## Premise 4 — the prior probability that a locus is a hidden duplication (then the posterior and FDR) *(proposed)*

**What it does, in words.** The per-locus score
(`paralog_log_likelihood_ratio`) says which explanation fits the data
better, but on its own it is not a probability — a locus with a slightly
positive score is not necessarily *probably* a paralog, because paralogs
are rare to begin with. To get a real probability we combine the score
with the **prior probability that a locus is a hidden duplication**
(`prior_probability`) — i.e. how common such duplications are across the
genome: the rarer they are, the more positive a score has to be before we
believe it.

We don't know that rate up front, so we estimate it **from the data
itself**. This is the iterative loop (the EM algorithm): start from a
guess, use it to assign every locus a paralog probability, average those
probabilities to get a better rate estimate, and repeat until the estimate
stops changing. With the rate in hand, each locus's score becomes a
posterior probability through the standard log-odds ↔ probability
conversion (`logit`/`sigmoid` in the glossary).

Finally we rank loci by that probability and assign each a **q-value** —
the fraction of real variants we'd wrongly remove if we cut at that locus.
That lets the operating point be set by false-discovery rate (default ≈
1%), which is exactly the introgression-preservation knob.

The whole thing is a pure function of the genome-wide vector of
`paralog_log_likelihood_ratio`s (faithful to `build_paralog_eb.py`):

```rust
/// The genome-wide prior probability that a locus is a hidden duplication
/// (present in the samples but collapsed to one copy in the reference),
/// estimated from the data. Applying it to a locus's score yields that
/// locus's posterior probability.
pub struct ParalogPrior {
    pub prior_probability: f64, // prior P(locus is a hidden duplication); the π of empirical Bayes
}

impl ParalogPrior {
    /// Estimate the prior from all loci's scores (the EM loop):
    /// posterior_i = sigmoid(ratio_i + logit(prior_probability));
    /// prior_probability ← mean(posterior); iterate to |Δ| < tol.
    pub fn estimate(paralog_log_likelihood_ratios: &[f64], cfg: &EmConfig) -> Self;
                                                 // start 0.03, tol=1e-9, max_iter=500
    /// Posterior probability the locus is a hidden duplication:
    /// `sigmoid(ratio + logit(prior_probability))`.
    pub fn paralog_posterior(&self, paralog_log_likelihood_ratio: f64) -> f64;
    /// The log-likelihood-ratio at which the posterior crosses 0.5:
    /// `log((1 − prior_probability) / prior_probability)` — positive, since
    /// hidden duplications are rare (the cut is NOT at ratio = 0).
    pub fn half_posterior_ratio(&self) -> f64;
}

/// Tail-FDR q-values: sort by posterior descending, take the running mean
/// of the local FDR (`1 − posterior`), re-index to input order.
pub fn q_values(paralog_posteriors: &[f64]) -> Vec<f64>;
```

The decision cut is **FDR-based and positive**, never `ratio = 0`
(spec §6): default high-confidence ≈ FDR 1%. The cut selection
(`q ≤ q_target`) is a thin helper, not state.

---

## Premise 5 — validating the production statistic *(proposed — revised)*

The Python prototype is a **reference draft**: it established the model and
let us validate the *mechanisms*. This review has deliberately moved the
production model past it — divergence dropped from the score (spec §3), a
folded site-frequency-spectrum prior on the allele frequency (Premise 2),
observed heterozygosity computed as the het rate with an inline `Hexp`
(Premise 3), the coverage likelihood kept Normal (the heavy coverage tail is
duplication signal, not noise — spec §7). So **bit-exact reproduction of the
prototype is neither the goal nor well-defined** (reproduce *which* draft?).
Validation is against **data**, not against the draft:

- **Primary — the Rust production code on the tomato2 data**, the same
  evidence base used throughout this review: π lands in the expected range,
  the flagged set carries the paralog profile (het excess + coverage excess),
  and the coverage model is *absolutely* calibrated — the single-copy peak
  sits at `relative_copy_number = 1.0`. That last check is the one a
  correlation cannot see (a constant rescale leaves correlation at 1.0 while
  moving the anchor the LR depends on), so it is checked directly: anchor the
  scale with the **mode** (which lands the single-copy peak at 1.0; the
  prototype's median mis-anchors it — measured 2.5% low typically, up to 17.5%
  on duplication-rich samples), and **fit σ₀ from the data** (≈ 0.28 on
  tomato2) rather than carrying over the prototype's constant.
- **Absolute rate — the ground truth** (read-level simulation / known
  duplications, §9), the one thing neither the prototype nor the tomato
  callset can supply.
- **Python↔Rust — a loose correlation sanity-check** (`examples/paralog_score_parity.rs`,
  mirroring the existing `dump_sample_summary.rs`): run both on tomato2 and
  confirm the per-locus LRs correlate highly, to catch gross porting bugs
  (a transposed argument, a sign error) — **not** an exact-parity gate.

The load-bearing point for this slice is unchanged: the statistics functions'
input types are exactly what the var-calling structs can produce, so the
later wiring (Premise 6) is a construction step, not a redesign.

---

## Premise 6 — the later var-calling wiring (deferred, sketched) *(not in this slice)*

For continuity. The integration slice will, in
`var_calling::pipeline::run_var_calling` (around the point all
`PspReader`s are open, before they are moved into the producer):

1. read each sample's metadata section → `SampleSummary`;
   `SingleCopyCoverageModel::fit` per sample; form each sample's
   `Hobs = n_het/callable`. Accumulate the cohort `Hexp = mean 2pq` in the
   per-locus genotype pass, then `F_s = 1 − Hobs_s/Hexp` (Premise 3).
   → a `ParalogPrePass` value.
2. at per-locus scoring (after `PosteriorRecord` is built, where
   `scalars` give per-sample `(alt_reads, total_reads)`), measure the locus
   tile's `(gc, mean_depth)` from the body, turn it into a
   `relative_copy_number` via the per-sample `SingleCopyCoverageModel`,
   build `LocusObservations`, call `score_locus_for_paralogy`. (MQDiff is
   still written to the VCF INFO as before — it is just not consumed by the
   score.)
3. **Score, then calibrate, then write — bounded RAM via a temporary spill.**
   The prior and FDR are *global* (π averages over every locus; the q-value
   ranks every locus), so no locus's verdict can be emitted until all loci
   are scored — there is no single-pass streaming form. To avoid holding a
   genome-wide per-locus vector in RAM, pass 1 scores each locus and **spills
   `(position, paralog_log_likelihood_ratio, and the fields the final VCF
   needs)` to a temporary binary file** (parquet or a `.psp`-style framing —
   columnar helps, since calibration reads only the LR column), while
   accumulating a **fixed-size histogram of the LRs** in RAM (a few KB). From
   that histogram alone — no re-read — `ParalogPrior::estimate` gives π (the
   EM is a weighted mean over bins) and the monotone `q(LR)` curve (the
   tail-FDR as a function of the LR). Pass 2 reads the spill once and writes
   the final VCF, stamping each locus's posterior / `FILTER` from `q(LR_i)` —
   the VCF is written once, correctly, with no in-place FILTER patching. RAM
   stays flat in both variant count and sample count; the only cost is a
   temporary, compressed, ~VCF-sized file on disk (the resource we *want* to
   spend for a memory-scaling tool).

   **This spill is an ephemeral scratch file — *not* a `.psp`-style artifact.**
   The `.psp` is a durable, first-class output built once per sample and
   reused across many callings. This file exists only to bound memory within
   a single var-calling run: it lives in scratch space and is deleted when the
   run finishes (and on failure). It is an internal implementation detail, not
   a format users see or that persists.

The point load-bearing now: **the statistics functions' input types are
exactly what the var-calling structs can produce**, so this wiring is a
construction step, not a redesign.

---

## Open items carried forward

- Wire the Python↔Rust LR-correlation sanity-check into
  `examples/paralog_score_parity.rs` (Premise 5) — a loose porting check, not
  a golden-fixture exact-parity gate.
- Producer-side: store the callable-position total alongside `HetCounts`
  (Premise 3); optional later `n_ambiguous`-weighted shrinkage of `Hobs`
  for low-coverage samples.
- Format for the temporary calibration spill (parquet vs a `.psp`-style
  framing) and its histogram bin resolution (Premise 6) — the two-pass +
  streaming-histogram design itself is settled.
