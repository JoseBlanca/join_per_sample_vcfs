# Filtering SNP/indel calls caused by hidden (reference-collapsed) paralogs

**Status:** draft, 2026-06-30 (rewritten for clarity; the model is unchanged
from the 2026-06-28 draft). Model-intent spec — *what* we compute to remove
false calls produced by duplications that are present in the sequenced genomes but
**collapsed into a single copy in the reference assembly**, and *why*. Grown one
agreed premise at a time during an exploratory study; the empirical work and all
figures live in [`benchmarks/tomato2/`](../../../benchmarks/tomato2/) (59-sample
*S. lycopersicum* WGS sliced to 160×200 kb random regions, ~6× per-sample depth).
The module/struct/signature shape (the *how*) and the production integration are
deferred to a later architecture doc; §8 sketches the intended architecture.

---

## 0. Vocabulary

This document and its architecture companion
([`hidden_paralog_locus_statistic.md`](../architecture/hidden_paralog_locus_statistic.md))
mix genetics, statistics, and a little software vocabulary. The genetics
you already know; the statistics terms are the ones worth pinning down,
so they get the fuller explanations. Skip whichever group you don't need.

### Genetics / model terms

- **Hidden (reference-collapsed) paralog.** A segment that is duplicated
  in the sequenced genome but assembled as a *single* copy in the
  reference. Reads from all copies pile onto the one reference locus.
  Where the copies differ in sequence, this manufactures a fake
  heterozygote and extra coverage — the artefact this whole filter exists
  to remove.
- **Introgression.** A genuinely single-copy but divergent allele,
  carried in from another population/species. It looks like a paralog at
  the variant level (a divergent ALT allele) but is *real* and must be
  kept. Telling the two apart is the entire problem (§2).
- **PSV — paralogous sequence variant.** A position where two collapsed
  copies differ. In the model it shows up as a single mutated copy in an
  array (so a *low* fraction of the reads carry the ALT — see VAF).
- **HWE — Hardy–Weinberg equilibrium.** The textbook relation between
  allele frequency `p` and genotype frequencies (`p²`, `2pq`, `q²`). We
  use the **inbreeding-adjusted** form (next entry).
- **Inbreeding coefficient `F`.** A number in `[0, 1]` measuring how
  self-fertilising / inbred an individual is. `F = 0` is fully outbred;
  `F → 1` is a near-pure selfer. It bends the HWE genotype proportions
  toward homozygotes: `P(het) = 2pq(1−F)`. A selfer almost never carries
  a real variant as a heterozygote, which is why a het call in a selfer
  is suspicious — a key signal the model uses.
- **Observed heterozygosity `Hobs`.** The fraction of an individual's
  **callable sites** that are heterozygous — the het *rate*, `n_het /
  callable_positions`. It is a *single-sample* quantity (needs no cohort)
  and places the sample on the selfer↔outbred axis — low `Hobs` ≈ high `F`.
  Note the denominator is *all callable positions*, **not** the variant
  sites: dividing het by variant sites (`n_het/(n_het+n_hom_alt)`) instead
  measures reference divergence, not heterozygosity, and inverts on
  reference-divergent samples (measured — §7). The `F` the model needs is
  derived from `Hobs` downstream (next entry).
- **Expected heterozygosity `Hexp`, and `F = 1 − Hobs/Hexp`.** `Hexp` is
  the heterozygosity a *randomly-mating* cohort would show (`mean 2pq` over
  sites). Dividing `Hobs` by it turns the raw het rate into the inbreeding
  coefficient `F`. `Hexp` is one cohort-wide number, accumulated **inside
  var-calling's existing per-locus genotype pass** (not a separate pre-pass)
  from the allele frequencies it already computes. So `Hobs` stays a stored
  single-sample summary and `F` is formed from it with one shared `Hexp`.
- **VAF — variant allele fraction.** Of the reads at a site, the fraction
  supporting the ALT allele. A real heterozygote sits near ½; a real
  homozygous-ALT near 1; a single mutated copy among `T` collapsed copies
  sits at `m/T` (low).
- **AD — allele depth.** The per-sample read counts behind each allele:
  `alt_reads` (supporting the ALT) out of `total_reads` (all alleles), with
  `alt_reads ≤ total_reads`. (Here "reads" are post-dedup,
  mate-overlap-resolved *fragments*.)
- **Copy number `T`, mutant copies `m`.** Under the paralog hypothesis a
  carrier has `T` total collapsed copies, `m` of them mutated → coverage
  `≈ T/2` (relative to single-copy) and VAF `= m/T`.
- **Relative copy number (GC-corrected relative coverage).** Windowed read
  depth after dividing out (a) the sample's overall sequencing depth and
  (b) its GC bias, so that **`1.0` = one copy**, `2.0` = two copies, and so
  on. The introgression-safe core of the filter (§4). (The tomato2
  prototype calls this `gc_rel`; the Rust code calls it
  `relative_copy_number`.)
- **Single-copy scale.** The depth level that corresponds to one copy in
  a given sample — the denominator that turns raw depth into the relative
  copy number.
- **MQDiff (mapping-quality difference).** Mean mapping quality of ALT
  reads minus that of REF reads. Negative means the ALT reads align worse
  (more divergent). It is **not part of the filter score** — divergence is
  a signal shared with introgression, so it cannot safely separate the two
  (§2, §3). We emit it as a per-locus VCF INFO field for inspection, and it
  is the natural hook for a future, clearly separate rescue of the
  diverging-paralog blind spot (§7).
- **H1 / H2.** The two competing explanations for a locus. **H1** = "this
  is a real, single-copy variant" (the keep hypothesis). **H2** = "this
  is a hidden paralog" (the remove hypothesis).

### Statistics terms

- **Likelihood (and log-likelihood).** The probability of the *observed
  data* under a given model — "if H2 were true, how probable are these
  read counts and this coverage?" We work with its logarithm
  (log-likelihood) because the numbers are tiny and logs turn products
  into sums.
- **Likelihood ratio (LR).** Here, `log P(data | H2) − log P(data | H1)`.
  Positive favours paralog, negative favours real variant. It is the core
  per-locus score. (The decision cut is *not* at `LR = 0` — see posterior
  / FDR.)
- **Prior / posterior.** The **prior** is what you believe before seeing
  this locus's data (e.g. "paralogs are rare"); the **posterior** is the
  updated belief after combining the prior with the likelihood. Bayes'
  rule is just "posterior ∝ prior × likelihood".
- **Flat (uniform) prior.** Declining to prefer any value of an unknown
  parameter — every value in range is weighted equally when we average
  over it. We use it for the paralog's carrier frequency `q`; for the
  real-variant allele frequency `p` we use the SFS prior instead (next).
- **Site-frequency spectrum (SFS) prior.** The population-genetic
  expectation that, at real single-copy variants, *most alleles are rare*.
  Instead of weighting every allele frequency `p` equally, we weight it by
  the neutral **folded SFS**, `∝ 1/(p(1−p))`, when averaging H1 over `p`
  (§5.2). "Folded" = symmetric in `p` and `1−p`, because the ALT allele is
  defined against the reference, not polarised to ancestral. It is bounded
  at `1/2N` and `1−1/2N` (`N` = samples) — the frequency of a single copy,
  the finest the panel can resolve — which removes the need for an
  arbitrary low-frequency cutoff.
- **Marginal likelihood / "marginalising over a parameter".** When a
  hypothesis has an unknown knob — for H1, the allele frequency `p`; for
  H2, the carrier frequency — we do **not** pick the single value that
  fits best. We average the likelihood over *all* plausible values of the
  knob (weighted by a prior). "Marginalise" is just the statisticians'
  word for "average it out". This matters because the more flexible
  hypothesis would always win if each were allowed to pick its
  best-fitting knobs; averaging instead charges a hypothesis for its
  flexibility — the built-in **Occam penalty** that stops us calling
  everything a paralog (§5.5).
- **Marginalising vs maximising.** Two ways to handle that averaging.
  *Maximising* takes the single best-fit value of the unknown (and over-
  rewards flexible models). *Marginalising* averages properly over all its
  values — here by a **grid average**: evaluate the likelihood on a grid of
  values and average them (weighted by the prior), done in logs with
  log-sum-exp so nothing underflows. We marginalise; we do not cherry-pick
  the best fit.
- **Empirical Bayes (EB).** A prior is usually chosen up front; in
  *empirical* Bayes we instead **estimate the prior from the data
  itself** — here, the genome-wide rate of paralogs (§6).
- **Prior fraction `π`, estimated by EM.** `π` is that data-estimated
  paralog rate (the prior probability a random locus is a paralog).
  **EM (expectation–maximisation)** is a standard iterative recipe for
  estimating such a hidden fraction: guess `π`, compute each locus's
  paralog probability, average them to update `π`, repeat until it stops
  moving.
- **Sigmoid / logit.** A matched pair of functions. **logit** turns a
  probability into a log-odds (`log(π/(1−π))`); **sigmoid** (its inverse)
  squashes a log-odds back into a probability in `(0, 1)`. The posterior
  is `sigmoid(LR + logit π)` — i.e. add the prior log-odds to the LR, then
  convert to a probability.
- **FDR / q-value / local FDR.** **FDR (false discovery rate)** is the
  expected fraction of your *flagged* loci that are actually fine. The
  **q-value** is the FDR analogue of a p-value: flag everything with
  `q ≤ 0.01` and ~1% of what you removed was a real variant. **Local
  FDR** is the same idea for a single locus (`1 − posterior`). This is
  exactly the introgression-preservation knob: "of the calls I throw out,
  what fraction am I throwing out by mistake?" (§6).
- **Normal (Gaussian) likelihood, `σ`.** Modelling a continuous
  measurement (here, the relative copy number) as scattering around an
  expected value with a bell-curve spread. `σ` (sigma) is that spread; `σ₀`
  is the single-copy coverage spread (≈ 0.26–0.28 at a 500 bp window, fit per
  sample). The Gaussian is kept deliberately — the heavy coverage tail is
  duplication signal, not noise (§7).
- **Binomial likelihood.** The probability of seeing `alt_reads` ALT reads
  out of `total_reads` if each read independently carries the ALT with some
  probability (the VAF). The natural model for allele read counts.
- **Log-sum-exp (LSE).** A numerical-safety trick: to add probabilities
  that are stored as logs (which would otherwise underflow to zero), you
  exponentiate, sum, and re-log in a way that avoids the underflow. It is
  *how* we compute the averages above, not a modelling choice.
- **Robust estimator / median / MAD.** "Robust" means not thrown off by a
  minority of outliers — important here because paralog windows are
  exactly the outliers we must not let distort the fit. The **median**
  (middle value) is robust where the mean is not; **MAD** (median
  absolute deviation, ×1.4826) is the robust counterpart of the standard
  deviation.

### Software terms

- **Pure function.** A function whose output depends *only* on its
  arguments and which changes nothing outside itself — no files read or
  written, no global or hidden state, no randomness. Same inputs always
  give the same output. We build the statistics this way so each piece
  can be tested and parity-checked in isolation, and so results are
  reproducible run-to-run.
- **Pre-pass.** A first sweep over the data that computes per-sample
  summaries (here: the coverage curve and `Hobs`) *before* the main
  per-locus scoring pass that uses them.
- **Window / tile.** A fixed-width stretch of the genome (default 500 bp)
  over which coverage is averaged. "Tiled" windows are non-overlapping and
  anchored at genome coordinates, so each locus maps to exactly one tile.
- **Parity (parity check).** When an algorithm is reimplemented — here, the
  validated Python prototype ported to Rust — "parity" means the new code
  produces the *same numbers* as the trusted reference for the same
  inputs, either exactly or within a tiny rounding tolerance. A parity
  check is the test that proves the port is faithful (like a regression
  test, but against a reference implementation rather than the code's own
  past output).
- **Fixture / golden fixture.** A small file of fixed test data used by a
  test. A *golden* fixture also stores the reference's expected outputs, so
  the test can assert the new code reproduces them — what the parity check
  feeds in.

---

## 1. The problem, in one statement

> A genomic segment that is duplicated in the sample but **single-copy in the
> reference** collects the reads from all of its copies onto that one reference
> locus. Where the copies differ in sequence, this manufactures **artefactual
> heterozygotes** (and excess coverage). We must remove these false calls
> **without removing introgressions** — divergent alleles that are genuinely
> single-copy.

At the variant level the two are indistinguishable — both appear as a divergent
ALT allele — so the whole task is to find a signal that tells them apart.

| | hidden paralog (REMOVE) | introgression (KEEP) |
|---|---|---|
| copy number at locus | **>1** (collapsed copies) | 1 (orthologous) |
| coverage | **excess** | normal |
| ALT allele | from a divergent collapsed copy | a real, divergent allele |
| genotype | pseudo-heterozygote (forced) | real (can be homozygous-ALT) |

## 2. Core principle — coverage is the one safe signal

**Coverage is the only signal that flags a hidden paralog without ever risking a
real variant.** It both discriminates *and* is safe, and that combination is what
singles it out — here is the reasoning.

The two classes differ on exactly one axis: copy number. Extra copies bring extra
reads, so a hidden paralog shows **excess coverage**; an introgression is
single-copy and so has **normal coverage**. Of the other signals, only **sequence
divergence** is genuinely shared by the two classes — a paralog and an introgression
both carry a divergent ALT allele — so it cannot separate them. **Heterozygote
excess** is different: it is a real paralog signature, *not* shared, because only
the collapsed copies force pseudo-heterozygotes, whereas a single-copy introgression
segregates as ordinary genotypes. It still cannot be the *primary* filter, though —
other legitimate high-heterozygosity real loci show the same excess — so the model
keeps it secondary (first corollary below).

The tomato2 data make this concrete
([`second_look.py`](../../../benchmarks/tomato2/src/second_look.py)). Among
high-heterozygosity loci, the normal-coverage "introgression-like" class is
actually **more** divergent (mean MQDiff −12.5) than the excess-coverage paralog
class (−7.3). A divergence-only filter would therefore be worse than useless: at
any MQDiff threshold it clears a *larger* fraction of introgressions than of
paralogs — it would delete the introgressions we want to keep faster than the
paralogs it is meant to remove.

An introgression has normal depth,
so any statistic that fires only on excess coverage leaves it untouched.

Two consequences shaped everything below:

- **Heterozygote excess stays a secondary, diagnostic signal, never a primary
  one.** Legitimately high-heterozygosity loci exist — balancing selection, the
  S-locus. Heterozygote excess is a tell-tale only *in
  combination with* coverage excess, never on its own.
- **Reference-based mappability is the wrong tool here.** A collapse is *unique in
  the reference* (the duplicate copy was never assembled), so a k-mer mappability
  score rates it as fully mappable — it is blind to exactly the loci we are
  hunting. Mappability is also entangled with the copy-number signal itself:
  normalising it away erases real copy-number variants (Janevski 2012, GENSENG).
  We do not use it. (Low-mapping-quality reads are already filtered upstream,
  which would make a mappability covariate inconsistent anyway.)

## 3. The signals, and what each one can and cannot do

The model draws on four signals. Only the first is a primary discriminator; the
rest are supporting.

- **Coverage excess** — the primary, introgression-safe signal. It must be measured
  **per sample, GC-corrected, and in windows** (§4). Per-base coverage at ~6× depth
  has no power to tell a 2× carrier (~12 reads) from a single-copy sample that
  happens to be reading high.

- **Per-sample observed heterozygosity (`Hobs`)** — how inbred / self-fertilising
  each individual is. The model uses it (§5) because in a selfer a real variant is
  almost always homozygous, so a heterozygous carrier is anomalous — *unless that
  individual is heterozygous across the board*. It must be **per sample**: a panel
  mixes near-homozygous inbreds with more heterozygous wild relatives. The pileup
  emits `Hobs` (a single-sample observable); the inbreeding coefficient the model
  needs, `F = 1 − Hobs/Hexp`, is formed downstream with one cohort `Hexp` (§0), no
  extra pass.

  `Hobs` is a single-sample quantity, estimated during the Stage-1 pileup from a
  **rough per-site genotype** — no cohort needed. The rough call is *not* a VAF
  threshold, which would be depth-blind (at ~6× a true heterozygote reads a VAF
  anywhere from 1/6 to 5/6 by sampling alone). Instead it is a **depth-aware
  binomial likelihood ratio** over three models for a site's allele counts
  (`alt_reads` out of `total_reads`):

  - homozygous-REF: `alt_reads ~ Binomial(total_reads, ε)` (alt reads are errors)
  - heterozygote:   `alt_reads ~ Binomial(total_reads, ½)`
  - homozygous-ALT: `alt_reads ~ Binomial(total_reads, 1−ε)` (ε = per-read error rate)

  in **two steps**, each using the confidence margin `M`. First a **variant gate**:
  a site counts as a variant only if it beats hom-REF —
  `max(LL_het, LL_homalt) − LL_homref > M` — so a lone error-level alt read at high
  depth stays hom-REF and is not counted (this keeps the variant set clean). Then,
  among admitted variants, the **het-vs-hom split** on
  `logLR = LL_het − LL_homalt =
  total_reads·ln½ − alt_reads·ln(1−ε) − (total_reads − alt_reads)·ln ε` (the binomial
  coefficient cancels, so it is a couple of multiplies). The confidence
  margin `M` splits them three ways: **confident heterozygote**
  (`logLR > +M`), **confident homozygous-ALT** (`logLR < −M`), and **ambiguous**
  (`|logLR| ≤ M`). The observed heterozygosity is then the het **rate**,
  `Hobs = n_het / callable_positions` — het count over *all* callable positions (the
  covered-position total the coverage side of the same walk already has), **not**
  over the variant sites. (Dividing by variant sites, `n_het/(n_het+n_hom_alt)`,
  tracks reference divergence and inverts on divergent samples — measured, §7.)
  `n_ambiguous` is kept as an **uncertainty signal** — it grows at low depth, telling
  the model to trust that sample's `Hobs` less. The counts (`n_het`, `n_hom_alt`,
  `n_ambiguous`, `n_variant`), the callable-position total, and `ε`, `M` are stored
  in the `.psp`. `M` is a settled threshold, frozen at pileup time — unlike the
  coverage histogram, whose model is still being calibrated (§6, §9).

- **Divergence (MQDiff)** — mean mapping quality of ALT reads minus REF reads;
  negative means the ALT reads align worse, i.e. more divergent. It is **not used
  in the score.** Divergence is shared with introgression — a paralog and a real
  divergent allele both carry it (§2 shows the introgression-like class is, if
  anything, *more* divergent) — so it cannot separate the two without doing harm.
  It is also unreliable for indels: indel-carrying reads get low mapping quality
  from the aligner regardless of paralogy. We keep it out of the filter and emit it
  as a VCF INFO field instead, leaving it available for a future, separate rescue of
  the diverging-paralog blind spot (§7, §9).

- **Allele balance and the homozygous-ALT veto** — a collapse keeps the two
  orthologous REF copies, so its VAF is capped below 1. A **confident
  homozygous-ALT** sample (VAF ≈ 1) is therefore almost impossible under the paralog
  hypothesis, and is strong evidence for a real variant (§5.4, the hom-alt veto).

## 4. Measuring coverage, step by step

Per-base read counts are too noisy at 6× to tell one copy from two (the count distributions overlap), so we average depth over 500 bp windows to beat the noise down, measure it against each sample's own one-copy level, and remove the GC artifact.
The measurement has two phases: **fit a coverage model for the
sample once** (from data the Stage-1 pileup already stored in the `.psp`), then
**apply it at each locus** to express that locus's expected coverage *in copies*.

Throughout, a **window** (or **tile**) is a fixed `W`-bp stretch of the genome
(`W = 500 bp` by default), non-overlapping and anchored at genome coordinates, so
each locus falls in exactly one. A window's depth is its mean fragment depth —
post-filter, deduplicated, mate-overlap-resolved coverage read from the `.psp` (the
right depth for copy number, and better than raw BAM depth).

### Phase 1 — fit the sample's coverage model (once per sample)

The Stage-1 pileup walks every position of the sample and stores a compact summary
in the `.psp`: a 2-D histogram of windows, each binned by its GC fraction and its
mean depth. From that histogram we fit three quantities.

1. **The single-copy depth scale.** In the **post-filter, covered-base** coverage
   distribution, single-copy windows form the **dominant peak** — even when much of
   the genome is repetitive. Mapping-quality filtering strips the multi-mapping reads
   of interspersed repeats, so those windows lose coverage or drop out of the
   histogram entirely (rather than showing excess), and collapsed duplications — the
   excess we hunt — are only a small high tail. So the **peak (mode)** of the
   covered-window depth distribution marks the one-copy level: estimate
   `single_copy_scale` as that mode (the densest depth bin, or a smoothed
   kernel-density / half-sample-mode peak). We use the mode deliberately, not the
   mean or the median: the mean is dragged by both tails, and the median assumes
   single-copy windows are an outright *majority* (which can fail on a
   duplication-rich or polyploid genome), whereas the mode only needs single-copy to
   be the *most common* state — and finds it without us having to know how much of
   the genome is multi-copy.

2. **The GC-bias curve.** Library prep and sequencing systematically under- or
   over-cover GC-poor and GC-rich regions, by a factor that has nothing to do with
   copy number. To measure it, first express each window's depth relative to the
   scale, `rel = window_depth / single_copy_scale` (so `rel ≈ 1` for a single-copy
   window). Then, for each GC level (1%-wide GC bins), take the median `rel` over the
   windows in that bin — using **only** windows in the single-copy band
   `0.4 < rel < 1.6`, so duplicated/deleted windows don't distort the estimate — and
   smooth across neighbouring GC bins. The result is `gc_bias_curve(GC)`: the
   multiplier by which coverage runs high or low at that GC (say ×0.85 where GC-poor,
   ×1.1 mid-GC). It is fit **per sample, not pooled** — GC bias is a PCR/library
   property and the curves fan out even within one project (in
   [`build_gc_normalization.py`](../../../benchmarks/tomato2/src/build_gc_normalization.py),
   a pooled curve leaves a residual the per-sample curve removes).

3. **The single-copy spread `σ₀`.** Record how much a true single-copy window's
   `rel` scatters around 1 — the robust standard deviation `σ₀` (`1.4826·MAD`,
   ≈ 0.26–0.28 at `W = 500 bp`, **fit per sample from the data** rather than
   hardcoded). The likelihood (§5) needs it to know how much coverage noise is
   normal.

### Phase 2 — the expected coverage at a locus (per sample)

With the model fitted, turning a locus into a copy number is a two-line computation.

4. **Find the window and measure it.** The locus falls in one window; measure that
   window's **mean depth** (from the `.psp` body) and its **GC fraction** (from the
   reference).

5. **Expected single-copy depth.** Scale the one-copy level by the GC multiplier for
   this window's GC:

   > `expected_single_copy_depth = single_copy_scale · gc_bias_curve(GC)`

   This is the depth we would expect here *if the locus were single-copy in this
   sample*.

6. **Relative copy number.** Divide the observed window depth by that expectation:

   > `relative_copy_number = observed_window_depth / expected_single_copy_depth`

   `1.0` means one copy, `2.0` two copies, `0.5` a heterozygous deletion. This one
   number per sample — `relative_copy_number`, used in §5 — is the coverage signal
   the likelihood uses.

**Why `W = 500 bp`.** Smaller windows dilute a partial duplication less but give a
noisier single-copy mode; the count threshold filters noise more effectively than the
window size does, and 500 bp sits at a good operating point on tomato2 (single-copy
mode robust SD ≈ 0.23; a 2× carrier separates cleanly per window). **Refinement
owed:** GC bias actually acts at the *fragment* scale (Benjamini & Speed 2012), so a
GC window the size of the insert predicts coverage a little better and would tighten
the single-copy mode further (§9).

## 5. The per-locus statistic — a likelihood ratio of two hypotheses

At each locus, every sample `s` brings three numbers:

- `relative_copy_number` — the locus's coverage in copies (§4), `1.0` = one copy;
- `alt_reads`, `total_reads` — the **allele counts** (the AD): ALT reads out of total
  reads at the locus;
- `F` — the sample's **inbreeding coefficient** (how selfing/inbred the individual
  is), a fixed per-sample value estimated upstream.

The coverage and the allele counts are **observed data**. The *true* copy number and
genotype are **hidden** — they are encoded inside the two hypotheses below and
averaged over, never read off directly. (So `relative_copy_number` is our noisy
*measurement* of copy number, not the true count.)

We ask which of two stories better explains *all* the samples' data at this locus:
**H1**, a real single-copy variant, or **H2**, a hidden paralog. The verdict is a
likelihood ratio — how probable the data are under H2 versus under H1. The two
likelihoods are built the same way; they differ only in which genotypes or
configurations they allow, and what coverage and VAF each one implies.

### 5.1 The two per-sample likelihood factors

Both hypotheses score each sample with the same two independent factors:

- **A coverage factor** — how well the sample's `relative_copy_number` matches the
  copy number the hypothesis expects, as a bell curve (Normal) around an expected
  value `μ` with spread `σ`: `Normal(relative_copy_number; μ, σ)`. Under H1 the sample
  is single-copy, so `μ = 1` and `σ = σ₀`. Under H2 a carrier of `T` copies has
  `μ = T/2` and a wider spread `σ = σ₀·√(T/2)` (more copies, more sampling scatter).
- **An allele factor** — how well `alt_reads` out of `total_reads` matches the variant
  allele fraction (VAF) the hypothesis expects, as a binomial:
  `VAF^{alt_reads} · (1−VAF)^{total_reads − alt_reads}`. (The binomial coefficient is
  identical under both hypotheses, so it cancels in the ratio and we drop it.)

A sample's likelihood under any one genotype or configuration is just **coverage
factor × allele factor**. What differs between the hypotheses is the menu of
genotypes/configurations, their prior weights, and which unknowns we have to average
over. Everything below is computed in logarithms, and every "average over an unknown"
is done with log-sum-exp so nothing underflows.

### 5.2 H1 — a real, single-copy variant, step by step

Under H1 the locus has one unknown: the population allele frequency `p`. Every sample
is single-copy; its genotype comes from inbreeding-adjusted Hardy–Weinberg, and the
genotype fixes its VAF.

1. **Genotype prior.** Given `p` and the sample's `F`, the three genotype
   probabilities are
   - homozygous-REF: `(1−p)² + F·p·(1−p)`
   - heterozygous:   `2p(1−p)·(1−F)`
   - homozygous-ALT: `p² + F·p·(1−p)`.

   Inbreeding (large `F`) shifts weight off the heterozygote onto the two
   homozygotes — which is why a het call in a selfer is surprising.
2. **VAF per genotype.** homozygous-REF → `ε`; heterozygous → `½`; homozygous-ALT →
   `1−ε`, where `ε ≈ 0.01` is a small error floor that keeps the allele factor away
   from `log 0`.
3. **Likelihood of each genotype.** For each of the three, multiply its allele factor
   (the binomial at that VAF) by the coverage factor `Normal(relative_copy_number; 1, σ₀)`. The
   coverage factor is the same for all three — under H1 the sample is single-copy
   whatever its genotype — so it is a shared multiplier.
4. **Average over the sample's genotype.** We don't know the genotype, so the sample's
   contribution is the genotype-prior-weighted sum of its three genotype likelihoods.
5. **Combine the samples.** Given `p` the samples are independent, so multiply their
   contributions (add the logs).
6. **Average over the allele frequency.** We don't know `p`, so we average the
   whole product over it — but *not* with a flat prior. Real variants follow the
   population-genetic **site-frequency spectrum**: most are rare. So we weight `p`
   by the neutral **folded SFS**, proportional to `1/(p(1−p))`, over
   `p ∈ [1/2N, 1−1/2N]` (`N` = number of samples). The bounds are the frequency of
   a *single copy* — the finest the panel can resolve — so the low end needs no
   arbitrary cutoff. "Folded" because the ALT allele is defined against the
   reference, not polarised to ancestral, so we do not assume ALT is the derived
   allele. This gives H1 proper credit for the rare-variant explanation (measured:
   it makes the paralog rate robust to the low-frequency bound, which a flat grid
   is not — §7).

The result is `P(data | H1)`.

### 5.3 H2 — a hidden paralog, step by step

Under H2 the locus is a duplication collapsed in the reference. It has two unknowns:
the **duplication configuration** — how many copies `T`, and how many of them carry
the variant, `m` — and the **carrier frequency** `q`, how common the duplication is
in the panel. Each sample is independently a carrier or a non-carrier.

1. **Enumerate the configurations.** A carrier with `T` total copies, `m` of them
   mutated, reads at coverage `T/2` and VAF `m/T` — so more copies imply a *lower*
   VAF (a single mutated copy in a tandem array shows up as a low-VAF ALT). We
   enumerate `T ∈ {3, 4, 6, 8}` (capping coverage at 4×) and, per `T`, keep only two
   cases: a single mutated copy (`m = 1`, low VAF) and a balanced split (`m ≈ T/2`,
   VAF ≈ ½). The crowded in-between states are not separately identifiable, so we
   don't fit them.
2. **Carrier vs non-carrier prior.** Given the carrier frequency `q` and the sample's
   `F`, the probability the sample is a **non-carrier** is the inbreeding-adjusted
   "homozygous-absent" term `(1−q)² + F·q·(1−q)`; otherwise it is a **carrier**.
   (Same inbreeding logic as H1 — a selfer carrier is usually homozygous for the
   duplication.)
3. **Non-carrier likelihood.** A non-carrier is plain single-copy REF: coverage factor
   at `μ = 1`, allele factor at `VAF = ε`.
4. **Carrier likelihood, per configuration.** A carrier of configuration `(T, m)`:
   coverage factor at `μ = T/2` with spread `σ₀·√(T/2)`, allele factor at
   `VAF = m/T`.
5. **Average over the sample's carrier status.** For a given configuration and `q`,
   the sample's contribution is
   `P(non-carrier)·[non-carrier likelihood] + P(carrier)·[carrier likelihood]`.
6. **Combine the samples.** Given the configuration and `q`, the samples are
   independent, so multiply their contributions.
7. **Average over the unknowns.** Average that product over both the configurations
   and a grid of carrier frequencies `q` (≈40 points, rare to common), under a flat
   prior.

The result is `P(data | H2)`.

### 5.4 The homozygous-ALT veto (falls out for free)

Because every H2 carrier has `VAF = m/T < 1`, a sample that is confidently
homozygous-ALT (VAF ≈ 1) has almost no likelihood under H2. Such a sample therefore
pushes the locus strongly toward H1, a real variant — there is no special-case rule,
it is just the model working.

### 5.5 We average over the unknowns, we don't maximise

Both hypotheses are **integrated** over their free parameters — the allele frequency
for H1 (weighted by the folded site-frequency-spectrum prior, §5.2), the
configuration and carrier frequency for H2 (flat priors) — the final averaging step
of each derivation above. Maximising instead (picking each hypothesis's best-fitting parameters)
would reward H2's extra flexibility and inflate the apparent paralog rate; averaging
charges a hypothesis for that flexibility — the Occam penalty. On tomato2 it moved
the estimated paralog fraction from 12.3% to 11.1%, so most of the rate is genuine
flagging, not flexibility (§7).

### 5.6 The score

`LR = log P(data | H2) − log P(data | H1)`. Positive favours a hidden paralog,
negative a real variant — but the decision cut is **not** at 0 (§6). Nothing is
added to this ratio: it is a **pure likelihood ratio**, which is exactly what lets
§6 treat it as a Bayes factor and turn it into a calibrated probability. Divergence
(MQDiff) is deliberately kept out — it is shared with introgression (§2, §3) — and
is emitted as a VCF INFO field instead of entering the score.

## 6. From score to probability: the prior rate of hidden paralogs, and the cut

The likelihood ratio says which explanation fits a locus better, but not how
*likely* the locus is to be a paralog. For that we also need the **prior
probability that any locus is a hidden paralog** — how common they are genome-wide.
Because the score is a (marginal) log-likelihood ratio, the per-locus posterior is
just the score shifted by the prior's log-odds and squashed into a probability:

> **P(paralog | data) = sigmoid( LR + log(π / (1−π)) )**

where **π is the prior paralog rate, estimated from the data itself** by EM over the
genome-wide scores
([`build_paralog_eb.py`](../../../benchmarks/tomato2/src/build_paralog_eb.py)).

This is what fixes the cut, and the cut is **not at LR = 0**. Paralogs are rare a
priori, so the point where the posterior reaches 0.5 sits at a **positive** score,
`log((1−π)/π)`; a locus at LR = 0 carries only the base-rate posterior (~π). We set
the operating point by **false-discovery rate (q-value)**: high-confidence removal
at about FDR 1% (posterior ≳ 0.9), aggressive removal at FDR 10–20% (posterior ≳
0.3–0.5). Given the goal of preserving introgressions, the default is
high-confidence.

## 7. Key findings and decisions (the *why*)

- **Coverage is the only introgression-safe discriminator** (§2); confirmed on
  tomato2, where the introgression-like class is *more* divergent than the paralog
  class.
- **Power requires windowed, per-sample, GC-corrected coverage** (§4); per-base
  depth at 6× cannot resolve a 2× carrier.
- **The coverage likelihood stays Normal — the heavy coverage tail is duplication
  *signal*, not measurement noise.** The single-copy coverage looked over-dispersed
  (excess kurtosis ≈ 1.2), which tempts a heavy-tailed (Student-t) likelihood. But
  fitting the tail-heaviness per sample from the coverage distribution shows the
  single-copy *bulk* is Gaussian (fitted `ν → 200` inside the single-copy band); the
  heaviness lives only in the far tail, which is *real elevated-coverage windows* —
  duplications, CNVs, GC/mapping artefacts — i.e. the very thing the filter detects.
  A Student-t treats those excursions as noise to be discounted, so it blunts the
  primary discriminator: on tomato2 even a mild `ν = 15` cut the flagged set by a
  third (Jaccard 0.66), and a synthetic-heavy-tail sim only "favoured" the t because
  it *assumed* heavy-tailed noise the data does not have. We keep the Normal
  deliberately (§9).
- **Per sample, not cohort-wide** — both the GC curve and observed heterozygosity
  vary per individual (PCR/library bias, and inbreeding: on tomato2 the
  observed-het-derived `F` ranged 0.00–0.97).
- **It handles both segregating and fixed duplications.** A *segregating*
  duplication is found by cross-sample coverage variation (carriers in excess,
  others normal); a duplication *fixed in the whole cohort* has no such variation,
  so it is found instead by **absolute** coverage (against the sample's own
  single-copy level) plus the heterozygote-excess / read-ratio signature —
  Hardy–Weinberg gives an absolute anchor that a coverage-vs-cohort comparison
  lacks. The model uses absolute per-sample coverage, so it covers both.
- **Single-sample detection needs confident coverage excess.** A lone heterozygote
  at a confidently 2× sample is flagged (the within-sample link between coverage and
  allele); a lone heterozygote at *normal* coverage is kept (introgression-safe). At
  marginal excess (1.2–1.5×) a singleton is not confidently separable from noise.
- **The diverging-paralog blind spot** — paralogs whose divergent reads were
  mapping-quality-filtered show only mild coverage excess, so coverage alone can
  miss them (worst for a single, private carrier). We accept this as a known recall
  limitation rather than patch it with a divergence term in the score: MQDiff is
  shared with introgression and would put the introgressions we must keep at risk.
  The MQDiff INFO field leaves the door open to a separate, opt-in rescue once
  read-level ground truth (§9) shows whether the blind spot matters in practice.
- **The cut is FDR-calibrated and positive** (§6), never LR = 0. With divergence out
  of the score (§3), the LR is a *pure* likelihood ratio, so the posterior and
  q-values are an internally-valid Bayesian FDR; the empirical-Bayes π-estimate
  (concave, so its fixed point is unique and start-independent) and the tail-FDR add
  no calibration debt of their own. The one remaining calibration gap is the
  *absolute* anchor (next bullet) — an empirical validation, not a flaw in the
  machinery.
- **The absolute paralog rate is not yet anchored.** With the folded-SFS allele-
  frequency prior the data-estimated rate is π ≈ 9% (a flat prior gives ≈ 11%; the
  SFS is the more defensible model and the shift is toward *keeping* rare variants).
  It is plausible for a duplication-rich, self-fertilising genome and the flagged set
  looks paralog-like, but only real ground truth (read-level simulation, or known
  duplications) can confirm the rate (§9).
- **The score's ranking is robust to the model's grid/prior choices.** A sensitivity
  sweep — varying the carrier-frequency range, the copy-number cap, the configuration
  set, and the allele-frequency prior — leaves the ranked list of most-paralog-like
  loci ≥ 99% unchanged. The one fragile knob was the *flat* allele-frequency prior's
  low-frequency cutoff: pulling it in mislabels rare real variants as paralogs. The
  folded-SFS prior (§5.2) removes that fragility — the rate becomes robust to the
  bound because the spectrum's shape, not the cutoff, carries the weight.
- **Samples are treated as independent given the locus's frequency** (§5.2 step 5,
  §5.3 step 6) — the standard Hardy–Weinberg assumption. Real cohorts are related
  (families, F2s, clonal lines), which breaks strict independence. Two checks show
  this is safe to assume:
  1. *On the tomato2 panel*, de-duplicating related lines (one per genotype-
     correlation cluster) moves π < 1 point and keeps ≥ 92% of the flagged set.
  2. *Under simulated extreme structure* — an F2 (every locus segregates 1:2:1) and a
     6-founder × 10-clone cohort — the **real-variant false-flag rate is unchanged**
     (~0.1% at FDR 1%, same as an unrelated panel). The only cost of the most clonal
     case is reduced paralog *recall* (87% vs 99.7%), from the smaller effective
     sample size — the *safe* direction: it keeps a few paralog artefacts, it does not
     delete real variants.

  The reason is structural: coverage — the safe signal — is measured per sample
  against the sample's *own* single-copy level, so it does not depend on the cohort's
  relatedness at all. We state the assumption and rely on it. (If a real, extremely
  related dataset ever showed a problem, the fix is to down-weight correlated samples
  toward an effective sample size — not needed on the evidence so far.)

## 8. Intended production architecture

This is **not** a VCF-only, per-locus filter. Its introgression-safe core — the
coverage signal — needs the GC-corrected windowed coverage model, and that **cannot
be derived from the VCF** (the VCF's per-site `DP` is the weak, noisy signal we
abandoned). The honest shape is a **pre-pass plus a per-locus scoring pass**:

```
PRE-PASS  (per sample, over the .psp + reference FASTA)
    - per-sample observed heterozygosity (overall het rate)
    - depth-vs-GC curve + single-copy coverage scale
    - per-W-bp window depths  ->  relative copy number
        |
        v
PER-LOCUS SCORING
    relative copy number (pre-pass)  +  genotype/AD (per locus, from .psp/VCF)
        ->  H1-vs-H2 marginal LR (a pure likelihood ratio)  ->  posterior  ->  FDR cut
```

- **The `.psp` is already in the pipeline** (`pileup → .psp → var-calling`), and
  var-calling already reads every sample's `.psp`. So the pre-pass adds computation,
  not a new data dependency — the same pre-pass pattern as SSR parameter estimation
  and [contamination estimation](contamination_estimation.md).
- **Better still, fold the coverage model into Stage 1 (the pileup).** The pileup
  already walks every position of every sample, so emitting the per-window depth,
  the per-sample GC curve, and the observed heterozygosity there is nearly free, and
  removes the separate pre-pass.
- **The per-locus terms are cheap and VCF-derivable.** Apart from aggregating the
  observed heterozygosity, the per-locus inputs are just allele balance, the
  homozygous-ALT flag, and genotypes. (MQDiff is still written to the VCF as an INFO
  field, but it is not an input to the score.)

## 9. Open items / owed work

- **Ground truth to anchor the absolute FDR** — read-level simulation through the
  real pileup, or known tomato segmental duplications. (The within-model
  [`paralog_simulation.py`](../../../benchmarks/tomato2/src/paralog_simulation.py)
  validated the *mechanisms*, not the absolute rate.)
- **Fragment-scale GC** (Benjamini & Speed) in place of the analysis-window GC.
- **The diverging-paralog blind spot** — whether to add a separate, opt-in rescue
  (e.g. mild coverage excess *and* strongly negative MQDiff), kept out of the
  false-discovery-controlled core so it never risks introgressions. Gate the
  decision on read-level ground truth first: it is only worth building if the blind
  spot proves common.
- **A cheap per-locus production filter** (relative copy number + observed
  heterozygosity), trained against this likelihood-ratio model as the "teacher",
  with the introgression-safe and aggressive modes selectable.
- **Indels** — covered by the coverage axis (which is allele-agnostic); they use the
  same score as SNPs now that divergence is out of it.
- **Coverage likelihood stays Normal — a heavy-tailed / negative-binomial term is
  deliberately declined** (see §7 for the evidence). The one measured noise-side
  effect is that the spread at true 2× grows slightly faster than the `√` law
  (1.66× vs `√2 = 1.41×`); a fitted scaling exponent is a possible future tweak, but
  it only affects paralog recall (the safe direction) and never touches the tail.
- **Optional, longer-term:** the joint latent-factor (gCNV-style) formulation as the
  principled coverage model.

## 10. Empirical testbed (where the evidence lives)

All of the empirical work lives under
[`benchmarks/tomato2/src/`](../../../benchmarks/tomato2/src/):

- **Callset:** `build_psp.sh`, `run_varcalling.sh` (a max-retention callset).
- **Per-site features:** `build_feature_table.py`.
- **Windowed GC-corrected coverage:** `build_window_coverage.sh`,
  `assemble_window_table.py`, `build_gc_normalization.py`, and
  `build_window_locus_features.py`.
- **The model:** `build_paralog_lr.py` (the likelihood ratio) and
  `build_paralog_eb.py` (the prior estimate, posterior, and FDR).
- **Synthetic validation:** `paralog_extended_test.py` and `paralog_simulation.py`.
- **Interactive exploration:** `dashboard.py` (a marimo notebook over everything).

The deep-research synthesis on coverage covariates (the mappability confound,
fragment-scale GC, PCR vs PCR-free) is summarised in the project notes.
