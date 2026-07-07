//! The observed-heterozygosity accumulator.
//!
//! Estimates a sample's observed heterozygosity from a **rough per-site
//! genotype** computed in the Stage-1 walk (architecture Premise 1b) —
//! no cohort needed, because Hobs is a single-sample observable.
//!
//! Each candidate site is fed as a [`SiteCounts`] (reference vs
//! non-reference fragment observations, post-dedup / mate-overlap), so
//! `alt ≤ total` holds by construction. It scores three binomial models
//! for the alt count,
//!
//! ```text
//! hom-ref :  alt ~ Binomial(n, ε̂)       (alt reads are errors)
//! het     :  alt ~ Binomial(n, ½)
//! hom-alt :  alt ~ Binomial(n, 1 − ε̂)
//! ```
//!
//! and classifies in two steps, both using one confidence margin `M`
//! (nats). The binomial coefficient is identical in all three models and
//! cancels in every comparison, so each log-likelihood is a couple of
//! multiplies.
//!
//! **The error rate `ε̂` is estimated per site from the reads' own
//! qualities**, not a flat genome-wide constant (research report §3.1). Each
//! site carries `log_error_sum` — the sum, over every read at the site, of
//! that read's log-probability of being an error (base + alignment + mapping
//! quality combined; the record's per-allele `q_sum`). Its geometric mean,
//! `ε̂ = exp(log_error_sum / n)`, is the site's effective error rate, clamped
//! to `[error_rate, MAX_SITE_ERROR_RATE]`. The floor (`error_rate`, the old
//! flat value) makes this a conservative first step: `ε̂` can only rise above
//! the old constant on low-quality sites — so the hom-ref model only ever gets
//! *more* forgiving there, demoting noisy false hets, and never *less*
//! forgiving on clean sites, so no new het is invented. The ceiling keeps `ε̂`
//! below ½ so hom-ref stays distinguishable from het: an error-dominated pile
//! (e.g. all-MQ0 mismapped reads) saturates at the ceiling and is effectively
//! uncallable → stays hom-ref.
//!
//! 1. **Variant gate.** A site is a variant unless hom-ref is the clear
//!    winner: it is admitted iff `max(LL_het, LL_hom_alt) − LL_hom_ref > M`.
//!    (A lone alt read at high depth stays hom-ref and is *not* counted.)
//! 2. **Het vs hom-alt.** Among admitted variants, `logLR = LL_het −
//!    LL_hom_alt`: **confident het** (`> +M`), **confident hom-alt**
//!    (`< −M`), **ambiguous** (`|logLR| ≤ M`).
//!
//! This is depth-aware by construction: a clean 50/50 site scores a large
//! positive `logLR` that grows with `n`, while a borderline `(n−1)/n` site
//! sits inside the margin at low `n` (ambiguous) and only resolves as the
//! depth rises. `Hobs = n_het / (n_het + n_hom_alt)` is formed downstream;
//! `n_ambiguous` is the per-sample uncertainty weight.

use super::HetCounts;

/// Ceiling on the per-site effective error rate `ε̂`. Kept below `½` so the
/// hom-ref and het models stay distinguishable (at `ε̂ = ½` they coincide). A
/// site whose reads are error-dominated — e.g. a pile of MQ0 mismapped reads,
/// whose `log_error_sum` is ~0 so `ε̂ → 1` — saturates here and becomes
/// effectively uncallable, so it stays hom-ref rather than producing a false
/// het. The exact value is a guardrail, not a tuning knob: any site with an
/// effective error rate this high is already uncallable.
pub const MAX_SITE_ERROR_RATE: f64 = 0.4;

/// Default strand-bias veto threshold (research report §3.2): a confident het
/// whose ALT allele's forward-strand fraction differs from the REF allele's by
/// more than this many standard errors (a two-proportion z) is demoted to
/// ambiguous as a likely strand artifact. `z ≈ 3` ⇒ a two-sided p ≈ 0.003.
/// Strand confinement — unlike low MAPQ — has no biological confound
/// (introgressed haplotypes are not strand-biased), so this veto targets pure
/// artifacts. `f64::INFINITY` disables it.
pub const DEFAULT_STRAND_BIAS_Z: f64 = 3.0;

// ---------------------------------------------------------------------
// Considered and rejected: an ALT-vs-REF MAPQ-diff veto (report §3.3)
//
// We prototyped and measured a second artifact veto: a one-sided Welch's t
// that demoted a confident het whose ALT allele's mean MAPQ was systematically
// lower than REF's — the collapsed-paralog / mismapping signature used by
// bcftools MQBZ, GATK MQRankSum and VarScan max-mapqual-diff. It was REMOVED
// after measurement (this comment is the record of why, for anyone tempted to
// re-add it; the per-allele `mapq_sum` / `mapq_sum_sq` fields it needed were
// removed from `AlleleGroupStats` with it — the record still carries them on
// `AlleleSupportStats`).
//
// On the 63-sample tomato1 cohort it removed ~8% of all confident hets, but
// UNIFORMLY: the cleanest, most-inbred samples lost 11–13% of their few, real
// hets — as much as or more than the noisiest samples — and it did NOT move
// the one strongly het-inflated sample it was meant to catch (that sample's
// reads are high-quality, well-mapped and strand-balanced). The reason is
// biological: low ALT MAPQ is not only a paralog signal — a genuine
// INTROGRESSED haplotype also maps with lower MAPQ, because divergent-donor
// reads mismatch the reference more. The veto therefore removes real
// heterozygosity broadly rather than targeting artifacts, biasing obs_het (and
// so F) downward across every sample. This reproduced the downstream caller's
// own experience, where the paralog filter abandoned MAPQ-diff for a
// COVERAGE-based signal (excess per-sample coverage flags collapsed paralogs
// without penalising divergence).
//
// The strand veto (below) has no such confound — introgressions are not
// strand-biased — and measured as cleanly targeted (noisy samples down, clean
// samples untouched), so it is kept. A het-inflated sample that survives the
// ε̂ and strand steps is not reachable by any cheap read-quality heuristic
// here; catching it is the coverage-based paralog signal's job, not this
// rough caller's.
// ---------------------------------------------------------------------

/// Aggregated support for one allele group at a site — the reference allele,
/// or all ALT alleles pooled. `obs` is the fragment count; the remaining
/// fields are the per-allele [`AlleleSupportStats`] moments summed over the
/// group. Splitting REF vs ALT this way gives the strand veto its REF-vs-ALT
/// contrast directly.
///
/// [`AlleleSupportStats`]: crate::pileup_record::AlleleSupportStats
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct AlleleGroupStats {
    /// Fragment observations supporting this group (post-dedup / mate-overlap).
    pub obs: u64,
    /// Forward-strand observations among `obs` (reverse = `obs − fwd`). The
    /// 2×2 `(ref.fwd, ref.rev, alt.fwd, alt.rev)` table drives the strand veto.
    pub fwd: u64,
    /// Σ per-read log P(error) over the group (the record's `q_sum`; `≤ 0`).
    /// REF + ALT sum to the site's `log_error_sum` for `ε̂`.
    pub log_error_sum: f64,
}

/// Per-site evidence for the het classifier: the reference allele group vs the
/// pooled ALT group. Splitting REF and ALT (rather than alt + total) makes
/// `alt ≤ total` hold by construction — there is no representable "alt exceeds
/// total" site — and supplies the vetoes their REF-vs-ALT contrast.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SiteCounts {
    /// The reference allele's support.
    pub reference: AlleleGroupStats,
    /// All ALT alleles pooled.
    pub alt: AlleleGroupStats,
}

impl SiteCounts {
    /// Total observations at the site (`ref + alt`), saturating at `u64::MAX`
    /// (unreachable for real fragment depths).
    pub fn total(&self) -> u64 {
        self.reference.obs.saturating_add(self.alt.obs)
    }

    /// Σ per-read log P(error) over **every** read at the site — the record's
    /// per-allele `q_sum` summed across REF + all ALT alleles. Always `≤ 0`.
    /// Divided by the total depth and exponentiated it yields the site's
    /// effective error rate `ε̂ = exp(log_error_sum / total)`, the
    /// quality-aware replacement for the flat `ε` in the hom models.
    pub fn log_error_sum(&self) -> f64 {
        self.reference.log_error_sum + self.alt.log_error_sum
    }
}

/// Parameters for the rough binomial-LR genotype call. The `pileup` CLI
/// supplies them (C1) and validates them at that boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HetClassifyParams {
    /// Minimum total fragment depth for a site to be scored at all.
    pub min_depth: u32,
    /// **Floor** on the per-site error rate `ε̂` for the hom models. The
    /// effective `ε̂` is estimated from the reads' qualities (see the module
    /// docs) and never allowed below this value — so a site is never judged
    /// cleaner than this genome-wide baseline. In the open interval
    /// `(0, MAX_SITE_ERROR_RATE)` so the clamp range is non-empty and
    /// `ln ε̂`, `ln(1 − ε̂)` are finite.
    pub error_rate: f64,
    /// Confidence margin `M` (nats) for both the variant gate and the
    /// het-vs-hom-alt split. `>= 0`.
    pub lr_margin: f64,
    /// Strand-bias veto threshold (see [`DEFAULT_STRAND_BIAS_Z`]). A confident
    /// het whose ALT strand fraction deviates from REF's by more than this
    /// many standard errors is demoted to ambiguous. `f64::INFINITY` disables
    /// the veto; must be `> 0` and not `NaN`.
    pub strand_bias_z: f64,
}

/// Accumulates the observed-heterozygosity counts across one sample's
/// candidate-variant sites. Build with [`new`], feed each site with
/// [`observe_site`], and reduce with [`finish`].
///
/// [`new`]: Self::new
/// [`observe_site`]: Self::observe_site
/// [`finish`]: Self::finish
#[derive(Debug, Clone)]
pub struct HetAccumulator {
    params: HetClassifyParams,
    // `ln ½`, the het model's per-read log-likelihood; precomputed once. The
    // hom models' `ln ε̂` / `ln(1 − ε̂)` are per-site (see `observe_site`), so
    // they cannot be precomputed here.
    ln_half: f64,
    n_het: u64,
    n_hom_alt: u64,
    n_ambiguous: u64,
}

impl HetAccumulator {
    /// Construct an empty accumulator.
    ///
    /// # Panics
    ///
    /// Panics (debug and release) if `error_rate` is not finite in
    /// `(0, 1)` or `lr_margin` is not finite and `>= 0`. These are
    /// CLI-validated config (C1); an invalid value here is a programmer
    /// error and fails loudly rather than producing `NaN`/`±∞` scores.
    pub fn new(params: HetClassifyParams) -> Self {
        assert!(
            params.error_rate.is_finite()
                && params.error_rate > 0.0
                && params.error_rate < MAX_SITE_ERROR_RATE,
            "error_rate (ε̂ floor) must be finite in (0, {MAX_SITE_ERROR_RATE}), got {}",
            params.error_rate,
        );
        assert!(
            params.lr_margin.is_finite() && params.lr_margin >= 0.0,
            "lr_margin must be finite and >= 0, got {}",
            params.lr_margin,
        );
        assert!(
            !params.strand_bias_z.is_nan() && params.strand_bias_z > 0.0,
            "strand_bias_z must be > 0 (or +inf to disable), got {}",
            params.strand_bias_z,
        );
        Self {
            ln_half: 0.5_f64.ln(),
            params,
            n_het: 0,
            n_hom_alt: 0,
            n_ambiguous: 0,
        }
    }

    /// Score one candidate site. Sites below `min_depth`, and sites the
    /// variant gate rejects as hom-ref, are not counted.
    pub fn observe_site(&mut self, site: SiteCounts) {
        let n_total = site.total();
        if n_total < u64::from(self.params.min_depth) {
            return;
        }

        // Note — min-alt-count (report §3.4, freebayes `minAltCount 2` /
        // VarScan `min-reads2 2`): considered and NOT added as an explicit
        // floor, because the depth-aware gate below already subsumes it. A lone
        // ALT read never passes the variant gate at `min_depth = 4` — the
        // thinnest such site, 3 ref / 1 alt at n = 4, scores 1.20 < ln 10 and
        // stays hom-ref (see `lone_alt_read_is_never_a_variant`). Adding an
        // explicit `alt >= 2` floor here was measured byte-identical on tomato.
        // Add one only if `min_depth` / the margin are ever loosened.

        // `alt ≤ total` by construction. The per-site `ε̂` is clamped to
        // `[error_rate, MAX_SITE_ERROR_RATE] ⊂ (0, ½)`, so `ln ε̂` and
        // `ln(1 − ε̂)` are finite and every LL below is finite — the `max` /
        // `<` / `>` comparisons are total and no NaN can arise.
        let k = site.alt.obs as f64;
        let n = n_total as f64;
        let ref_obs = site.reference.obs as f64;

        // Site-level effective error rate from the reads' own qualities
        // (research report §3.1, Option A). `log_error_sum ≤ 0`, so the
        // geometric mean `exp(log_error_sum / n) ∈ (0, 1]`; the clamp pins it
        // into the valid model range. `n > 0` here (min-depth guard passed).
        let eps = (site.log_error_sum() / n)
            .exp()
            .clamp(self.params.error_rate, MAX_SITE_ERROR_RATE);
        let ln_eps = eps.ln();
        let ln_1m_eps = (1.0 - eps).ln();

        let ll_hom_ref = k * ln_eps + ref_obs * ln_1m_eps;
        let ll_het = n * self.ln_half;
        let ll_hom_alt = k * ln_1m_eps + ref_obs * ln_eps;

        // Variant gate: admit unless hom-ref is the clear winner.
        let best_variant = ll_het.max(ll_hom_alt);
        if best_variant - ll_hom_ref <= self.params.lr_margin {
            return; // hom-ref / not a confident variant
        }

        // Het vs hom-alt among admitted variants. The counters are bounded
        // by the number of candidate sites in a genome (<< u64::MAX), so
        // `+= 1` and the `finish` sum cannot overflow.
        let log_lr = ll_het - ll_hom_alt;
        if log_lr > self.params.lr_margin {
            // Tentative confident het. The strand veto runs here — a het whose
            // ALT allele is strand-confined relative to REF (§3.2) is demoted
            // to ambiguous (not counted as het), the cheap false-het rejection
            // that `obs_het` cares about. Only hets are vetoed; hom-alt and
            // ambiguous are untouched. (A MAPQ-diff veto was tried here and
            // rejected — see the note above `AlleleGroupStats`.)
            if strand_biased(&site, self.params.strand_bias_z) {
                self.n_ambiguous += 1;
            } else {
                self.n_het += 1;
            }
        } else if log_lr < -self.params.lr_margin {
            self.n_hom_alt += 1;
        } else {
            self.n_ambiguous += 1;
        }
    }

    /// Reduce to the per-sample [`HetCounts`].
    pub fn finish(self) -> HetCounts {
        let n_variant = self.n_het + self.n_hom_alt + self.n_ambiguous;
        HetCounts {
            n_het_sites: self.n_het,
            n_hom_alt_sites: self.n_hom_alt,
            n_ambiguous_sites: self.n_ambiguous,
            n_variant_sites: n_variant,
            min_depth: self.params.min_depth,
            error_rate: self.params.error_rate,
            lr_margin: self.params.lr_margin,
        }
    }
}

/// Strand-bias veto (research report §3.2). Tests whether the ALT allele's
/// forward-strand fraction differs from the REF allele's, via a two-proportion
/// z on the `(ref_fwd, ref_rev, alt_fwd, alt_rev)` table. Returns `true` (→
/// demote the het) when `|z| > threshold`.
///
/// Comparing ALT *to REF* (rather than to a flat ½) accounts for sites whose
/// coverage is itself strand-skewed. The test only fires when every **margin**
/// is `≥ 2` (bcftools `bam2bcf.c:1361`): it needs reads on both strands and
/// both alleles to be meaningful, and this stops it firing on thin counts. The
/// z is inherently depth-aware — a strand-confined ALT at low depth may be
/// chance and is not vetoed until the deviation is significant.
fn strand_biased(site: &SiteCounts, threshold: f64) -> bool {
    if threshold.is_infinite() {
        return false; // veto disabled
    }
    let (nr, na) = (site.reference.obs, site.alt.obs);
    let (rf, af) = (site.reference.fwd, site.alt.fwd);
    // Reverse-strand cells; `fwd ≤ obs` by construction so these don't wrap.
    let (rr, ar) = (nr - rf, na - af);
    let fwd_total = rf + af;
    let rev_total = rr + ar;
    // Margin guard: both alleles present, both strands present.
    if nr < 2 || na < 2 || fwd_total < 2 || rev_total < 2 {
        return false;
    }
    let (nr, na) = (nr as f64, na as f64);
    let p_alt = af as f64 / na; // ALT forward fraction
    let p_ref = rf as f64 / nr; // REF forward fraction
    let p_pooled = fwd_total as f64 / (nr + na);
    let variance = p_pooled * (1.0 - p_pooled) * (1.0 / na + 1.0 / nr);
    if variance <= 0.0 {
        // Both alleles fully on one (the same) strand: no ALT-vs-REF strand
        // contrast to speak of — not a strand artifact of the ALT allele.
        return false;
    }
    let z = (p_alt - p_ref) / variance.sqrt();
    z.abs() > threshold
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::LN_10;

    fn params() -> HetClassifyParams {
        HetClassifyParams {
            min_depth: 4,
            error_rate: 0.02,
            lr_margin: LN_10, // 10:1 odds
            strand_bias_z: DEFAULT_STRAND_BIAS_Z,
        }
    }

    /// One allele group at `per_read_eps` quality with **balanced strands**
    /// (`fwd = obs/2`) — so the strand veto never fires unless a test sets the
    /// strands itself.
    fn group(obs: u64, per_read_eps: f64) -> AlleleGroupStats {
        AlleleGroupStats {
            obs,
            fwd: obs / 2,
            log_error_sum: obs as f64 * per_read_eps.ln(),
        }
    }

    /// `(ref_obs, alt_obs)` site builder for readability. Encodes a per-read
    /// error rate of exactly `0.02` (the classifier's floor) with balanced
    /// strands, so `ε̂` clamps to `0.02` and no veto fires — these tests pin the
    /// *same* classifications the old flat-`ε` model produced. The
    /// quality-aware path is exercised by [`site_q`], the strand path by
    /// [`site_strand`].
    fn site(ref_obs: u64, alt_obs: u64) -> SiteCounts {
        site_q(ref_obs, alt_obs, 0.02)
    }

    /// Site builder with an explicit per-read error rate `per_read_eps` (every
    /// read contributes `ln(per_read_eps)` ⇒ `ε̂ = clamp(per_read_eps, ..)`),
    /// balanced strands.
    fn site_q(ref_obs: u64, alt_obs: u64, per_read_eps: f64) -> SiteCounts {
        SiteCounts {
            reference: group(ref_obs, per_read_eps),
            alt: group(alt_obs, per_read_eps),
        }
    }

    /// Site builder with explicit forward-strand counts per allele group (for
    /// the strand veto), clean quality (`ε̂ = 0.02`) and uniform MAPQ.
    fn site_strand(ref_obs: u64, ref_fwd: u64, alt_obs: u64, alt_fwd: u64) -> SiteCounts {
        SiteCounts {
            reference: AlleleGroupStats {
                fwd: ref_fwd,
                ..group(ref_obs, 0.02)
            },
            alt: AlleleGroupStats {
                fwd: alt_fwd,
                ..group(alt_obs, 0.02)
            },
        }
    }

    /// A clean 50/50 site at good depth is a confident het.
    #[test]
    fn clean_balanced_site_is_het() {
        let mut acc = HetAccumulator::new(params());
        acc.observe_site(site(15, 15)); // n = 30
        let h = acc.finish();
        assert_eq!(
            (h.n_het_sites, h.n_hom_alt_sites, h.n_ambiguous_sites),
            (1, 0, 0)
        );
        assert_eq!(h.n_variant_sites, 1);
    }

    /// An all-alt site is a confident hom-alt.
    #[test]
    fn all_alt_site_is_hom_alt() {
        let mut acc = HetAccumulator::new(params());
        acc.observe_site(site(0, 30));
        let h = acc.finish();
        assert_eq!(
            (h.n_het_sites, h.n_hom_alt_sites, h.n_ambiguous_sites),
            (0, 1, 0)
        );
    }

    /// A site with no alt (or only error-level alt) is hom-ref — not a
    /// variant, not counted.
    #[test]
    fn hom_ref_sites_are_not_counted() {
        let mut acc = HetAccumulator::new(params());
        acc.observe_site(site(30, 0)); // no alt
        acc.observe_site(site(29, 1)); // one alt — consistent with ε error
        let h = acc.finish();
        assert_eq!(h.n_variant_sites, 0);
    }

    /// Depth-awareness: the *same* near-hom-alt VAF (n−1 of n alt) is
    /// ambiguous at low depth and resolves to hom-alt as depth rises.
    #[test]
    fn near_hom_alt_is_ambiguous_at_low_depth_resolves_high() {
        let mut low = HetAccumulator::new(params());
        low.observe_site(site(1, 5)); // 5/6 alt at 6x
        let lh = low.finish();
        assert_eq!(lh.n_ambiguous_sites, 1, "5/6 at 6x should be ambiguous");
        assert_eq!(lh.n_het_sites + lh.n_hom_alt_sites, 0);

        let mut high = HetAccumulator::new(params());
        high.observe_site(site(1, 59)); // 59/60 alt at 60x
        let hh = high.finish();
        assert_eq!(hh.n_hom_alt_sites, 1, "59/60 at 60x should be hom-alt");
    }

    /// Sites below `min_depth` are skipped entirely.
    #[test]
    fn below_min_depth_is_skipped() {
        let mut acc = HetAccumulator::new(params()); // min_depth 4
        acc.observe_site(site(2, 1)); // n = 3
        acc.observe_site(site(0, 2)); // n = 2
        let h = acc.finish();
        assert_eq!(h.n_variant_sites, 0);
    }

    /// The depth-aware gate inherently enforces min-alt-count (§3.4): a single
    /// ALT read is never a variant — at the thinnest admitted depth (3 ref /
    /// 1 alt at `n = 4`) and at high depth (29 / 1, 99 / 1). This is why no
    /// explicit `alt >= 2` floor is needed (see the note in `observe_site`).
    #[test]
    fn lone_alt_read_is_never_a_variant() {
        for (r, a) in [(3, 1), (29, 1), (99, 1)] {
            let mut acc = HetAccumulator::new(params());
            acc.observe_site(site(r, a));
            assert_eq!(
                acc.finish().n_variant_sites,
                0,
                "{r} ref / {a} alt must not be a variant"
            );
        }
    }

    /// A site at exactly `min_depth` is admitted to scoring (the guard is
    /// `n_total < min_depth`, so equality passes).
    #[test]
    fn site_at_exactly_min_depth_is_scored() {
        let mut acc = HetAccumulator::new(params()); // min_depth 4
        acc.observe_site(site(2, 2)); // n == 4, balanced -> het
        assert_eq!(acc.finish().n_variant_sites, 1);
    }

    /// A pure-alt site at exactly `min_depth` classifies as hom-alt
    /// (`ref_obs == 0` arithmetic edge; the regime where collapsed-paralog
    /// artefacts live).
    #[test]
    fn pure_hom_alt_at_min_depth_is_hom_alt() {
        let mut acc = HetAccumulator::new(params());
        acc.observe_site(site(0, 4)); // n == 4, all alt
        let h = acc.finish();
        assert_eq!(
            (h.n_het_sites, h.n_hom_alt_sites, h.n_ambiguous_sites),
            (0, 1, 0)
        );
    }

    // ---- Quality-aware ε̂ (research report §3.1, Option A) ----

    /// A site whose counts are a confident het at clean quality is **demoted
    /// to hom-ref (not counted)** when the reads are junky enough that `ε̂`
    /// rises well above the flat floor. Same `(ref, alt)`, only the quality
    /// differs — this is the core false-het rejection lever.
    #[test]
    fn junky_reads_demote_het_to_hom_ref() {
        // 20 ref / 10 alt (33% VAF) is a confident het at ε̂ = 0.02 ...
        let mut clean = HetAccumulator::new(params());
        clean.observe_site(site_q(20, 10, 0.02));
        assert_eq!(clean.finish().n_het_sites, 1, "clean 20/10 must be het");

        // ... but at ε̂ = 0.20 (e.g. low-MAPQ mismapped reads) the same counts
        // are consistent with error and the variant gate rejects them.
        let mut junky = HetAccumulator::new(params());
        junky.observe_site(site_q(20, 10, 0.20));
        assert_eq!(
            junky.finish().n_variant_sites,
            0,
            "junky 20/10 must be demoted to hom-ref"
        );
    }

    /// The demotion is monotone in read quality: as `ε̂` climbs, a fixed het
    /// site only ever loses confidence (het → … → not counted), never regains
    /// it. Swept on the 20/10 site whose flip lies between 0.10 and 0.20.
    #[test]
    fn het_demotion_is_monotone_in_error_rate() {
        let het_at = |eps: f64| {
            let mut acc = HetAccumulator::new(params());
            acc.observe_site(site_q(20, 10, eps));
            acc.finish().n_het_sites
        };
        assert_eq!(het_at(0.02), 1, "0.02 het");
        assert_eq!(het_at(0.10), 1, "0.10 still het");
        assert_eq!(het_at(0.20), 0, "0.20 demoted");
        assert_eq!(het_at(0.35), 0, "0.35 demoted");
    }

    /// The `ε̂` floor: reads cleaner than the floor (`error_rate = 0.02`) are
    /// pinned at the floor, so an ultra-high-quality site classifies exactly
    /// as the old flat-`ε` model would — the change never *invents* a het.
    #[test]
    fn ultra_clean_reads_floor_to_baseline() {
        let mut floored = HetAccumulator::new(params());
        floored.observe_site(site_q(15, 15, 1e-9)); // ε̂ clamps up to 0.02
        let mut baseline = HetAccumulator::new(params());
        baseline.observe_site(site_q(15, 15, 0.02));
        assert_eq!(floored.finish().n_het_sites, baseline.finish().n_het_sites);
    }

    /// An error-dominated pile — every read MAPQ0, so `log_error_sum ≈ 0` and
    /// the raw `ε̂ → 1` — saturates at [`MAX_SITE_ERROR_RATE`] and stays
    /// hom-ref rather than producing a false het.
    #[test]
    fn all_mq0_site_saturates_and_stays_hom_ref() {
        let mut acc = HetAccumulator::new(params());
        // per_read_eps = 1.0 ⇒ log_error_sum = 0 ⇒ raw ε̂ = exp(0) = 1 ⇒
        // clamped to 0.4.
        acc.observe_site(site_q(20, 10, 1.0));
        assert_eq!(
            acc.finish().n_variant_sites,
            0,
            "MQ0-saturated site must stay hom-ref"
        );
    }

    // ---- Strand-bias veto (research report §3.2) ----

    /// A confident het whose ALT allele is strand-confined (all forward) while
    /// REF sits on both strands is demoted from het to ambiguous.
    #[test]
    fn strand_confined_alt_demotes_het_to_ambiguous() {
        let mut acc = HetAccumulator::new(params());
        // 15/15 (a clean het) but every ALT read forward, REF balanced.
        acc.observe_site(site_strand(15, 7, 15, 15));
        let h = acc.finish();
        assert_eq!(
            (h.n_het_sites, h.n_hom_alt_sites, h.n_ambiguous_sites),
            (0, 0, 1),
            "strand-confined ALT must demote het → ambiguous"
        );
    }

    /// The same balanced-count het with ALT on **both** strands (matching REF)
    /// is kept — the veto only targets strand skew, not heterozygosity.
    #[test]
    fn strand_balanced_het_is_kept() {
        let mut acc = HetAccumulator::new(params());
        acc.observe_site(site_strand(15, 7, 15, 7));
        assert_eq!(acc.finish().n_het_sites, 1, "balanced-strand het must stay");
    }

    /// Disabling the veto (`strand_bias_z = +∞`) keeps even a strand-confined
    /// het — the veto is the only thing that would have demoted it.
    #[test]
    fn strand_veto_disabled_keeps_confined_het() {
        let mut p = params();
        p.strand_bias_z = f64::INFINITY;
        let mut acc = HetAccumulator::new(p);
        acc.observe_site(site_strand(15, 7, 15, 15));
        assert_eq!(acc.finish().n_het_sites, 1, "disabled veto must keep het");
    }

    /// The margin guard: a strand-confined het too thin on one strand (here no
    /// reverse reads at all → `rev_total < 2`) is **not** vetoed — the test
    /// needs both strands present to be meaningful.
    #[test]
    fn strand_veto_skips_when_a_strand_margin_is_thin() {
        let mut acc = HetAccumulator::new(params());
        // All reads forward (rev_total = 0): can't judge strand bias → kept.
        acc.observe_site(site_strand(15, 15, 15, 15));
        assert_eq!(
            acc.finish().n_het_sites,
            1,
            "no reverse reads ⇒ guard skips veto ⇒ het kept"
        );
    }

    /// The variant gate uses a strict `>`: a site whose gate value equals
    /// the margin exactly is rejected (not counted). A balanced site has
    /// `gate == split logLR`, so setting the margin to that value rejects
    /// at the gate.
    #[test]
    fn site_at_exact_gate_margin_is_rejected() {
        let (ref_obs, alt_obs) = (15.0_f64, 15.0_f64);
        let eps = 0.02_f64;
        let n = ref_obs + alt_obs;
        let ll_het = n * 0.5_f64.ln();
        let ll_hom_alt = alt_obs * (1.0 - eps).ln() + ref_obs * eps.ln();
        let ll_hom_ref = alt_obs * eps.ln() + ref_obs * (1.0 - eps).ln();
        let gate = ll_het.max(ll_hom_alt) - ll_hom_ref;

        let mut p = params();
        p.error_rate = eps;
        p.lr_margin = gate; // gate == M -> `> M` false -> reject
        let mut acc = HetAccumulator::new(p);
        acc.observe_site(site(15, 15));
        assert_eq!(acc.finish().n_variant_sites, 0, "gate == M must reject");
    }

    /// The het/hom split uses a strict `>`: an admitted site whose
    /// `logLR == +M` exactly is ambiguous, not confident het. A mildly
    /// alt-biased site has `gate > split`, so the margin can equal the
    /// split while the gate still admits.
    #[test]
    fn site_at_exact_split_margin_is_ambiguous() {
        let (ref_obs, alt_obs) = (10.0_f64, 20.0_f64);
        let eps = 0.02_f64;
        let n = ref_obs + alt_obs;
        let ll_het = n * 0.5_f64.ln();
        let ll_hom_alt = alt_obs * (1.0 - eps).ln() + ref_obs * eps.ln();
        let split = ll_het - ll_hom_alt; // positive, gate is larger

        let mut p = params();
        p.error_rate = eps;
        p.lr_margin = split; // logLR == M -> `> M` false -> ambiguous
        let mut acc = HetAccumulator::new(p);
        acc.observe_site(site(10, 20));
        let h = acc.finish();
        assert_eq!(
            (h.n_het_sites, h.n_hom_alt_sites, h.n_ambiguous_sites),
            (0, 0, 1),
            "logLR == +M must be ambiguous, not confident het"
        );
    }

    /// The finished counts satisfy the B1 `HetCounts` invariants
    /// (sum identity, error-rate range) via `SampleSummary::validate`.
    #[test]
    fn finished_counts_validate() {
        let mut acc = HetAccumulator::new(params());
        acc.observe_site(site(15, 15)); // het
        acc.observe_site(site(0, 30)); // hom-alt
        acc.observe_site(site(1, 5)); // ambiguous
        let het = acc.finish();
        assert_eq!(het.n_variant_sites, 3);

        let summary = super::super::SampleSummary {
            version: super::super::SAMPLE_SUMMARY_VERSION,
            coverage_by_gc: super::super::CoverageByGcHistogram {
                window_bp: 500,
                gc_bins: 1,
                depth_bin_width: 1.0,
                depth_bins: 1,
                n_positions: 0,
                n_skipped_tiles: 0,
                callable_positions: 0,
                counts: vec![0, 0],
            },
            heterozygosity: het,
        };
        summary
            .validate()
            .expect("het accumulator output must validate");
    }

    /// An empty accumulator finishes to all-zero counts.
    #[test]
    fn empty_accumulator_finishes_zero() {
        let h = HetAccumulator::new(params()).finish();
        assert_eq!(h.n_variant_sites, 0);
        assert_eq!(
            (h.n_het_sites, h.n_hom_alt_sites, h.n_ambiguous_sites),
            (0, 0, 0)
        );
    }

    /// A non-finite / out-of-range error rate is a programmer error and
    /// panics at construction.
    #[test]
    #[should_panic(expected = "error_rate (ε̂ floor) must be finite")]
    fn invalid_error_rate_panics() {
        let mut p = params();
        p.error_rate = 0.0;
        let _ = HetAccumulator::new(p);
    }

    proptest::proptest! {
        /// Every observed site increments exactly one of the three class
        /// counts (or none, if skipped/hom-ref); the four-count identity
        /// `n_het + n_hom_alt + n_ambiguous == n_variant` holds for
        /// arbitrary `(ref_obs, alt_obs)` sequences. `SiteCounts` makes
        /// `alt <= total` unrepresentable, so no input is invalid.
        #[test]
        fn accumulator_maintains_count_identity(
            sites in proptest::collection::vec((0u64..200, 0u64..200), 0..100),
        ) {
            let mut acc = HetAccumulator::new(params());
            for (ref_obs, alt_obs) in sites {
                acc.observe_site(site(ref_obs, alt_obs));
            }
            let h = acc.finish();
            proptest::prop_assert_eq!(
                h.n_het_sites + h.n_hom_alt_sites + h.n_ambiguous_sites,
                h.n_variant_sites
            );
        }
    }
}
