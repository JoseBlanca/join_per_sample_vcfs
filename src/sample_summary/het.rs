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
//! hom-ref :  alt ~ Binomial(n, ε)       (alt reads are errors)
//! het     :  alt ~ Binomial(n, ½)
//! hom-alt :  alt ~ Binomial(n, 1 − ε)
//! ```
//!
//! and classifies in two steps, both using one confidence margin `M`
//! (nats). The binomial coefficient is identical in all three models and
//! cancels in every comparison, so each log-likelihood is a couple of
//! multiplies.
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

/// Per-site fragment counts for the het classifier: reference vs
/// non-reference (alt) observations at a candidate site, post-dedup /
/// mate-overlap. Taking ref and alt **separately** (rather than alt +
/// total) makes `alt ≤ total` hold by construction — there is no
/// representable "alt exceeds total" site — and removes the
/// argument-transposition risk of two same-typed counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SiteCounts {
    /// Reference-allele observations.
    pub ref_obs: u64,
    /// Non-reference (alt) observations, summed over all ALT alleles.
    pub alt_obs: u64,
}

impl SiteCounts {
    /// Total observations at the site (`ref_obs + alt_obs`), saturating at
    /// `u64::MAX` (unreachable for real fragment depths).
    pub fn total(&self) -> u64 {
        self.ref_obs.saturating_add(self.alt_obs)
    }
}

/// Parameters for the rough binomial-LR genotype call. The `pileup` CLI
/// supplies them (C1) and validates them at that boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HetClassifyParams {
    /// Minimum total fragment depth for a site to be scored at all.
    pub min_depth: u32,
    /// Per-read error rate `ε` for the hom models. In the open interval
    /// `(0, 1)` so `ln ε` and `ln(1 − ε)` are finite.
    pub error_rate: f64,
    /// Confidence margin `M` (nats) for both the variant gate and the
    /// het-vs-hom-alt split. `>= 0`.
    pub lr_margin: f64,
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
    // Precomputed logs of the binomial success probabilities.
    ln_eps: f64,
    ln_1m_eps: f64,
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
            params.error_rate.is_finite() && params.error_rate > 0.0 && params.error_rate < 1.0,
            "error_rate must be finite in (0, 1), got {}",
            params.error_rate,
        );
        assert!(
            params.lr_margin.is_finite() && params.lr_margin >= 0.0,
            "lr_margin must be finite and >= 0, got {}",
            params.lr_margin,
        );
        Self {
            ln_eps: params.error_rate.ln(),
            ln_1m_eps: (1.0 - params.error_rate).ln(),
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

        // `alt_obs <= total` by construction and the precomputed logs are
        // finite (constructor asserts), so every LL below is finite and
        // the `max` / `<` / `>` comparisons are total — no NaN can arise.
        let k = site.alt_obs as f64;
        let n = n_total as f64;
        let ref_obs = site.ref_obs as f64;

        let ll_hom_ref = k * self.ln_eps + ref_obs * self.ln_1m_eps;
        let ll_het = n * self.ln_half;
        let ll_hom_alt = k * self.ln_1m_eps + ref_obs * self.ln_eps;

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
            self.n_het += 1;
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::LN_10;

    fn params() -> HetClassifyParams {
        HetClassifyParams {
            min_depth: 4,
            error_rate: 0.02,
            lr_margin: LN_10, // 10:1 odds
        }
    }

    /// `(ref_obs, alt_obs)` site builder for readability.
    fn site(ref_obs: u64, alt_obs: u64) -> SiteCounts {
        SiteCounts { ref_obs, alt_obs }
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
    #[should_panic(expected = "error_rate must be finite in (0, 1)")]
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
