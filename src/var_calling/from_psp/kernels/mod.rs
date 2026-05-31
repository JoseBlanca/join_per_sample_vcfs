//! Native column-native kernels — Phase A.1 replacement for the
//! Phase A.0 worker adapter that pipes column data through the
//! existing row-shape
//! [`PerGroupMerger`](crate::var_calling::per_group_merger) and
//! [`PosteriorEngine`](crate::var_calling::posterior_engine) kernels.
//!
//! The kernels split into four layers (each ported in its own
//! sub-step):
//!
//! 1. [`unify_alleles`] — column-native allele unification:
//!    per-position projection, chain-anchored compound detection,
//!    max-alleles cap. Phase A.1 layer 1.
//! 2. [`project_scalars`] — per-(sample, allele) `AlleleSupportStats`
//!    projection. Phase A.1 layer 2.
//! 3. [`compute_log_likelihoods`] — per-(sample, genotype) closed-form
//!    log-likelihood plus chain-broken fallback. Phase A.1 layer 3.
//! 4. *EM* (layer 4; queued).
//!
//! Each layer is byte-identity-tested against the existing kernel's
//! output on the same input as it lands.

pub mod compute_log_likelihoods;
pub mod project_scalars;
pub mod unify_alleles;
