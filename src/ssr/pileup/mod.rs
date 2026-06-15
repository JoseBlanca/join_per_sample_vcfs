//! Stage 1 — `ssr-pileup`: per-sample evidence extraction.
//!
//! Turns one sample's BAM/CRAM + the catalog into one `.ssr.psp` evidence file
//! — a sparse, columnar, per-locus summary of what the reads say about each
//! locus's repeat length, with stutter deliberately left out (it is learned and
//! applied in Stage 2). Design: `doc/devel/architecture/ssr_pileup.md` and the
//! implementation sketch `doc/devel/implementation_plans/ssr_pileup.md`.
//!
//! **Build status.** Built bottom-up — the leaf scorers/counters first, the I/O
//! fetcher and driver last — so each layer is testable before it has a consumer.
//! Present so far: [`count_repeats`] (the fast-path motif counter).

pub mod count_repeats;
