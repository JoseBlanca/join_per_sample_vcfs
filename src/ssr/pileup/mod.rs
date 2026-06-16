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
//! Present so far: [`count_repeats`] (the fast-path motif counter),
//! [`pair_hmm`] (the slow-path forward scorer + per-read candidate scoring),
//! [`candidate_generation`] (on-ladder rungs + off-ladder candidates), and
//! [`triage`] (the content pre-probe so far; read classification next).

pub mod candidate_generation;
pub mod count_repeats;
pub mod driver;
pub mod fetch_reads;
pub mod locus_record;
pub mod pair_hmm;
pub mod read_analysis;
pub mod triage;
