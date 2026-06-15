//! Stage 0 — `ssr-catalog`: build the per-genome SSR locus catalog.
//!
//! The catalog is the SSR pipeline's first artefact: one self-describing,
//! bgzip-wrapped TSV listing every short-tandem-repeat locus in a reference,
//! each row carrying the tract coordinates, the repeat motif, a recomputed
//! purity, and the embedded local reference (`ref_seq` + `ref_seq_start`). It
//! is the *only* reference-bearing input the downstream stages need
//! (architecture [`ssr_catalog.md`](../../../doc/devel/architecture/ssr_catalog.md)).
//!
//! ```text
//! reference FASTA ─► [TRF-mod detect] ─► [post-process] ─► [embed ref_seq] ─► catalog.ssr_catalog.bed.gz (+ index)
//! ```
//!
//! **Build status (incremental).** The format I/O layer ([`io`]) lands first —
//! it is the cross-stage contract Stage 1's `fetch_reads` reader consumes, and
//! it is buildable and testable without the external `trf-mod` binary. The
//! detection/post-processing front-end follows:
//!
//! - [`io`] — **built**: [`io::CatalogHeader`], [`io::CatalogWriter`],
//!   [`io::CatalogReader`], and the `Locus` ⇄ row serialisation + round-trip.
//! - [`postprocess`] — **built**: [`postprocess::build_loci`] — the period≤6 →
//!   drop-compound → drop-bundle → end-trim → recompute-purity → embed-`ref_seq`
//!   pipeline (faithful GangSTR `minimal_trim`/`remove_bundles` port).
//! - [`trf`] — **partial**: the parsed [`trf::TrfRecord`] type (what
//!   `postprocess` consumes). The locate / spawn / BED-parse functions are the
//!   next increment (the `trf-mod` binary is now in the dev container).
//! - `run()` orchestrator + the `ssr-catalog` CLI subcommand — *pending*.

pub mod io;
pub mod postprocess;
pub mod trf;

/// Errors building or reading an SSR catalog (Stage 0).
///
/// `#[non_exhaustive]`: the detection / post-processing front-end will add
/// `trf-mod`-spawn and FASTA variants in later increments; only the format-I/O
/// variants exist today.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub(crate) enum CatalogError {
    /// An underlying read/write/bgzf failure. `context` names the operation.
    #[error("catalog I/O failed ({context})")]
    Io {
        context: &'static str,
        #[source]
        source: std::io::Error,
    },

    /// The `##` metadata header is missing a required key or carries an
    /// unparsable value.
    #[error("malformed catalog header: {reason}")]
    HeaderParse { reason: String },

    /// A data row does not have the expected tab-separated column shape, or a
    /// numeric field failed to parse.
    #[error("malformed catalog row at line {line}: {reason}")]
    RowParse { line: usize, reason: String },

    /// A row's `motif` field is not a valid SSR period.
    #[error("invalid motif at line {line}")]
    InvalidMotif {
        line: usize,
        #[source]
        source: crate::ssr::types::MotifError,
    },

    /// A row's coordinates / purity violate the [`Locus`](crate::ssr::types::Locus)
    /// invariants.
    #[error("invalid locus at line {line}")]
    InvalidLocus {
        line: usize,
        #[source]
        source: crate::ssr::types::LocusError,
    },
}
