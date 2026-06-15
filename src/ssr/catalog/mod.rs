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
//! - [`trf`] — **built**: [`trf::locate_trf_mod`], [`trf::version`],
//!   [`trf::run_on_contig`] (temp-file spawn, no pipes), and
//!   [`trf::parse_bed_line`] (the 10-column BED) → [`trf::TrfRecord`].
//! - `run()` orchestrator + the `ssr-catalog` CLI subcommand — *pending* (the
//!   last Stage-0 piece: per-contig fan-out/collect + header build + CSI index).

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

    /// `trf-mod` could not be located: no usable override, not found beside our
    /// own executable, and not on `PATH` (architecture §2.4 layered discovery).
    #[error("trf-mod binary not found (checked: override, exe-sibling, PATH)")]
    TrfModNotFound,

    /// Spawning or waiting on `trf-mod` failed at the OS level.
    #[error("failed to run trf-mod ({context})")]
    TrfSpawn {
        context: &'static str,
        #[source]
        source: std::io::Error,
    },

    /// `trf-mod` exited unsuccessfully while processing a contig.
    #[error("trf-mod failed on contig {contig:?} ({status})")]
    TrfRun { contig: String, status: String },

    /// `trf-mod -v` produced no recognisable version line.
    #[error("could not parse a version from `trf-mod -v`")]
    TrfVersion,

    /// A `trf-mod` BED line did not match the expected 10-column layout, a field
    /// failed to parse, or the contig-name column disagreed.
    #[error("malformed trf-mod BED at line {line}: {reason}")]
    TrfParse { line: usize, reason: String },
}
