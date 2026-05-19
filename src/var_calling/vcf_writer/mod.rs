//! Cohort VCF writer (Stage 6 sink).
//!
//! Streams [`PosteriorRecord`] items from the posterior engine and
//! emits a multi-sample VCF (plain text or bgzipped). The writer is
//! built around three boundaries:
//!
//! * [`CohortMetadata`] carries everything pinned at construction
//!   time — sample names (in the order matching
//!   `PosteriorRecord.posteriors` rows), the contig table sourced
//!   from the upstream merger's chromosome slice, and the provenance
//!   strings that go into `##source` / `##commandline` header lines.
//! * [`WriterConfig`] carries per-run knobs — the output path
//!   (suffix selects plain vs bgzf, matched case-insensitively), and
//!   the off-by-default [`WriterConfig::emit_gp`] flag. Construct via
//!   [`WriterConfig::new`]; there is no `impl Default` because
//!   `output` has no sensible default.
//! * [`CohortVcfWriter`] is the runtime state — opens
//!   `<output>.tmp`, writes the header at construction, accepts
//!   per-record `write_record` calls, and atomically renames into
//!   place on [`CohortVcfWriter::finish`].
//!
//! Errors flow through [`VcfWriteError`]; the enum is
//! `#[non_exhaustive]`, so callers writing `match` against it must
//! include a wildcard arm. A forgotten `finish` leaves `<output>.tmp`
//! on disk and no `<output>` — this is the intended loud-failure
//! mode: callers do not see a half-written VCF.
//!
//! Plan: `doc/devel/implementation_plans/cohort_vcf_writer.md`.
//!
//! [`PosteriorRecord`]: crate::var_calling::posterior_engine::PosteriorRecord

use std::path::PathBuf;

mod errors;
mod header;
mod record_encode;
mod sink;
mod writer;

pub use errors::VcfWriteError;
pub use header::CohortMetadata;
pub use writer::CohortVcfWriter;

/// Default for [`WriterConfig::emit_gp`]. Off — `GP` is `Number=G`
/// and most consumers don't read it (size grows as
/// `(ploidy + n_alleles - 1) choose ploidy`; 21 floats per sample
/// at ploidy=2, n_alleles=6). Opt in when the caller specifically
/// wants posteriors on disk.
pub const DEFAULT_EMIT_GP: bool = false;

/// Knobs that apply to one writer run. Narrow today; structured so
/// future flags slot in without changing the public constructor
/// signature.
///
/// Constructed via [`WriterConfig::new`]. There is no `Default` impl
/// because `output` has no sensible default — an empty `PathBuf`
/// would silently produce `.tmp` files in the process cwd and fail
/// at rename time with a confusing `io::Error`. Forcing the caller
/// to pass `output` explicitly makes the missing-field bug a compile
/// error.
///
/// **v1 FILTER policy:** every emitted record carries `PASS`. The
/// `default_filter_pass` knob that used to live on this struct has
/// been dropped — a future filter slice will re-introduce filter
/// expressions through a different surface (typed filter rules,
/// not a binary flag).
#[derive(Debug, Clone)]
pub struct WriterConfig {
    /// Output path. Suffix selects the sink kind (matched
    /// case-insensitively):
    /// * `.vcf.gz` or `.vcf.bgz` → bgzf
    /// * anything else           → plain text
    pub output: PathBuf,

    /// When true, declare `##FORMAT=<ID=GP,Number=G,Type=Float,…>`
    /// in the header and emit a `GP` cell per sample. Off by default
    /// ([`DEFAULT_EMIT_GP`]) — see the constant's doc for the
    /// rationale.
    pub emit_gp: bool,
}

impl WriterConfig {
    /// Build a config with [`DEFAULT_EMIT_GP`].
    pub fn new(output: PathBuf) -> Self {
        Self {
            output,
            emit_gp: DEFAULT_EMIT_GP,
        }
    }
}
