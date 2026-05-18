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
//!   (suffix selects plain vs bgzf), and the off-by-default
//!   [`WriterConfig::emit_gp`] flag.
//! * [`CohortVcfWriter`] is the runtime state — opens
//!   `<output>.tmp`, writes the header at construction, accepts
//!   per-record `write_record` calls, and atomically renames into
//!   place on [`CohortVcfWriter::finish`].
//!
//! Errors flow through [`VcfWriteError`]. A forgotten `finish` leaves
//! `<output>.tmp` on disk and no `<output>` — this is the intended
//! loud-failure mode: callers do not see a half-written VCF.
//!
//! Plan:
//! [doc/devel/implementation_plans/cohort_vcf_writer.md](https://example.invalid).

use std::path::PathBuf;

mod errors;
mod header;
mod record_encode;
mod sink;
mod writer;

pub use errors::VcfWriteError;
pub use header::CohortMetadata;
pub use writer::CohortVcfWriter;

/// Knobs that apply to one writer run. Narrow today; structured so
/// future flags (`--no-info-dp`, `--strict-monomorphic-as-ref`, etc.)
/// slot in without changing the public constructor signature.
#[derive(Debug, Clone)]
pub struct WriterConfig {
    /// Output path. Suffix selects the sink kind:
    /// * `.vcf.gz` or `.vcf.bgz` → bgzf
    /// * anything else           → plain text
    pub output: PathBuf,

    /// Always emit `PASS` in the FILTER column. Reserved for a future
    /// filter slice; v1 has no filter expressions.
    pub default_filter_pass: bool,

    /// When true, declare `##FORMAT=<ID=GP,Number=G,Type=Float,…>` in
    /// the header and emit a `GP` cell per sample. Off by default —
    /// `GP` is `Number=G` (size grows as `(ploidy + n_alleles - 1)
    /// choose ploidy`; 21 floats per sample at ploidy=2, n_alleles=6)
    /// and most downstream consumers don't read it. Opt in when the
    /// caller specifically wants posteriors on disk.
    pub emit_gp: bool,
}

impl Default for WriterConfig {
    fn default() -> Self {
        Self {
            output: PathBuf::new(),
            default_filter_pass: true,
            emit_gp: false,
        }
    }
}
