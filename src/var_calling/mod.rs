//! Variant-calling stages (Stages 3–6).
//!
//! Sits between Stage 2 (the per-sample `.psp` reader) and the
//! posterior engine: DUST filter, grouping, per-group merger,
//! posterior. Used in cohort mode across many samples, but the same
//! code path applies to a single sample. See
//! `doc/devel/specs/calling_pipeline_architecture.md` for the full
//! stage breakdown.
//!
//! First occupant: [`per_position_merger`], the linear-scan k-way
//! merge over per-sample `.psp` record streams. Stage 4 lands in
//! [`variant_grouping`], the streaming overlap bundler that turns
//! the merger's output into `OverlappingVariantGroup`s for Stage 5.
//! Later stages land in sibling modules as they are implemented.

pub mod contamination_estimation;
pub mod dust_filter;
pub mod per_group_merger;
pub mod per_position_merger;
pub mod posterior_engine;
pub mod variant_grouping;

// Cohort-level chunk-at-a-time materialisation (formerly the `from_psp`
// submodule, merged up after the direct from-BAM path was removed). At each
// step of the chunk loop the driver loads one genomic chunk × N samples into
// memory (variant positions only), partitions it into windows whose
// boundaries fall in `max_group_span`-wide safe gaps, and runs the worker
// pipeline on each window to produce final variant records.
//
// - [`columns`] — [`SampleColumns`] (per-sample CSR columnar storage) and
//   [`MaterialisedChunk`] (one chunk × N samples).
// - [`loader`] — the batch [`load_chunk_from_iters`] (equivalence-test oracle)
//   and the production streaming [`StreamingBlockLoader`].
// M8 (review finding): the chunk-driver internals are implementation detail
// of this crate's cohort `var-calling`, not a public library surface —
// nothing outside the crate consumes them (the only external reach is the
// root `DEFAULT_*` consts below and `driver::record_fails_mapq_diff_t_for_test`,
// used by an integration test). They are kept `pub` for now: tightening to
// `pub(crate)` un-masks a layer of dead / test-only items that `pub` hid
// from dead-code analysis (unused re-exports + ~15 never-used methods/types
// — e.g. `WorkerUnifyError`, `clone_from_columns`, `materialise_record`,
// `into_reader`, `project_per_position_alleles_columnar`). Internalising the
// surface therefore needs a dead-code removal pass first; tracked as a
// follow-up rather than bundled into the review-fix run.
pub mod column_span_reader;
pub mod columns;
pub mod driver;
pub mod kernels;
pub mod loader;
pub mod partition;
#[cfg(test)]
pub(crate) mod test_helpers;
pub mod two_pass;
pub mod worker;

pub use column_span_reader::ColumnSpanReader;
pub use columns::{MaterialisedChunk, SampleColumns};
pub use driver::{
    ChunkDriverError, ChunkDriverParams, ChunkDriverStats, ChunkSizingParams,
    DEFAULT_CHUNK_GENOMIC_SPAN, DownstreamFilterParams, drive_cohort_chunked,
};
pub use loader::{
    ChunkLoadError, ChunkLoadExtent, ChunkLoadScratch, ChunkLoadStats, SpanColumnSource,
    StreamingBlockLoader, load_chunk_from_iters,
};
pub use partition::{PartitionError, PartitionScratch, WindowPartition, partition_window};
pub use worker::{WindowRunStats, WorkerScratch, run_window};

// ---------------------------------------------------------------------
// Cohort-pipeline downstream-filter defaults
// ---------------------------------------------------------------------
//
// These three CLI-default constants describe the per-record downstream
// filters applied post-EM (or, for `min_alt_obs`, at the merger→EM
// boundary). They live here, at the `var_calling` root, because they are
// consumed both by the CLI args layer (`cli::shared_args`) and by the
// chunk driver's `DownstreamFilterParams`.

/// Default minimum site-level `QUAL` (phred) required to emit a record.
/// Matches GATK HaplotypeCaller's
/// `--standard-min-confidence-threshold-for-calling` (30.0) as a hygiene
/// gate.
///
/// **Caveat.** Our site-level QUAL is currently
/// `-10·log10(Π_s P(GT_s = hom-ref))` — a *product* over samples, not a
/// marginalised cohort-AF posterior the way GATK (`P(AC = 0 | data)`,
/// integrated over a shared AF) and FreeBayes (`1 - P(K = 1 | reads)`)
/// compute it. Under our formula QUAL grows roughly linearly with cohort
/// size even when nothing about the per-sample evidence changes, so a
/// fixed threshold has different bite at different N. A fix to the QUAL
/// semantics is tracked separately; until then `30.0` is a coarse low-end
/// hygiene cutoff (drops the obvious tail) rather than a calibrated
/// quality bar.
pub const DEFAULT_MIN_QUAL_PHRED: f64 = 30.0;

/// Default minimum per-allele alt-read support. Records where no ALT
/// allele has `max(num_obs across samples) >= 2` are dropped *before*
/// reaching the EM. Empirically (multichrom 10-duplicate cohort,
/// 2026-05-20) `2` cuts FPs vs GATK by 70 % while losing < 1 pp recall;
/// 89.9 % of pileup-level candidates are single-alt-read positions that
/// come from base / alignment error rather than real variation.
///
/// `0` and `1` are both effective no-ops (every candidate has at least
/// one supporting alt read by construction); they're allowed as escape
/// hatches for tests / debugging.
///
/// Closest GATK analogue: `--min-pruning = 2` on the De Bruijn assembly,
/// which prunes assembly paths supported by fewer than two reads in the
/// cohort. Our filter is per-allele max across samples rather than
/// per-path, but the effect on candidate-set size is comparable.
pub const DEFAULT_MIN_ALT_OBS_PER_SAMPLE: u32 = 2;

/// Default Welch's-t threshold for the multi-mapper MAPQ-difference
/// filter. Empirically validated against GATK's cohort-level `MQRankSum`
/// annotation on the 10-duplicate synthetic cohort: catches 25 % of
/// GATK-extreme sites at a 2.0 % false-flag rate on GATK-clean sites. See
/// the per_allele_mapq_tracking implementation plan for the
/// cross-validation table.
pub const DEFAULT_MIN_MAPQ_DIFF_T: f32 = -3.0;
