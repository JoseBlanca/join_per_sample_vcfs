//! Map a [`PosteriorRecord`] onto a `noodles_vcf::variant::RecordBuf`
//! ready for the writer to serialise.
//!
//! This module owns the field-by-field translation; everything
//! upstream is project types and everything downstream is noodles
//! types. Both directions are pure functions — no I/O, no state.
//!
//! The mapping is documented in
//! `doc/devel/implementation_plans/cohort_vcf_writer.md` §"Record
//! encoding (`record_encode.rs`)". Highlights:
//!
//! * `QUAL` is clamped to `[0, 9999]` (see [`QUAL_MAX`]).
//! * `GT` is unphased — Stage 5 does not currently emit phase — and
//!   the multiset of allele indices is sorted ascending (`0/1` not
//!   `1/0`) per the VCF spec.
//! * `GQ` is clamped to `[0, GQ_MAX]` ([`GQ_MAX`]).
//! * `GP` is only assembled when `config.emit_gp` is true.
//! * Per-record dependent vectors are validated up-front against
//!   the record's declared shape; a malformed `PosteriorRecord`
//!   surfaces as `VcfWriteError::InconsistentRecord` instead of
//!   panicking on slice indexing.

use std::fmt::Write as _;

use noodles_core::Position;
use noodles_vcf::variant::RecordBuf;
use noodles_vcf::variant::record_buf::{
    AlternateBases, Filters, Ids, Info, Samples, info::field::Value as InfoValue,
    info::field::value::Array as InfoArray, samples::Keys, samples::sample::Value as SampleValue,
};

use super::WriterConfig;
use super::errors::VcfWriteError;
use super::writable::VcfWritable;
use crate::psp::header::ParsedChromosome;

/// Hard cap on the QUAL column.
///
/// `PosteriorRecord.qual_phred` is `f64::INFINITY` when every sample
/// is certainly variant. VCF has no sentinel for "infinite quality";
/// we cap at the same value GATK and freebayes use so downstream
/// tooling parses it as a number rather than choking. Search the
/// codebase for `9999` to find this comment.
pub(super) const QUAL_MAX: f32 = 9999.0;

/// Hard cap on per-sample GQ.
///
/// GATK and bcftools convention is to cap GQ at 99 (single-byte
/// Phred). The writer cap must be at least as high as the engine's
/// `PosteriorEngineConfig::max_gq_phred`; the `gq_writer_cap_is_at_least_engine_cap`
/// test latches that drift.
pub(super) const GQ_MAX: f32 = 99.0;

/// Build the four-or-five-key FORMAT-keys list, once per writer
/// instance. Stored on `CohortVcfWriter` and cloned into each
/// `Samples` so we don't rebuild it per record.
pub(super) fn build_format_keys(config: &WriterConfig) -> Keys {
    let mut keys: Vec<String> = vec!["GT".into(), "GQ".into(), "DP".into(), "AD".into()];
    if config.emit_gp {
        keys.push("GP".into());
    }
    keys.into_iter().collect()
}

/// Encode a single posterior record into a noodles `RecordBuf`.
///
/// `contigs` is the same slice the writer holds — used to map
/// `record.locus.chrom_id` back to the contig name. `format_keys`
/// is pre-built via [`build_format_keys`] and shared across records.
/// `table` is the canonical `genotype_order(ploidy, n_alleles)`
/// table for this record's allele cardinality; the writer caches it
/// per `(ploidy, n_alleles)` so we don't rebuild per record.
pub(super) fn encode<R: VcfWritable>(
    record: &R,
    contigs: &[ParsedChromosome],
    config: &WriterConfig,
    format_keys: &Keys,
    expected_samples: usize,
    table: &[Vec<u8>],
) -> Result<RecordBuf, VcfWriteError> {
    validate_record_shape(record, expected_samples, table)?;

    let chrom_name = contigs
        .get(record.chrom_id() as usize)
        .map(|c| c.name.clone())
        .ok_or(VcfWriteError::UnknownChromId {
            chrom_id: record.chrom_id(),
            pos: record.pos_1based(),
            n_contigs: contigs.len(),
        })?;

    let position = build_position(record)?;

    let ref_bases = allele_to_string(record.allele_seq(0), record.chrom_id(), record.pos_1based())?;
    let alt_strings: Vec<String> = (1..record.n_alleles())
        .map(|i| allele_to_string(record.allele_seq(i), record.chrom_id(), record.pos_1based()))
        .collect::<Result<_, _>>()?;
    let alternate_bases = AlternateBases::from(alt_strings);

    // Refine QUAL: deflate it at sites whose alt-allele support looks
    // like a systematic artifact (allele-balance + strand/position bias).
    // Affects the QUAL column only; genotypes/GQ/AF are unchanged.
    let qual = clamp_qual(super::qual_refine::refine_qual(
        record,
        table,
        record.qual_phred(),
    ));

    // FILTER policy:
    //   * `PASS`     — EM converged within `max_iterations`.
    //   * `EMNoConv` — EM hit the iteration cap; the record is still
    //                  emitted with allele frequencies + posteriors
    //                  from a final E-step under the un-converged
    //                  p̂/f̂. Downstream consumers can prune via
    //                  `bcftools view -f PASS`.
    // A future filter slice will re-introduce filter *expressions*
    // through a typed-rule surface; the per-record bit here is the
    // engine's diagnostics signal, not the rule surface.
    let filters = if record.converged() {
        Filters::pass()
    } else {
        [String::from(super::header::EM_NO_CONV_FILTER_ID)]
            .into_iter()
            .collect()
    };

    let n_alleles = record.n_alleles();

    let info = build_info(record, n_alleles, table)?;
    let samples = build_samples(record, format_keys, config, table)?;

    Ok(RecordBuf::builder()
        .set_reference_sequence_name(chrom_name)
        .set_variant_start(position)
        .set_ids(Ids::default())
        .set_reference_bases(ref_bases)
        .set_alternate_bases(alternate_bases)
        .set_quality_score(qual)
        .set_filters(filters)
        .set_info(info)
        .set_samples(samples)
        .build())
}

/// M6: validate every per-record vector length against the record's
/// declared shape before any indexing. A malformed upstream becomes
/// a typed `InconsistentRecord` or `SampleCountMismatch` instead of
/// a panic at the Stage 6 sink.
fn validate_record_shape<R: VcfWritable>(
    record: &R,
    expected_samples: usize,
    table: &[Vec<u8>],
) -> Result<(), VcfWriteError> {
    if record.n_samples() != expected_samples {
        return Err(VcfWriteError::SampleCountMismatch {
            chrom_id: record.chrom_id(),
            pos: record.pos_1based(),
            expected_samples,
            got_samples: record.n_samples(),
        });
    }
    let n_alleles = record.n_alleles();
    if n_alleles == 0 {
        return Err(VcfWriteError::InconsistentRecord {
            chrom_id: record.chrom_id(),
            pos: record.pos_1based(),
            field: "alleles",
            expected: 1,
            actual: 0,
        });
    }
    // The table passed in must match the record's `n_genotypes`. The
    // writer is responsible for building it from `(ploidy, n_alleles)`;
    // a mismatch is an internal bug rather than upstream malformation,
    // but we surface it as `InconsistentRecord` so the error path is
    // uniform.
    if table.len() != record.n_genotypes() {
        return Err(VcfWriteError::InconsistentRecord {
            chrom_id: record.chrom_id(),
            pos: record.pos_1based(),
            field: "n_genotypes/genotype_order table",
            expected: record.n_genotypes(),
            actual: table.len(),
        });
    }
    let n_samples = record.n_samples();
    for &(field, expected, actual) in &[
        ("best_genotype", n_samples, record.best_genotype().len()),
        ("gq_phred", n_samples, record.gq_phred().len()),
        (
            "allele_frequencies",
            n_alleles,
            record.allele_frequencies().len(),
        ),
        (
            "compound_frequencies",
            n_alleles,
            record.compound_frequencies().len(),
        ),
        ("scalars", n_samples * n_alleles, record.scalars_len()),
        (
            "posteriors",
            n_samples * record.n_genotypes(),
            record.posteriors_len(),
        ),
        (
            "chain_anchor_flags",
            n_samples * n_alleles,
            record.chain_anchor_flags_len(),
        ),
    ] {
        if expected != actual {
            return Err(VcfWriteError::InconsistentRecord {
                chrom_id: record.chrom_id(),
                pos: record.pos_1based(),
                field,
                expected,
                actual,
            });
        }
    }
    Ok(())
}

/// Convert the record's 1-based `start` to a noodles `Position`.
/// Carries operation context through `VcfWriteError::Encode`.
fn build_position<R: VcfWritable>(record: &R) -> Result<Position, VcfWriteError> {
    let pos_usize = usize::try_from(record.pos_1based()).map_err(|e| VcfWriteError::Encode {
        operation: "1-based position to usize",
        source: Box::new(e),
    })?;
    Position::try_from(pos_usize).map_err(|e| VcfWriteError::Encode {
        operation: "usize to noodles Position",
        source: Box::new(e),
    })
}

/// Convert a `Vec<u8>` allele sequence (guaranteed by Stage 1 to be
/// uppercase `{A,C,G,T,N}`) into an owned `String` for the VCF
/// encoder. An unexpected byte surfaces as
/// `VcfWriteError::Encode { operation: "allele bytes UTF-8", ... }`.
fn allele_to_string(seq: &[u8], _chrom_id: u32, _pos: u32) -> Result<String, VcfWriteError> {
    // Convert to a String, but on UTF-8 failure preserve the original
    // bytes' provenance by wrapping the typed `Utf8Error` rather than
    // formatting it away.
    std::str::from_utf8(seq)
        .map(|s| s.to_string())
        .map_err(|e| VcfWriteError::Encode {
            operation: "allele bytes UTF-8",
            source: Box::new(e),
        })
}

/// Clamp QUAL to `[0, QUAL_MAX]`. `INFINITY` maps to `QUAL_MAX` (every
/// sample certainly variant); `NaN` maps to 0 defensively — the EM
/// math should never produce NaN for an emitted record, so this is a
/// guard not an expected path.
fn clamp_qual(qual_phred: f64) -> f32 {
    if qual_phred.is_nan() {
        return 0.0;
    }
    let q = qual_phred as f32;
    q.clamp(0.0, QUAL_MAX)
}

fn build_info<R: VcfWritable>(
    record: &R,
    n_alleles: usize,
    table: &[Vec<u8>],
) -> Result<Info, VcfWriteError> {
    let mut info = Info::default();
    let n_alts = n_alleles.saturating_sub(1);

    // AF — allele frequencies for ALT alleles only (REF is implied).
    // `Number=A`: empty list when no ALTs.
    if n_alts > 0 {
        let af_values: Vec<Option<f32>> = record.allele_frequencies()[1..]
            .iter()
            .map(|f| Some(*f as f32))
            .collect();
        info.as_mut().insert(
            "AF".into(),
            Some(InfoValue::Array(InfoArray::Float(af_values))),
        );
    }

    // AC — per-ALT count of allele copies in argmax genotypes summed
    // across samples; AN — total called allele number.
    let (ac_per_alt, an_total) = tally_called_alleles(record, n_alleles, table)?;
    if !ac_per_alt.is_empty() {
        let ac_values: Vec<Option<i32>> = ac_per_alt.iter().map(|&c| Some(c as i32)).collect();
        info.as_mut().insert(
            "AC".into(),
            Some(InfoValue::Array(InfoArray::Integer(ac_values))),
        );
    }
    info.as_mut()
        .insert("AN".into(), Some(InfoValue::Integer(an_total as i32)));

    // DP — total depth across samples = sum of all per-allele
    // observation counts.
    let dp_total: u64 = (0..record.n_samples())
        .map(|s| {
            record
                .scalars_row(s)
                .iter()
                .map(|stats| u64::from(stats.num_obs))
                .sum::<u64>()
        })
        .sum();
    let dp_total_i32 = i32::try_from(dp_total).map_err(|_| VcfWriteError::DepthOverflow {
        chrom_id: record.chrom_id(),
        pos: record.pos_1based(),
        sample_idx: None,
        depth: dp_total,
    })?;
    info.as_mut()
        .insert("DP".into(), Some(InfoValue::Integer(dp_total_i32)));

    // CA — flag set when any sample carries a chain-anchor-broken
    // call. `chain_anchor_flags` is a flat row-major `Vec<bool>` of
    // length `n_samples * n_alleles`; iterate row-by-row through the
    // trait's `chain_anchor_flags_row(sample_idx)` accessor.
    if (0..record.n_samples()).any(|s| record.chain_anchor_flags_row(s).iter().any(|&b| b)) {
        info.as_mut().insert("CA".into(), Some(InfoValue::Flag));
    }

    // MQRef / MQAlt / MQDiff / MQDiffT — cohort-pooled MAPQ stats.
    // Number=1 for MQRef (REF is single); Number=A for MQAlt /
    // MQDiff / MQDiffT (one entry per ALT). All four emit `.` when
    // the underlying counts are insufficient (variance undefined or
    // no reads of the relevant class).
    emit_mapq_info(record, n_alleles, &mut info);

    Ok(info)
}

/// Per-allele MAPQ summary across the cohort. For each allele, sum
/// `num_obs`, `mapq_sum`, `mapq_sum_sq` across samples; from the
/// totals derive mean and Welch's-t-statistic for ALT vs REF.
///
/// Storage shape: returned vector indexed by allele (0 = REF).
struct CohortMapqStats {
    n: u32,
    mean: f64,
    variance: f64,
}
impl CohortMapqStats {
    fn from_pool(n: u32, sum: u64, sum_sq: u128) -> Option<Self> {
        if n == 0 {
            return None;
        }
        let n_f = n as f64;
        let mean = (sum as f64) / n_f;
        // Sample variance via the identity
        //   var = (sum_sq − sum²/n) / (n − 1)
        // Requires n ≥ 2; for n == 1 we report variance = 0 (single
        // observation — Welch's t can't be computed regardless).
        let variance = if n >= 2 {
            let sum_sq_f = sum_sq as f64;
            ((sum_sq_f - (sum as f64) * mean) / ((n - 1) as f64)).max(0.0)
        } else {
            0.0
        };
        Some(Self { n, mean, variance })
    }
}

fn pool_allele_mapq<R: VcfWritable>(record: &R, allele_idx: usize) -> Option<CohortMapqStats> {
    let mut n: u64 = 0;
    let mut sum: u64 = 0;
    let mut sum_sq: u128 = 0;
    for sample_idx in 0..record.n_samples() {
        let row = record.scalars_row(sample_idx);
        let stats = &row[allele_idx];
        n += u64::from(stats.num_obs);
        sum += u64::from(stats.mapq_sum);
        sum_sq += stats.mapq_sum_sq as u128;
    }
    // Cohort `num_obs` overflows u32 only at depths well beyond the
    // pipeline's per-sample max — but the per-allele stat downstream
    // only ever fits the per-sample-summed u64 sum into mean/var
    // computations, so we keep n as u32 by saturating.
    let n_u32 = u32::try_from(n).unwrap_or(u32::MAX);
    CohortMapqStats::from_pool(n_u32, sum, sum_sq)
}

fn welch_t(alt: &CohortMapqStats, ref_: &CohortMapqStats) -> Option<f32> {
    // Need ≥2 observations on each side to estimate variance.
    if alt.n < 2 || ref_.n < 2 {
        return None;
    }
    let se2 = alt.variance / (alt.n as f64) + ref_.variance / (ref_.n as f64);
    if se2 <= 0.0 {
        // Degenerate variance: if means are equal report t = 0;
        // otherwise the statistic is undefined (a perfectly uniform
        // group means the t-test can't speak; let the filter pass
        // these — they aren't the multi-mapper pattern we target).
        if (alt.mean - ref_.mean).abs() < f64::EPSILON {
            return Some(0.0);
        }
        return None;
    }
    Some(((alt.mean - ref_.mean) / se2.sqrt()) as f32)
}

fn emit_mapq_info<R: VcfWritable>(record: &R, n_alleles: usize, info: &mut Info) {
    if n_alleles < 2 {
        // REF-only sites never reach the VCF writer (cohort_driver
        // drops them upstream), but guard defensively.
        return;
    }
    let ref_stats = pool_allele_mapq(record, 0);
    // MQRef — Number=1 Type=Float. Emit only when defined; omit
    // entirely if no REF reads anywhere in the cohort.
    if let Some(r) = ref_stats.as_ref() {
        info.as_mut()
            .insert("MQRef".into(), Some(InfoValue::Float(r.mean as f32)));
    }
    // Per-ALT MQAlt, MQDiff, MQDiffT
    let n_alts = n_alleles - 1;
    let mut mq_alt_vec: Vec<Option<f32>> = Vec::with_capacity(n_alts);
    let mut mq_diff_vec: Vec<Option<f32>> = Vec::with_capacity(n_alts);
    let mut mq_diff_t_vec: Vec<Option<f32>> = Vec::with_capacity(n_alts);
    for alt_idx in 1..n_alleles {
        let alt_stats = pool_allele_mapq(record, alt_idx);
        let alt_mean = alt_stats.as_ref().map(|s| s.mean as f32);
        mq_alt_vec.push(alt_mean);
        let diff = match (alt_stats.as_ref(), ref_stats.as_ref()) {
            (Some(a), Some(r)) => Some((a.mean - r.mean) as f32),
            _ => None,
        };
        mq_diff_vec.push(diff);
        let t = match (alt_stats.as_ref(), ref_stats.as_ref()) {
            (Some(a), Some(r)) => welch_t(a, r),
            _ => None,
        };
        mq_diff_t_vec.push(t);
    }
    info.as_mut().insert(
        "MQAlt".into(),
        Some(InfoValue::Array(InfoArray::Float(mq_alt_vec))),
    );
    info.as_mut().insert(
        "MQDiff".into(),
        Some(InfoValue::Array(InfoArray::Float(mq_diff_vec))),
    );
    info.as_mut().insert(
        "MQDiffT".into(),
        Some(InfoValue::Array(InfoArray::Float(mq_diff_t_vec))),
    );
}

/// Walk each sample's argmax genotype, decode it through the
/// `genotype_order` table, and tally per-ALT counts plus the total
/// called allele number.
///
/// B1: an allele index past the record's `alleles.len() - 1`
/// surfaces as `VcfWriteError::AlleleIndexOutOfBounds` — *not* a
/// silent drop. The previous shape (`if alt_pos < n_alts { tally
/// += 1 }`) violated the VCF spec invariant `sum(AC) + REF == AN`
/// when the upstream produced an inconsistent genotype.
fn tally_called_alleles<R: VcfWritable>(
    record: &R,
    n_alleles: usize,
    table: &[Vec<u8>],
) -> Result<(Vec<u32>, u32), VcfWriteError> {
    let n_alts = n_alleles.saturating_sub(1);
    let mut ac_per_alt = vec![0u32; n_alts];
    let mut an_total: u32 = 0;
    for (sample_idx, &gt_idx) in record.best_genotype().iter().enumerate() {
        let gt = lookup_genotype(table, record, sample_idx, gt_idx)?;
        for &allele_idx in gt {
            if (allele_idx as usize) >= n_alleles {
                return Err(VcfWriteError::AlleleIndexOutOfBounds {
                    chrom_id: record.chrom_id(),
                    pos: record.pos_1based(),
                    sample_idx,
                    allele_idx,
                    n_alleles,
                });
            }
            an_total += 1;
            if allele_idx == 0 {
                continue; // REF
            }
            let alt_pos = (allele_idx as usize) - 1;
            ac_per_alt[alt_pos] += 1;
        }
    }
    Ok((ac_per_alt, an_total))
}

fn build_samples<R: VcfWritable>(
    record: &R,
    format_keys: &Keys,
    config: &WriterConfig,
    table: &[Vec<u8>],
) -> Result<Samples, VcfWriteError> {
    let n_alleles = record.n_alleles();
    let mut rows: Vec<Vec<Option<SampleValue>>> = Vec::with_capacity(record.n_samples());
    for sample_idx in 0..record.n_samples() {
        let gt_idx = record.best_genotype()[sample_idx];
        let gt = lookup_genotype(table, record, sample_idx, gt_idx)?;
        // Mi8: drop the gt_buf clear/reuse pattern — the trailing
        // .clone() inside the loop allocated a fresh String per row
        // anyway, so the buffer reuse saved zero allocations and
        // hid the per-iteration shape from the reader. A short
        // String::new() per iteration is clearer and the same cost.
        let mut gt_buf = String::with_capacity(2 * gt.len());
        format_gt_unphased(&mut gt_buf, gt);

        let gq = record.gq_phred()[sample_idx]
            .round()
            .clamp(0.0, GQ_MAX as f64) as i32;

        let scalars = record.scalars_row(sample_idx);
        let dp_sample: u64 = scalars.iter().map(|stats| u64::from(stats.num_obs)).sum();
        let dp_value = i32::try_from(dp_sample).map_err(|_| VcfWriteError::DepthOverflow {
            chrom_id: record.chrom_id(),
            pos: record.pos_1based(),
            sample_idx: Some(sample_idx),
            depth: dp_sample,
        })?;

        // M8: replace the silent `num_obs as i32` wrap with a typed
        // overflow error matched to the DP overflow shape.
        let ad_values: Vec<Option<i32>> = scalars
            .iter()
            .take(n_alleles)
            .map(|stats| {
                i32::try_from(stats.num_obs)
                    .map(Some)
                    .map_err(|_| VcfWriteError::DepthOverflow {
                        chrom_id: record.chrom_id(),
                        pos: record.pos_1based(),
                        sample_idx: Some(sample_idx),
                        depth: u64::from(stats.num_obs),
                    })
            })
            .collect::<Result<_, _>>()?;

        let mut row: Vec<Option<SampleValue>> = Vec::with_capacity(format_keys.as_ref().len());
        row.push(Some(SampleValue::String(gt_buf)));
        row.push(Some(SampleValue::Integer(gq)));
        row.push(Some(SampleValue::Integer(dp_value)));
        row.push(Some(SampleValue::from(ad_values)));

        if config.emit_gp {
            let gp_values: Vec<Option<f32>> = record
                .posteriors_row(sample_idx)
                .iter()
                .map(|p| Some(*p as f32))
                .collect();
            row.push(Some(SampleValue::from(gp_values)));
        }

        rows.push(row);
    }

    Ok(Samples::new(format_keys.clone(), rows))
}

/// Mi7: shared genotype-table lookup. Both `tally_called_alleles`
/// and `build_samples` walk `record.best_genotype` and need to
/// resolve `gt_idx` against the table; centralising the lookup keeps
/// the five-field `GenotypeIndexOutOfBounds` constructor in one
/// place.
fn lookup_genotype<'a, R: VcfWritable>(
    table: &'a [Vec<u8>],
    record: &R,
    sample_idx: usize,
    gt_idx: usize,
) -> Result<&'a [u8], VcfWriteError> {
    table
        .get(gt_idx)
        .map(Vec::as_slice)
        .ok_or(VcfWriteError::GenotypeIndexOutOfBounds {
            chrom_id: record.chrom_id(),
            pos: record.pos_1based(),
            sample_idx,
            got: gt_idx,
            n_genotypes: table.len(),
        })
}

/// Format a multiset of allele indices as an unphased `a/b[/c…]`
/// genotype string. The input is assumed already sorted ascending
/// (the `genotype_order` table builds them in canonical order, so a
/// re-sort would be redundant for now — kept as a function so a
/// future phased-emission flag is a localised change).
fn format_gt_unphased(out: &mut String, alleles: &[u8]) {
    for (i, &a) in alleles.iter().enumerate() {
        if i > 0 {
            out.push('/');
        }
        // `write!` into a String never fails (the only path that
        // could is `fmt::Error`, which `String`'s `Write` impl never
        // produces); name that explicitly via `.expect` so the
        // discarded-Result rule isn't violated.
        write!(out, "{a}").expect("write to String is infallible");
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;
    use crate::pileup_record::AlleleSupportStats;
    use crate::var_calling::per_group_merger::{MergedAllele, genotype_order};
    use crate::var_calling::posterior_engine::{EmDiagnostics, PosteriorRecord, RecordLocus};

    fn ref_allele(seq: &[u8]) -> MergedAllele {
        MergedAllele {
            seq: seq.to_vec(),
            is_compound: false,
            constituents: Vec::new(),
        }
    }

    fn alt_allele(seq: &[u8]) -> MergedAllele {
        ref_allele(seq)
    }

    fn support(num_obs: u32) -> AlleleSupportStats {
        AlleleSupportStats {
            num_obs,
            q_sum: 0.0,
            fwd: 0,
            placed_left: 0,
            placed_start: 0,

            mapq_sum: 0,
            mapq_sum_sq: 0,
        }
    }

    /// 2-sample biallelic SNP A→T; best_genotype = [0/0, 0/1].
    /// Genotype enumeration at ploidy=2, n_alleles=2 is
    /// [[0,0], [0,1], [1,1]], so indices 0 → 0/0 and 1 → 0/1.
    fn biallelic_two_samples() -> PosteriorRecord {
        let alleles = vec![ref_allele(b"A"), alt_allele(b"T")];
        let scalars = vec![
            support(20), // S0, REF
            support(0),  // S0, ALT
            support(10), // S1, REF
            support(10), // S1, ALT
        ];
        let posteriors = vec![
            0.98, 0.01, 0.01, // S0
            0.05, 0.90, 0.05, // S1
        ];
        let best_genotype = vec![0usize, 1usize];
        PosteriorRecord {
            locus: RecordLocus {
                chrom_id: 0,
                start: 1234,
                end: 1234,
            },
            alleles,
            ploidy: 2,
            n_samples: 2,
            n_genotypes: 3,
            allele_frequencies: vec![0.75, 0.25],
            compound_frequencies: vec![None, None],
            posteriors,
            best_genotype,
            gq_phred: vec![60.0, 40.0],
            qual_phred: 200.0,
            scalars,
            other_scalars: vec![],
            chain_anchor_flags: vec![false; 4],
            diagnostics: EmDiagnostics {
                iterations: 5,
                final_max_delta_p: 1e-6,
                converged: true,
            },
        }
    }

    fn fixture_contigs() -> Vec<ParsedChromosome> {
        vec![ParsedChromosome {
            name: "chr1".into(),
            length: 1_000_000,
            md5: "00000000000000000000000000000001".into(),
        }]
    }

    fn cfg_default() -> WriterConfig {
        WriterConfig::new(PathBuf::from("/dev/null"))
    }

    fn cfg_emit_gp_on() -> WriterConfig {
        WriterConfig::new(PathBuf::from("/dev/null")).with_emit_gp(true)
    }

    fn table_for(record: &PosteriorRecord) -> Vec<Vec<u8>> {
        genotype_order(record.ploidy, record.alleles.len())
    }

    /// Drive a `RecordBuf` through the writer to a String; lets tests
    /// inspect the actual rendered bytes.
    fn record_to_line(record: &RecordBuf, header: &noodles_vcf::Header) -> String {
        use noodles_vcf::variant::io::Write as _;
        let mut writer = noodles_vcf::io::Writer::new(Vec::new());
        writer.write_variant_record(header, record).unwrap();
        String::from_utf8(writer.into_inner()).unwrap()
    }

    fn fixture_header(config: &WriterConfig) -> noodles_vcf::Header {
        use super::super::header::{CohortMetadata, build_vcf_header};
        let metadata = CohortMetadata {
            sample_names: vec!["S0".into(), "S1".into()],
            contigs: fixture_contigs(),
            tool_string: "test 0".into(),
            command_line: String::new(),
        };
        build_vcf_header(&metadata, config).unwrap()
    }

    #[test]
    fn biallelic_snp_default_config_renders_expected_line() {
        let record = biallelic_two_samples();
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);
        let table = table_for(&record);

        let buf = encode(&record, &contigs, &config, &keys, 2, &table).unwrap();
        let line = record_to_line(&buf, &header);

        let line = line.trim_end();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 11, "unexpected field count: {line}");
        assert_eq!(fields[0], "chr1");
        assert_eq!(fields[1], "1234");
        assert_eq!(fields[2], ".");
        assert_eq!(fields[3], "A");
        assert_eq!(fields[4], "T");
        assert!(fields[5].starts_with("200"), "QUAL = {}", fields[5]);
        assert_eq!(fields[6], "PASS");
        let info = fields[7];
        assert!(info.contains("AF=0.25"), "INFO missing AF: {info}");
        assert!(info.contains("AC=1"), "INFO missing AC: {info}");
        assert!(info.contains("AN=4"), "INFO missing AN: {info}");
        assert!(info.contains("DP=40"), "INFO missing DP: {info}");
        assert!(!info.contains("CA"), "INFO carried unexpected CA: {info}");
        assert_eq!(fields[8], "GT:GQ:DP:AD");
        assert_eq!(fields[9], "0/0:60:20:20,0", "S0 cell: {}", fields[9]);
        assert_eq!(fields[10], "0/1:40:20:10,10", "S1 cell: {}", fields[10]);
    }

    #[test]
    fn unconverged_record_renders_with_emnoconv_filter() {
        let mut record = biallelic_two_samples();
        record.diagnostics.converged = false;
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);
        let table = table_for(&record);

        let buf = encode(&record, &contigs, &config, &keys, 2, &table).unwrap();
        let line = record_to_line(&buf, &header);
        let fields: Vec<&str> = line.trim_end().split('\t').collect();

        // Same record shape as the happy-path test — only FILTER changes.
        assert_eq!(
            fields[6], "EMNoConv",
            "non-converging records must carry the EMNoConv FILTER, got {} (full line: {line})",
            fields[6],
        );

        // And the header must declare the filter so a strict
        // consumer accepts the file. `write_header` is an inherent
        // method on `noodles_vcf::io::Writer` — no trait import.
        let header_text = {
            let mut writer = noodles_vcf::io::Writer::new(Vec::new());
            writer.write_header(&header).unwrap();
            String::from_utf8(writer.into_inner()).unwrap()
        };
        assert!(
            header_text.contains("##FILTER=<ID=EMNoConv"),
            "header missing ##FILTER=<ID=EMNoConv> declaration: {header_text}",
        );
    }

    #[test]
    fn converged_record_renders_with_pass_filter() {
        // Lower-stakes regression check — the converged path must
        // still emit PASS after the FILTER-policy change.
        let mut record = biallelic_two_samples();
        record.diagnostics.converged = true;
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);
        let table = table_for(&record);

        let buf = encode(&record, &contigs, &config, &keys, 2, &table).unwrap();
        let line = record_to_line(&buf, &header);
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(fields[6], "PASS");
    }

    #[test]
    fn biallelic_snp_with_emit_gp_adds_gp_column() {
        let record = biallelic_two_samples();
        let contigs = fixture_contigs();
        let config = cfg_emit_gp_on();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);
        let table = table_for(&record);

        let buf = encode(&record, &contigs, &config, &keys, 2, &table).unwrap();
        let line = record_to_line(&buf, &header).trim_end().to_string();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[8], "GT:GQ:DP:AD:GP");
        let s0 = fields[9];
        let s0_parts: Vec<&str> = s0.split(':').collect();
        assert_eq!(s0_parts.len(), 5, "S0 with GP: {s0}");
        let gp_floats: Vec<f32> = s0_parts[4].split(',').map(|s| s.parse().unwrap()).collect();
        assert_eq!(gp_floats.len(), 3);
        assert!((gp_floats.iter().sum::<f32>() - 1.0).abs() < 1e-3);
    }

    #[test]
    fn triallelic_snp_emits_two_alts_with_number_a_lists() {
        let mut record = biallelic_two_samples();
        record.alleles.push(alt_allele(b"C"));
        record.n_genotypes = 6;
        record.allele_frequencies = vec![0.6, 0.25, 0.15];
        record.compound_frequencies = vec![None, None, None];
        record.posteriors = vec![
            0.7, 0.1, 0.05, 0.1, 0.025, 0.025, // S0
            0.05, 0.05, 0.05, 0.7, 0.1, 0.05, // S1 → best = 3 (0/2)
        ];
        record.best_genotype = vec![0, 3];
        record.scalars = vec![
            support(20),
            support(0),
            support(0), // S0
            support(10),
            support(0),
            support(10), // S1
        ];
        record.chain_anchor_flags = vec![false; 6];

        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);
        let table = table_for(&record);

        let buf = encode(&record, &contigs, &config, &keys, 2, &table).unwrap();
        let line = record_to_line(&buf, &header).trim_end().to_string();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[3], "A");
        assert_eq!(fields[4], "T,C");
        let info = fields[7];
        assert!(info.contains("AF=0.25,0.15"), "INFO: {info}");
        assert!(info.contains("AC=0,1"), "INFO: {info}");
    }

    #[test]
    fn chain_anchor_flag_surfaces_as_ca_info_flag() {
        let mut record = biallelic_two_samples();
        record.chain_anchor_flags[3] = true;
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);
        let table = table_for(&record);

        let buf = encode(&record, &contigs, &config, &keys, 2, &table).unwrap();
        let line = record_to_line(&buf, &header).trim_end().to_string();
        let fields: Vec<&str> = line.split('\t').collect();
        let info = fields[7];
        assert!(info.contains("CA"), "INFO should carry CA: {info}");
    }

    #[test]
    fn infinite_qual_renders_as_9999_cap() {
        let mut record = biallelic_two_samples();
        record.qual_phred = f64::INFINITY;
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);
        let table = table_for(&record);

        let buf = encode(&record, &contigs, &config, &keys, 2, &table).unwrap();
        let line = record_to_line(&buf, &header).trim_end().to_string();
        let fields: Vec<&str> = line.split('\t').collect();
        let qual: f32 = fields[5].parse().unwrap();
        assert_eq!(qual, QUAL_MAX, "QUAL should be capped at 9999, got {qual}");
    }

    #[test]
    fn nan_qual_renders_as_zero() {
        let mut record = biallelic_two_samples();
        record.qual_phred = f64::NAN;
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);
        let table = table_for(&record);

        let buf = encode(&record, &contigs, &config, &keys, 2, &table).unwrap();
        let line = record_to_line(&buf, &header).trim_end().to_string();
        let fields: Vec<&str> = line.split('\t').collect();
        let qual: f32 = fields[5].parse().unwrap();
        assert_eq!(qual, 0.0);
    }

    #[test]
    fn unknown_chrom_id_is_an_error() {
        let mut record = biallelic_two_samples();
        record.locus.chrom_id = 99;
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let table = table_for(&record);
        let err = encode(&record, &contigs, &config, &keys, 2, &table).unwrap_err();
        assert!(matches!(
            err,
            VcfWriteError::UnknownChromId { chrom_id: 99, .. }
        ));
    }

    #[test]
    fn sample_count_mismatch_is_an_error() {
        let record = biallelic_two_samples();
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let table = table_for(&record);
        // Cohort metadata names 3 samples but the record has 2.
        let err = encode(&record, &contigs, &config, &keys, 3, &table).unwrap_err();
        assert!(matches!(
            err,
            VcfWriteError::SampleCountMismatch {
                expected_samples: 3,
                got_samples: 2,
                ..
            }
        ));
    }

    /// M5: SampleCountMismatch's rendered message names the cohort
    /// metadata count first ("expected"), the posterior arrays' count
    /// second ("got"). Locks the wording against future field renames.
    #[test]
    fn sample_count_mismatch_message_names_cohort_first() {
        let err = VcfWriteError::SampleCountMismatch {
            chrom_id: 0,
            pos: 1,
            expected_samples: 3,
            got_samples: 2,
        };
        let s = err.to_string();
        assert!(s.contains("cohort metadata names 3"), "{s}");
        assert!(s.contains("posterior arrays carry 2"), "{s}");
    }

    #[test]
    fn format_keys_track_emit_gp() {
        let off = build_format_keys(&cfg_default());
        let on = build_format_keys(&cfg_emit_gp_on());

        let off_vec: Vec<&str> = off.as_ref().iter().map(String::as_str).collect();
        let on_vec: Vec<&str> = on.as_ref().iter().map(String::as_str).collect();
        assert_eq!(off_vec, vec!["GT", "GQ", "DP", "AD"]);
        assert_eq!(on_vec, vec!["GT", "GQ", "DP", "AD", "GP"]);
    }

    #[test]
    fn format_gt_unphased_joins_with_slashes() {
        let mut s = String::new();
        format_gt_unphased(&mut s, &[0, 1]);
        assert_eq!(s, "0/1");
        s.clear();
        format_gt_unphased(&mut s, &[1, 1, 2]);
        assert_eq!(s, "1/1/2");
    }

    // ---------- new regression tests for review fixes ----------

    /// B1: out-of-bounds allele index in a decoded genotype surfaces
    /// as the typed `AlleleIndexOutOfBounds` error instead of
    /// silently dropping the index.
    #[test]
    fn tally_called_alleles_errors_on_out_of_bounds_allele_index() {
        let record = biallelic_two_samples();
        // Pretend the writer was handed a 3-allele genotype table
        // while the record only carries 2 alleles; the entry
        // `[0, 2]` references allele 2 which doesn't exist on the
        // record.
        let table3 = genotype_order(2, 3); // [[0,0],[0,1],[1,1],[0,2],[1,2],[2,2]]
        let mut bad_record = record.clone();
        bad_record.best_genotype = vec![0, 3]; // index 3 is [0, 2] — allele 2 is out of bounds for n_alleles=2
        let err = tally_called_alleles(&bad_record, /*n_alleles=*/ 2, &table3).unwrap_err();
        match err {
            VcfWriteError::AlleleIndexOutOfBounds {
                sample_idx,
                allele_idx,
                n_alleles,
                ..
            } => {
                assert_eq!(sample_idx, 1);
                assert_eq!(allele_idx, 2);
                assert_eq!(n_alleles, 2);
            }
            other => panic!("expected AlleleIndexOutOfBounds, got {other:?}"),
        }
    }

    /// M6: empty allele set is rejected with `InconsistentRecord`
    /// instead of panicking on `record.alleles[0]`.
    #[test]
    fn encode_errors_on_empty_alleles() {
        let mut record = biallelic_two_samples();
        record.alleles.clear();
        record.allele_frequencies.clear();
        record.compound_frequencies.clear();
        record.scalars.clear();
        record.chain_anchor_flags.clear();
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let table = genotype_order(record.ploidy, 0);
        let err = encode(&record, &contigs, &config, &keys, 2, &table).unwrap_err();
        match err {
            VcfWriteError::InconsistentRecord {
                field,
                expected,
                actual,
                ..
            } => {
                assert_eq!(field, "alleles");
                assert_eq!(expected, 1);
                assert_eq!(actual, 0);
            }
            other => panic!("expected InconsistentRecord, got {other:?}"),
        }
    }

    /// M6: per-sample vectors must match the declared `n_samples`.
    #[test]
    fn encode_errors_on_undersized_best_genotype() {
        let mut record = biallelic_two_samples();
        record.best_genotype = vec![0]; // n_samples == 2 but len == 1
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let table = table_for(&record);
        let err = encode(&record, &contigs, &config, &keys, 2, &table).unwrap_err();
        match err {
            VcfWriteError::InconsistentRecord {
                field,
                expected,
                actual,
                ..
            } => {
                assert_eq!(field, "best_genotype");
                assert_eq!(expected, 2);
                assert_eq!(actual, 1);
            }
            other => panic!("expected InconsistentRecord, got {other:?}"),
        }
    }

    /// M7: per-sample DP overflowing `i32` surfaces as the typed
    /// `DepthOverflow` error rather than silently saturating.
    #[test]
    fn encode_errors_when_per_sample_depth_overflows_i32() {
        let mut record = biallelic_two_samples();
        // S0's two allele observations: each u32::MAX → sum > i32::MAX.
        record.scalars[0] = support(u32::MAX);
        record.scalars[1] = support(u32::MAX);
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let table = table_for(&record);
        let err = encode(&record, &contigs, &config, &keys, 2, &table).unwrap_err();
        match err {
            VcfWriteError::DepthOverflow {
                sample_idx, depth, ..
            } => {
                // Either the cohort-total (None) or the per-sample
                // path can fire first depending on iteration order;
                // both are valid signals of the same bug. Accept
                // both, assert non-zero depth.
                assert!(depth >= u32::MAX as u64);
                assert!(sample_idx.is_none() || sample_idx == Some(0));
            }
            other => panic!("expected DepthOverflow, got {other:?}"),
        }
    }

    /// M8: a single allele's `num_obs > i32::MAX` triggers
    /// `DepthOverflow` for the per-sample/per-allele AD cell rather
    /// than wrapping to a negative `i32`.
    #[test]
    fn encode_errors_when_single_allele_depth_overflows_i32() {
        let mut record = biallelic_two_samples();
        // One allele at u32::MAX (which is > i32::MAX), the other 0
        // — per-sample DP also overflows but the AD cell would
        // wrap. Either DP or AD path can fire; both are the same
        // typed variant.
        record.scalars[0] = support(u32::MAX);
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let table = table_for(&record);
        let err = encode(&record, &contigs, &config, &keys, 2, &table).unwrap_err();
        assert!(matches!(err, VcfWriteError::DepthOverflow { .. }));
    }

    /// M9: the writer's GQ cap must not be lower than the engine's
    /// `max_gq_phred`, or the writer will silently truncate
    /// confidence values the engine produces.
    #[test]
    fn gq_writer_cap_is_at_least_engine_cap() {
        use crate::var_calling::posterior_engine::PosteriorEngineConfig;
        let engine_cap = PosteriorEngineConfig::default().max_gq_phred;
        assert!(
            (GQ_MAX as f64) >= engine_cap,
            "writer GQ_MAX = {GQ_MAX} but engine max_gq_phred = {engine_cap}; \
             bumping the engine cap requires bumping GQ_MAX (and the matching \
             ##FORMAT description) in lockstep."
        );
    }

    /// Mi18: non-UTF-8 allele bytes surface as the typed encode error
    /// with the operation tag set.
    #[test]
    fn invalid_utf8_allele_surfaces_encode_error() {
        let mut record = biallelic_two_samples();
        record.alleles[1].seq = vec![0xff, 0xfe];
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let table = table_for(&record);
        let err = encode(&record, &contigs, &config, &keys, 2, &table).unwrap_err();
        match err {
            VcfWriteError::Encode { operation, .. } => {
                assert_eq!(operation, "allele bytes UTF-8");
            }
            other => panic!("expected Encode, got {other:?}"),
        }
    }

    /// Mi18: a `best_genotype` value past the genotype-order table
    /// length surfaces as `GenotypeIndexOutOfBounds` with the right
    /// sample and index recorded.
    #[test]
    fn out_of_range_best_genotype_surfaces_typed_error() {
        let mut record = biallelic_two_samples();
        record.best_genotype = vec![0, 99]; // table has length 3
        let contigs = fixture_contigs();
        let config = cfg_default();
        let keys = build_format_keys(&config);
        let table = table_for(&record);
        let err = encode(&record, &contigs, &config, &keys, 2, &table).unwrap_err();
        match err {
            VcfWriteError::GenotypeIndexOutOfBounds {
                sample_idx, got, ..
            } => {
                assert_eq!(sample_idx, 1);
                assert_eq!(got, 99);
            }
            other => panic!("expected GenotypeIndexOutOfBounds, got {other:?}"),
        }
    }
}
