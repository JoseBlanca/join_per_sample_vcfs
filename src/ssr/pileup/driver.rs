//! The Stage-1 `ssr-pileup` driver (Mark-2) — turns the sorted catalog + one
//! sample's alignment files into a per-locus `.ssr.psp` evidence file.
//!
//! The work for **one** locus is a single self-contained unit, [`process_locus`]
//! — fetch the locus's reads, delimit each (Viterbi+traceback) + quality-gate it,
//! tally the observed repeat-region sequences. It is pure over its shared `&`
//! inputs, so [`run`] processes the catalog in parallel batches (warm decode cache
//! per chunk), writing each batch in catalog order. Output is **byte-identical for
//! any thread count**: each locus's reservoir is seeded by `(chrom, start)`, scored
//! atomically on one thread, the delimiter tie-break is fixed, and the observed
//! sequences are stored sorted by bytes.
//!
//! This identity holds **across runs and thread counts on a fixed target +
//! toolchain**. It is *not* a cross-platform guarantee: the delimiter scores are
//! `f64` sums of transcendental log-probabilities, so a 1-ULP difference (a
//! different libm, fp-contraction, or reassociation on another target) could flip
//! a near-tie in the traceback and change the delimited bytes. Comparing or merging
//! `.ssr.psp` files produced on heterogeneous machines is therefore out of scope of
//! the guarantee.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use crate::bam::alignment_input::{PileupInputs, build_fasta_repository, load_pileup_inputs};
use crate::bam::errors::AlignmentInputError;
use crate::bam::segment_reader::{AlignmentFile, SegmentReadFilter, WorkerReader};
use crate::fasta::ContigList;
use crate::pop_var_caller::common::{
    DEFAULT_BUFFERED_IO_CAPACITY, basename, current_command_line, format_md5_hex, rfc3339_now,
};
use crate::psp::errors::PspWriteError;
use crate::psp::header::{ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance};
use crate::psp::registry_ssr::{SsrKind, SsrLocusRecord as ContainerRecord};
use crate::psp::writer::PspWriter;
use crate::ssr::catalog::CatalogError;
use crate::ssr::catalog::io::{CatalogHeader, CatalogReader};
use crate::ssr::types::Locus;

use rayon::prelude::*;

use super::alignment::{
    Delimited, HmmModel, MIN_REGION_Q1, ViterbiScratch, delimit_read, passes_quality_gate,
};
use super::fetch_reads::{LocusReads, fetch_locus_reads};
use super::footprint::{MIN_FLANK_BP, extract_region, read_footprint};
use super::locus_tally::{QcCounts, ReadObs, SsrLocusObs, tally};

/// Loci processed per parallel batch — bounds in-flight memory while staying large
/// enough to amortize scheduling.
const LOCUS_BATCH: usize = 8192;

/// Parallel-chunk granularity: each batch is split into contiguous chunks
/// (`par_chunks`), each owning one warm decode cache. The chunk **count** scales
/// with the worker count so a CRAM container is decoded ~once per chunk, bounding
/// the decode count as `--threads` grows.
const CHUNKS_PER_THREAD: usize = 4;
/// Floor on chunk size — never split below this (cold-cache boundaries for no gain).
const MIN_FETCH_CHUNK: usize = 64;

/// Errors from the Stage-1 driver.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub(crate) enum SsrPileupError {
    #[error("reading alignment input failed")]
    Read(#[from] AlignmentInputError),
    #[error("reading the catalog failed")]
    Catalog(#[from] CatalogError),
    #[error("writing the .ssr.psp failed")]
    Write(#[from] PspWriteError),
    #[error("failed to open the catalog {path:?}")]
    OpenCatalog {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("failed to create the temporary output {path:?}")]
    CreateOutput {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("failed to flush the .ssr.psp output")]
    FlushOutput {
        #[source]
        source: std::io::Error,
    },
    #[error("failed to sync the output {path:?}")]
    SyncOutput {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("failed to rename {tmp:?} into {dest:?}")]
    RenameOutput {
        tmp: PathBuf,
        dest: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("locus contig {name:?} is not in the reference / alignment header")]
    ContigNotInReference { name: String },
    #[error("contig {name:?} has no md5 in the reference (required for the .psp header)")]
    MissingMd5 { name: String },
    #[error("contig {name:?} length {length} exceeds the .psp u32 limit")]
    ContigLengthOverflow { name: String, length: u64 },
    #[error("could not parse the run timestamp {value:?}")]
    TimestampFormat {
        value: String,
        #[source]
        source: toml::value::DatetimeParseError,
    },
    #[error("failed to build the worker thread pool")]
    ThreadPool(#[from] rayon::ThreadPoolBuildError),
}

/// Per-worker scratch reused across loci — the Viterbi DP workspace the delimiter
/// writes through. **One per worker thread** (the parallel driver hands each
/// thread its own).
pub(crate) struct LocusScratch {
    viterbi: ViterbiScratch,
}

impl LocusScratch {
    pub(crate) fn new() -> Self {
        Self {
            viterbi: ViterbiScratch::new(),
        }
    }
}

/// Assemble the locus's QC scalars from the fetch pass. `depth` = usable primary
/// reads considered; `n_filtered` = the gate's quality drops + duplicates +
/// length-malformed reads; `mapped_reads` = the dup-free coverage denominator.
/// (Decided in the Mark-1 driver doc §4/§7.3, carried over.) Counts are u64 →
/// u32; the adds saturate and the casts clamp so a pathologically deep locus
/// records `u32::MAX` rather than wrapping silently.
pub(crate) fn qc_counts(fetched: &LocusReads) -> QcCounts {
    let d = &fetched.filtered;
    let primary_quality_drops = d
        .qc_fail
        .saturating_add(d.low_mapq)
        .saturating_add(d.too_short);
    let n_filtered = primary_quality_drops
        .saturating_add(d.duplicate)
        .saturating_add(fetched.malformed);
    QcCounts {
        depth: u32::try_from(fetched.yielded).unwrap_or(u32::MAX),
        n_filtered: u32::try_from(n_filtered).unwrap_or(u32::MAX),
        mapped_reads: u32::try_from(fetched.yielded.saturating_add(primary_quality_drops))
            .unwrap_or(u32::MAX),
    }
}

/// Process one locus end to end: fetch + depth-cap its reads, delimit each
/// (Viterbi+traceback) and quality-gate the repeat region, and tally the observed
/// sequences into the in-memory (chrom-**name**-keyed) [`SsrLocusObs`].
/// Concurrency-safe for different loci: the only `&mut` is the caller-owned
/// `scratch`.
pub(crate) fn process_locus(
    readers: &mut [WorkerReader<'_>],
    locus: &Locus,
    cap: usize,
    model: &HmmModel,
    scratch: &mut LocusScratch,
) -> Result<SsrLocusObs, SsrPileupError> {
    let fetched = fetch_locus_reads(readers, locus, cap)?;
    let mut outcomes = Vec::with_capacity(fetched.reads.len());
    for read in &fetched.reads {
        let fp = read_footprint(&read.cigar, read.pos);
        let region = extract_region(&read.cigar, fp, read.seq.len(), locus);
        let region_seq = &read.seq[region.clone()];
        let region_qual = &read.qual[region];
        let outcome =
            match delimit_read(region_seq, region_qual, locus, model, &mut scratch.viterbi) {
                Delimited::Region(r) => {
                    if passes_quality_gate(
                        &region_qual[r.clone()],
                        MIN_REGION_Q1,
                        &mut scratch.viterbi,
                    ) {
                        ReadObs::Sequence(region_seq[r].to_vec().into_boxed_slice())
                    } else {
                        ReadObs::LowQuality
                    }
                }
                Delimited::BorderOffEnd => ReadObs::BorderOffEnd,
            };
        outcomes.push(outcome);
    }
    Ok(tally(locus, &outcomes, qc_counts(&fetched)))
}

/// Adapt the in-memory (chrom-**name**-keyed, 0-based) [`SsrLocusObs`] to the
/// container's (chrom-**id**-keyed, 1-based) record. Both coordinate bounds shift
/// by `+1`; `n_obs` is derived from `observed.len()` by the container, so it is
/// not carried here.
pub(crate) fn to_container_record(
    record: SsrLocusObs,
    name_to_id: &HashMap<&str, u32>,
) -> Result<ContainerRecord, SsrPileupError> {
    let chrom_id = *name_to_id.get(record.chrom.as_ref()).ok_or_else(|| {
        SsrPileupError::ContigNotInReference {
            name: record.chrom.to_string(),
        }
    })?;
    Ok(ContainerRecord {
        chrom_id,
        start: record.start + 1, // 0-based inclusive -> 1-based inclusive
        end: record.end + 1,     // 0-based exclusive -> 1-based exclusive
        depth: record.depth,
        n_filtered: record.n_filtered,
        mapped_reads: record.mapped_reads,
        n_low_quality: record.n_low_quality,
        n_border_off_end: record.n_border_off_end,
        observed: record.observed,
    })
}

/// Inputs for one `ssr-pileup` run.
pub(crate) struct SsrPileupConfig {
    /// One or more coordinate-sorted, indexed BAM/CRAM files for a single sample.
    pub alignment_files: Vec<PathBuf>,
    /// Reference FASTA (with sibling `.fai`).
    pub reference: PathBuf,
    /// The sorted `.ssr.catalog` to genotype against.
    pub catalog: PathBuf,
    /// Destination `.ssr.psp`.
    pub output: PathBuf,
    /// The cheap read filter (MAPQ / dup / qc-fail / length).
    pub filter: SegmentReadFilter,
    /// Per-locus reservoir depth cap.
    pub cap: usize,
    /// Build a missing alignment index in place rather than erroring.
    pub build_index_if_missing: bool,
    /// Sample name override; defaults to the one the inputs cross-validate.
    pub sample: Option<String>,
    /// Worker threads for the per-locus pool. Output is identical for any value.
    pub threads: usize,
}

/// Build the `.ssr.psp` writer header — one [`ChromosomeEntry`] per reference
/// contig, `subcommand = "ssr-pileup"`, and the run parameters (including the
/// catalog binding: its `reference_md5` and `flank_bp`).
fn build_ssr_writer_header(
    sample: &str,
    cfg: &SsrPileupConfig,
    contigs: &ContigList,
    catalog: &CatalogHeader,
) -> Result<WriterHeader, SsrPileupError> {
    let mut chromosomes = Vec::with_capacity(contigs.entries.len());
    for entry in &contigs.entries {
        let length =
            u32::try_from(entry.length).map_err(|_| SsrPileupError::ContigLengthOverflow {
                name: entry.name.clone(),
                length: entry.length,
            })?;
        let md5 = entry
            .md5
            .map(format_md5_hex)
            .ok_or_else(|| SsrPileupError::MissingMd5 {
                name: entry.name.clone(),
            })?;
        chromosomes.push(ChromosomeEntry {
            name: entry.name.clone(),
            length,
            md5,
        });
    }

    let now = rfc3339_now();
    let created: toml::value::Datetime =
        now.parse()
            .map_err(|source| SsrPileupError::TimestampFormat {
                value: now.clone(),
                source,
            })?;

    let input_alignment_files: Vec<String> =
        cfg.alignment_files.iter().map(|p| basename(p)).collect();
    let input_fasta = basename(&cfg.reference);

    let mut parameters = BTreeMap::new();
    parameters.insert(
        "quality_q1_threshold".to_string(),
        ParameterValue::Integer(i64::from(MIN_REGION_Q1)),
    );
    parameters.insert(
        "reservoir_cap".to_string(),
        ParameterValue::Integer(i64::try_from(cfg.cap).unwrap_or(i64::MAX)),
    );
    parameters.insert(
        "reach_min_flank_bp".to_string(),
        ParameterValue::Integer(i64::from(MIN_FLANK_BP)),
    );
    parameters.insert(
        "flank_bp".to_string(),
        ParameterValue::Integer(i64::from(catalog.params.flank_bp)),
    );
    parameters.insert(
        "catalog_reference_md5".to_string(),
        ParameterValue::String(catalog.reference_md5.clone()),
    );
    if let Some(mapq) = cfg.filter.min_mapq {
        parameters.insert(
            "min_mapq".to_string(),
            ParameterValue::Integer(i64::from(mapq)),
        );
    }
    if let Some(len) = cfg.filter.min_read_length {
        parameters.insert(
            "min_read_length".to_string(),
            ParameterValue::Integer(i64::from(len)),
        );
    }
    parameters.insert(
        "drop_duplicate".to_string(),
        ParameterValue::Boolean(cfg.filter.drop_duplicate),
    );
    parameters.insert(
        "drop_qc_fail".to_string(),
        ParameterValue::Boolean(cfg.filter.drop_qc_fail),
    );

    Ok(WriterHeader {
        format_version: (1, 0),
        sample: sample.to_string(),
        reference: input_fasta.clone(),
        created,
        chromosomes,
        writer: WriterProvenance {
            tool: "pop_var_caller".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            subcommand: "ssr-pileup".to_string(),
            input_crams: input_alignment_files,
            input_fasta,
            command_line: current_command_line(),
            parameters,
        },
    })
}

/// `.tmp` sibling of `output`, written then atomically renamed into place.
fn tmp_path(output: &Path) -> PathBuf {
    let mut s = output.to_path_buf().into_os_string();
    s.push(".tmp");
    PathBuf::from(s)
}

/// Run Stage 1: walk the sorted catalog, genotype each locus from the sample's
/// alignment files, and write the per-locus evidence to a `.ssr.psp`. Output is
/// written to a temp file and atomically renamed.
pub(crate) fn run(cfg: &SsrPileupConfig) -> Result<(), SsrPileupError> {
    let PileupInputs {
        sample_name,
        contigs,
        headers,
        indexes,
    } = load_pileup_inputs(
        &cfg.alignment_files,
        &cfg.reference,
        cfg.build_index_if_missing,
    )?;
    let repo = build_fasta_repository(&cfg.reference)?;

    // One pooled reader per input file (built once; shared across all loci).
    let files: Vec<AlignmentFile> = headers
        .into_iter()
        .zip(indexes)
        .enumerate()
        .map(|(i, (header, index))| {
            AlignmentFile::from_input(
                cfg.alignment_files[i].clone(),
                header,
                index,
                Some(repo.clone()),
                cfg.filter,
                i,
            )
        })
        .collect::<Result<_, _>>()?;

    let sample = cfg.sample.clone().unwrap_or(sample_name);
    let catalog_file = File::open(&cfg.catalog).map_err(|source| SsrPileupError::OpenCatalog {
        path: cfg.catalog.clone(),
        source,
    })?;
    let mut catalog = CatalogReader::new(catalog_file)?;
    let header = build_ssr_writer_header(&sample, cfg, &contigs, catalog.header())?;

    let name_to_id: HashMap<&str, u32> = contigs
        .entries
        .iter()
        .enumerate()
        .map(|(i, e)| (e.name.as_str(), i as u32))
        .collect();

    let tmp = tmp_path(&cfg.output);
    let output_file = File::create(&tmp).map_err(|source| SsrPileupError::CreateOutput {
        path: tmp.clone(),
        source,
    })?;
    let writer_sink = BufWriter::with_capacity(DEFAULT_BUFFERED_IO_CAPACITY, output_file);
    let mut writer = PspWriter::<_, SsrKind>::new_ssr(writer_sink, header)?;

    let model = HmmModel::new();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cfg.threads)
        .build()?;

    let mut batch: Vec<Locus> = Vec::with_capacity(LOCUS_BATCH);
    loop {
        batch.clear();
        while batch.len() < LOCUS_BATCH {
            match catalog.read_locus() {
                Some(locus) => batch.push(locus?),
                None => break,
            }
        }
        if batch.is_empty() {
            break;
        }
        let n_chunks = (cfg.threads.max(1) * CHUNKS_PER_THREAD).max(1);
        let chunk_size = batch.len().div_ceil(n_chunks).max(MIN_FETCH_CHUNK);
        let chunks: Vec<Vec<ContainerRecord>> = pool.install(|| {
            batch
                .par_chunks(chunk_size)
                .map(|chunk| -> Result<Vec<ContainerRecord>, SsrPileupError> {
                    let mut readers: Vec<WorkerReader> =
                        files.iter().map(|f| f.worker_reader()).collect();
                    let mut scratch = LocusScratch::new();
                    let mut out = Vec::with_capacity(chunk.len());
                    for locus in chunk {
                        let record =
                            process_locus(&mut readers, locus, cfg.cap, &model, &mut scratch)?;
                        out.push(to_container_record(record, &name_to_id)?);
                    }
                    Ok(out)
                })
                .collect::<Result<Vec<_>, SsrPileupError>>()
        })?;
        // Chunks are in catalog order (`par_chunks` is indexed); within a chunk
        // records are sequential — so this writes in catalog order.
        for chunk_records in &chunks {
            for record in chunk_records {
                writer.write_locus(record)?;
            }
        }
    }

    let file = writer
        .finish()?
        .into_inner()
        .map_err(|e| SsrPileupError::FlushOutput {
            source: e.into_error(),
        })?;
    file.sync_all()
        .map_err(|source| SsrPileupError::SyncOutput {
            path: tmp.clone(),
            source,
        })?;
    std::fs::rename(&tmp, &cfg.output).map_err(|source| SsrPileupError::RenameOutput {
        tmp: tmp.clone(),
        dest: cfg.output.clone(),
        source,
    })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::alignment_input::FilterCounts;

    fn locus_reads(yielded: u64, malformed: u64, filtered: FilterCounts) -> LocusReads {
        LocusReads {
            reads: Vec::new(),
            yielded,
            malformed,
            filtered,
        }
    }

    #[test]
    fn qc_counts_excludes_dups_from_coverage_but_keeps_them_filtered() {
        // Every field named (no `..default()`) so a new FilterCounts bucket
        // forces a decision here rather than being silently absorbed.
        let filtered = FilterCounts {
            unmapped: 2,
            secondary: 5,
            supplementary: 4,
            qc_fail: 2,
            duplicate: 3,
            low_mapq: 1,
            too_short: 1,
            high_mismatch_fraction: 0,
            bad_cigar: 0,
            baq_rejected: 0,
        };
        let qc = qc_counts(&locus_reads(10, 0, filtered));
        assert_eq!(qc.depth, 10);
        assert_eq!(qc.n_filtered, 4 + 3); // quality drops + duplicates
        assert_eq!(qc.mapped_reads, 10 + 4); // depth + quality drops (dup-free)
    }

    #[test]
    fn qc_counts_folds_malformed_reads_into_n_filtered() {
        // Malformed (length-inconsistent) reads are counted in n_filtered, not
        // in depth/mapped_reads (B2 — they never enter the analyzed set).
        let qc = qc_counts(&locus_reads(8, 3, FilterCounts::default()));
        assert_eq!(qc.depth, 8);
        assert_eq!(qc.n_filtered, 3); // only the malformed reads
        assert_eq!(qc.mapped_reads, 8);
    }

    fn in_memory_obs(chrom: &str, start: u32, end: u32) -> SsrLocusObs {
        SsrLocusObs {
            chrom: chrom.into(),
            start,
            end,
            depth: 5,
            n_filtered: 3,
            mapped_reads: 8,
            n_low_quality: 1,
            n_border_off_end: 2,
            observed: vec![(b"CACACA".to_vec().into_boxed_slice(), 4)],
        }
    }

    #[test]
    fn adapter_maps_name_to_id_and_shifts_coords_to_1_based() {
        let name_to_id = HashMap::from([("chr1", 0u32), ("chr2", 1)]);
        let c = to_container_record(in_memory_obs("chr2", 16, 22), &name_to_id).unwrap();
        assert_eq!(c.chrom_id, 1);
        assert_eq!((c.start, c.end), (17, 23)); // 0-based [16,22) -> 1-based [17,23)
        assert_eq!(
            (
                c.depth,
                c.n_filtered,
                c.mapped_reads,
                c.n_low_quality,
                c.n_border_off_end
            ),
            (5, 3, 8, 1, 2)
        );
        assert_eq!(c.observed, vec![(b"CACACA".to_vec().into_boxed_slice(), 4)]);
    }

    #[test]
    fn adapter_errors_on_a_contig_absent_from_the_reference() {
        let name_to_id = HashMap::from([("chr1", 0u32)]);
        let err = to_container_record(in_memory_obs("chrX", 1, 5), &name_to_id).unwrap_err();
        assert!(matches!(
            err,
            SsrPileupError::ContigNotInReference { name } if name == "chrX"
        ));
    }

    #[test]
    fn build_ssr_writer_header_errors_on_overflow_and_missing_md5() {
        use crate::fasta::{ContigEntry, ContigList};
        let cfg = SsrPileupConfig {
            alignment_files: vec![PathBuf::from("s.bam")],
            reference: PathBuf::from("ref.fa"),
            catalog: PathBuf::from("c.ssr.catalog"),
            output: PathBuf::from("o.ssr.psp"),
            filter: SegmentReadFilter::default(),
            cap: 1000,
            build_index_if_missing: false,
            sample: None,
            threads: 1,
        };
        let cat = catalog_header();

        // (a) a contig longer than the .psp u32 limit -> ContigLengthOverflow.
        let too_long = u64::from(u32::MAX) + 1;
        let big = ContigList {
            entries: vec![ContigEntry {
                name: "chr1".into(),
                length: too_long,
                md5: Some([0u8; 16]),
            }],
        };
        let err = build_ssr_writer_header("S", &cfg, &big, &cat).unwrap_err();
        assert!(matches!(
            err,
            SsrPileupError::ContigLengthOverflow { ref name, length }
                if name == "chr1" && length == too_long
        ));

        // (b) a contig without an md5 -> MissingMd5.
        let no_md5 = ContigList {
            entries: vec![ContigEntry {
                name: "chr2".into(),
                length: 100,
                md5: None,
            }],
        };
        let err = build_ssr_writer_header("S", &cfg, &no_md5, &cat).unwrap_err();
        assert!(matches!(err, SsrPileupError::MissingMd5 { name } if name == "chr2"));

        // (c) a valid contig records the run parameters (incl. reach_min_flank_bp).
        let ok = ContigList {
            entries: vec![ContigEntry {
                name: "chr1".into(),
                length: 100,
                md5: Some([1u8; 16]),
            }],
        };
        let header = build_ssr_writer_header("S", &cfg, &ok, &cat).unwrap();
        assert_eq!(header.chromosomes.len(), 1);
        for key in [
            "quality_q1_threshold",
            "reservoir_cap",
            "reach_min_flank_bp",
        ] {
            assert!(
                header.writer.parameters.contains_key(key),
                "header parameters should record {key}"
            );
        }
    }

    // --- end-to-end: catalog + reference + BAM -> .ssr.psp ------------

    use crate::ssr::catalog::CatalogParams;
    use std::num::NonZero;
    use tempfile::TempDir;

    fn catalog_header() -> CatalogHeader {
        CatalogHeader {
            tool_version: "0.1.0".into(),
            reference: "ref.fa".into(),
            reference_md5: "a".repeat(32),
            trf_mod_version: "1.0".into(),
            params: CatalogParams::default(),
            date: "2026-06-17".into(),
        }
    }

    use noodles_sam::alignment::RecordBuf;

    /// One clean CA(3) spanning read for the locus at `start`. Aligned at ref
    /// [start-11, start+19) (pos = start-10, M30): the window ref [start-6,
    /// start+12) maps to read[5..23] = 6 bp G-flank + CACACA + 6 bp T-flank.
    /// `quals` is the 30-base quality string (lets a caller dim the tract).
    fn spanning_read(start: u32, qname: &str, quals: Vec<u8>) -> RecordBuf {
        use noodles_core::Position;
        use noodles_sam::alignment::record::cigar::Op;
        use noodles_sam::alignment::record::cigar::op::Kind;
        use noodles_sam::alignment::record::{Flags, MappingQuality};
        use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AAAAA");
        seq.extend_from_slice(b"GGGGGGCACACATTTTTT");
        seq.extend_from_slice(b"AAAAAAA");
        RecordBuf::builder()
            .set_name(qname.as_bytes())
            .set_reference_sequence_id(0)
            .set_flags(Flags::default())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_alignment_start(Position::try_from((start - 10) as usize).unwrap())
            .set_cigar([Op::new(Kind::Match, 30)].into_iter().collect())
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(quals))
            .build()
    }

    /// Write a coordinate-sorted single-contig BAM (`@RG SM`, a trusted `@SQ M5`
    /// of length `contig_len`) holding `reads` (caller supplies them in
    /// coordinate order).
    fn write_bam(path: &std::path::Path, contig_len: usize, reads: &[RecordBuf]) {
        use bstr::BString;
        use noodles_bam as bam;
        use noodles_sam as sam;
        use sam::alignment::io::Write as _;
        use sam::header::record::value::Map;
        use sam::header::record::value::map::header::Version;
        use sam::header::record::value::map::header::tag::SORT_ORDER;
        use sam::header::record::value::map::read_group::tag::SAMPLE;
        use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
        use sam::header::record::value::map::{Header as HeaderMap, ReadGroup, ReferenceSequence};

        let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
        hd.other_fields_mut()
            .insert(SORT_ORDER, BString::from("coordinate"));
        let mut sq = Map::<ReferenceSequence>::new(NonZero::new(contig_len).unwrap());
        sq.other_fields_mut()
            .insert(MD5_CHECKSUM, BString::from("0".repeat(32)));
        let mut rg = Map::<ReadGroup>::default();
        rg.other_fields_mut()
            .insert(SAMPLE, BString::from("sample0"));
        let header = sam::Header::builder()
            .set_header(hd)
            .add_reference_sequence(b"chr1".to_vec(), sq)
            .add_read_group(b"rg0".to_vec(), rg)
            .build();

        let mut writer = bam::io::Writer::new(File::create(path).unwrap());
        writer.write_header(&header).unwrap();
        for read in reads {
            writer.write_alignment_record(&header, read).unwrap();
        }
        writer.try_finish().unwrap();
    }

    /// Build a fixture with one CA(3) catalog locus per (ascending) tract start
    /// and `reps` identical clean spanning reads per locus. `reps = 1` is the
    /// one-read-per-locus case.
    fn stage1_fixture_reps(
        tract_starts: &[u32],
        reps: usize,
    ) -> (TempDir, TempDir, PathBuf, PathBuf, PathBuf) {
        use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};
        use crate::ssr::catalog::io::CatalogWriter;
        use crate::ssr::types::Motif;

        let contig_len = tract_starts.iter().copied().max().unwrap_or(0) as usize + 20;
        let (fa_dir, fasta) = build_fasta(&[ContigSpec {
            name: "chr1".into(),
            length: contig_len as u64,
        }])
        .unwrap();
        let dir = TempDir::new().unwrap();
        let bam = dir.path().join("s.bam");
        let mut reads = Vec::new();
        for &start in tract_starts {
            for i in 0..reps {
                reads.push(spanning_read(
                    start,
                    &format!("r{start}_{i}"),
                    vec![40u8; 30],
                ));
            }
        }
        write_bam(&bam, contig_len, &reads);

        let catalog = dir.path().join("c.ssr.catalog");
        let mut w = CatalogWriter::new(File::create(&catalog).unwrap(), &catalog_header()).unwrap();
        for &start in tract_starts {
            let locus = Locus::new(
                "chr1".into(),
                start,
                start + 6,
                Motif::new(b"CA").unwrap(),
                1.0,
                (*b"GGGGGGCACACATTTTTT").into(),
                start - 6,
            )
            .unwrap();
            w.write_locus(&locus).unwrap();
        }
        w.finish().unwrap();
        (fa_dir, dir, bam, fasta, catalog)
    }

    fn stage1_fixture(tract_starts: &[u32]) -> (TempDir, TempDir, PathBuf, PathBuf, PathBuf) {
        stage1_fixture_reps(tract_starts, 1)
    }

    fn read_records(path: &std::path::Path) -> Vec<ContainerRecord> {
        use crate::psp::PspReader;
        let mut reader = PspReader::new(File::open(path).unwrap()).unwrap();
        assert_eq!(reader.header().kind, "ssr");
        reader.records_of::<SsrKind>().map(|r| r.unwrap()).collect()
    }

    fn run_config_with_cap(
        bam: &std::path::Path,
        fasta: &std::path::Path,
        catalog: &std::path::Path,
        output: &std::path::Path,
        threads: usize,
        cap: usize,
    ) -> SsrPileupConfig {
        SsrPileupConfig {
            alignment_files: vec![bam.to_path_buf()],
            reference: fasta.to_path_buf(),
            catalog: catalog.to_path_buf(),
            output: output.to_path_buf(),
            filter: SegmentReadFilter::default(),
            cap,
            build_index_if_missing: true,
            sample: Some("S1".into()),
            threads,
        }
    }

    fn run_config(
        bam: &std::path::Path,
        fasta: &std::path::Path,
        catalog: &std::path::Path,
        output: &std::path::Path,
        threads: usize,
    ) -> SsrPileupConfig {
        use crate::ssr::pileup::fetch_reads::MAX_READS_PER_LOCUS;
        run_config_with_cap(bam, fasta, catalog, output, threads, MAX_READS_PER_LOCUS)
    }

    #[test]
    fn run_genotypes_one_locus_end_to_end() {
        let (_fa, dir, bam, fasta, catalog) = stage1_fixture(&[16]);
        let output = dir.path().join("out.ssr.psp");
        run(&run_config(&bam, &fasta, &catalog, &output, 1)).expect("run");

        let records = read_records(&output);
        assert_eq!(records.len(), 1);
        let rec = &records[0];
        assert_eq!(rec.chrom_id, 0);
        assert_eq!((rec.start, rec.end), (17, 23)); // 0-based [16,22) -> 1-based [17,23)
        assert_eq!(rec.depth, 1);
        assert_eq!(rec.mapped_reads, 1);
        assert_eq!(rec.n_filtered, 0);
        // One spanning read showing the clean reference tract CACACA, count 1.
        assert_eq!(
            rec.observed,
            vec![(b"CACACA".to_vec().into_boxed_slice(), 1)]
        );
    }

    #[test]
    fn run_output_is_thread_count_invariant() {
        let (_fa, dir, bam, fasta, catalog) = stage1_fixture(&[16, 66, 116]);
        let out1 = dir.path().join("t1.ssr.psp");
        let out3 = dir.path().join("t3.ssr.psp");
        run(&run_config(&bam, &fasta, &catalog, &out1, 1)).expect("run @ 1 thread");
        run(&run_config(&bam, &fasta, &catalog, &out3, 3)).expect("run @ 3 threads");

        let r1 = read_records(&out1);
        let r3 = read_records(&out3);
        assert_eq!(r1.len(), 3);
        assert_eq!(r1, r3); // identical records + order regardless of thread count
    }

    #[test]
    fn run_output_is_thread_invariant_across_multiple_chunks() {
        // > MIN_FETCH_CHUNK (64) loci so a batch splits into >1 `par_chunks`
        // chunk; the ordered collect/write must preserve catalog order at any
        // thread count (B3 — the 3-locus fixture is always a single chunk).
        let starts: Vec<u32> = (0..80u32).map(|i| 16 + 30 * i).collect();
        let (_fa, dir, bam, fasta, catalog) = stage1_fixture(&starts);
        let out1 = dir.path().join("t1.ssr.psp");
        let out8 = dir.path().join("t8.ssr.psp");
        run(&run_config(&bam, &fasta, &catalog, &out1, 1)).expect("run @ 1 thread");
        run(&run_config(&bam, &fasta, &catalog, &out8, 8)).expect("run @ 8 threads");

        let r1 = read_records(&out1);
        let r8 = read_records(&out8);
        assert_eq!(r1.len(), 80);
        assert_eq!(r1, r8); // identical records + order across the chunk split
        // Strictly ascending start => chunks were emitted in catalog order.
        assert!(r1.windows(2).all(|w| w[0].start < w[1].start));
    }

    #[test]
    fn run_output_is_thread_invariant_when_the_reservoir_cap_bites() {
        // 50 spanning reads at one locus, cap 4: the reservoir eviction branch
        // runs on the end-to-end path; depth counts every read, observed caps at
        // the cap, and the record is identical at 1 vs 4 threads (B3).
        let (_fa, dir, bam, fasta, catalog) = stage1_fixture_reps(&[16], 50);
        let out1 = dir.path().join("t1.ssr.psp");
        let out4 = dir.path().join("t4.ssr.psp");
        run(&run_config_with_cap(&bam, &fasta, &catalog, &out1, 1, 4)).expect("run @ 1 thread");
        run(&run_config_with_cap(&bam, &fasta, &catalog, &out4, 4, 4)).expect("run @ 4 threads");

        let r1 = read_records(&out1);
        assert_eq!(r1, read_records(&out4)); // cap evicted identically across threads
        assert_eq!(r1.len(), 1);
        assert_eq!(r1[0].depth, 50);
        assert_eq!(
            r1[0].observed,
            vec![(b"CACACA".to_vec().into_boxed_slice(), 4)]
        );
    }

    #[test]
    fn process_locus_routes_sequence_and_low_quality_outcomes() {
        // One locus, two spanning reads: a clean one (Q40) and one whose tract
        // bases are below the Phred-15 gate. The clean read tallies as an
        // observed sequence; the dim one routes to n_low_quality (M9).
        use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};
        use crate::ssr::catalog::io::CatalogWriter;
        use crate::ssr::types::Motif;

        let contig_len: usize = 36;
        let (_fa, fasta) = build_fasta(&[ContigSpec {
            name: "chr1".into(),
            length: contig_len as u64,
        }])
        .unwrap();
        let dir = TempDir::new().unwrap();
        let bam = dir.path().join("s.bam");
        // read[11..17] is the tract CACACA; dim those bases below Phred 15.
        let mut low_q = vec![40u8; 30];
        for q in &mut low_q[11..17] {
            *q = 10;
        }
        let reads = vec![
            spanning_read(16, "clean", vec![40u8; 30]),
            spanning_read(16, "lowq", low_q),
        ];
        write_bam(&bam, contig_len, &reads);

        let catalog = dir.path().join("c.ssr.catalog");
        let mut w = CatalogWriter::new(File::create(&catalog).unwrap(), &catalog_header()).unwrap();
        w.write_locus(
            &Locus::new(
                "chr1".into(),
                16,
                22,
                Motif::new(b"CA").unwrap(),
                1.0,
                (*b"GGGGGGCACACATTTTTT").into(),
                10,
            )
            .unwrap(),
        )
        .unwrap();
        w.finish().unwrap();

        let output = dir.path().join("out.ssr.psp");
        run(&run_config(&bam, &fasta, &catalog, &output, 1)).expect("run");
        let records = read_records(&output);
        assert_eq!(records.len(), 1);
        let rec = &records[0];
        assert_eq!(rec.depth, 2);
        assert_eq!(rec.n_low_quality, 1);
        assert_eq!(rec.n_border_off_end, 0);
        assert_eq!(
            rec.observed,
            vec![(b"CACACA".to_vec().into_boxed_slice(), 1)]
        );
    }
}
