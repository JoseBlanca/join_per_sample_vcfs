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
use super::footprint::{extract_region, read_footprint};
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
    #[error("I/O error")]
    Io(#[from] std::io::Error),
    #[error("locus contig {name:?} is not in the reference / alignment header")]
    ContigNotInReference { name: String },
    #[error("contig {name:?} has no md5 in the reference (required for the .psp header)")]
    MissingMd5 { name: String },
    #[error("contig {name:?} length {length} exceeds the .psp u32 limit")]
    ContigLengthOverflow { name: String, length: u64 },
    #[error("could not format the run timestamp")]
    TimestampFormat,
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
/// reads considered; `n_filtered` = the gate's quality drops + duplicates;
/// `mapped_reads` = the dup-free coverage denominator. (Decided in the Mark-1
/// driver doc §4/§7.3, carried over.)
pub(crate) fn qc_counts(fetched: &LocusReads) -> QcCounts {
    let d = &fetched.filtered;
    let primary_quality_drops = d.qc_fail + d.low_mapq + d.too_short;
    QcCounts {
        depth: fetched.yielded as u32,
        n_filtered: (primary_quality_drops + d.duplicate) as u32,
        mapped_reads: (fetched.yielded + primary_quality_drops) as u32,
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

    let created: toml::value::Datetime = rfc3339_now()
        .parse()
        .map_err(|_| SsrPileupError::TimestampFormat)?;

    let input_crams: Vec<String> = cfg.alignment_files.iter().map(|p| basename(p)).collect();
    let input_fasta = basename(&cfg.reference);

    let mut parameters = BTreeMap::new();
    parameters.insert(
        "quality_q1_threshold".to_string(),
        ParameterValue::Integer(i64::from(MIN_REGION_Q1)),
    );
    parameters.insert(
        "reservoir_cap".to_string(),
        ParameterValue::Integer(cfg.cap as i64),
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
            input_crams,
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
    let mut catalog = CatalogReader::new(File::open(&cfg.catalog)?)?;
    let header = build_ssr_writer_header(&sample, cfg, &contigs, catalog.header())?;

    let name_to_id: HashMap<&str, u32> = contigs
        .entries
        .iter()
        .enumerate()
        .map(|(i, e)| (e.name.as_str(), i as u32))
        .collect();

    let tmp = tmp_path(&cfg.output);
    let writer_sink = BufWriter::with_capacity(DEFAULT_BUFFERED_IO_CAPACITY, File::create(&tmp)?);
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
        .map_err(|e| SsrPileupError::Io(e.into_error()))?;
    file.sync_all()?;
    std::fs::rename(&tmp, &cfg.output)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::alignment_input::FilterCounts;

    fn locus_reads(yielded: u64, filtered: FilterCounts) -> LocusReads {
        LocusReads {
            reads: Vec::new(),
            yielded,
            filtered,
        }
    }

    #[test]
    fn qc_counts_excludes_dups_from_coverage_but_keeps_them_filtered() {
        let filtered = FilterCounts {
            qc_fail: 2,
            low_mapq: 1,
            too_short: 1,
            duplicate: 3,
            secondary: 5,
            supplementary: 4,
            unmapped: 2,
            ..FilterCounts::default()
        };
        let qc = qc_counts(&locus_reads(10, filtered));
        assert_eq!(qc.depth, 10);
        assert_eq!(qc.n_filtered, 4 + 3); // quality drops + duplicates
        assert_eq!(qc.mapped_reads, 10 + 4); // depth + quality drops (dup-free)
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

    /// A coordinate-sorted BAM (one contig, `@RG SM`, a trusted `@SQ M5`) with one
    /// clean spanning read per CA(3) locus at the given ascending tract starts.
    fn write_bam_for_loci(path: &std::path::Path, tract_starts: &[u32]) {
        use bstr::BString;
        use noodles_bam as bam;
        use noodles_sam as sam;
        use sam::alignment::RecordBuf;
        use sam::alignment::io::Write as _;
        use sam::alignment::record::cigar::Op;
        use sam::alignment::record::cigar::op::Kind;
        use sam::alignment::record::{Flags, MappingQuality};
        use sam::alignment::record_buf::{QualityScores, Sequence};
        use sam::header::record::value::Map;
        use sam::header::record::value::map::header::Version;
        use sam::header::record::value::map::header::tag::SORT_ORDER;
        use sam::header::record::value::map::read_group::tag::SAMPLE;
        use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
        use sam::header::record::value::map::{Header as HeaderMap, ReadGroup, ReferenceSequence};

        let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
        hd.other_fields_mut()
            .insert(SORT_ORDER, BString::from("coordinate"));
        let mut sq = Map::<ReferenceSequence>::new(NonZero::new(200usize).unwrap());
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
        for &start in tract_starts {
            // Tract [start, start+6); window ref [start-6, start+12). Read aligned
            // at ref [start-11, start+19) (pos = start-10, M30): window maps to
            // read[5..23] = 6 bp G-flank + CACACA + 6 bp T-flank.
            let mut seq = Vec::new();
            seq.extend_from_slice(b"AAAAA");
            seq.extend_from_slice(b"GGGGGGCACACATTTTTT");
            seq.extend_from_slice(b"AAAAAAA");
            let read = RecordBuf::builder()
                .set_name(format!("r{start}").as_bytes())
                .set_reference_sequence_id(0)
                .set_flags(Flags::default())
                .set_mapping_quality(MappingQuality::new(60).unwrap())
                .set_alignment_start(
                    noodles_core::Position::try_from((start - 10) as usize).unwrap(),
                )
                .set_cigar([Op::new(Kind::Match, 30)].into_iter().collect())
                .set_sequence(Sequence::from(seq))
                .set_quality_scores(QualityScores::from(vec![40u8; 30]))
                .build();
            writer.write_alignment_record(&header, &read).unwrap();
        }
        writer.try_finish().unwrap();
    }

    fn stage1_fixture(tract_starts: &[u32]) -> (TempDir, TempDir, PathBuf, PathBuf, PathBuf) {
        use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};
        use crate::ssr::catalog::io::CatalogWriter;
        use crate::ssr::types::Motif;

        let (fa_dir, fasta) = build_fasta(&[ContigSpec {
            name: "chr1".into(),
            length: 200,
        }])
        .unwrap();
        let dir = TempDir::new().unwrap();
        let bam = dir.path().join("s.bam");
        write_bam_for_loci(&bam, tract_starts);

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

    fn read_records(path: &std::path::Path) -> Vec<ContainerRecord> {
        use crate::psp::PspReader;
        let mut reader = PspReader::new(File::open(path).unwrap()).unwrap();
        assert_eq!(reader.header().kind, "ssr");
        reader.records_of::<SsrKind>().map(|r| r.unwrap()).collect()
    }

    fn run_config(
        bam: &std::path::Path,
        fasta: &std::path::Path,
        catalog: &std::path::Path,
        output: &std::path::Path,
        threads: usize,
    ) -> SsrPileupConfig {
        use crate::ssr::pileup::fetch_reads::MAX_READS_PER_LOCUS;
        SsrPileupConfig {
            alignment_files: vec![bam.to_path_buf()],
            reference: fasta.to_path_buf(),
            catalog: catalog.to_path_buf(),
            output: output.to_path_buf(),
            filter: SegmentReadFilter::default(),
            cap: MAX_READS_PER_LOCUS,
            build_index_if_missing: true,
            sample: Some("S1".into()),
            threads,
        }
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
}
