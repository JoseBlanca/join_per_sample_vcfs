//! The Stage-1 `ssr-pileup` driver — turns the sorted catalog + one sample's
//! alignment files into a per-locus `.ssr.psp` evidence file.
//!
//! The work for **one** locus is a single self-contained unit, [`process_locus`]
//! — fetch the locus's reads, analyze each, fold them into the locus record. It
//! is pure over its shared `&` inputs (and [`AlignmentFile`] is `Sync`), so the
//! eventual parallel step is just `for` → `par_iter().map_init(…)`: each worker
//! reuses its own [`LocusScratch`], nothing inside the unit changes (arch §8.4).
//!
//! **Build status:** Stage 1 is runnable end to end — the per-locus unit
//! ([`process_locus`], [`LocusScratch`], [`qc_counts`]), the name→chrom_id
//! container adapter ([`to_container_record`]), the writer-header build, and the
//! single-threaded [`run`] loop, wired to the `ssr-pileup` CLI subcommand and
//! covered by a catalog→reference→BAM→`.ssr.psp` round-trip test. The
//! worker-pool parallelization (turn the `run` loop's `for` into a
//! `par_iter().map_init`) is the deferred §8.4 step.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use crate::bam::alignment_input::{PileupInputs, build_fasta_repository, load_pileup_inputs};
use crate::bam::errors::AlignmentInputError;
use crate::bam::segment_reader::{AlignmentFile, SegmentReadFilter};
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

use super::candidate_generation::CandidateAllele;
use super::fetch_reads::{LocusReads, fetch_locus_reads};
use super::locus_record::{QcCounts, SsrLocusRecord, aggregate};
use super::pair_hmm::{HmmModel, PairHmmScratch};
use super::read_analysis::analyze_read;

/// Default `analyze_read` candidate half-width (rungs): the pair-HMM scores
/// `observed_count ± DEFAULT_WINDOW` on-ladder lengths per spanning read. A
/// **calibration** placeholder (arch §14), like `MAX_READS_PER_LOCUS` /
/// `MIN_FLANK_BP` — wide enough to bracket genuine stutter/length variation
/// around the content pre-probe's estimate without inflating per-read work.
pub(crate) const DEFAULT_WINDOW: u16 = 10;

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
}

/// Per-locus scratch reused across loci to avoid allocation churn — the
/// candidate-allele buffer and the pair-HMM forward-matrix workspace
/// [`analyze_read`] writes through. **One per worker thread:** the serial driver
/// keeps a single instance; the parallel driver hands each thread its own via
/// `rayon`'s `map_init` (so the reuse survives parallelization).
pub(crate) struct LocusScratch {
    cands: Vec<CandidateAllele>,
    hmm: PairHmmScratch,
}

impl LocusScratch {
    pub(crate) fn new() -> Self {
        Self {
            cands: Vec::new(),
            hmm: PairHmmScratch::new(),
        }
    }
}

/// Assemble the locus's QC scalars from the fetch pass. All three describe
/// **independent primary alignments** at the locus:
///
/// - `depth` — usable primary reads considered (passed the reader's cheap
///   filter and overlapped the window; uncapped, pre reach-gate/reservoir).
/// - `n_filtered` — primary reads the gate dropped: quality (QC-fail / low-MAPQ
///   / too-short) **plus** duplicates (auditable as filtered).
/// - `mapped_reads` — the dup-free primary coverage denominator (`depth` + the
///   quality drops, *excluding* duplicates).
///
/// Secondary / supplementary (non-independent) and unmapped (not at the locus)
/// reads are excluded everywhere. So `mapped_reads ≥ depth`, and `n_filtered`
/// adds duplicates on top — it is deliberately *not* `mapped_reads − depth`.
/// (Decided in `ssr_pileup_driver.md` §4/§7.3.)
pub(crate) fn qc_counts(fetched: &LocusReads) -> QcCounts {
    let d = &fetched.filtered;
    let primary_quality_drops = d.qc_fail + d.low_mapq + d.too_short;
    QcCounts {
        depth: fetched.yielded as u32,
        n_filtered: (primary_quality_drops + d.duplicate) as u32,
        mapped_reads: (fetched.yielded + primary_quality_drops) as u32,
    }
}

/// Process one locus end to end: fetch + depth-cap its reads, analyze each
/// (triage → realign spanning reads → `Qᵣ`), and fold the outcomes + QC into the
/// in-memory (chrom-**name**-keyed) [`SsrLocusRecord`]. The name→chrom_id
/// container adapter is applied by the writer stage (next increment).
///
/// `window` is the `analyze_read` candidate half-width; `cap` is the per-locus
/// reservoir depth cap. Concurrency-safe for different loci: the only `&mut` is
/// the caller-owned `scratch`.
pub(crate) fn process_locus(
    files: &[AlignmentFile],
    locus: &Locus,
    window: u16,
    cap: usize,
    model: &HmmModel,
    scratch: &mut LocusScratch,
) -> Result<SsrLocusRecord, SsrPileupError> {
    let fetched = fetch_locus_reads(files, locus, cap)?;
    let outcomes: Vec<_> = fetched
        .reads
        .iter()
        .map(|read| {
            analyze_read(
                read,
                locus,
                window,
                model,
                &mut scratch.cands,
                &mut scratch.hmm,
            )
        })
        .collect();
    Ok(aggregate(locus, &outcomes, qc_counts(&fetched)))
}

/// Adapt the in-memory (chrom-**name**-keyed, 0-based) [`SsrLocusRecord`] to the
/// container's (chrom-**id**-keyed, 1-based) record — the SNP-symmetric name→id
/// boundary that keeps the two `SsrLocusRecord` types distinct by design.
///
/// **Coordinate shift:** the in-memory record carries the catalog locus's
/// 0-based half-open `[start, end)`; the container is 1-based half-open (its
/// `start` is "1-based", and an end-of-contig locus stores `end == length + 1`),
/// so both bounds shift by `+1`. The container derives `n_spanning` from
/// `spanning.len()`, so that field is dropped here.
pub(crate) fn to_container_record(
    record: SsrLocusRecord,
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
        n_flanking: record.n_flanking,
        n_frr: record.n_frr,
        n_filtered: record.n_filtered,
        mapped_reads: record.mapped_reads,
        spanning: record.spanning,
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
    /// `analyze_read` candidate half-width (rungs).
    pub window: u16,
    /// Per-locus reservoir depth cap.
    pub cap: usize,
    /// Build a missing alignment index in place rather than erroring.
    pub build_index_if_missing: bool,
    /// Sample name override; defaults to the one the inputs cross-validate.
    pub sample: Option<String>,
}

/// Build the `.ssr.psp` writer header — one [`ChromosomeEntry`] per reference
/// contig (mirrors the SNP `build_writer_header`), `subcommand = "ssr-pileup"`,
/// and the run parameters (including the catalog binding: its `reference_md5`
/// and `flank_bp`).
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
        "window".to_string(),
        ParameterValue::Integer(i64::from(cfg.window)),
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
/// alignment files, and write the per-locus evidence to a `.ssr.psp`.
///
/// Single-threaded: the per-locus body is [`process_locus`] (the parallelization
/// seam — see the module docs). Output is written to a temp file and atomically
/// renamed, mirroring the SNP `pileup` path.
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

    let model = HmmModel::default();
    let mut scratch = LocusScratch::new();

    while let Some(locus) = catalog.read_locus() {
        let locus = locus?;
        let record = process_locus(&files, &locus, cfg.window, cfg.cap, &model, &mut scratch)?;
        writer.write_locus(&to_container_record(record, &name_to_id)?)?;
    }

    // Atomic finish: flush + fsync the temp file, then rename into place.
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
        // 10 usable; quality drops 2+1+1=4; 3 dups; plus non-independent /
        // not-at-locus reads that must be ignored everywhere.
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

        assert_eq!(qc.depth, 10); // yielded
        assert_eq!(qc.n_filtered, 4 + 3); // quality drops + duplicates
        assert_eq!(qc.mapped_reads, 10 + 4); // depth + quality drops (dup-free)
        // Secondary / supplementary / unmapped influenced nothing.
    }

    #[test]
    fn qc_counts_coverage_is_at_least_depth_and_dup_free() {
        // Only duplicates filtered → coverage == depth (dups don't count), but
        // n_filtered still reports them.
        let filtered = FilterCounts {
            duplicate: 7,
            ..FilterCounts::default()
        };
        let qc = qc_counts(&locus_reads(20, filtered));

        assert_eq!(qc.depth, 20);
        assert_eq!(qc.mapped_reads, 20); // dups excluded from coverage
        assert_eq!(qc.n_filtered, 7); // but visible as filtered
    }

    #[test]
    fn qc_counts_all_zero_for_an_empty_locus() {
        let qc = qc_counts(&locus_reads(0, FilterCounts::default()));
        assert_eq!((qc.depth, qc.n_filtered, qc.mapped_reads), (0, 0, 0));
    }

    fn in_memory_record(chrom: &str, start: u32, end: u32) -> SsrLocusRecord {
        SsrLocusRecord {
            chrom: chrom.into(),
            start,
            end,
            depth: 5,
            n_spanning: 2,
            n_flanking: 1,
            n_frr: 0,
            n_filtered: 3,
            mapped_reads: 8,
            spanning: vec![vec![(10u16, -0.5f32)]],
        }
    }

    #[test]
    fn adapter_maps_name_to_id_and_shifts_coords_to_1_based() {
        let name_to_id = HashMap::from([("chr1", 0u32), ("chr2", 1)]);
        let c = to_container_record(in_memory_record("chr2", 16, 22), &name_to_id).unwrap();

        assert_eq!(c.chrom_id, 1);
        assert_eq!((c.start, c.end), (17, 23)); // 0-based [16,22) -> 1-based [17,23)
        // Other QC/profile fields pass through; n_spanning is derived from
        // spanning.len() by the container, so it is not carried here.
        assert_eq!(
            (c.depth, c.n_flanking, c.n_frr, c.n_filtered, c.mapped_reads),
            (5, 1, 0, 3, 8)
        );
        assert_eq!(c.spanning, vec![vec![(10u16, -0.5f32)]]);
    }

    #[test]
    fn adapter_errors_on_a_contig_absent_from_the_reference() {
        let name_to_id = HashMap::from([("chr1", 0u32)]);
        let err = to_container_record(in_memory_record("chrX", 1, 5), &name_to_id).unwrap_err();
        assert!(matches!(
            err,
            SsrPileupError::ContigNotInReference { name } if name == "chrX"
        ));
    }

    // --- header build -------------------------------------------------

    use crate::fasta::ContigEntry;
    use crate::ssr::catalog::CatalogParams;

    fn catalog_header() -> CatalogHeader {
        CatalogHeader {
            tool_version: "0.1.0".into(),
            reference: "ref.fa".into(),
            reference_md5: "a".repeat(32),
            trf_mod_version: "1.0".into(),
            params: CatalogParams::default(),
            date: "2026-06-16".into(),
        }
    }

    fn config_for(files: Vec<PathBuf>) -> SsrPileupConfig {
        SsrPileupConfig {
            alignment_files: files,
            reference: PathBuf::from("/data/ref.fa"),
            catalog: PathBuf::from("/data/cat.ssr.catalog"),
            output: PathBuf::from("/data/out.ssr.psp"),
            filter: SegmentReadFilter::default(),
            window: 8,
            cap: 1000,
            build_index_if_missing: false,
            sample: None,
        }
    }

    #[test]
    fn header_build_maps_contigs_records_params_and_subcommand() {
        let contigs = ContigList {
            entries: vec![
                ContigEntry {
                    name: "chr1".into(),
                    length: 1000,
                    md5: Some([1u8; 16]),
                },
                ContigEntry {
                    name: "chr2".into(),
                    length: 2000,
                    md5: Some([2u8; 16]),
                },
            ],
        };
        let cfg = config_for(vec![
            PathBuf::from("/data/s1.bam"),
            PathBuf::from("/data/s2.cram"),
        ]);
        let h = build_ssr_writer_header("SAMPLE", &cfg, &contigs, &catalog_header()).unwrap();

        assert_eq!(h.sample, "SAMPLE");
        assert_eq!(h.writer.subcommand, "ssr-pileup");
        assert_eq!(
            h.writer.input_crams,
            vec!["s1.bam".to_string(), "s2.cram".to_string()]
        );
        assert_eq!(h.writer.input_fasta, "ref.fa");
        assert_eq!(h.chromosomes.len(), 2);
        assert_eq!(h.chromosomes[1].name, "chr2");
        assert_eq!(h.chromosomes[1].length, 2000);
        assert_eq!(h.chromosomes[0].md5, format_md5_hex([1u8; 16]));
        // Params, incl. the catalog binding (flank_bp + reference_md5).
        assert!(matches!(
            h.writer.parameters.get("window"),
            Some(ParameterValue::Integer(8))
        ));
        assert!(matches!(
            h.writer.parameters.get("flank_bp"),
            Some(ParameterValue::Integer(50))
        ));
        assert!(matches!(
            h.writer.parameters.get("catalog_reference_md5"),
            Some(ParameterValue::String(s)) if *s == "a".repeat(32)
        ));
    }

    #[test]
    fn header_build_errors_when_a_contig_lacks_md5() {
        let contigs = ContigList {
            entries: vec![ContigEntry {
                name: "chr1".into(),
                length: 1000,
                md5: None,
            }],
        };
        let cfg = config_for(vec![PathBuf::from("/data/s1.bam")]);
        let err = build_ssr_writer_header("S", &cfg, &contigs, &catalog_header()).unwrap_err();
        assert!(matches!(err, SsrPileupError::MissingMd5 { name } if name == "chr1"));
    }

    // --- end-to-end: catalog + reference + BAM -> .ssr.psp ------------

    /// A coordinate-sorted BAM (one contig, `@RG SM`, a trusted `@SQ M5`) with a
    /// single clean spanning read for a CA(3) locus. `run()` builds the index
    /// (`build_index_if_missing`), so no `.csi` is written here.
    fn write_spanning_bam(path: &std::path::Path) {
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
        use std::num::NonZero;

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

        // Aligned ref [5,35) (pos 6, M30). The locus window is ref [10,28),
        // which maps to read[5..23] — laid out as the locus's embedded ref:
        // 6 bp G-flank + CACACA + 6 bp T-flank, padded to 30 bp.
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AAAAA");
        seq.extend_from_slice(b"GGGGGGCACACATTTTTT");
        seq.extend_from_slice(b"AAAAAAA");
        assert_eq!(seq.len(), 30);
        let read = RecordBuf::builder()
            .set_name(b"r1")
            .set_reference_sequence_id(0)
            .set_flags(Flags::default())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_alignment_start(noodles_core::Position::try_from(6).unwrap())
            .set_cigar([Op::new(Kind::Match, 30)].into_iter().collect())
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(vec![40u8; 30]))
            .build();

        let mut writer = bam::io::Writer::new(File::create(path).unwrap());
        writer.write_header(&header).unwrap();
        writer.write_alignment_record(&header, &read).unwrap();
        writer.try_finish().unwrap();
    }

    #[test]
    fn run_genotypes_one_locus_end_to_end() {
        use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};
        use crate::psp::PspReader;
        use crate::psp::registry_ssr::SsrKind;
        use crate::ssr::catalog::io::CatalogWriter;
        use crate::ssr::pileup::fetch_reads::MAX_READS_PER_LOCUS;
        use crate::ssr::types::{Locus, Motif};
        use tempfile::TempDir;

        let dir = TempDir::new().unwrap();
        let (_fa_dir, fasta) = build_fasta(&[ContigSpec {
            name: "chr1".into(),
            length: 200,
        }])
        .unwrap();

        // BAM with one spanning read.
        let bam = dir.path().join("s.bam");
        write_spanning_bam(&bam);

        // One-locus catalog: CA tract ref [16,22), 6 bp flanks embedded.
        let locus = Locus::new(
            "chr1".into(),
            16,
            22,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"GGGGGGCACACATTTTTT").into(),
            10,
        )
        .unwrap();
        let catalog = dir.path().join("c.ssr.catalog");
        {
            let mut w =
                CatalogWriter::new(File::create(&catalog).unwrap(), &catalog_header()).unwrap();
            w.write_locus(&locus).unwrap();
            w.finish().unwrap();
        }

        let output = dir.path().join("out.ssr.psp");
        let cfg = SsrPileupConfig {
            alignment_files: vec![bam],
            reference: fasta,
            catalog,
            output: output.clone(),
            filter: SegmentReadFilter::default(),
            window: DEFAULT_WINDOW,
            cap: MAX_READS_PER_LOCUS,
            build_index_if_missing: true,
            sample: Some("S1".into()),
        };
        run(&cfg).expect("run");

        let mut reader = PspReader::new(File::open(&output).unwrap()).unwrap();
        assert_eq!(reader.header().kind, "ssr");
        let records: Vec<_> = reader.records_of::<SsrKind>().map(|r| r.unwrap()).collect();

        assert_eq!(records.len(), 1);
        let rec = &records[0];
        assert_eq!(rec.chrom_id, 0);
        assert_eq!((rec.start, rec.end), (17, 23)); // 0-based [16,22) -> 1-based [17,23)
        assert_eq!(rec.depth, 1);
        assert_eq!(rec.mapped_reads, 1);
        assert_eq!(rec.n_filtered, 0);
        assert_eq!(rec.spanning.len(), 1); // one spanning read -> one Qᵣ profile
        // The profile's most-likely length is the true allele: 3 CA copies.
        let best = rec.spanning[0]
            .iter()
            .copied()
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap();
        assert_eq!(best.0, 3);
    }
}
