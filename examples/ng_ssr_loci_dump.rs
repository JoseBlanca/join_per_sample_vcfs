//! **The STR locus generator's dump tool** — `locus_generation_ssr.md` Milestone E, the
//! acceptance anchor. Runs the real pipeline over a reference + a sample's reads and prints, per
//! microsatellite tract, what that sample's reads showed: the observed tract sequences and their
//! support. This generator emits no calls, so "done" must not decay into "compiles" — the output
//! is inspectable and, on a committed fixture, asserted (spec §9).
//!
//! ```text
//! ng_ssr_loci_dump <reference.fa> <sample.bam|cram> [contig]
//! ```
//!
//! The reference's sibling `<reference.fa>.fai` is used (created if absent); reads go through ng's
//! real ingestion (`SampleReads`, the filtered reads a caller sees). An optional third argument
//! restricts the walk to one contig by name. The pipeline is: region typing over the reference →
//! one `SsrSegment` per tract → the `SsrGenerator` turns each into one locus → its observations
//! become rows.
//!
//! Output: a `#`-prefixed `key=value` counts header, a bare TSV column line, then one tab-separated
//! row per observed sequence. `depth` sums the reads behind **complete** observations only (the ones
//! that pinned the tract); partial reads appear on their own rows, tagged `partial:left` /
//! `partial:right`, because conflating a lower bound with an exact length is the mistake spec §3
//! exists to prevent.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::ng::locus_generation::ssr::{SsrGenerator, SsrGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    LocusGenerator, LocusKind, ReadCoverage, SampleLocusObservations,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::segment_criteria::SsrSegment;
use pop_var_caller::ng::region_typing::{RegionKind, TypedRegionConfig, TypedRegionIterator};
use pop_var_caller::ng::types::{Bp, ContigId};

/// One TSV row: an observed tract sequence at a locus, with its support.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ObservationRow {
    contig: String,
    /// 1-based tract coordinates (inclusive).
    start: u64,
    end: u64,
    motif: Vec<u8>,
    ref_tract: Vec<u8>,
    /// The locus's complete-observation depth (sum of the complete rows' read counts) — shown on
    /// every row of the locus, so a partial row is read against the depth that actually pinned it.
    depth: u32,
    read_coverage: &'static str,
    observed: Vec<u8>,
    reads: u32,
}

/// The whole dump: the run-level counts plus the rows. Separated from rendering so the tests can
/// assert on the numbers directly, not by parsing text.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct DumpReport {
    /// Loci emitted — one per `SsrSegment`, including uncovered ones.
    ssr_loci: u64,
    /// Loci no read reached at all (present-but-empty).
    zero_coverage: u64,
    /// Reads the depth cap discarded, across all loci.
    reads_capped: u64,
    /// Reads that reached the aligner and yielded nothing (by any reason).
    reads_without_observation: u64,
    /// Complete / partial observation totals (reads that became each kind).
    obs_complete: u64,
    obs_partial: u64,
    /// Reads fetched over all loci — **per fetch event, not per distinct read** (a read
    /// overlapping two loci' query spans is fetched, and counted, at each). This is the left side
    /// of the accounting identity (`= obs_complete + obs_partial + reads_capped +
    /// reads_without_observation`), which holds because every fetched read folds into exactly one
    /// of those run-level counters.
    reads_fetched: u64,
    /// The typed-region walk's own `ssr_loci` — `ssr_loci` must equal it (one locus per segment).
    walk_ssr_loci: u64,
    rows: Vec<ObservationRow>,
}

impl DumpReport {
    /// Fold one locus into the report: count it, note zero-coverage, and emit a row per observation.
    fn push_locus(&mut self, locus: &SampleLocusObservations, segment: &SsrSegment) {
        self.ssr_loci += 1;
        // Zero coverage is "no read reached the tract" — distinct from "reads reached it and said
        // nothing" (that is `reads_without_observation`) and from "the cap dropped them all".
        if locus.observed_sequences.is_empty()
            && locus.reads_without_observation == 0
            && locus.reads_discarded_by_cap == 0
        {
            self.zero_coverage += 1;
        }
        let depth: u32 = locus.complete_observations().map(|obs| obs.num_obs).sum();
        let motif = match &locus.kind {
            LocusKind::Ssr(detail) => detail.motif.as_bytes().to_vec(),
            _ => Vec::new(),
        };
        for obs in &locus.observed_sequences {
            self.rows.push(ObservationRow {
                contig: segment.chrom().to_string(),
                start: locus.region.start.get(),
                end: locus.region.end.get(),
                motif: motif.clone(),
                ref_tract: locus.reference_bases.to_vec(),
                depth,
                read_coverage: coverage_label(obs.read_coverage),
                observed: obs.bases.to_vec(),
                reads: obs.num_obs,
            });
        }
    }

    /// The dump text: the two `#` header lines, the TSV column line, then the rows.
    fn render(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        // Writing into a `String` via `fmt::Write` never fails, so the `write!` results are ignored.
        let _ = writeln!(
            out,
            "# ssr_loci={} zero_coverage={} reads_capped={} reads_without_observation={}",
            self.ssr_loci, self.zero_coverage, self.reads_capped, self.reads_without_observation
        );
        let _ = writeln!(
            out,
            "# obs_complete={} obs_partial={}",
            self.obs_complete, self.obs_partial
        );
        out.push_str(
            "contig\tstart\tend\tmotif\tref_tract\tdepth\tread_coverage\tobserved\treads\n",
        );
        for row in &self.rows {
            let _ = writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                row.contig,
                row.start,
                row.end,
                String::from_utf8_lossy(&row.motif),
                String::from_utf8_lossy(&row.ref_tract),
                row.depth,
                row.read_coverage,
                String::from_utf8_lossy(&row.observed),
                row.reads,
            );
        }
        out
    }
}

/// The tag a read-coverage carries in the `read_coverage` column.
fn coverage_label(coverage: ReadCoverage) -> &'static str {
    match coverage {
        ReadCoverage::Complete => "complete",
        ReadCoverage::PartialLeft(_) => "partial:left",
        ReadCoverage::PartialRight(_) => "partial:right",
    }
}

/// Run the whole pipeline over `fasta` + `bams`, optionally restricted to `contig_filter` (a contig
/// name), building the [`DumpReport`]. `gen_config` is the STR generator's config (the tests vary
/// its cap). The reference's `.fai` is created if absent; the BAM index likewise (`SampleReads`).
fn run_dump(
    fasta: &Path,
    bams: &[PathBuf],
    contig_filter: Option<&str>,
    gen_config: SsrGeneratorConfig,
) -> Result<DumpReport, Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(&cache, fasta.to_path_buf())?;
    let contigs = info.contig_list();

    let sample = SampleReads::open(bams, &info, ReadFilterConfig::default(), true)?;

    let walk_config = TypedRegionConfig::default();
    // The generator holds its own reference (margin fetch) and a factory (the per-file read query),
    // both windowed over the same FASTA — the reference seam the STR generator's doc calls the Arc
    // gap. Both are cheap: a path plus the contig table, nothing resident until a fetch.
    let mut generator = SsrGenerator::new(
        WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone()),
        {
            let fasta = fasta.to_path_buf();
            let contigs = contigs.clone();
            move || WindowedRefSeq::new(fasta.clone(), contigs.clone())
        },
        gen_config,
        Bp(walk_config.criteria.bundle_threshold),
    )?;

    let mut report = DumpReport::default();
    for (index, entry) in contigs.entries.iter().enumerate() {
        if contig_filter.is_some_and(|name| entry.name != name) {
            continue;
        }
        // A fresh windowed reference per contig for the walk (it takes the reference by value).
        let walk_reference = WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone());
        let mut walk = TypedRegionIterator::over_contig(
            walk_reference,
            ContigId(index as u32),
            walk_config.clone(),
        )?;
        for region in walk.by_ref() {
            let region = region?;
            if let RegionKind::SsrSegment(segment) = &region.kind {
                generator.begin_segment(region.region);
                while let Some(locus) = generator.next_locus(segment, &sample)? {
                    report.push_locus(&locus, segment);
                }
            }
        }
        report.walk_ssr_loci += walk.counts().ssr_loci;
    }

    // The background `.fai` verification (only when a pre-existing `.fai` was used) is joined at the
    // end, so a stale index is a failure rather than a silently-wrong walk.
    if let Some(handle) = verify {
        handle.join()?;
    }

    let counts = generator.counts();
    report.reads_fetched = counts.reads_fetched;
    report.reads_capped = counts.reads_discarded_by_cap;
    report.obs_complete = counts.observations_complete;
    report.obs_partial = counts.observations_partial;
    report.reads_without_observation =
        counts.no_border_anchored + counts.low_quality + counts.window_truncated;
    Ok(report)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: ng_ssr_loci_dump <reference.fa> <sample.bam|cram> [contig]\n\
             dumps, per microsatellite tract, the observed tract sequences one sample's reads showed."
        );
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&args[1]);
    let bam = PathBuf::from(&args[2]);
    let contig_filter = args.get(3).map(String::as_str);

    match run_dump(&fasta, &[bam], contig_filter, SsrGeneratorConfig::default()) {
        Ok(report) => {
            print!("{}", report.render());
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use tempfile::TempDir;

    // chr1 carries TWO clean (AC)×12 tracts, 40 bp apart: LEFT + tract1 + MID + tract2 + RIGHT.
    // Each tract (period 2, 12 copies, purity 1.0) is past the copy-number floor; every flank
    // (LEFT/MID/RIGHT) exceeds the 30 bp bundle threshold, and the tracts are 40 bp apart (>30), so
    // neither is bundled and both sit ≥30 bp from the contig ends — region typing emits both as
    // clean `SsrSegment`s. The reads only reach tract1 (all end ≤ position 74, clear of tract2's
    // query span 75..158), so tract2 is a genuine **zero-coverage** locus. The flanks are aperiodic
    // and their joins with the AC tracts ("…CTT|AC…", "…AC|CAG…", "…AC|GG…") do not extend the run.
    const LEFT: &[u8] = b"GATCTTGCAAGCTGGAATCCGTTACGATCGGATCAAGCTT";
    const TRACT: &[u8] = b"ACACACACACACACACACACACAC"; // (AC)×12 = 24 bp; tract1 = 41..64, tract2 = 105..128
    const MID: &[u8] = b"CAGTTGCACGATCCTAAGGCTTGACCATGGATCCAAGTTG";
    const RIGHT: &[u8] = b"GGTTCAAGATCCGGATCTTGCAATCGGATCAAGCTTGACT";

    fn contig() -> Vec<u8> {
        assert_eq!(LEFT.len(), 40);
        assert_eq!(TRACT.len(), 24);
        assert_eq!(MID.len(), 40);
        assert_eq!(RIGHT.len(), 40);
        [LEFT, TRACT, MID, TRACT, RIGHT].concat()
    }

    /// A reference-allele read at 1-based `start` spanning `len` bases — its sequence is the contig
    /// slice it maps to, so it carries zero mismatches and clears the mismatch filter.
    fn read(
        contig: &[u8],
        qname: &str,
        start: usize,
        len: usize,
    ) -> noodles_sam::alignment::RecordBuf {
        use noodles_core::Position;
        use noodles_sam::alignment::record::cigar::Op;
        use noodles_sam::alignment::record::cigar::op::Kind;
        use noodles_sam::alignment::record::{Flags, MappingQuality};
        use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
        let seq = contig[start - 1..start - 1 + len].to_vec();
        noodles_sam::alignment::RecordBuf::builder()
            .set_name(qname.as_bytes())
            .set_reference_sequence_id(0)
            .set_flags(Flags::empty())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar([Op::new(Kind::Match, len)].into_iter().collect())
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(vec![40u8; len]))
            .build()
    }

    /// Write `chr1` as a one-line FASTA (its `.fai` is created by `run_dump`).
    fn write_fasta(path: &Path, contig: &[u8]) {
        use std::io::Write as _;
        let mut file = File::create(path).unwrap();
        writeln!(file, ">chr1").unwrap();
        file.write_all(contig).unwrap();
        writeln!(file).unwrap();
    }

    /// Write a coordinate-sorted single-contig BAM (`@RG SM`, no `@SQ M5`) holding `reads`.
    fn write_bam(path: &Path, contig_len: usize, reads: &[noodles_sam::alignment::RecordBuf]) {
        use bstr::BString;
        use noodles_bam as bam;
        use noodles_sam as sam;
        use sam::alignment::io::Write as _;
        use sam::header::record::value::Map;
        use sam::header::record::value::map::header::Version;
        use sam::header::record::value::map::header::tag::SORT_ORDER;
        use sam::header::record::value::map::read_group::tag::SAMPLE;
        use sam::header::record::value::map::{Header as HeaderMap, ReadGroup, ReferenceSequence};
        use std::num::NonZero;

        let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
        hd.other_fields_mut()
            .insert(SORT_ORDER, BString::from("coordinate"));
        let sq = Map::<ReferenceSequence>::new(NonZero::new(contig_len).unwrap());
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
        for record in reads {
            writer.write_alignment_record(&header, record).unwrap();
        }
        writer.try_finish().unwrap();
    }

    /// The committed fixture: the reference and a BAM whose reads produce, at **tract1**, four
    /// complete observations, one left-partial and one right-partial — while **tract2** is left
    /// uncovered (zero coverage). Every read ends at or before position 74, clear of tract2's query
    /// span (75..158). Reads are in coordinate order.
    fn fixture() -> (TempDir, PathBuf, PathBuf) {
        let contig = contig();
        let dir = TempDir::new().unwrap();
        let fasta = dir.path().join("ref.fa");
        let bam = dir.path().join("sample.bam");
        write_fasta(&fasta, &contig);

        let reads = vec![
            // Four identical reads spanning tract1 + both its flanks (21..74) → complete "AC"×12.
            read(&contig, "c0", 21, 54),
            read(&contig, "c1", 21, 54),
            read(&contig, "c2", 21, 54),
            read(&contig, "c3", 21, 54),
            // Anchors tract1's left flank, runs off inside the tract (21..60) → partial:left.
            read(&contig, "pl", 21, 40),
            // Begins inside tract1, anchors its right flank (45..74) → partial:right.
            read(&contig, "pr", 45, 30),
        ];
        write_bam(&bam, contig.len(), &reads);
        (dir, fasta, bam)
    }

    fn dump(fasta: &Path, bam: &Path, config: SsrGeneratorConfig) -> DumpReport {
        run_dump(
            fasta,
            std::slice::from_ref(&bam.to_path_buf()),
            None,
            config,
        )
        .expect("the fixture dumps")
    }

    /// Both AC tracts are emitted as loci (one per `SsrSegment`, including the uncovered tract2),
    /// and every fetched read is accounted for — a complete observation, a partial one, a cap
    /// discard, or a no-observation (spec §9.1, §9.2).
    #[test]
    fn every_segment_is_one_locus_and_every_read_is_accounted() {
        let (_dir, fasta, bam) = fixture();
        let report = dump(&fasta, &bam, SsrGeneratorConfig::default());

        assert_eq!(report.walk_ssr_loci, 2, "both AC tracts are detected");
        assert_eq!(
            report.ssr_loci, report.walk_ssr_loci,
            "one emitted locus per typed SsrSegment (spec §9.1)"
        );
        assert_eq!(
            report.zero_coverage, 1,
            "tract2 is present but uncovered (spec §9.1) — 'looked and saw nothing' ≠ 'never looked'"
        );
        assert_eq!(
            report.reads_fetched,
            report.obs_complete
                + report.obs_partial
                + report.reads_capped
                + report.reads_without_observation,
            "every fetched read is accounted for (spec §9.2)"
        );
        // tract1's six reads land as four complete + two partial observations; tract2 gets none.
        assert_eq!(report.obs_complete, 4);
        assert_eq!(report.obs_partial, 2);
        assert_eq!(report.reads_without_observation, 0);
        assert_eq!(report.reads_capped, 0);
    }

    /// The rendered text is exactly the spec §9 shape: two `#` header lines, the TSV column line,
    /// then a tab-separated row per observation. This pins the **format** — column order, the tabs,
    /// the header keys — which the structured assertions (reading `report.rows`) cannot catch.
    #[test]
    fn render_emits_the_spec_9_header_and_tsv_rows() {
        let (_dir, fasta, bam) = fixture();
        let report = dump(&fasta, &bam, SsrGeneratorConfig::default());
        let text = report.render();
        let lines: Vec<&str> = text.lines().collect();

        assert_eq!(
            lines[0],
            "# ssr_loci=2 zero_coverage=1 reads_capped=0 reads_without_observation=0"
        );
        assert_eq!(lines[1], "# obs_complete=4 obs_partial=2");
        assert_eq!(
            lines[2],
            "contig\tstart\tend\tmotif\tref_tract\tdepth\tread_coverage\tobserved\treads"
        );
        // The complete row, built from the report's own coordinates so the assertion pins the TSV
        // format (tabs, column order, values) without hard-coding region typing's tract bounds.
        let complete = report
            .rows
            .iter()
            .find(|row| row.read_coverage == "complete")
            .expect("a complete row");
        let expected = format!(
            "chr1\t{}\t{}\tAC\tACACACACACACACACACACACAC\t4\tcomplete\tACACACACACACACACACACACAC\t4",
            complete.start, complete.end
        );
        assert!(
            text.lines().any(|line| line == expected),
            "the complete row should render as:\n{expected}\ngot:\n{text}"
        );
    }

    /// Partial observations exist — which proves the relevance gate admitted the partially-covering
    /// reads (a spanning-only gate would have dropped them, spec §9.4). One left, one right, on
    /// their own rows, tagged and distinct from the complete rows.
    #[test]
    fn partial_observations_are_present_and_tagged() {
        let (_dir, fasta, bam) = fixture();
        let report = dump(&fasta, &bam, SsrGeneratorConfig::default());

        assert!(
            report.obs_partial >= 1,
            "the relevance gate admitted partials"
        );
        let left = report
            .rows
            .iter()
            .find(|row| row.read_coverage == "partial:left")
            .expect("a left partial row");
        let right = report
            .rows
            .iter()
            .find(|row| row.read_coverage == "partial:right")
            .expect("a right partial row");
        // A partial is a lower bound — shorter than the complete tract it sits under.
        assert!(left.observed.len() < left.ref_tract.len());
        assert!(right.observed.len() < right.ref_tract.len());
        // Depth on a partial row is the locus's *complete* depth (the four spanning reads).
        assert_eq!(left.depth, 4);
        assert_eq!(right.depth, 4);
    }

    /// The complete row is the reference tract, at depth 4, and the row for it carries the tract as
    /// both `ref_tract` and `observed`.
    #[test]
    fn the_complete_observation_is_the_reference_tract() {
        let (_dir, fasta, bam) = fixture();
        let report = dump(&fasta, &bam, SsrGeneratorConfig::default());
        let complete = report
            .rows
            .iter()
            .find(|row| row.read_coverage == "complete")
            .expect("a complete row");
        assert_eq!(complete.observed, TRACT);
        assert_eq!(complete.ref_tract, TRACT);
        assert_eq!(complete.motif, b"AC");
        assert_eq!(complete.reads, 4);
        assert_eq!(complete.depth, 4);
    }

    /// The output is byte-identical across repeated runs (spec §9.5).
    #[test]
    fn output_is_deterministic_across_runs() {
        let (_dir, fasta, bam) = fixture();
        let first = dump(&fasta, &bam, SsrGeneratorConfig::default()).render();
        let second = dump(&fasta, &bam, SsrGeneratorConfig::default()).render();
        assert_eq!(first, second);
    }

    /// Raising the cap above the deepest locus (6 reads) leaves the output unchanged; a cap *below*
    /// it changes the output — that is what a cap does, not a determinism failure (spec §9.5).
    #[test]
    fn a_cap_above_the_depth_is_invisible_and_a_cap_below_it_bites() {
        let (_dir, fasta, bam) = fixture();
        let uncapped = dump(
            &fasta,
            &bam,
            SsrGeneratorConfig {
                flank_bp: Bp(30),
                max_reads_per_locus: None,
            },
        )
        .render();
        let cap_above = dump(
            &fasta,
            &bam,
            SsrGeneratorConfig {
                flank_bp: Bp(30),
                max_reads_per_locus: Some(100),
            },
        )
        .render();
        assert_eq!(uncapped, cap_above, "a cap above the depth changes nothing");

        let cap_below = dump(
            &fasta,
            &bam,
            SsrGeneratorConfig {
                flank_bp: Bp(30),
                max_reads_per_locus: Some(2),
            },
        );
        assert_ne!(
            uncapped,
            cap_below.render(),
            "a cap below the depth changes the output"
        );
        assert!(cap_below.reads_capped > 0, "the cap discarded reads");
    }
}
