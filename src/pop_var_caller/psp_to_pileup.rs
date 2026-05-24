//! `pop_var_caller psp-to-pileup` — stream a `.psp` artefact as a
//! samtools-mpileup-style text dump with an extra trailing column for
//! PSP-specific per-allele aggregates.
//!
//! Read-side surface: re-opens the file with [`PspReader`] and walks
//! either every record (`records()`) or a coordinate-clamped window
//! (`region_records()`). Each record is encoded into a 7-column
//! tab-delimited line and streamed to the sink.
//!
//! Plan: `ia/feature_implementation_plans/psp_to_pileup.md`.

use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Args;
use thiserror::Error;

use crate::per_sample_pileup::psp::{PspReadError, PspReader, RecordsIter};
use crate::pileup_record::{AlleleObservation, PileupRecord};
use crate::pop_var_caller::common::DEFAULT_BUFFERED_IO_CAPACITY;

// ---------------------------------------------------------------------
// Clap surface
// ---------------------------------------------------------------------

/// Arguments accepted by the `psp-to-pileup` subcommand.
#[derive(Debug, Args, Clone)]
pub struct PspToPileupArgs {
    /// Input .psp file.
    #[arg(long)]
    pub input: PathBuf,

    /// Output path. Use `-` (the default) to stream to stdout.
    #[arg(long, default_value = "-")]
    pub output: PathBuf,

    /// Restrict output to a region: `chrom` or `chrom:start-end`
    /// (1-based, inclusive on both ends). When omitted the whole
    /// file is streamed in genomic order.
    #[arg(long)]
    pub region: Option<String>,

    /// Append the per-allele `chain_ids` list (semicolon-separated)
    /// to each allele record in the trailing ALLELE_DETAILS column.
    #[arg(long)]
    pub show_chain_ids: bool,
}

// ---------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------

#[derive(Debug, Error)]
pub enum PspToPileupError {
    #[error("PSP read: {0}")]
    Read(#[from] PspReadError),
    #[error("io: {0}")]
    Io(#[from] io::Error),
    #[error("region: {0}")]
    Region(String),
    #[error("chromosome '{0}' not in input .psp")]
    UnknownChromosome(String),
}

// ---------------------------------------------------------------------
// Region spec
// ---------------------------------------------------------------------

/// A parsed `--region` value. `start` / `end` are 1-based and
/// inclusive on both ends, with `None` meaning "from the beginning"
/// or "to the end" respectively.
#[derive(Debug, Clone, PartialEq, Eq)]
struct RegionSpec {
    chrom: String,
    start: Option<u32>,
    end: Option<u32>,
}

fn parse_region(s: &str) -> Result<RegionSpec, PspToPileupError> {
    if s.is_empty() {
        return Err(PspToPileupError::Region("empty region".into()));
    }
    // Split on the first colon. Anything after is the position range;
    // the chromosome name may itself contain no colons by spec.
    let (chrom, range) = match s.split_once(':') {
        None => {
            return Ok(RegionSpec {
                chrom: s.into(),
                start: None,
                end: None,
            });
        }
        Some((c, r)) => (c, r),
    };
    if chrom.is_empty() {
        return Err(PspToPileupError::Region(format!(
            "missing chromosome in '{s}'"
        )));
    }
    // `range` is either "<start>" or "<start>-<end>".
    let (start_str, end_str) = match range.split_once('-') {
        None => (range, None),
        Some((a, b)) => (a, Some(b)),
    };
    let start: u32 = start_str.parse().map_err(|_| {
        PspToPileupError::Region(format!("invalid start position '{start_str}' in '{s}'"))
    })?;
    if start == 0 {
        return Err(PspToPileupError::Region(format!(
            "start must be 1-based (got 0) in '{s}'"
        )));
    }
    let end: Option<u32> = match end_str {
        None => None,
        Some(e) => {
            let v: u32 = e.parse().map_err(|_| {
                PspToPileupError::Region(format!("invalid end position '{e}' in '{s}'"))
            })?;
            if v < start {
                return Err(PspToPileupError::Region(format!(
                    "end ({v}) is before start ({start}) in '{s}'"
                )));
            }
            Some(v)
        }
    };
    Ok(RegionSpec {
        chrom: chrom.into(),
        start: Some(start),
        end,
    })
}

// ---------------------------------------------------------------------
// Top-level driver
// ---------------------------------------------------------------------

/// Open the `.psp`, build the output sink, and stream every record
/// (or every record inside `--region`) as a 7-column tab-delimited
/// line.
pub fn run_psp_to_pileup(args: &PspToPileupArgs) -> Result<(), PspToPileupError> {
    let file = File::open(&args.input)?;
    let mut reader = PspReader::new(BufReader::with_capacity(DEFAULT_BUFFERED_IO_CAPACITY, file))?;

    // Snapshot the chrom_id → name map up-front. The reader's parsed
    // header sits behind an immutable borrow, which would conflict
    // with the iterator we are about to mint; cloning avoids both
    // sides fighting over the same &mut.
    let chrom_names: Vec<String> = reader
        .header()
        .chromosomes
        .iter()
        .map(|c| c.name.clone())
        .collect();

    let region: Option<(u32, u32, u32)> = match &args.region {
        None => None,
        Some(raw) => {
            let spec = parse_region(raw)?;
            let chrom_id = chrom_names
                .iter()
                .position(|n| n == &spec.chrom)
                .ok_or_else(|| PspToPileupError::UnknownChromosome(spec.chrom.clone()))?
                as u32;
            Some((
                chrom_id,
                spec.start.unwrap_or(1),
                spec.end.unwrap_or(u32::MAX),
            ))
        }
    };

    let mut sink: Box<dyn Write> = open_sink(&args.output)?;
    let mut line = String::with_capacity(256);
    let show_chains = args.show_chain_ids;

    // Local helper to drive a `RecordsIter` to exhaustion. The
    // closure form keeps the duplication between the
    // whole-file and region paths to one place — both branches
    // hand the iterator to the same drain loop.
    let drain = |it: RecordsIter<'_, _>,
                 sink: &mut Box<dyn Write>,
                 line: &mut String|
     -> Result<(), PspToPileupError> {
        for record in it {
            let record = record?;
            line.clear();
            let chrom_name = chrom_name_for(&chrom_names, record.chrom_id);
            emit_line(line, &record, chrom_name, show_chains);
            line.push('\n');
            sink.write_all(line.as_bytes())?;
        }
        Ok(())
    };

    match region {
        None => drain(reader.records(), &mut sink, &mut line)?,
        Some((chrom_id, start, end)) => drain(
            reader.region_records(chrom_id, start, end),
            &mut sink,
            &mut line,
        )?,
    }

    sink.flush()?;
    Ok(())
}

fn open_sink(path: &PathBuf) -> Result<Box<dyn Write>, io::Error> {
    if path.as_os_str() == "-" {
        Ok(Box::new(io::stdout().lock()))
    } else {
        let file = File::create(path)?;
        Ok(Box::new(BufWriter::with_capacity(
            DEFAULT_BUFFERED_IO_CAPACITY,
            file,
        )))
    }
}

fn chrom_name_for(chrom_names: &[String], chrom_id: u32) -> &str {
    chrom_names
        .get(chrom_id as usize)
        .map(String::as_str)
        // The PSP writer validates `chrom_id < chromosomes.len()` at
        // record-emit time, so a record carrying an out-of-range
        // chrom_id means the file violated its own invariant. Fall
        // back to a literal that is obviously wrong rather than
        // silently misnaming, so any downstream eyeball catches it.
        .unwrap_or("<unknown>")
}

// ---------------------------------------------------------------------
// Line encoder
// ---------------------------------------------------------------------

/// Format a single `PileupRecord` into 7 tab-separated columns and
/// append the result to `out` (without a trailing newline). The
/// caller controls newline placement and buffer reuse.
fn emit_line(out: &mut String, record: &PileupRecord, chrom_name: &str, show_chains: bool) {
    let ref_seq = &record.alleles[0].seq;
    let ref_span = ref_seq.len();
    let ref_anchor = ref_seq[0]; // alleles[0].seq is non-empty by walker invariant
    let depth: u32 = record.alleles.iter().map(|a| a.support.num_obs).sum();

    // Columns 1..=3
    out.push_str(chrom_name);
    out.push('\t');
    let _ = write!(out, "{}", record.pos);
    out.push('\t');
    out.push(ref_anchor as char);
    out.push('\t');

    // Column 4: depth
    let _ = write!(out, "{}", depth);
    out.push('\t');

    // Column 5: read bases — emit per-allele runs in PSP order.
    for (idx, allele) in record.alleles.iter().enumerate() {
        emit_bases_for_allele(out, idx, allele, ref_anchor, ref_span, ref_seq);
    }
    out.push('\t');

    // Column 6: quality placeholder. `!` (Phred 0) × depth.
    for _ in 0..depth {
        out.push('!');
    }
    out.push('\t');

    // Column 7: per-allele details, comma-separated.
    let mut first = true;
    for allele in &record.alleles {
        if !first {
            out.push(',');
        }
        first = false;
        emit_allele_detail(out, allele, show_chains);
    }
}

fn emit_bases_for_allele(
    out: &mut String,
    allele_idx: usize,
    allele: &AlleleObservation,
    ref_anchor: u8,
    ref_span: usize,
    ref_seq: &[u8],
) {
    let num_obs = allele.support.num_obs;
    let fwd = allele.support.fwd.min(num_obs);
    let rev = num_obs - fwd;
    let anchor_byte = *allele.seq.first().unwrap_or(&ref_anchor);
    let anchor_matches = allele_idx == 0 || anchor_byte == ref_anchor;

    match allele.seq.len().cmp(&ref_span) {
        std::cmp::Ordering::Equal => {
            // SNP, MNP, or REF — emit one character per supporting read.
            if anchor_matches {
                for _ in 0..fwd {
                    out.push('.');
                }
                for _ in 0..rev {
                    out.push(',');
                }
            } else {
                let up = anchor_byte.to_ascii_uppercase() as char;
                let lo = anchor_byte.to_ascii_lowercase() as char;
                for _ in 0..fwd {
                    out.push(up);
                }
                for _ in 0..rev {
                    out.push(lo);
                }
            }
        }
        std::cmp::Ordering::Greater => {
            // Insertion: anchor + `+N<inserted>` per read.
            let inserted = &allele.seq[ref_span..];
            let fwd_marker = format_insertion_marker(inserted, true);
            let rev_marker = format_insertion_marker(inserted, false);
            emit_indel_runs(
                out,
                anchor_byte,
                anchor_matches,
                fwd,
                rev,
                &fwd_marker,
                &rev_marker,
            );
        }
        std::cmp::Ordering::Less => {
            // Deletion: anchor + `-N<deleted-ref-bases>` per read.
            // For single-event deletions at offset 0 the deleted bases
            // are exactly the suffix of REF that allele.seq stops
            // short of. Multi-event alleles aren't represented
            // unambiguously by seq alone; that limitation is
            // documented in the plan.
            let deleted = &ref_seq[allele.seq.len()..];
            let fwd_marker = format_deletion_marker(deleted, true);
            let rev_marker = format_deletion_marker(deleted, false);
            emit_indel_runs(
                out,
                anchor_byte,
                anchor_matches,
                fwd,
                rev,
                &fwd_marker,
                &rev_marker,
            );
        }
    }
}

fn emit_indel_runs(
    out: &mut String,
    anchor_byte: u8,
    anchor_matches: bool,
    fwd: u32,
    rev: u32,
    fwd_marker: &str,
    rev_marker: &str,
) {
    let (fwd_anchor, rev_anchor) = if anchor_matches {
        ('.', ',')
    } else {
        (
            anchor_byte.to_ascii_uppercase() as char,
            anchor_byte.to_ascii_lowercase() as char,
        )
    };
    for _ in 0..fwd {
        out.push(fwd_anchor);
        out.push_str(fwd_marker);
    }
    for _ in 0..rev {
        out.push(rev_anchor);
        out.push_str(rev_marker);
    }
}

fn format_insertion_marker(inserted: &[u8], upper: bool) -> String {
    let mut s = String::with_capacity(2 + inserted.len() + 4);
    s.push('+');
    let _ = write!(s, "{}", inserted.len());
    for &b in inserted {
        s.push(if upper {
            b.to_ascii_uppercase() as char
        } else {
            b.to_ascii_lowercase() as char
        });
    }
    s
}

fn format_deletion_marker(deleted: &[u8], upper: bool) -> String {
    let mut s = String::with_capacity(2 + deleted.len() + 4);
    s.push('-');
    let _ = write!(s, "{}", deleted.len());
    for &b in deleted {
        s.push(if upper {
            b.to_ascii_uppercase() as char
        } else {
            b.to_ascii_lowercase() as char
        });
    }
    s
}

fn emit_allele_detail(out: &mut String, allele: &AlleleObservation, show_chains: bool) {
    // seq field: render bytes as UTF-8 if possible (always {A,C,G,T,N}
    // by spec), or a `*` if the allele is empty (defensive — the
    // walker shouldn't produce empty alleles today).
    if allele.seq.is_empty() {
        out.push('*');
    } else {
        for &b in &allele.seq {
            out.push(b as char);
        }
    }
    let _ = write!(
        out,
        ":{}:{}:{}:{}:{:.3}",
        allele.support.num_obs,
        allele.support.fwd,
        allele.support.placed_left,
        allele.support.placed_start,
        allele.support.q_sum,
    );
    if show_chains {
        out.push(':');
        let mut first = true;
        for chain_id in &allele.chain_ids {
            if !first {
                out.push(';');
            }
            first = false;
            let _ = write!(out, "{}", chain_id);
        }
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};

    fn supp(num_obs: u32, fwd: u32) -> AlleleSupportStats {
        AlleleSupportStats::new(num_obs, -1.234, fwd, 0, 0, 0, 0)
    }

    fn allele(seq: &[u8], support: AlleleSupportStats) -> AlleleObservation {
        AlleleObservation::new(seq.to_vec(), support, Vec::new())
    }

    fn allele_with_chains(
        seq: &[u8],
        support: AlleleSupportStats,
        chains: Vec<u64>,
    ) -> AlleleObservation {
        AlleleObservation::new(seq.to_vec(), support, chains)
    }

    // ---------- parse_region ----------

    #[test]
    fn parse_region_chrom_only() {
        let r = parse_region("chr1").unwrap();
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, None);
        assert_eq!(r.end, None);
    }

    #[test]
    fn parse_region_chrom_start() {
        let r = parse_region("chr1:1000").unwrap();
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, Some(1000));
        assert_eq!(r.end, None);
    }

    #[test]
    fn parse_region_chrom_start_end() {
        let r = parse_region("chr1:1000-2000").unwrap();
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, Some(1000));
        assert_eq!(r.end, Some(2000));
    }

    #[test]
    fn parse_region_rejects_empty() {
        assert!(parse_region("").is_err());
    }

    #[test]
    fn parse_region_rejects_missing_chrom() {
        assert!(parse_region(":1000").is_err());
    }

    #[test]
    fn parse_region_rejects_nonnumeric_start() {
        assert!(parse_region("chr1:abc").is_err());
    }

    #[test]
    fn parse_region_rejects_nonnumeric_end() {
        assert!(parse_region("chr1:1000-xyz").is_err());
    }

    #[test]
    fn parse_region_rejects_end_before_start() {
        assert!(parse_region("chr1:2000-1000").is_err());
    }

    #[test]
    fn parse_region_rejects_zero_start() {
        assert!(parse_region("chr1:0-1000").is_err());
    }

    // ---------- emit_line: REF-only positions ----------

    fn ref_only_record(pos: u32, ref_base: u8, num_obs: u32, fwd: u32) -> PileupRecord {
        PileupRecord::new(0, pos, vec![allele(&[ref_base], supp(num_obs, fwd))])
    }

    #[test]
    fn emit_line_ref_only_all_forward() {
        let mut out = String::new();
        let r = ref_only_record(100, b'A', 5, 5);
        emit_line(&mut out, &r, "chr1", false);
        let cols: Vec<&str> = out.split('\t').collect();
        assert_eq!(cols[0], "chr1");
        assert_eq!(cols[1], "100");
        assert_eq!(cols[2], "A");
        assert_eq!(cols[3], "5");
        assert_eq!(cols[4], ".....");
        assert_eq!(cols[5], "!!!!!");
        assert_eq!(cols[6], "A:5:5:0:0:-1.234");
    }

    #[test]
    fn emit_line_ref_only_mixed_strand() {
        let mut out = String::new();
        let r = ref_only_record(100, b'C', 6, 4);
        emit_line(&mut out, &r, "chr1", false);
        let cols: Vec<&str> = out.split('\t').collect();
        assert_eq!(cols[3], "6");
        assert_eq!(cols[4], "....,,");
        assert_eq!(cols[5], "!!!!!!");
    }

    // ---------- emit_line: SNP ----------

    #[test]
    fn emit_line_snp_alt_mixed_strand() {
        let mut out = String::new();
        let record = PileupRecord::new(
            0,
            200,
            vec![
                allele(b"A", supp(3, 2)), // REF, fwd=2 rev=1
                allele(b"G", supp(4, 1)), // ALT G, fwd=1 rev=3
            ],
        );
        emit_line(&mut out, &record, "chr1", false);
        let cols: Vec<&str> = out.split('\t').collect();
        assert_eq!(cols[2], "A");
        assert_eq!(cols[3], "7");
        // REF block: `..,` (2 fwd, 1 rev). ALT block: `Gggg` (1 fwd, 3 rev).
        assert_eq!(cols[4], "..,Gggg");
        assert_eq!(cols[5], "!!!!!!!");
        let allele_field = cols[6];
        assert!(allele_field.starts_with("A:3:2:0:0:"));
        assert!(allele_field.contains(",G:4:1:0:0:"));
    }

    // ---------- emit_line: insertion ----------

    #[test]
    fn emit_line_insertion_emits_plus_marker_per_read() {
        let mut out = String::new();
        // REF "A" at pos 500; insertion-bearing allele "ACG" (anchor
        // A + inserted CG). 2 forward + 1 reverse supports.
        let record = PileupRecord::new(
            0,
            500,
            vec![
                allele(b"A", supp(0, 0)),   // REF, no supports
                allele(b"ACG", supp(3, 2)), // INS allele
            ],
        );
        emit_line(&mut out, &record, "chr1", false);
        let cols: Vec<&str> = out.split('\t').collect();
        assert_eq!(cols[3], "3");
        assert_eq!(cols[4], ".+2CG.+2CG,+2cg");
        assert_eq!(cols[5], "!!!");
    }

    // ---------- emit_line: deletion ----------

    #[test]
    fn emit_line_deletion_emits_minus_marker_per_read() {
        let mut out = String::new();
        // REF span 3 ("ACG"); DEL allele "A" (anchor preserved,
        // CG deleted). 1 fwd + 1 rev.
        let record = PileupRecord::new(
            0,
            500,
            vec![
                allele(b"ACG", supp(0, 0)), // REF (multi-base span)
                allele(b"A", supp(2, 1)),   // DEL
            ],
        );
        emit_line(&mut out, &record, "chr1", false);
        let cols: Vec<&str> = out.split('\t').collect();
        assert_eq!(cols[2], "A");
        assert_eq!(cols[3], "2");
        assert_eq!(cols[4], ".-2CG,-2cg");
        assert_eq!(cols[5], "!!");
    }

    // ---------- emit_line: chain ids toggle ----------

    #[test]
    fn emit_line_chain_ids_field_off_by_default() {
        let mut out = String::new();
        let record = PileupRecord::new(
            0,
            10,
            vec![allele_with_chains(b"A", supp(2, 1), vec![7, 9])],
        );
        emit_line(&mut out, &record, "chr1", false);
        let cols: Vec<&str> = out.split('\t').collect();
        assert_eq!(cols[6], "A:2:1:0:0:-1.234");
    }

    #[test]
    fn emit_line_chain_ids_field_on_when_flag_set() {
        let mut out = String::new();
        let record = PileupRecord::new(
            0,
            10,
            vec![allele_with_chains(b"A", supp(2, 1), vec![7, 9])],
        );
        emit_line(&mut out, &record, "chr1", true);
        let cols: Vec<&str> = out.split('\t').collect();
        assert_eq!(cols[6], "A:2:1:0:0:-1.234:7;9");
    }

    #[test]
    fn emit_line_chain_ids_empty_list_when_no_chains() {
        let mut out = String::new();
        let record = PileupRecord::new(0, 10, vec![allele_with_chains(b"A", supp(2, 1), vec![])]);
        emit_line(&mut out, &record, "chr1", true);
        let cols: Vec<&str> = out.split('\t').collect();
        // Chains field is present but empty (trailing colon).
        assert_eq!(cols[6], "A:2:1:0:0:-1.234:");
    }

    // ---------- placed_left / placed_start surface through ----------

    #[test]
    fn allele_detail_includes_placed_counters() {
        let mut out = String::new();
        let support = AlleleSupportStats::new(5, -2.0, 3, 1, 2, 0, 0);
        let record = PileupRecord::new(0, 1, vec![allele(b"A", support)]);
        emit_line(&mut out, &record, "chr1", false);
        let cols: Vec<&str> = out.split('\t').collect();
        // num_obs:fwd:placed_left:placed_start:q_sum
        assert_eq!(cols[6], "A:5:3:1:2:-2.000");
    }
}
