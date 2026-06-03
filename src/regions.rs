//! Genomic analysis regions.
//!
//! The pipeline always operates on a [`RegionSet`] — a sorted,
//! non-overlapping list of spans to analyze. When the user supplies a
//! BED file, that is the set ([`RegionSet::from_bed_reader`]); when
//! they do not, the set is one full-length span per contig
//! ([`RegionSet::whole_contigs`]). "Whole genome" is not a special
//! case — it is the region set whose every span covers an entire
//! contig. Both stages (`pileup`, `var-calling`) build a `RegionSet`
//! and seek to its spans rather than walking the genome.
//!
//! **Coordinate convention.** Spans are **1-based inclusive**
//! `[start, end]`, matching the PSP reader's region API
//! ([`crate::psp::PspReader::region_records`]) and the pileup walker's
//! position convention. A BED file is 0-based half-open `[start, end)`;
//! the parser converts a BED span `[b_start, b_end)` to the 1-based
//! inclusive span `[b_start + 1, b_end]` at the single boundary where
//! BED text is read.
//!
//! This module is a top-level peer consumed by both pipeline stages,
//! so it deliberately does not depend on any stage's contig type
//! (`fasta::ContigList`, `psp::header::ChromosomeEntry`). Callers
//! adapt their own contig list into a slice of the neutral
//! [`ContigBounds`].

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

use thiserror::Error;

/// One analysis span on a single contig, in **1-based inclusive**
/// coordinates. Invariants held by every `Region` a [`RegionSet`]
/// hands out: `1 <= start <= end`, and `chrom_id` indexes the contig
/// slice the set was built from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Region {
    /// Zero-based index of the contig in the [`ContigBounds`] slice
    /// the set was built from — the `chrom_id` the PSP reader and
    /// pileup walker use.
    pub chrom_id: u32,
    /// 1-based inclusive lower bound.
    pub start: u32,
    /// 1-based inclusive upper bound (`>= start`).
    pub end: u32,
}

/// A reference contig's identity for region resolution: the name
/// matched against a BED `chrom` column, and the length spans are
/// clamped to. The position in the slice handed to a [`RegionSet`]
/// constructor is the resulting `chrom_id`, so callers must pass
/// contigs in the same order the downstream reader uses (the
/// reference/`.fai` order for pileup, the PSP header order for
/// var-calling).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ContigBounds<'a> {
    /// Contig name, matched verbatim against the BED `chrom` column.
    pub name: &'a str,
    /// Contig length in bases. A span's inclusive end is clamped to
    /// this; a span that starts past it is rejected.
    pub length: u32,
}

/// A sorted, non-overlapping set of [`Region`]s to analyze.
///
/// The internal `Vec<Region>` is ordered by `(chrom_id, start)` and
/// carries no two spans that overlap or touch — directly-adjacent
/// spans (a zero-base gap) are coalesced, since they select a
/// contiguous run of bases and merging them changes which bases are
/// selected not at all, only the span count. Iterating the set yields
/// spans in genomic order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionSet {
    regions: Vec<Region>,
}

impl RegionSet {
    /// One full-length span per contig: `[1, length]` for each contig
    /// in `contigs`, in slice order. Contigs of length 0 contribute no
    /// span (there is no 1-based position to cover). This is the
    /// default set used when the user supplies no BED file.
    pub fn whole_contigs(contigs: &[ContigBounds]) -> RegionSet {
        let regions = contigs
            .iter()
            .enumerate()
            .filter(|(_, c)| c.length > 0)
            .map(|(idx, c)| Region {
                chrom_id: idx as u32,
                start: 1,
                end: c.length,
            })
            .collect();
        // Already sorted by chrom_id and non-overlapping by
        // construction (one span per contig); no normalization needed.
        RegionSet { regions }
    }

    /// Parse a BED file from `reader` and resolve it against `contigs`.
    ///
    /// Each data line contributes a span; the result is sorted by
    /// `(chrom_id, start)` with overlapping and adjacent spans merged.
    /// Blank lines, `#` comments, and UCSC `track`/`browser` header
    /// lines are ignored. Columns past the first three (start/end) are
    /// ignored.
    ///
    /// # Errors
    ///
    /// See [`BedError`]: a line with fewer than three columns or
    /// non-numeric coordinates ([`BedError::Parse`]), `end <= start`
    /// ([`BedError::InvalidInterval`]), a `chrom` not in `contigs`
    /// ([`BedError::UnknownContig`]), or a span whose first base lies
    /// past the contig's end ([`BedError::IntervalBeyondContig`]). A
    /// span that merely overhangs the contig end is clamped to the
    /// contig length rather than rejected.
    pub fn from_bed_reader<R: BufRead>(
        reader: R,
        contigs: &[ContigBounds],
    ) -> Result<RegionSet, BedError> {
        let by_name: HashMap<&str, (u32, u32)> = contigs
            .iter()
            .enumerate()
            .map(|(idx, c)| (c.name, (idx as u32, c.length)))
            .collect();

        let mut regions: Vec<Region> = Vec::new();
        for (line_idx, line) in reader.lines().enumerate() {
            let line_number = line_idx + 1;
            let line = line.map_err(BedError::Io)?;
            if let Some(region) = parse_bed_line(&line, line_number, &by_name)? {
                regions.push(region);
            }
        }

        sort_and_merge(&mut regions);
        Ok(RegionSet { regions })
    }

    /// Open `path` and parse it as a BED file via
    /// [`from_bed_reader`](Self::from_bed_reader). The opened path is
    /// attached to any I/O error.
    pub fn from_bed_path(path: &Path, contigs: &[ContigBounds]) -> Result<RegionSet, BedError> {
        let file = File::open(path).map_err(|source| BedError::Open {
            path: path.to_path_buf(),
            source,
        })?;
        Self::from_bed_reader(BufReader::new(file), contigs)
    }

    /// The spans, in genomic order (`(chrom_id, start)`).
    pub fn regions(&self) -> &[Region] {
        &self.regions
    }

    /// Iterate the spans in genomic order.
    pub fn iter(&self) -> std::slice::Iter<'_, Region> {
        self.regions.iter()
    }

    /// Number of spans.
    pub fn len(&self) -> usize {
        self.regions.len()
    }

    /// Whether the set selects no spans at all.
    pub fn is_empty(&self) -> bool {
        self.regions.is_empty()
    }
}

/// Parse one BED line into a [`Region`], or `None` for lines that
/// carry no span (blank, comment, or UCSC header). `by_name` resolves
/// the `chrom` column to its `(chrom_id, length)`.
fn parse_bed_line(
    line: &str,
    line_number: usize,
    by_name: &HashMap<&str, (u32, u32)>,
) -> Result<Option<Region>, BedError> {
    let trimmed = line.trim();
    if trimmed.is_empty()
        || trimmed.starts_with('#')
        || trimmed.starts_with("track")
        || trimmed.starts_with("browser")
    {
        return Ok(None);
    }

    let mut fields = trimmed.split_whitespace();
    let chrom = fields.next().ok_or_else(|| BedError::Parse {
        line_number,
        reason: "empty line after trimming".to_string(),
    })?;
    let (start_raw, end_raw) = match (fields.next(), fields.next()) {
        (Some(s), Some(e)) => (s, e),
        _ => {
            return Err(BedError::Parse {
                line_number,
                reason: "expected at least three columns (chrom, start, end)".to_string(),
            });
        }
    };

    let bed_start: u64 = start_raw.parse().map_err(|_| BedError::Parse {
        line_number,
        reason: format!("start column `{start_raw}` is not a non-negative integer"),
    })?;
    let bed_end: u64 = end_raw.parse().map_err(|_| BedError::Parse {
        line_number,
        reason: format!("end column `{end_raw}` is not a non-negative integer"),
    })?;

    // BED is 0-based half-open; an empty or inverted span carries no
    // bases and is almost certainly a mistake.
    if bed_end <= bed_start {
        return Err(BedError::InvalidInterval {
            line_number,
            start: bed_start,
            end: bed_end,
        });
    }

    let (chrom_id, length) =
        by_name
            .get(chrom)
            .copied()
            .ok_or_else(|| BedError::UnknownContig {
                line_number,
                name: chrom.to_string(),
            })?;

    // Convert to 1-based inclusive: [b_start, b_end) → [b_start+1, b_end].
    let start_1based = bed_start + 1;
    if start_1based > length as u64 {
        return Err(BedError::IntervalBeyondContig {
            line_number,
            name: chrom.to_string(),
            start: bed_start,
            contig_length: length,
        });
    }
    // A span that overhangs the contig end is clamped rather than
    // rejected (bedtools-style leniency for slightly-long BEDs).
    let end_1based = bed_end.min(length as u64);

    Ok(Some(Region {
        chrom_id,
        start: start_1based as u32,
        end: end_1based as u32,
    }))
}

/// Sort `regions` by `(chrom_id, start)` and merge spans that overlap
/// or directly abut. In-place; leaves a sorted, non-overlapping list.
fn sort_and_merge(regions: &mut Vec<Region>) {
    regions.sort_unstable_by_key(|r| (r.chrom_id, r.start, r.end));

    let mut merged: Vec<Region> = Vec::with_capacity(regions.len());
    for region in regions.drain(..) {
        match merged.last_mut() {
            // Same contig and the next span starts no later than one
            // base past the current end → they describe a contiguous
            // run; extend. `cur.end + 1` is computed in u64 so an end
            // at u32::MAX cannot overflow.
            Some(cur)
                if cur.chrom_id == region.chrom_id
                    && (region.start as u64) <= cur.end as u64 + 1 =>
            {
                cur.end = cur.end.max(region.end);
            }
            _ => merged.push(region),
        }
    }
    *regions = merged;
}

/// Failure modes of BED parsing and resolution. Every variant that
/// originates from a line carries its 1-based `line_number` for a
/// pointable error message.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum BedError {
    /// The BED file could not be opened.
    #[error("could not open BED file {path}: {source}")]
    Open { path: PathBuf, source: io::Error },

    /// An I/O error occurred while reading a line.
    #[error("reading BED file: {0}")]
    Io(io::Error),

    /// A data line was structurally malformed: fewer than three
    /// columns, or non-numeric coordinates.
    #[error("BED line {line_number}: {reason}")]
    Parse { line_number: usize, reason: String },

    /// `end <= start` — a BED span must cover at least one base
    /// (`start < end`, half-open).
    #[error(
        "BED line {line_number}: empty or inverted interval \
         (start={start}, end={end}); BED requires start < end"
    )]
    InvalidInterval {
        line_number: usize,
        start: u64,
        end: u64,
    },

    /// The `chrom` column names a contig absent from the reference /
    /// PSP contig list.
    #[error("BED line {line_number}: contig `{name}` is not in the reference contig list")]
    UnknownContig { line_number: usize, name: String },

    /// The span's first base lies entirely past the contig's end —
    /// likely a BED built against a different reference.
    #[error(
        "BED line {line_number}: interval start {start} on contig `{name}` \
         is past the contig length {contig_length}"
    )]
    IntervalBeyondContig {
        line_number: usize,
        name: String,
        start: u64,
        contig_length: u32,
    },
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn contigs() -> Vec<ContigBounds<'static>> {
        vec![
            ContigBounds {
                name: "chr1",
                length: 1000,
            },
            ContigBounds {
                name: "chr2",
                length: 500,
            },
        ]
    }

    fn parse(bed: &str) -> Result<RegionSet, BedError> {
        RegionSet::from_bed_reader(bed.as_bytes(), &contigs())
    }

    // --- whole_contigs --------------------------------------------------

    #[test]
    fn whole_contigs_spans_each_contig_one_based_inclusive() {
        let set = RegionSet::whole_contigs(&contigs());
        assert_eq!(
            set.regions(),
            &[
                Region {
                    chrom_id: 0,
                    start: 1,
                    end: 1000
                },
                Region {
                    chrom_id: 1,
                    start: 1,
                    end: 500
                },
            ]
        );
    }

    #[test]
    fn whole_contigs_skips_zero_length_contigs() {
        let cs = [
            ContigBounds {
                name: "empty",
                length: 0,
            },
            ContigBounds {
                name: "chr1",
                length: 10,
            },
        ];
        let set = RegionSet::whole_contigs(&cs);
        // chrom_id is the slice index, so the surviving span keeps id 1.
        assert_eq!(
            set.regions(),
            &[Region {
                chrom_id: 1,
                start: 1,
                end: 10
            }]
        );
    }

    // --- coordinate conversion -----------------------------------------

    #[test]
    fn bed_zero_based_half_open_becomes_one_based_inclusive() {
        // BED [0, 100) covers 0-based bases 0..=99 = 1-based 1..=100.
        let set = parse("chr1\t0\t100\n").unwrap();
        assert_eq!(
            set.regions(),
            &[Region {
                chrom_id: 0,
                start: 1,
                end: 100
            }]
        );
    }

    #[test]
    fn single_base_bed_span_maps_to_single_position() {
        // BED [4, 5) is just 0-based base 4 = 1-based position 5.
        let set = parse("chr1\t4\t5\n").unwrap();
        assert_eq!(
            set.regions(),
            &[Region {
                chrom_id: 0,
                start: 5,
                end: 5
            }]
        );
    }

    // --- sorting + merging ---------------------------------------------

    #[test]
    fn spans_are_sorted_by_contig_then_start() {
        let set = parse("chr2\t10\t20\nchr1\t50\t60\nchr1\t0\t10\n").unwrap();
        let ids_starts: Vec<(u32, u32)> = set.iter().map(|r| (r.chrom_id, r.start)).collect();
        assert_eq!(ids_starts, vec![(0, 1), (0, 51), (1, 11)]);
    }

    #[test]
    fn overlapping_spans_merge() {
        // [0,100) and [50,150) overlap → one span [1,150].
        let set = parse("chr1\t0\t100\nchr1\t50\t150\n").unwrap();
        assert_eq!(
            set.regions(),
            &[Region {
                chrom_id: 0,
                start: 1,
                end: 150
            }]
        );
    }

    #[test]
    fn adjacent_spans_merge() {
        // [0,100) → 1-based [1,100]; [100,200) → 1-based [101,200].
        // They abut (gap of zero bases) and coalesce to [1,200].
        let set = parse("chr1\t0\t100\nchr1\t100\t200\n").unwrap();
        assert_eq!(
            set.regions(),
            &[Region {
                chrom_id: 0,
                start: 1,
                end: 200
            }]
        );
    }

    #[test]
    fn spans_with_a_gap_stay_separate() {
        // [0,100) → [1,100]; [200,300) → [201,300]. One base gap at
        // position 101..=200 (well, a real gap) keeps them apart.
        let set = parse("chr1\t0\t100\nchr1\t200\t300\n").unwrap();
        assert_eq!(
            set.regions(),
            &[
                Region {
                    chrom_id: 0,
                    start: 1,
                    end: 100
                },
                Region {
                    chrom_id: 0,
                    start: 201,
                    end: 300
                },
            ]
        );
    }

    #[test]
    fn merge_does_not_cross_contig_boundaries() {
        let set = parse("chr1\t0\t100\nchr2\t0\t100\n").unwrap();
        assert_eq!(set.len(), 2);
        assert_eq!(set.regions()[0].chrom_id, 0);
        assert_eq!(set.regions()[1].chrom_id, 1);
    }

    // --- line skipping + extra columns ---------------------------------

    #[test]
    fn comment_blank_and_header_lines_are_ignored() {
        let bed = "# a comment\n\
                   browser position chr1\n\
                   track name=foo\n\
                   \n\
                   chr1\t0\t10\n";
        let set = parse(bed).unwrap();
        assert_eq!(set.len(), 1);
        assert_eq!(
            set.regions()[0],
            Region {
                chrom_id: 0,
                start: 1,
                end: 10
            }
        );
    }

    #[test]
    fn extra_columns_are_ignored() {
        let set = parse("chr1\t0\t10\tname\t0\t+\n").unwrap();
        assert_eq!(
            set.regions()[0],
            Region {
                chrom_id: 0,
                start: 1,
                end: 10
            }
        );
    }

    #[test]
    fn space_separated_columns_are_accepted() {
        let set = parse("chr1 0 10\n").unwrap();
        assert_eq!(
            set.regions()[0],
            Region {
                chrom_id: 0,
                start: 1,
                end: 10
            }
        );
    }

    // --- clamping + errors ---------------------------------------------

    #[test]
    fn span_overhanging_contig_end_is_clamped() {
        // chr2 has length 500; BED [400, 600) → 1-based [401, 600],
        // clamped to [401, 500].
        let set = parse("chr2\t400\t600\n").unwrap();
        assert_eq!(
            set.regions(),
            &[Region {
                chrom_id: 1,
                start: 401,
                end: 500
            }]
        );
    }

    #[test]
    fn span_exactly_at_contig_end_is_kept() {
        // chr2 length 500; BED [499, 500) → 1-based [500, 500].
        let set = parse("chr2\t499\t500\n").unwrap();
        assert_eq!(
            set.regions(),
            &[Region {
                chrom_id: 1,
                start: 500,
                end: 500
            }]
        );
    }

    #[test]
    fn span_starting_past_contig_end_errors() {
        // chr2 length 500; 0-based start 500 → 1-based 501 > 500.
        let err = parse("chr2\t500\t600\n").unwrap_err();
        assert!(matches!(
            err,
            BedError::IntervalBeyondContig { ref name, start: 500, contig_length: 500, .. } if name == "chr2"
        ));
    }

    #[test]
    fn unknown_contig_errors() {
        let err = parse("chrX\t0\t10\n").unwrap_err();
        assert!(
            matches!(err, BedError::UnknownContig { ref name, line_number: 1 } if name == "chrX")
        );
    }

    #[test]
    fn empty_interval_errors() {
        let err = parse("chr1\t10\t10\n").unwrap_err();
        assert!(matches!(
            err,
            BedError::InvalidInterval {
                start: 10,
                end: 10,
                line_number: 1
            }
        ));
    }

    #[test]
    fn inverted_interval_errors() {
        let err = parse("chr1\t20\t10\n").unwrap_err();
        assert!(matches!(
            err,
            BedError::InvalidInterval {
                start: 20,
                end: 10,
                ..
            }
        ));
    }

    #[test]
    fn too_few_columns_errors() {
        let err = parse("chr1\t10\n").unwrap_err();
        assert!(matches!(err, BedError::Parse { line_number: 1, .. }));
    }

    #[test]
    fn non_numeric_start_errors() {
        let err = parse("chr1\tfoo\t10\n").unwrap_err();
        assert!(matches!(err, BedError::Parse { line_number: 1, .. }));
    }

    #[test]
    fn line_numbers_point_at_the_offending_line() {
        // Two good lines, then a bad one on line 3.
        let err = parse("chr1\t0\t10\nchr1\t20\t30\nchrX\t0\t10\n").unwrap_err();
        assert!(matches!(
            err,
            BedError::UnknownContig { line_number: 3, .. }
        ));
    }

    #[test]
    fn empty_bed_yields_empty_region_set() {
        let set = parse("# only a comment\n").unwrap();
        assert!(set.is_empty());
        assert_eq!(set.len(), 0);
    }
}
