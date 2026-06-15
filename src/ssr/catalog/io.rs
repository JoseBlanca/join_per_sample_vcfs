//! The catalog wire format: a self-describing, bgzip-wrapped TSV
//! (architecture [`ssr_catalog.md`](../../../doc/devel/architecture/ssr_catalog.md) §7).
//!
//! Layout, in order:
//! - a `##`-prefixed metadata block (`## key: value` lines) — the reference +
//!   its md5, the `trf-mod` version, the post-process parameters, `flank_bp`,
//!   the tool version, and the build date;
//! - a single `#`-prefixed column-header line;
//! - one tab-separated data row per locus:
//!   `chrom  start  end  motif  purity_fraction  ref_seq_start  ref_seq`.
//!
//! Everything else (`period`, copy number, perfect/imperfect, motif class) is
//! *derived* from these columns, never stored. The whole stream is bgzip-framed
//! so Stages 1–2 can region-query it through a coordinate index; comment lines
//! (`##` / `#`) are skipped by tabix/CSI.
//!
//! This module owns serialisation in one place — [`CatalogWriter`] (write side,
//! Stage 0) and [`CatalogReader`] (read side, Stages 1–2) share the row and
//! header codecs, so the two can never drift. The coordinate index
//! (`write_index` / `CatalogReader::query`) is a later increment; the
//! sequential `read_locus` path is here now.

use std::io::{BufRead, BufReader, Read, Write};

use noodles_bgzf as bgzf;

use super::CatalogError;
use crate::ssr::types::{Locus, Motif};

/// The column-header line (without its trailing newline). Written verbatim and
/// validated on read, so an upstream column reshuffle fails loudly.
const COLUMN_HEADER: &str = "#chrom\tstart\tend\tmotif\tpurity_fraction\tref_seq_start\tref_seq";

/// The post-process / build parameters recorded in the catalog header, so a
/// reader can see exactly how the catalog was constructed (and Stage 1 can read
/// back `flank_bp`, which sizes its read-anchoring band).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CatalogParams {
    /// Flank margin (bp) embedded each side of the tract in `ref_seq`.
    pub flank_bp: u32,
    /// Purity floor applied after recomputation (a degeneracy cutoff).
    pub min_purity: f32,
    /// Early TRF-score accept-gate.
    pub min_score: i32,
    /// Bundle-drop radius (bp). `>= flank_bp` guarantees clean flanks.
    pub bundle_threshold: u32,
}

/// The catalog's `##` metadata block.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CatalogHeader {
    /// Building tool version (e.g. `"0.1.0"`).
    pub tool_version: String,
    /// Reference path, as supplied to the builder.
    pub reference: String,
    /// Reference md5 (32 lowercase hex chars; upper-cased-content convention).
    pub reference_md5: String,
    /// The `trf-mod` version string, pinned and recorded.
    pub trf_mod_version: String,
    /// Build parameters.
    pub params: CatalogParams,
    /// Build date (ISO `YYYY-MM-DD`), supplied by the caller — this module
    /// never reads the clock.
    pub date: String,
}

impl CatalogHeader {
    /// Serialise the `##` metadata block (including the trailing column-header
    /// line), each line newline-terminated. Key order is fixed for
    /// determinism.
    fn write_to(&self, out: &mut impl Write) -> std::io::Result<()> {
        writeln!(out, "## tool: ssr-catalog")?;
        writeln!(out, "## tool_version: {}", self.tool_version)?;
        writeln!(out, "## reference: {}", self.reference)?;
        writeln!(out, "## reference_md5: {}", self.reference_md5)?;
        writeln!(out, "## trf_mod_version: {}", self.trf_mod_version)?;
        writeln!(out, "## flank_bp: {}", self.params.flank_bp)?;
        writeln!(out, "## min_purity: {}", self.params.min_purity)?;
        writeln!(out, "## min_score: {}", self.params.min_score)?;
        writeln!(out, "## bundle_threshold: {}", self.params.bundle_threshold)?;
        writeln!(out, "## date: {}", self.date)?;
        writeln!(out, "{COLUMN_HEADER}")?;
        Ok(())
    }

    /// Parse the `##` metadata block + the column-header line from `lines`,
    /// consuming exactly through the column header. `*line_no` is advanced to
    /// the count of lines consumed.
    fn parse_from(
        reader: &mut impl BufRead,
        line_no: &mut usize,
    ) -> Result<CatalogHeader, CatalogError> {
        let mut meta: std::collections::BTreeMap<String, String> =
            std::collections::BTreeMap::new();
        let mut buf = String::new();
        loop {
            buf.clear();
            let n = reader
                .read_line(&mut buf)
                .map_err(|source| CatalogError::Io {
                    context: "read catalog header",
                    source,
                })?;
            if n == 0 {
                return Err(CatalogError::HeaderParse {
                    reason: "unexpected end of stream before the column header".to_string(),
                });
            }
            *line_no += 1;
            let line = buf.trim_end_matches(['\n', '\r']);
            if let Some(rest) = line.strip_prefix("## ") {
                let (key, value) =
                    rest.split_once(": ")
                        .ok_or_else(|| CatalogError::HeaderParse {
                            reason: format!("metadata line is not `## key: value`: {line:?}"),
                        })?;
                meta.insert(key.to_string(), value.to_string());
            } else if line == COLUMN_HEADER {
                break;
            } else {
                return Err(CatalogError::HeaderParse {
                    reason: format!(
                        "expected a `## ` metadata line or the column header, got {line:?}"
                    ),
                });
            }
        }

        let get = |key: &str| -> Result<String, CatalogError> {
            meta.get(key)
                .cloned()
                .ok_or_else(|| CatalogError::HeaderParse {
                    reason: format!("missing required header key {key:?}"),
                })
        };
        let parse_u32 = |key: &str, raw: &str| -> Result<u32, CatalogError> {
            raw.parse::<u32>().map_err(|_| CatalogError::HeaderParse {
                reason: format!("header key {key:?} = {raw:?} is not a u32"),
            })
        };

        let flank_raw = get("flank_bp")?;
        let purity_raw = get("min_purity")?;
        let score_raw = get("min_score")?;
        let bundle_raw = get("bundle_threshold")?;
        let params = CatalogParams {
            flank_bp: parse_u32("flank_bp", &flank_raw)?,
            min_purity: purity_raw
                .parse::<f32>()
                .map_err(|_| CatalogError::HeaderParse {
                    reason: format!("header key \"min_purity\" = {purity_raw:?} is not an f32"),
                })?,
            min_score: score_raw
                .parse::<i32>()
                .map_err(|_| CatalogError::HeaderParse {
                    reason: format!("header key \"min_score\" = {score_raw:?} is not an i32"),
                })?,
            bundle_threshold: parse_u32("bundle_threshold", &bundle_raw)?,
        };
        Ok(CatalogHeader {
            tool_version: get("tool_version")?,
            reference: get("reference")?,
            reference_md5: get("reference_md5")?,
            trf_mod_version: get("trf_mod_version")?,
            params,
            date: get("date")?,
        })
    }
}

/// Serialise one locus to a catalog data row (no trailing newline). The motif
/// and `ref_seq` bytes are ASCII (`A/C/G/T/N`), so they render as text. `f32`
/// purity is formatted with the shortest round-tripping representation.
fn locus_to_row(locus: &Locus) -> String {
    // Motif / ref_seq are validated ASCII; lossless as UTF-8. `motif()`
    // returns by value, so bind it before borrowing its bytes.
    let motif = locus.motif();
    let motif_str =
        std::str::from_utf8(motif.as_bytes()).expect("motif bytes are ASCII by construction");
    let ref_seq =
        std::str::from_utf8(locus.ref_bytes()).expect("ref_seq bytes are ASCII by construction");
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        locus.chrom(),
        locus.start(),
        locus.end(),
        motif_str,
        locus.purity_fraction(),
        locus.ref_bytes_start(),
        ref_seq,
    )
}

/// Parse one catalog data row back into a [`Locus`], validating column count,
/// numeric fields, and the motif / locus invariants. `line` is the 1-based
/// source line for diagnostics.
fn row_to_locus(row: &str, line: usize) -> Result<Locus, CatalogError> {
    let mut fields = row.split('\t');
    let mut next = |name: &str| -> Result<&str, CatalogError> {
        fields.next().ok_or_else(|| CatalogError::RowParse {
            line,
            reason: format!("missing column {name:?}"),
        })
    };
    let chrom = next("chrom")?.to_string();
    let start_raw = next("start")?;
    let end_raw = next("end")?;
    let motif_raw = next("motif")?;
    let purity_raw = next("purity_fraction")?;
    let ref_start_raw = next("ref_seq_start")?;
    let ref_seq_raw = next("ref_seq")?;
    if fields.next().is_some() {
        return Err(CatalogError::RowParse {
            line,
            reason: "too many columns (expected 7)".to_string(),
        });
    }

    let parse_u32 = |name: &str, raw: &str| -> Result<u32, CatalogError> {
        raw.parse::<u32>().map_err(|_| CatalogError::RowParse {
            line,
            reason: format!("column {name:?} = {raw:?} is not a u32"),
        })
    };
    let start = parse_u32("start", start_raw)?;
    let end = parse_u32("end", end_raw)?;
    let ref_seq_start = parse_u32("ref_seq_start", ref_start_raw)?;
    let purity_fraction = purity_raw
        .parse::<f32>()
        .map_err(|_| CatalogError::RowParse {
            line,
            reason: format!("column \"purity_fraction\" = {purity_raw:?} is not an f32"),
        })?;

    let motif = Motif::new(motif_raw.as_bytes())
        .map_err(|source| CatalogError::InvalidMotif { line, source })?;
    Locus::new(
        chrom.into_boxed_str(),
        start,
        end,
        motif,
        purity_fraction,
        ref_seq_raw.as_bytes().to_vec().into_boxed_slice(),
        ref_seq_start,
    )
    .map_err(|source| CatalogError::InvalidLocus { line, source })
}

/// Streaming catalog writer over any sink. Wraps the sink in a bgzip frame,
/// emits the `##` header + column header on construction, then one row per
/// [`CatalogWriter::write_locus`]. Rows must arrive already `(chrom, start)`
/// sorted (the coordinate-index requirement); the writer does not sort.
pub(crate) struct CatalogWriter<W: Write> {
    inner: bgzf::io::Writer<W>,
}

impl<W: Write> CatalogWriter<W> {
    /// Frame `sink` as bgzip and write the header block.
    pub(crate) fn new(sink: W, header: &CatalogHeader) -> Result<Self, CatalogError> {
        let mut inner = bgzf::io::Writer::new(sink);
        header
            .write_to(&mut inner)
            .map_err(|source| CatalogError::Io {
                context: "write catalog header",
                source,
            })?;
        Ok(Self { inner })
    }

    /// Append one locus as a data row.
    pub(crate) fn write_locus(&mut self, locus: &Locus) -> Result<(), CatalogError> {
        let row = locus_to_row(locus);
        writeln!(self.inner, "{row}").map_err(|source| CatalogError::Io {
            context: "write catalog row",
            source,
        })
    }

    /// Finalise the bgzip stream (flush + EOF block) and return the sink.
    pub(crate) fn finish(self) -> Result<W, CatalogError> {
        self.inner.finish().map_err(|source| CatalogError::Io {
            context: "finish catalog bgzf stream",
            source,
        })
    }
}

/// Streaming catalog reader over any bgzip source. Parses the header on
/// construction; [`CatalogReader::read_locus`] yields the next data row as a
/// [`Locus`]. (Region-query via a coordinate index is a later increment.)
pub(crate) struct CatalogReader<R: Read> {
    inner: BufReader<bgzf::io::Reader<R>>,
    header: CatalogHeader,
    line_no: usize,
}

impl<R: Read> CatalogReader<R> {
    /// Decode the bgzip frame and parse the header block.
    pub(crate) fn new(source: R) -> Result<Self, CatalogError> {
        let mut inner = BufReader::new(bgzf::io::Reader::new(source));
        let mut line_no = 0usize;
        let header = CatalogHeader::parse_from(&mut inner, &mut line_no)?;
        Ok(Self {
            inner,
            header,
            line_no,
        })
    }

    /// The catalog's metadata header (reference / md5 / params / `flank_bp`).
    pub(crate) fn header(&self) -> &CatalogHeader {
        &self.header
    }

    /// Read the next locus, or `None` at end of stream.
    pub(crate) fn read_locus(&mut self) -> Option<Result<Locus, CatalogError>> {
        let mut buf = String::new();
        loop {
            buf.clear();
            match self.inner.read_line(&mut buf) {
                Ok(0) => return None,
                Ok(_) => {
                    self.line_no += 1;
                    let line = buf.trim_end_matches(['\n', '\r']);
                    if line.is_empty() {
                        continue; // tolerate trailing/blank lines
                    }
                    return Some(row_to_locus(line, self.line_no));
                }
                Err(source) => {
                    return Some(Err(CatalogError::Io {
                        context: "read catalog row",
                        source,
                    }));
                }
            }
        }
    }

    /// Collect every remaining locus (convenience for whole-catalog consumers).
    pub(crate) fn read_all(&mut self) -> Result<Vec<Locus>, CatalogError> {
        let mut out = Vec::new();
        while let Some(locus) = self.read_locus() {
            out.push(locus?);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn sample_header() -> CatalogHeader {
        CatalogHeader {
            tool_version: "0.1.0".to_string(),
            reference: "ref.fa".to_string(),
            reference_md5: "0".repeat(32),
            trf_mod_version: "TRF-mod 0.0.0".to_string(),
            params: CatalogParams {
                flank_bp: 50,
                min_purity: 0.8,
                min_score: 20,
                bundle_threshold: 50,
            },
            date: "2026-06-15".to_string(),
        }
    }

    /// `chrom` len-12 tract `AT`×6 at [100,112) on a 0-based contig, with a
    /// 50 bp flank each side embedded (here a short hand-built ref_seq).
    fn sample_locus(chrom: &str, start: u32) -> Locus {
        let motif = Motif::new(b"AT").unwrap();
        let tract = b"ATATATATATAT"; // 12 bp = 6 copies of AT
        let flank = b"CGCG"; // 4 bp flanks for the fixture
        let mut ref_seq = Vec::new();
        ref_seq.extend_from_slice(flank);
        ref_seq.extend_from_slice(tract);
        ref_seq.extend_from_slice(flank);
        let ref_bytes_start = start - flank.len() as u32;
        Locus::new(
            chrom.to_string().into_boxed_str(),
            start,
            start + tract.len() as u32,
            motif,
            1.0,
            ref_seq.into_boxed_slice(),
            ref_bytes_start,
        )
        .unwrap()
    }

    #[test]
    fn locus_row_round_trips() {
        let locus = sample_locus("chr1", 100);
        let row = locus_to_row(&locus);
        // Column shape pinned: 7 tab fields, no embedded tabs.
        assert_eq!(row.split('\t').count(), 7);
        let back = row_to_locus(&row, 1).unwrap();
        assert_eq!(back, locus);
    }

    #[test]
    fn imperfect_purity_round_trips_exactly() {
        // f32 must survive the text round-trip bit-for-bit.
        let motif = Motif::new(b"AT").unwrap();
        let locus = Locus::new(
            "chr2".to_string().into_boxed_str(),
            10,
            22,
            motif,
            0.8333333,
            b"NNATATATATATATNN".to_vec().into_boxed_slice(),
            8,
        )
        .unwrap();
        let back = row_to_locus(&locus_to_row(&locus), 1).unwrap();
        assert_eq!(back.purity_fraction(), locus.purity_fraction());
        assert_eq!(back, locus);
    }

    #[test]
    fn writer_reader_round_trip_through_bgzf() {
        let header = sample_header();
        let loci = vec![
            sample_locus("chr1", 100),
            sample_locus("chr1", 500),
            sample_locus("chr2", 40),
        ];

        let mut writer = CatalogWriter::new(Cursor::new(Vec::<u8>::new()), &header).unwrap();
        for l in &loci {
            writer.write_locus(l).unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();

        let mut reader = CatalogReader::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.header(), &header);
        let read_back = reader.read_all().unwrap();
        assert_eq!(read_back, loci);
    }

    #[test]
    fn header_missing_required_key_is_rejected() {
        // Drop the `reference_md5` line from a valid header block.
        let header = sample_header();
        let mut raw = Vec::new();
        header.write_to(&mut raw).unwrap();
        let text = String::from_utf8(raw).unwrap();
        let mutated: String = text
            .lines()
            .filter(|l| !l.starts_with("## reference_md5"))
            .map(|l| format!("{l}\n"))
            .collect();
        // Re-frame as bgzf so the reader's decode path is exercised.
        let mut w = bgzf::io::Writer::new(Cursor::new(Vec::<u8>::new()));
        w.write_all(mutated.as_bytes()).unwrap();
        let bytes = w.finish().unwrap().into_inner();
        let err = match CatalogReader::new(Cursor::new(bytes)) {
            Ok(_) => panic!("missing md5 must fail"),
            Err(e) => e,
        };
        assert!(
            matches!(err, CatalogError::HeaderParse { ref reason } if reason.contains("reference_md5")),
            "expected HeaderParse naming reference_md5, got {err:?}"
        );
    }

    #[test]
    fn row_with_wrong_column_count_is_rejected() {
        let err = row_to_locus("chr1\t100\t112\tAT", 7).unwrap_err();
        assert!(
            matches!(err, CatalogError::RowParse { line: 7, .. }),
            "expected RowParse at line 7, got {err:?}"
        );
    }

    #[test]
    fn row_with_bad_motif_is_rejected() {
        // Empty motif is outside the SSR period range.
        let row = "chr1\t100\t112\t\t1\t96\tNNATATATATATATNN";
        let err = row_to_locus(row, 3).unwrap_err();
        assert!(
            matches!(err, CatalogError::InvalidMotif { line: 3, .. }),
            "expected InvalidMotif at line 3, got {err:?}"
        );
    }
}
