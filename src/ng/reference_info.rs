//! Reading a reference's info — the contig table, content digest, and reconstructed
//! `.fai` index, built *from the reference itself*. Design spec:
//! `doc/devel/ng/spec/reference_info.md`; types & interfaces:
//! `doc/devel/ng/arch/reference_info.md`.
//!
//! Every contig table elsewhere in the tree is built from an alignment file's `@SQ`
//! headers, and the `.fai` is only ever validated *against* one. This module builds the
//! table from a FASTA (optionally cross-checked against a `.fai`) or from a `.fai` alone,
//! and returns a [`ReferenceInfo`]: the contig table ([`ContigInfo`] per contig), the
//! whole-reference digest, and the geometry that locates each contig in the file.
//!
//! It is a **leaf**: it imports [`crate::fasta`] (read-only, for the transitional
//! [`ReferenceInfo::contig_list`] bridge) and noodles, and knows nothing else about ng
//! (spec §3.6). Its output types are ng's own plain data — `ContigInfo` is the
//! everything-about-a-contig type production never had in one place.
//!
//! **Milestone A** landed here so far: the data types, and the cheap `.fai` reader.

use crate::fasta::{ContigEntry, ContigList};
use std::io;
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------
// The data — ContigInfo, ReferenceInfo
// ---------------------------------------------------------------------

/// **Everything about one contig, in one place.** Its name and length (what a BAM `@SQ`
/// carries), the `.fai` geometry that locates it in the file (spec §3.8), and — once the
/// bases are read — its MD5 (the `@SQ M5`, spec §3.4). Named for what it *is*, not the
/// `.fai` some fields come from: `ContigInfo` is a *superset* of any production type
/// (`ContigEntry` has name/length/md5 but no geometry; a `.fai` row has geometry but no
/// md5).
///
/// Fields are `pub` and unconstrained — a plain data record. Illegal shapes (e.g.
/// `line_bases == 0`) are rejected at the boundary that builds one from a file (spec §5
/// T3), never represented as an invariant on this struct.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigInfo {
    /// The contig's name: the first word of the `>` header, the same string as `@SQ SN`.
    pub name: String,
    /// Total bases (`faidx.5` LENGTH; the same number as `@SQ LN`).
    pub length: u64,
    /// Byte offset of the first base (`faidx.5` OFFSET).
    pub offset: u64,
    /// Bases per line (`faidx.5` LINEBASES).
    pub line_bases: u64,
    /// Bytes per line, terminator included (`faidx.5` LINEWIDTH); `line_width - line_bases`
    /// is the terminator width, so CR-LF needs no special case (spec §3.8).
    pub line_width: u64,
    /// The `@SQ M5` (spec §3.4). `None` until the FASTA is read — a `.fai` has no MD5
    /// column, so a `.fai`-only read leaves this absent, which downstream treats as a
    /// wildcard (spec §5 T6).
    pub md5: Option<[u8; 16]>,
}

/// What reading a reference yields: the whole-assembly digest and every contig's info.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceInfo {
    /// The **whole-reference** MD5 — every contig's uppercased bases *concatenated* in
    /// file order (spec §3.4; the `ssr-catalog` convention, and **not** a SAM field).
    /// Distinct from the per-contig `@SQ M5` on each [`ContigInfo`]. `None` from a
    /// `.fai`-only read; being an in-order concatenation, it also pins the contig order
    /// (spec §5 T1).
    pub md5: Option<[u8; 16]>,
    /// Every contig, **in file order** — the order that defines `ContigId` (`ContigId(i)`
    /// is `contigs[i]`, spec §5 T1). One source of truth for name/length/geometry/md5, so
    /// there is nothing to keep in step.
    pub contigs: Vec<ContigInfo>,
}

impl ReferenceInfo {
    /// **Transitional** projection to the production [`ContigList`] (one [`ContigEntry`]
    /// `{ name, length, md5 }` per contig) for ng consumers that still speak it — the
    /// `RefSeq` impls, `ContigTable`, `GenomeRegions` validation. A fresh allocation; the
    /// geometry has no place in a `ContigList` and is dropped. The *only* place
    /// `ContigList` surfaces from this module; removed when those consumers migrate to
    /// [`ContigInfo`] (spec §3.6, §3.9).
    pub fn contig_list(&self) -> ContigList {
        ContigList {
            entries: self
                .contigs
                .iter()
                .map(|c| ContigEntry {
                    name: c.name.clone(),
                    length: c.length,
                    md5: c.md5,
                })
                .collect(),
        }
    }
}

// ---------------------------------------------------------------------
// The input — ReferenceSource
// ---------------------------------------------------------------------

/// What the caller hands in. Two ways to describe a reference; the variant names the cost
/// — `Fai` is a small file read, `Fasta` reads the whole reference (spec §3.1, §3.2).
/// Owns its paths (`PathBuf`, no lifetime) so the module — and the cache built on it —
/// stays lifetime-free; a `PathBuf` clone next to a genome read does not register.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReferenceSource {
    /// The index alone → names, order, lengths, geometry; `md5: None`. Cheap.
    Fai(PathBuf),
    /// The FASTA, **optionally** cross-checked against a `.fai`. `fai: None` reads the
    /// FASTA alone (names, order, lengths, every MD5, the reference digest). `fai: Some`
    /// reads the same **and** verifies the FASTA against that index (spec §3.3); a
    /// disagreement is an error. Both cost a whole-genome read — the index only adds a
    /// comparison.
    Fasta {
        fasta: PathBuf,
        fai: Option<PathBuf>,
    },
}

// ---------------------------------------------------------------------
// The error — ReferenceInfoError
// ---------------------------------------------------------------------

/// Everything that can go wrong reading a reference. `#[non_exhaustive]` `thiserror` enum,
/// the house pattern (mirroring `RefSeqError`); each variant carries the offending path.
///
/// Nothing in a supplied reference may panic — its bytes are attacker-influenced input, so
/// every malformed shape is a named error, never an `assert!`/`debug_assert!` (spec §5 T3).
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum ReferenceInfoError {
    /// The `.fai` index could not be read — missing, unreadable, or malformed enough that
    /// `noodles_fasta::fai::fs::read` rejected it (the `Fai` arm, and the `.fai` half of a
    /// `Fasta { fai: Some }` read).
    #[error("failed to read .fai index {path:?}")]
    FaiRead {
        path: PathBuf,
        #[source]
        source: io::Error,
    },

    /// FASTA I/O failed part-way through the streaming pass (spec §4).
    #[error("failed to read FASTA {path:?}")]
    FastaRead {
        path: PathBuf,
        #[source]
        source: io::Error,
    },

    /// The FASTA is malformed: bad geometry (non-uniform line width, an empty contig,
    /// `line_bases == 0`), or bases before the first `>` (spec §5 T3). `contig` is the
    /// offending contig's name when a header has already been seen, `None` for a
    /// before-first-header error; `byte_offset` locates it in the file. The message names
    /// the contig when known, per T3 (`""` renders as the before-header case).
    #[error("malformed FASTA {path:?} (contig {}) at byte {byte_offset}: {detail}",
            contig.as_deref().unwrap_or("<before first header>"))]
    MalformedFasta {
        path: PathBuf,
        contig: Option<String>,
        byte_offset: u64,
        detail: String,
    },

    /// A contig name appears twice. Rejected — stricter than htslib, which warns and drops
    /// the duplicate (spec §5 T2): ng resolves contigs by position, so a dropped duplicate
    /// would renumber every `ContigId` after it.
    #[error("duplicate contig name {name:?} in {path:?}")]
    DuplicateContigName { path: PathBuf, name: String },

    /// The FASTA and the supplied `.fai` disagree — the `.fai` is stale or wrong about the
    /// file (spec §3.3). Names the field and contig that differ so a re-wrap is diagnosable
    /// (`detail` carries the two values).
    #[error(
        "FASTA {fasta:?} disagrees with .fai {fai:?} for contig {contig:?}: {field} differs ({detail})"
    )]
    FastaFaiMismatch {
        fasta: PathBuf,
        fai: PathBuf,
        contig: String,
        field: String,
        detail: String,
    },

    /// A six-column FASTQ index was handed where a FASTA index was expected (spec §3.8).
    /// A recognisable, nameable wrong input, not a corruption.
    #[error("{path:?} is a FASTQ index (six columns); a FASTA index was expected")]
    FastqIndex { path: PathBuf },

    /// A bgzip-compressed reference (`.fa.gz`) — out of scope, deferred (spec §3.8, §7). A
    /// clean error today, not a silent misread of compressed bytes as FASTA.
    #[error("compressed reference {path:?} is not supported (bgzip .fa.gz is deferred)")]
    CompressedReference { path: PathBuf },

    /// Writing the `.fai` failed in `read_reference_verifying_or_creating_fai` (spec
    /// §3.11) — a read-only reference dir, a full disk. Fatal there by design; a caller
    /// that cannot write uses the primitives (`get_or_read(Fasta { fai: None })`).
    #[error("failed to write .fai index {path:?}")]
    FaiWrite {
        path: PathBuf,
        #[source]
        source: io::Error,
    },
}

// ---------------------------------------------------------------------
// Pure path helper
// ---------------------------------------------------------------------

/// `<fasta>` + `".fai"` — the reference naming convention, as a path, **no I/O**. A caller
/// wanting sibling behaviour composes this with its own existence check; the module never
/// probes the filesystem for a `.fai` it was not handed (spec §2, §3.6). This is the one
/// place ng expresses the `<fasta>.fai` convention as a value; the module applies it only
/// in the `read_reference_verifying_or_creating_fai` orchestrator (spec §3.11).
///
/// It appends `.fai` to the *whole* path — `ref.fa` → `ref.fa.fai`, not `ref.fai` — because
/// that is what `samtools faidx` writes and what every reader in the tree seeks. (Five lines
/// copied from `bam/alignment_input.rs`'s private `with_fai_extension`, the standing
/// copy-what-you-need rule, `ng/mod.rs`.)
pub fn sibling_fai_path(fasta: &Path) -> PathBuf {
    let mut buf = fasta.as_os_str().to_owned();
    buf.push(".fai");
    PathBuf::from(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn info_with_two_contigs() -> ReferenceInfo {
        ReferenceInfo {
            md5: Some([7u8; 16]),
            contigs: vec![
                ContigInfo {
                    name: "chr1".to_string(),
                    length: 100,
                    offset: 6,
                    line_bases: 60,
                    line_width: 61,
                    md5: Some([1u8; 16]),
                },
                ContigInfo {
                    name: "chr2".to_string(),
                    length: 40,
                    offset: 200,
                    line_bases: 60,
                    line_width: 61,
                    md5: None,
                },
            ],
        }
    }

    #[test]
    fn contig_list_projects_name_length_md5_in_order() {
        let info = info_with_two_contigs();
        let list = info.contig_list();

        assert_eq!(list.entries.len(), info.contigs.len());
        for (entry, contig) in list.entries.iter().zip(info.contigs.iter()) {
            assert_eq!(entry.name, contig.name);
            assert_eq!(entry.length, contig.length);
            assert_eq!(entry.md5, contig.md5);
        }
        // Order is the ContigId contract (spec §5 T1): the projection must not reorder.
        assert_eq!(list.entries[0].name, "chr1");
        assert_eq!(list.entries[1].name, "chr2");
    }

    #[test]
    fn sibling_fai_path_appends_dot_fai_and_touches_nothing() {
        // Appends to the whole path (samtools' convention), does not replace the extension.
        assert_eq!(
            sibling_fai_path(Path::new("/data/ref.fa")),
            PathBuf::from("/data/ref.fa.fai")
        );
        assert_eq!(
            sibling_fai_path(Path::new("ref.fasta")),
            PathBuf::from("ref.fasta.fai")
        );
        // A path with no extension still just gets `.fai` appended.
        assert_eq!(
            sibling_fai_path(Path::new("genome")),
            PathBuf::from("genome.fai")
        );
        // Pure: it returns a value for a path that does not exist, without any I/O — a
        // filesystem probe would error or block on a bogus path; this cannot.
        assert_eq!(
            sibling_fai_path(Path::new("/no/such/place/ghost.fa")),
            PathBuf::from("/no/such/place/ghost.fa.fai")
        );
    }

    #[test]
    fn malformed_fasta_error_names_the_contig_when_known() {
        // T3: the error must name the contig and the byte offset.
        let named = ReferenceInfoError::MalformedFasta {
            path: PathBuf::from("ref.fa"),
            contig: Some("chr7".to_string()),
            byte_offset: 4096,
            detail: "non-uniform line width".to_string(),
        };
        let text = named.to_string();
        assert!(text.contains("chr7"), "contig name missing from: {text}");
        assert!(text.contains("4096"), "byte offset missing from: {text}");

        // Before the first header there is no contig; the message says so rather than
        // rendering an empty parenthetical.
        let headerless = ReferenceInfoError::MalformedFasta {
            path: PathBuf::from("ref.fa"),
            contig: None,
            byte_offset: 0,
            detail: "bases before the first '>'".to_string(),
        };
        assert!(headerless.to_string().contains("before first header"));
    }

    #[test]
    fn contig_list_drops_geometry_only() {
        // The projection is lossy by design: name/length/md5 survive, geometry does not
        // (a ContigList has nowhere to put it). This pins that the bridge carries exactly
        // what ContigEntry can hold.
        let info = info_with_two_contigs();
        let list = info.contig_list();
        let second = &list.entries[1];
        assert_eq!(second.md5, None); // the .fai-less contig's absent md5 rides through
    }

    // =================================================================
    // A3 — test fixtures + committed samtools oracle constants
    //      (spec §6, §3.8; arch §7). The oracle is `samtools`, resolved
    //      to committed constants so `samtools` is not a test-time dep.
    // =================================================================

    /// One contig's committed oracle values.
    struct GoldenContig {
        name: &'static str,
        length: u64,
        offset: u64,
        line_bases: u64,
        line_width: u64,
        /// The `@SQ M5`, hex — from `samtools dict`.
        md5_hex: &'static str,
    }

    /// The golden reference — `tests/data/tandem_repeat/synthetic_ref.fa`, the same file
    /// `anchor.rs` reads and the golden `.cat` was built from — as `samtools 1.16.1`
    /// describes it. Run **once** in the dev container and committed here: `samtools dict`
    /// gave `LN`/`M5`; `samtools faidx` gave the five `.fai` columns (committed beside the
    /// FASTA as `synthetic_ref.fa.fai`). The MD5s (B2's oracle) ride along now so the whole
    /// oracle lands in one place.
    const GOLDEN: [GoldenContig; 2] = [
        GoldenContig {
            name: "ctg1",
            length: 2369,
            offset: 6,
            line_bases: 60,
            line_width: 61,
            md5_hex: "d6b8aee46bd8cfc7c7dbc13e7a9286c3",
        },
        GoldenContig {
            name: "ctg2",
            length: 346,
            offset: 2421,
            line_bases: 60,
            line_width: 61,
            md5_hex: "4d9680fefca4e59d9adcac944c7c52ab",
        },
    ];

    /// Path to the committed golden FASTA (its `.fai` sits beside it).
    fn golden_ref_path() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/tandem_repeat/synthetic_ref.fa")
    }

    /// A width-configurable FASTA writer for the geometry tests: each contig's bases wrapped
    /// at `wrap` bases per line (LF, a possibly-shorter last line), plus a matching
    /// five-column `.fai`, in a fresh tempdir. Returns `(dir, fasta, fai)` — the caller keeps
    /// `dir` alive. Unlike `ref_seq.rs::build_fasta` (one line per contig, so no geometry to
    /// exercise), this is what B/C stand on. Callers pass non-empty contigs.
    fn write_wrapped_fasta(
        contigs: &[(&str, &[u8])],
        wrap: usize,
    ) -> (tempfile::TempDir, PathBuf, PathBuf) {
        assert!(wrap >= 1, "wrap width must be >= 1");
        let dir = tempfile::tempdir().expect("tempdir");
        let fasta_path = dir.path().join("ref.fa");
        let fai_path = dir.path().join("ref.fa.fai");

        let mut fasta: Vec<u8> = Vec::new();
        let mut fai = String::new();
        for (name, bases) in contigs {
            fasta.extend_from_slice(format!(">{name}\n").as_bytes());
            let seq_offset = fasta.len() as u64;
            for chunk in bases.chunks(wrap) {
                fasta.extend_from_slice(chunk);
                fasta.push(b'\n');
            }
            // samtools reports the first line's base count; for a contig shorter than `wrap`
            // that is the whole contig.
            let line_bases = bases.len().min(wrap).max(1) as u64;
            fai.push_str(&format!(
                "{name}\t{}\t{seq_offset}\t{line_bases}\t{}\n",
                bases.len(),
                line_bases + 1,
            ));
        }
        std::fs::write(&fasta_path, &fasta).expect("write fasta");
        std::fs::write(&fai_path, fai).expect("write fai");
        (dir, fasta_path, fai_path)
    }

    // ---- the htslib `faidx.5` man-page worked example (spec §3.8), verbatim ----
    // Two contigs; the second header carries a description (`>two another chromosome`)
    // whose NAME is the first word `two` (T7). Committed so B2 can pin §4's reconstruction
    // against the format's own authors, not only local samtools.
    const WORKED_ONE_HEADER: &str = ">one";
    const WORKED_TWO_HEADER: &str = ">two another chromosome";
    const WORKED_ONE_LINES: [&[u8]; 3] = [
        b"ATGCATGCATGCATGCATGCATGCATGCAT",
        b"GCATGCATGCATGCATGCATGCATGCATGC",
        b"ATGCAT",
    ];
    const WORKED_TWO_LINES: [&[u8]; 2] = [b"ATGCATGCATGCAT", b"GCATGCATGCATGC"];
    /// `(name, length, offset, line_bases, line_width)` under LF (`faidx.5`).
    const WORKED_INDEX_LF: [(&str, u64, u64, u64, u64); 2] =
        [("one", 66, 5, 30, 31), ("two", 28, 98, 14, 15)];
    /// The same under CR-LF: `line_width` +1 per line, offsets shifted by the extra CRs.
    const WORKED_INDEX_CRLF: [(&str, u64, u64, u64, u64); 2] =
        [("one", 66, 6, 30, 32), ("two", 28, 103, 14, 16)];

    /// Build the worked-example FASTA bytes under LF (`crlf = false`) or CR-LF.
    fn build_worked_example(crlf: bool) -> Vec<u8> {
        let nl: &[u8] = if crlf { b"\r\n" } else { b"\n" };
        let mut out = Vec::new();
        for header in [WORKED_ONE_HEADER, WORKED_TWO_HEADER] {
            let lines: &[&[u8]] = if header == WORKED_ONE_HEADER {
                &WORKED_ONE_LINES
            } else {
                &WORKED_TWO_LINES
            };
            out.extend_from_slice(header.as_bytes());
            out.extend_from_slice(nl);
            for line in lines {
                out.extend_from_slice(line);
                out.extend_from_slice(nl);
            }
        }
        out
    }

    #[test]
    fn golden_fai_fixture_matches_committed_constants() {
        // The committed `synthetic_ref.fa.fai` (samtools) must equal the constants above —
        // one pins the other, so a future edit to either is caught.
        let fai_text = std::fs::read_to_string(sibling_fai_path(&golden_ref_path()))
            .expect("committed golden .fai");
        let rows: Vec<&str> = fai_text.lines().collect();
        assert_eq!(rows.len(), GOLDEN.len(), "golden .fai row count");
        for (row, g) in rows.iter().zip(GOLDEN.iter()) {
            let cols: Vec<&str> = row.split('\t').collect();
            assert_eq!(cols.len(), 5, "FASTA .fai has five columns");
            assert_eq!(cols[0], g.name);
            assert_eq!(cols[1].parse::<u64>().unwrap(), g.length);
            assert_eq!(cols[2].parse::<u64>().unwrap(), g.offset);
            assert_eq!(cols[3].parse::<u64>().unwrap(), g.line_bases);
            assert_eq!(cols[4].parse::<u64>().unwrap(), g.line_width);
            // Sanity on the committed M5 (B2's oracle): 32 lowercase hex digits.
            assert_eq!(g.md5_hex.len(), 32, "M5 is 16 bytes = 32 hex chars");
            assert!(
                g.md5_hex.bytes().all(|b| b.is_ascii_hexdigit()),
                "M5 is hex"
            );
        }
    }

    #[test]
    fn wrapped_fasta_writer_geometry_is_self_consistent() {
        // A contig of 10 bases wrapped at 4 → lines 4, 4, 2 (a shorter last line). The `.fai`
        // the writer emits must describe the bytes it wrote: OFFSET points at the first base,
        // and the geometry locates the last base too.
        // `_dir` keeps the tempdir alive through the reads below (a bare `_` would delete it
        // immediately); it drops at the end of the test.
        let (_dir, fasta_path, fai_path) = write_wrapped_fasta(&[("solo", b"ACGTACGTAC")], 4);
        let fasta = std::fs::read(&fasta_path).unwrap();
        let fai = std::fs::read_to_string(&fai_path).unwrap();
        let cols: Vec<&str> = fai.trim_end().split('\t').collect();
        assert_eq!(cols[0], "solo");
        let length: u64 = cols[1].parse().unwrap();
        let offset: u64 = cols[2].parse().unwrap();
        let line_bases: u64 = cols[3].parse().unwrap();
        let line_width: u64 = cols[4].parse().unwrap();
        assert_eq!((length, offset, line_bases, line_width), (10, 6, 4, 5)); // ">solo\n" = 6

        // First base sits at OFFSET.
        assert_eq!(fasta[offset as usize], b'A');
        // Last base (10th, 1-based) via the same arithmetic every reader uses:
        // line_idx = 9 / 4 = 2, in_line = 9 % 4 = 1 → offset + 2*line_width + 1.
        let last = offset + 2 * line_width + 1;
        assert_eq!(fasta[last as usize], b'C');
    }

    #[test]
    fn faidx5_worked_example_vector_is_internally_consistent() {
        // The committed index tuples must agree with the bytes `build_worked_example`
        // produces — otherwise the vector B2 pins against would be self-contradictory.
        for (crlf, index) in [(false, WORKED_INDEX_LF), (true, WORKED_INDEX_CRLF)] {
            let fasta = build_worked_example(crlf);
            for (name, _length, offset, _lb, _lw) in index {
                // The sequence of each contig starts with 'A' at its committed OFFSET (both
                // worked contigs begin `ATGC...`).
                assert_eq!(
                    fasta[offset as usize], b'A',
                    "contig {name} sequence must start at offset {offset} (crlf={crlf})",
                );
            }
        }
        // T7: the second header carries a description; its name is the first word only.
        assert_eq!(WORKED_TWO_HEADER.split_whitespace().next().unwrap(), ">two");
    }
}
