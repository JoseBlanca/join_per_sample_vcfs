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
use std::path::PathBuf;

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
}
