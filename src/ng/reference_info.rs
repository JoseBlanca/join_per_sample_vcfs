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
//! Landed so far: the data types and the cheap `.fai` reader (Milestone A), and the
//! from-byte-zero FASTA streaming pass (Milestone B) — `Fasta { fai: None }` reconstructs
//! geometry + per-contig and whole-reference MD5 in one buffer; `Fasta { fai: Some }` reads
//! the same and proves the supplied `.fai` describes the same genome (spec §3.3). And
//! `write_fai` writes a `.fai` byte-identical to `samtools faidx`, atomically (Milestone C).
//! `ReferenceInfoCache` is the compute-once, single-flight cache over the reader (Milestone D).
//! `read_fai_verify_in_background` and `read_reference_verifying_or_creating_fai` are the two
//! composing entry points — background verify, and verify-or-index (Milestone E). Module complete.

use crate::fasta::{ContigEntry, ContigList};
use md5::{Digest, Md5};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::SystemTime;

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
    /// The FASTA this was read from, when it was read from one.
    ///
    /// `None` for a `.fai`-only read, and that absence is **meaningful rather
    /// than incidental**: a `.fai` describes a genome's geometry but holds no
    /// bases, so a consumer that needs the sequence itself — CRAM decoding is
    /// the one that exists — cannot proceed from it. Carrying the path lets
    /// such a consumer build a `fasta::Repository`, and lets it fail with a
    /// clear message when it cannot (owner, 2026-07-20).
    ///
    /// Not used for anything this module does; it is provenance passed through
    /// for consumers that need to go back to the bases.
    pub fasta_path: Option<PathBuf>,
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
    /// The `.fai` index could not be read: missing, unreadable, malformed enough that
    /// `noodles_fasta::fai::fs::read` rejected it, or it parsed but failed ng's field
    /// guards (`line_bases > 0`, `line_width >= line_bases` — the `ContigFai::validate`
    /// checks, copied). The `source` carries the specifics (contig + reason for a guard
    /// failure). The `Fai` arm, and the `.fai` half of a `Fasta { fai: Some }` read.
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
// The reader
// ---------------------------------------------------------------------

/// Read a reference. **Pure** — no shared state, no memoisation, no filesystem probing
/// beyond opening the paths it was handed; every call does the full work the source names
/// (spec §3.1). A single-call consumer calls this directly; one that reads a reference
/// repeatedly or from several threads holds a `ReferenceInfoCache` (Milestone D) instead.
///
/// Milestone A implements the cheap `Fai` arm; the `Fasta` streaming pass lands in
/// Milestone B (spec §4).
pub fn read_reference_info(source: ReferenceSource) -> Result<ReferenceInfo, ReferenceInfoError> {
    match source {
        ReferenceSource::Fai(path) => read_fai(&path),
        // Read the FASTA alone: names, order, lengths, every MD5, the reference digest.
        ReferenceSource::Fasta { fasta, fai: None } => read_fasta(&fasta),
        // Read the FASTA and prove the supplied `.fai` describes the same genome (spec §3.3).
        ReferenceSource::Fasta {
            fasta,
            fai: Some(fai),
        } => read_fasta_verifying(&fasta, &fai),
    }
}

/// The `Fasta { fai: Some }` arm: run the pass, then **prove** the supplied `.fai` describes
/// the same genome (spec §3.3). "They match" means the whole index — all five `faidx.5`
/// columns, reconstructed by the pass, compared field-for-field against the on-disk `.fai`.
/// Names and lengths alone are not enough: a `fold -w` re-wrap keeps both and breaks the
/// geometry every reader seeks by, so the check is exactly what a fetch depends on
/// (`offset`, `length`, `line_bases`, `line_width`). A reordering is caught here too — the
/// comparison is position-wise, so a permuted `.fai` disagrees on `name` (T1).
fn read_fasta_verifying(fasta: &Path, fai: &Path) -> Result<ReferenceInfo, ReferenceInfoError> {
    // Read the cheap `.fai` first so a missing or malformed index (T4, T2, FASTQ, a bad
    // field) fails fast, before the expensive whole-genome pass.
    let index = read_fai(fai)?;
    let info = read_fasta(fasta)?;

    if info.contigs.len() != index.contigs.len() {
        return Err(ReferenceInfoError::FastaFaiMismatch {
            fasta: fasta.to_path_buf(),
            fai: fai.to_path_buf(),
            contig: "(whole file)".to_string(),
            field: "contig count".to_string(),
            detail: format!(
                "FASTA has {}, .fai has {}",
                info.contigs.len(),
                index.contigs.len()
            ),
        });
    }
    for (from_fasta, from_fai) in info.contigs.iter().zip(index.contigs.iter()) {
        if let Some((field, detail)) = first_fai_field_disagreement(from_fasta, from_fai) {
            return Err(ReferenceInfoError::FastaFaiMismatch {
                fasta: fasta.to_path_buf(),
                fai: fai.to_path_buf(),
                contig: from_fasta.name.clone(),
                field: field.to_string(),
                detail,
            });
        }
    }
    Ok(info)
}

/// The first of the five `faidx.5` fields on which the FASTA's reconstruction and the `.fai`
/// disagree, with a `fasta … / .fai …` detail — or `None` when they agree. `md5` is ignored
/// (a `.fai` has none). Order matters: `name` first so a reordering reports `name`, then the
/// geometry a re-wrap breaks, so the counterexample reports `line_bases` (spec §3.3).
fn first_fai_field_disagreement(
    from_fasta: &ContigInfo,
    from_fai: &ContigInfo,
) -> Option<(&'static str, String)> {
    if from_fasta.name != from_fai.name {
        return Some((
            "name",
            format!("fasta {:?}, .fai {:?}", from_fasta.name, from_fai.name),
        ));
    }
    let fields: [(&'static str, u64, u64); 4] = [
        ("length", from_fasta.length, from_fai.length),
        ("offset", from_fasta.offset, from_fai.offset),
        ("line_bases", from_fasta.line_bases, from_fai.line_bases),
        ("line_width", from_fasta.line_width, from_fai.line_width),
    ];
    for (name, fasta_value, fai_value) in fields {
        if fasta_value != fai_value {
            return Some((name, format!("fasta {fasta_value}, .fai {fai_value}")));
        }
    }
    None
}

/// The `Fai` arm: parse a `.fai` into `ContigInfo`s (`md5: None` — a `.fai` has no MD5
/// column). Names, order, lengths and geometry only, at the cost of one small file read
/// (spec §3.1). Rejects a FASTQ index (six columns, §3.8), a duplicate contig name (T2),
/// and a `.fai` failing the field guards (T3).
fn read_fai(path: &Path) -> Result<ReferenceInfo, ReferenceInfoError> {
    // Detect a FASTQ index (six columns) *before* the numeric parse. noodles splits each
    // line into five fields (`splitn(5, '\t')`), so a sixth column is folded into
    // `line_width` and rejected as a generic parse error — losing the diagnosis. htslib
    // distinguishes the two indices the same way, by column count (spec §3.8).
    let text = std::fs::read_to_string(path).map_err(|source| ReferenceInfoError::FaiRead {
        path: path.to_path_buf(),
        source,
    })?;
    for line in text.lines().filter(|l| !l.is_empty()) {
        if line.split('\t').count() >= 6 {
            return Err(ReferenceInfoError::FastqIndex {
                path: path.to_path_buf(),
            });
        }
    }

    // Reuse noodles for the authoritative numeric parse — the same call the rest of the
    // tree makes (`raw_chrom_reader.rs`, `alignment_input.rs`).
    let index =
        noodles_fasta::fai::fs::read(path).map_err(|source| ReferenceInfoError::FaiRead {
            path: path.to_path_buf(),
            source,
        })?;

    let records: &[noodles_fasta::fai::Record] = index.as_ref();
    let mut contigs = Vec::with_capacity(records.len());
    let mut seen: HashSet<&[u8]> = HashSet::with_capacity(records.len());
    for record in records {
        let name_bytes = AsRef::<[u8]>::as_ref(record.name());
        // T2: a duplicate contig name is rejected, not warned-and-dropped like htslib —
        // ng resolves contigs by position, so dropping one renumbers every ContigId.
        if !seen.insert(name_bytes) {
            return Err(ReferenceInfoError::DuplicateContigName {
                path: path.to_path_buf(),
                name: String::from_utf8_lossy(name_bytes).into_owned(),
            });
        }
        let name = String::from_utf8_lossy(name_bytes).into_owned();
        let line_bases = record.line_bases();
        let line_width = record.line_width();
        // The field guards `ContigFai::validate` enforces (copied, spec §3.8): a
        // `line_bases` of 0 divides by zero in every reader's offset arithmetic, and
        // `line_width < line_bases` cannot hold (the width includes the terminator).
        if line_bases == 0 {
            return Err(fai_field_error(
                path,
                &name,
                "line_bases = 0 (would divide-by-zero in offset arithmetic)".to_string(),
            ));
        }
        if line_width < line_bases {
            return Err(fai_field_error(
                path,
                &name,
                format!(
                    "line_width ({line_width}) < line_bases ({line_bases}) — line_width must \
                     include the trailing newline"
                ),
            ));
        }
        contigs.push(ContigInfo {
            name,
            length: record.length(),
            offset: record.offset(),
            line_bases,
            line_width,
            md5: None,
        });
    }
    Ok(ReferenceInfo {
        md5: None,
        contigs,
        // A `.fai` holds no bases, so there is no FASTA to hand on.
        fasta_path: None,
    })
}

/// A `.fai` field-guard failure, as a `FaiRead` carrying a synthesised `InvalidData` error
/// (the shape `RawChromReader` wraps `ContigFai::validate`'s error in).
fn fai_field_error(path: &Path, contig: &str, detail: String) -> ReferenceInfoError {
    ReferenceInfoError::FaiRead {
        path: path.to_path_buf(),
        source: io::Error::new(
            io::ErrorKind::InvalidData,
            format!("malformed .fai for contig {contig}: {detail}"),
        ),
    }
}

// ---------------------------------------------------------------------
// The FASTA streaming pass (spec §4) — the heart
// ---------------------------------------------------------------------

/// One read window. Resident memory is one of these regardless of contig size — the pass
/// never holds a whole contig (spec §4). The same size `compute_contig_md5_streaming` uses.
const FASTA_PASS_BUFFER_SIZE: usize = 64 * 1024;

/// The `Fasta { fai: None }` arm: read the FASTA from **byte zero**, in one streaming pass,
/// reconstructing every contig's geometry and MD5 and the whole-reference digest (spec §4).
/// From byte zero, not seeking by an index, because a `.fai`-driven reader can only confirm
/// the index agrees with itself (spec §3.3 — the circular check). One buffer, never a whole
/// contig.
fn read_fasta(path: &Path) -> Result<ReferenceInfo, ReferenceInfoError> {
    let mut info = read_fasta_bases(path)?;
    // The bases came from here, so a consumer that needs them (CRAM decoding)
    // can get back to them.
    info.fasta_path = Some(path.to_path_buf());
    Ok(info)
}

fn read_fasta_bases(path: &Path) -> Result<ReferenceInfo, ReferenceInfoError> {
    let mut file = File::open(path).map_err(|source| ReferenceInfoError::FastaRead {
        path: path.to_path_buf(),
        source,
    })?;
    let mut buf = [0u8; FASTA_PASS_BUFFER_SIZE];
    let mut pass = FastaPass::new(path);
    let mut first_window = true;
    loop {
        let n = file
            .read(&mut buf)
            .map_err(|source| ReferenceInfoError::FastaRead {
                path: path.to_path_buf(),
                source,
            })?;
        if n == 0 {
            break;
        }
        if first_window {
            first_window = false;
            // A bgzip `.fa.gz` starts with the gzip magic; its `.fai` offsets index the
            // *uncompressed* stream, so parsing the raw bytes as FASTA is nonsense. Reject
            // it cleanly rather than misread (spec §3.8; the reader is deferred, §7).
            if n >= 2 && buf[0] == 0x1f && buf[1] == 0x8b {
                return Err(ReferenceInfoError::CompressedReference {
                    path: path.to_path_buf(),
                });
            }
        }
        for &b in &buf[..n] {
            pass.push_byte(b)?;
        }
    }
    pass.finish()
}

/// Where the pass is in the file's grammar.
enum PassMode {
    /// Before the first `>` header.
    Start,
    /// Inside a definition line (after `>`, before its terminator).
    Header,
    /// Reading a contig's sequence lines.
    Sequence,
}

/// The from-byte-zero streaming state. Holds only bounded state — the current contig's
/// name, the running line counters, the two MD5 states, and one flush buffer — never a
/// whole contig. `offset`/geometry are reconstructed the way htslib's own indexer does, so
/// a `.fai` comparison (B3) is like-for-like (spec §4).
struct FastaPass<'a> {
    path: &'a Path,
    /// Byte offset of the current byte being processed (0-based).
    pos: u64,
    mode: PassMode,
    /// True at file start and immediately after a `\n` — the only place a `>` opens a header.
    at_line_start: bool,

    /// The whole-reference digest: every contig's uppercased bases, concatenated in file
    /// order (spec §3.4).
    reference_md5: Md5,
    contigs: Vec<ContigInfo>,
    seen_names: HashSet<String>,

    // ---- current contig (meaningful in Header / Sequence) ----
    name_bytes: Vec<u8>,
    /// Set once the name's first word ends (at whitespace); the rest of the header is skipped.
    name_done: bool,
    /// Byte offset of the contig's first base (`faidx.5` OFFSET).
    seq_offset: u64,
    length: u64,
    contig_md5: Md5,
    /// The first sequence line's bases; fixes the geometry every later line must match.
    /// `None` until the first line completes.
    line_bases0: Option<u64>,
    line_width0: u64,
    cur_line_bases: u64,
    cur_line_width: u64,
    /// A line shorter than `line_bases0` was seen; only the *last* line may be short, so a
    /// further sequence line after one is a non-uniform-geometry error (spec §3.8, T3).
    short_line_seen: bool,

    /// Uppercased bases batched before feeding the MD5s (one `update` per full window
    /// instead of per byte). Reused across the whole pass.
    upper: Vec<u8>,
}

impl<'a> FastaPass<'a> {
    fn new(path: &'a Path) -> Self {
        FastaPass {
            path,
            pos: 0,
            mode: PassMode::Start,
            at_line_start: true,
            reference_md5: Md5::new(),
            contigs: Vec::new(),
            seen_names: HashSet::new(),
            name_bytes: Vec::new(),
            name_done: false,
            seq_offset: 0,
            length: 0,
            contig_md5: Md5::new(),
            line_bases0: None,
            line_width0: 0,
            cur_line_bases: 0,
            cur_line_width: 0,
            short_line_seen: false,
            upper: Vec::with_capacity(FASTA_PASS_BUFFER_SIZE),
        }
    }

    fn push_byte(&mut self, b: u8) -> Result<(), ReferenceInfoError> {
        match self.mode {
            PassMode::Start => {
                if b == b'>' {
                    self.begin_header();
                } else if b <= 0x20 {
                    // Tolerate leading blank lines / whitespace before the first header.
                    if b == b'\n' {
                        self.at_line_start = true;
                    }
                } else {
                    return Err(self.malformed("sequence data before the first '>' header"));
                }
            }
            PassMode::Header => {
                if b == b'\n' {
                    // A definition with an empty name is malformed (spec §5 T3).
                    if self.name_bytes.is_empty() {
                        return Err(self.malformed("definition line has an empty contig name"));
                    }
                    // The definition line ends; the next byte is the contig's first base.
                    self.seq_offset = self.pos + 1;
                    self.mode = PassMode::Sequence;
                    self.at_line_start = true;
                    self.cur_line_bases = 0;
                    self.cur_line_width = 0;
                } else if !self.name_done {
                    // The name is the first word (spec §5 T7): stop at whitespace.
                    if b == b' ' || b == b'\t' || b == b'\r' {
                        self.name_done = true;
                    } else {
                        self.name_bytes.push(b);
                    }
                }
                // else: past the name, inside the description — skipped until '\n'.
            }
            PassMode::Sequence => {
                if self.at_line_start && b == b'>' {
                    // This contig ends and a new one begins.
                    self.finalize_contig()?;
                    self.begin_header();
                } else if b == b'\n' {
                    self.cur_line_width += 1;
                    self.complete_line(false)?;
                    self.at_line_start = true;
                } else {
                    self.at_line_start = false;
                    self.cur_line_width += 1;
                    // The single predicate (spec §3.4/§4): a base is `[0x21, 0x7E]`
                    // (printable, non-space). `\r`, spaces and tabs fall out of the count
                    // and out of the MD5 for free.
                    if (0x21..=0x7e).contains(&b) {
                        self.cur_line_bases += 1;
                        self.length += 1;
                        self.push_base(b);
                    }
                }
            }
        }
        self.pos += 1;
        Ok(())
    }

    /// Start a new contig's header (the `>` has just been seen but is not part of the name).
    fn begin_header(&mut self) {
        self.mode = PassMode::Header;
        self.at_line_start = false;
        self.name_bytes.clear();
        self.name_done = false;
        self.length = 0;
        self.line_bases0 = None;
        self.line_width0 = 0;
        self.cur_line_bases = 0;
        self.cur_line_width = 0;
        self.short_line_seen = false;
        self.contig_md5 = Md5::new();
    }

    /// Feed one base byte (already known to be `[0x21, 0x7E]`), uppercased, to both digests
    /// — batched so the MD5s see one `update` per full window, not one per byte.
    fn push_base(&mut self, b: u8) {
        self.upper.push(b.to_ascii_uppercase());
        if self.upper.len() == self.upper.capacity() {
            self.flush_md5();
        }
    }

    fn flush_md5(&mut self) {
        if !self.upper.is_empty() {
            self.contig_md5.update(&self.upper);
            self.reference_md5.update(&self.upper);
            self.upper.clear();
        }
    }

    /// A sequence line has ended (`\n`, or `is_last` at EOF). Fix or check the geometry: the
    /// first line sets it; every later line must match it, except the last, which may be
    /// shorter (`faidx.c` — "Different line length in sequence", spec §3.8/§4).
    fn complete_line(&mut self, is_last: bool) -> Result<(), ReferenceInfoError> {
        match self.line_bases0 {
            None => {
                self.line_bases0 = Some(self.cur_line_bases);
                self.line_width0 = self.cur_line_width;
            }
            Some(lb0) => {
                if self.cur_line_bases > lb0 {
                    return Err(self.malformed(&format!(
                        "line has {} bases, more than the first line's {lb0}",
                        self.cur_line_bases
                    )));
                }
                if !is_last {
                    if self.short_line_seen {
                        return Err(self.malformed(
                            "a line shorter than the first was not the last line of the contig",
                        ));
                    }
                    if self.cur_line_bases < lb0 {
                        self.short_line_seen = true;
                    } else if self.cur_line_width != self.line_width0 {
                        // A full interior line whose width changed = a mid-contig terminator
                        // change (e.g. LF→CR-LF); caught here, not silently reconstructed.
                        return Err(self.malformed(&format!(
                            "line width {} differs from the first line's {}",
                            self.cur_line_width, self.line_width0
                        )));
                    }
                }
            }
        }
        self.cur_line_bases = 0;
        self.cur_line_width = 0;
        Ok(())
    }

    /// Emit the current contig. Flushes its bases, rejects an empty contig (T3) and a
    /// duplicate name (T2), and records name/length/geometry/md5 in file order.
    fn finalize_contig(&mut self) -> Result<(), ReferenceInfoError> {
        self.flush_md5();
        if self.length == 0 {
            return Err(self.malformed("empty contig (no sequence bases)"));
        }
        // A contig with bases always has a completed line (the `\n` before the next `>`, or
        // the last-line completion at EOF), so this is an internal invariant — but the
        // module never panics on a supplied reference (spec §5 T3), so it surfaces as an
        // error rather than an `expect`.
        let Some(line_bases) = self.line_bases0 else {
            return Err(self.malformed("contig has bases but no completed sequence line"));
        };
        let digest: [u8; 16] = std::mem::replace(&mut self.contig_md5, Md5::new())
            .finalize()
            .into();
        let name = String::from_utf8_lossy(&self.name_bytes).into_owned();
        if !self.seen_names.insert(name.clone()) {
            return Err(ReferenceInfoError::DuplicateContigName {
                path: self.path.to_path_buf(),
                name,
            });
        }
        self.contigs.push(ContigInfo {
            name,
            length: self.length,
            offset: self.seq_offset,
            line_bases,
            line_width: self.line_width0,
            md5: Some(digest),
        });
        Ok(())
    }

    fn finish(mut self) -> Result<ReferenceInfo, ReferenceInfoError> {
        match self.mode {
            // An empty file (or whitespace only): no contigs, the digest of nothing.
            PassMode::Start => {}
            PassMode::Header => {
                return Err(self.malformed("contig header not followed by any sequence"));
            }
            PassMode::Sequence => {
                if self.cur_line_width > 0 {
                    // A last line with no trailing newline.
                    self.complete_line(true)?;
                }
                self.finalize_contig()?;
            }
        }
        let digest: [u8; 16] = self.reference_md5.finalize().into();
        Ok(ReferenceInfo {
            md5: Some(digest),
            contigs: self.contigs,
            // Filled in by the callers that know which FASTA was read.
            fasta_path: None,
        })
    }

    fn malformed(&self, detail: &str) -> ReferenceInfoError {
        ReferenceInfoError::MalformedFasta {
            path: self.path.to_path_buf(),
            contig: (!self.name_bytes.is_empty())
                .then(|| String::from_utf8_lossy(&self.name_bytes).into_owned()),
            byte_offset: self.pos,
            detail: detail.to_string(),
        }
    }
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

// ---------------------------------------------------------------------
// The `.fai` writer
// ---------------------------------------------------------------------

/// Write a `.fai` from contig info — the five-column `faidx.5` form (`md5` ignored, a `.fai`
/// has no MD5 column), **byte-identical to `samtools faidx`**: five TAB-separated columns,
/// one `\n`-terminated line per contig, no trailing blank line. The real use is indexing a
/// FASTA that lacks a `.fai` (spec §3.9).
///
/// Written **atomically** — a sibling temp file, then a rename over `out` — because a
/// half-written `.fai` is silently valid (no header, no checksum, nothing a reader would
/// reject), so every later run would trust a truncated index. The temp file is a sibling so
/// the rename stays on one filesystem (a cross-device rename is not atomic), and on any
/// failure it is removed, leaving `out` untouched.
pub fn write_fai(contigs: &[ContigInfo], out: &Path) -> io::Result<()> {
    use std::fmt::Write as _;
    let mut text = String::new();
    for c in contigs {
        // name \t length \t offset \t line_bases \t line_width — the five faidx.5 columns.
        // `writeln!` emits `\n` (never `\r\n`), matching samtools on every platform.
        writeln!(
            text,
            "{}\t{}\t{}\t{}\t{}",
            c.name, c.length, c.offset, c.line_bases, c.line_width
        )
        .expect("formatting into a String is infallible");
    }

    let dir = out
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let mut tmp = tempfile::NamedTempFile::new_in(dir)?;
    tmp.write_all(text.as_bytes())?;
    tmp.flush()?;
    tmp.persist(out).map_err(|e| e.error)?;
    Ok(())
}

// ---------------------------------------------------------------------
// The cache — compute-once, single-flight (spec §3.7)
// ---------------------------------------------------------------------

/// A file's identity for cache keying: `(path as given, size, mtime)`. Stat-based, which is
/// sound because a reference is write-once during a run (spec §5 T8). Paths are used as
/// given (not canonicalised), so two spellings of one file compute twice — harmless.
type FileStat = (PathBuf, u64, SystemTime);

/// The cache key. The **source discriminant is part of the key** because the three shapes
/// return *different* results for the same bytes — `Fai` yields `md5: None`; `FastaAlone`
/// the table with MD5s but unverified; `FastaWithFai` the same verified against that index —
/// so one must never satisfy another (spec §3.7).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum CacheKey {
    Fai(FileStat),
    FastaAlone(FileStat),
    FastaWithFai(FileStat, FileStat),
}

/// One key's single-flight slot: `None` until the read finishes, then the shared result. The
/// slot lock is held across the read (the single-flight), so same-key callers wait then hit.
type Slot = Arc<Mutex<Option<Arc<ReferenceInfo>>>>;

/// Stat one path into a `FileStat`, or `None` if it cannot be stat'd or its mtime is
/// unavailable — in which case the caller bypasses the cache rather than fail (spec §3.7).
fn file_stat(path: &Path) -> Option<FileStat> {
    let meta = std::fs::metadata(path).ok()?;
    let mtime = meta.modified().ok()?;
    Some((path.to_path_buf(), meta.len(), mtime))
}

/// Build the cache key for a source, or `None` when any involved file cannot be stat'd (a
/// missing file, or a platform without mtime) — the signal to bypass the cache.
fn cache_key(source: &ReferenceSource) -> Option<CacheKey> {
    Some(match source {
        ReferenceSource::Fai(p) => CacheKey::Fai(file_stat(p)?),
        ReferenceSource::Fasta { fasta, fai: None } => CacheKey::FastaAlone(file_stat(fasta)?),
        ReferenceSource::Fasta {
            fasta,
            fai: Some(f),
        } => CacheKey::FastaWithFai(file_stat(fasta)?, file_stat(f)?),
    })
}

/// A thread-safe, compute-once, single-flight cache over [`read_reference_info`]. Caller-held
/// (not a global static): construct one per run, share `&cache` across worker threads, and
/// every thread asking for a reference already being read **waits** rather than re-reading
/// it (spec §3.7). `Send + Sync`.
///
/// The consumer that makes it matter is parallel per-sample readers validating each sample's
/// alignment file against the *same* reference: N workers would otherwise each read the whole
/// genome to build the same table; the cache turns that into one read the others wait on.
#[derive(Default)]
pub struct ReferenceInfoCache {
    /// Two-level lock: the map lock is held only to get-or-create a key's slot; the slot lock
    /// is held across the read, which *is* the single-flight (same-key callers wait, then
    /// hit). A read on one key therefore never blocks a read on another.
    map: Mutex<HashMap<CacheKey, Slot>>,
}

impl ReferenceInfoCache {
    pub fn new() -> Self {
        Self::default()
    }

    /// Read a reference through the cache. The result is `Arc`-shared, so a hit is a pointer
    /// clone. Keyed on `(source discriminant, per-file (path, size, mtime))` — `Fai`,
    /// `Fasta { None }` and `Fasta { Some }` are distinct keys. **Successes only** are cached:
    /// a transient I/O error is not, so it can retry. When a file cannot be stat'd (missing,
    /// or no mtime), the cache is bypassed and the read runs uncached — the cache is an
    /// optimisation, never a correctness gate (spec §3.7).
    pub fn get_or_read(
        &self,
        source: ReferenceSource,
    ) -> Result<Arc<ReferenceInfo>, ReferenceInfoError> {
        let key = cache_key(&source);
        self.get_or_read_with(key, move || read_reference_info(source))
    }

    /// The single-flight core, parameterised over the actual read so tests can drive it with
    /// a counting / barrier-blocked read (the public `get_or_read` passes the real reader).
    fn get_or_read_with<F>(
        &self,
        key: Option<CacheKey>,
        read: F,
    ) -> Result<Arc<ReferenceInfo>, ReferenceInfoError>
    where
        F: FnOnce() -> Result<ReferenceInfo, ReferenceInfoError>,
    {
        // No key (un-stattable file / no mtime): bypass, compute directly, do not cache. The
        // read still produces the proper error for a missing file.
        let Some(key) = key else {
            return read().map(Arc::new);
        };

        // Hold the map lock only long enough to get-or-create this key's slot.
        let slot = {
            let mut map = self.map.lock().expect("reference-info cache map poisoned");
            Arc::clone(map.entry(key).or_default())
        };

        // Hold the slot lock across the read — this is the single-flight: a second caller on
        // the same key blocks here until the first finishes, then falls through to the hit.
        let mut guard = slot.lock().expect("reference-info cache slot poisoned");
        if let Some(existing) = guard.as_ref() {
            return Ok(Arc::clone(existing));
        }
        // Miss: run the read with the slot lock still held. Errors are NOT cached — the slot
        // stays `None`, so a transient failure can retry.
        let info = Arc::new(read()?);
        *guard = Some(Arc::clone(&info));
        Ok(info)
    }
}

// ---------------------------------------------------------------------
// Background verification (spec §3.10)
// ---------------------------------------------------------------------

const VERIFY_HANDLE_DROP_WARNING: &str = "warning: ReferenceInfo verification handle dropped without join(); a possibly-pending \
     reference-verification error went unobserved";

/// The pending result of a background FASTA verification (spec §3.10). **`#[must_use]`, and it
/// must be joined before the caller commits any output** — an unverified reference is a
/// possibly-wrong run, and this handle is the *only* place that error surfaces.
///
/// The verified `ReferenceInfo` (MD5s filled) arrives by **return** from [`Self::join`], never
/// by in-place mutation: `ContigInfo` stays a plain immutable struct, and immutable data plus
/// message-passing (the returned value *is* the message) avoids shared-mutable-across-threads.
#[must_use]
#[derive(Debug)]
pub struct VerificationHandle {
    /// `Some` until `join` takes it; a dedicated OS thread's handle whose `join` yields the
    /// verify's `Result` (a `JoinHandle` *is* the result channel — no separate channel needed;
    /// a dedicated thread, not a rayon worker, so a seconds-long read cannot starve the pool).
    handle: Option<std::thread::JoinHandle<Result<Arc<ReferenceInfo>, ReferenceInfoError>>>,
}

impl VerificationHandle {
    /// Non-blocking: `true` once verification has finished, so a caller can defer `join` until
    /// it will not block.
    pub fn is_finished(&self) -> bool {
        self.handle.as_ref().is_none_or(|h| h.is_finished())
    }

    /// Block until verification finishes and take its result — the **verified**
    /// `ReferenceInfo` (MD5s filled) on success, the geometry or read error on failure.
    pub fn join(mut self) -> Result<Arc<ReferenceInfo>, ReferenceInfoError> {
        let handle = self
            .handle
            .take()
            .expect("VerificationHandle already joined");
        match handle.join() {
            Ok(result) => result,
            // The verify thread panicked (read_reference_info does not, so this is a bug, not
            // bad input): re-raise rather than swallow it.
            Err(panic) => std::panic::resume_unwind(panic),
        }
    }
}

impl Drop for VerificationHandle {
    fn drop(&mut self) {
        // Abandoning the run's only verification hides a possibly-pending failure — warn
        // rather than swallow it (but stay quiet if we are already unwinding a panic).
        if self.handle.is_some() && !std::thread::panicking() {
            eprintln!("{VERIFY_HANDLE_DROP_WARNING}");
        }
    }
}

/// Read the `.fai` **now** (fast, `md5: None`, returned immediately) and verify the FASTA
/// against it on a **background thread** — both reads going **through the cache**, so the
/// verify populates it and a concurrent reader of the same key coordinates on the
/// single-flight slot instead of reading again (spec §3.10). Takes `&Arc<ReferenceInfoCache>`
/// because the background thread outlives the call and holds a clone.
///
/// The initial `.fai` read is synchronous, so a missing or malformed `.fai` errors here,
/// before any thread is spawned.
pub fn read_fai_verify_in_background(
    cache: &Arc<ReferenceInfoCache>,
    fasta: PathBuf,
    fai: PathBuf,
) -> Result<(Arc<ReferenceInfo>, VerificationHandle), ReferenceInfoError> {
    // Foreground: the cheap `.fai` read, available immediately (names/lengths/geometry, no
    // digests). A missing or malformed `.fai` surfaces synchronously.
    let cached = cache.get_or_read(ReferenceSource::Fai(fai.clone()))?;

    // The cached value is keyed on `Fai(..)` and so carries no FASTA path — correctly, since
    // a caller who asked for a `.fai` alone has no FASTA. *This* caller does, and a consumer
    // that needs the bases (CRAM decoding) would otherwise be told to "supply the reference
    // FASTA" when it just had. Fill the field in on a clone rather than in the cache, so the
    // pure-`.fai` entry stays truthful for whoever shares it.
    let info = if cached.fasta_path.is_none() {
        Arc::new(ReferenceInfo {
            fasta_path: Some(fasta.clone()),
            md5: cached.md5,
            contigs: cached.contigs.clone(),
        })
    } else {
        cached
    };

    // Background: verify the FASTA against the `.fai` through the cache (§3.3 + §3.4), so the
    // verified result — and its error — reach the caller via `join`, and the cache is
    // populated for any later or concurrent reader of the same key.
    let cache = Arc::clone(cache);
    let handle = std::thread::spawn(move || {
        cache.get_or_read(ReferenceSource::Fasta {
            fasta,
            fai: Some(fai),
        })
    });

    Ok((
        info,
        VerificationHandle {
            handle: Some(handle),
        },
    ))
}

// ---------------------------------------------------------------------
// The batteries-included orchestrator (spec §3.11)
// ---------------------------------------------------------------------

/// Read a reference the convenient way: derive the sibling `<fasta>.fai` and do the right
/// thing (spec §3.11).
/// - **`.fai` present** → verify in the background ([`read_fai_verify_in_background`]);
///   returns immediately with `(info, Some(handle))` the caller must `join` before committing
///   output.
/// - **`.fai` absent** → scan the FASTA now (verified, MD5s), **write** the `.fai` beside it,
///   and return `(verified info, None)`.
///
/// This is the *one* place the module adopts the `<fasta>.fai` convention (an exception to the
/// §3.6 modularity line) and the *one* place a read can write a file. A `.fai`-write failure is
/// **fatal** (`FaiWrite`): a caller on a read-only reference dir uses the primitives instead
/// (`cache.get_or_read(Fasta { fai: None })`, which reads and caches but writes nothing).
///
/// The return is asymmetric on purpose: `Some(handle)` = verification pending (join to
/// verify); `None` = already verified, nothing to await. A caller wanting uniformly-verified
/// info writes `if let Some(h) = handle { h.join()?; }`.
pub fn read_reference_verifying_or_creating_fai(
    cache: &Arc<ReferenceInfoCache>,
    fasta: PathBuf,
) -> Result<(Arc<ReferenceInfo>, Option<VerificationHandle>), ReferenceInfoError> {
    let fai = sibling_fai_path(&fasta);
    if fai.exists() {
        let (info, handle) = read_fai_verify_in_background(cache, fasta, fai)?;
        Ok((info, Some(handle)))
    } else {
        // No `.fai`: scan now (verified, with the MD5s), then index the reference ourselves.
        let info = cache.get_or_read(ReferenceSource::Fasta { fasta, fai: None })?;
        write_fai(&info.contigs, &fai)
            .map_err(|source| ReferenceInfoError::FaiWrite { path: fai, source })?;
        Ok((info, None))
    }
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
            fasta_path: None,
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

    // =================================================================
    // A4 — read_reference_info, the `Fai` arm (spec §4 Fai arm, §3.8,
    //      §5 T2/T3; arch §5).
    // =================================================================

    /// The FASTQ index from the `faidx.5` worked example — six columns (spec §3.8).
    const WORKED_FASTQ_FAI: &str = "fastq1\t66\t8\t30\t31\t79\nfastq2\t28\t156\t14\t15\t188\n";

    /// Write raw `.fai` text into a fresh tempdir; returns `(dir, fai_path)`.
    fn write_fai_text(text: &str) -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("ref.fa.fai");
        std::fs::write(&path, text).expect("write .fai");
        (dir, path)
    }

    #[test]
    fn fai_arm_parses_the_committed_golden_fai() {
        let info = read_reference_info(ReferenceSource::Fai(sibling_fai_path(&golden_ref_path())))
            .expect("golden .fai parses");

        // A `.fai`-only read carries no digests.
        assert_eq!(info.md5, None);
        assert_eq!(info.contigs.len(), GOLDEN.len());
        for (got, g) in info.contigs.iter().zip(GOLDEN.iter()) {
            assert_eq!(got.name, g.name);
            assert_eq!(got.length, g.length);
            assert_eq!(got.offset, g.offset);
            assert_eq!(got.line_bases, g.line_bases);
            assert_eq!(got.line_width, g.line_width);
            assert_eq!(got.md5, None, "the Fai arm has no MD5 column to read");
        }
    }

    #[test]
    fn fai_arm_rejects_a_duplicate_contig_name() {
        let (_dir, path) = write_fai_text("dup\t10\t5\t10\t11\ndup\t10\t22\t10\t11\n");
        match read_reference_info(ReferenceSource::Fai(path)).unwrap_err() {
            ReferenceInfoError::DuplicateContigName { name, .. } => assert_eq!(name, "dup"),
            other => panic!("expected DuplicateContigName, got {other:?}"),
        }
    }

    #[test]
    fn fai_arm_rejects_a_fastq_index() {
        let (_dir, path) = write_fai_text(WORKED_FASTQ_FAI);
        assert!(
            matches!(
                read_reference_info(ReferenceSource::Fai(path)).unwrap_err(),
                ReferenceInfoError::FastqIndex { .. }
            ),
            "a six-column FASTQ index must be named, not misparsed",
        );
    }

    #[test]
    fn fai_arm_rejects_line_bases_zero() {
        // `line_bases == 0` divides by zero downstream; the guard catches it (T3).
        let (_dir, path) = write_fai_text("c\t10\t3\t0\t1\n");
        match read_reference_info(ReferenceSource::Fai(path)).unwrap_err() {
            ReferenceInfoError::FaiRead { source, .. } => {
                assert!(
                    source.to_string().contains("line_bases = 0"),
                    "source detail: {source}"
                );
            }
            other => panic!("expected FaiRead, got {other:?}"),
        }
    }

    #[test]
    fn fai_arm_rejects_line_width_below_line_bases() {
        // `line_width` must include the terminator, so it cannot be < `line_bases` (T3).
        let (_dir, path) = write_fai_text("c\t10\t3\t10\t5\n");
        match read_reference_info(ReferenceSource::Fai(path)).unwrap_err() {
            ReferenceInfoError::FaiRead { source, .. } => {
                assert!(
                    source
                        .to_string()
                        .contains("line_width (5) < line_bases (10)"),
                    "source detail: {source}"
                );
            }
            other => panic!("expected FaiRead, got {other:?}"),
        }
    }

    #[test]
    fn fai_arm_errors_on_a_missing_index() {
        // T4: a `.fai` that does not open is an error, not a silent fall-back.
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("nope.fa.fai");
        assert!(matches!(
            read_reference_info(ReferenceSource::Fai(missing)).unwrap_err(),
            ReferenceInfoError::FaiRead { .. }
        ));
    }

    // =================================================================
    // B (B1+B2) — the FASTA streaming pass (`Fasta { fai: None }`), and
    //      its samtools/`.cat`/`faidx.5` oracle (spec §4, §3.4, §3.8).
    // =================================================================

    fn hex_to_md5(hex: &str) -> [u8; 16] {
        let mut out = [0u8; 16];
        for (i, byte) in out.iter_mut().enumerate() {
            *byte = u8::from_str_radix(&hex[i * 2..i * 2 + 2], 16).expect("hex digit pair");
        }
        out
    }

    fn write_bytes_fasta(bytes: &[u8]) -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("ref.fa");
        std::fs::write(&path, bytes).expect("write fasta");
        (dir, path)
    }

    fn read_fasta_none(path: PathBuf) -> Result<ReferenceInfo, ReferenceInfoError> {
        read_reference_info(ReferenceSource::Fasta {
            fasta: path,
            fai: None,
        })
    }

    #[test]
    fn fasta_pass_reconstructs_the_golden_index_and_per_contig_md5() {
        let info = read_fasta_none(golden_ref_path()).expect("golden FASTA parses");
        assert_eq!(info.contigs.len(), GOLDEN.len());
        for (got, g) in info.contigs.iter().zip(GOLDEN.iter()) {
            assert_eq!(got.name, g.name);
            assert_eq!(got.length, g.length, "LN for {}", g.name);
            assert_eq!(got.offset, g.offset, "offset for {}", g.name);
            assert_eq!(got.line_bases, g.line_bases);
            assert_eq!(got.line_width, g.line_width);
            // The per-contig MD5 equals `samtools dict`'s @SQ M5 (spec §3.4).
            assert_eq!(got.md5, Some(hex_to_md5(g.md5_hex)), "M5 for {}", g.name);
        }
    }

    #[test]
    fn fasta_pass_reference_md5_matches_the_golden_cat_header() {
        use crate::ssr::catalog::io::CatalogReader;
        let info = read_fasta_none(golden_ref_path()).unwrap();
        let cat_path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/tandem_repeat/golden.ssr_catalog.bed.gz");
        let reader = CatalogReader::new(std::fs::File::open(cat_path).unwrap()).unwrap();
        let expected = reader.header().reference_md5.clone();
        let got = crate::pop_var_caller::common::format_md5_hex(info.md5.unwrap());
        assert_eq!(
            got, expected,
            "whole-reference digest vs the golden .cat header"
        );
    }

    #[test]
    fn fasta_pass_reproduces_the_faidx5_worked_example() {
        // The strongest geometry check: reconstruct the man-page vector byte-for-byte, under
        // both line endings (spec §3.8). CR-LF falls out of `line_width - line_bases`.
        for (crlf, index) in [(false, WORKED_INDEX_LF), (true, WORKED_INDEX_CRLF)] {
            let (_dir, path) = write_bytes_fasta(&build_worked_example(crlf));
            let info = read_fasta_none(path).unwrap();
            assert_eq!(info.contigs.len(), 2, "crlf={crlf}");
            for (got, (name, length, offset, lb, lw)) in info.contigs.iter().zip(index) {
                assert_eq!(got.name, name, "crlf={crlf}");
                assert_eq!(
                    (got.length, got.offset, got.line_bases, got.line_width),
                    (length, offset, lb, lw),
                    "contig {name} crlf={crlf}"
                );
            }
        }
    }

    #[test]
    fn fasta_pass_per_contig_md5_matches_one_shot() {
        // Streaming MD5 == the one-shot `Md5::digest` of the filtered, uppercased bases —
        // production's own streaming-matches-one-shot check, over the golden reference.
        let info = read_fasta_none(golden_ref_path()).unwrap();
        let file = std::fs::File::open(golden_ref_path()).unwrap();
        let mut reader = noodles_fasta::io::Reader::new(std::io::BufReader::new(file));
        for (rec, got) in reader.records().zip(info.contigs.iter()) {
            let rec = rec.unwrap();
            let upper: Vec<u8> = rec
                .sequence()
                .as_ref()
                .iter()
                .filter(|&&b| (0x21..=0x7e).contains(&b))
                .map(|b| b.to_ascii_uppercase())
                .collect();
            let one_shot: [u8; 16] = Md5::digest(&upper).into();
            assert_eq!(
                got.md5,
                Some(one_shot),
                "one-shot vs streaming for {}",
                got.name
            );
        }
    }

    #[test]
    fn fasta_pass_space_and_tab_hash_as_absent() {
        // The predicate edge that distinguishes our rule from production's (spec §3.4/§6):
        // a sequence line carrying a space and a tab hashes as if they were not there, and
        // they count toward neither `line_bases` nor `length`.
        let (_d1, with_ws) = write_bytes_fasta(b">sp\nACG T\tACG\n");
        let (_d2, control) = write_bytes_fasta(b">ctl\nACGTACG\n");
        let a = read_fasta_none(with_ws).unwrap();
        let b = read_fasta_none(control).unwrap();
        assert_eq!(a.contigs[0].length, 7, "space and tab are not bases");
        assert_eq!(a.contigs[0].line_bases, 7);
        assert_eq!(
            a.contigs[0].md5, b.contigs[0].md5,
            "the space/tab line hashes identically to the clean one"
        );
    }

    #[test]
    fn fasta_pass_rejects_bases_before_the_first_header() {
        let (_dir, path) = write_bytes_fasta(b"ACGT\n>c\nACGT\n");
        match read_fasta_none(path).unwrap_err() {
            ReferenceInfoError::MalformedFasta { detail, contig, .. } => {
                assert!(detail.contains("before the first"), "detail: {detail}");
                assert_eq!(contig, None, "no contig header seen yet");
            }
            other => panic!("expected MalformedFasta, got {other:?}"),
        }
    }

    #[test]
    fn fasta_pass_rejects_a_non_uniform_interior_line() {
        // A short interior line followed by another line — the short one was not last (T3).
        let (_dir, path) = write_bytes_fasta(b">c\nACGTAC\nACG\nACGTAC\n");
        assert!(matches!(
            read_fasta_none(path).unwrap_err(),
            ReferenceInfoError::MalformedFasta { .. }
        ));
    }

    #[test]
    fn fasta_pass_rejects_an_empty_contig() {
        let (_dir, path) = write_bytes_fasta(b">a\n>b\nACGT\n");
        match read_fasta_none(path).unwrap_err() {
            ReferenceInfoError::MalformedFasta { detail, contig, .. } => {
                assert!(detail.contains("empty contig"), "detail: {detail}");
                assert_eq!(contig.as_deref(), Some("a"));
            }
            other => panic!("expected MalformedFasta, got {other:?}"),
        }
    }

    #[test]
    fn fasta_pass_rejects_a_compressed_reference() {
        // gzip magic → a clean CompressedReference, not compressed bytes misread as FASTA.
        let (_dir, path) = write_bytes_fasta(&[0x1f, 0x8b, 0x08, 0x00, 0x00]);
        assert!(matches!(
            read_fasta_none(path).unwrap_err(),
            ReferenceInfoError::CompressedReference { .. }
        ));
    }

    #[test]
    fn fasta_pass_rejects_a_duplicate_contig_name() {
        let (_dir, path) = write_bytes_fasta(b">dup\nACGT\n>dup\nTTTT\n");
        match read_fasta_none(path).unwrap_err() {
            ReferenceInfoError::DuplicateContigName { name, .. } => assert_eq!(name, "dup"),
            other => panic!("expected DuplicateContigName, got {other:?}"),
        }
    }

    #[test]
    fn fasta_pass_rejects_an_empty_contig_name() {
        // T3: a definition line with no name (here `>` then a space) is malformed.
        let (_dir, path) = write_bytes_fasta(b"> desc only\nACGT\n");
        match read_fasta_none(path).unwrap_err() {
            ReferenceInfoError::MalformedFasta { detail, .. } => {
                assert!(detail.contains("empty contig name"), "detail: {detail}");
            }
            other => panic!("expected MalformedFasta, got {other:?}"),
        }
    }

    #[test]
    fn fasta_pass_handles_a_last_line_without_a_trailing_newline() {
        // A full last line at EOF (no `\n`) must still index correctly (length/geometry).
        let (_dir, path) = write_bytes_fasta(b">c\nACGT\nACGT\nAC");
        let info = read_fasta_none(path).unwrap();
        assert_eq!(info.contigs.len(), 1);
        assert_eq!(info.contigs[0].length, 10);
        assert_eq!(info.contigs[0].line_bases, 4);
        assert_eq!(info.contigs[0].line_width, 5);
    }

    // =================================================================
    // B3 — the fasta-vs-`.fai` check (`Fasta { fai: Some }`, spec §3.3,
    //      §5 T1).
    // =================================================================

    /// `n` aperiodic-ish bases (`ACGT…` cycled) — content is irrelevant to geometry.
    fn bases(n: usize) -> Vec<u8> {
        (0..n).map(|i| b"ACGT"[i % 4]).collect()
    }

    fn read_fasta_with_fai(
        fasta: PathBuf,
        fai: PathBuf,
    ) -> Result<ReferenceInfo, ReferenceInfoError> {
        read_reference_info(ReferenceSource::Fasta {
            fasta,
            fai: Some(fai),
        })
    }

    #[test]
    fn fasta_fai_check_passes_on_a_matching_index() {
        // The `.fai` `write_wrapped_fasta` emits describes the FASTA it wrote; the pass must
        // reconstruct the same five columns and agree.
        let c1 = bases(150);
        let c2 = bases(80);
        let (_dir, fasta, fai) = write_wrapped_fasta(&[("c1", &c1), ("c2", &c2)], 60);
        let info = read_fasta_with_fai(fasta, fai).expect("matching .fai verifies clean");
        assert_eq!(info.contigs.len(), 2);
        assert!(info.md5.is_some()); // the FASTA read still fills the digests
    }

    #[test]
    fn fasta_fai_check_catches_a_rewrap_naming_line_bases() {
        // The spec's counterexample (§3.3): re-wrapping keeps names and lengths but breaks
        // the geometry. Same contig, indexed at width 40, FASTA written at width 60.
        let c = bases(150);
        let (dir60, fasta60, _fai60) = write_wrapped_fasta(&[("c", &c)], 60);
        let (_dir40, _fasta40, fai40) = write_wrapped_fasta(&[("c", &c)], 40);

        // Mutation-verify against a names-only check: names AND lengths match, so a cheaper
        // check would pass — only the geometry gives it away.
        let from_fasta = &read_fasta_none(fasta60.clone()).unwrap().contigs[0];
        let from_fai = &read_fai(&fai40).unwrap().contigs[0];
        assert_eq!(from_fasta.name, from_fai.name, "a names check would pass");
        assert_eq!(
            from_fasta.length, from_fai.length,
            "a lengths check would pass"
        );
        assert_ne!(
            from_fasta.line_bases, from_fai.line_bases,
            "only geometry differs"
        );

        match read_fasta_with_fai(fasta60, fai40).unwrap_err() {
            ReferenceInfoError::FastaFaiMismatch { field, contig, .. } => {
                assert_eq!(
                    field, "line_bases",
                    "the re-wrap must be named at line_bases"
                );
                assert_eq!(contig, "c");
            }
            other => panic!("expected FastaFaiMismatch, got {other:?}"),
        }
        drop(dir60);
    }

    #[test]
    fn fasta_fai_check_catches_a_single_contig_rewrap_where_offset_is_unchanged() {
        // On a single contig the header is unchanged, so OFFSET matches even after a re-wrap
        // (§3.3) — the case an offsets-only check would miss. `line_bases` still catches it.
        let c = bases(200);
        let (_d60, fasta60, _f60) = write_wrapped_fasta(&[("solo", &c)], 60);
        let (_d50, _f50fa, fai50) = write_wrapped_fasta(&[("solo", &c)], 50);

        let from_fasta = &read_fasta_none(fasta60.clone()).unwrap().contigs[0];
        let from_fai = &read_fai(&fai50).unwrap().contigs[0];
        assert_eq!(
            from_fasta.offset, from_fai.offset,
            "single-contig offset is unchanged"
        );

        match read_fasta_with_fai(fasta60, fai50).unwrap_err() {
            ReferenceInfoError::FastaFaiMismatch { field, .. } => assert_eq!(field, "line_bases"),
            other => panic!("expected FastaFaiMismatch, got {other:?}"),
        }
    }

    #[test]
    fn fasta_fai_check_catches_a_reordering_and_the_digest_differs() {
        // A permuted FASTA against its original `.fai` errors on `name` (position-wise), and
        // the reference digest is order-dependent (T1).
        let c1 = bases(120);
        let c2 = bases(90);
        let (_da, fasta_a, fai_a) = write_wrapped_fasta(&[("c1", &c1), ("c2", &c2)], 60);
        let (_db, fasta_b, _fai_b) = write_wrapped_fasta(&[("c2", &c2), ("c1", &c1)], 60);

        match read_fasta_with_fai(fasta_b.clone(), fai_a).unwrap_err() {
            ReferenceInfoError::FastaFaiMismatch { field, contig, .. } => {
                assert_eq!(field, "name");
                assert_eq!(contig, "c2", "the FASTA's first contig after the swap");
            }
            other => panic!("expected FastaFaiMismatch, got {other:?}"),
        }

        // The digest pins the order all by itself.
        let md5_a = read_fasta_none(fasta_a).unwrap().md5;
        let md5_b = read_fasta_none(fasta_b).unwrap().md5;
        assert_ne!(
            md5_a, md5_b,
            "a permuted reference has a different digest (T1)"
        );
    }

    #[test]
    fn fasta_fai_check_catches_a_contig_count_mismatch() {
        let c1 = bases(120);
        let c2 = bases(90);
        let (_dir, fasta2, _fai2) = write_wrapped_fasta(&[("c1", &c1), ("c2", &c2)], 60);
        let (_dir1, _fasta1, fai1) = write_wrapped_fasta(&[("c1", &c1)], 60);
        match read_fasta_with_fai(fasta2, fai1).unwrap_err() {
            ReferenceInfoError::FastaFaiMismatch { field, .. } => assert_eq!(field, "contig count"),
            other => panic!("expected FastaFaiMismatch, got {other:?}"),
        }
    }

    #[test]
    fn fasta_fai_check_errors_when_the_supplied_fai_is_missing() {
        // T4: a supplied `.fai` that does not open errors — it does not degrade to reading
        // the FASTA alone.
        let c = bases(80);
        let (dir, fasta, fai) = write_wrapped_fasta(&[("c", &c)], 60);
        std::fs::remove_file(&fai).unwrap();
        assert!(matches!(
            read_fasta_with_fai(fasta, fai).unwrap_err(),
            ReferenceInfoError::FaiRead { .. }
        ));
        drop(dir);
    }

    // =================================================================
    // C1 — write_fai (spec §3.9): byte-identical to samtools, atomic.
    // =================================================================

    #[test]
    fn write_fai_is_byte_identical_to_samtools() {
        // The oracle: a `.fai` we write is indistinguishable from the committed
        // `samtools faidx` one for the golden reference.
        let info = read_fasta_none(golden_ref_path()).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.fai");
        write_fai(&info.contigs, &out).unwrap();
        let ours = std::fs::read(&out).unwrap();
        let samtools = std::fs::read(sibling_fai_path(&golden_ref_path())).unwrap();
        assert_eq!(
            ours, samtools,
            "write_fai must equal samtools faidx byte-for-byte"
        );
    }

    #[test]
    fn write_fai_round_trip_is_a_fixpoint() {
        // read → write → read: the `.fai` we write verifies clean against its own FASTA.
        let c1 = bases(150);
        let c2 = bases(85);
        let (_dir, fasta, _fai) = write_wrapped_fasta(&[("c1", &c1), ("c2", &c2)], 60);
        let info = read_fasta_none(fasta.clone()).unwrap();
        let out_dir = tempfile::tempdir().unwrap();
        let out = out_dir.path().join("ref.fa.fai");
        write_fai(&info.contigs, &out).unwrap();
        read_fasta_with_fai(fasta, out).expect("write_fai output verifies against its FASTA");
    }

    #[test]
    fn write_fai_is_atomic_no_partial_file_on_failure() {
        // The atomic writer only ever creates `out` by renaming a fully-written temp file,
        // so a failure at any point leaves no partial `.fai`. Force the rename to fail by
        // making `out` an existing directory (renaming a file over a directory is rejected
        // even for root, so this holds in the container too).
        let info = read_fasta_none(golden_ref_path()).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("occupied");
        std::fs::create_dir(&out).unwrap();

        let result = write_fai(&info.contigs, &out);
        assert!(result.is_err(), "renaming over a directory must fail");
        assert!(
            out.is_dir(),
            "out is untouched — no partial .fai replaced it"
        );
        // The temp file was cleaned up on failure — only the directory remains.
        let entries = std::fs::read_dir(dir.path()).unwrap().count();
        assert_eq!(entries, 1, "no leftover temp file after a failed write");
    }

    // =================================================================
    // D1 — ReferenceInfoCache: compute-once, single-flight (spec §3.7).
    // =================================================================

    use std::sync::atomic::{AtomicUsize, Ordering};

    fn dummy_info() -> ReferenceInfo {
        ReferenceInfo {
            md5: None,
            contigs: Vec::new(),
            fasta_path: None,
        }
    }

    /// A fixed cache key for the seam tests (no real file needed).
    fn test_key(tag: &str) -> CacheKey {
        CacheKey::Fai((PathBuf::from(tag), 0, SystemTime::UNIX_EPOCH))
    }

    #[test]
    fn cache_hit_computes_the_read_once() {
        let cache = ReferenceInfoCache::new();
        let key = test_key("ref");
        let n = AtomicUsize::new(0);
        let a = cache
            .get_or_read_with(Some(key.clone()), || {
                n.fetch_add(1, Ordering::SeqCst);
                Ok(dummy_info())
            })
            .unwrap();
        let b = cache
            .get_or_read_with(Some(key), || {
                n.fetch_add(1, Ordering::SeqCst);
                Ok(dummy_info())
            })
            .unwrap();
        assert_eq!(n.load(Ordering::SeqCst), 1, "the second call is a hit");
        assert!(Arc::ptr_eq(&a, &b), "a hit returns the same Arc");
    }

    #[test]
    fn cache_single_flight_reads_once_under_threads() {
        let cache = Arc::new(ReferenceInfoCache::new());
        let key = test_key("shared");
        let n = Arc::new(AtomicUsize::new(0));
        let handles: Vec<_> = (0..8)
            .map(|_| {
                let cache = Arc::clone(&cache);
                let key = key.clone();
                let n = Arc::clone(&n);
                std::thread::spawn(move || {
                    cache
                        .get_or_read_with(Some(key), || {
                            n.fetch_add(1, Ordering::SeqCst);
                            // Hold the slot long enough that the other threads pile up on it.
                            std::thread::sleep(std::time::Duration::from_millis(100));
                            Ok(dummy_info())
                        })
                        .unwrap()
                })
            })
            .collect();
        let results: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();
        assert_eq!(
            n.load(Ordering::SeqCst),
            1,
            "the real read ran exactly once"
        );
        assert!(
            results.iter().all(|r| Arc::ptr_eq(r, &results[0])),
            "every thread got the same Arc"
        );
    }

    #[test]
    fn cache_different_keys_do_not_serialize() {
        // The map lock is not held across a read: a read blocked on one key must not wedge a
        // read on another. Thread A's read blocks until released; thread B (a different key)
        // must complete *before* A is released.
        let cache = Arc::new(ReferenceInfoCache::new());
        let (release, blocked_read) = std::sync::mpsc::channel::<()>();

        let a = {
            let cache = Arc::clone(&cache);
            std::thread::spawn(move || {
                cache
                    .get_or_read_with(Some(test_key("a")), move || {
                        blocked_read.recv().unwrap(); // block until the test releases us
                        Ok(dummy_info())
                    })
                    .unwrap()
            })
        };
        let b = {
            let cache = Arc::clone(&cache);
            std::thread::spawn(move || {
                cache
                    .get_or_read_with(Some(test_key("b")), || Ok(dummy_info()))
                    .unwrap()
            })
        };

        b.join().unwrap(); // completes while A is still blocked — no serialization
        release.send(()).unwrap();
        a.join().unwrap();
    }

    #[test]
    fn cache_distinguishes_fai_from_fasta_of_the_same_file() {
        // The source discriminant is in the key: a cheap `Fai` read must not poison the entry
        // a `Fasta` read needs (T6-adjacent).
        let cache = ReferenceInfoCache::new();
        let from_fai = cache
            .get_or_read(ReferenceSource::Fai(sibling_fai_path(&golden_ref_path())))
            .unwrap();
        let from_fasta = cache
            .get_or_read(ReferenceSource::Fasta {
                fasta: golden_ref_path(),
                fai: None,
            })
            .unwrap();
        assert_eq!(from_fai.md5, None, "the Fai read has no digest");
        assert!(from_fasta.md5.is_some(), "the Fasta read has the digest");
        assert!(
            !Arc::ptr_eq(&from_fai, &from_fasta),
            "distinct keys, distinct entries"
        );
    }

    #[test]
    fn cache_does_not_cache_errors() {
        let cache = ReferenceInfoCache::new();
        let key = test_key("flaky");
        let first = cache.get_or_read_with(Some(key.clone()), || {
            Err(ReferenceInfoError::FastaRead {
                path: PathBuf::from("x"),
                source: io::Error::other("transient"),
            })
        });
        assert!(first.is_err());

        // The error was not cached, so the slot is retried and this read actually runs.
        let n = AtomicUsize::new(0);
        cache
            .get_or_read_with(Some(key), || {
                n.fetch_add(1, Ordering::SeqCst);
                Ok(dummy_info())
            })
            .unwrap();
        assert_eq!(n.load(Ordering::SeqCst), 1, "the retry ran the read");
    }

    #[test]
    fn cache_bypasses_when_there_is_no_key() {
        // No key (un-stattable file / no mtime) => not cached => every call reads.
        let cache = ReferenceInfoCache::new();
        let n = AtomicUsize::new(0);
        let a = cache
            .get_or_read_with(None, || {
                n.fetch_add(1, Ordering::SeqCst);
                Ok(dummy_info())
            })
            .unwrap();
        let b = cache
            .get_or_read_with(None, || {
                n.fetch_add(1, Ordering::SeqCst);
                Ok(dummy_info())
            })
            .unwrap();
        assert_eq!(n.load(Ordering::SeqCst), 2, "no key => both calls read");
        assert!(!Arc::ptr_eq(&a, &b), "the bypass does not share an Arc");

        // A real missing file has no key, so it bypasses to the proper read error.
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("nope.fa.fai");
        assert!(matches!(
            cache
                .get_or_read(ReferenceSource::Fai(missing))
                .unwrap_err(),
            ReferenceInfoError::FaiRead { .. }
        ));
    }

    #[test]
    fn cache_is_send_and_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<ReferenceInfoCache>();
    }

    // =================================================================
    // E1 — read_fai_verify_in_background (spec §3.10).
    // =================================================================

    /// A FASTA + matching `.fai` written into a fresh tempdir, both at wrap `w`.
    fn golden_like_pair(w: usize) -> (tempfile::TempDir, PathBuf, PathBuf) {
        let c1 = bases(150);
        let c2 = bases(95);
        write_wrapped_fasta(&[("c1", &c1), ("c2", &c2)], w)
    }

    #[test]
    fn background_verify_info_available_before_join_and_join_upgrades() {
        let cache = Arc::new(ReferenceInfoCache::new());
        let (_dir, fasta, fai) = golden_like_pair(60);
        let (info, handle) = read_fai_verify_in_background(&cache, fasta, fai).unwrap();
        // The foreground `.fai` read has no digests.
        assert_eq!(info.md5, None);
        assert_eq!(info.contigs.len(), 2);
        // `join` upgrades to the verified info, with the digests filled.
        let verified = handle.join().unwrap();
        assert!(verified.md5.is_some(), "join fills the reference digest");
        assert!(verified.contigs.iter().all(|c| c.md5.is_some()));
    }

    #[test]
    fn background_verify_stale_fai_errors_on_join_not_before() {
        // A `.fai` from a different wrap of the same contigs: names/lengths match (foreground
        // read succeeds) but the geometry is stale (the verify fails on join).
        let cache = Arc::new(ReferenceInfoCache::new());
        let (_dir60, fasta60, _fai60) = golden_like_pair(60);
        let (_dir40, _fasta40, fai40) = golden_like_pair(40);

        let (info, handle) = read_fai_verify_in_background(&cache, fasta60, fai40).unwrap();
        // The foreground read succeeded — the error is NOT surfaced before join.
        assert_eq!(info.md5, None);
        // Mutation-verify: a swallowed error would make this `Ok` and fail the test.
        match handle.join() {
            Err(ReferenceInfoError::FastaFaiMismatch { field, .. }) => {
                assert_eq!(field, "line_bases")
            }
            other => panic!("expected a stale-.fai FastaFaiMismatch on join, got {other:?}"),
        }
    }

    #[test]
    fn background_verify_populates_the_cache() {
        let cache = Arc::new(ReferenceInfoCache::new());
        let (_dir, fasta, fai) = golden_like_pair(60);
        let (_info, handle) =
            read_fai_verify_in_background(&cache, fasta.clone(), fai.clone()).unwrap();
        let verified = handle.join().unwrap();
        // A later read of the same `Fasta { Some }` key is a hit — the verify populated it.
        let again = cache
            .get_or_read(ReferenceSource::Fasta {
                fasta,
                fai: Some(fai),
            })
            .unwrap();
        assert!(
            Arc::ptr_eq(&verified, &again),
            "the background verify populated the cache"
        );
    }

    #[test]
    fn background_verify_concurrent_reader_coordinates_through_the_cache() {
        // A reader of the same `Fasta { Some }` key, started while the verify runs, coordinates
        // on the single-flight slot and gets the *same* Arc — not a second genome read.
        let cache = Arc::new(ReferenceInfoCache::new());
        let (_dir, fasta, fai) = golden_like_pair(60);
        let (_info, handle) =
            read_fai_verify_in_background(&cache, fasta.clone(), fai.clone()).unwrap();
        let cache2 = Arc::clone(&cache);
        let other = std::thread::spawn(move || {
            cache2
                .get_or_read(ReferenceSource::Fasta {
                    fasta,
                    fai: Some(fai),
                })
                .unwrap()
        });
        let verified = handle.join().unwrap();
        let concurrent = other.join().unwrap();
        assert!(
            Arc::ptr_eq(&verified, &concurrent),
            "the concurrent reader shared the verify's Arc through the cache"
        );
    }

    #[test]
    fn background_verify_failure_does_not_poison_the_key() {
        let cache = Arc::new(ReferenceInfoCache::new());
        let (_dir60, fasta60, _f60) = golden_like_pair(60);
        let (_dir40, _f40, fai40) = golden_like_pair(40);
        let (_info, handle) =
            read_fai_verify_in_background(&cache, fasta60.clone(), fai40.clone()).unwrap();
        assert!(handle.join().is_err());
        // The failed verify was not cached, so a later read re-runs (and re-fails).
        assert!(matches!(
            cache
                .get_or_read(ReferenceSource::Fasta {
                    fasta: fasta60,
                    fai: Some(fai40),
                })
                .unwrap_err(),
            ReferenceInfoError::FastaFaiMismatch { .. }
        ));
    }

    #[test]
    fn background_verify_missing_fai_errors_before_spawning() {
        // The foreground `.fai` read is synchronous: a missing `.fai` errors here, not on join.
        let cache = Arc::new(ReferenceInfoCache::new());
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("ref.fa");
        std::fs::write(&fasta, b">c\nACGT\n").unwrap();
        let fai = dir.path().join("ref.fa.fai");
        assert!(matches!(
            read_fai_verify_in_background(&cache, fasta, fai).unwrap_err(),
            ReferenceInfoError::FaiRead { .. }
        ));
    }

    #[test]
    fn abandoning_the_handle_does_not_panic_and_the_warning_is_pinned() {
        // The Drop warning is the guard against a lost error; pin its text, and prove that
        // dropping an un-joined handle runs cleanly (it warns to stderr — visible under
        // --nocapture — rather than swallowing silently).
        assert!(VERIFY_HANDLE_DROP_WARNING.contains("without join"));
        let cache = Arc::new(ReferenceInfoCache::new());
        let (_dir, fasta, fai) = golden_like_pair(60);
        let (_info, handle) = read_fai_verify_in_background(&cache, fasta, fai).unwrap();
        drop(handle); // must not panic or hang
    }

    // =================================================================
    // E2 — read_reference_verifying_or_creating_fai (spec §3.11).
    // =================================================================

    /// Copy the golden FASTA into a fresh dir with **no** sibling `.fai`.
    fn golden_fasta_without_fai() -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("ref.fa");
        std::fs::copy(golden_ref_path(), &fasta).unwrap();
        assert!(!sibling_fai_path(&fasta).exists());
        (dir, fasta)
    }

    #[test]
    fn orchestrator_fai_present_verifies_in_background_and_writes_nothing() {
        let cache = Arc::new(ReferenceInfoCache::new());
        let (_dir, fasta, fai) = golden_like_pair(60);
        let fai_before = std::fs::read(&fai).unwrap();

        let (info, handle) = read_reference_verifying_or_creating_fai(&cache, fasta).unwrap();
        assert_eq!(info.md5, None, "immediate info is the .fai read");
        let handle = handle.expect(".fai present => a handle to join");
        assert_eq!(
            std::fs::read(&fai).unwrap(),
            fai_before,
            "nothing is written when the .fai is present"
        );
        let verified = handle.join().unwrap();
        assert!(verified.md5.is_some(), "join upgrades to verified info");
    }

    #[test]
    fn orchestrator_fai_absent_scans_and_writes_a_samtools_identical_fai() {
        let cache = Arc::new(ReferenceInfoCache::new());
        let (_dir, fasta) = golden_fasta_without_fai();
        let fai = sibling_fai_path(&fasta);

        let (info, handle) =
            read_reference_verifying_or_creating_fai(&cache, fasta.clone()).unwrap();
        assert!(info.md5.is_some(), "the scan path returns verified info");
        assert!(handle.is_none(), "nothing to await on the scan path");
        assert!(fai.exists(), "the .fai was written beside the FASTA");
        // Byte-identical to samtools faidx — the golden ref's committed `.fai` is the oracle.
        let written = std::fs::read(&fai).unwrap();
        let samtools = std::fs::read(sibling_fai_path(&golden_ref_path())).unwrap();
        assert_eq!(written, samtools, "written .fai == samtools faidx");
    }

    #[test]
    fn orchestrator_escape_hatch_reads_and_writes_nothing() {
        // A read-only-dir caller drops to the primitive, which writes no `.fai`.
        let cache = Arc::new(ReferenceInfoCache::new());
        let (_dir, fasta) = golden_fasta_without_fai();
        let fai = sibling_fai_path(&fasta);
        let info = cache
            .get_or_read(ReferenceSource::Fasta { fasta, fai: None })
            .unwrap();
        assert!(info.md5.is_some());
        assert!(!fai.exists(), "the primitive writes nothing");
    }

    /// Whether the filesystem actually blocks writes into a read-only directory. `false` when
    /// running as root (the dev container), where the permission bits are bypassed — the
    /// write-failure test then skips rather than pass silently. Host runs (non-root) exercise
    /// it for real.
    #[cfg(unix)]
    fn readonly_dirs_are_enforced() -> bool {
        use std::os::unix::fs::PermissionsExt;
        let dir = tempfile::tempdir().unwrap();
        std::fs::set_permissions(dir.path(), std::fs::Permissions::from_mode(0o555)).unwrap();
        let blocked = std::fs::File::create(dir.path().join("probe")).is_err();
        let _ = std::fs::set_permissions(dir.path(), std::fs::Permissions::from_mode(0o755));
        blocked
    }

    #[cfg(unix)]
    #[test]
    fn orchestrator_fai_write_failure_is_fatal() {
        use std::os::unix::fs::PermissionsExt;
        if !readonly_dirs_are_enforced() {
            eprintln!(
                "skipping orchestrator_fai_write_failure_is_fatal: read-only dirs not enforced \
                 (running as root, e.g. the container) — validated on the host instead"
            );
            return;
        }
        let cache = Arc::new(ReferenceInfoCache::new());
        let (dir, fasta) = golden_fasta_without_fai();

        // Make the reference dir read-only: the scan still reads the FASTA, but writing the
        // sibling `.fai` fails — and that failure is fatal.
        std::fs::set_permissions(dir.path(), std::fs::Permissions::from_mode(0o555)).unwrap();
        let result = read_reference_verifying_or_creating_fai(&cache, fasta.clone());
        std::fs::set_permissions(dir.path(), std::fs::Permissions::from_mode(0o755)).unwrap();
        assert!(
            matches!(result, Err(ReferenceInfoError::FaiWrite { .. })),
            "a .fai-write failure is fatal, got {result:?}"
        );

        // And the escape hatch: the primitive succeeds (a cache hit from the scan) and wrote
        // no `.fai`.
        let info = cache
            .get_or_read(ReferenceSource::Fasta {
                fasta: fasta.clone(),
                fai: None,
            })
            .unwrap();
        assert!(info.md5.is_some());
        assert!(
            !sibling_fai_path(&fasta).exists(),
            "the escape hatch wrote no .fai"
        );
    }
}
