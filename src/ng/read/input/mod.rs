//! Step 1's **input edge**: turning alignment files on disk into ordered,
//! filtered streams of reads for a region.
//!
//! Two layers, split by subject. [`open_bam`] and [`region_query`] own the
//! reader logic whose subject is **one alignment file** — the validate-on-open
//! gate, the reader pool, the BAM/CRAM region queries and the order guard
//! (`doc/devel/ng/spec/alignment_file.md`). [`merge`] and this module's
//! `SampleReads` own what is true only of **the sample** — the cross-file
//! checks and the k-way merge (`doc/devel/ng/spec/sample_reads.md`). Both
//! layers' error types live here, where the arch docs put them.
//!
//! The layering is strict: the sample layer consumes what the file layer
//! produces and contains no reader logic of its own. It is also where the
//! module's two guarantees come from — every file's `ref_id` was proved equal
//! to its `ContigId` at open, and every per-file stream was proved
//! coordinate-monotonic while streaming — which is what makes it sound for the
//! merge to compare positions *across* files without re-checking anything.
//!
//! This is step 1's input edge rather than a new step, so it lives under
//! `read/` beside `filtering.rs` and reuses that module's
//! [`RecordSource`](super::RecordSource)/[`RawRecord`](super::RawRecord) seam
//! (`doc/devel/ng/arch/module_layout.md` principle 1, note b).

pub mod merge;
pub mod open_bam;
pub mod region_query;

#[cfg(test)]
pub(crate) mod test_fixtures;

use std::io;
use std::path::{Path, PathBuf};

use crate::bam::errors::AlignmentIndexError;
use crate::ng::read::filtering::ReadFilterError;
use crate::ng::reference_info::ReferenceInfo;
use crate::ng::types::{GenomePosition, GenomeRegion};
use crate::pop_var_caller::common::format_md5_hex;

/// Everything that can go wrong for **one** alignment file.
///
/// Split by *when* it fires. The first five are returned by
/// `AlignmentFile::open` before any read flows, so a handle never exists in an
/// unvalidated state. [`Self::OutOfOrderRead`] and [`Self::Filter`] arrive in
/// the item stream instead, after which the iterator fuses — the first `Err` is
/// yielded once and then `None`, so a caller writing `let read = read?;` cannot
/// mistake a fatal condition for a clean end of input.
///
/// Design: `doc/devel/ng/spec/alignment_file.md` §4,
/// `doc/devel/ng/arch/alignment_file.md` §2.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum AlignmentFileError {
    /// The file could not be opened, or its header could not be parsed.
    #[error("opening alignment file '{path}' failed")]
    Open {
        path: PathBuf,
        #[source]
        source: io::Error,
    },

    /// `@HD SO` is absent or is not `coordinate` (spec §3.1, check 1).
    ///
    /// Cheap and checked first, but **not** a substitute for the streaming
    /// order guard: this catches a file that does not *claim* to be sorted,
    /// while [`Self::OutOfOrderRead`] catches the file that claims it and lies.
    ///
    /// `sort_order` carries what the header actually said — `queryname` and a
    /// missing tag are the same rejection but not the same diagnosis, and a
    /// user whose file says `SO:queryname` should be told that rather than
    /// left to guess. `None` means the tag (or the whole `@HD` line) was
    /// absent. (Arch §2 sketched this variant with only `path`; signatures
    /// there are illustrative, and production's equivalent carries the observed
    /// value too.)
    #[error(
        "alignment file '{path}' is not coordinate-sorted (@HD SO is {})",
        match sort_order {
            Some(value) => format!("'{value}'"),
            None => "missing".to_string(),
        }
    )]
    NotCoordinateSorted {
        path: PathBuf,
        sort_order: Option<String>,
    },

    /// The file's `@SQ` list is not the reference's contig table (spec §3.1,
    /// check 2). `detail` comes from `ContigList::first_disagreement` and names
    /// the first differing field and index.
    ///
    /// This is the variant that closes the **permutation hole**: a file whose
    /// `@SQ` list is a re-ordering of the reference's passes a resolves-only
    /// check and then fetches the wrong contig for every read.
    #[error("alignment file '{path}' does not match the reference contig table: {detail}")]
    ContigReconcile { path: PathBuf, detail: String },

    /// An `@SQ` line carries an `M5` tag that is not a 32-character hex digest.
    ///
    /// A *missing* `M5` is ordinary input and never an error (spec §3.1) — the
    /// digest check is opportunistic, and refusing a file for not offering a
    /// bonus would punish the common case. A **malformed** one is different: it
    /// is a header its writer got wrong, and reading it as absent would hide a
    /// real defect. Worse, it would hide it *silently and completely* — a file
    /// whose digests are all malformed would report as untagged, and
    /// wrong-assembly detection would switch off with nobody told. Errors are
    /// not passed silently, so this is fatal at open.
    #[error("alignment file '{path}': @SQ '{contig}' has a malformed M5 tag ({detail})")]
    MalformedMd5 {
        path: PathBuf,
        contig: String,
        detail: String,
    },

    /// The index is missing or unparseable (spec §3.1, check 3).
    #[error("loading the index for alignment file '{path}' failed")]
    Index {
        path: PathBuf,
        #[source]
        source: AlignmentIndexError,
    },

    /// The `@RG` records name more than one sample (spec §3.1, check 4).
    /// Agreement *across* files is not a property of one file and is checked
    /// one layer up, as [`IngestError::SampleNameMismatch`].
    #[error("alignment file '{path}' names more than one sample: {}", names.join(", "))]
    MultipleSampleNames { path: PathBuf, names: Vec<String> },

    /// The `@RG` records name **no** sample — either the file has none at all,
    /// or one carries no `SM` tag (spec §3.1, check 4).
    ///
    /// The other half of "exactly one sample". Arch §2's variant list covered
    /// only the *more than one* case, but a file that cannot say whose reads it
    /// holds fails the same check for the same reason, and the sample layer
    /// above has nothing to compare. Production rejects these files too.
    ///
    /// `read_group` distinguishes the two causes: `Some(id)` names the record
    /// missing its tag, `None` means the file has no `@RG` lines.
    #[error(
        "alignment file '{path}' names no sample: {}",
        match read_group {
            Some(id) => format!("@RG '{id}' has no SM tag"),
            None => "the file has no @RG records".to_string(),
        }
    )]
    MissingSampleName {
        path: PathBuf,
        read_group: Option<String>,
    },

    /// A record regressed in genome position: the file claims `SO:coordinate`
    /// and is not sorted (spec §3.2). Carries **both** keys so the message can
    /// say where the file breaks rather than only that it does.
    #[error(
        "alignment file '{path}' is not sorted: a read at contig {} position {} \
         follows one at contig {} position {}",
        current.contig.get(), current.position.get(),
        previous.contig.get(), previous.position.get()
    )]
    OutOfOrderRead {
        path: PathBuf,
        previous: GenomePosition,
        current: GenomePosition,
    },

    /// Step 1 hit a fatal condition — a failed source read, a failed decode, or
    /// filter #8's reference fetch.
    #[error("read filtering failed")]
    Filter(#[source] ReadFilterError),

    /// The requested region is invalid, or names a contig absent from the
    /// reference.
    #[error(
        "invalid region: contig {} [{}, {}]",
        region.contig.get(), region.start.get(), region.end.get()
    )]
    Region { region: GenomeRegion },

    /// A CRAM was opened against a reference that carries no FASTA — a
    /// `.fai`-only `ReferenceInfo`.
    ///
    /// CRAM stores its sequences as differences from the reference, so decoding
    /// a record *requires the bases themselves*. A `.fai` describes a genome's
    /// geometry and holds no bases, so there is nothing to decode against. BAM
    /// is unaffected: it stores its own sequences.
    ///
    /// Raised at open rather than at the first query, so the run fails while
    /// the cause is still obvious (owner, 2026-07-20).
    #[error(
        "'{path}' is a CRAM, which can only be decoded against the reference \
         sequence itself, but the reference was read from a .fai index (which \
         holds no bases) — supply the reference FASTA"
    )]
    CramNeedsReferenceFasta { path: PathBuf },

    /// Querying the parsed index for a region's chunks failed.
    ///
    /// Distinct from [`Self::Region`]: the region has already been checked
    /// against the contig table by the time this can fire, so a failure here is
    /// the *index* being unusable — corrupt or truncated — not the caller
    /// asking for something silly. Reporting it as a bad region would name the
    /// wrong culprit and throw away the cause.
    #[error(
        "querying the index of '{path}' for contig {} [{}, {}] failed",
        region.contig.get(), region.start.get(), region.end.get()
    )]
    IndexQuery {
        path: PathBuf,
        region: GenomeRegion,
        #[source]
        source: io::Error,
    },
}

/// What [`check_assembly`] found: how many contigs could be compared, out of
/// how many there are.
///
/// Returned as **information, not a complaint**. `compared == 0` is ordinary —
/// it means no contig carried a digest on both sides — and a caller that
/// ignores this is behaving correctly. It exists so a caller that wants to say
/// "assembly verified for 18 of 24 contigs" can, rather than having to imply a
/// guarantee it does not have.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AssemblyCheck {
    /// Contigs where **both** sides carried a digest and the two agreed.
    pub compared: usize,
    /// Contigs in the reference — the denominator, so `compared` can be read
    /// as a fraction of the genome actually verified.
    pub total: usize,
}

/// A contig whose `@SQ M5` disagrees with the reference's digest: the file is
/// aligned to a **different assembly** with the same contig names and lengths.
///
/// A separate type rather than an [`AlignmentFileError`] variant, because
/// [`check_assembly`] runs long after the streams are gone — so this can never
/// appear as a stream item, and a variant of that enum would be unreachable
/// there (spec §4).
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error(
    "alignment file '{path}' is aligned to a different assembly: contig \
     '{contig}' has M5 {} but the reference's is {}",
    format_md5_hex(*observed),
    format_md5_hex(*expected)
)]
pub struct AssemblyMismatch {
    pub path: PathBuf,
    pub contig: String,
    pub expected: [u8; 16],
    pub observed: [u8; 16],
}

/// Compare the `@SQ M5` tags captured at open against a **verified**
/// `ReferenceInfo` — the one carrying real per-contig digests, available once
/// the caller has joined `reference_info`'s background verification.
///
/// **Why this is deferred rather than done at open.** With a `.fai` present,
/// `reference_info` hands back the cheap table immediately and verifies the
/// genome on a background thread, so the `ReferenceInfo` a file is opened
/// against usually carries no digests at all. Waiting for them would mean
/// blocking startup on a whole-genome read. Instead the tags are captured for
/// free at open and compared here, at the point the caller joins that
/// verification — before anything is published. A wrong-assembly run therefore
/// wastes its work before aborting, which is the right trade for something
/// catastrophic and rare (spec §3.1).
///
/// **A contig is compared only when both sides carry a digest**, and a contig
/// that carries none on either side is skipped, not failed. Many BAM/CRAMs have
/// no `M5` at all, and a `.fai`-only reference has none either; refusing or
/// nagging would punish the common case for missing an opportunistic extra
/// check. (A *malformed* tag is a different matter and was already rejected at
/// open, so everything reaching here is absent or valid.)
///
/// `observed` is indexed by `ContigId`, parallel to `verified.contigs` — which
/// the open gate guarantees, having proved the file's `@SQ` list *is* the
/// reference's contig table.
pub fn check_assembly(
    path: &Path,
    observed: &[Option<[u8; 16]>],
    verified: &ReferenceInfo,
) -> Result<AssemblyCheck, AssemblyMismatch> {
    // The open gate proves the file's `@SQ` list *is* the reference's contig
    // table, so these are parallel. Checked because this is `pub` and takes two
    // loose slices: if they ever diverged, pair `i` would compare the file's
    // contig against a *different* reference contig — missing a real mismatch,
    // or naming the wrong contig in the error.
    debug_assert_eq!(
        observed.len(),
        verified.contigs.len(),
        "observed digests must be parallel to the reference's contigs"
    );

    let mut compared = 0;

    for (observed_digest, contig) in observed.iter().zip(verified.contigs.iter()) {
        let (Some(observed_digest), Some(expected)) = (observed_digest, contig.md5) else {
            continue;
        };

        if *observed_digest != expected {
            return Err(AssemblyMismatch {
                path: path.to_path_buf(),
                contig: contig.name.clone(),
                expected,
                observed: *observed_digest,
            });
        }
        compared += 1;
    }

    Ok(AssemblyCheck {
        compared,
        // The reference's count, not the number of pairs walked: if the two
        // slices disagree in length — which the open gate rules out — the
        // shortfall shows up as `compared < total` rather than vanishing.
        total: verified.contigs.len(),
    })
}

/// Sample-level failures — the layer above [`AlignmentFileError`].
///
/// The two own variants are the checks that can only be made *across* files.
/// Everything else arrives from one file and is **wrapped** with the index of
/// the file that raised it, rather than flattened: flattening would duplicate
/// every [`AlignmentFileError`] variant and lose which file was at fault
/// (`doc/devel/ng/arch/sample_reads.md` §2).
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum IngestError {
    /// The k files do not name the same sample. Fires at `SampleReads::open`,
    /// before any read (spec `sample_reads.md` §3.1).
    ///
    /// This guard earns its keep because a sample's files are usually separate
    /// *experiments* gathered by hand from different projects, where picking up
    /// a foreign file is the realistic failure mode.
    #[error(
        "the files of one sample name different samples: {}",
        files.iter().zip(names.iter())
            .map(|(file, name)| format!("'{}' names '{name}'", file.display()))
            .collect::<Vec<_>>().join(", ")
    )]
    SampleNameMismatch {
        files: Vec<PathBuf>,
        names: Vec<String>,
    },

    /// The identical read surfaced from two files, which means the caller
    /// passed the same file twice (or two files with overlapping content).
    ///
    /// Not a deduplication: across files there is no such thing as a legitimate
    /// duplicate, because reads from different experiments are different reads.
    /// So this is an **input-sanity error**, never a silent drop (spec §3.2).
    #[error(
        "read '{}' appeared in both input file #{} and input file #{} \
         at contig {} position {}",
        String::from_utf8_lossy(qname), files.0, files.1,
        key.contig.get(), key.position.get()
    )]
    DuplicateReadAcrossFiles {
        qname: Vec<u8>,
        key: GenomePosition,
        files: (usize, usize),
    },

    /// Anything one file raised, at open or mid-stream, tagged with which one.
    #[error("alignment file #{source_file_index} failed")]
    File {
        source_file_index: usize,
        #[source]
        source: AlignmentFileError,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::ng::types::{ContigId, Position};

    fn genome_position(contig: u32, position: u64) -> GenomePosition {
        GenomePosition {
            contig: ContigId(contig),
            position: Position(position),
        }
    }

    /// The error messages are the user-facing surface of this module, and these
    /// two exist precisely to name *where* a fault is — which read follows
    /// which, and which two files collided.
    ///
    /// Asserted against the whole rendered string rather than by substring,
    /// because the mistake worth catching is **transposing `previous` and
    /// `current`**: that turns "where the file breaks" into a lie while leaving
    /// every substring assertion green.
    #[test]
    fn out_of_order_message_says_which_read_follows_which() {
        let message = AlignmentFileError::OutOfOrderRead {
            path: PathBuf::from("/data/sample.bam"),
            previous: genome_position(2, 5000),
            current: genome_position(2, 120),
        }
        .to_string();

        assert_eq!(
            message,
            "alignment file '/data/sample.bam' is not sorted: \
             a read at contig 2 position 120 follows one at contig 2 position 5000"
        );
    }

    // -----------------------------------------------------------------
    // check_assembly (D1) — T2b
    // -----------------------------------------------------------------
    use crate::ng::reference_info::ContigInfo;

    fn contig(name: &str, md5: Option<[u8; 16]>) -> ContigInfo {
        ContigInfo {
            name: name.to_string(),
            length: 100,
            offset: 0,
            line_bases: 60,
            line_width: 61,
            md5,
        }
    }

    fn verified_reference(digests: &[Option<[u8; 16]>]) -> ReferenceInfo {
        ReferenceInfo {
            md5: None,
            contigs: digests
                .iter()
                .enumerate()
                .map(|(i, md5)| contig(&format!("chr{}", i + 1), *md5))
                .collect(),
            fasta_path: None,
        }
    }

    const A: [u8; 16] = [0xaa; 16];
    const B: [u8; 16] = [0xbb; 16];

    /// **T2b, the fault it exists to catch.** A file aligned to a different
    /// assembly with the same contig names and lengths passes every check at
    /// open — that is exactly why this check is worth having, and why
    /// production, which trusts the `@SQ M5`, cannot make it.
    #[test]
    fn t2b_a_disagreeing_digest_names_the_contig_and_both_values() {
        let error = check_assembly(
            Path::new("/data/sample.bam"),
            &[Some(A), Some(B)],
            &verified_reference(&[Some(A), Some(A)]),
        )
        .expect_err("the second contig disagrees");

        assert_eq!(error.contig, "chr2");
        assert_eq!(error.expected, A);
        assert_eq!(error.observed, B);

        // Asserted whole, because the mistake worth catching is transposing
        // `expected` and `observed` in the message: that inverts the diagnosis
        // and sends the user after the wrong artefact, while leaving every
        // "contains" assertion green.
        assert_eq!(
            error.to_string(),
            format!(
                "alignment file '/data/sample.bam' is aligned to a different \
                 assembly: contig 'chr2' has M5 {} but the reference's is {}",
                "bb".repeat(16),
                "aa".repeat(16)
            )
        );
    }

    /// **T2b, the ordinary case.** A file with no `M5` tags at all is not a
    /// fault — it simply cannot be checked. `compared == 0` says so plainly
    /// rather than implying a guarantee that was never obtained.
    #[test]
    fn t2b_a_file_without_digests_is_not_an_error_and_reports_nothing_compared() {
        let check = check_assembly(
            Path::new("/data/sample.bam"),
            &[None, None, None],
            &verified_reference(&[Some(A), Some(A), Some(A)]),
        )
        .expect("no tags is ordinary input");

        assert_eq!(
            check,
            AssemblyCheck {
                compared: 0,
                total: 3
            }
        );
    }

    /// **T2b, the partial case.** `compared` counts exactly the contigs that
    /// could be checked — so "verified 2 of 3" is a statement a caller can make
    /// truthfully.
    #[test]
    fn t2b_compared_counts_exactly_the_contigs_both_sides_carry() {
        let check = check_assembly(
            Path::new("/data/sample.bam"),
            &[Some(A), None, Some(A)],
            &verified_reference(&[Some(A), Some(A), Some(A)]),
        )
        .expect("the tagged contigs agree");
        assert_eq!(
            check,
            AssemblyCheck {
                compared: 2,
                total: 3
            }
        );

        // And symmetrically: a digest on the file's side is useless if the
        // reference has none, which is the whole `.fai`-only situation.
        let check = check_assembly(
            Path::new("/data/sample.bam"),
            &[Some(A), Some(A), Some(A)],
            &verified_reference(&[Some(A), None, None]),
        )
        .expect("only the first contig is comparable");
        assert_eq!(
            check,
            AssemblyCheck {
                compared: 1,
                total: 3
            }
        );
    }

    /// The headline case a caller wants to report: everything tagged, all
    /// agreeing, so `compared == total`.
    #[test]
    fn a_fully_tagged_agreeing_file_reports_every_contig_verified() {
        let check = check_assembly(
            Path::new("/data/sample.bam"),
            &[Some(A), Some(A), Some(A)],
            &verified_reference(&[Some(A), Some(A), Some(A)]),
        )
        .expect("all agree");
        assert_eq!(
            check,
            AssemblyCheck {
                compared: 3,
                total: 3
            }
        );
    }

    /// A `.fai`-only reference carries no digests at all, so nothing is
    /// comparable and nothing is claimed — the case that makes this check
    /// *deferred* rather than redundant with the one at open.
    #[test]
    fn nothing_is_comparable_against_a_fai_only_reference() {
        let check = check_assembly(
            Path::new("/data/sample.bam"),
            &[Some(A), Some(B)],
            &verified_reference(&[None, None]),
        )
        .expect("a .fai has no digests to disagree with");
        assert_eq!(
            check,
            AssemblyCheck {
                compared: 0,
                total: 2
            }
        );
    }

    /// The first disagreement wins, so the message names the earliest contig
    /// rather than an arbitrary one.
    #[test]
    fn the_first_disagreeing_contig_is_the_one_reported() {
        let error = check_assembly(
            Path::new("/data/sample.bam"),
            &[Some(B), Some(B)],
            &verified_reference(&[Some(A), Some(A)]),
        )
        .expect_err("both disagree");
        assert_eq!(error.contig, "chr1");
    }

    /// The colliding read's name is what lets a user grep the inputs and
    /// confirm the diagnosis, and the file indices have to read as *files*
    /// rather than as another coordinate pair sitting beside `contig 0
    /// position 99`.
    #[test]
    fn duplicate_read_message_names_the_read_and_both_files() {
        let message = IngestError::DuplicateReadAcrossFiles {
            qname: b"read-1".to_vec(),
            key: genome_position(0, 99),
            files: (0, 1),
        }
        .to_string();

        assert_eq!(
            message,
            "read 'read-1' appeared in both input file #0 and input file #1 \
             at contig 0 position 99"
        );
    }

    /// The variant carries `files` as well as `names`, and the whole point of
    /// the check is that someone picked up a stray file — so the message must
    /// name the *path*, not just repeat the sample names the user already
    /// knows disagree.
    #[test]
    fn sample_name_mismatch_message_names_the_offending_paths() {
        let message = IngestError::SampleNameMismatch {
            files: vec![PathBuf::from("/data/a.bam"), PathBuf::from("/data/b.bam")],
            names: vec!["NA12878".to_string(), "NA12892".to_string()],
        }
        .to_string();

        assert_eq!(
            message,
            "the files of one sample name different samples: \
             '/data/a.bam' names 'NA12878', '/data/b.bam' names 'NA12892'"
        );
    }
}
