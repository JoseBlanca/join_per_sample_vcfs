//! `AlignmentFile::open` — the validate-on-open gate, and the handle it
//! produces.
//!
//! Every check lives inside the function that opens the file, so a file is
//! either opened *and* validated or it is an `Err`: there is no window in which
//! an unvalidated handle exists. The checks run fail-fast in this order —
//! `@HD SO`, `@SQ`↔reference, the index, `@RG SM`
//! (`doc/devel/ng/spec/alignment_file.md` §3.1).
//!
//! **The invariant this establishes** is what every later layer leans on: once
//! the gate passes, `ref_id == ContigId` holds by construction for every record
//! the file can yield. That is what makes it sound for the merge one layer up
//! to compare positions across files without re-checking anything.
//!
//! ng reads `@HD SO`, `@SQ` and `@RG SM` off the noodles `sam::Header` itself
//! rather than reusing production's extractors: those are module-private, and
//! making them visible would be a production edit the ng freeze forbids
//! (`doc/devel/ng/arch/alignment_file.md` §5).
//!
//! The name says BAM but the module opens CRAM too — "BAM" in the everyday
//! sense of "the alignment file" (spec §6).

use std::fs::File;
use std::path::{Path, PathBuf};

use noodles_bam as bam;
use noodles_cram as cram;
use noodles_sam as sam;

use crate::bam::index_preflight::{
    AlignmentFileKind, AlignmentIndex, load_alignment_index, preflight_alignment_indexes,
};
use crate::fasta::{ContigEntry, ContigList};
use crate::ng::read::filtering::ReadFilterConfig;
use crate::ng::reference_info::ReferenceInfo;

use super::AlignmentFileError;

/// One opened, **validated** alignment file.
///
/// The gate of [`AlignmentFile::open`] has passed, so `ref_id == ContigId`
/// holds for every record this file will ever yield — the invariant the region
/// query, the order guard and the cross-file merge all lean on without
/// re-checking. There is no way to hold one of these unvalidated: `open` either
/// returns a validated handle or an error.
///
/// The parsed index is owned here and lives for the whole run: a region query
/// is an in-memory lookup plus a seek, never a file open or an index parse
/// (spec §3.3). Not `Clone` — it owns an index, and will own a reader pool;
/// sharing is by reference.
#[derive(Debug)]
pub struct AlignmentFile {
    path: PathBuf,
    /// Kept for the region query, which resolves a contig to a `ref_id` and
    /// hands the header to the record source. Read from C2 onwards.
    #[expect(dead_code, reason = "the region query (C2) reads this")]
    header: sam::Header,
    /// Parsed once, at open — never re-read per query (spec §3.3). Queried from
    /// C2 onwards, which is the guarantee the whole per-query cost model rests
    /// on: a query is an in-memory lookup plus a seek.
    #[expect(dead_code, reason = "the region query (C2) queries this")]
    index: AlignmentIndex,
    /// From `@RG SM`. The k-file agreement check belongs to `SampleReads`.
    sample_name: String,
    /// The `@SQ M5` tags, indexed by `ContigId` — which is sound precisely
    /// because the gate just proved this file's `@SQ` order *is* the
    /// reference's. `None` where the file carries no usable `M5`.
    ///
    /// Captured at open and compared much later, by `check_assembly` (D1),
    /// once the caller has joined `reference_info`'s background verification.
    sq_md5s: Vec<Option<[u8; 16]>>,
    /// Handed to the per-query `ReadFilter` from C4 onwards. Held on the file
    /// rather than passed per query because the filtering policy is the file's
    /// for the whole run, not the caller's per region.
    #[expect(dead_code, reason = "reads_in_region (C4) builds the filter from this")]
    filter_config: ReadFilterConfig,
}

impl AlignmentFile {
    /// Open one indexed BAM/CRAM and **validate it or fail**.
    ///
    /// Four checks, fail-fast **in this order** (spec §3.1): `@HD SO` is
    /// `coordinate`; the `@SQ` list equals `reference.contig_list()` exactly,
    /// order included; the index loads; the `@RG` records name exactly one
    /// sample. The order is cheapest-and-most-fundamental first, so a file that
    /// is wrong in several ways reports the most basic fault rather than
    /// whichever check happens to run.
    ///
    /// **The `@SQ` check is the permutation fix.** Comparing *in order* against
    /// the reference is what distinguishes this from the resolves-only check it
    /// replaces: a file whose contig list is a re-ordering of the reference's
    /// resolves every index and then fetches the wrong contig for every read.
    ///
    /// With `build_index_if_missing`, an absent index is built next to the file
    /// rather than rejected; that is a caller policy, not this module's.
    pub fn open(
        path: &Path,
        reference: &ReferenceInfo,
        filter_config: ReadFilterConfig,
        build_index_if_missing: bool,
    ) -> Result<Self, AlignmentFileError> {
        let header = read_header(path)?;

        // 1. @HD SO.
        let observed_sort_order = sort_order(&header);
        if observed_sort_order.as_deref() != Some("coordinate") {
            return Err(AlignmentFileError::NotCoordinateSorted {
                path: path.to_path_buf(),
                sort_order: observed_sort_order,
            });
        }

        // 2. @SQ vs the reference's contig table — names, lengths, order, and
        //    digests where both sides carry one. `first_disagreement` applies
        //    the absent-digest-is-a-wildcard rule itself (`fasta/mod.rs`, its
        //    `if let (Some, Some)` arm) rather than going through
        //    `ContigEntry`'s `PartialEq`, which encodes the same rule
        //    separately — so a `.fai`-only reference compares on name and
        //    length alone with no branch needed here.
        //
        //    The file is the left operand so the message reads "file value vs
        //    reference value", the direction a user needs: the reference is the
        //    authority and the file is the thing that is wrong.
        let file_contigs = contig_list(&header);
        file_contigs
            .first_disagreement(&reference.contig_list())
            .map_err(|detail| AlignmentFileError::ContigReconcile {
                path: path.to_path_buf(),
                detail,
            })?;

        // 3. The index — parsed here, once, and held for the life of the run.
        if build_index_if_missing {
            preflight_alignment_indexes(&[path.to_path_buf()], true).map_err(|source| {
                AlignmentFileError::Index {
                    path: path.to_path_buf(),
                    source,
                }
            })?;
        }
        let index = load_alignment_index(path).map_err(|source| AlignmentFileError::Index {
            path: path.to_path_buf(),
            source,
        })?;

        // 4. Exactly one @RG SM.
        let sample_name = match sample_names(&header) {
            SampleNames::One(name) => name,
            SampleNames::Several(names) => {
                return Err(AlignmentFileError::MultipleSampleNames {
                    path: path.to_path_buf(),
                    names,
                });
            }
            SampleNames::MissingTag { read_group_id } => {
                return Err(AlignmentFileError::MissingSampleName {
                    path: path.to_path_buf(),
                    read_group: Some(read_group_id),
                });
            }
            SampleNames::NoReadGroups => {
                return Err(AlignmentFileError::MissingSampleName {
                    path: path.to_path_buf(),
                    read_group: None,
                });
            }
        };

        // Free, and only sound now: check 2 proved this order is the
        // reference's, so position i really is `ContigId(i)`.
        let sq_md5s = file_contigs.entries.iter().map(|entry| entry.md5).collect();

        Ok(Self {
            path: path.to_path_buf(),
            header,
            index,
            sample_name,
            sq_md5s,
            filter_config,
        })
    }

    /// The single sample this file's `@RG` records name.
    pub fn sample_name(&self) -> &str {
        &self.sample_name
    }

    /// The `@SQ M5` tags, indexed by `ContigId`, for the deferred assembly
    /// check. `None` where the file carried no usable digest.
    pub fn sq_md5s(&self) -> &[Option<[u8; 16]>] {
        &self.sq_md5s
    }

    pub fn path(&self) -> &Path {
        &self.path
    }
}

/// Read just the SAM header, dispatching on the file's extension.
///
/// The reader is opened and dropped here: the gate needs the header, and the
/// readers a query will use come from the pool (step C1).
fn read_header(path: &Path) -> Result<sam::Header, AlignmentFileError> {
    let open_error = |source: std::io::Error| AlignmentFileError::Open {
        path: path.to_path_buf(),
        source,
    };

    let file = File::open(path).map_err(open_error)?;

    match AlignmentFileKind::from_path(path) {
        Some(AlignmentFileKind::Bam) => bam::io::Reader::new(file).read_header(),
        Some(AlignmentFileKind::Cram) => cram::io::Reader::new(file).read_header(),
        // `load_alignment_index` rejects this too, but the header read comes
        // first, so the extension has to be understood here.
        None => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "expected a '.bam' or '.cram' extension",
        )),
    }
    .map_err(open_error)
}

/// The `@HD SO` value, or `None` when the header carries no `@HD` line or no
/// sort-order tag at all.
///
/// The value is returned rather than compared here so the caller owns the
/// policy: "missing" and "queryname" are the same rejection, but they are not
/// the same *diagnosis*, and only the caller knows whether its error can say
/// which it was.
pub(crate) fn sort_order(header: &sam::Header) -> Option<String> {
    use sam::header::record::value::map::header::tag::SORT_ORDER;

    let raw = header.header()?.other_fields().get(&SORT_ORDER)?;
    Some(String::from_utf8_lossy(raw.as_ref()).into_owned())
}

/// The header's `@SQ` list as a [`ContigList`] — name, length, and the `M5`
/// digest where the file carries one.
///
/// The order is the header's own, which is the whole point: the gate compares
/// this list against the reference's **including order**, and that is what
/// catches a permutation. Building it in header order is therefore load-bearing
/// rather than incidental.
///
/// **A malformed `M5` is read as absent, not as an error — this module's
/// call, not the spec's.** Spec §3.1 rules only on a *missing* `M5`: never an
/// error, never a warning, because refusing or nagging would punish an ordinary
/// file for not offering a bonus check. A tag that is present but not 32 hex
/// characters is not quite that case — it is a header its writer got wrong —
/// but it cannot yield a digest either, so treating it the same keeps a stray
/// bad tag from being more fatal than the opportunistic check is worth. The
/// file is still fully validated on name, length and order.
///
/// The cost is that nothing downstream can tell the two apart: `Option` has no
/// room for "present but unreadable", so a file whose digests are *all*
/// malformed silently reports as untagged and wrong-assembly detection switches
/// off unremarked. Production instead rejects the file
/// (`AlignmentInputError::MalformedMd5`). If that blind spot matters, the fix
/// belongs with `check_assembly` (step D1): count the unreadable tags and let
/// `AssemblyCheck` report "verified 17 of 18, 1 unreadable" — informative
/// without becoming fatal.
pub(crate) fn contig_list(header: &sam::Header) -> ContigList {
    let entries = header
        .reference_sequences()
        .iter()
        .map(|(name, reference_sequence)| ContigEntry {
            name: String::from_utf8_lossy(name.as_ref()).into_owned(),
            length: usize::from(reference_sequence.length()) as u64,
            md5: md5_tag(reference_sequence),
        })
        .collect();

    ContigList { entries }
}

/// The `M5` tag of one `@SQ` entry, decoded from hex. `None` when the tag is
/// absent *or* unparseable — see [`contig_list`] for why those are one case.
fn md5_tag(
    reference_sequence: &sam::header::record::value::Map<
        sam::header::record::value::map::ReferenceSequence,
    >,
) -> Option<[u8; 16]> {
    use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;

    let hex = reference_sequence.other_fields().get(&MD5_CHECKSUM)?;
    decode_md5_hex(hex.as_ref())
}

/// Decode 32 hex characters into the 16 digest bytes they spell. `None` for any
/// other length or a non-hex character.
fn decode_md5_hex(hex: &[u8]) -> Option<[u8; 16]> {
    let hex: &[u8; 32] = hex.try_into().ok()?;

    let mut digest = [0u8; 16];
    for (byte, pair) in digest.iter_mut().zip(hex.chunks_exact(2)) {
        *byte = (hex_nibble(pair[0])? << 4) | hex_nibble(pair[1])?;
    }
    Some(digest)
}

fn hex_nibble(character: u8) -> Option<u8> {
    match character {
        b'0'..=b'9' => Some(character - b'0'),
        b'a'..=b'f' => Some(character - b'a' + 10),
        b'A'..=b'F' => Some(character - b'A' + 10),
        _ => None,
    }
}

/// What the `@RG` records say the file's sample is.
///
/// The gate requires **exactly one** (spec §3.1, check 4), so each way of
/// failing that is a separate variant carrying the *facts* — never a rendered
/// message. Phrasing belongs to the caller's error type, which is also the only
/// layer that knows the file's path; and cross-file agreement is a different
/// check again, one layer further up, which needs the name rather than an
/// error.
pub(crate) enum SampleNames {
    /// Every `@RG` carries an `SM`, and they all agree.
    One(String),
    /// Two or more distinct `SM` values, in first-seen order.
    Several(Vec<String>),
    /// An `@RG` record carries no `SM` tag.
    MissingTag { read_group_id: String },
    /// The file has no `@RG` records at all.
    NoReadGroups,
}

/// Collect the distinct `@RG SM` values, in first-seen order.
///
/// First-seen rather than sorted, so a message lists the names the way the file
/// does and a reader can find them.
///
/// **When a file is wrong in two ways at once, [`SampleNames::Several`] wins.**
/// A header can both name two samples and contain an `@RG` with no `SM`, and
/// the precedence is a deliberate choice rather than whichever the loop meets
/// first: "this file holds two samples" is the more consequential diagnosis, so
/// scanning finishes before deciding. Production reports whichever anomaly
/// comes first in header order
/// ([`extract_single_sample_name`](crate::bam::alignment_input)), which can hide
/// the two-sample fact behind a missing tag on a later read group.
pub(crate) fn sample_names(header: &sam::Header) -> SampleNames {
    use sam::header::record::value::map::read_group::tag::SAMPLE;

    let mut distinct: Vec<String> = Vec::new();
    let mut first_untagged: Option<String> = None;

    for (read_group_id, read_group) in header.read_groups() {
        match read_group.other_fields().get(&SAMPLE) {
            Some(raw) => {
                let sample = String::from_utf8_lossy(raw.as_ref()).into_owned();
                if !distinct.contains(&sample) {
                    distinct.push(sample);
                }
            }
            None => {
                first_untagged.get_or_insert_with(|| {
                    String::from_utf8_lossy(read_group_id.as_ref()).into_owned()
                });
            }
        }
    }

    if distinct.len() > 1 {
        return SampleNames::Several(distinct);
    }
    if let Some(read_group_id) = first_untagged {
        return SampleNames::MissingTag { read_group_id };
    }
    match distinct.into_iter().next() {
        Some(name) => SampleNames::One(name),
        None => SampleNames::NoReadGroups,
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use sam::header::record::value::Map;
    use sam::header::record::value::map::{ReadGroup, ReferenceSequence};

    use super::*;

    /// A header builder taking `(sort order, contigs, read groups)`, so each
    /// test states only the field it is about.
    fn header(
        sort_order_value: Option<&str>,
        contigs: &[(&str, usize, Option<&str>)],
        read_groups: &[(&str, Option<&str>)],
    ) -> sam::Header {
        use sam::header::record::value::map::header::tag::SORT_ORDER;
        use sam::header::record::value::map::read_group::tag::SAMPLE;
        use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;

        let mut hd = Map::<sam::header::record::value::map::Header>::default();
        if let Some(value) = sort_order_value {
            hd.other_fields_mut()
                .insert(SORT_ORDER, value.as_bytes().into());
        }

        let mut builder = sam::Header::builder().set_header(hd);

        for (name, length, md5) in contigs {
            let mut sq = Map::<ReferenceSequence>::new(NonZero::new(*length).unwrap());
            if let Some(md5) = md5 {
                sq.other_fields_mut()
                    .insert(MD5_CHECKSUM, md5.as_bytes().into());
            }
            builder = builder.add_reference_sequence(*name, sq);
        }

        for (id, sample) in read_groups {
            let mut rg = Map::<ReadGroup>::default();
            if let Some(sample) = sample {
                rg.other_fields_mut()
                    .insert(SAMPLE, sample.as_bytes().into());
            }
            builder = builder.add_read_group(*id, rg);
        }

        builder.build()
    }

    // --- @HD SO ---

    #[test]
    fn sort_order_reads_the_tag_when_present() {
        assert_eq!(
            sort_order(&header(Some("coordinate"), &[], &[])).as_deref(),
            Some("coordinate")
        );
        assert_eq!(
            sort_order(&header(Some("queryname"), &[], &[])).as_deref(),
            Some("queryname"),
            "a wrong value is reported, not silently normalised — the gate \
             wants to name what it found"
        );
    }

    #[test]
    fn sort_order_is_none_when_the_tag_or_the_hd_line_is_absent() {
        assert_eq!(sort_order(&header(None, &[], &[])), None, "@HD with no SO");
        assert_eq!(
            sort_order(&sam::Header::default()),
            None,
            "no @HD line at all"
        );
    }

    // --- @SQ ---

    #[test]
    fn contig_list_preserves_header_order_with_names_and_lengths() {
        let contigs = contig_list(&header(
            Some("coordinate"),
            &[("chr2", 200, None), ("chr1", 100, None)],
            &[],
        ));

        let names: Vec<&str> = contigs.entries.iter().map(|e| e.name.as_str()).collect();
        assert_eq!(
            names,
            vec!["chr2", "chr1"],
            "header order is kept as-is — comparing it against the reference \
             *in order* is what catches a permutation"
        );
        assert_eq!(contigs.entries[0].length, 200);
        assert_eq!(contigs.entries[1].length, 100);
    }

    #[test]
    fn contig_list_decodes_an_m5_tag_and_tolerates_its_absence() {
        let contigs = contig_list(&header(
            Some("coordinate"),
            &[
                ("chr1", 100, Some("0123456789abcdef0123456789ABCDEF")),
                ("chr2", 200, None),
            ],
            &[],
        ));

        assert_eq!(
            contigs.entries[0].md5,
            Some([
                0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0x01, 0x23, 0x45, 0x67, 0x89, 0xab,
                0xcd, 0xef
            ]),
            "decoded from hex, and upper case is accepted alongside lower"
        );
        assert_eq!(
            contigs.entries[1].md5, None,
            "a contig without M5 is ordinary input, not a fault"
        );
    }

    /// A malformed digest is deliberately read as *absent* rather than
    /// rejected, so the contig is skipped by the assembly check exactly as an
    /// untagged one is. Documented on `contig_list`; asserted here so the
    /// choice is visible rather than implicit.
    #[test]
    fn a_malformed_m5_is_read_as_absent() {
        for bad in [
            "",
            "abc",
            "zzzz56789abcdef0123456789abcdef0",
            &"a".repeat(33),
        ] {
            let contigs = contig_list(&header(
                Some("coordinate"),
                &[("chr1", 100, Some(bad))],
                &[],
            ));
            assert_eq!(
                contigs.entries[0].md5, None,
                "malformed M5 {bad:?} should read as absent"
            );
        }
    }

    // --- @RG SM ---

    #[test]
    fn one_sample_name_across_several_read_groups_is_one_name() {
        let names = sample_names(&header(
            Some("coordinate"),
            &[],
            &[("rg1", Some("NA12878")), ("rg2", Some("NA12878"))],
        ));
        assert!(matches!(names, SampleNames::One(name) if name == "NA12878"));
    }

    /// The single-file half of the sample check: two `@RG` records naming
    /// different samples is a file that cannot say whose reads it holds.
    #[test]
    fn two_distinct_sample_names_are_reported_in_first_seen_order() {
        let names = sample_names(&header(
            Some("coordinate"),
            &[],
            &[
                ("rg1", Some("NA12892")),
                ("rg2", Some("NA12878")),
                ("rg3", Some("NA12892")),
            ],
        ));
        match names {
            SampleNames::Several(names) => assert_eq!(
                names,
                vec!["NA12892", "NA12878"],
                "distinct, de-duplicated, and in the order the file lists them"
            ),
            _ => panic!("expected two distinct names"),
        }
    }

    #[test]
    fn a_file_with_no_read_groups_names_no_sample() {
        let names = sample_names(&header(Some("coordinate"), &[("chr1", 100, None)], &[]));
        assert!(matches!(names, SampleNames::NoReadGroups));
    }

    #[test]
    fn a_read_group_without_an_sm_tag_is_reported_by_id() {
        let names = sample_names(&header(
            Some("coordinate"),
            &[],
            &[("rg1", Some("NA12878")), ("rg2", None)],
        ));
        assert!(matches!(
            names,
            SampleNames::MissingTag { read_group_id } if read_group_id == "rg2"
        ));
    }

    /// Two faults at once, and the precedence is deliberate: naming two samples
    /// is the more consequential diagnosis, so it must not be hidden behind a
    /// missing tag on a later read group. This is where ng diverges from
    /// production's first-anomaly-wins, so it is pinned rather than left to
    /// whichever branch the loop reaches first.
    #[test]
    fn two_distinct_names_outrank_a_later_missing_sm_tag() {
        let names = sample_names(&header(
            Some("coordinate"),
            &[],
            &[("rg1", Some("A")), ("rg2", Some("B")), ("rg3", None)],
        ));
        match names {
            SampleNames::Several(names) => assert_eq!(names, vec!["A", "B"]),
            _ => panic!("the two-sample fact must win over the missing tag"),
        }
    }

    #[test]
    fn contig_list_is_empty_for_a_header_with_no_reference_sequences() {
        // B2 compares this against the reference, so an empty list is rejected
        // on length rather than on anything subtler.
        assert!(
            contig_list(&header(Some("coordinate"), &[], &[]))
                .entries
                .is_empty()
        );
    }

    // --- the hex decoder's edges ---

    /// The four malformed cases above put the bad character early or vary the
    /// length; none proves the **last** pair is inspected. A 32-character string
    /// that is valid hex except for its final character catches a decoder that
    /// stops one pair short — the off-by-one the `zip`/`chunks_exact` pairing
    /// could plausibly hide.
    #[test]
    fn decode_md5_hex_rejects_a_non_hex_character_in_the_final_position() {
        let mut hex = [b'0'; 32];
        assert!(decode_md5_hex(&hex).is_some(), "all-zeros is valid hex");

        hex[31] = b'z';
        assert!(decode_md5_hex(&hex).is_none(), "last nibble is inspected");

        let mut hex = [b'0'; 32];
        hex[30] = b'z';
        assert!(
            decode_md5_hex(&hex).is_none(),
            "high nibble of the last byte is inspected too"
        );
    }

    /// Range patterns fail at their edges, and a `z`-only test cannot see it.
    /// These six characters each sit one step outside a valid range.
    // -----------------------------------------------------------------
    // The gate — `AlignmentFile::open` (T1, T2, T3, T12a)
    // -----------------------------------------------------------------
    use std::fs::File;

    use noodles_bam as bam;
    use noodles_core::Position as RecordPosition;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::MappingQuality;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
    use tempfile::TempDir;

    use crate::ng::reference_info::{ReferenceSource, read_reference_info};
    use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};

    /// The fixture reference: two contigs of different lengths, so a
    /// permutation of the `@SQ` list is detectable on *name* and a re-labelling
    /// on *length*.
    const FIXTURE_CONTIGS: [(&str, usize); 2] = [("chr1", 100), ("chr2", 200)];

    fn fixture_reference(with_digests: bool) -> (TempDir, ReferenceInfo) {
        let specs: Vec<ContigSpec> = FIXTURE_CONTIGS
            .iter()
            .map(|(name, length)| ContigSpec {
                name: (*name).to_string(),
                length: *length as u64,
            })
            .collect();
        let (dir, fasta) = build_fasta(&specs).expect("build fasta");

        // The `Fasta` arm reads the genome and carries real per-contig digests;
        // the `Fai` arm cannot, so its digests are `None` and the MD5 half of
        // reconciliation is a no-op. T2 needs both.
        let source = if with_digests {
            ReferenceSource::Fasta {
                fasta: fasta.clone(),
                fai: None,
            }
        } else {
            ReferenceSource::Fai(crate::ng::reference_info::sibling_fai_path(&fasta))
        };
        let reference = read_reference_info(source).expect("read reference");
        (dir, reference)
    }

    /// A header whose `@SQ` list is `contigs`, with `SO:coordinate` and one
    /// read group naming `NA12878` unless overridden.
    fn bam_header(contigs: &[(&str, usize, Option<&str>)]) -> sam::Header {
        header(Some("coordinate"), contigs, &[("rg1", Some("NA12878"))])
    }

    fn matching_contigs() -> Vec<(&'static str, usize, Option<&'static str>)> {
        FIXTURE_CONTIGS
            .iter()
            .map(|(name, length)| (*name, *length, None))
            .collect()
    }

    fn a_record(reference_sequence_id: usize, start: usize) -> RecordBuf {
        RecordBuf::builder()
            .set_name(b"read-1")
            .set_reference_sequence_id(reference_sequence_id)
            .set_mapping_quality(MappingQuality::new(60).expect("mapq in range"))
            .set_alignment_start(RecordPosition::try_from(start).unwrap())
            .set_cigar([Op::new(Kind::Match, 10)].into_iter().collect())
            .set_sequence(Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(QualityScores::from(vec![30u8; 10]))
            .build()
    }

    /// Write a BAM with `header` and build its index next to it, so
    /// `AlignmentFile::open` finds a real indexed file on disk.
    fn indexed_bam(header: &sam::Header) -> (TempDir, PathBuf) {
        let dir = TempDir::new().expect("tempdir");
        let path = dir.path().join("sample.bam");

        let mut writer = bam::io::Writer::new(File::create(&path).expect("create bam"));
        writer.write_header(header).expect("write header");
        if !header.reference_sequences().is_empty() {
            writer
                .write_alignment_record(header, &a_record(0, 1))
                .expect("write record");
        }
        writer.try_finish().expect("finish");

        preflight_alignment_indexes(std::slice::from_ref(&path), true).expect("build index");
        (dir, path)
    }

    /// An opened fixture **plus the temp dirs its files live in**.
    ///
    /// The dirs are returned rather than dropped at helper exit because they
    /// own the BAM and the reference on disk: drop them and the handle points
    /// at deleted files. Today's assertions would survive that (the index is
    /// already parsed in memory), but the first test that actually *queries*
    /// the file would fail somewhere far from the cause.
    struct OpenedFixture {
        file: Result<AlignmentFile, AlignmentFileError>,
        _reference_dir: TempDir,
        _bam_dir: TempDir,
    }

    fn open_fixture(
        contigs: &[(&str, usize, Option<&str>)],
        reference_has_digests: bool,
    ) -> OpenedFixture {
        open_fixture_with_header(&bam_header(contigs), reference_has_digests)
    }

    fn open_fixture_with_header(
        header: &sam::Header,
        reference_has_digests: bool,
    ) -> OpenedFixture {
        let (reference_dir, reference) = fixture_reference(reference_has_digests);
        let (bam_dir, path) = indexed_bam(header);
        OpenedFixture {
            file: AlignmentFile::open(&path, &reference, ReadFilterConfig::default(), false),
            _reference_dir: reference_dir,
            _bam_dir: bam_dir,
        }
    }

    /// The `ContigReconcile` detail for a header that should fail check 2.
    fn reconcile_detail(
        contigs: &[(&str, usize, Option<&str>)],
        reference_has_digests: bool,
    ) -> String {
        match open_fixture(contigs, reference_has_digests)
            .file
            .expect_err("must not open")
        {
            AlignmentFileError::ContigReconcile { detail, .. } => detail,
            other => panic!("expected ContigReconcile, got {other:?}"),
        }
    }

    #[test]
    fn a_matching_file_opens_and_exposes_its_sample_and_digests() {
        let fixture = open_fixture(&matching_contigs(), false);
        let file = fixture.file.as_ref().expect("the file matches");

        assert_eq!(file.sample_name(), "NA12878");
        assert_eq!(
            file.sq_md5s().len(),
            2,
            "one slot per contig, indexed by ContigId"
        );
        assert!(
            file.sq_md5s().iter().all(Option::is_none),
            "this fixture's @SQ lines carry no M5, which is ordinary input"
        );
    }

    /// **T1 — the permutation hole, closed.**
    ///
    /// `@SQ` lists the reference's contigs in the wrong order. Every `ref_id`
    /// still *resolves*, which is why the check this replaces let it through
    /// and then fetched the wrong contig for every read. Order-aware equality
    /// catches it on the first transposed name.
    ///
    /// Mutation-verified below by
    /// `a_resolves_only_check_accepts_the_permutation_this_gate_rejects`.
    #[test]
    fn t1_a_permuted_sq_list_is_rejected_naming_the_first_transposed_contig() {
        let permuted = vec![("chr2", 200, None), ("chr1", 100, None)];

        let detail = reconcile_detail(&permuted, false);

        assert_eq!(
            detail, "name disagreement at index 0 ('chr2' vs 'chr1')",
            "the first transposed position, file value before reference value"
        );
    }

    /// **Why T1 needs an ordered comparison — the superseded probe, recorded.**
    ///
    /// `ReadFilter::new`'s check asks only "does every `@SQ` index resolve in
    /// the reference?", fetching each contig and accepting if the fetch
    /// succeeds. Run here against the permuted header, it **passes** — order is
    /// never consulted, so the file goes on to fetch the wrong contig for every
    /// read.
    ///
    /// **This is documentation, not the mutation guard.** It does not call
    /// `AlignmentFile::open`, so gutting the gate's comparison would leave it
    /// green; T1 is what catches that. The real mutation was performed against
    /// the gate itself during B2 — the `first_disagreement` call was replaced
    /// with a resolves-only length check and T1 duly failed — and this test is
    /// the standing record of *what* the old check did, so the reason the gate
    /// is stricter does not have to be taken on trust.
    #[test]
    fn the_superseded_resolves_only_probe_cannot_see_order() {
        use crate::ng::ref_seq::{InMemoryRefSeq, RefSeq};
        use crate::ng::types::ContigId;

        let (_reference_dir, reference) = fixture_reference(false);
        let permuted = bam_header(&[("chr2", 200, None), ("chr1", 100, None)]);

        // The superseded check: fetch every contig the file's `@SQ` list can
        // name, and accept if each one resolves. Order is never consulted.
        let reference_bases = InMemoryRefSeq::from_contigs(
            FIXTURE_CONTIGS
                .iter()
                .map(|(_, length)| vec![b'A'; *length])
                .collect(),
        );
        let mut probe = Vec::new();
        let every_contig_resolves = (0..permuted.reference_sequences().len()).all(|index| {
            reference_bases
                .fetch_into(ContigId(index as u32), 1, 0, &mut probe)
                .is_ok()
        });
        assert!(
            every_contig_resolves,
            "every @SQ index of the permuted file resolves — so the old probe \
             accepts it, and then every read fetches the wrong contig"
        );

        // The gate's own comparison rejects that very same header.
        assert!(
            contig_list(&permuted)
                .first_disagreement(&reference.contig_list())
                .is_err(),
            "the ordered comparison must reject what the probe accepted"
        );
    }

    /// **T2 — name, length and count mismatches**, each naming the right field
    /// and index.
    #[test]
    fn t2_name_length_and_count_mismatches_are_each_named() {
        let cases = [
            (
                vec![("chrX", 100, None), ("chr2", 200, None)],
                "name disagreement at index 0",
            ),
            (
                vec![("chr1", 100, None), ("chr2", 999, None)],
                "length disagreement at index 1",
            ),
            (
                vec![("chr1", 100, None)],
                "@SQ list length differs (1 vs 2)",
            ),
        ];

        for (contigs, expected) in cases {
            let detail = reconcile_detail(&contigs, false);
            assert!(
                detail.contains(expected),
                "expected {expected:?}, got {detail:?}"
            );
        }
    }

    /// **T2, the digest half.** A wrong `@SQ M5` is caught at open **only**
    /// when the `ReferenceInfo` in hand carries digests — the `Fasta` arm. Under
    /// a `.fai`-only table the comparison is a wildcard and the same file opens
    /// cleanly, which is the behaviour production ships and what makes the
    /// deferred `check_assembly` (D1) necessary rather than redundant.
    #[test]
    fn t2_a_wrong_m5_is_caught_only_when_the_reference_carries_digests() {
        let wrong_digest = "ffffffffffffffffffffffffffffffff";
        let contigs = vec![
            ("chr1", 100, Some(wrong_digest)),
            ("chr2", 200, Some(wrong_digest)),
        ];

        let detail = reconcile_detail(&contigs, true);
        assert!(
            detail.contains("md5 disagreement at index 0"),
            "got: {detail}"
        );

        assert!(
            open_fixture(&contigs, false).file.is_ok(),
            "against a .fai-only reference the digest is a wildcard, so the \
             same file opens — the gap check_assembly closes later"
        );
    }

    /// **T3 — `SO` wrong or missing**, rejected before anything else, and the
    /// message says which it was.
    #[test]
    fn t3_a_file_that_is_not_coordinate_sorted_is_rejected_at_open() {
        let contigs = matching_contigs();

        let queryname = header(Some("queryname"), &contigs, &[("rg1", Some("NA12878"))]);
        let error = open_fixture_with_header(&queryname, false)
            .file
            .expect_err("must not open");
        match &error {
            AlignmentFileError::NotCoordinateSorted { sort_order, .. } => assert_eq!(
                sort_order.as_deref(),
                Some("queryname"),
                "carries the observed value as a fact, not pre-quoted prose"
            ),
            other => panic!("expected NotCoordinateSorted, got {other:?}"),
        }
        assert!(
            error.to_string().contains("@HD SO is 'queryname'"),
            "and the rendered message quotes it: {error}"
        );

        let no_sort_order = header(None, &contigs, &[("rg1", Some("NA12878"))]);
        let error = open_fixture_with_header(&no_sort_order, false)
            .file
            .expect_err("must not open");
        match &error {
            AlignmentFileError::NotCoordinateSorted { sort_order, .. } => {
                assert_eq!(sort_order.as_deref(), None)
            }
            other => panic!("expected NotCoordinateSorted, got {other:?}"),
        }
        assert!(error.to_string().contains("@HD SO is missing"), "{error}");
    }

    /// The gate's fail-fast order is spec §3.1's — `SO` → `@SQ` → index → `SM`
    /// — and it is only observable through a file that fails two checks at
    /// once. Each boundary gets its own case, so re-ordering any adjacent pair
    /// fails a test rather than quietly changing which fault a user is told
    /// about.
    #[test]
    fn the_four_checks_run_in_the_specified_order() {
        // SO before @SQ.
        let bad_sort_and_bad_contigs = header(
            Some("queryname"),
            &[("chrX", 1, None)],
            &[("rg1", Some("NA12878"))],
        );
        assert!(matches!(
            open_fixture_with_header(&bad_sort_and_bad_contigs, false).file,
            Err(AlignmentFileError::NotCoordinateSorted { .. })
        ));

        // @SQ before the index: bad contigs on a BAM with no index at all.
        let bad_contigs = bam_header(&[("chrX", 1, None)]);
        assert!(matches!(
            open_unindexed_with_header(&bad_contigs).1,
            Err(AlignmentFileError::ContigReconcile { .. })
        ));

        // The index before @RG SM: no index and no sample.
        let no_sample = header(Some("coordinate"), &matching_contigs(), &[]);
        assert!(matches!(
            open_unindexed_with_header(&no_sample).1,
            Err(AlignmentFileError::Index { .. })
        ));
    }

    /// Write a BAM with **no** index beside it and try to open it. Returns the
    /// temp dirs so the file outlives the call.
    fn open_unindexed_with_header(
        header: &sam::Header,
    ) -> (
        (TempDir, TempDir),
        Result<AlignmentFile, AlignmentFileError>,
    ) {
        let (reference_dir, reference) = fixture_reference(false);
        let dir = TempDir::new().expect("tempdir");
        let path = dir.path().join("unindexed.bam");

        let mut writer = bam::io::Writer::new(File::create(&path).expect("create"));
        writer.write_header(header).expect("write header");
        if !header.reference_sequences().is_empty() {
            writer
                .write_alignment_record(header, &a_record(0, 1))
                .expect("write record");
        }
        writer.try_finish().expect("finish");

        let opened = AlignmentFile::open(&path, &reference, ReadFilterConfig::default(), false);
        ((reference_dir, dir), opened)
    }

    /// **T12a — one file naming two samples.** The cross-file half (two files
    /// naming different samples) is `SampleReads`' T12b.
    #[test]
    fn t12a_a_file_whose_read_groups_name_two_samples_is_rejected_at_open() {
        let contigs = matching_contigs();
        let two_samples = header(
            Some("coordinate"),
            &contigs,
            &[("rg1", Some("NA12878")), ("rg2", Some("NA12892"))],
        );

        match open_fixture_with_header(&two_samples, false)
            .file
            .expect_err("must not open")
        {
            AlignmentFileError::MultipleSampleNames { names, .. } => {
                assert_eq!(names, vec!["NA12878", "NA12892"])
            }
            other => panic!("expected MultipleSampleNames, got {other:?}"),
        }
    }

    /// The gate's mapping of `SampleNames::MissingTag`, and the `Some` branch
    /// of `MissingSampleName`'s message — the only non-trivial formatting in
    /// the new variant, and until now exercised by nothing.
    #[test]
    fn a_read_group_with_no_sm_tag_is_rejected_naming_that_read_group() {
        let contigs = matching_contigs();
        let untagged = header(Some("coordinate"), &contigs, &[("rg1", None)]);

        let error = open_fixture_with_header(&untagged, false)
            .file
            .expect_err("must not open");
        match &error {
            AlignmentFileError::MissingSampleName { read_group, .. } => {
                assert_eq!(read_group.as_deref(), Some("rg1"))
            }
            other => panic!("expected MissingSampleName, got {other:?}"),
        }
        assert!(
            error.to_string().contains("@RG 'rg1' has no SM tag"),
            "{error}"
        );
    }

    /// The other half of "exactly one sample": a file with no `@RG` at all
    /// cannot say whose reads it holds, so the sample layer would have nothing
    /// to agree on.
    #[test]
    fn a_file_naming_no_sample_is_rejected_at_open() {
        let contigs = matching_contigs();
        let no_read_groups = header(Some("coordinate"), &contigs, &[]);

        match open_fixture_with_header(&no_read_groups, false)
            .file
            .expect_err("must not open")
        {
            AlignmentFileError::MissingSampleName { read_group, .. } => {
                assert_eq!(read_group, None)
            }
            other => panic!("expected MissingSampleName, got {other:?}"),
        }
    }

    /// A missing index is an error, not a silent whole-file scan — and with the
    /// build flag it is repaired instead.
    #[test]
    fn a_missing_index_is_rejected_unless_the_caller_asks_for_one() {
        let header = bam_header(&matching_contigs());
        let (dirs, opened) = open_unindexed_with_header(&header);

        assert!(
            matches!(opened, Err(AlignmentFileError::Index { .. })),
            "an unindexed file is an error, never a silent whole-file scan"
        );

        // The same file in the same place — only the caller's policy differs.
        let (_reference_dir, bam_dir) = &dirs;
        let (_fresh_reference_dir, reference) = fixture_reference(false);
        let path = bam_dir.path().join("unindexed.bam");
        assert!(
            AlignmentFile::open(&path, &reference, ReadFilterConfig::default(), true).is_ok(),
            "with build_index_if_missing the index is created next to the file"
        );
    }

    #[test]
    fn hex_nibble_rejects_the_characters_bordering_each_valid_range() {
        for character in [b'/', b':', b'`', b'g', b'@', b'G'] {
            assert_eq!(
                hex_nibble(character),
                None,
                "{:?} borders a valid range but is not hex",
                character as char
            );
        }
        assert_eq!(hex_nibble(b'0'), Some(0));
        assert_eq!(hex_nibble(b'9'), Some(9));
        assert_eq!(hex_nibble(b'a'), Some(10));
        assert_eq!(hex_nibble(b'f'), Some(15));
        assert_eq!(hex_nibble(b'A'), Some(10));
        assert_eq!(hex_nibble(b'F'), Some(15));
    }
}
