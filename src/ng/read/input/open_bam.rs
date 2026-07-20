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

// The header readers below are the gate's parts, built (step B1) one commit
// before the gate that calls them (step B2). They have tests but no production
// caller yet, which `-D warnings` reads as dead code.
//
// `expect` rather than `allow` on purpose: the moment B2 calls these, the
// expectation goes unfulfilled and the build fails *naming this line*, so
// removing it is enforced by the compiler rather than remembered.
//
// `not(test)`, because under `cfg(test)` the tests below are callers and the
// lint would not fire — an unconditional `expect` is itself unfulfilled there.
#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "B2 (AlignmentFile::open) consumes these; remove this attribute with it"
    )
)]

use noodles_sam as sam;

use crate::fasta::{ContigEntry, ContigList};

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
