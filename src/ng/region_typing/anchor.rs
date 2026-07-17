//! ng step 3 — the typed-region generator, end to end (impl plan E3, spec §8).
//!
//! **The port anchor.** Everything else in step 3 is a unit test standing next to the code
//! it checks, driving an `InMemoryRefSeq` and reaching in at `partition_resident` or
//! `partition_windowed`. This module drives the shipping stack instead: a **real
//! multi-contig FASTA on disk**, read through the reference implementation that ships
//! (`WindowedRefSeq` — file-backed, `.fai`-indexed, evicting), through the public surface
//! (`TypedRegionIterator`), and nothing else.
//!
//! **`#[cfg(test)]` in-crate, not `tests/`, for two reasons.** The golden catalog's reader
//! is `pub(crate)` and production is frozen (spec Revision), so an out-of-crate test
//! cannot open the oracle at all; and step 3 must touch no file outside `src/ng/`. It is
//! the `scanner_parity` shape exactly — an ng test that needs production as a yardstick,
//! living in ng (B2, commit `d097ebf`). It still consumes only what a caller could.
//!
//! It exists to fail when the pieces stop fitting together, which no unit test can:
//! `RawChromReader`'s windowed reads, the contig table's provenance, eviction, the
//! iterator's ownership, the scan/emit split, and the walk itself, all at once, on
//! sequence nobody wrote to make a point.
//!
//! What it asserts (spec §8):
//!
//! - **`.cat` parity** — the walk reproduces the committed trf-mod-built golden catalog;
//! - the **partition invariant** — contiguous, non-overlapping, complete, maximal;
//! - **window-invariance** — `window_bp` is a memory knob;
//! - **BED-invariance** — a BED chooses what you are shown, never what things are;
//! - the **edge cases** §8 lists: a tract at position 1, a repeat-free contig, and a tract
//!   at one contig's end abutting one at the next's start.

use crate::fasta::{ContigEntry, ContigList};
use crate::ng::WindowedRefSeq;
use crate::ng::region_typing::{
    GenomeRegions, RegionKind, TypedRegion, TypedRegionConfig, TypedRegionIterator,
};
use crate::ng::types::{Bp, ContigId, Position};
use crate::regions::ContigBounds;
use std::io::Write;
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------
// Fixtures on disk
// ---------------------------------------------------------------------

/// Write `contigs` as a FASTA plus its `.fai`, **one line per contig** so the index is
/// arithmetic rather than a guess: a record's bases start right after its header and run
/// unbroken, so `offset`, `linebases` and `linewidth` fall out of the lengths.
///
/// Returns the FASTA path and the contig table. The `TempDir` must outlive them, which is
/// why it comes back too.
fn write_fasta(contigs: &[(&str, Vec<u8>)]) -> (tempfile::TempDir, PathBuf, ContigList) {
    std::fs::create_dir_all("tmp").expect("project-local scratch (CLAUDE.md: never /tmp)");
    let dir = tempfile::tempdir_in("tmp").unwrap();
    let fa = dir.path().join("ref.fa");
    let fai = dir.path().join("ref.fa.fai");

    let mut fasta = std::fs::File::create(&fa).unwrap();
    let mut index = std::fs::File::create(&fai).unwrap();
    let mut offset = 0u64;
    for (name, bases) in contigs {
        let header = format!(">{name}\n");
        fasta.write_all(header.as_bytes()).unwrap();
        fasta.write_all(bases).unwrap();
        fasta.write_all(b"\n").unwrap();

        let bases_at = offset + header.len() as u64;
        writeln!(
            index,
            "{}\t{}\t{}\t{}\t{}",
            name,
            bases.len(),
            bases_at,
            bases.len(),
            bases.len() + 1
        )
        .unwrap();
        offset = bases_at + bases.len() as u64 + 1;
    }

    let table = ContigList {
        entries: contigs
            .iter()
            .map(|(name, bases)| ContigEntry {
                name: (*name).to_string(),
                length: bases.len() as u64,
                md5: None,
            })
            .collect(),
    };
    (dir, fa, table)
}

fn bounds(table: &ContigList) -> Vec<ContigBounds<'_>> {
    table
        .entries
        .iter()
        .map(|e| ContigBounds {
            name: &e.name,
            length: e.length as u32,
        })
        .collect()
}

/// The committed golden reference, as `(name, bases)` — real sequence, and the same file
/// the golden `.cat` catalog was built from.
fn golden_contigs() -> Vec<(String, Vec<u8>)> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("tandem_repeat")
        .join("synthetic_ref.fa");
    let file = std::fs::File::open(path).unwrap();
    let mut reader = noodles_fasta::io::Reader::new(std::io::BufReader::new(file));
    reader
        .records()
        .map(|r| {
            let rec = r.unwrap();
            (
                String::from_utf8_lossy(rec.name()).into_owned(),
                rec.sequence().as_ref().to_vec(),
            )
        })
        .collect()
}

/// Walk a reference **through the shipping stack**: a file-backed `WindowedRefSeq` and the
/// public iterator. Every assertion in this file goes through here.
fn walk(
    fasta: &Path,
    table: &ContigList,
    spans: GenomeRegions,
    config: TypedRegionConfig,
) -> Vec<TypedRegion> {
    let reference = WindowedRefSeq::new(fasta.to_path_buf(), table.clone());
    TypedRegionIterator::over_regions(reference, spans, config)
        .expect("the spans name contigs this reference has")
        .collect::<Result<_, _>>()
        .expect("the reference reads cleanly")
}

// ---------------------------------------------------------------------
// The invariant, stated here rather than borrowed
// ---------------------------------------------------------------------

/// **The partition invariant** (spec §2.3), over one contig's slice of the output:
/// contiguous, non-overlapping, complete over `[1, contig_len]`, and maximal.
///
/// Written out again rather than shared with the unit tests, deliberately: this file is
/// the outside view, and an anchor that imports the thing it is anchoring proves less. One
/// property — *concatenating the regions reconstructs what was asked for, exactly.*
#[track_caller]
fn assert_partitions(regions: &[TypedRegion], contig: ContigId, contig_len: u64, case: &str) {
    assert!(
        !regions.is_empty(),
        "{case}: a non-empty contig has regions"
    );
    let mut expect = 1u64;
    let mut prev: Option<std::mem::Discriminant<RegionKind>> = None;
    for r in regions {
        assert_eq!(r.region.contig, contig, "{case}: contig");
        assert_eq!(
            r.region.start.get(),
            expect,
            "{case}: gap or overlap at {}",
            r.region.start.get()
        );
        assert!(r.region.end >= r.region.start, "{case}: empty region");
        let kind = std::mem::discriminant(&r.kind);
        assert_ne!(
            Some(kind),
            prev,
            "{case}: two consecutive regions share a kind at {} — MAXIMALITY. For Generic \
             that is a correctness bug: the pileup mints loci INSIDE a generic region, so a \
             split run makes an indel across the join callable by neither half",
            r.region.start.get()
        );
        prev = Some(kind);
        expect = r.region.end.get() + 1;
    }
    assert_eq!(
        expect - 1,
        contig_len,
        "{case}: the partition must cover exactly [1, {contig_len}] — COMPLETENESS"
    );
}

fn kinds(regions: &[TypedRegion]) -> Vec<&'static str> {
    regions
        .iter()
        .map(|r| match &r.kind {
            RegionKind::SsrLocus(_) => "locus",
            RegionKind::SsrBundle { .. } => "bundle",
            RegionKind::Generic => "generic",
            RegionKind::Satellite => "satellite",
        })
        .collect()
}

fn on_contig(regions: &[TypedRegion], contig: ContigId) -> Vec<TypedRegion> {
    regions
        .iter()
        .filter(|r| r.region.contig == contig)
        .cloned()
        .collect()
}

// ---------------------------------------------------------------------
// The real reference, end to end
// ---------------------------------------------------------------------

/// **`.cat` parity through the whole stack** (spec §8.1) — the anchor's headline.
///
/// The walk at the catalog's settings must reproduce the catalog: every golden locus is
/// present, **or** absent *and* inside a satellite run. A strict subset, and that shape is
/// earned by spec §2.4's ordering — the cap applies to the *cleaned* coverage, after
/// admission, so the difference can only go one way.
///
/// The oracle is the committed **trf-mod-built** golden catalog: a different detector, a
/// different code path, nothing ng touched. Overlap matching, inherited from
/// `scanner_parity`, because the detector difference is characterised (±1–2 bp of
/// boundary/phase wobble) and is a yardstick rather than a confound.
///
/// D1 asserts this too, resident and in-process. What this adds is everything between:
/// the FASTA on disk, `RawChromReader`'s windowed reads, the contig table, eviction, and
/// the iterator.
#[test]
fn the_walk_reproduces_the_golden_catalog_through_the_shipping_stack() {
    use crate::ssr::catalog::io::CatalogReader;

    let fixture = |name: &str| {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("data")
            .join("tandem_repeat")
            .join(name)
    };
    let mut reader =
        CatalogReader::new(std::fs::File::open(fixture("golden.ssr_catalog.bed.gz")).unwrap())
            .unwrap();
    let cat_params = reader.header().params.clone();
    let golden = reader.read_all().unwrap();
    assert!(!golden.is_empty(), "the golden catalog must have loci");

    let contigs = golden_contigs();
    assert!(
        contigs.len() > 1,
        "the anchor walks a MULTI-contig reference"
    );
    let owned: Vec<(&str, Vec<u8>)> = contigs
        .iter()
        .map(|(n, b)| (n.as_str(), b.clone()))
        .collect();
    let (_dir, fa, table) = write_fasta(&owned);

    // The catalog's settings, pinned explicitly rather than to whatever `Default` is —
    // or this starts failing the first time someone moves a floor and reads a *result*
    // as a bug.
    let config = TypedRegionConfig {
        admission: crate::ng::region_typing::admission::SsrAdmissionParams {
            min_purity: cat_params.min_purity,
            min_score: cat_params.min_score,
            flank_bp: u64::from(cat_params.flank_bp),
            ..Default::default()
        },
        ..TypedRegionConfig::default()
    };

    let regions = walk(
        &fa,
        &table,
        GenomeRegions::whole_contigs(&bounds(&table)),
        config,
    );

    // The partition holds on every contig of real sequence.
    for (idx, (name, bases)) in contigs.iter().enumerate() {
        let contig = ContigId(idx as u32);
        assert_partitions(
            &on_contig(&regions, contig),
            contig,
            bases.len() as u64,
            &format!("golden contig {name}"),
        );
    }

    // Present, or absent AND inside a satellite.
    let name_of = |id: ContigId| contigs[id.get() as usize].0.clone();
    let ours: Vec<(String, u64, u64)> = regions
        .iter()
        .filter_map(|r| match &r.kind {
            RegionKind::SsrLocus(l) => Some((name_of(r.region.contig), l.start(), l.end())),
            _ => None,
        })
        .collect();
    let satellites: Vec<(String, u64, u64)> = regions
        .iter()
        .filter(|r| matches!(r.kind, RegionKind::Satellite))
        .map(|r| {
            (
                name_of(r.region.contig),
                r.region.start.get(),
                r.region.end.get(),
            )
        })
        .collect();
    assert!(!ours.is_empty(), "the walk must find loci");

    // Production's `Locus` is 0-based half-open, ng's 1-based inclusive: `[s, e)` is
    // `[s + 1, e]` (spec §4).
    let overlaps =
        |a: &(String, u64, u64), b: &(String, u64, u64)| a.0 == b.0 && a.1 <= b.2 && b.1 <= a.2;
    let mut missed = Vec::new();
    for g in &golden {
        let g1 = (
            g.chrom().to_string(),
            u64::from(g.start()) + 1,
            u64::from(g.end()),
        );
        if !ours.iter().any(|o| overlaps(&g1, o)) && !satellites.iter().any(|s| overlaps(&g1, s)) {
            missed.push(format!("{}:{}-{}", g1.0, g1.1, g1.2));
        }
    }
    assert!(
        missed.is_empty(),
        "every golden locus must be present, or absent AND inside a satellite run. At the \
         catalog's settings, a locus missing for any other reason is a machinery bug. \
         Missing: {missed:#?}"
    );
}

/// **Window-invariance through the shipping stack** (spec §2.3): `window_bp` is a memory
/// knob and must not move the output — including through a *file-backed, evicting*
/// reference, where a window boundary is also a buffer boundary and a re-read.
///
/// The unit tests prove this against an in-memory reference, where a window edge costs
/// nothing. Here the reference slides and evicts underneath the walk, which is the one
/// place a window edge could genuinely change an answer.
#[test]
fn window_bp_does_not_change_the_output_through_a_file_backed_reference() {
    let contigs = golden_contigs();
    let owned: Vec<(&str, Vec<u8>)> = contigs
        .iter()
        .map(|(n, b)| (n.as_str(), b.clone()))
        .collect();
    let (_dir, fa, table) = write_fasta(&owned);
    let spans = || GenomeRegions::whole_contigs(&bounds(&table));

    let baseline = walk(&fa, &table, spans(), TypedRegionConfig::default());
    assert!(!baseline.is_empty());

    for window_bp in [200u64, 700, 100_000] {
        let config = TypedRegionConfig {
            window_bp: Bp(window_bp),
            ..TypedRegionConfig::default()
        };
        assert_eq!(
            walk(&fa, &table, spans(), config),
            baseline,
            "window_bp = {window_bp} changed the output"
        );
    }
}

/// **BED-invariance through the shipping stack** (spec §2.5): the BED chooses what you are
/// shown, never what things are.
#[test]
fn a_bed_does_not_change_what_things_are_through_the_shipping_stack() {
    let contigs = golden_contigs();
    let owned: Vec<(&str, Vec<u8>)> = contigs
        .iter()
        .map(|(n, b)| (n.as_str(), b.clone()))
        .collect();
    let (_dir, fa, table) = write_fasta(&owned);

    let whole = walk(
        &fa,
        &table,
        GenomeRegions::whole_contigs(&bounds(&table)),
        TypedRegionConfig::default(),
    );

    // A BED over the middle of the first contig, written 0-based half-open as BED is.
    let dir = tempfile::tempdir_in("tmp").unwrap();
    let bed = dir.path().join("r.bed");
    let (first, len) = (&contigs[0].0, contigs[0].1.len() as u64);
    assert!(len > 800, "the fixture contig must be worth subsetting");
    writeln!(std::fs::File::create(&bed).unwrap(), "{first}\t400\t800").unwrap();

    let spans = GenomeRegions::from_bed_path(&bed, &bounds(&table)).expect("valid bed");
    let got = walk(&fa, &table, spans, TypedRegionConfig::default());
    assert!(!got.is_empty());

    for r in &got {
        for pos in [r.region.start.get(), r.region.end.get()] {
            // An object straddling the edge is emitted whole, so only ask about the
            // bases actually requested (1-based [401, 800]).
            if !(401..=800).contains(&pos) {
                continue;
            }
            let truth = whole
                .iter()
                .find(|w| w.region.contig == r.region.contig && w.region.contains(Position(pos)))
                .expect("the whole-genome run covers every base");
            assert_eq!(
                std::mem::discriminant(&r.kind),
                std::mem::discriminant(&truth.kind),
                "base {pos} is {:?} with a BED and {:?} without one",
                r.kind,
                truth.kind
            );
        }
    }
}

// ---------------------------------------------------------------------
// The edge cases spec §8 lists
// ---------------------------------------------------------------------

/// Aperiodic filler — **not** a homopolymer, which would be a period-1 tract and make
/// these fixtures pass for the wrong reason (the lesson D1 learned the hard way).
fn filler(n: usize) -> Vec<u8> {
    b"ACGTTGCAAGCTTGCA"
        .iter()
        .copied()
        .cycle()
        .take(n)
        .collect()
}

fn tract(copies: usize) -> Vec<u8> {
    b"AT".iter().copied().cycle().take(copies * 2).collect()
}

/// **The three edge cases spec §8 names, on one multi-contig reference** — and the third
/// is the one only a multi-contig walk can reach.
///
/// | contig | shape | what it pins |
/// |---|---|---|
/// | `at_start` | a tract at base 1, then filler | no left flank at the CONTIG's start → `Generic`, and the partition still starts at 1 |
/// | `empty_of_repeats` | aperiodic filler only | **one** `Generic`, not a run of them (maximality) |
/// | `ends_with_tract` / `starts_with_tract` | a tract at one contig's end abutting one at the next's start | contigs are independent: neither tract borrows the other's flank, and nothing carries across the seam |
/// | `flanked` | **the control**: the same tract, flanks both sides | it IS a locus |
///
/// **The control is what makes the rest mean anything.** Three of these assert *"not a
/// locus"*, and a tract that was never admissible would satisfy all three for free — which
/// is exactly how D1's satellite test passed for the wrong reason. `flanked` is the same
/// `tract(10)` built by the same helper at the same settings: it comes back a locus, so the
/// absences above are the contig **ends** doing their job and not the tract being
/// unremarkable.
#[test]
fn the_edge_cases_hold_on_a_multi_contig_reference() {
    let mut at_start = tract(10);
    at_start.extend(filler(200));

    let repeat_free = filler(300);

    let mut ends_with = filler(200);
    ends_with.extend(tract(10));

    let mut starts_with = tract(10);
    starts_with.extend(filler(200));

    let mut flanked = filler(200);
    flanked.extend(tract(10));
    flanked.extend(filler(200));

    let contigs: Vec<(&str, Vec<u8>)> = vec![
        ("at_start", at_start.clone()),
        ("empty_of_repeats", repeat_free.clone()),
        ("ends_with_tract", ends_with.clone()),
        ("starts_with_tract", starts_with.clone()),
        ("flanked", flanked.clone()),
    ];
    let (_dir, fa, table) = write_fasta(&contigs);
    let regions = walk(
        &fa,
        &table,
        GenomeRegions::whole_contigs(&bounds(&table)),
        TypedRegionConfig::default(),
    );

    for (idx, (name, bases)) in contigs.iter().enumerate() {
        assert_partitions(
            &on_contig(&regions, ContigId(idx as u32)),
            ContigId(idx as u32),
            bases.len() as u64,
            name,
        );
    }

    // A tract at base 1 has no left flank to anchor against, so it is not a locus — and
    // the partition still tiles from base 1.
    let first = on_contig(&regions, ContigId(0));
    assert_eq!(first[0].region.start, Position(1));
    assert!(
        !kinds(&first).contains(&"locus"),
        "a tract at base 1 has no left flank: {:?}",
        kinds(&first)
    );

    // A repeat-free contig is exactly ONE Generic region. Not many.
    let second = on_contig(&regions, ContigId(1));
    assert_eq!(kinds(&second), vec!["generic"]);
    assert_eq!(second[0].region.end, Position(repeat_free.len() as u64));

    // **The seam.** A tract at one contig's end and one at the next's start are 20 bp
    // apart *in the file* and on different chromosomes: neither may borrow the other's
    // flank, bundle with it, or carry a coverage run across. Both are dropped for the
    // same reason — each abuts its own contig's end — and the two partitions are
    // independent.
    let third = on_contig(&regions, ContigId(2));
    let fourth = on_contig(&regions, ContigId(3));
    assert!(
        !kinds(&third).contains(&"locus"),
        "the tract at this contig's END has no right flank: {:?}",
        kinds(&third)
    );
    assert!(
        !kinds(&fourth).contains(&"locus"),
        "the tract at this contig's START has no left flank: {:?}",
        kinds(&fourth)
    );
    assert!(
        !kinds(&third).contains(&"bundle") && !kinds(&fourth).contains(&"bundle"),
        "and they must not bundle with each other ACROSS the contig seam"
    );
    assert_eq!(
        third.last().unwrap().region.end,
        Position(ends_with.len() as u64),
        "the third contig's partition ends at its own last base"
    );
    assert_eq!(
        fourth[0].region.start,
        Position(1),
        "and the fourth's starts at its own first base"
    );

    // **The control.** The same tract, the same settings, flanks either side: a locus. So
    // every "not a locus" above is the contig's end doing its job, and not a tract that was
    // never going to be admitted (D1's satellite test taught this the hard way).
    let fifth = on_contig(&regions, ContigId(4));
    assert_eq!(
        kinds(&fifth),
        vec!["generic", "locus", "generic"],
        "the SAME tract, given flanks, is a locus"
    );
}

/// The running tally describes what came out, on the real reference — spec §3.1's "no
/// silent caps", checked against the regions themselves rather than against literals.
#[test]
fn the_counts_describe_the_walk_of_a_real_reference() {
    let contigs = golden_contigs();
    let owned: Vec<(&str, Vec<u8>)> = contigs
        .iter()
        .map(|(n, b)| (n.as_str(), b.clone()))
        .collect();
    let (_dir, fa, table) = write_fasta(&owned);

    let reference = WindowedRefSeq::new(fa.clone(), table.clone());
    let mut iter = TypedRegionIterator::over_regions(
        reference,
        GenomeRegions::whole_contigs(&bounds(&table)),
        TypedRegionConfig::default(),
    )
    .expect("valid spans");

    let mut regions = Vec::new();
    for r in iter.by_ref() {
        regions.push(r.expect("the reference reads cleanly"));
    }
    let counts = iter.counts();

    assert_eq!(counts.spans, contigs.len() as u64, "one span per contig");
    let count = |k: &str| kinds(&regions).iter().filter(|x| **x == k).count() as u64;
    assert_eq!(counts.ssr_loci, count("locus"));
    assert_eq!(counts.ssr_bundles, count("bundle"));
    assert_eq!(counts.generic, count("generic"));
    assert_eq!(counts.satellites, count("satellite"));
    assert!(counts.ssr_loci > 0, "real sequence, real loci");

    // Repeat coverage that yielded no locus is a *subset* of the repeat coverage, so it
    // cannot exceed the bases the walk typed as repeat at all.
    let repeat_bp: u64 = regions
        .iter()
        .filter(|r| !matches!(r.kind, RegionKind::Generic))
        .map(|r| r.region.len())
        .sum();
    assert!(
        counts.repeat_bp_with_no_locus <= repeat_bp + counts.ssr_loci,
        "the no-locus gap ({}) cannot exceed the repeat coverage it is part of ({repeat_bp})",
        counts.repeat_bp_with_no_locus
    );
}
