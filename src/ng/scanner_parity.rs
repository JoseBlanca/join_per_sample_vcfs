//! **Golden-catalog parity for ng's tandem-repeat scanner.**
//!
//! [`find_tandem_repeats`](crate::ng::tandem_repeat::find_tandem_repeats) must
//! reproduce the **final catalog** the external `trf-mod` produces, run through
//! production's *unchanged* post-filter (`ssr::catalog::postprocess`).
//!
//! The oracle is a golden catalog snapshot committed under
//! `tests/data/tandem_repeat/` — built once by the current `trf-mod` →
//! `postprocess` path (`ssr-catalog`) on a synthetic STR-diversity reference. This
//! test bypasses `trf-mod`: it runs the scanner on the same reference, bridges the
//! intervals into the post-filter's input via the test-only `TrfRecord::for_test`,
//! and asserts the resulting `Locus` set reproduces the golden one.
//! **Production stays on `trf-mod`** — this only *validates* ng's scanner; nothing
//! here changes the shipping path, and it reads `src/ssr/` without writing to it.
//!
//! ## Why this lives in `src/ng/`, having been written in `src/ssr/catalog/`
//!
//! It is **ng's test**: its subject is ng's scanner, and production is only the
//! yardstick. It began in `src/ssr/catalog/` because that is where the golden path
//! was, and that made `src/ssr/` — which is **frozen** (owner, 2026-07-16: *"leave
//! production as is"*) — the **only** part of the tree that `use`d `crate::ng`.
//!
//! That dependency pointed the wrong way, and it bit: typed-regions Milestone B2
//! widens ng's coordinates to `u64`, and `RepeatInterval`'s `u32`-ness was baked
//! into this file (four sites, most bindingly `TrfRecord::for_test(start: u32,
//! end: u32, …)`) — so an **ng-internal** change could not compile without editing
//! frozen production. A freeze that ng can break is not a freeze; the whole point
//! is that an experiment cannot destabilise the yardstick it is scored against.
//!
//! Moving it (owner-approved) severs that: **production now depends on nothing in
//! ng**, permanently, and this test follows ng's types as they move. The move is
//! test-only — no shipped behaviour changed, and the parity bar is the same one.
//!
//! ## The catalog pre-filter (a validated integration finding)
//!
//! `trf-mod` hands `postprocess` a **clean** candidate set — statistically
//! significant repeats, redundancy already eliminated. The raw scanner is
//! deliberately permissive (`min_copies = 2`), so in aperiodic sequence it emits
//! many low-copy noise repeats and every period-multiple of a real tract. Fed
//! straight to `build_loci`, that noise triggers `drop_bundles` (which runs
//! *before* the copy-number floor) and cascades away the real loci. So a consumer
//! must apply, **before** `build_loci`, the same two cleanups `trf-mod` bakes in —
//! the per-period copy floor and period-multiple redundancy elimination (its
//! `IsRedundant`).
//!
//! [`catalog_prefilter`] is that policy, kept here **verbatim as trf-mod's
//! shape**. Note ng has since ported it properly, with the floors as a real knob,
//! to [`crate::ng::region_typing::segment_criteria::prefilter`] — this copy stays frozen
//! against the golden path rather than tracking ng's, so that the two oracles stay
//! independent: if ng's prefilter drifts, `classification`'s own differential says so,
//! and this one keeps measuring the scanner against trf-mod.

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use crate::ng::tandem_repeat::{PeriodRange, RepeatInterval, ScanParams, find_tandem_repeats};
use crate::ssr::catalog::CatalogParams;
use crate::ssr::catalog::io::CatalogReader;
use crate::ssr::catalog::postprocess::build_loci;
use crate::ssr::catalog::trf::TrfRecord;
use crate::ssr::types::Locus;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("tandem_repeat")
        .join(name)
}

/// The per-period minimum copy count `postprocess` applies (as its `copy_number_floor`),
/// used here as a **pre-`build_loci`** filter — the catalog's job, since `drop_bundles` runs
/// before the post-filter's own floor.
fn copy_floor(period: u8) -> u32 {
    match period {
        2 => 5,
        3 => 4,
        _ => 3,
    }
}

/// The catalog's pre-`build_loci` cleanup of raw scanner intervals (see the module docs):
/// keep period-2..6 intervals clearing the per-period copy floor, then drop period-multiple
/// redundancy (a higher-period interval overlapping a divisor-period one is the same tract).
fn catalog_prefilter(intervals: &[RepeatInterval]) -> Vec<RepeatInterval> {
    let mut floored: Vec<RepeatInterval> = intervals
        .iter()
        .copied()
        .filter(|iv| {
            iv.period >= 2
                && (iv.end - iv.start) / u64::from(iv.period) >= u64::from(copy_floor(iv.period))
        })
        .collect();
    // Process low periods first so a fundamental tract is kept and its multiples dropped.
    floored.sort_by_key(|iv| (iv.period, iv.start));
    let mut kept: Vec<RepeatInterval> = Vec::new();
    for iv in floored {
        let redundant = kept.iter().any(|a| {
            a.period < iv.period
                && iv.period % a.period == 0
                && iv.start < a.end
                && a.start < iv.end
        });
        if !redundant {
            kept.push(iv);
        }
    }
    kept
}

/// Run the scanner path (`find_tandem_repeats` → catalog pre-filter → `build_loci`) over the
/// reference the golden was built from, and return its loci.
fn scanner_loci(reference: &Path, params: &CatalogParams) -> Vec<Locus> {
    let file = File::open(reference).unwrap();
    let mut reader = noodles_fasta::io::Reader::new(BufReader::new(file));
    let mut out = Vec::new();
    for result in reader.records() {
        let rec = result.unwrap();
        let name = String::from_utf8_lossy(rec.name()).into_owned();
        let seq = rec.sequence().as_ref();
        // The catalog scans period 1..=6; the pre-filter and post-filter drop period 1.
        let intervals =
            find_tandem_repeats(seq, PeriodRange::new(1, 6).unwrap(), &ScanParams::default());
        // **The narrowing that made this file block B2, now on ng's side of the
        // fence.** `TrfRecord` is trf-mod's parse shape and is `u32`; ng's
        // `RepeatInterval` is `u64` (spec §4). While this test lived in
        // `src/ssr/catalog/`, widening ng meant editing frozen production — which is
        // why it moved here. Now the conversion is ng's own, at ng's own boundary,
        // and production never notices.
        //
        // PANIC-FREE: the synthetic parity fixture is a few kb; a `u32` holds 4.29 Gb.
        // `expect` over `as` because a truncated coordinate would silently compare the
        // *wrong tract* against the golden catalog — a green test asserting nothing.
        let recs: Vec<TrfRecord> = catalog_prefilter(&intervals)
            .iter()
            .map(|iv| {
                TrfRecord::for_test(
                    u32::try_from(iv.start).expect("parity fixture coordinates fit u32"),
                    u32::try_from(iv.end).expect("parity fixture coordinates fit u32"),
                    u16::from(iv.period),
                    iv.score,
                    b"",
                )
            })
            .collect();
        out.extend(build_loci(recs, &name, seq, params));
    }
    out
}

/// Two loci are the *same tract* if they are on the same contig and their `[start, end)`
/// spans overlap — the parity-relevant match, tolerant of the small boundary/phase wobble
/// two different detectors place on the same repeat (spec §6.1). `(chrom, start, end)` for a
/// readable diff.
fn overlaps(a: &Locus, b: &Locus) -> bool {
    a.chrom() == b.chrom() && a.start() < b.end() && b.start() < a.end()
}

fn describe(l: &Locus) -> String {
    format!(
        "{}:{}-{} {}",
        l.chrom(),
        l.start(),
        l.end(),
        String::from_utf8_lossy(l.motif().as_bytes())
    )
}

#[test]
fn scanner_reproduces_the_trf_mod_golden_catalog() {
    // The golden catalog (trf-mod → postprocess) and its build parameters.
    let mut golden_reader =
        CatalogReader::new(File::open(fixture("golden.ssr_catalog.bed.gz")).unwrap()).unwrap();
    let params = golden_reader.header().params.clone();
    let golden = golden_reader.read_all().unwrap();

    // The scanner path over the same reference, with the same post-filter parameters.
    let scanner = scanner_loci(&fixture("synthetic_ref.fa"), &params);

    // Recall: every golden tract must be reproduced (overlap match). Track exact-coordinate
    // hits vs boundary/phase wobble, and scanner-only loci, so the diff stays reviewed.
    let mut exact = 0usize;
    let mut wobble = Vec::new();
    let mut missed = Vec::new();
    for g in &golden {
        match scanner.iter().find(|s| overlaps(g, s)) {
            Some(s) if (s.start(), s.end(), s.motif()) == (g.start(), g.end(), g.motif()) => {
                exact += 1;
            }
            Some(s) => wobble.push(format!("{}  ~  {}", describe(g), describe(s))),
            None => missed.push(describe(g)),
        }
    }
    let extra: Vec<String> = scanner
        .iter()
        .filter(|s| !golden.iter().any(|g| overlaps(g, s)))
        .map(describe)
        .collect();

    let recovered = golden.len() - missed.len();
    let recall = recovered as f64 / golden.len() as f64;

    assert!(
        recall >= 0.99,
        "scanner recall {recall:.4} ({recovered}/{}) below the 0.99 parity bar (spec §9.2).\n\
         missed golden loci: {missed:#?}",
        golden.len(),
    );

    // Reviewed diff (spec §6.1): a handful of boundary/phase wobbles and scanner-only loci are
    // expected between two different detectors. Printed (visible with `--nocapture`), and the
    // scanner-only loci must be genuine STRs — pinned by count so the diff can't silently grow.
    eprintln!(
        "parity: {}/{} recall — {exact} exact, {} boundary/phase wobble, {} scanner-only.\n\
         wobble: {wobble:#?}\n  scanner-only (genuine STRs trf-mod's significance model rejected): {extra:#?}",
        recovered,
        golden.len(),
        wobble.len(),
        extra.len(),
    );
    assert!(
        extra.len() <= 1,
        "more scanner-only loci ({}) than the reviewed baseline of 1 (the ctg2 CCG tract): {extra:#?}",
        extra.len()
    );
}
