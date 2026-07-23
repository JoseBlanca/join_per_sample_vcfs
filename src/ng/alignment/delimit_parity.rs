//! **Byte-parity of ng's tract delimiter against production's — the port anchor.**
//!
//! Algorithm 3 ([`SsrFlatGapAligner`](super::ssr_best_path_flat_gap::SsrFlatGapAligner)) is
//! a port of production's `delimit_read`, and it is the **only aligner in this module with a
//! parity oracle**: algorithms 2, 4, 5 and 6 are *measured*, not verified (spec §10.3, arch
//! §Test & bench shape). So this test is what the whole port rests on, and it is what makes
//! Milestone C's banding a *performance* change rather than a behaviour change — banding
//! must leave every assertion here untouched.
//!
//! # What parity means here, exactly
//!
//! Production answers with a two-case `Delimited`: a `Region(range)` or a side-blind
//! `BorderOffEnd`. ng answers with the four-case [`RepeatSpan`]. The two are **not** the same
//! shape, because the widening is the point of step B2 — so parity is asserted where the two
//! are comparable and *counted* where they are not:
//!
//! - **Every read production calls `Region` must measure the same bytes here.** Not the same
//!   *classification* — the same **offsets**, which is the measurement itself. This is the
//!   hard assertion.
//! - **Every read production calls `BorderOffEnd` must land in a one-flank or no-flank
//!   case.** It cannot be checked against production more finely than that, because
//!   production discarded the information: `BorderOffEnd` does not say which flank was
//!   missing. So the reclassification is checked **by count and by shape**, which is all the
//!   oracle can support (plan step B3).
//!
//! # Why this lives in `src/ng/` and reads `src/ssr/`
//!
//! It is **ng's test**: its subject is ng's aligner, and production is only the yardstick.
//! `src/ssr/` is frozen (owner, 2026-07-16) and is read here, never written — the same shape
//! as [`scanner_parity`](crate::ng::scanner_parity), and the same shape as the `MAX_SLIP`
//! cross-check in [`stutter`](super::stutter). It is `#[cfg(test)]`, so **shipping ng code
//! still depends on nothing in production**.

use super::ssr_best_path_flat_gap::{SsrFlatGapAligner, ViterbiScratch};
use super::{
    BestPathAligner, PerQualityEmission, ReadBases, RepeatContext, RepeatGeometry, RepeatSpan,
};
use super::{StutterModel, stutter};
use crate::ng::types::Bp;
use crate::ng::types::Motif as NgMotif;
use crate::ssr::pileup::alignment::{
    Delimited, HmmModel, ViterbiScratch as ProdScratch, delimit_read,
};
use crate::ssr::types::{Locus, Motif as ProdMotif};

/// A deterministic PRNG, so a failure is reproducible from its seed alone.
///
/// SplitMix64 — the same generator ng already ported for the locus reservoir. Written out
/// rather than pulled from a crate because a parity fixture must be reproducible from the
/// source in front of you, not from a dependency's version.
struct SplitMix64(u64);

impl SplitMix64 {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn below(&mut self, bound: usize) -> usize {
        (self.next_u64() % bound as u64) as usize
    }

    fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4)]
    }
}

/// One generated case: a locus frame and a read to delimit against it.
struct Case {
    left_flank: Vec<u8>,
    tract: Vec<u8>,
    right_flank: Vec<u8>,
    motif: Vec<u8>,
    read: Vec<u8>,
    quality: Vec<u8>,
}

impl Case {
    fn reference(&self) -> Vec<u8> {
        let mut reference = self.left_flank.clone();
        reference.extend_from_slice(&self.tract);
        reference.extend_from_slice(&self.right_flank);
        reference
    }
}

/// Build a case: a random motif tiled into a tract, random unique flanks, and a read derived
/// from the reference by a random length change, substitutions and end truncations.
///
/// The generator deliberately produces the four things the plan's fixture must carry — a
/// clean repeat, a long allele, an interrupted repeat, and reads running off each end — plus
/// the awkward middles between them.
fn generate(rng: &mut SplitMix64) -> Case {
    let period = 1 + rng.below(6);
    let motif: Vec<u8> = (0..period).map(|_| rng.base()).collect();
    let units = 2 + rng.below(10);
    let tract: Vec<u8> = motif.iter().cycle().take(period * units).copied().collect();

    // Flanks are random sequence, and their two lengths are drawn independently — equal
    // flank lengths hid a transposition in this module once already.
    //
    // **Zero-length flanks are generated on purpose**, about one case in eight per side. They
    // are the geometry at a contig edge, and they are the one place ng *deliberately*
    // disagrees with production — so a fixture that never produced them would leave the
    // widening's most consequential divergence untested, and would let the harness look
    // greener than it is.
    let left_len = if rng.below(8) == 0 {
        0
    } else {
        1 + rng.below(12)
    };
    let right_len = if rng.below(8) == 0 {
        0
    } else {
        1 + rng.below(12)
    };
    let left_flank: Vec<u8> = (0..left_len).map(|_| rng.base()).collect();
    let right_flank: Vec<u8> = (0..right_len).map(|_| rng.base()).collect();

    // The read starts as the reference, then gains or loses tract units, picks up
    // substitutions, and may be truncated at either end.
    let mut read = left_flank.clone();
    let unit_change = rng.below(9) as isize - 4;
    let read_units = (units as isize + unit_change).max(0) as usize;
    read.extend(motif.iter().cycle().take(period * read_units).copied());
    read.extend_from_slice(&right_flank);

    // An out-of-frame indel, sometimes — an interrupted or ragged repeat.
    if !read.is_empty() && rng.below(4) == 0 {
        let at = rng.below(read.len());
        if rng.below(2) == 0 {
            read.insert(at, rng.base());
        } else {
            read.remove(at);
        }
    }
    // Substitutions.
    for _ in 0..rng.below(3) {
        if !read.is_empty() {
            let at = rng.below(read.len());
            read[at] = rng.base();
        }
    }
    // Truncations, which is how reads run off their own ends inside a repeat.
    if rng.below(3) == 0 {
        let cut = rng.below(read.len().max(1) / 2 + 1);
        read.drain(..cut.min(read.len()));
    }
    if rng.below(3) == 0 {
        let keep = read
            .len()
            .saturating_sub(rng.below(read.len().max(1) / 2 + 1));
        read.truncate(keep);
    }

    let quality: Vec<u8> = (0..read.len()).map(|_| (rng.below(41)) as u8).collect();

    Case {
        left_flank,
        tract,
        right_flank,
        motif,
        read,
        quality,
    }
}

/// Run production's delimiter on a case.
fn production_answer(case: &Case) -> Delimited {
    let reference = case.reference();
    let tract_start = case.left_flank.len() as u32;
    let tract_end = tract_start + case.tract.len() as u32;
    let locus = Locus::new(
        "parity".into(),
        tract_start,
        tract_end,
        ProdMotif::new(&case.motif).expect("a valid motif"),
        1.0,
        reference.into_boxed_slice(),
        0,
    )
    .expect("a valid locus");

    let model = HmmModel::new();
    let mut scratch = ProdScratch::new();
    delimit_read(&case.read, &case.quality, &locus, &model, &mut scratch)
}

/// Run ng's aligner on the same case.
///
/// **The scratch is passed in, not created here.** Hoisting it out of the loop is what makes
/// this the one randomized, size-varying driver that exercises `ViterbiScratch`'s
/// grow-but-never-clear reuse — the hazard `ssr_best_path_flat_gap` names as *the* Milestone-C
/// risk, since banding will stop the fill being exhaustive. A fresh scratch per case would be
/// structurally blind to it.
fn ng_answer(case: &Case, scratch: &mut ViterbiScratch) -> RepeatSpan {
    let reference = case.reference();
    let geometry = RepeatGeometry {
        left_flank_len: Bp(case.left_flank.len() as u64),
        right_flank_len: Bp(case.right_flank.len() as u64),
        motif: NgMotif::new(&case.motif).expect("a valid motif"),
    };
    let stutter_model = StutterModel::hipstr_shipped();
    let context = RepeatContext {
        geometry: &geometry,
        stutter: &stutter_model,
    };
    let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
    let bases = ReadBases::try_new(&case.read, &case.quality).expect("generated in step");
    aligner.align(bases, &reference, context, scratch)
}

/// Check ng's answer against **itself** — that the case it reports is the one its own offsets
/// imply.
///
/// This exists because the parity oracle cannot see a transposed side: production's
/// `BorderOffEnd` discarded which flank was missing, so swapping `left_anchored` and
/// `right_anchored` leaves every parity assertion satisfied. It is not a tautology, because
/// the side is independently derivable from the offsets — a read that ran off its **3′** end
/// has its span reaching the read's end, and one that began inside the repeat has its span
/// starting at zero.
///
/// **Both implications are conditional on the flank existing**, and the first version of this
/// function was wrong for omitting that — it asserted the unconditional form and its own test
/// caught it. A side reports "not anchored" for *two* different reasons: the read ran past
/// that junction, **or** there is no flank on that side to anchor against (step B2). Only the
/// first constrains the offsets; with a zero-length flank the tract begins at reference
/// column 0 and the read's tract offset can still be non-zero, from a leading insertion.
fn assert_self_consistent(span: &RepeatSpan, case: &Case, context: &str) {
    let read_len = case.read.len() as u64;
    match span {
        RepeatSpan::FromLeft(observed) if !case.right_flank.is_empty() => assert_eq!(
            observed.end, read_len,
            "{context}: the right flank exists but did not anchor, so the read ran off its \
             3′ end and the span must reach the read's end"
        ),
        RepeatSpan::FromRight(observed) if !case.left_flank.is_empty() => assert_eq!(
            observed.start, 0,
            "{context}: the left flank exists but did not anchor, so the read began inside \
             the repeat and the span must start at 0"
        ),
        _ => {}
    }
}

/// How many generated cases each seed contributes. Large enough that the awkward
/// combinations (a long allele *and* a truncation, an out-of-frame indel at a junction) are
/// hit many times over, small enough to stay a unit test.
const CASES_PER_SEED: usize = 3_000;

/// Cases per seed, overridable by `PVC_PARITY_CASES` so a soak run is one command away.
///
/// The default is a **sentinel**, not a campaign: every recurrence mutation tried against
/// this harness died within the first 40 cases, so 3,000 × 4 catches drift promptly and
/// cheaply. But the port was originally justified by a far larger sweep, and the ability to
/// repeat that should not need editing the source — `PVC_PARITY_CASES=50000 cargo test
/// ng_measures_the_same_bytes` reproduces it.
fn cases_per_seed() -> usize {
    std::env::var("PVC_PARITY_CASES")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(CASES_PER_SEED)
}

/// The seeds. Several, because one generator stream can systematically miss a corner.
const SEEDS: [u64; 4] = [
    0x5EED_0001,
    0xC0FF_EE42,
    0x1234_5678_9ABC_DEF0,
    0xDEAD_BEEF_CAFE,
];

/// **The port anchor.** Every read production measures, ng measures identically — byte for
/// byte, not approximately — and every read production rejects lands in a case that says so.
///
/// A failure here means the port has drifted from the oracle it was justified by, and the
/// seed plus index in the message is enough to replay the exact case.
#[test]
fn ng_measures_the_same_bytes_as_production() {
    let mut measured = 0usize;
    let mut off_end = 0usize;
    let mut flankless_divergences = 0usize;
    let mut off_end_shapes = [0usize; 3]; // FromLeft, FromRight, Unanchored
    // One scratch for the whole run, so its grow-but-never-clear reuse is exercised across
    // thousands of differently-sized reads and frames.
    let mut scratch = ViterbiScratch::new();

    for seed in SEEDS {
        let mut rng = SplitMix64(seed);
        for index in 0..cases_per_seed() {
            let case = generate(&mut rng);
            if case.read.is_empty() {
                continue;
            }
            let where_ = format!("seed {seed:#x} case {index}");
            let production = production_answer(&case);
            let ours = ng_answer(&case, &mut scratch);
            assert_self_consistent(&ours, &case, &where_);

            let flankless = case.left_flank.is_empty() || case.right_flank.is_empty();

            match production {
                // THE HARD ASSERTION: the same offsets, not merely the same shape.
                Delimited::Region(ref range) => {
                    // **The one place ng deliberately disagrees.** With a zero-length flank,
                    // production still answers `Region`, because its `left_off`/`right_off`
                    // tests are guarded on the flank existing. ng answers a *bound*, because
                    // a flank that does not exist cannot anchor (step B2): at a contig edge a
                    // read that ended because the tract ended and one that ran out mid-tract
                    // produce an identical readout, so a measurement claim would over-claim
                    // for one of them.
                    //
                    // The divergence is asserted here rather than tripped over. It is in the
                    // *classification* only — where ng still reports a span, the bytes must
                    // match production exactly, which is what the parity claim actually means.
                    if flankless {
                        flankless_divergences += 1;
                        assert!(
                            !ours.is_measurement(),
                            "{where_}: a flankless locus must not yield a measurement ({ours:?})"
                        );
                        if let Some(observed) = ours.observed_span() {
                            assert_eq!(
                                (observed.start, observed.end),
                                (range.start as u64, range.end as u64),
                                "{where_}: flankless bytes diverged from production"
                            );
                        }
                        continue;
                    }

                    measured += 1;
                    let observed = ours.observed_span().unwrap_or_else(|| {
                        panic!(
                            "{where_}: production measured {range:?}, ng reported no span \
                             ({ours:?})"
                        )
                    });
                    assert_eq!(
                        (observed.start, observed.end),
                        (range.start as u64, range.end as u64),
                        "{where_}: measured bytes diverged"
                    );
                    assert!(
                        ours.is_measurement(),
                        "{where_}: production measured it and both flanks exist, so ng must \
                         call it a measurement ({ours:?})"
                    );
                }
                // Production discarded which flank was missing, so this is all the oracle
                // can support: ng must not claim a measurement where production would not.
                Delimited::BorderOffEnd => {
                    off_end += 1;
                    assert!(
                        !ours.is_measurement(),
                        "{where_}: production rejected the read as off-end, ng called it a \
                         measurement ({ours:?})"
                    );
                    match ours {
                        RepeatSpan::FromLeft(_) => off_end_shapes[0] += 1,
                        RepeatSpan::FromRight(_) => off_end_shapes[1] += 1,
                        RepeatSpan::Unanchored => off_end_shapes[2] += 1,
                        RepeatSpan::Between(_) => unreachable!("excluded by the assertion above"),
                    }
                }
            }
        }
    }

    // The fixture has to actually exercise both halves, or the assertions above are vacuous.
    assert!(
        measured > 1_000,
        "too few measured cases to be meaningful: {measured}"
    );
    assert!(
        off_end > 100,
        "too few off-end cases to be meaningful: {off_end}"
    );
    // The flankless divergence must actually be reached, or the assertion encoding it is
    // dead and the harness has silently gone back to being blind to it.
    assert!(
        flankless_divergences > 50,
        "too few flankless cases to exercise the deliberate divergence: {flankless_divergences}"
    );
    // And the reclassification must reach **all three** off-end shapes — the measured
    // distribution produces every one, so anything less means the generator has drifted.
    assert!(
        off_end_shapes.iter().all(|&count| count > 0),
        "the off-end reclassification did not reach every shape: {off_end_shapes:?}"
    );
}

/// The four cases the plan names for the shared fixture, spelled out as fixed inputs rather
/// than left to the generator — so a failure names the biological case directly, and so the
/// fixture is readable by someone checking it by hand.
#[test]
fn the_named_fixture_cases_all_agree_with_production() {
    let motif = b"CAG".to_vec();
    let left_flank = b"ACGTACGT".to_vec();
    let right_flank = b"TTGGTTGGAT".to_vec();
    let tract = b"CAGCAGCAGCAG".to_vec();

    let case_for = |read: &[u8]| Case {
        left_flank: left_flank.clone(),
        tract: tract.clone(),
        right_flank: right_flank.clone(),
        motif: motif.clone(),
        read: read.to_vec(),
        quality: vec![35u8; read.len()],
    };

    let named: [(&str, &[u8]); 6] = [
        ("a clean repeat", b"ACGTACGTCAGCAGCAGCAGTTGGTTGGAT"),
        ("a long allele", b"ACGTACGTCAGCAGCAGCAGCAGCAGCAGTTGGTTGGAT"),
        ("a short allele", b"ACGTACGTCAGCAGTTGGTTGGAT"),
        ("an interrupted repeat", b"ACGTACGTCAGCAGCTGCAGTTGGTTGGAT"),
        ("a read off its 3' end", b"ACGTACGTCAGCAGCAG"),
        ("a read off its 5' end", b"CAGCAGTTGGTTGGAT"),
    ];

    let mut scratch = ViterbiScratch::new();
    for (name, read) in named {
        let case = case_for(read);
        let production = production_answer(&case);
        let ours = ng_answer(&case, &mut scratch);
        assert_self_consistent(&ours, &case, name);

        match production {
            Delimited::Region(ref range) => {
                let observed = ours
                    .observed_span()
                    .unwrap_or_else(|| panic!("{name}: production measured {range:?}, ng did not"));
                assert_eq!(
                    (observed.start, observed.end),
                    (range.start as u64, range.end as u64),
                    "{name}: measured bytes diverged"
                );
            }
            Delimited::BorderOffEnd => assert!(
                !ours.is_measurement(),
                "{name}: production rejected it, ng called it a measurement ({ours:?})"
            ),
        }
    }
}

/// ng's copied slip cutoff must not drift from production's — asserted where the stutter
/// model lives, and re-stated here because this file is where "ng agrees with production" is
/// the subject.
#[test]
fn the_shared_constants_still_agree() {
    assert_eq!(
        stutter::MAX_SLIP as usize,
        crate::ssr::cohort::param_estimation::MAX_SLIP
    );
}
