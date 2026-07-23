//! Milestone D2 — the algorithm-3-vs-4 delimiter comparison on synthetic reads.
//!
//! Both delimiters *measure* a read's repeat length; this simulates reads whose true length is
//! known and scores each aligner against it, so the two can be compared where it matters
//! (spec §10.3). Run: `cargo run --example ssr_delimiter_comparison`.
//!
//! # What is measured, and why "calibration" means bias here
//!
//! Algorithm 3 (`SsrFlatGapAligner`) and algorithm 4 (`SsrUnitSlipAligner`) are **best-path**
//! delimiters: each returns one measured length, not a probability. So the marginal sense of
//! "calibration" (does 90%-confident mean right 90% of the time) does not apply. For a
//! point-estimate *ruler*, the analog is its **bias** — does it systematically pull the
//! measurement away from the truth, toward a tidy in-frame length? That is precisely the §4.2
//! risk ("making tidy in-frame lengths cheaper pulls the measurement toward the answer the
//! model already expects, quietly rounding interruptions away"). So three numbers per cell:
//!
//! - **accuracy** — the fraction measured exactly right;
//! - **bias** — the mean *signed* error `measured − true` (calibration: a pull, if non-zero);
//! - **spread** — the mean *absolute* error.
//!
//! # The splits the spec makes mandatory (§10.3)
//!
//! Period 1 **and** period 2-and-above, because they exercise different parts of the model.
//! At period 1, an indel **of the repeat's own base** (genuine stutter) is scored apart from
//! an indel of **a different base** (an interruption the arithmetic in-frame test mis-routes
//! as stutter) — averaged together the two effects cancel (spec §4.2).

use pop_var_caller::ng::alignment::ssr_best_path_flat_gap::{SsrFlatGapAligner, ViterbiScratch};
use pop_var_caller::ng::alignment::ssr_best_path_unit_slip::{SsrUnitSlipAligner, UnitSlipScratch};
use pop_var_caller::ng::alignment::{
    BestPathAligner, PerQualityEmission, ReadBases, RepeatContext, RepeatGeometry, RepeatSpan,
    StutterModel, StutterRates,
};
use pop_var_caller::ng::types::{Bp, Motif};

/// Deterministic PRNG (SplitMix64) so the whole comparison is reproducible from its seed.
struct Rng(u64);
impl Rng {
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
    fn chance(&mut self, p: f64) -> bool {
        let unit = ((self.next_u64() >> 11) as f64) / ((1u64 << 53) as f64);
        unit < p
    }
    /// A base different from `b`.
    fn other_base(&mut self, b: u8) -> u8 {
        let choice = self.below(3);
        *b"ACGT".iter().filter(|&&x| x != b).nth(choice).unwrap()
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4)]
    }
}

/// A canonical motif of a given period (distinct bases so a homopolymer is period 1 only).
fn motif_bytes(period: usize) -> &'static [u8] {
    match period {
        1 => b"A",
        2 => b"CA",
        3 => b"CAG",
        4 => b"CAGT",
        _ => b"CAGTA",
    }
}

/// The scenario a simulated read is drawn from.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Scenario {
    /// A pure tiling of the motif at some unit count — the baseline.
    Clean,
    /// One interior base substituted for a different one — an interruption, no length change.
    Substitution,
    /// A single base inserted or deleted that is **not** a whole unit — out of frame (period ≥ 2).
    OutOfFrameIndel,
    /// Period 1 only: a copy of the homopolymer's own base inserted or deleted — genuine stutter.
    OwnBaseIndel,
    /// Period 1 only: a **foreign** base inserted — arithmetically in-frame, composition wrong.
    DifferentBaseIndel,
}

impl Scenario {
    fn label(self) -> &'static str {
        match self {
            Scenario::Clean => "clean",
            Scenario::Substitution => "substitution",
            Scenario::OutOfFrameIndel => "out-of-frame indel",
            Scenario::OwnBaseIndel => "own-base indel (p1)",
            Scenario::DifferentBaseIndel => "different-base indel (p1)",
        }
    }
}

/// A simulated read plus the ground truth its tract length should measure.
struct Sample {
    read: Vec<u8>,
    quality: Vec<u8>,
    true_tract_bp: i64,
}

/// The fixed flanks — unique sequence that anchors the tract **without a tiling collision at
/// either junction**, which would let the boundary slide and confound the comparison. The
/// slide condition is a *tiling continuation*, not merely a shared base: the right flank must
/// not start with `motif[0]` (the base that would continue the tract rightward), and the left
/// flank must not end with `motif[period − 1]` (the base that would continue it leftward). The
/// motifs tested have `motif[0] ∈ {A, C}` and `motif[period − 1] ∈ {A, G, T}`, so the right
/// flank starts with `G` (neither `A` nor `C`) and the left flank ends with `C` (none of
/// `A`/`G`/`T`) — safe for every period. (An earlier pick had the right flank start with `C`,
/// which *is* period 3's `motif[0]`, and period-3 clean accuracy collapsed to 83% from the
/// boundary sliding a whole unit into the flank — a good reminder that the collision is about
/// the tiling, not a coincidental base.)
const LEFT_FLANK: &[u8] = b"TGCATGAC";
const RIGHT_FLANK: &[u8] = b"GATCGATCGT";

/// Build the reference frame for a locus of `period`, with a `ref_units`-unit tract.
fn reference(period: usize, ref_units: usize) -> (Vec<u8>, RepeatGeometry) {
    let motif = motif_bytes(period);
    let tract: Vec<u8> = motif
        .iter()
        .cycle()
        .take(period * ref_units)
        .copied()
        .collect();
    let mut reference = LEFT_FLANK.to_vec();
    reference.extend_from_slice(&tract);
    reference.extend_from_slice(RIGHT_FLANK);
    let geometry = RepeatGeometry {
        left_flank_len: Bp(LEFT_FLANK.len() as u64),
        right_flank_len: Bp(RIGHT_FLANK.len() as u64),
        motif: Motif::new(motif).expect("valid motif"),
    };
    (reference, geometry)
}

/// Simulate one read of the given scenario. The true tract length is what the read actually
/// carries after the scenario's edit — that is what a delimiter should measure.
fn simulate(rng: &mut Rng, period: usize, ref_units: usize, scenario: Scenario) -> Sample {
    let motif = motif_bytes(period);

    // The read's true unit count: the reference, jittered by a small stutter draw so accuracy
    // is measured across a spread of lengths, not just the reference length.
    let jitter = rng.below(5) as i64 - 2; // −2..=2 units
    let read_units = (ref_units as i64 + jitter).max(1) as usize;
    let mut tract: Vec<u8> = motif
        .iter()
        .cycle()
        .take(period * read_units)
        .copied()
        .collect();

    match scenario {
        Scenario::Clean => {}
        Scenario::Substitution => {
            if !tract.is_empty() {
                let at = rng.below(tract.len());
                tract[at] = rng.other_base(tract[at]);
            }
        }
        Scenario::OutOfFrameIndel => {
            // Insert or delete a single base — not a whole unit (period ≥ 2 here).
            if rng.chance(0.5) && !tract.is_empty() {
                let at = rng.below(tract.len());
                tract.remove(at);
            } else {
                let at = rng.below(tract.len() + 1);
                tract.insert(at, rng.base());
            }
        }
        Scenario::OwnBaseIndel => {
            // Period 1: insert or delete a copy of the homopolymer's own base — genuine stutter.
            let own = motif[0];
            if rng.chance(0.5) && tract.len() > 1 {
                tract.pop();
            } else {
                tract.push(own);
            }
        }
        Scenario::DifferentBaseIndel => {
            // Period 1: insert a foreign base. Arithmetically in-frame (period 1 divides any
            // length), composition wrong — the mis-routing the spec §4.2 works through.
            let at = rng.below(tract.len() + 1);
            tract.insert(at, rng.other_base(motif[0]));
        }
    }

    let true_tract_bp = tract.len() as i64;

    // Assemble the read and add a modest per-base sequencing substitution error.
    let mut read = LEFT_FLANK.to_vec();
    read.extend_from_slice(&tract);
    read.extend_from_slice(RIGHT_FLANK);
    for b in &mut read {
        if rng.chance(0.01) {
            *b = rng.other_base(*b);
        }
    }
    let quality = vec![30u8; read.len()];

    Sample {
        read,
        quality,
        true_tract_bp,
    }
}

/// Running tallies for one (aligner, period, scenario) cell.
#[derive(Default, Clone, Copy)]
struct Tally {
    n: u64,
    exact: u64,
    signed_error: i64,
    abs_error: i64,
    /// Reads where the tract could not be measured (a flank ran off the end): excluded from the
    /// error stats, counted here so a scenario that produces many is not silently averaged.
    unmeasured: u64,
}

impl Tally {
    fn record(&mut self, span: &RepeatSpan, true_bp: i64) {
        self.n += 1;
        match span.measured_length() {
            Some(measured) => {
                let error = measured as i64 - true_bp;
                if error == 0 {
                    self.exact += 1;
                }
                self.signed_error += error;
                self.abs_error += error.abs();
            }
            None => self.unmeasured += 1,
        }
    }
    fn measured(&self) -> u64 {
        self.n - self.unmeasured
    }
    fn accuracy(&self) -> f64 {
        if self.measured() == 0 {
            return f64::NAN;
        }
        self.exact as f64 / self.measured() as f64
    }
    fn bias(&self) -> f64 {
        if self.measured() == 0 {
            return f64::NAN;
        }
        self.signed_error as f64 / self.measured() as f64
    }
    fn spread(&self) -> f64 {
        if self.measured() == 0 {
            return f64::NAN;
        }
        self.abs_error as f64 / self.measured() as f64
    }
}

fn main() {
    // Contraction-biased parameters — HipSTR's fitted values are, and it is the direction
    // asymmetry that algorithm 4 can express and algorithm 3 cannot.
    let model = StutterModel::new(StutterRates {
        in_up: 0.03,
        in_down: 0.07,
        in_geom: 0.9,
        out_up: 0.004,
        out_down: 0.012,
        out_geom: 0.8,
    });
    let emission = PerQualityEmission::new();
    let flat = SsrFlatGapAligner::new(emission);
    let slip = SsrUnitSlipAligner::new(emission);
    let mut flat_scratch = ViterbiScratch::new();
    let mut slip_scratch = UnitSlipScratch::new();

    let trials = 20_000usize;
    let ref_units = 6usize;

    // (period, scenario) -> (flat tally, slip tally)
    let mut cells: Vec<(usize, Scenario, Tally, Tally)> = Vec::new();
    // How often the two aligners disagree on the measured length, per (period, scenario).
    let mut disagreements: Vec<(usize, Scenario, u64, u64)> = Vec::new();

    let periods_and_scenarios: &[(usize, &[Scenario])] = &[
        (
            1,
            &[
                Scenario::Clean,
                Scenario::OwnBaseIndel,
                Scenario::DifferentBaseIndel,
            ],
        ),
        (
            2,
            &[
                Scenario::Clean,
                Scenario::Substitution,
                Scenario::OutOfFrameIndel,
            ],
        ),
        (
            3,
            &[
                Scenario::Clean,
                Scenario::Substitution,
                Scenario::OutOfFrameIndel,
            ],
        ),
        (
            4,
            &[
                Scenario::Clean,
                Scenario::Substitution,
                Scenario::OutOfFrameIndel,
            ],
        ),
    ];

    let mut rng = Rng(0x5EED_D2C0_FFEE);
    for &(period, scenarios) in periods_and_scenarios {
        let (reference, geometry) = reference(period, ref_units);
        for &scenario in scenarios {
            let mut flat_tally = Tally::default();
            let mut slip_tally = Tally::default();
            let mut disagree = 0u64;
            let mut both_measured = 0u64;

            for _ in 0..trials {
                let sample = simulate(&mut rng, period, ref_units, scenario);
                let context = RepeatContext {
                    geometry: &geometry,
                    stutter: &model,
                };

                let bases = ReadBases::try_new(&sample.read, &sample.quality).expect("lengths");
                let flat_span = flat.align(bases, &reference, context, &mut flat_scratch);
                let bases = ReadBases::try_new(&sample.read, &sample.quality).expect("lengths");
                let slip_span = slip.align(bases, &reference, context, &mut slip_scratch);

                flat_tally.record(&flat_span, sample.true_tract_bp);
                slip_tally.record(&slip_span, sample.true_tract_bp);

                if let (Some(a), Some(b)) =
                    (flat_span.measured_length(), slip_span.measured_length())
                {
                    both_measured += 1;
                    if a != b {
                        disagree += 1;
                    }
                }
            }

            cells.push((period, scenario, flat_tally, slip_tally));
            disagreements.push((period, scenario, disagree, both_measured));
        }
    }

    // ---- Report ----
    println!("# Algorithm 3 vs 4 — synthetic delimiter comparison");
    println!(
        "\n{trials} reads per cell, reference tract {ref_units} units, contraction-biased model.\n\
         accuracy = fraction measured exactly right; bias = mean signed error (measured − true);\n\
         spread = mean absolute error. All in base pairs. Higher accuracy and lower |bias|/spread\n\
         are better.\n"
    );
    println!(
        "{:<28} {:>4} | {:>18} | {:>18} | {:>10}",
        "scenario", "per", "algo 3 (flat gap)", "algo 4 (unit slip)", "disagree"
    );
    println!(
        "{:<28} {:>4} | {:>6} {:>5} {:>5} | {:>6} {:>5} {:>5} | {:>10}",
        "", "", "acc", "bias", "sprd", "acc", "bias", "sprd", "rate"
    );
    println!("{}", "-".repeat(96));
    for ((period, scenario, flat_tally, slip_tally), (_, _, disagree, both)) in
        cells.iter().zip(disagreements.iter())
    {
        let disagree_rate = if *both > 0 {
            *disagree as f64 / *both as f64
        } else {
            f64::NAN
        };
        println!(
            "{:<28} {:>4} | {:>5.1}% {:>+5.2} {:>5.2} | {:>5.1}% {:>+5.2} {:>5.2} | {:>9.2}%",
            scenario.label(),
            period,
            100.0 * flat_tally.accuracy(),
            flat_tally.bias(),
            flat_tally.spread(),
            100.0 * slip_tally.accuracy(),
            slip_tally.bias(),
            slip_tally.spread(),
            100.0 * disagree_rate,
        );
    }

    // Unmeasured counts, so a scenario dominated by off-end reads is not read as agreement.
    println!("\nUnmeasured (a flank ran off the end), per cell:");
    for (period, scenario, flat_tally, slip_tally) in &cells {
        if flat_tally.unmeasured > 0 || slip_tally.unmeasured > 0 {
            println!(
                "  p{period} {:<26} algo3 {} / algo4 {} of {}",
                scenario.label(),
                flat_tally.unmeasured,
                slip_tally.unmeasured,
                flat_tally.n
            );
        }
    }
}
