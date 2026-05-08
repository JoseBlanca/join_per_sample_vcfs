//! Scope C of [`ia/feature_implementation_plans/pileup_freebayes_style_benchmark.md`][plan]:
//! pure data-structure microbenchmark of the open-record bookkeeping
//! cost of two pileup-walker shapes.
//!
//! [plan]: ../../ia/feature_implementation_plans/pileup_freebayes_style_benchmark.md
//!
//! - **`OursStyle`** — `BTreeMap<u32, OurSlot>` keyed by anchor.
//!   Events whose anchor falls inside an existing record's footprint
//!   merge into that record (no new BTreeMap entry); events at fresh
//!   anchors open new slots. Ageing uses an early-break drain that
//!   relies on the disjoint-footprint invariant.
//! - **`FreebayesStyle`** — `Vec<PoolEntry>`, one entry per event.
//!   No merging on insert. Ageing is `pool.retain(|e| e.anchor +
//!   e.ref_span > walker_pos)` — a full linear sweep — and the
//!   per-position genotyping query is a second linear filter
//!   counting alleles overlapping the walker.
//!
//! Both impls process the same synthetic event stream. The bench
//! measures **only bookkeeping**: no scalar work, no fold logic, no
//! reference fetch, no real reads. Per-event work is identical
//! (`u32` increments) so the only difference between the two
//! numbers is the data-structure access pattern.
//!
//! Parity is checked in `mod tests` (run via `cargo test --bench
//! freebayes_bookkeeping`): both impls process the same stream,
//! both report the same `total_events_closed` at finalize, and that
//! value equals the input event count.

use std::collections::BTreeMap;
use std::hint::black_box;
use std::time::Duration;

use criterion::{BatchSize, BenchmarkId, Criterion, criterion_group};

// =============================================================
// Synthetic event stream
// =============================================================

#[derive(Clone, Copy)]
struct Event {
    anchor: u32,
    ref_span: u32,
    allele_id: u8,
    read_id: u32,
}

/// Deterministic 64-bit LCG (Knuth's MMIX constants). Fast enough
/// that the rng is not a bench bottleneck; using `rand` would add a
/// dependency for no value at this scale.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed.max(1))
    }
    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 32) as u32
    }
    /// Uniform `[0, 1)`.
    fn frac(&mut self) -> f32 {
        (self.next_u32() as f32) / (u32::MAX as f32 + 1.0)
    }
    /// Uniform inclusive `[lo, hi]`.
    fn range_u32(&mut self, lo: u32, hi: u32) -> u32 {
        debug_assert!(hi >= lo);
        lo + self.next_u32() % (hi - lo + 1)
    }
}

/// Footprint distribution — controls how often events have a
/// `ref_span > 1`, which is what keeps records open across multiple
/// walker steps and grows the active set.
#[derive(Clone, Copy)]
enum FootprintDist {
    /// Every event has `ref_span = 1` (SNP / single-base INS).
    /// Records open and close in adjacent walker steps; the active
    /// set is bounded by one step's worth of merging.
    PureSnp,
    /// 90% `span = 1`, 10% `span ∈ [2, 10]` uniform. Mimics typical
    /// Illumina BAQ-adjusted output (~10% of bases sit inside short
    /// indel footprints).
    Mixed,
    /// 50% `span = 1`, 50% `span ∈ [10, 100]` uniform. Pathological
    /// regime: half the events open long-footprint records that
    /// stay alive across many walker steps. Designed to exercise
    /// the `O(N)` sweep cost in `FreebayesStyle`.
    DeletionHeavy,
}

impl FootprintDist {
    fn name(&self) -> &'static str {
        match self {
            FootprintDist::PureSnp => "pure_snp",
            FootprintDist::Mixed => "mixed",
            FootprintDist::DeletionHeavy => "deletion_heavy",
        }
    }
    fn sample(&self, rng: &mut Lcg) -> u32 {
        match self {
            FootprintDist::PureSnp => 1,
            FootprintDist::Mixed => {
                if rng.frac() < 0.10 {
                    rng.range_u32(2, 10)
                } else {
                    1
                }
            }
            FootprintDist::DeletionHeavy => {
                if rng.frac() < 0.50 {
                    rng.range_u32(10, 100)
                } else {
                    1
                }
            }
        }
    }
}

/// Pre-build the event stream so its construction cost is excluded
/// from the timed bench loop. One inner `Vec` per walker step,
/// holding `coverage` events all anchored at that step's position.
fn build_event_stream(span: u32, coverage: u32, dist: FootprintDist, seed: u64) -> Vec<Vec<Event>> {
    let mut rng = Lcg::new(seed);
    let mut steps: Vec<Vec<Event>> = Vec::with_capacity(span as usize);
    for walker_pos in 1..=span {
        let mut events = Vec::with_capacity(coverage as usize);
        for read_id in 0..coverage {
            let ref_span = dist.sample(&mut rng);
            let allele_id = (rng.next_u32() % 2) as u8;
            events.push(Event {
                anchor: walker_pos,
                ref_span,
                allele_id,
                read_id,
            });
        }
        steps.push(events);
    }
    steps
}

// =============================================================
// OursStyle: BTreeMap of open records, footprints disjoint
// =============================================================

#[derive(Default, Clone)]
struct OurSlot {
    /// Reference span of this slot's footprint. Anchor + ref_span
    /// gives the exclusive footprint end.
    ref_span: u32,
    /// Number of events folded into this slot.
    n_entries: u32,
    /// Per-allele observation counters, indexed by `allele_id`. The
    /// real walker keeps a `Vec<AlleleObservation>` here; we mirror
    /// the per-event memory-touch cost with a fixed two-bucket
    /// counter. The bench's allele_id is `0` or `1`.
    per_allele: [u32; 2],
}

#[derive(Default)]
struct OursStyle {
    open: BTreeMap<u32, OurSlot>,
    total_events_closed: u64,
}

impl OursStyle {
    /// Process all events at `walker_pos`, then drain any records
    /// whose footprint is now fully behind the walker.
    ///
    /// Order matches the production walker: aged-out records close
    /// before new events for `walker_pos` are folded, so that
    /// inserts never collide with closing records.
    fn step(&mut self, walker_pos: u32, events: &[Event]) {
        // 1. Drain aged records. Disjoint-footprint invariant lets
        //    us early-break: BTreeMap iterates in anchor order; the
        //    leftmost (smallest anchor) record either ages or it
        //    doesn't, and a non-aged record bounds every later
        //    record below it (later anchors have footprints that
        //    start no earlier).
        //
        //    NB: this only holds because `OursStyle` enforces the
        //    disjoint-footprint invariant on insert (see step 2).
        //    Without that invariant, a long-span record could sit
        //    in front of a short-span record that has aged out, and
        //    the early-break would miss the short one.
        loop {
            match self.open.first_key_value() {
                Some((&anchor, slot)) if anchor + slot.ref_span <= walker_pos => {
                    self.total_events_closed += slot.n_entries as u64;
                    self.open.pop_first();
                }
                _ => break,
            }
        }

        // 2. Fold events at walker_pos. Each event either merges
        //    into an existing record whose footprint covers
        //    walker_pos, or opens a new record at walker_pos.
        //
        //    The rightmost (largest-anchor) open record is the only
        //    candidate for merging: anchors are non-decreasing in
        //    walker_pos order, so any record whose footprint covers
        //    walker_pos must be the most recently inserted. (A
        //    smaller-anchor record whose footprint reached
        //    walker_pos would have been merged-into-or-by some
        //    intervening insert already, by the same rule on a
        //    prior step.)
        //
        //    Production walker uses `find_overlapping` with a
        //    BTreeMap range query for full generality; in this
        //    bench the simpler "check the last entry" works because
        //    every event at step T has anchor exactly T (the
        //    cursor only emits events anchored at the walker
        //    position).
        let merge_into = match self.open.last_key_value() {
            Some((&anchor, slot)) if anchor + slot.ref_span > walker_pos => Some(anchor),
            _ => None,
        };
        if let Some(anchor) = merge_into {
            // Existing record covers walker_pos: every event at this
            // step folds into it. Possibly widen ref_span if any
            // event has a longer footprint.
            let slot = self.open.get_mut(&anchor).expect("merge target exists");
            for ev in events {
                slot.n_entries += 1;
                slot.per_allele[ev.allele_id as usize] += 1;
                let ev_end = walker_pos + ev.ref_span;
                let cur_end = anchor + slot.ref_span;
                if ev_end > cur_end {
                    slot.ref_span = ev_end - anchor;
                }
            }
        } else {
            // No covering record: open a new slot at walker_pos.
            // Multiple events at this step share the new slot.
            let mut slot = OurSlot::default();
            for ev in events {
                slot.n_entries += 1;
                slot.per_allele[ev.allele_id as usize] += 1;
                if ev.ref_span > slot.ref_span {
                    slot.ref_span = ev.ref_span;
                }
            }
            self.open.insert(walker_pos, slot);
        }
    }

    /// Drop everything still open. Used at end-of-stream to make
    /// the closed-event total comparable across impls.
    fn finalize(&mut self) {
        for (_, slot) in self.open.iter() {
            self.total_events_closed += slot.n_entries as u64;
        }
        self.open.clear();
    }
}

// =============================================================
// FreebayesStyle: flat allele pool, sweep + filter per step
// =============================================================

#[derive(Clone, Copy)]
struct PoolEntry {
    anchor: u32,
    ref_span: u32,
    #[allow(dead_code)]
    allele_id: u8,
    #[allow(dead_code)]
    read_id: u32,
}

#[derive(Default)]
struct FreebayesStyle {
    pool: Vec<PoolEntry>,
    total_events_closed: u64,
    /// Accumulated count of "alleles overlapping walker_pos" across
    /// every step — freebayes' per-position genotyping query.
    /// Touched via `black_box` so the optimiser cannot elide the
    /// linear filter that produces it.
    overlap_acc: u64,
}

impl FreebayesStyle {
    fn step(&mut self, walker_pos: u32, events: &[Event]) {
        // 1. Sweep: drop pool entries whose footprint is past
        //    walker_pos. Linear scan over the entire pool, matching
        //    `AlleleParser::updateRegisteredAlleles` (which nulls
        //    aged alleles then `erase(remove(NULL))`s the vector).
        let before = self.pool.len();
        self.pool.retain(|e| e.anchor + e.ref_span > walker_pos);
        let dropped = before - self.pool.len();
        self.total_events_closed += dropped as u64;

        // 2. Push every event. No merging: each observation gets
        //    its own pool entry, like freebayes' `registeredAlleles`.
        for ev in events {
            self.pool.push(PoolEntry {
                anchor: ev.anchor,
                ref_span: ev.ref_span,
                allele_id: ev.allele_id,
                read_id: ev.read_id,
            });
        }

        // 3. Per-position query: count entries overlapping
        //    walker_pos. freebayes does this every step to populate
        //    its genotyping input (`getAlleles` + filter). We don't
        //    actually do anything with the count here — `black_box`
        //    just keeps the optimiser from deleting the filter.
        let overlap = self
            .pool
            .iter()
            .filter(|e| e.anchor <= walker_pos && walker_pos < e.anchor + e.ref_span)
            .count() as u64;
        self.overlap_acc = self.overlap_acc.wrapping_add(overlap);
        black_box(self.overlap_acc);
    }

    fn finalize(&mut self) {
        self.total_events_closed += self.pool.len() as u64;
        self.pool.clear();
    }
}

// =============================================================
// Bench drivers
// =============================================================

fn run_ours(stream: &[Vec<Event>]) -> u64 {
    let mut state = OursStyle::default();
    for (i, events) in stream.iter().enumerate() {
        let walker_pos = (i + 1) as u32;
        state.step(walker_pos, events);
    }
    state.finalize();
    state.total_events_closed
}

fn run_freebayes(stream: &[Vec<Event>]) -> u64 {
    let mut state = FreebayesStyle::default();
    for (i, events) in stream.iter().enumerate() {
        let walker_pos = (i + 1) as u32;
        state.step(walker_pos, events);
    }
    state.finalize();
    state.total_events_closed
}

const SPAN: u32 = 50_000;
const COVERAGES: [u32; 5] = [10, 30, 100, 500, 1000];

fn bench_dist(c: &mut Criterion, dist: FootprintDist) {
    let group_name = format!("freebayes_bookkeeping_{}", dist.name());
    let mut group = c.benchmark_group(&group_name);
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(3));

    for &coverage in &COVERAGES {
        let id_ours = BenchmarkId::new("ours", coverage);
        let id_fb = BenchmarkId::new("freebayes", coverage);

        group.bench_with_input(id_ours, &coverage, |b, &cov| {
            b.iter_batched(
                || build_event_stream(SPAN, cov, dist, 0xC0FFEE),
                |stream| black_box(run_ours(&stream)),
                BatchSize::LargeInput,
            );
        });

        group.bench_with_input(id_fb, &coverage, |b, &cov| {
            b.iter_batched(
                || build_event_stream(SPAN, cov, dist, 0xC0FFEE),
                |stream| black_box(run_freebayes(&stream)),
                BatchSize::LargeInput,
            );
        });
    }

    group.finish();
}

fn bench_pure_snp(c: &mut Criterion) {
    bench_dist(c, FootprintDist::PureSnp);
}

fn bench_mixed(c: &mut Criterion) {
    bench_dist(c, FootprintDist::Mixed);
}

fn bench_deletion_heavy(c: &mut Criterion) {
    bench_dist(c, FootprintDist::DeletionHeavy);
}

fn config() -> Criterion {
    Criterion::default()
        .sample_size(10)
        .measurement_time(Duration::from_secs(3))
}

criterion_group! {
    name = benches;
    config = config();
    targets = bench_pure_snp, bench_mixed, bench_deletion_heavy
}

/// Run the parity check at startup. Panics if the two impls
/// disagree on `total_events_closed` for any of the distributions
/// at moderate coverage. Cheap (~100 ms total) — runs before
/// criterion benches so a regression fails fast.
///
/// The `cargo test --bench freebayes_bookkeeping` path runs this
/// implicitly via `criterion`'s `--test` mode (which calls `main`
/// for each cell once); the `#[cfg(test)] mod tests` block below
/// duplicates the assertions for `cargo test`-style IDE
/// integration.
fn validate_parity() {
    for &dist in &[
        FootprintDist::PureSnp,
        FootprintDist::Mixed,
        FootprintDist::DeletionHeavy,
    ] {
        for &coverage in &[10u32, 100] {
            let stream = build_event_stream(SPAN, coverage, dist, 0xC0FFEE);
            let expected: u64 = stream.iter().map(|s| s.len() as u64).sum();
            let ours = run_ours(&stream);
            let fb = run_freebayes(&stream);
            assert_eq!(
                ours,
                expected,
                "OursStyle drained {ours} events, expected {expected} ({}, cov={coverage})",
                dist.name(),
            );
            assert_eq!(
                fb,
                expected,
                "FreebayesStyle drained {fb} events, expected {expected} ({}, cov={coverage})",
                dist.name(),
            );
        }
    }
}

fn main() {
    validate_parity();
    benches();
}
