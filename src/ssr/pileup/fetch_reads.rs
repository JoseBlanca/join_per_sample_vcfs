//! The single-threaded read fetcher (arch §8) — walks the sorted catalog,
//! index-queries each locus's reads, depth-caps them, and (eventually) hands
//! per-locus bundles to the worker pool.
//!
//! **Built so far: the per-locus depth cap — reservoir sampling (§8.3).** This is
//! the net-new, byte-identity-critical piece: the SNP side has only a *column*
//! depth cap (`MPLP_MAX_DEPTH`), not a read-level reservoir. The catalog walk +
//! [`crate::bam::alignment_input::AlignmentMergedReader::query`] driver (mirroring
//! the SNP `run_pileup` shape — load handles once, share the FASTA repository,
//! `clear()` per contig transition), the cheap coordinate-reach admission gate
//! (reusing triage's footprint), the bundle handoff, and the fetcher-thread /
//! bounded-queue / worker-pool / ordered-collector topology all land with the
//! driver, once the container writer exists for output. Single-threaded semantics
//! are built first (trivially deterministic); the pool is a later,
//! determinism-gated step (§8.4).

/// Per-locus cap on **admitted** reads (arch §8.3/§10). A hypervariable,
/// high-depth locus is reservoir-sampled down to this many reads so one locus
/// can't make one worker's task (or one bundle) huge. A **calibration**
/// placeholder (arch §14) — a few hundred spanning reads already genotype a
/// locus; the cap only bites at pathologically deep loci.
pub(crate) const MAX_READS_PER_LOCUS: usize = 1000;

/// A tiny deterministic PRNG (SplitMix64) — seeded per locus so the depth-cap
/// subsample is reproducible and `--threads`-invariant (§8.4), with no external
/// RNG dependency whose stream could shift under us.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e3779b97f4a7615);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
}

/// Deterministic per-locus reservoir seed from `(chrom, start)` (arch §8.3):
/// FNV-1a over the contig name folded with the tract start. Stable for a given
/// locus across runs and thread counts — never derived from wall-clock or
/// thread id. Recorded in the `.ssr.psp` header so the subsample is reproducible.
pub(crate) fn locus_seed(chrom: &str, start: u32) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf29ce484222325;
    const FNV_PRIME: u64 = 0x100000001b3;
    let mut h = FNV_OFFSET;
    for &b in chrom.as_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h ^= start as u64;
    h.wrapping_mul(FNV_PRIME)
}

/// Reservoir sampler (Algorithm R) — an unbiased uniform sample of up to
/// `capacity` items from a stream of unknown length, in one pass with `O(capacity)`
/// memory (arch §8.3). The caller `offer`s each **admitted** read in a fixed
/// total order (`AlignmentMergedReader`'s `(ref_id, pos)` → source-file → record
/// order); with the deterministic per-locus seed, the kept set is identical on
/// every run and at every `--threads`. The caller must not reorder the stream
/// (§8.4).
pub(crate) struct Reservoir<T> {
    capacity: usize,
    held: Vec<T>,
    /// Admitted items offered so far (the `i` of Algorithm R).
    seen: u64,
    rng: SplitMix64,
}

impl<T> Reservoir<T> {
    pub(crate) fn new(capacity: usize, seed: u64) -> Self {
        Self {
            capacity,
            held: Vec::with_capacity(capacity),
            seen: 0,
            rng: SplitMix64::new(seed),
        }
    }

    /// Offer one admitted item. Keeps the first `capacity`; for the `i`-th item
    /// (`i > capacity`) keeps it with probability `capacity / i`, evicting one
    /// held item uniformly at random if kept.
    pub(crate) fn offer(&mut self, item: T) {
        self.seen += 1;
        if self.held.len() < self.capacity {
            self.held.push(item);
        } else {
            // j uniform in [0, seen); replace held[j] when it lands in-window.
            let j = (self.rng.next_u64() % self.seen) as usize;
            if j < self.capacity {
                self.held[j] = item;
            }
        }
    }

    /// The admitted depth (total items offered) — the reservoir sees only
    /// admitted reads, so this is `n_adm`, not the locus's raw depth.
    pub(crate) fn seen(&self) -> u64 {
        self.seen
    }

    /// Consume the reservoir, yielding the sampled items (≤ `capacity`).
    pub(crate) fn into_held(self) -> Vec<T> {
        self.held
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn keeps_everything_when_offered_at_most_capacity() {
        let mut r = Reservoir::new(5, locus_seed("chr1", 100));
        for x in [10u32, 20, 30] {
            r.offer(x);
        }
        assert_eq!(r.seen(), 3);
        assert_eq!(r.into_held(), vec![10, 20, 30]); // first-K kept, in order
    }

    #[test]
    fn caps_at_capacity_and_counts_all_offers() {
        let mut r = Reservoir::new(10, locus_seed("chr1", 100));
        for x in 1..=100u32 {
            r.offer(x);
        }
        assert_eq!(r.seen(), 100);
        assert_eq!(r.into_held().len(), 10);
    }

    #[test]
    fn is_deterministic_for_a_fixed_seed_and_order() {
        let run = || {
            let mut r = Reservoir::new(10, locus_seed("chr7", 4242));
            for x in 1..=100u32 {
                r.offer(x);
            }
            r.into_held()
        };
        assert_eq!(run(), run());
    }

    #[test]
    fn different_loci_sample_differently() {
        let sample = |chrom, start| {
            let mut r = Reservoir::new(10, locus_seed(chrom, start));
            for x in 1..=100u32 {
                r.offer(x);
            }
            r.into_held()
        };
        // Two distinct loci over the same stream pick (overwhelmingly likely)
        // different subsets — the seed actually drives the sampling.
        assert_ne!(sample("chr1", 100), sample("chr1", 101));
        assert_ne!(sample("chr1", 100), sample("chr2", 100));
    }

    #[test]
    fn every_item_can_be_selected_no_structural_exclusion() {
        // Over many seeds, the union of sampled items covers the whole stream —
        // the reservoir reaches every position, not a fixed prefix/suffix.
        let mut covered = HashSet::new();
        for seed in 0..500u64 {
            let mut r = Reservoir::new(10, seed);
            for x in 1..=100u32 {
                r.offer(x);
            }
            covered.extend(r.into_held());
        }
        assert_eq!(covered.len(), 100);
    }

    #[test]
    fn locus_seed_is_deterministic_and_distinguishes_loci() {
        assert_eq!(locus_seed("chr1", 100), locus_seed("chr1", 100));
        assert_ne!(locus_seed("chr1", 100), locus_seed("chr1", 101));
        assert_ne!(locus_seed("chr1", 100), locus_seed("chr2", 100));
    }
}
