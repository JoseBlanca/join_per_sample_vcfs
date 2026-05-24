//! Reusable buffers for `probaln_glocal`. `BaqEngine` owns one
//! `ProbalnScratch` and threads it through the HMM on every read so a
//! per-read call does not allocate. The Phred→P_err lookup `Q2P` is
//! process-wide shared via `LazyLock` so each rayon worker reads from
//! one copy rather than carrying its own.

use std::sync::LazyLock;

/// Phred → P_err lookup. Initialised once per process from
/// `10^(-q/10)` for `q in 0..=255`. Mirrors htslib's `g_qual2prob`
/// (`htslib/probaln.c:46, 122-127`).
//
// PARITY: computed at first use because `f32::powf` is not yet
// `const fn`. The static `ENCODE_BASE` LUT in `engine.rs` is built
// at compile time because its initialiser is `const`-evaluable.
pub static Q2P: LazyLock<[f32; 256]> = LazyLock::new(|| {
    let mut t = [0f32; 256];
    for (i, slot) in t.iter_mut().enumerate() {
        *slot = 10f32.powf(-(i as f32) / 10.0);
    }
    t
});

/// Per-call working storage for `probaln_glocal`: the forward (`f`) and
/// backward (`b`) banded DP tables, the per-row scale factors (`s`),
/// the per-base error-probability working vector (`qual`), and the HMM
/// outputs (`state`, `q`) that `probaln_glocal` writes for each query
/// base.
///
/// Sized lazily by `resize_for`. The buffers grow to the largest
/// problem seen so far and stay at that capacity, so a long read
/// followed by a short one pays the long read's allocation only once.
pub struct ProbalnScratch {
    pub(super) f: Vec<f64>,
    pub(super) b: Vec<f64>,
    pub(super) s: Vec<f64>,
    pub(super) qual: Vec<f32>,
    /// HMM output: alignment state for each query base.
    /// `state[i] >> 2` is the 0-based ref position the base was
    /// aligned to; `state[i] & 3` is 0 (match) or 1 (insertion).
    pub(super) state: Vec<i32>,
    /// HMM output: Phred-scaled posterior probability that `state[i]`
    /// is wrong, capped at 99.
    pub(super) q: Vec<u8>,
}

impl Default for ProbalnScratch {
    fn default() -> Self {
        Self::new()
    }
}

impl ProbalnScratch {
    pub fn new() -> Self {
        Self {
            f: Vec::new(),
            b: Vec::new(),
            s: Vec::new(),
            qual: Vec::new(),
            state: Vec::new(),
            q: Vec::new(),
        }
    }

    /// HMM output: the alignment state for each query base.
    /// `state[i] >> 2` is the 0-based ref position the base was
    /// aligned to; `state[i] & 3` is 0 (match) or 1 (insertion).
    pub fn state(&self) -> &[i32] {
        &self.state
    }

    /// HMM output: Phred-scaled posterior probability that `state[i]`
    /// is wrong, capped at 99.
    pub fn q(&self) -> &[u8] {
        &self.q
    }

    /// Resize all DP buffers to fit an `(l_query, i_dim)`-sized problem
    /// and zero every cell. Capacity is retained across calls;
    /// re-zeroing is a single `memset` on each buffer.
    pub(super) fn resize_for(&mut self, l_query: usize, i_dim: usize) {
        let cells = (l_query + 1) * i_dim;
        self.f.clear();
        self.f.resize(cells, 0.0);
        self.b.clear();
        self.b.resize(cells, 0.0);
        self.s.clear();
        self.s.resize(l_query + 2, 0.0);
        self.qual.clear();
        self.qual.resize(l_query, 0.0);
        // state / q: resize without a prior clear() — only the growth
        // above the current length pays a zero-fill. The HMM's
        // MAP-decode loop writes every entry of `state` and `q`
        // before they leave the engine, so the residual values from
        // the previous read are harmless.
        self.state.resize(l_query, 0);
        self.q.resize(l_query, 0);
    }
}
