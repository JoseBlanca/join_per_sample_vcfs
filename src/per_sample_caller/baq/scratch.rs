//! Reusable buffers for `probaln_glocal`. `BaqEngine` owns one
//! `Scratch` and threads it through the HMM on every read so a per-read
//! call does not allocate.

/// Per-call working storage for `probaln_glocal`: the forward (`f`) and
/// backward (`b`) banded DP tables, the per-row scale factors (`s`),
/// the per-base error-probability working vector (`qual`), and the
/// shared Phred → P_err lookup (`q2p`).
///
/// Sized lazily by `resize_for`. The buffers grow to the largest
/// problem seen so far and stay at that capacity, so a long read
/// followed by a short one pays the long read's allocation only once.
pub(super) struct Scratch {
    pub(super) f: Vec<f64>,
    pub(super) b: Vec<f64>,
    pub(super) s: Vec<f64>,
    pub(super) qual: Vec<f32>,
    pub(super) q2p: [f32; 256],
}

impl Scratch {
    pub(super) fn new() -> Self {
        let mut q2p = [0f32; 256];
        for (i, slot) in q2p.iter_mut().enumerate() {
            *slot = 10f32.powf(-(i as f32) / 10.0);
        }
        Self {
            f: Vec::new(),
            b: Vec::new(),
            s: Vec::new(),
            qual: Vec::new(),
            q2p,
        }
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
    }
}
