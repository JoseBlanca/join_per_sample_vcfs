//! Pure HMM-Glocal forward/backward + posterior decoding. Port of
//! htslib's `probaln_glocal`
//! ([htslib/probaln.c:77](../../../htslib/probaln.c#L77)). See Heng Li
//! 2011 (doi:10.1093/bioinformatics/btr076) for the algorithm.
//!
//! The Rust port reads like htslib's PROBALN_ORIG branch
//! (straightforward indexing, no pointer arithmetic) for everything
//! except the forward I-state update, which is written in the OPTIMIZED
//! branch's order because htslib's production build uses OPTIMIZED and
//! the two differ by float-associativity rounding on that one
//! expression. The other forward and the backward expressions are
//! bit-identical between the two branches. Floating-point types match
//! htslib byte-for-byte: `f64` for the DP tables, `f32` for the per-Q
//! lookup and `BaqConfig::gap_open_prob` / `::gap_extend_prob`.

use super::BaqConfig;
use super::errors::ProbalnError;
use super::scratch::{ProbalnScratch, Q2P};

const EI: f64 = 0.25;
// PARITY: htslib uses the decimal literal `.33333333333` (11 digits), not
// 1/3. The two f64 representations differ by ~3.3e-12, which is large
// enough to change the HMM's posterior at mismatch positions enough to
// flip a Phred rounding. Match the literal exactly for parity — do not
// promote to `1.0 / 3.0` or `f64::consts::FRAC_1_3`.
#[allow(clippy::approx_constant)]
pub(super) const EM: f64 = 0.33333333333;

// PARITY: 10/ln(10), 4-digit truncation. Mirrors htslib's literal in
// `probaln.c:285, 388, 421` — keep the truncation; do **not** promote
// to `10.0 / f64::ln(10.0)` or `f64::consts::LOG10_E.recip()`. The
// rounding boundary on `q[i]` depends on this exact value.
pub(super) const PHRED_PER_NAT: f64 = 4.343;

// Rounding bias for `as i32` truncation: `(x + ROUND_HALF_BIAS) as i32`
// rounds half-up for positive `x`. Mirrors htslib.
const ROUND_HALF_BIAS: f64 = 0.499;

// Forward/backward rescale threshold (probaln.c:283). Below this the
// running product `p` is too small to keep in f64; the loop renormalises.
const RESCALE_THRESHOLD: f64 = 1e-100;

// Posterior-Q overflow probe. htslib's `if (qv > 100) qv = 99`.
const PHRED_CAP_OVERFLOW: i32 = 100;
const PHRED_CAP: u8 = 99;

/// Offset of HMM cell `(i, k)` inside its banded row, in 3-cell strides
/// (one stride per HMM cell, holding M/I/D state probabilities). Direct
/// port of htslib's `set_u` macro
/// ([probaln.c:48](../../../htslib/probaln.c#L48)).
///
/// Returns the raw `int` arithmetic htslib does — for callers that pass
/// `(i, k)` outside the band (the `f[l_query+1]`, `b[l_query]`, `b[0]`
/// boundary passes), the result can be negative. Those callers wrap
/// the value in a `< 3 || >= i_dim` bounds check before indexing;
/// negative offsets wrap to a huge `usize` and are caught by the
/// `>= i_dim` half. In-band callers (the main forward / backward loop
/// over `k in beg..=end`) are guaranteed non-negative by construction.
#[inline(always)]
pub(super) fn set_u(bw: i32, i: i32, k: i32) -> usize {
    let x = (i - bw).max(0);
    ((k - x + 1) * 3) as usize
}

/// Posterior-probability local realignment. Direct port of htslib's
/// [`probaln_glocal`](../../../htslib/probaln.c#L77).
///
/// `ref_seq` and `query` are encoded as 0/1/2/3/4 (A/C/G/T/N). `iqual`
/// is the Phred 0-93 per-base quality of `query`, one byte per query
/// base. On success, fills `scratch.state` and `scratch.q` and returns
/// the Phred-scaled alignment likelihood (the value htslib's `int`
/// return carries).
///
/// `scratch.state[i]` encodes the alignment of query base `i`:
/// - `state[i] >> 2` — 0-based reference position the base was aligned
///   to (relative to the start of `ref_seq`).
/// - `state[i] & 3` — 0 if match, 1 if insertion.
///
/// `scratch.q[i]` is the Phred-scaled posterior probability that
/// `state[i]` is wrong, capped at 99 (htslib's same cap at
/// [probaln.c:418](../../../htslib/probaln.c#L418)).
pub(super) fn probaln_glocal(
    scratch: &mut ProbalnScratch,
    ref_seq: &[u8],
    query: &[u8],
    iqual: &[u8],
    cfg: &BaqConfig,
) -> Result<i32, ProbalnError> {
    if iqual.len() != query.len() {
        return Err(ProbalnError::SliceLengthMismatch);
    }
    let l_ref = i32::try_from(ref_seq.len()).map_err(|_| ProbalnError::SequenceTooLong)?;
    let l_query = i32::try_from(query.len()).map_err(|_| ProbalnError::SequenceTooLong)?;
    if l_query >= i32::MAX - 2 {
        return Err(ProbalnError::SequenceTooLong);
    }
    if cfg.band_half_width < 0 {
        return Err(ProbalnError::InvalidBandwidth);
    }
    if l_ref == 0 || l_query == 0 {
        // Resize state/q to the (possibly zero) query length so callers
        // that read `scratch.state` / `scratch.q` after a zero-length
        // call see consistent slice lengths.
        scratch.state.clear();
        scratch.q.clear();
        return Ok(0);
    }

    // Bandwidth selection (probaln.c:93-97).
    let mut bw = l_ref.max(l_query);
    if bw > cfg.band_half_width {
        bw = cfg.band_half_width;
    }
    let len_diff = (l_ref - l_query).abs();
    if bw < len_diff {
        bw = len_diff;
    }
    // PARITY/SAFETY: `bw * 2 + 1` and `bw2 * 3 + 6` in plain i32 can
    // wrap for pathological bandwidths. Surface as AllocationOverflow
    // rather than silently corrupt the matrix shape.
    let bw2 = bw
        .checked_mul(2)
        .and_then(|v| v.checked_add(1))
        .ok_or(ProbalnError::AllocationOverflow)?;
    let i_dim_i32 = if bw2 < l_ref {
        bw2.checked_mul(3)
            .and_then(|v| v.checked_add(6))
            .ok_or(ProbalnError::AllocationOverflow)?
    } else {
        l_ref
            .checked_mul(3)
            .and_then(|v| v.checked_add(6))
            .ok_or(ProbalnError::AllocationOverflow)?
    };
    let i_dim = i_dim_i32 as usize;
    let l_q = l_query as usize;

    // Allocation-size guard (probaln.c:104-107). Done before
    // `resize_for` so an overflow returns `Err` instead of panicking
    // inside Vec.
    let rows = l_q.checked_add(1).ok_or(ProbalnError::AllocationOverflow)?;
    let cells = rows
        .checked_mul(i_dim)
        .ok_or(ProbalnError::AllocationOverflow)?;
    let bytes = cells
        .checked_mul(std::mem::size_of::<f64>())
        .ok_or(ProbalnError::AllocationOverflow)?;
    if bytes > isize::MAX as usize {
        return Err(ProbalnError::AllocationOverflow);
    }

    scratch.resize_for(l_q, i_dim);
    // Borrow-split the disjoint scratch fields so the HMM body can
    // index them by short names without re-walking `scratch.`:
    //   f     — forward DP table, banded, 3 cells per (i, k) for M/I/D.
    //   b     — backward DP table, same shape as `f`.
    //   s     — per-row scale factor s[i], used to renormalise each row.
    //   qual  — per-query-base P_err, populated from `iqual` via Q2P.
    //   state — HMM output: alignment state per query base (out-param).
    //   q     — HMM output: Phred-scaled posterior error (out-param).
    let f = &mut scratch.f;
    let b = &mut scratch.b;
    let s = &mut scratch.s;
    let qual = &mut scratch.qual;
    let state = &mut scratch.state;
    let q = &mut scratch.q;
    // Process-wide shared lookup — see `scratch::Q2P`.
    let q2p: &[f32; 256] = &Q2P;

    // Populate qual from iqual using the shared Phred → P_err lookup.
    for i in 0..l_q {
        qual[i] = q2p[iqual[i] as usize];
    }

    // Transition probabilities (probaln.c:130-141).
    let sm = 1.0 / (2.0 * l_query as f64 + 2.0);
    let si = sm;
    let d = cfg.gap_open_prob as f64;
    let e = cfg.gap_extend_prob as f64;
    // PARITY: trans[from*3 + to] for state 0=M, 1=I, 2=D. Rows in
    // declaration order: M→M, M→I, M→D, I→M, I→I, I→D, D→M, D→I, D→D.
    // Renamed from htslib's `m` to disambiguate from `m_scale` and the
    // M-state index encoding.
    let trans = [
        (1.0 - d - d) * (1.0 - sm), // M→M
        d * (1.0 - sm),             // M→I
        d * (1.0 - sm),             // M→D
        (1.0 - e) * (1.0 - si),     // I→M
        e * (1.0 - si),             // I→I
        0.0,                        // I→D (disallowed)
        1.0 - e,                    // D→M
        0.0,                        // D→I (disallowed)
        e,                          // D→D
    ];
    let bm = (1.0 - d) / l_ref as f64;
    let bi = d / l_ref as f64;

    // === Forward ===
    // f[0] — entry state lives at the (i=0, k=0) padding cell.
    {
        let k = set_u(bw, 0, 0);
        f[k] = 1.0;
        s[0] = 1.0;
    }
    // f[1] — entry into the first query base from the bM/bI priors.
    {
        let beg: i32 = 1;
        let end = l_ref.min(bw + 1);
        let mut sum = 0.0;
        let base = i_dim;
        for k in beg..=end {
            let r = ref_seq[(k - 1) as usize];
            let qb = query[0];
            let qe = qual[0] as f64;
            let e_val = if r > 3 || qb > 3 {
                1.0
            } else if r == qb {
                1.0 - qe
            } else {
                qe * EM
            };
            let u = set_u(bw, 1, k);
            f[base + u] = e_val * bm;
            f[base + u + 1] = EI * bi;
            sum += f[base + u] + f[base + u + 1];
        }
        s[1] = sum;
    }
    // f[2..=l_query] — main forward recursion.
    //
    // PARITY: the xi[1] (I-state) update is written in htslib's
    // OPTIMIZED-branch form, not the PROBALN_ORIG form, because the
    // htslib production build uses OPTIMIZED. The two are mathematically
    // equal but differ by float-associativity rounding (ORIG:
    // `EI * (a + b)`; OPT: `xm3*y3 + xm4*y4` with `xm3 = EI*M*trans[1]`
    // precomputed). Do **not** reassociate these expressions. The xi[0]
    // and xi[2] updates are bit-identical between the two branches.
    for i in 2..=l_query {
        let i_us = i as usize;
        let qli = qual[i_us - 1] as f64;
        let qyi = query[i_us - 1];
        let beg = 1.max(i - bw);
        let end = l_ref.min(i + bw);
        // E[(any_N << 1) | match] — emission probability per base.
        let e_table = [qli * EM, 1.0 - qli, 1.0, 1.0];
        let m_scale = 1.0 / s[i_us - 1];
        // Precomputed xm[3..=4] from htslib's OPT branch.
        let xm3 = EI * m_scale * trans[1];
        let xm4 = EI * m_scale * trans[4];
        let mut sum = 0.0;
        let base_i = i_us * i_dim;
        let base_i1 = (i_us - 1) * i_dim;
        // L8: row-invariant `set_u` terms hoisted out of the inner k-loop.
        // `set_u(bw, i, k) = (k - x_i + 1) * 3`,
        // `set_u(bw, i-1, k) = (k - x_im1 + 1) * 3`, etc.
        let x_i = (i - bw).max(0);
        let x_im1 = ((i - 1) - bw).max(0);
        // UNREACHABLE: i_dim >= bw2*3+6 and every in-band cell offset
        // (k - x_i + 1) * 3 + 2 <= bw2*3 + 2 < i_dim by construction.
        // A single per-row assert hoists all per-cell bounds checks.
        assert!(base_i + i_dim <= f.len());
        assert!(base_i1 + i_dim <= f.len());
        for k in beg..=end {
            let r = ref_seq[(k - 1) as usize];
            let cond = ((r > 3 || qyi > 3) as usize) * 2 + ((r == qyi) as usize);
            let e_val = e_table[cond];
            let u = ((k - x_i + 1) * 3) as usize;
            let v11 = ((k - 1 - x_im1 + 1) * 3) as usize;
            let v10 = ((k - x_im1 + 1) * 3) as usize;
            let v01 = ((k - 1 - x_i + 1) * 3) as usize;
            f[base_i + u] = e_val
                * (trans[0] * m_scale * f[base_i1 + v11]
                    + trans[3] * m_scale * f[base_i1 + v11 + 1]
                    + trans[6] * m_scale * f[base_i1 + v11 + 2]);
            f[base_i + u + 1] = xm3 * f[base_i1 + v10] + xm4 * f[base_i1 + v10 + 1];
            f[base_i + u + 2] = trans[2] * f[base_i + v01] + trans[8] * f[base_i + v01 + 2];
            sum += f[base_i + u] + f[base_i + u + 1] + f[base_i + u + 2];
        }
        s[i_us] = sum;
    }
    // f[l_query+1] — exit. Per probaln.c:263-267 the endpoint must be
    // checked against `i_dim` (a 1.8-1.17 bug used `i_dim - 3`).
    {
        let m_scale = 1.0 / s[l_q];
        let base = l_q * i_dim;
        let mut sum = 0.0;
        for k in 1..=l_ref {
            let u = set_u(bw, l_query, k);
            if u < 3 || u >= i_dim {
                continue;
            }
            sum += m_scale * f[base + u] * sm + m_scale * f[base + u + 1] * si;
        }
        s[l_q + 1] = sum;
    }

    // Likelihood (probaln.c:277-287). Captured for parity with
    // htslib's return value; the BAQ-capping consumer ignores it.
    let pr_score: i32;
    {
        let mut p = 1.0;
        let mut pr1 = 0.0;
        for &scale in &s[..=(l_q + 1)] {
            p *= scale;
            if p < RESCALE_THRESHOLD {
                pr1 += -PHRED_PER_NAT * p.ln();
                p = 1.0;
            }
        }
        pr1 += -PHRED_PER_NAT * (p * l_ref as f64 * l_query as f64).ln();
        // PARITY: `f64 → i32 as` saturates in Rust (post-1.45); NaN
        // becomes 0. htslib has the same behaviour on x86-64 builds.
        debug_assert!(pr1.is_finite(), "pr_score input must be finite");
        pr_score = (pr1 + ROUND_HALF_BIAS) as i32;
    }

    // === Backward ===
    // b[l_query] — exit boundary.
    {
        let base = l_q * i_dim;
        for k in 1..=l_ref {
            let u = set_u(bw, l_query, k);
            if u < 3 || u >= i_dim {
                continue;
            }
            b[base + u] = sm / s[l_q] / s[l_q + 1];
            b[base + u + 1] = si / s[l_q] / s[l_q + 1];
        }
    }
    // b[l_query-1..=1] — backward recursion.
    for i in (1..=(l_query - 1)).rev() {
        let i_us = i as usize;
        let y_d = if i > 1 { 1.0 } else { 0.0 };
        let qli1 = qual[i_us] as f64;
        let qyi1 = query[i_us];
        let beg = 1.max(i - bw);
        let end = l_ref.min(i + bw);
        let e_table = [qli1 * EM, 1.0 - qli1, 1.0, 1.0];
        let base_i = i_us * i_dim;
        let base_i1 = (i_us + 1) * i_dim;
        // L8: row-invariant set_u terms hoisted (uses i+1 for the backward
        // diagonal predecessor rather than i-1 as forward does).
        let x_i = (i - bw).max(0);
        let x_ip1 = ((i + 1) - bw).max(0);
        // UNREACHABLE: same band-width invariant as the forward loop;
        // `base_i1 > base_i` so one assert covers both rows.
        assert!(base_i1 + i_dim <= b.len());
        for k in (beg..=end).rev() {
            let u = ((k - x_i + 1) * 3) as usize;
            let v11 = ((k + 1 - x_ip1 + 1) * 3) as usize;
            let v10 = ((k - x_ip1 + 1) * 3) as usize;
            let v01 = ((k + 1 - x_i + 1) * 3) as usize;
            let e_val = if k >= l_ref {
                0.0
            } else {
                let r = ref_seq[k as usize];
                let cond = ((r > 3 || qyi1 > 3) as usize) * 2 + ((r == qyi1) as usize);
                e_table[cond] * b[base_i1 + v11]
            };
            b[base_i + u] = e_val * trans[0]
                + EI * trans[1] * b[base_i1 + v10 + 1]
                + trans[2] * b[base_i + v01 + 2];
            b[base_i + u + 1] = e_val * trans[3] + EI * trans[4] * b[base_i1 + v10 + 1];
            b[base_i + u + 2] = (e_val * trans[6] + trans[8] * b[base_i + v01 + 2]) * y_d;
        }
        // Rescale this row (probaln.c:327-330).
        let scale = 1.0 / s[i_us];
        let beg_off = ((beg - x_i + 1) * 3) as usize;
        let end_off = ((end - x_i + 1) * 3) as usize + 2;
        for k_off in beg_off..=end_off {
            b[base_i + k_off] *= scale;
        }
    }
    // b[0] — sanity check (probaln.c:366-378). Not used downstream;
    // computing it preserves htslib's full pass for parity-debug logs.
    {
        let beg: i32 = 1;
        let end = l_ref.min(bw + 1);
        let mut sum = 0.0;
        for k in (beg..=end).rev() {
            let r = ref_seq[(k - 1) as usize];
            let qb = query[0];
            let qe = qual[0] as f64;
            let u = set_u(bw, 1, k);
            if u < 3 || u >= i_dim {
                continue;
            }
            let e_val = if r > 3 || qb > 3 {
                1.0
            } else if r == qb {
                1.0 - qe
            } else {
                qe * EM
            };
            sum += e_val * b[i_dim + u] * bm + EI * b[i_dim + u + 1] * bi;
        }
        let k = set_u(bw, 0, 0);
        b[k] = sum / s[0];
    }

    // === MAP / posterior decoding (probaln.c:380-425) ===
    // PARITY: the strict `>` tie-break controls which `k` is selected
    // when two cells share the maximum probability. Do not rewrite as
    // branchless argmax — the parity tests pin the chosen `state[i]`.
    for i in 1..=l_query {
        let i_us = i as usize;
        let beg = 1.max(i - bw);
        let end = l_ref.min(i + bw);
        let m_scale = 1.0 / s[i_us];
        let mut max_val = 0.0;
        let mut max_k: i32 = -1;
        let mut sum = 0.0;
        let base = i_us * i_dim;
        // L8: hoist the row-invariant `(i - bw).max(0)` from `set_u`.
        let x_i = (i - bw).max(0);
        // UNREACHABLE: band-width invariant (see forward loop).
        assert!(base + i_dim <= f.len());
        assert!(base + i_dim <= b.len());
        for k in beg..=end {
            let u = ((k - x_i + 1) * 3) as usize;
            let z_m = m_scale * f[base + u] * b[base + u];
            if z_m > max_val {
                max_val = z_m;
                max_k = (k - 1) << 2;
            }
            sum += z_m;
            let z_i = m_scale * f[base + u + 1] * b[base + u + 1];
            if z_i > max_val {
                max_val = z_i;
                max_k = ((k - 1) << 2) | 1;
            }
            sum += z_i;
        }
        max_val /= sum;
        state[i_us - 1] = max_k;
        // PARITY: htslib's `(int)(-4.343 * log(1.0 - max_val) + 0.499)`.
        // NaN-on-input is impossible here by construction (max_val is
        // a probability in [0, 1]); debug-assert as a regression guard.
        let nat = -PHRED_PER_NAT * (1.0 - max_val).ln();
        debug_assert!(nat.is_finite() || max_val >= 1.0, "qv input non-finite");
        let qv = (nat + ROUND_HALF_BIAS) as i32;
        q[i_us - 1] = if qv > PHRED_CAP_OVERFLOW {
            PHRED_CAP
        } else {
            qv.max(0) as u8
        };
    }

    Ok(pr_score)
}
