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
//! lookup and `BaqConfig.d`/`.e`.

use super::BaqConfig;
use super::errors::BaqOverflow;

const EI: f64 = 0.25;
// htslib uses the decimal literal `.33333333333` (11 digits), not 1/3.
// The two f64 representations differ by ~3.3e-12, which is large
// enough to change the HMM's posterior at mismatch positions enough
// to flip a Phred rounding. Match the literal exactly for parity.
const EM: f64 = 0.33333333333;

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
fn set_u(bw: i32, i: i32, k: i32) -> usize {
    let x = (i - bw).max(0);
    ((k - x + 1) * 3) as usize
}

/// Posterior-probability local realignment. Direct port of htslib's
/// [`probaln_glocal`](../../../htslib/probaln.c#L77).
///
/// `ref_seq` and `query` are encoded as 0/1/2/3/4 (A/C/G/T/N). `iqual`
/// is the Phred 0-93 per-base quality of `query`, one byte per query
/// base. On success, fills `state` and `q` and returns the
/// Phred-scaled alignment likelihood (the value htslib's `int` return
/// carries).
///
/// `state[i]` encodes the alignment of query base `i`:
/// - `state[i] >> 2` — 0-based reference position the base was aligned
///   to (relative to the start of `ref_seq`).
/// - `state[i] & 3` — 0 if match, 1 if insertion.
///
/// `q[i]` is the Phred-scaled posterior probability that `state[i]` is
/// wrong, capped at 99 (htslib's same cap at
/// [probaln.c:418](../../../htslib/probaln.c#L418)).
pub(super) fn probaln_glocal(
    ref_seq: &[u8],
    query: &[u8],
    iqual: &[u8],
    cfg: &BaqConfig,
    state: &mut [i32],
    q: &mut [u8],
) -> Result<i32, BaqOverflow> {
    if iqual.len() != query.len() || state.len() != query.len() || q.len() != query.len() {
        return Err(BaqOverflow::InvalidInput);
    }
    let l_ref = i32::try_from(ref_seq.len()).map_err(|_| BaqOverflow::InvalidInput)?;
    let l_query = i32::try_from(query.len()).map_err(|_| BaqOverflow::InvalidInput)?;
    if l_query >= i32::MAX - 2 {
        return Err(BaqOverflow::InvalidInput);
    }
    if l_ref == 0 || l_query == 0 {
        return Ok(0);
    }

    // Bandwidth selection (probaln.c:93-97).
    let mut bw = l_ref.max(l_query);
    if bw > cfg.bw {
        bw = cfg.bw;
    }
    let len_diff = (l_ref - l_query).abs();
    if bw < len_diff {
        bw = len_diff;
    }
    let bw2 = bw * 2 + 1;
    let i_dim = (if bw2 < l_ref {
        bw2 * 3 + 6
    } else {
        l_ref * 3 + 6
    }) as usize;
    let l_q = l_query as usize;

    // Allocation-size guard (probaln.c:104-107).
    let rows = l_q.checked_add(1).ok_or(BaqOverflow::AllocationOverflow)?;
    let cells = rows
        .checked_mul(i_dim)
        .ok_or(BaqOverflow::AllocationOverflow)?;
    let bytes = cells
        .checked_mul(std::mem::size_of::<f64>())
        .ok_or(BaqOverflow::AllocationOverflow)?;
    if bytes > isize::MAX as usize {
        return Err(BaqOverflow::AllocationOverflow);
    }

    let mut f: Vec<f64> = vec![0.0; cells];
    let mut b: Vec<f64> = vec![0.0; cells];
    let mut s: Vec<f64> = vec![0.0; l_q + 2];

    // Per-base error prob lookup (probaln.c:46, 122-127).
    let mut q2p = [0f32; 256];
    for (i, slot) in q2p.iter_mut().enumerate() {
        *slot = 10f32.powf(-(i as f32) / 10.0);
    }
    let mut qual: Vec<f32> = vec![0.0; l_q];
    for i in 0..l_q {
        qual[i] = q2p[iqual[i] as usize];
    }

    // Transition probabilities (probaln.c:130-141).
    let sm = 1.0 / (2.0 * l_query as f64 + 2.0);
    let si = sm;
    let d = cfg.d as f64;
    let e = cfg.e as f64;
    // m[from*3 + to] for state 0=M, 1=I, 2=D.
    let m = [
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
    // The xi[1] (I-state) update is written in htslib's OPTIMIZED-branch
    // form, not the PROBALN_ORIG form, because the htslib production
    // build uses OPTIMIZED. The two are mathematically equal but differ
    // by float-associativity rounding (ORIG: `EI * (a + b)`; OPT:
    // `xm3*y3 + xm4*y4` with `xm3 = EI*M*m[1]` precomputed). The xi[0]
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
        let xm3 = EI * m_scale * m[1];
        let xm4 = EI * m_scale * m[4];
        let mut sum = 0.0;
        let base_i = i_us * i_dim;
        let base_i1 = (i_us - 1) * i_dim;
        for k in beg..=end {
            let r = ref_seq[(k - 1) as usize];
            let cond = ((r > 3 || qyi > 3) as usize) * 2 + ((r == qyi) as usize);
            let e_val = e_table[cond];
            let u = set_u(bw, i, k);
            let v11 = set_u(bw, i - 1, k - 1);
            let v10 = set_u(bw, i - 1, k);
            let v01 = set_u(bw, i, k - 1);
            f[base_i + u] = e_val
                * (m[0] * m_scale * f[base_i1 + v11]
                    + m[3] * m_scale * f[base_i1 + v11 + 1]
                    + m[6] * m_scale * f[base_i1 + v11 + 2]);
            f[base_i + u + 1] = xm3 * f[base_i1 + v10] + xm4 * f[base_i1 + v10 + 1];
            f[base_i + u + 2] = m[2] * f[base_i + v01] + m[8] * f[base_i + v01 + 2];
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
        for i in 0..=(l_q + 1) {
            p *= s[i];
            if p < 1e-100 {
                pr1 += -4.343 * p.ln();
                p = 1.0;
            }
        }
        pr1 += -4.343 * (p * l_ref as f64 * l_query as f64).ln();
        pr_score = (pr1 + 0.499) as i32;
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
        for k in (beg..=end).rev() {
            let u = set_u(bw, i, k);
            let v11 = set_u(bw, i + 1, k + 1);
            let v10 = set_u(bw, i + 1, k);
            let v01 = set_u(bw, i, k + 1);
            let e_val = if k >= l_ref {
                0.0
            } else {
                let r = ref_seq[k as usize];
                let cond = ((r > 3 || qyi1 > 3) as usize) * 2 + ((r == qyi1) as usize);
                e_table[cond] * b[base_i1 + v11]
            };
            b[base_i + u] =
                e_val * m[0] + EI * m[1] * b[base_i1 + v10 + 1] + m[2] * b[base_i + v01 + 2];
            b[base_i + u + 1] = e_val * m[3] + EI * m[4] * b[base_i1 + v10 + 1];
            b[base_i + u + 2] = (e_val * m[6] + m[8] * b[base_i + v01 + 2]) * y_d;
        }
        // Rescale this row (probaln.c:327-330).
        let scale = 1.0 / s[i_us];
        let beg_off = set_u(bw, i, beg);
        let end_off = set_u(bw, i, end) + 2;
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
    for i in 1..=l_query {
        let i_us = i as usize;
        let beg = 1.max(i - bw);
        let end = l_ref.min(i + bw);
        let m_scale = 1.0 / s[i_us];
        let mut max_val = 0.0;
        let mut max_k: i32 = -1;
        let mut sum = 0.0;
        let base = i_us * i_dim;
        for k in beg..=end {
            let u = set_u(bw, i, k);
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
        let qv = (-4.343 * (1.0 - max_val).ln() + 0.499) as i32;
        q[i_us - 1] = if qv > 100 { 99 } else { qv.max(0) as u8 };
    }

    Ok(pr_score)
}
