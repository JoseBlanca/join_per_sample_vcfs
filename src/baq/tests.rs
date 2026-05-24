//! Algorithm-only tests for the BAQ core: `probaln_glocal` input
//! validation and output invariants, PARITY-constant pins, `set_u`
//! boundary contract, and the `Q2P` lookup. Pileup-glue tests (the
//! htslib `realn` parity fixtures, engine driver, rayon stream
//! adapter) live in
//! `src/per_sample_pileup/baq_tests.rs`.

use super::BaqConfig;
use super::errors::ProbalnError;
use super::probaln::probaln_glocal;
use super::scratch::{ProbalnScratch, Q2P};

// ---------------------------------------------------------------------
// probaln_glocal — input validation paths (B2)
// ---------------------------------------------------------------------

#[test]
fn probaln_glocal_returns_slice_length_mismatch_on_unequal_iqual() {
    let mut scratch = ProbalnScratch::new();
    let err = probaln_glocal(
        &mut scratch,
        &[0u8; 5],
        &[0u8; 3],
        &[40u8; 2], // iqual shorter than query
        &BaqConfig::default(),
    )
    .unwrap_err();
    assert_eq!(err, ProbalnError::SliceLengthMismatch);
}

#[test]
fn probaln_glocal_returns_ok_zero_on_empty_query() {
    let mut scratch = ProbalnScratch::new();
    let r = probaln_glocal(&mut scratch, &[0u8; 4], &[], &[], &BaqConfig::default()).unwrap();
    assert_eq!(r, 0);
}

#[test]
fn probaln_glocal_returns_ok_zero_on_empty_ref() {
    let mut scratch = ProbalnScratch::new();
    let r = probaln_glocal(
        &mut scratch,
        &[],
        &[0u8; 3],
        &[40u8; 3],
        &BaqConfig::default(),
    )
    .unwrap();
    assert_eq!(r, 0);
}

#[test]
fn probaln_glocal_returns_invalid_bandwidth_on_negative_bw() {
    let mut scratch = ProbalnScratch::new();
    let cfg = BaqConfig {
        band_half_width: -1,
        ..BaqConfig::default()
    };
    let err = probaln_glocal(&mut scratch, &[0u8; 6], &[0u8; 4], &[40u8; 4], &cfg).unwrap_err();
    assert_eq!(err, ProbalnError::InvalidBandwidth);
}

// M16 note: a runtime test for the `bw * 2 + 1` overflow guard would
// require `max(l_ref, l_query) > i32::MAX / 2` — `probaln_glocal`
// clamps the working bw down to `max(l_ref, l_query)` regardless of
// `cfg.band_half_width`. That means triggering the overflow needs
// multi-GB input vectors and cannot be exercised in a unit test on
// realistic hardware. The `checked_mul` chain in `probaln.rs` is
// still load-bearing as defensive arithmetic; the practical
// regression-test is `parity_realn02`, which exercises a non-trivial
// (l_ref, l_query, bw) triple end-to-end.

// ---------------------------------------------------------------------
// probaln_glocal — output invariants (M17)
// ---------------------------------------------------------------------

#[test]
fn probaln_glocal_caps_q_at_99_and_state_within_ref() {
    // All-A reference vs all-T query, qual=1: the HMM saturates the cap.
    let mut scratch = ProbalnScratch::new();
    let ref_seq = vec![0u8; 8];
    let query = vec![3u8; 8];
    let iqual = vec![1u8; 8];
    probaln_glocal(
        &mut scratch,
        &ref_seq,
        &query,
        &iqual,
        &BaqConfig::default(),
    )
    .unwrap();
    let l_ref = ref_seq.len() as i32;
    for (i, (&qv, &st)) in scratch.q.iter().zip(scratch.state.iter()).enumerate() {
        assert!(qv <= 99, "q[{i}] = {qv} exceeds 99 cap");
        let k = st >> 2;
        // state[i] >> 2 < l_ref is the doc invariant. May be -1 if no
        // cell exceeded zero, but should not exceed l_ref - 1.
        assert!(k < l_ref, "state[{i}] >> 2 = {k} exceeds l_ref = {l_ref}");
    }
}

// ---------------------------------------------------------------------
// PARITY constant pins (M8)
// ---------------------------------------------------------------------

#[test]
fn em_constant_matches_htslib_literal() {
    // PARITY: `EM = 0.33333333333` (11 digits) must not drift to
    // `1.0/3.0` (~0.3333333333333333). A unit-test pin catches the
    // drift before the parity fixture does.
    const _: () = {
        assert!(super::probaln::EM > 0.333333333329);
        assert!(super::probaln::EM < 0.333333333331);
    };
}

#[test]
fn phred_per_nat_constant_matches_htslib_literal() {
    // PARITY: 10/ln(10) truncated to 4 digits = 4.343. Catches a
    // "use exact 10.0 / ln(10)" refactor.
    assert_eq!(super::probaln::PHRED_PER_NAT, 4.343);
}

// ---------------------------------------------------------------------
// set_u — boundary contract (M23)
// ---------------------------------------------------------------------

#[test]
fn set_u_returns_three_at_in_band_origin() {
    // set_u(7, 0, 0) = ((0 - max(0-7,0) + 1) * 3) = ((0 - 0 + 1) * 3) = 3
    assert_eq!(super::probaln::set_u(7, 0, 0), 3);
}

#[test]
fn set_u_returns_three_at_in_band_step() {
    // set_u(7, 1, 1) = ((1 - max(1-7,0) + 1) * 3) = ((1 - 0 + 1) * 3) = 6
    assert_eq!(super::probaln::set_u(7, 1, 1), 6);
}

#[test]
fn set_u_wraps_out_of_band_to_huge_usize() {
    // PARITY: for callers that pass (i, k) outside the band,
    // `k - x + 1` may be negative; `as usize` wraps to near usize::MAX.
    // Each caller bounds-checks with `< 3 || >= i_dim` before indexing.
    // This test pins the wrap contract — a refactor that returns
    // Option<usize> would break the boundary loops.
    // set_u(7, 5, -3): x = max(5-7, 0) = 0; (-3 - 0 + 1)*3 = -6 → huge usize.
    let u = super::probaln::set_u(7, 5, -3);
    assert!(u > usize::MAX / 2, "expected wrap-to-huge; got {u}");
}

// ---------------------------------------------------------------------
// Q2P — Phred → P_err lookup invariants
// ---------------------------------------------------------------------

#[test]
fn q2p_returns_one_for_phred_zero_and_decreases_monotonically() {
    let q2p: &[f32; 256] = &Q2P;
    assert!((q2p[0] - 1.0).abs() < 1e-6);
    for w in q2p.windows(2).take(93) {
        assert!(w[0] >= w[1], "Q2P should be monotonic non-increasing");
    }
}
