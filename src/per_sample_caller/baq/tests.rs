//! Parity test for `probaln_glocal` against htslib's shipped realn
//! fixtures (`tests/data/baq/realn0N.{sam,fa}` and
//! `realn0N_exp.sam`). Decodes the `BQ:Z:` tag htslib emits for
//! `BAQ_APPLY` mode and asserts our post-cap quality matches byte for
//! byte.

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    alignment::{
        record::{cigar::op::Kind, data::field::Tag},
        record_buf::data::field::Value,
    },
};

use super::BaqConfig;
use super::probaln::probaln_glocal;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("baq")
        .join(name)
}

fn read_single_fasta(path: &Path) -> Vec<u8> {
    let mut reader = fasta::io::Reader::new(BufReader::new(File::open(path).unwrap()));
    let mut records = reader.records();
    let rec = records
        .next()
        .expect("fixture FASTA must contain a record")
        .unwrap();
    assert!(
        records.next().is_none(),
        "fixture FASTA must contain a single record"
    );
    rec.sequence().as_ref().to_vec()
}

fn encode_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// First CIGAR walk — locate the match-only span (`xb..xe` ref,
/// `yb..ye` query) and the indel-adjusted bandwidth. Mirrors
/// [htslib/realn.c:178-198](../../../htslib/realn.c#L178-L198).
fn locate_alignment_span(
    cigar: &sam::alignment::record_buf::Cigar,
    pos_0: i32,
    cfg_bw: i32,
) -> Option<(i32, i32, i32, i32, i32)> {
    let mut x = pos_0;
    let mut y: i32 = 0;
    let mut yb: i32 = -1;
    let mut ye: i32 = -1;
    let mut xb: i32 = -1;
    let mut xe: i32 = -1;
    for op in cigar.as_ref() {
        let l = op.len() as i32;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if yb < 0 {
                    yb = y;
                }
                if xb < 0 {
                    xb = x;
                }
                ye = y + l;
                xe = x + l;
                x += l;
                y += l;
            }
            Kind::SoftClip | Kind::Insertion => y += l,
            Kind::Deletion => x += l,
            Kind::Skip => return None,
            _ => {}
        }
    }
    if xb < 0 {
        return None;
    }
    let mut bw = cfg_bw;
    let diff = ((xe - xb) - (ye - yb)).abs();
    if diff > bw {
        bw = diff + 3;
    }
    Some((xb, xe, yb, ye, bw))
}

/// Apply the ref-window extension, clamp, and asymmetric rebalance from
/// [realn.c:200-203](../../../htslib/realn.c#L200-L203). The rebalance
/// is faithful to htslib's comma-operator order-of-eval, where the
/// second delta is computed against the *already-updated* `xb`.
fn extend_ref_window(
    mut xb: i32,
    mut xe: i32,
    yb: i32,
    ye: i32,
    l_qseq: i32,
    bw: i32,
    ref_len: i32,
) -> (i32, i32) {
    xb -= yb + bw / 2;
    if xb < 0 {
        xb = 0;
    }
    xe += l_qseq - ye + bw / 2;
    if xe - xb - l_qseq > bw {
        let dx_xb = (xe - xb - l_qseq - bw) / 2;
        xb += dx_xb;
        let dx_xe = (xe - xb - l_qseq - bw) / 2;
        xe -= dx_xe;
    }
    if xe > ref_len {
        xe = ref_len;
    }
    (xb, xe)
}

/// Walk CIGAR a second time and cap match-position qualities with the
/// HMM's `q`. Mirrors [realn.c:244-265](../../../htslib/realn.c#L244-L265).
fn apply_baq_cap(
    cigar: &sam::alignment::record_buf::Cigar,
    qual: &[u8],
    state: &[i32],
    q: &[u8],
    pos_0: i32,
    xb: i32,
) -> Vec<u8> {
    let mut capped = qual.to_vec();
    let mut x = pos_0;
    let mut y: usize = 0;
    let l_qseq = qual.len();
    for op in cigar.as_ref() {
        let l = op.len();
        if l == 0 {
            continue;
        }
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let l_clamped = l.min(l_qseq - y);
                for i in y..(y + l_clamped) {
                    let expected_k = (x - xb) + (i as i32 - y as i32);
                    let agree = state[i] & 3 == 0 && state[i] >> 2 == expected_k;
                    capped[i] = if agree { qual[i].min(q[i]) } else { 0 };
                }
                x += l as i32;
                y += l_clamped;
            }
            Kind::SoftClip | Kind::Insertion => {
                y += l.min(l_qseq - y);
            }
            Kind::Deletion => {
                x += l as i32;
            }
            _ => {}
        }
    }
    capped
}

fn check_parity(stem: &str) {
    let ref_seq = read_single_fasta(&fixture(&format!("{stem}.fa")));
    let ref_len = ref_seq.len() as i32;

    let mut sam_reader = sam::io::Reader::new(BufReader::new(
        File::open(fixture(&format!("{stem}.sam"))).unwrap(),
    ));
    let sam_header = sam_reader.read_header().unwrap();
    let in_records: Vec<_> = sam_reader
        .record_bufs(&sam_header)
        .map(|r| r.unwrap())
        .collect();

    let mut exp_reader = sam::io::Reader::new(BufReader::new(
        File::open(fixture(&format!("{stem}_exp.sam"))).unwrap(),
    ));
    let exp_header = exp_reader.read_header().unwrap();
    let exp_records: Vec<_> = exp_reader
        .record_bufs(&exp_header)
        .map(|r| r.unwrap())
        .collect();

    assert_eq!(
        in_records.len(),
        exp_records.len(),
        "{stem}: input/expected record count mismatch",
    );

    let cfg = BaqConfig::default();

    for (idx, (in_rec, exp_rec)) in in_records.iter().zip(exp_records.iter()).enumerate() {
        let bq_bytes: Vec<u8> = match exp_rec.data().get(&Tag::BASE_ALIGNMENT_QUALITY_OFFSETS) {
            Some(Value::String(bs)) => {
                let slice: &[u8] = bs.as_ref();
                slice.to_vec()
            }
            Some(other) => panic!("{stem} record {idx}: BQ tag has wrong type: {other:?}"),
            None => continue, // htslib produced no BAQ for this read (e.g. no-match CIGAR)
        };

        let pos_0 = in_rec.alignment_start().unwrap().get() as i32 - 1;
        let qual: Vec<u8> = in_rec.quality_scores().as_ref().to_vec();
        let seq_ascii: Vec<u8> = in_rec.sequence().as_ref().to_vec();
        let cigar = in_rec.cigar();
        let l_qseq = seq_ascii.len() as i32;

        let (xb0, xe0, yb, ye, bw) = match locate_alignment_span(cigar, pos_0, cfg.bw) {
            Some(span) => span,
            None => panic!(
                "{stem} record {idx}: BQ tag present but no realignable span (CREF_SKIP or empty match)",
            ),
        };
        let (xb, xe) = extend_ref_window(xb0, xe0, yb, ye, l_qseq, bw, ref_len);

        let tref: Vec<u8> = ref_seq[xb as usize..xe as usize]
            .iter()
            .map(|&b| encode_base(b))
            .collect();
        let tseq: Vec<u8> = seq_ascii.iter().map(|&b| encode_base(b)).collect();

        let mut bw_cfg = cfg;
        bw_cfg.bw = bw;

        let mut state = vec![0i32; seq_ascii.len()];
        let mut q = vec![0u8; seq_ascii.len()];
        probaln_glocal(&tref, &tseq, &qual, &bw_cfg, &mut state, &mut q)
            .expect("probaln_glocal should not fail on fixture inputs");

        let capped = apply_baq_cap(cigar, &qual, &state, &q, pos_0, xb);

        assert_eq!(
            bq_bytes.len(),
            seq_ascii.len(),
            "{stem} record {idx}: BQ tag length mismatch",
        );
        for i in 0..seq_ascii.len() {
            let delta = bq_bytes[i] as i32 - 64;
            assert!(
                delta >= 0,
                "{stem} record {idx} pos {i}: negative BQ delta ({delta})",
            );
            let expected = (qual[i] as i32 - delta) as u8;
            assert_eq!(
                capped[i], expected,
                "{stem} record {idx} pos {i}: capped qual mismatch (ours={}, htslib={}, raw qual={}, delta={})",
                capped[i], expected, qual[i], delta,
            );
        }
    }
}

#[test]
fn parity_realn01() {
    check_parity("realn01");
}

#[test]
fn parity_realn02() {
    check_parity("realn02");
}

// realn03 is intentionally not used as a parity fixture: its expected
// SAM is generated by htslib's `test_realn -e` (extended BAQ), which
// runs a left/right running-max post-pass on top of the HMM output and
// is out of scope for this plan. The realn03 data files remain in
// `tests/data/baq/` so the corpus is easy to extend if/when extended
// BAQ lands. See `tests/data/baq/README.md`.
