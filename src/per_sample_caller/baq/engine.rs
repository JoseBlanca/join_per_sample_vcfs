//! Per-read BAQ driver. Wraps `probaln_glocal` with the CIGAR walks,
//! ref-window extension, and `MappedRead` → `PreparedRead` conversion
//! that match htslib's `sam_prob_realn(b, ref, len, BAQ_APPLY)`. The
//! pileup walker consumes the resulting `PreparedRead.bq_baq` as opaque
//! per-base BQ; this module is the only place capping happens.

use std::sync::Arc;

use crate::per_sample_caller::cram_input::{
    CigarOp, FLAG_FIRST_OF_PAIR, FLAG_PAIRED, FLAG_REVERSE_STRAND, FLAG_UNMAPPED, MappedRead,
};
use crate::per_sample_caller::pileup::{MateRole, PreparedRead, RefSeqFetcher};

use super::BaqConfig;
use super::probaln::probaln_glocal;
use super::scratch::Scratch;

/// Per-read result of [`BaqEngine::process`]: either a fully-built
/// `PreparedRead` with `bq_baq = min(BQ, BAQ)`, or a structured skip
/// reason that the pipeline integration layer counts into
/// `FilterCounts.baq_rejected`.
#[derive(Debug)]
pub enum BaqOutcome {
    Capped(PreparedRead),
    Skipped(BaqSkipReason),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaqSkipReason {
    /// SAM flag `0x4` set. Upstream filtering already drops unmapped
    /// reads, but `BaqEngine` guards defensively in case it is called
    /// from a less-strict caller.
    Unmapped,
    /// `MappedRead.seq` is empty — no bases to score.
    EmptyQuery,
    /// `MappedRead.qual` is empty or its length disagrees with `seq`.
    /// Mirrors htslib's `qual[0] == 0xff` early-return at
    /// [realn.c:134](../../../htslib/realn.c#L134).
    QualAbsent,
    /// CIGAR is entirely insertions / soft-clips / pads — no M/=/X
    /// operation means no reference position to anchor the HMM to.
    /// Mirrors realn.c:192.
    NoMatchInCigar,
    /// CIGAR contains a `BAM_CREF_SKIP` (`N`) operation. RNA-seq
    /// splice gaps are out of scope; reject defensively. Mirrors
    /// realn.c:190.
    ContainsRefSkip,
    /// `probaln_glocal` returned [`BaqOverflow`](super::errors::BaqOverflow)
    /// — practically only triggered by pathological read lengths.
    HmmOverflow,
    /// Reference window extension reached past the chromosome end and
    /// the fetcher refused. Matches htslib's `xe = i; break` clamp at
    /// realn.c:231-234 in spirit — we skip the read rather than try to
    /// run the HMM on a truncated window with no clear interpretation
    /// of the trailing query bases.
    RefWindowPastChromEnd,
}

/// Per-read BAQ driver. Owns reusable scratch buffers so a batched
/// per-thread pass over reads does not allocate.
///
/// `bq_buf` is the in-flight `bq_baq` working buffer that gets handed
/// off to `PreparedRead` via `mem::take` on the success path; the
/// engine starts every new read with a fresh empty `Vec` whose
/// capacity is paid for once (on the first capping that allocates and
/// every time the next read is longer than any prior).
pub struct BaqEngine {
    cfg: BaqConfig,
    scratch: Scratch,
    state: Vec<i32>,
    q: Vec<u8>,
    tref: Vec<u8>,
    tquery: Vec<u8>,
    bq_buf: Vec<u8>,
}

impl BaqEngine {
    pub fn new(cfg: BaqConfig) -> Self {
        Self {
            cfg,
            scratch: Scratch::new(),
            state: Vec::new(),
            q: Vec::new(),
            tref: Vec::new(),
            tquery: Vec::new(),
            bq_buf: Vec::new(),
        }
    }

    /// Run BAQ on one read. See the [`BaqOutcome`] / [`BaqSkipReason`]
    /// docs for the success and skip shapes.
    ///
    /// Takes `MappedRead` by value so the success path can **move**
    /// `cigar`, `seq`, and `qname` straight into the resulting
    /// `PreparedRead` instead of cloning them. The skip paths drop
    /// the `MappedRead` immediately, same as the previous borrowing
    /// API.
    pub fn process(&mut self, read: MappedRead, ref_fetcher: &dyn RefSeqFetcher) -> BaqOutcome {
        if read.flag & FLAG_UNMAPPED != 0 {
            return BaqOutcome::Skipped(BaqSkipReason::Unmapped);
        }
        if read.seq.is_empty() {
            return BaqOutcome::Skipped(BaqSkipReason::EmptyQuery);
        }
        if read.qual.is_empty() || read.qual.len() != read.seq.len() {
            return BaqOutcome::Skipped(BaqSkipReason::QualAbsent);
        }

        let pos_0 = read.pos as i32 - 1;
        let (xb0, xe0, yb, ye, bw) = match locate_alignment_span(&read.cigar, pos_0, self.cfg.bw) {
            Ok(span) => span,
            Err(reason) => return BaqOutcome::Skipped(reason),
        };
        let l_qseq = read.seq.len() as i32;
        let (xb, xe) = extend_ref_window(xb0, xe0, yb, ye, l_qseq, bw);
        let length = (xe - xb) as u32;
        if length == 0 {
            return BaqOutcome::Skipped(BaqSkipReason::RefWindowPastChromEnd);
        }

        let ref_bytes = match ref_fetcher.fetch(read.ref_id as u32, (xb + 1) as u32, length) {
            Ok(bytes) if !bytes.is_empty() => bytes,
            _ => return BaqOutcome::Skipped(BaqSkipReason::RefWindowPastChromEnd),
        };

        self.tref.clear();
        self.tref.extend(ref_bytes.iter().map(|&b| encode_base(b)));
        self.tquery.clear();
        self.tquery.extend(read.seq.iter().map(|&b| encode_base(b)));

        let hmm_cfg = BaqConfig { bw, ..self.cfg };

        // resize without a prior clear(): only the growth above the
        // current length pays a zero-fill (matters for variable-length
        // workloads; uniform L=150 pays the fill exactly once on the
        // first read). The HMM's MAP-decode loop writes every entry
        // of `state` and `q` before they leave the engine, so the
        // residual values from the previous read are harmless.
        self.state.resize(read.seq.len(), 0);
        self.q.resize(read.seq.len(), 0);

        if probaln_glocal(
            &mut self.scratch,
            &self.tref,
            &self.tquery,
            &read.qual,
            &hmm_cfg,
            &mut self.state,
            &mut self.q,
        )
        .is_err()
        {
            return BaqOutcome::Skipped(BaqSkipReason::HmmOverflow);
        }

        apply_baq_cap_into(
            &mut self.bq_buf,
            &read.cigar,
            &read.qual,
            &self.state,
            &self.q,
            pos_0,
            xb,
        );
        let bq_baq = std::mem::take(&mut self.bq_buf);

        BaqOutcome::Capped(mapped_to_prepared(read, bq_baq))
    }
}

// ---------------------------------------------------------------------
// CIGAR / base helpers — shared between `BaqEngine::process` and the
// parity-test driver via `pub(super)` visibility.
// ---------------------------------------------------------------------

/// Lookup table for [`encode_base`]: A/C/G/T → 0/1/2/3 (case-insensitive),
/// anything else → 4. Replacing the previous 5-way match with a
/// straight-line gather lets `Vec::extend(iter.map(encode_base))`
/// vectorize on AVX2 / NEON.
static ENCODE_BASE: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

/// Encode an ASCII base byte (case-insensitive A/C/G/T) into the 0..=3
/// values `probaln_glocal` consumes; anything else (including `N`)
/// becomes 4 (ambiguous). Mirrors htslib's `seq_nt16_int[seq_nt16_table[ch]]`.
#[inline]
pub(super) fn encode_base(b: u8) -> u8 {
    ENCODE_BASE[b as usize]
}

/// First CIGAR walk — find the match span `(xb..xe ref, yb..ye query)`
/// and an indel-widened bandwidth. Mirrors realn.c:178-198.
fn locate_alignment_span(
    cigar: &[CigarOp],
    pos_0: i32,
    cfg_bw: i32,
) -> Result<(i32, i32, i32, i32, i32), BaqSkipReason> {
    let mut x = pos_0;
    let mut y: i32 = 0;
    let mut yb: i32 = -1;
    let mut ye: i32 = -1;
    let mut xb: i32 = -1;
    let mut xe: i32 = -1;
    for op in cigar {
        match *op {
            CigarOp::Match(l) | CigarOp::SeqMatch(l) | CigarOp::SeqMismatch(l) => {
                let l = l as i32;
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
            CigarOp::SoftClip(l) | CigarOp::Insertion(l) => y += l as i32,
            CigarOp::Deletion(l) => x += l as i32,
            CigarOp::Skip(_) => return Err(BaqSkipReason::ContainsRefSkip),
            CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
        }
    }
    if xb < 0 {
        return Err(BaqSkipReason::NoMatchInCigar);
    }
    let mut bw = cfg_bw;
    let diff = ((xe - xb) - (ye - yb)).abs();
    if diff > bw {
        bw = diff + 3;
    }
    Ok((xb, xe, yb, ye, bw))
}

/// Realn.c:200-203 ref-window extension and rebalance. The asymmetric
/// rebalance preserves htslib's comma-operator order-of-eval (the
/// second delta uses the already-updated `xb`).
fn extend_ref_window(
    mut xb: i32,
    mut xe: i32,
    yb: i32,
    ye: i32,
    l_qseq: i32,
    bw: i32,
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
    (xb, xe)
}

/// Second CIGAR walk: cap match-position qualities by the HMM's `q`,
/// writing into `out`. Mirrors the `!extend_baq` branch at
/// realn.c:244-265. `out` is cleared first; the engine reuses one
/// buffer across reads via `mem::take` on success.
fn apply_baq_cap_into(
    out: &mut Vec<u8>,
    cigar: &[CigarOp],
    qual: &[u8],
    state: &[i32],
    q: &[u8],
    pos_0: i32,
    xb: i32,
) {
    out.clear();
    out.extend_from_slice(qual);
    let mut x = pos_0;
    let mut y: usize = 0;
    let l_qseq = qual.len();
    for op in cigar {
        let l = op_len(*op);
        if l == 0 {
            continue;
        }
        match *op {
            CigarOp::Match(_) | CigarOp::SeqMatch(_) | CigarOp::SeqMismatch(_) => {
                let l_clamped = l.min(l_qseq - y);
                for i in y..(y + l_clamped) {
                    let expected_k = (x - xb) + (i as i32 - y as i32);
                    let agree = state[i] & 3 == 0 && state[i] >> 2 == expected_k;
                    out[i] = if agree { qual[i].min(q[i]) } else { 0 };
                }
                x += l as i32;
                y += l_clamped;
            }
            CigarOp::SoftClip(_) | CigarOp::Insertion(_) => {
                y += l.min(l_qseq - y);
            }
            CigarOp::Deletion(_) => {
                x += l as i32;
            }
            _ => {}
        }
    }
}

fn op_len(op: CigarOp) -> usize {
    match op {
        CigarOp::Match(l)
        | CigarOp::Insertion(l)
        | CigarOp::Deletion(l)
        | CigarOp::Skip(l)
        | CigarOp::SoftClip(l)
        | CigarOp::HardClip(l)
        | CigarOp::Padding(l)
        | CigarOp::SeqMatch(l)
        | CigarOp::SeqMismatch(l) => l as usize,
    }
}

fn mapped_to_prepared(read: MappedRead, bq_baq: Vec<u8>) -> PreparedRead {
    let alignment_start = read.pos as u32;
    let alignment_end = compute_alignment_end(alignment_start, &read.cigar);
    let qname = qname_to_arc(&read.qname);
    let mate_role = derive_mate_role(read.flag);
    let is_reverse_strand = read.flag & FLAG_REVERSE_STRAND != 0;
    let mq_log_err = phred_to_ln_perr(read.mapq);
    PreparedRead {
        chrom_id: read.ref_id as u32,
        alignment_start,
        alignment_end,
        cigar: read.cigar,
        seq: read.seq,
        bq_baq,
        mq_log_err,
        is_reverse_strand,
        qname,
        mate_role,
        adaptor_boundary: read.adaptor_boundary,
    }
}

/// SAM-spec qnames are printable ASCII, so the happy path is one
/// `Arc::<str>::from(&str)` allocation. The fallback round-trips
/// through `from_utf8_lossy` only for the rare bad-byte case — without
/// the fast path that's an extra `String` allocation on every read.
fn qname_to_arc(qname: &[u8]) -> Arc<str> {
    match std::str::from_utf8(qname) {
        Ok(s) => Arc::<str>::from(s),
        Err(_) => Arc::<str>::from(String::from_utf8_lossy(qname).as_ref()),
    }
}

fn compute_alignment_end(start_1: u32, cigar: &[CigarOp]) -> u32 {
    let mut ref_span: u32 = 0;
    for op in cigar {
        match *op {
            CigarOp::Match(l)
            | CigarOp::SeqMatch(l)
            | CigarOp::SeqMismatch(l)
            | CigarOp::Deletion(l)
            | CigarOp::Skip(l) => ref_span = ref_span.saturating_add(l),
            _ => {}
        }
    }
    start_1.saturating_add(ref_span).saturating_sub(1)
}

fn derive_mate_role(flag: u16) -> MateRole {
    if flag & FLAG_PAIRED == 0 {
        MateRole::Solo
    } else if flag & FLAG_FIRST_OF_PAIR != 0 {
        MateRole::FirstOfPair
    } else {
        MateRole::SecondOfPair
    }
}

/// `Q -> ln(P_err)` where `P_err = 10^(-Q/10)`. Mirrors the private
/// `phred_to_ln_perr` helper in the pileup walker's `open_record.rs`;
/// kept local here so the BAQ stage does not depend on that module's
/// internal visibility.
fn phred_to_ln_perr(q: u8) -> f64 {
    if q == 0 {
        return 0.0;
    }
    -(q as f64) * std::f64::consts::LN_10 / 10.0
}
