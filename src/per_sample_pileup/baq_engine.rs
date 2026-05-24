//! Per-read BAQ driver. Wraps `probaln_glocal` with the CIGAR walks,
//! ref-window extension, and `MappedRead` → `PreparedRead` conversion
//! that match htslib's `sam_prob_realn(b, ref, len, BAQ_APPLY)`. The
//! pileup walker consumes the resulting `PreparedRead.bq_baq` as opaque
//! per-base BQ; this module is the only place capping happens.

use std::sync::Arc;

use crate::fasta::ManualEvictChromRefFetcher;
use crate::per_sample_pileup::cram_input::{
    CigarOp, FLAG_FIRST_OF_PAIR, FLAG_PAIRED, FLAG_REVERSE_STRAND, FLAG_UNMAPPED, MappedRead,
};
use crate::per_sample_pileup::pileup::{MateRole, PreparedRead};

use crate::baq::{BaqConfig, ProbalnScratch, probaln_glocal};

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
    /// [realn.c:134](../../htslib/realn.c#L134).
    QualAbsent,
    /// CIGAR is entirely insertions / soft-clips / pads — no M/=/X
    /// operation means no reference position to anchor the HMM to.
    /// Mirrors realn.c:192.
    NoMatchInCigar,
    /// CIGAR contains a `BAM_CREF_SKIP` (`N`) operation. RNA-seq
    /// splice gaps are out of scope; reject defensively. Mirrors
    /// realn.c:190.
    ContainsRefSkip,
    /// `probaln_glocal` returned [`ProbalnError`](crate::baq::ProbalnError)
    /// — practically only triggered by pathological read lengths.
    HmmOverflow,
    /// Reference window extension reached past the chromosome end and
    /// the fetcher refused. Matches htslib's `xe = i; break` clamp at
    /// realn.c:231-234 in spirit — we skip the read rather than try to
    /// run the HMM on a truncated window with no clear interpretation
    /// of the trailing query bases.
    RefWindowPastChromEnd,
    /// `MappedRead.pos` exceeds `i32::MAX` (positions on chromosomes
    /// longer than 2 Gbp). The HMM uses `i32` for ref-axis arithmetic;
    /// narrowing would corrupt the alignment.
    PosOutOfRange,
    /// `MappedRead.seq.len()` exceeds `i32::MAX`. Same i32-axis reason.
    ReadTooLong,
    /// `MappedRead.ref_id` exceeds `u32::MAX`. The fetcher API uses
    /// `u32` for chromosome ids.
    ChromIdOutOfRange,
}

/// Per-read BAQ driver. Owns reusable scratch buffers so a batched
/// per-thread pass over reads does not allocate.
///
/// `bq_baq_buf` is the in-flight `bq_baq` working buffer that gets
/// handed off to `PreparedRead` via `mem::take` on the success path;
/// the engine starts every new read with a fresh empty `Vec` whose
/// capacity is paid for once (on the first capping that allocates and
/// every time the next read is longer than any prior).
pub struct BaqEngine {
    cfg: BaqConfig,
    scratch: ProbalnScratch,
    encoded_ref: Vec<u8>,
    encoded_query: Vec<u8>,
    bq_baq_buf: Vec<u8>,
}

impl BaqEngine {
    /// Construct a per-thread BAQ engine. Cheap — allocates only the
    /// empty scratch `Vec`s; DP buffers grow lazily on the first
    /// `process` call. Use one engine per worker thread (`BaqStream`
    /// already does this via `map_init`); the engine is not `Sync`.
    pub fn new(cfg: BaqConfig) -> Self {
        Self {
            cfg,
            scratch: ProbalnScratch::new(),
            encoded_ref: Vec::new(),
            encoded_query: Vec::new(),
            bq_baq_buf: Vec::new(),
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
    pub fn process(
        &mut self,
        read: MappedRead,
        ref_fetcher: &mut ManualEvictChromRefFetcher,
    ) -> BaqOutcome {
        if read.flag & FLAG_UNMAPPED != 0 {
            return BaqOutcome::Skipped(BaqSkipReason::Unmapped);
        }
        if read.seq.is_empty() {
            return BaqOutcome::Skipped(BaqSkipReason::EmptyQuery);
        }
        if read.qual.is_empty() || read.qual.len() != read.seq.len() {
            return BaqOutcome::Skipped(BaqSkipReason::QualAbsent);
        }

        // Checked narrowing — `MappedRead` is decoded from external
        // CRAM/SAM and is not type-bounded to the narrower types the
        // BAQ pipeline uses. Plain `as` casts wrap silently in release.
        let pos_0 = match i32::try_from(read.pos).ok().and_then(|p| p.checked_sub(1)) {
            Some(v) => v,
            None => return BaqOutcome::Skipped(BaqSkipReason::PosOutOfRange),
        };
        let l_qseq = match i32::try_from(read.seq.len()) {
            Ok(v) => v,
            Err(_) => return BaqOutcome::Skipped(BaqSkipReason::ReadTooLong),
        };
        // `chrom_id` is no longer passed to the fetcher (the new
        // API binds it to one contig at construction — BAQ's
        // per-worker fetcher is built for the chunk's chrom), but
        // it's still needed downstream by `mapped_to_prepared` for
        // the resulting `PreparedRead`. Range-checking here catches
        // negative/invalid SAM rows the same as before.
        let chrom_id = match u32::try_from(read.ref_id) {
            Ok(v) => v,
            Err(_) => return BaqOutcome::Skipped(BaqSkipReason::ChromIdOutOfRange),
        };

        // (xb..xe): ref span; (yb..ye): query span. Both 0-based,
        // half-open. `bw` may be wider than `cfg.band_half_width` if a
        // single indel forces it.
        let (xb0, xe0, yb, ye, bw) =
            match locate_alignment_span(&read.cigar, pos_0, self.cfg.band_half_width) {
                Ok(span) => span,
                Err(reason) => return BaqOutcome::Skipped(reason),
            };
        let (xb, xe) = extend_ref_window(xb0, xe0, yb, ye, l_qseq, bw);
        if xe <= xb {
            return BaqOutcome::Skipped(BaqSkipReason::RefWindowPastChromEnd);
        }
        let length = (xe - xb) as u32;

        // Fetch into `self.encoded_ref` directly from the borrowed
        // slice — no allocation, no `Vec<u8>` round-trip. The
        // `encoded_ref` `clear`+`extend` happens here (not after the
        // fetch landing) because the borrow on `ref_fetcher` ends as
        // soon as we finish iterating the slice, and that needs to
        // be before any further fetcher mutation.
        self.encoded_ref.clear();
        match ref_fetcher.fetch((xb + 1) as u32, length) {
            Ok(bytes) if !bytes.is_empty() => {
                self.encoded_ref
                    .extend(bytes.iter().map(|&b| encode_base(b)));
            }
            Ok(_empty) => return BaqOutcome::Skipped(BaqSkipReason::RefWindowPastChromEnd),
            Err(e) => {
                // Surface genuine I/O failures distinct from past-chrom-end.
                // The pileup-walker side has no way to tell from skip
                // telemetry alone whether a transient fetch failure
                // occurred, so log it here and still drop the read.
                eprintln!(
                    "warning: BAQ ref fetch failed (ref_id={}, xb={}, length={}): {}; \
                     skipping read",
                    read.ref_id, xb, length, e,
                );
                return BaqOutcome::Skipped(BaqSkipReason::RefWindowPastChromEnd);
            }
        }
        self.encoded_query.clear();
        self.encoded_query
            .extend(read.seq.iter().map(|&b| encode_base(b)));

        // Spell every BaqConfig field — adding a new field becomes a
        // compile error at this site rather than silently inheriting
        // from `self.cfg` (refactor-safety: parity drift hazard).
        let hmm_cfg = BaqConfig {
            gap_open_prob: self.cfg.gap_open_prob,
            gap_extend_prob: self.cfg.gap_extend_prob,
            band_half_width: bw,
        };

        if probaln_glocal(
            &mut self.scratch,
            &self.encoded_ref,
            &self.encoded_query,
            &read.qual,
            &hmm_cfg,
        )
        .is_err()
        {
            return BaqOutcome::Skipped(BaqSkipReason::HmmOverflow);
        }

        apply_baq_cap_into(
            &mut self.bq_baq_buf,
            &read.cigar,
            &read.qual,
            self.scratch.state(),
            self.scratch.q(),
            pos_0,
            xb,
        );
        let bq_baq = std::mem::take(&mut self.bq_baq_buf);

        BaqOutcome::Capped(mapped_to_prepared(read, chrom_id, bq_baq))
    }
}

// ---------------------------------------------------------------------
// CIGAR / base helpers — `pub(super)` so the test module can call them
// directly rather than maintaining a duplicate copy.
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

/// Encode an ASCII base byte into the 0..=3 values `probaln_glocal`
/// consumes: A/a → 0, C/c → 1, G/g → 2, T/t → 3. Anything else
/// (including `N` and every IUPAC ambiguity code) becomes 4
/// ("ambiguous"). htslib's `seq_nt16_int[seq_nt16_table[ch]]` maps the
/// IUPAC ambiguity bases to specific 0..=3 slots; we collapse them
/// all to 4 because the BAQ parity fixtures contain no ambiguity
/// codes and the BAQ caller's downstream consumer cannot use the
/// finer distinction.
#[inline]
pub(super) fn encode_base(b: u8) -> u8 {
    ENCODE_BASE[b as usize]
}

/// First CIGAR walk — find the match span `(xb..xe ref, yb..ye query)`
/// and an indel-widened bandwidth. Mirrors realn.c:178-198.
pub(super) fn locate_alignment_span(
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
///
/// Note: this function does **not** clamp `xe` to the actual reference
/// length. The caller relies on `ChromRefFetcher::fetch` to reject
/// over-end windows; the test driver's separate copy used to clamp
/// against a known `ref_len` for unit-test convenience.
pub(super) fn extend_ref_window(
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
pub(super) fn apply_baq_cap_into(
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
            // Refactor-safety: spell out the no-op variants so a new
            // `CigarOp` variant is a compile error here, not a silent
            // misaligned-cursor bug. `Skip` is rejected upstream by
            // `locate_alignment_span` but spelled out anyway.
            CigarOp::HardClip(_) | CigarOp::Padding(_) | CigarOp::Skip(_) => {}
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

/// Build a [`PreparedRead`] from a [`MappedRead`] *without* running
/// the BAQ HMM — `bq_baq` is a clone of the raw `qual` from the CRAM
/// record. Used by `pop_var_caller pileup --no-baq` to bypass BAQ
/// while keeping the rest of the pipeline unchanged.
///
/// Reuses the same field-wiring as the BAQ-on path
/// ([`mapped_to_prepared`]) so the two branches stay in lockstep on
/// alignment-end computation, mate-role derivation, and adaptor
/// boundary propagation.
pub fn prepare_passthrough(read: MappedRead, chrom_id: u32) -> PreparedRead {
    let bq_baq = read.qual.clone();
    mapped_to_prepared(read, chrom_id, bq_baq)
}

fn mapped_to_prepared(read: MappedRead, chrom_id: u32, bq_baq: Vec<u8>) -> PreparedRead {
    let alignment_start = read.pos as u32;
    let alignment_end = compute_alignment_end(alignment_start, &read.cigar);
    let qname = qname_to_arc(&read.qname);
    let mate_role = derive_mate_role(read.flag);
    let is_reverse_strand = read.flag & FLAG_REVERSE_STRAND != 0;
    let mq_log_err = phred_to_ln_perr(read.mapq);
    let mapq = read.mapq;
    PreparedRead {
        chrom_id,
        alignment_start,
        alignment_end,
        cigar: read.cigar,
        seq: read.seq,
        bq_baq,
        mq_log_err,
        mapq,
        is_reverse_strand,
        qname,
        mate_role,
        adaptor_boundary: read.adaptor_boundary,
    }
}

/// SAM-spec qnames are printable ASCII, so the happy path is one
/// `Arc::<str>::from(&str)` allocation. The fallback consumes the
/// `Cow<str>` from `from_utf8_lossy` via `Arc<str>: From<Cow<'_, str>>`,
/// which steals the `Owned` arm's allocation rather than copying.
fn qname_to_arc(qname: &[u8]) -> Arc<str> {
    match std::str::from_utf8(qname) {
        Ok(s) => Arc::<str>::from(s),
        Err(_) => Arc::<str>::from(String::from_utf8_lossy(qname)),
    }
}

pub(super) fn compute_alignment_end(start_1: u32, cigar: &[CigarOp]) -> u32 {
    let mut ref_span: u32 = 0;
    for op in cigar {
        match *op {
            CigarOp::Match(l)
            | CigarOp::SeqMatch(l)
            | CigarOp::SeqMismatch(l)
            | CigarOp::Deletion(l)
            | CigarOp::Skip(l) => ref_span = ref_span.saturating_add(l),
            // Refactor-safety: explicit no-op variants so a new
            // `CigarOp` variant cannot silently produce a wrong
            // `alignment_end` for the pileup walker.
            CigarOp::Insertion(_)
            | CigarOp::SoftClip(_)
            | CigarOp::HardClip(_)
            | CigarOp::Padding(_) => {}
        }
    }
    start_1.saturating_add(ref_span).saturating_sub(1)
}

/// Derive [`MateRole`] from SAM flags.
///
/// # Defaults / fallbacks
///
/// - `FLAG_PAIRED` unset → `Solo`.
/// - `FLAG_PAIRED` set & `FLAG_FIRST_OF_PAIR` set → `FirstOfPair`.
/// - `FLAG_PAIRED` set & `FLAG_FIRST_OF_PAIR` unset → `SecondOfPair`.
///   This treats `FLAG_SECOND_OF_PAIR` (0x80) as redundant — any
///   paired read that does not declare itself first-of-pair is taken
///   to be second-of-pair. Malformed records (both bits unset, or
///   both set) fall into this branch silently; upstream CRAM
///   validation is expected to reject them.
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
/// internal visibility. Candidate for promotion to a shared
/// `per_sample_pileup/phred.rs` module.
fn phred_to_ln_perr(q: u8) -> f64 {
    if q == 0 {
        return 0.0;
    }
    -(q as f64) * std::f64::consts::LN_10 / 10.0
}
