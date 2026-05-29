//! Indel left-alignment (normalization) of a read's CIGAR.
//!
//! A short-read aligner places an indel inside a tandem repeat or
//! homopolymer at an arbitrary offset among equally-scoring positions.
//! Two reads carrying the *same* biological indel can therefore present
//! it at different CIGAR offsets, so recorded verbatim they become
//! distinct `(anchor, REF, ALT)` triples that the cohort merge buckets
//! apart, fragmenting support. Left-aligning every indel to its leftmost
//! equivalent position canonicalises the representation (the
//! `bcftools norm` / `vt normalize` / GATK form GIAB truth uses) so
//! identical indels consolidate onto one allele.
//!
//! This is a port of GATK `AlignmentUtils.leftAlignIndels` +
//! `normalizeAlleles` (`gatk/.../utils/read/AlignmentUtils.java`),
//! cross-checked against freebayes `LeftAlign.cpp`. Both reference
//! callers left-align by rewriting the read's whole CIGAR — traversing
//! right-to-left, trimming for parsimony, then shifting each indel left
//! across the preceding alignment block as long as the moved bases match
//! **both the reference and the read**, bounded so one indel cannot
//! cross the previous (a collision merges them). See
//! `doc/devel/implementation_plans/indel_normalization.md` and the
//! architecture spec §"Indel normalization (left-alignment)".

use super::CigarOp;

/// Result of left-aligning a read's CIGAR.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct LeftAlignResult {
    /// The normalized CIGAR. Synthesized alignment runs are emitted as
    /// [`CigarOp::Match`]; the original op variant is preserved for the
    /// un-shifted remainder of each alignment block.
    pub(crate) cigar: Vec<CigarOp>,
    /// Reference bases removed because a deletion left-aligned all the
    /// way to the read's start (a leading deletion). The caller bumps
    /// the read's `alignment_start` right by this much. Mirrors GATK's
    /// `leadingDeletionBasesRemoved`. Always 0 when the routine is called
    /// with `remove_deletions_at_ends = false`.
    pub(crate) leading_deletion_bases_removed: u32,
}

// --- CigarOp classification helpers -------------------------------------

#[inline]
fn op_len(op: CigarOp) -> u32 {
    match op {
        CigarOp::Match(n)
        | CigarOp::Insertion(n)
        | CigarOp::Deletion(n)
        | CigarOp::Skip(n)
        | CigarOp::SoftClip(n)
        | CigarOp::HardClip(n)
        | CigarOp::Padding(n)
        | CigarOp::SeqMatch(n)
        | CigarOp::SeqMismatch(n) => n,
    }
}

#[inline]
fn with_len(op: CigarOp, n: u32) -> CigarOp {
    match op {
        CigarOp::Match(_) => CigarOp::Match(n),
        CigarOp::Insertion(_) => CigarOp::Insertion(n),
        CigarOp::Deletion(_) => CigarOp::Deletion(n),
        CigarOp::Skip(_) => CigarOp::Skip(n),
        CigarOp::SoftClip(_) => CigarOp::SoftClip(n),
        CigarOp::HardClip(_) => CigarOp::HardClip(n),
        CigarOp::Padding(_) => CigarOp::Padding(n),
        CigarOp::SeqMatch(_) => CigarOp::SeqMatch(n),
        CigarOp::SeqMismatch(_) => CigarOp::SeqMismatch(n),
    }
}

#[inline]
fn is_indel(op: CigarOp) -> bool {
    matches!(op, CigarOp::Insertion(_) | CigarOp::Deletion(_))
}

#[inline]
fn is_alignment(op: CigarOp) -> bool {
    matches!(
        op,
        CigarOp::Match(_) | CigarOp::SeqMatch(_) | CigarOp::SeqMismatch(_)
    )
}

#[inline]
fn consumes_ref(op: CigarOp) -> bool {
    matches!(
        op,
        CigarOp::Match(_)
            | CigarOp::Deletion(_)
            | CigarOp::Skip(_)
            | CigarOp::SeqMatch(_)
            | CigarOp::SeqMismatch(_)
    )
}

#[inline]
fn consumes_read(op: CigarOp) -> bool {
    matches!(
        op,
        CigarOp::Match(_)
            | CigarOp::Insertion(_)
            | CigarOp::SoftClip(_)
            | CigarOp::SeqMatch(_)
            | CigarOp::SeqMismatch(_)
    )
}

#[inline]
fn length_on_ref(op: CigarOp) -> usize {
    if consumes_ref(op) {
        op_len(op) as usize
    } else {
        0
    }
}

#[inline]
fn length_on_read(op: CigarOp) -> usize {
    if consumes_read(op) {
        op_len(op) as usize
    } else {
        0
    }
}

// --- Index range (half-open [start, end)) -------------------------------
//
// Mirrors GATK's `IndexRange`. Signed fields so the trim/shift bookkeeping
// can momentarily express the right-shift trimming case (start_shift < 0)
// without underflow; indices are always in-bounds when used to read bases.

#[derive(Debug, Clone, Copy)]
struct Range {
    start: i64,
    end: i64,
}

impl Range {
    #[inline]
    fn size(self) -> usize {
        // The range is always well-formed (`start <= end`) when `size` is
        // read for op-length emission; the assert catches a future edit
        // that inverts it rather than letting `as usize` wrap to ~2^64.
        debug_assert!(self.end >= self.start, "inverted Range: {self:?}");
        (self.end - self.start) as usize
    }
    #[inline]
    fn shift(&mut self, n: i64) {
        self.start += n;
        self.end += n;
    }
    #[inline]
    fn shift_left(&mut self, n: i64) {
        self.shift(-n);
    }
    #[inline]
    fn shift_start(&mut self, n: i64) {
        self.start += n;
    }
    #[inline]
    fn shift_start_left(&mut self, n: i64) {
        self.start -= n;
    }
    #[inline]
    fn shift_end_left(&mut self, n: i64) {
        self.end -= n;
    }
}

// --- Allele normalization (the shift core) ------------------------------

#[inline]
fn last_base_on_right_is_same(seqs: &[&[u8]], bounds: &[Range]) -> bool {
    let first = seqs[0][(bounds[0].end - 1) as usize];
    seqs.iter()
        .zip(bounds)
        .all(|(s, b)| s[(b.end - 1) as usize] == first)
}

#[inline]
fn first_base_on_left_is_same(seqs: &[&[u8]], bounds: &[Range]) -> bool {
    let first = seqs[0][bounds[0].start as usize];
    seqs.iter()
        .zip(bounds)
        .all(|(s, b)| s[b.start as usize] == first)
}

#[inline]
fn next_base_on_left_is_same(seqs: &[&[u8]], bounds: &[Range]) -> bool {
    let first = seqs[0][(bounds[0].start - 1) as usize];
    seqs.iter()
        .zip(bounds)
        .all(|(s, b)| s[(b.start - 1) as usize] == first)
}

/// Trim shared flanking bases (parsimony), then shift the alleles left as
/// far as the moved bases stay equal across **all** `seqs`, bounded by
/// `max_shift`. Mutates `bounds` in place and returns
/// `(start_shift, end_shift)` — the number of bases the allele start and
/// end moved left (start may be negative when trimming forced a right
/// shift). Port of GATK `normalizeAlleles`.
fn normalize_alleles(
    seqs: &[&[u8]],
    bounds: &mut [Range],
    max_shift: i64,
    trim: bool,
) -> (i64, i64) {
    let mut start_shift: i64 = 0;
    let mut end_shift: i64 = 0;

    let mut min_size = bounds.iter().map(|b| b.size()).min().unwrap_or(0);

    // Consume redundant shared bases at the end of the alleles.
    while trim && min_size > 0 && last_base_on_right_is_same(seqs, bounds) {
        for b in bounds.iter_mut() {
            b.shift_end_left(1);
        }
        min_size -= 1;
        end_shift += 1;
    }

    // Consume redundant shared bases at the start of the alleles.
    while trim && min_size > 0 && first_base_on_left_is_same(seqs, bounds) {
        for b in bounds.iter_mut() {
            b.shift_start(1);
        }
        min_size -= 1;
        start_shift -= 1;
    }

    // Shift left while the next base on the left is equal across all
    // sequences and the last base on the right is equal. For an empty
    // range (e.g. the reference relative to an insertion) the last base
    // on the right is the next base on the left.
    while start_shift < max_shift
        && next_base_on_left_is_same(seqs, bounds)
        && last_base_on_right_is_same(seqs, bounds)
    {
        for b in bounds.iter_mut() {
            b.shift_left(1);
        }
        start_shift += 1;
        end_shift += 1;
    }

    (start_shift, end_shift)
}

// --- CIGAR builder ------------------------------------------------------

/// Assemble a CIGAR from a left-to-right element list, handling the
/// complications a naive rewrite ignores (mirrors GATK's `CigarBuilder`):
///
/// 1. drop zero-length elements,
/// 2. order a deletion before an adjacent insertion (canonical form),
/// 3. merge consecutive identical operators,
/// 4. when `remove_deletions_at_ends`, strip leading and trailing
///    deletions, recording the leading deletion's reference length so the
///    caller can bump the read start. With it `false`, a leading/trailing
///    deletion is left in place (the walker's cursor rejects first/last-op
///    indels anyway, and keeping `alignment_start` fixed preserves the
///    coordinate-sort invariant the read stream relies on).
fn build_cigar(elements: &[CigarOp], remove_deletions_at_ends: bool) -> LeftAlignResult {
    // (1) drop zeros.
    let mut nonzero: Vec<CigarOp> = elements
        .iter()
        .copied()
        .filter(|op| op_len(*op) > 0)
        .collect();

    // (2) within each maximal run of consecutive indels, place the
    // (single) deletion before the (single) insertion.
    {
        let mut i = 0;
        while i < nonzero.len() {
            if is_indel(nonzero[i]) {
                let mut j = i;
                while j < nonzero.len() && is_indel(nonzero[j]) {
                    j += 1;
                }
                let mut del_len = 0u32;
                let mut ins_len = 0u32;
                for op in &nonzero[i..j] {
                    match op {
                        CigarOp::Deletion(n) => del_len += n,
                        CigarOp::Insertion(n) => ins_len += n,
                        _ => unreachable!("run is all indels"),
                    }
                }
                let mut replacement = Vec::with_capacity(2);
                if del_len > 0 {
                    replacement.push(CigarOp::Deletion(del_len));
                }
                if ins_len > 0 {
                    replacement.push(CigarOp::Insertion(ins_len));
                }
                let rep_len = replacement.len();
                nonzero.splice(i..j, replacement);
                i += rep_len;
            } else {
                i += 1;
            }
        }
    }

    // (3) merge consecutive identical operators.
    let mut merged: Vec<CigarOp> = Vec::with_capacity(nonzero.len());
    for op in nonzero {
        if let Some(last) = merged.last_mut()
            && std::mem::discriminant(last) == std::mem::discriminant(&op)
        {
            *last = with_len(*last, op_len(*last) + op_len(op));
        } else {
            merged.push(op);
        }
    }

    // (4) optionally strip leading / trailing deletions. Clips may sit
    // outside the deletion (e.g. `S D M`), so skip past leading/trailing
    // clips when locating the boundary deletion.
    let mut leading_deletion_bases_removed = 0u32;
    if remove_deletions_at_ends {
        let is_clip = |op: CigarOp| matches!(op, CigarOp::SoftClip(_) | CigarOp::HardClip(_));

        // Leading: first non-clip element.
        let lead_idx = merged.iter().position(|op| !is_clip(*op));
        if let Some(idx) = lead_idx
            && let CigarOp::Deletion(n) = merged[idx]
        {
            leading_deletion_bases_removed = n;
            merged.remove(idx);
        }

        // Trailing: last non-clip element.
        let trail_idx = merged.iter().rposition(|op| !is_clip(*op));
        if let Some(idx) = trail_idx
            && matches!(merged[idx], CigarOp::Deletion(_))
        {
            merged.remove(idx);
        }
    }

    LeftAlignResult {
        cigar: merged,
        leading_deletion_bases_removed,
    }
}

// --- Public entry point -------------------------------------------------

/// Left-align every indel in `cigar` to its leftmost equivalent position.
///
/// `ref_bases` are the reference bases spanning the read, with the read's
/// first aligned base at index `read_start`; `read_seq` is the full read
/// sequence (including soft-clipped bases). The returned CIGAR places
/// each indel at the leftmost position reachable without introducing a
/// mismatch on either the reference or the read, with adjacent indels
/// merged. If `cigar` has no indel the input is returned unchanged.
///
/// `remove_deletions_at_ends` controls whether a deletion that left-aligns
/// to a read end is stripped (bumping `alignment_start` via
/// `leading_deletion_bases_removed`) or kept in place. The per-sample read
/// prep passes `false` so `alignment_start` stays fixed and the read
/// stream stays coordinate-sorted; the cursor rejects the resulting
/// first-op deletion anyway.
///
/// Port of GATK `AlignmentUtils.leftAlignIndels`.
pub(crate) fn left_align_cigar(
    cigar: &[CigarOp],
    ref_bases: &[u8],
    read_seq: &[u8],
    read_start: usize,
    remove_deletions_at_ends: bool,
) -> LeftAlignResult {
    if !cigar.iter().copied().any(is_indel) {
        return LeftAlignResult {
            cigar: cigar.to_vec(),
            leading_deletion_bases_removed: 0,
        };
    }

    // Fail safe on malformed or under-provisioned input by leaving the
    // CIGAR untouched (the downstream walker's `read.length()` check then
    // rejects a genuinely malformed read, exactly as without normalization):
    //
    // - the reference window must cover the read's footprint
    //   (`read_start + ref_length`) — should always hold, since the caller
    //   fetches that footprint;
    // - the CIGAR's read-consumption must equal `read_seq.len()`. This is an
    //   *untrusted-input* guard: a CIGAR from a corrupt/adversarial CRAM/BAM
    //   whose read-consuming ops disagree with `seq` would otherwise drive
    //   the right-to-left walk to a non-zero residual and emit a
    //   wrong-length CIGAR in release (where the `read_indel_range.start`
    //   debug-assert below is compiled out).
    let ref_length: usize = cigar.iter().copied().map(length_on_ref).sum();
    let read_length: usize = cigar.iter().copied().map(length_on_read).sum();
    if read_start + ref_length > ref_bases.len()
        || read_length != read_seq.len()
        || read_seq.is_empty()
    {
        return LeftAlignResult {
            cigar: cigar.to_vec(),
            leading_deletion_bases_removed: 0,
        };
    }

    let seqs: [&[u8]; 2] = [ref_bases, read_seq];

    // Traverse right to left. `ref_indel_range` / `read_indel_range`
    // accumulate the bases consumed by a pending indel until we reach the
    // alignment block (or read start) it left-aligns into.
    let mut result_right_to_left: Vec<CigarOp> = Vec::with_capacity(cigar.len() + 4);
    let ref_end = (read_start + ref_length) as i64;
    let read_end = read_seq.len() as i64;
    let mut ref_indel_range = Range {
        start: ref_end,
        end: ref_end,
    };
    let mut read_indel_range = Range {
        start: read_end,
        end: read_end,
    };

    for n in (0..cigar.len()).rev() {
        let element = cigar[n];
        if is_indel(element) {
            // Accumulate; don't shift until we hit an alignment block.
            ref_indel_range.shift_start_left(length_on_ref(element) as i64);
            read_indel_range.shift_start_left(length_on_read(element) as i64);
        } else if ref_indel_range.size() == 0 && read_indel_range.size() == 0 {
            // No pending indel — pass the element through.
            result_right_to_left.push(element);
            ref_indel_range.shift_left(length_on_ref(element) as i64);
            read_indel_range.shift_left(length_on_read(element) as i64);
        } else {
            // A pending indel left-aligns into this element. We may shift
            // into alignment blocks but not into clips.
            let max_shift: i64 = if is_alignment(element) {
                op_len(element) as i64
            } else {
                0
            };
            let mut bounds = [ref_indel_range, read_indel_range];
            let (start_shift, end_shift) = normalize_alleles(&seqs, &mut bounds, max_shift, true);
            ref_indel_range = bounds[0];
            read_indel_range = bounds[1];

            // New match alignments on the right due to left-alignment.
            result_right_to_left.push(CigarOp::Match(end_shift as u32));

            // Emit the indel here unless it can keep shifting into the
            // next block to the left.
            let emit_indel = n == 0 || start_shift < max_shift || !is_alignment(element);
            let new_match_on_left_due_to_trimming = if start_shift < 0 { -start_shift } else { 0 };
            let remaining_bases_on_left = if start_shift < 0 {
                op_len(element) as i64
            } else {
                op_len(element) as i64 - start_shift
            };

            if emit_indel {
                result_right_to_left.push(CigarOp::Deletion(ref_indel_range.size() as u32));
                result_right_to_left.push(CigarOp::Insertion(read_indel_range.size() as u32));
                // Ranges now empty, pointing at the start of the
                // left-aligned indel.
                ref_indel_range.shift_end_left(ref_indel_range.size() as i64);
                read_indel_range.shift_end_left(read_indel_range.size() as i64);

                ref_indel_range.shift_left(
                    new_match_on_left_due_to_trimming
                        + if consumes_ref(element) {
                            remaining_bases_on_left
                        } else {
                            0
                        },
                );
                read_indel_range.shift_left(
                    new_match_on_left_due_to_trimming
                        + if consumes_read(element) {
                            remaining_bases_on_left
                        } else {
                            0
                        },
                );
            }
            result_right_to_left.push(CigarOp::Match(new_match_on_left_due_to_trimming as u32));
            result_right_to_left.push(with_len(element, remaining_bases_on_left as u32));
        }
    }

    // Any indel still pending at the read start (no alignment block to its
    // left) is emitted here; `build_cigar` strips a resulting leading
    // deletion (and records its length) when `remove_deletions_at_ends`.
    result_right_to_left.push(CigarOp::Deletion(ref_indel_range.size() as u32));
    result_right_to_left.push(CigarOp::Insertion(read_indel_range.size() as u32));

    debug_assert_eq!(
        read_indel_range.start, 0,
        "left-aligned cigar does not account for all read bases"
    );

    result_right_to_left.reverse();
    build_cigar(&result_right_to_left, remove_deletions_at_ends)
}

/// Left-align every indel in a read's `cigar` to its leftmost equivalent
/// position, in place. `ref_seq` is the reference slice covering the read's
/// aligned footprint, with `ref_seq[0]` the base at the read's first aligned
/// position (`read_start = 0`); `seq` is the read sequence. The CIGAR is
/// rewritten to the canonical (leftmost) form, leaving the read's start
/// position untouched (`remove_deletions_at_ends = false`, so a deletion
/// that rolls to a read edge stays a first/last-op deletion the cursor
/// rejects). No-op for reads with no indel.
///
/// This is the drop-in replacement for the alignment-input cascade's F3
/// pass: same `(&mut cigar, seq, ref_seq)` shape, but backed by the GATK
/// `leftAlignIndels` port instead of the prior single-forward-pass shifter.
///
/// A debug-only invariant checks that left-alignment — a lossless
/// re-placement of the same indels — leaves the read's reference-mismatch
/// count unchanged, guarding the port against silent corruption.
pub(crate) fn left_align_indels(cigar: &mut Vec<CigarOp>, seq: &[u8], ref_seq: &[u8]) {
    if !cigar.iter().copied().any(is_indel) {
        return;
    }
    let result = left_align_cigar(cigar, ref_seq, seq, 0, false);
    #[cfg(debug_assertions)]
    {
        let before = count_mismatches(cigar, ref_seq, seq, 0);
        let after = count_mismatches(&result.cigar, ref_seq, seq, 0);
        debug_assert_eq!(
            before, after,
            "indel left-alignment changed the read's reference-mismatch count",
        );
    }
    *cigar = result.cigar;
}

/// Count read/reference mismatches across `cigar`'s alignment blocks,
/// with the read's first aligned base at `ref_bases[read_start]`. Used as
/// a debug invariant: left-alignment is a lossless re-placement of the
/// same indels, so it must leave the mismatch count unchanged. Mirrors
/// freebayes' `countMismatches`. Debug-only — the sole caller is the
/// `debug_assert_eq!` in [`left_align_prepared`].
#[cfg(debug_assertions)]
fn count_mismatches(
    cigar: &[CigarOp],
    ref_bases: &[u8],
    read_seq: &[u8],
    read_start: usize,
) -> usize {
    let mut sp = read_start;
    let mut rp = 0usize;
    let mut mismatches = 0usize;
    for &op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                for _ in 0..n {
                    if sp >= ref_bases.len() || rp >= read_seq.len() {
                        return mismatches;
                    }
                    if read_seq[rp] != ref_bases[sp] {
                        mismatches += 1;
                    }
                    sp += 1;
                    rp += 1;
                }
            }
            CigarOp::Deletion(n) | CigarOp::Skip(n) => sp += n as usize,
            CigarOp::Insertion(n) | CigarOp::SoftClip(n) => rp += n as usize,
            CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
        }
    }
    mismatches
}

#[cfg(test)]
mod tests;
