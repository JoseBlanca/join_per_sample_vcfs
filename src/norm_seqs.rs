//! Shared sequence-normalization kernel — the representation-neutral core of
//! indel left-alignment, operating on raw byte sequences and index ranges.
//!
//! Given a set of parallel byte sequences (e.g. a reference allele and an
//! alternate allele, or a reference window and a read) and a half-open index
//! range into each marking the variant span, [`normalize_alleles`] trims shared
//! flanking bases for parsimony and shifts the spans as far left as the moved
//! bases stay equal across **all** sequences. That is the canonical
//! left-aligned, minimally-trimmed form `bcftools norm` / `vt normalize` / GATK
//! produce — the form GIAB truth uses, so identical variants consolidate onto
//! one allele instead of fragmenting across representations.
//!
//! This is a port of GATK `AlignmentUtils.normalizeAlleles`
//! (`gatk/.../utils/read/AlignmentUtils.java`), cross-checked against freebayes
//! `LeftAlign.cpp`. It is **representation-neutral**: it knows nothing about
//! CIGARs or alleles, only sequences and ranges. Two users build on it:
//!
//! - the SNP caller's CIGAR left-alignment ([`crate::pileup::walker::indel_norm`])
//!   wraps it in a right-to-left CIGAR walk;
//! - the SSR caller's off-ladder path (Stage 1 `candidate_generation`) wraps it
//!   in a thin `(seqs, bounds)` adapter to canonicalize an off-ladder tract.
//!
//! Keeping one kernel means a fix to the subtle left-alignment logic reaches
//! both callers; a hand-rolled copy in either would be a divergence-bug farm.
//! See `doc/devel/implementation_plans/indel_normalization.md` and
//! `doc/devel/architecture/ssr_shared_types.md` §4.

// --- Index range (half-open [start, end)) -------------------------------
//
// Mirrors GATK's `IndexRange`. Signed fields so the trim/shift bookkeeping
// can momentarily express the right-shift trimming case (start_shift < 0)
// without underflow; indices are always in-bounds when used to read bases.

/// A half-open `[start, end)` index range into a byte sequence, with signed
/// bounds so trimming can momentarily express a right shift (`start < 0` mid
/// step) without unsigned underflow. Bounds are always in-range when used to
/// read bases. Mirrors GATK's `IndexRange`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct IndexRange {
    pub(crate) start: i64,
    pub(crate) end: i64,
}

impl IndexRange {
    #[inline]
    pub(crate) fn size(self) -> usize {
        // The range is always well-formed (`start <= end`) when `size` is
        // read for op-length emission; the assert catches a future edit
        // that inverts it rather than letting `as usize` wrap to ~2^64.
        debug_assert!(self.end >= self.start, "inverted IndexRange: {self:?}");
        (self.end - self.start) as usize
    }
    #[inline]
    fn shift(&mut self, n: i64) {
        self.start += n;
        self.end += n;
    }
    #[inline]
    pub(crate) fn shift_left(&mut self, n: i64) {
        self.shift(-n);
    }
    #[inline]
    fn shift_start(&mut self, n: i64) {
        self.start += n;
    }
    #[inline]
    pub(crate) fn shift_start_left(&mut self, n: i64) {
        self.start -= n;
    }
    #[inline]
    pub(crate) fn shift_end_left(&mut self, n: i64) {
        self.end -= n;
    }
}

// --- Allele normalization (the shift core) ------------------------------

#[inline]
fn last_base_on_right_is_same(seqs: &[&[u8]], bounds: &[IndexRange]) -> bool {
    let first = seqs[0][(bounds[0].end - 1) as usize];
    seqs.iter()
        .zip(bounds)
        .all(|(s, b)| s[(b.end - 1) as usize] == first)
}

#[inline]
fn first_base_on_left_is_same(seqs: &[&[u8]], bounds: &[IndexRange]) -> bool {
    let first = seqs[0][bounds[0].start as usize];
    seqs.iter()
        .zip(bounds)
        .all(|(s, b)| s[b.start as usize] == first)
}

#[inline]
fn next_base_on_left_is_same(seqs: &[&[u8]], bounds: &[IndexRange]) -> bool {
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
pub(crate) fn normalize_alleles(
    seqs: &[&[u8]],
    bounds: &mut [IndexRange],
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trims_shared_prefix_and_suffix_to_the_minimal_variant() {
        // REF=TCAG, ALT=TGAG differ only at index 1 (C vs G). Parsimony trim
        // should peel the shared "AG" suffix and shared "T" prefix, leaving the
        // single-base variant at [1, 2) in both sequences.
        let r: &[u8] = b"TCAG";
        let a: &[u8] = b"TGAG";
        let seqs: [&[u8]; 2] = [r, a];
        let mut bounds = [
            IndexRange { start: 0, end: 4 },
            IndexRange { start: 0, end: 4 },
        ];
        let (start_shift, end_shift) = normalize_alleles(&seqs, &mut bounds, 0, true);
        // End trimmed twice ("G", "A"); start trimmed once ("T", a right shift).
        assert_eq!((start_shift, end_shift), (-1, 2));
        for b in &bounds {
            assert_eq!((b.start, b.end), (1, 2));
        }
    }

    #[test]
    fn shifts_a_span_left_through_a_homopolymer_run() {
        // Identical sequences with an A-run; a span sitting at the right of the
        // run slides left one base at a time while the moved base stays equal,
        // stopping at max_shift.
        let s: &[u8] = b"GAAAAC";
        let seqs: [&[u8]; 2] = [s, s];
        let mut bounds = [
            IndexRange { start: 4, end: 5 },
            IndexRange { start: 4, end: 5 },
        ];
        // max_shift = 3 stops the slide at [1, 2) before it could read index -1.
        let (start_shift, end_shift) = normalize_alleles(&seqs, &mut bounds, 3, false);
        assert_eq!((start_shift, end_shift), (3, 3));
        for b in &bounds {
            assert_eq!((b.start, b.end), (1, 2));
        }
    }

    #[test]
    fn no_trim_no_shift_leaves_bounds_untouched() {
        // trim=false, max_shift=0: a pure no-op that returns zero shifts.
        let r: &[u8] = b"ACGT";
        let a: &[u8] = b"ATGT";
        let seqs: [&[u8]; 2] = [r, a];
        let mut bounds = [
            IndexRange { start: 1, end: 2 },
            IndexRange { start: 1, end: 2 },
        ];
        let (start_shift, end_shift) = normalize_alleles(&seqs, &mut bounds, 0, false);
        assert_eq!((start_shift, end_shift), (0, 0));
        for b in &bounds {
            assert_eq!((b.start, b.end), (1, 2));
        }
    }
}
