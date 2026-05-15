// Scenarios for the pileup walker animation.
//
// Each scenario is a small, hand-crafted setup: a stretch of reference,
// a handful of reads with explicit CIGAR + sequence + base qualities, and
// a free-text description that explains *what concept* the scenario
// demonstrates. The simulator (simulator.js) runs the same algorithm the
// production walker runs and records a step-by-step trace; the viewer
// (viewer.js) renders that trace.
//
// Reads must be sorted by alignment_start (the production walker enforces
// this as a coordinate-order invariant). Sequences are uppercase
// {A,C,G,T,N}. Base qualities are Phred 0-93. CIGAR is given as an array
// of {op, len} pairs over op letters {M, I, D, S, H, N, =, X}; for these
// scenarios we only use M, I, D.
//
// MAPQ is 60 throughout (mq_log_err = ln(10^-6) ≈ -13.815). is_first_mate
// and has_mate are SAM flags. adaptor_boundary is null when the scenario
// has no adaptor read-through to demonstrate.

const SCENARIOS = {
  // ----------------------------------------------------------------- 1
  plain_snp: {
    name: "Scenario 1 — Plain SNP, two reads",
    summary:
      "Two reads overlap on the reference. Both observe the same SNP " +
      "at one position. Demonstrates the walker loop, the cursor " +
      "stepping over M-ops, lazy record opening, the per-step " +
      "closure rule, and how an alt allele's per-allele scalars " +
      "accumulate as additional reads contribute.",
    teaches: [
      "the per-step phase order: admit → process → expire → close → advance",
      "staggered admission: r2 enters one walker-step after r1",
      "the cursor's events_at(walker_pos) emits a Match event per M base",
      "open records are created lazily on the first event at a position",
      "the closure rule: a record at pos with ref_span s closes when walker_pos ≥ pos + s",
      "depth accumulation: two reads agreeing on the alt allele double its num_obs and sum its q",
    ],
    reference: { start: 1, seq: "ACGTACGTAC" },
    reads: [
      {
        qname: "r1",
        chrom_id: 0,
        alignment_start: 2,
        alignment_end: 6,
        cigar: [{ op: "M", len: 5 }],
        // ref at pos 2..6 is C G T A C; SNP at pos 5 (A → G).
        seq:    "CGTGC",
        bq:     [30, 30, 30, 30, 30],
        mq_log_err: -13.815,
        is_reverse_strand: false,
        is_first_mate: true,
        has_mate: false,
        adaptor_boundary: null,
      },
      {
        qname: "r2",
        chrom_id: 0,
        alignment_start: 3,
        alignment_end: 7,
        cigar: [{ op: "M", len: 5 }],
        // ref at pos 3..7 is G T A C G; same SNP at pos 5 (A → G).
        seq:    "GTGCG",
        bq:     [30, 30, 30, 30, 30],
        mq_log_err: -13.815,
        is_reverse_strand: false,
        is_first_mate: true,
        has_mate: false,
        adaptor_boundary: null,
      },
    ],
  },

  // ----------------------------------------------------------------- 2
  compound_del: {
    name: "Scenario 2 — Deletion extends REF span, two reads",
    summary:
      "Two reads cover the same region. r1 carries a 3-base deletion " +
      "(CIGAR 3M3D4M); r2 has no deletion (10M, all-REF). The deletion " +
      "is anchored at position 3 (the last matched position before the " +
      "deletion, per VCF convention) and reaches to position 7 — so " +
      "the record at position 3 stays open for four walker steps. At " +
      "that record both alleles end up populated: REF=AGGG (observed " +
      "by r2) and DEL=A (observed by r1). This is the trickiest moment " +
      "in the walker.",
    teaches: [
      "indel anchoring at the M-position before the indel (VCF convention)",
      "REF-span widening: open_record.widen() extends the REF allele " +
        "and every alt's seq when an event reaches past the current span",
      "why closure is delayed: the long REF span keeps pos + ref_span " +
        "large, so the walker has to advance further before the record " +
        "is safe to close",
      "two reads contribute different alleles to the same record — " +
        "REF and DEL coexist with one observation each",
      "records can close out of coordinate order: short-span records " +
        "at pos 4, 5, 6 close before the long-span record at pos 3",
    ],
    reference: { start: 1, seq: "AAAGGGCCCC" },
    reads: [
      {
        qname: "r1",
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 10, // last reference position the read covers
        cigar: [
          { op: "M", len: 3 }, // ref pos 1-3 (AAA), read pos 0-2 (AAA)
          { op: "D", len: 3 }, // ref pos 4-6 (GGG) deleted
          { op: "M", len: 4 }, // ref pos 7-10 (CCCC), read pos 3-6 (CCCC)
        ],
        seq: "AAACCCC",
        bq: [30, 30, 30, 30, 30, 30, 30],
        mq_log_err: -13.815,
        is_reverse_strand: false,
        is_first_mate: true,
        has_mate: false,
        adaptor_boundary: null,
      },
      {
        qname: "r2",
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 10,
        // ref AAAGGGCCCC across pos 1-10; r2 sees every base, no indel.
        cigar: [{ op: "M", len: 10 }],
        seq: "AAAGGGCCCC",
        bq: [30, 30, 30, 30, 30, 30, 30, 30, 30, 30],
        mq_log_err: -13.815,
        is_reverse_strand: false,
        is_first_mate: true,
        has_mate: false,
        adaptor_boundary: null,
      },
    ],
  },

  // ----------------------------------------------------------------- 3
  mate_overlap: {
    name: "Scenario 3 — Mate-pair overlap on a SNP",
    summary:
      "Two reads from the same fragment (same qname, has_mate set) " +
      "overlap on a SNP. The slot allocator gives them the same " +
      "phase-chain SlotId — the only way the walker can later detect " +
      "the overlap is the shared slot. When both mates contribute at " +
      "the SNP position, mate-overlap resolution applies samtools' BQ " +
      "math (sum on agreement, scale-down on disagreement) so the " +
      "molecule isn't double-counted.",
    teaches: [
      "slot allocator pairs mates by qname (refcount pre-bumped to 2 " +
        "on first mate)",
      "second-mate admission reuses the first mate's slot",
      "resolve_mate_overlap_at_pos detects the overlap by " +
        "grouping contributors by slot_id",
      "samtools BQ math: BQ ← min(BQ₁ + BQ₂, 200) when bases agree",
    ],
    reference: { start: 1, seq: "ACGTACGTAC" },
    reads: [
      {
        qname: "pair1",
        chrom_id: 0,
        alignment_start: 2,
        alignment_end: 6,
        cigar: [{ op: "M", len: 5 }],
        // ref CGTAC pos 2-6; SNP at pos 5 (A → G).
        seq:    "CGTGC",
        bq:     [30, 30, 30, 30, 30],
        mq_log_err: -13.815,
        is_reverse_strand: false,
        is_first_mate: true,
        has_mate: true,
        adaptor_boundary: null,
      },
      {
        qname: "pair1",
        chrom_id: 0,
        alignment_start: 4,
        alignment_end: 8,
        cigar: [{ op: "M", len: 5 }],
        // ref TACGT pos 4-8; SNP at pos 5 (A → G), agreeing with mate.
        seq:    "TGCGT",
        bq:     [25, 25, 25, 25, 25],
        mq_log_err: -13.815,
        is_reverse_strand: true,
        is_first_mate: false,
        has_mate: true,
        adaptor_boundary: null,
      },
    ],
  },
};

// Top-level `const` in a non-module script doesn't attach to `window`;
// expose the table explicitly so viewer.js can read it.
window.SCENARIOS = SCENARIOS;
