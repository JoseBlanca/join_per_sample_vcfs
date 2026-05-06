# Pileup walker — Stage 1 per-sample observation aggregation

**Status:** Draft, 2026-05-06. Implementation-ready specification of
the sequential walker that turns a coordinate-sorted stream of
prepared (filtered + BAQ-adjusted) reads into the per-position
records consumed by the Stage 2 encoder. Slots underneath
[per_sample_caller.md](per_sample_caller.md), which already commits
to the high-level Stage 1 shape, the input/output contract, the
five-scalar layout, the phase-chain semantics, the mate-overlap
rule, and the per-read filter cascade. This document fills in the
runtime mechanics those commitments leave open: the in-memory
data structures, the CIGAR-decomposition, the open-record
formation and merging algorithm, the position emission and
closure protocol, the phase-chain slot allocator, and the
mate-overlap implementation.

Open design questions are tagged **[DECISION]** (a choice we have
to make before implementation) or **[QUESTION]** (a sub-question
or detail to confirm).

## Purpose and scope

The pileup walker is the **single-threaded core** of Stage 1
([per_sample_caller.md §"Parallelism"](per_sample_caller.md)).
Its job is to walk the stream of prepared reads, accumulate
per-allele evidence at every reference position covered by at
least one read, and emit a sequence of per-position records to
the Stage 2 encoder.

Inside this scope:
- Active-set bookkeeping (which reads cover the current walker
  position).
- CIGAR-driven event decomposition (turning a read into a list
  of position-anchored allele observations).
- Open-record formation, including merging when overlapping
  events force a wider REF span — freebayes' haplotype-allele
  model (see [calling_pipeline_architecture.md §"Overlapping
  events extend the anchor REF"](calling_pipeline_architecture.md)).
- Five-scalar accumulation per (open record, allele) pair.
- Phase-chain slot allocation, recycling, and lifecycle markers.
- Mate-overlap resolution at per-position granularity.
- Open-record closure (eager — as soon as the record's
  footprint is fully behind the walker) and emission to Stage 2
  in coordinate order. `MAX_RECORD_SPAN` bounds per-record memory by capping
  read span upstream; it is not the closure trigger. See
  §"Closure rule".

Outside this scope:
- CRAM decoding and multi-file merge — handled in
  [`cram_input`](../../src/per_sample_caller/cram_input.rs).
- Per-read filtering (flags, MAPQ, etc.) — runs upstream of the
  walker as part of the parallel filter+BAQ stage
  ([per_sample_caller.md §"Read filters"](per_sample_caller.md)).
- BAQ HMM — runs upstream alongside filtering. The walker only
  ever sees BAQ-capped base qualities. BAQ is its own spec.
- Byte-level encoding to disk — Stage 2 (`psf_writer`).

## Input interface

The walker consumes a stream of `PreparedRead` values produced
by the upstream filter+BAQ stage. Each value is fully
self-describing:

| Field | Type | Notes |
|---|---|---|
| `chrom_id` | `u32` | Index into the file header's chromosome table. |
| `alignment_start` | `u32` | 1-based reference position of the first reference base the read covers. |
| `alignment_end` | `u32` | 1-based reference position of the last reference base the read covers (inclusive). Cached at decode time so the walker does not re-walk the CIGAR. |
| `cigar` | `Vec<CigarOp>` | Owned list of CIGAR ops. |
| `seq` | `Vec<u8>` | Read bases, uppercase ASCII over `{A,C,G,T,N}`. Length matches the CIGAR's read-consuming op count. |
| `bq_baq` | `Vec<u8>` | Per-base BAQ-capped base quality, Phred. Length matches `seq`. |
| `mq_log_err` | `f64` | `ln(P(misalignment))` derived from MAPQ once at filter time. |
| `is_reverse_strand` | `bool` | From the BAM flag. |
| `qname` | `Arc<str>` | Cheap to clone; used by mate-overlap and phase-chain logic. |
| `is_first_mate` | `bool` | From flag `0x40`; tiebreaker on equal-BQ overlaps. |
| `has_mate` | `bool` | From flag `0x1`; controls whether the QNAME is registered for second-mate lookup. |

`CigarOp` is the project's own enum (per [design_principles.md
principle 5](design_principles.md), no noodles types in our
internal interfaces — they're wrapped at the cram_input boundary).

The stream invariants the walker relies on:
- Monotonically non-decreasing `(chrom_id, alignment_start)`.
- Every read has already been through the filter cascade (no
  unmapped, secondary, supplementary, QC-fail, duplicate, or
  low-MAPQ reads survive). The walker does not re-check.
- BAQ has already capped each base quality. The walker does not
  re-apply BAQ.
- `alignment_end ≥ alignment_start`. A read decoding zero
  reference bases is rejected upstream.
- Reads whose reference span exceeds `MAX_RECORD_SPAN` are also rejected
  upstream (counted in the run summary, per
  [calling_pipeline_architecture.md §"Anchor window and
  closure"](calling_pipeline_architecture.md)).

A violation of the order invariant is a hard error — the walker
halts with a message naming the offending read's `qname` and the
out-of-order coordinate, per [design_principles.md principle
3](design_principles.md).

## Output interface

The walker emits a sequence of `PileupRecord` values to the
Stage 2 encoder running on a separate thread, one per closed
record, in coordinate order. Hand-off is via a bounded
channel — see §"Output channel" for the protocol; the walker's
`run` signature takes the channel's `Sender<PileupRecord>` end
and pushes a record through it the moment the record closes.

The shape of `PileupRecord`:

| Field | Type | Notes |
|---|---|---|
| `chrom_id` | `u32` | |
| `pos` | `u32` | Anchor position, 1-based. Stage 2 deltas these into `delta_pos`. |
| `new_chains` | `SmallVec<[SlotId; 4]>` | Slot ids that started since the previous emitted position (in any order; Stage 2 may reorder). |
| `expired_chains` | `SmallVec<[SlotId; 4]>` | Slot ids that ended since the previous emitted position. |
| `alleles` | `Vec<AlleleObservation>` | At least one entry; `alleles[0]` is always REF. |

Each `AlleleObservation`:

| Field | Type | Notes |
|---|---|---|
| `seq` | `SmallVec<[u8; 8]>` | Literal allele over `{A,C,G,T,N}`. For SNP/MNP/DEL alleles the length equals the record's REF-span; INS-bearing alleles are longer. |
| `num_obs` | `u32` | Reads supporting this allele. |
| `q_sum` | `f64` | `Σ max(ln_BQ_BAQ, ln_MQ)` over supporting reads. |
| `fwd` | `u32` | Forward-strand reads in `num_obs`. |
| `placed_left` | `u32` | Reads whose mapped 5′ end is to the left of `pos`. |
| `placed_start` | `u32` | Reads whose mapped 5′ end *is* `pos`. |
| `chain_slots` | `SmallVec<[SlotId; 4]>` | Distinct slot ids that contributed to this allele. |

The record's reference span is not a stored field — it is
**derivable as `alleles[0].seq.len()`**, since `alleles[0]` is
always REF and REF's literal sequence covers the full
reference stretch. Storing it separately would duplicate the
information and risk desync; consumers needing the span call
`record.ref_span()` (a one-line accessor) or read it directly
off `alleles[0]`.

Allele equality (so two reads contribute to the same entry) is
**byte-for-byte identity of `seq`** under a fixed REF-span.
That is the consequence of the haplotype-allele convention:
when the open record's REF-span widens during processing,
every existing allele's `seq` is rewritten against the wider
span before new contributions are folded in, so that
identical-haplotype reads converge on the same `seq`
regardless of arrival order.

The first entry is always REF: literal reference bases over
`[pos, pos + ref_span)` taken from the FASTA. REF is created
when the open record is created, even if no read directly
observes the unmodified reference at that span (reads
supporting any non-REF allele are *not* counted toward REF —
the allele list is partition-style, one bucket per supporting
read). REF entries with `num_obs = 0` are valid and emitted: a
record with only non-REF alleles and a zero-count REF means
no read covering this record's full span matched the
reference.

A zero-observation REF entry carries no data (zero scalars,
empty `chain_slots`), but the entry itself is kept so that
`alleles[0]` is unconditionally REF — the rule that lets us
derive `ref_span` from `alleles[0].seq.len()` and keeps Stage
5's allele-iteration uniform. The byte cost on disk is the
Stage 2 encoder's concern, not the in-memory record's: the
encoder sets a per-allele flag bit indicating
observed/unobserved, and writes the scalar block only for the
observed case. Zero-obs REF entries therefore cost a single
flag bit on disk (plus REF's own `seq`, which is needed for
the span anyway) — no scalar bytes. Non-REF entries always
have `num_obs ≥ 1` under the current walker (allele buckets
are created lazily on first event), so the flag-bit path is in
practice REF-specific.

## Data model

The walker maintains four coupled structures. Each is small in
absolute size: per-record memory is bounded by `MAX_RECORD_SPAN` (the
upstream read-span filter caps any open record's `ref_span` at
`MAX_RECORD_SPAN`), and per-position memory usage is bounded by
typical coverage depth.

### Active read table

`ActiveReads = Vec<ActiveRead>` — a flat vector of reads
currently overlapping the walker's position. A plain `Vec` (vs
a `SmallVec` with an inline buffer) is fine here: there is one
`ActiveReads` per walker, allocated once and reused for the
whole run, so amortised heap-allocation cost is negligible —
and `Vec` carries no fixed cap, so the structure handles
arbitrarily-deep pileups without us having to defend a
particular inline-size choice. The shape:

| Field | Type | Purpose |
|---|---|---|
| `read_id` | `u32` | Local monotonically-increasing id, allocated on entry. Used as the key in mate-overlap lookups and chain bookkeeping. |
| `read` | `PreparedRead` | The owned read content, moved in when added. |
| `events` | `Vec<ReadEvent>` | Pre-computed CIGAR decomposition (see §"Read decomposition"). |
| `event_cursor` | `u32` | Index into `events` of the next unprocessed event — advanced as the walker passes each event's reference position. |
| `chain_slot_id` | `SlotId` | The phase-chain slot this read belongs to. |
| `is_first_mate` | `bool` | Cached from `read` for fast access. |
| `mate_read_id` | `Option<u32>` | When the read's mate is also active, points at it. Filled when the second mate enters. |

`SlotId = u16`, with a hard cap on concurrently active chains
of `MAX_ACTIVE_SLOTS = 4096` (12 of the 16 bits used). `u16`
matches the defensive choice from
[per_sample_caller.md §"Per-record encoding"](per_sample_caller.md)
and gives ~16× headroom over the cap; the per-element memory
cost over `u8` is one byte per slot id, negligible at the
collection sizes involved (`SmallVec<[SlotId; 4]>` carries 8
inline bytes vs 4). The cap is set high enough that typical
30× coverage does not approach it; long-tail repeat regions
that do hit it produce a hard error rather than silent slot
reuse.

The active set is iterated in `read_id` order in a few places
(deterministic record-formation tiebreaks, deterministic
emission order). A flat `Vec` rather than a `HashMap` is
chosen as the *primary* container because iteration is
frequent and Vec iteration is cache-friendly.

For the *lookup* path — given a `read_id`, find the matching
`ActiveRead` (used by mate-overlap fold-time pairing) — the
walker keeps a secondary index:

```
active_by_read_id: AHashMap<u32, usize>   // read_id → index into active_reads
```

Why a secondary index even at small sizes: lookup cost grows
quadratically with coverage. At 5× the active set is ~10
entries and a linear scan is trivially fast; at 30× it is
~60–90 entries and lookups happen on every second-mate read,
so the linear-scan path becomes a measurable fraction of
walker time. The secondary index removes the coverage-scaling
penalty for ~one extra hash insert/remove per read admission
or expiry — well under the cost of the corresponding
CIGAR decomposition and event folding.

The two structures are kept in sync through encapsulation:
admission inserts into both, expiry uses `Vec::swap_remove` on
`active_reads` and updates `active_by_read_id` for both the
removed entry and the entry that took its place. The pattern
is standard and contained in `active_set.rs`; tests assert
`active_by_read_id[r] == active_reads.iter().position(|x|
x.read_id == r)` after each admission and expiry.

### Pending mates map

`pending_mates: AHashMap<Arc<str>, PendingMate>` where
`PendingMate { chain_slot_id: SlotId, first_mate_read_id: u32,
seen_at: u32 }`. The map's purpose is short-term bookkeeping
for paired reads: when the first mate of a pair is admitted
and the second has not yet arrived, the entry holds the slot
allocation and the first mate's `read_id` so that the second
mate (when it arrives) can be linked to the same phase chain
slot and to its partner's `ActiveRead`.

**Unpaired reads are not inserted.** An entry is created
*only* when the read's `has_mate` flag is set (BAM flag
`0x1`). Single-end reads, or paired-end reads whose mate has
been filtered out upstream, never enter the map — they're
allocated a fresh slot directly and treated as solo (the slot
will release when their `ActiveRead` exits, no partner
linkage). This keeps the map sized to genuinely-pending pairs,
not to all in-flight reads.

An entry is dropped as soon as the second mate arrives
(consumed in `allocate_for_read`'s "second-mate path") or
when the defensive window times out (below).

`MATE_LOOKUP_WINDOW`: default `10_000` bp. A read whose mate
has not arrived within 10 kbp past the first-mate's position
is treated as solo — its `chain_slot_id` is finalised, the
map entry is released, and a counter is bumped in the run
summary. The bound is sized for typical Illumina paired-end
fragments (200–500 bp insert; outer-mate distance under a
few kb even on jumping libraries), with comfortable headroom
for unusual but valid cases. It exists to bound map size on
malformed inputs (a stray first-mate-without-partner cannot
grow the map without limit), not to constrain typical
biology. Long-read paired protocols with mate gaps beyond
10 kb would need this raised; the value lives in a named
constant so a CLI override is a one-line addition if real
data demands it.

### Open record table

`OpenPileupRecordTable = BTreeMap<u32, OpenPileupRecord>` keyed by 1-based
anchor position. The shape of `OpenPileupRecord`:

| Field | Type | Purpose |
|---|---|---|
| `pos` | `u32` | Anchor position. Redundant with the map key, kept for ergonomics. |
| `ref_seq` | `SmallVec<[u8; 8]>` | The reference bases over `[pos, pos + ref_seq.len())`, taken from the FASTA on first open and re-fetched / re-checked when the span widens. The current REF-span is `ref_seq.len()` — no separate `ref_span` field is stored. |
| `alleles` | `SmallVec<[OpenAllele; 4]>` | One entry per distinct allele. `alleles[0]` is always REF (`seq == ref_seq`). |
| `new_chain_emit` | `SmallVec<[SlotId; 4]>` | Slots that started since the previous emitted record (filled by the slot allocator, drained on emission). |

`OpenAllele`:

| Field | Type | Purpose |
|---|---|---|
| `seq` | `SmallVec<[u8; 8]>` | Length = `ref_seq.len()` for SNP/MNP/DEL alleles; longer for INS-bearing alleles. |
| `scalars` | `FiveScalars` | The five running per-allele scalars. |
| `chain_slots` | `SmallVec<[SlotId; 4]>` | Sorted, deduplicated. |

`BTreeMap` over `HashMap` because the walker frequently needs to
**find open records whose REF span overlaps a given event
span**, and that query is naturally a small range scan: iterate
`open_records.range(..=event_anchor_pos)` and check the last
few entries for span overlap. Open records live for at most
their own footprint length (see §"Closure rule"), so the map
is shallow — most entries are short SNP records (`ref_span =
1`) that close one walker step after they open. The worst-case
upper bound is the longest open record's `ref_span`, itself
bounded by `MAX_RECORD_SPAN`.

### Output channel

The walker does not buffer closed records internally. As soon
as a record closes (its footprint is fully behind the walker),
it is pushed directly to a **bounded MPSC-shaped channel**
(`crossbeam_channel::bounded` or `std::sync::mpsc::sync_channel`)
whose receiving end is owned by the Stage 2 encoder thread.
The walker's emission step is one `tx.send(record)` per closed
record.

The channel itself supplies all the buffering needed:

- **Bursty closures handled by the buffer.** When the walker
  advances over a stretch where several short records age out
  together (or jumps a long gap of uncovered positions
  finalising many records at once), the walker pushes each
  one in turn; the channel buffer absorbs the burst.
- **Backpressure on a slow encoder.** If Stage 2 falls behind,
  the channel buffer fills, and the next `send` blocks until
  Stage 2 drains an entry. The walker naturally throttles
  itself to the encoder's pace — no extra coordination
  primitive needed.
- **Coordinate order preserved.** Records close in
  monotonically-non-decreasing position order (the front of
  `OpenPileupRecordTable` is closed first), and the channel
  preserves send order, so Stage 2 receives records in
  coordinate order without sorting.

**Channel buffer size.** Default `64` slots. Sized for the
worst-case burst of simultaneous closures (a long uncovered
gap closes every open record at once) plus enough headroom
for a few subsequent steps to proceed before backpressure
kicks in. Small enough that excess memory pressure is
negligible.

**Cross-thread requirement.** Per
[per_sample_caller.md §"Parallelism"](per_sample_caller.md),
the Stage 2 encoder runs on a separate thread from the walker
so that zstd compression overlaps with walker work. The
channel is the only synchronisation point between them; the
walker has no other thread-shared state with Stage 2.
`PileupRecord` is owned (no borrowed data), so it crosses the
thread boundary cleanly.

## Read decomposition

When a read enters the active set, the walker walks its CIGAR
once and produces `events: Vec<ReadEvent>` — the read's ordered
list of per-reference-position contributions. Doing this eagerly
once on entry is simpler than walking lazily as the walker
position advances, because the same data is consulted twice
(once when the walker reaches each event, once during
mate-overlap lookup).

`ReadEvent` is a small enum. All variants carry the 1-based
reference position they apply to (the walker sorts events by
this position).

| Variant | Fields | Meaning |
|---|---|---|
| `Match { ref_pos, base, bq_baq }` | `u32, u8, u8` | The read has `base` aligned to `ref_pos` with BAQ-capped quality `bq_baq`. CIGAR ops `M`, `=`, `X` all produce this. |
| `Insertion { anchor_ref_pos, seq, bq_proxy }` | `u32, SmallVec<[u8;8]>, u8` | An insertion of `seq` is anchored at `anchor_ref_pos` (the reference position immediately before the inserted run). |
| `Deletion { anchor_ref_pos, deleted_len, bq_proxy }` | `u32, u8, u8` | A deletion of `deleted_len` bases starting at `anchor_ref_pos + 1`, anchored at `anchor_ref_pos`. |

Soft-clipped (`S`) and hard-clipped (`H`) bases produce no
events. `N` ops produce no events either — the read does not
observe those reference positions, so the walker simply has no
contribution from this read at those positions (it is *not*
that the read disagrees; it is that the read does not speak to
those positions). The position-emission rule below is what
makes this work correctly.

`P` (padding) ops produce no events.

### CIGAR walk algorithm

```
fn decompose(read: &PreparedRead, fasta: &FastaSlice) -> Vec<ReadEvent> {
    let mut events = Vec::with_capacity(read.cigar.len());
    let mut ref_pos = read.alignment_start;              // 1-based
    let mut read_pos = 0_usize;                          // 0-based into seq

    for (op_idx, op) in read.cigar.iter().enumerate() {
        let is_first = op_idx == 0;
        let is_last  = op_idx + 1 == read.cigar.len();

        match op.kind {
            M | EQ | X => {
                for k in 0..op.len {
                    events.push(Match {
                        ref_pos: ref_pos + k as u32,
                        base:    read.seq[read_pos + k],
                        bq_baq:  read.bq_baq[read_pos + k],
                    });
                }
                ref_pos  += op.len as u32;
                read_pos += op.len as usize;
            }
            I => {
                // Drop the indel only when:
                //   (a) it is the first or last CIGAR op (no flanking evidence — freebayes rule), OR
                //   (b) its anchor falls before position 1 (off the chromosome's start edge).
                if !is_first && !is_last && ref_pos > 1 {
                    let anchor_pos = ref_pos - 1;            // VCF anchor: position before insertion
                    let seq        = read.seq[read_pos .. read_pos + op.len as usize].to_smallvec();
                    let bq         = indel_bq_proxy(&read.bq_baq, read_pos, op.len, op.kind);
                    events.push(Insertion { anchor_ref_pos: anchor_pos, seq, bq_proxy: bq });
                }
                read_pos += op.len as usize;
                // ref_pos unchanged: insertions do not consume reference.
            }
            D => {
                if !is_first && !is_last && ref_pos > 1 {
                    let anchor_pos = ref_pos - 1;
                    let bq         = indel_bq_proxy(&read.bq_baq, read_pos, op.len, op.kind);
                    events.push(Deletion { anchor_ref_pos: anchor_pos, deleted_len: op.len, bq_proxy: bq });
                }
                ref_pos += op.len as u32;
                // read_pos unchanged: deletions do not consume read bases.
            }
            N => { ref_pos += op.len as u32; }            // splice — no event, both ends produce events independently
            S => { read_pos += op.len as usize; }         // soft clip — drop bases
            H | P => { /* nothing */ }
        }
    }
    events
}
```

Edge cases worth pinning:

- **Indel anchored before position 1 (chromosome start edge).**
  An indel is anchored at `ref_pos - 1`, which can fall off
  the chromosome if `ref_pos = 1` when the I/D op is reached.
  Two CIGAR shapes hit this in practice:
  - **Leading I/D (e.g. `2I3M…` at `alignment_start = 1`)** —
    already filtered by the first/last-op rule (the indel is
    op index 0).
  - **Soft-clip-then-indel (e.g. `1S2I5M` at `alignment_start
    = 1`)** — *not* caught by the first/last-op rule (the
    `1S` is op 0, the `2I` is op 1), but the soft clip
    doesn't advance `ref_pos`, so the I still anchors at
    `1 - 1 = 0`. The CIGAR walk above adds the bound check
    `ref_pos > 1` to drop the indel in this case. Only the
    indel event is dropped — subsequent `M` ops on the same
    read still produce normal `Match` events at valid ref
    positions.

  This is an extremely rare case in real data — most aligners
  either soft-clip more aggressively or move `alignment_start`
  past the inserted bases — but the bound is cheap and
  prevents an internal underflow on `ref_pos - 1`.

- **A read with CIGAR `5M2D5M` ending at the deletion's last
  base.** The `D` is op-index 1, not last (the trailing `M` is
  last), so it survives. The trailing `M` provides flanking
  evidence as required.

- **Consecutive indels: CIGAR `3M2I3D2M`.** Two adjacent indel
  events. Decompose each as its own event. Anchor position for
  the insertion: end of first `M` block - 1. Anchor position
  for the deletion: same position (no reference bases between
  them). Both anchored at the same reference position.
  Open-record merging (next §) folds them.

### Indel BQ proxy

Per [calling_pipeline_architecture.md §"Indel BQ
proxy"](calling_pipeline_architecture.md), the per-read
quality contribution to an indel allele is the **minimum**
BAQ-adjusted base quality over a window of `l + 2` bases
centred on the indel, edge-clamped at the read ends.
Confirmed against freebayes: same window length and same `min`
aggregation, see
[AlleleParser.cpp:1638](../../freebayes/src/AlleleParser.cpp#L1638)
(deletions, `int L = l + 2`),
[AlleleParser.cpp:1659](../../freebayes/src/AlleleParser.cpp#L1659)
(`minQuality(qualstr)` when `useMinIndelQuality` is set, the
default), and
[AlleleParser.cpp:1714](../../freebayes/src/AlleleParser.cpp#L1714)
(insertions). Freebayes exposes a `-H --harmonic-indel-quality`
flag to swap `min` for a sum-quality alternative; we do not
expose this — the architecture commits to `min`.

Window centring differs slightly between insertion and deletion
in freebayes (insertion centres on `rp - 1`, deletion on `rp -
L/2`). Both edge-clamp at read ends. We follow the same shape
in `indel_bq_proxy` below.

```
const INDEL_BQ_PROXY_PAD: usize = 1;   // bases on each side of the indel run

fn indel_bq_proxy(bq_baq: &[u8], read_pos: usize, op_len: u32, op: CigarKind) -> u8 {
    let (lo, hi) = match op {
        CigarKind::I => {
            // inserted bases occupy bq_baq[read_pos .. read_pos + op_len]
            let lo = read_pos.saturating_sub(INDEL_BQ_PROXY_PAD);
            let hi = (read_pos + op_len as usize + INDEL_BQ_PROXY_PAD).min(bq_baq.len());
            (lo, hi)
        }
        CigarKind::D => {
            // deleted bases have no read bytes; pad ±1 around the deletion's read-position cursor
            let lo = read_pos.saturating_sub(INDEL_BQ_PROXY_PAD);
            let hi = (read_pos + INDEL_BQ_PROXY_PAD).min(bq_baq.len());
            (lo, hi.max(lo + 1))   // at least one base; should always hold given filters
        }
        _ => unreachable!(),
    };
    bq_baq[lo..hi].iter().copied().min().unwrap()
}
```

`INDEL_BQ_PROXY_PAD = 1` (one flanking base on each side) is what
yields freebayes' `l + 2` window. The constant lives at module
scope with a doc comment noting the source (architecture spec +
freebayes `AlleleParser.cpp:1626`).

## The walker loop

The loop runs single-threaded and is the sequential backbone of
Stage 1. The walker takes a `Sender<PileupRecord>` — one end of
the bounded channel described in §"Output channel" — and pushes
each closed record through it as it closes. Pseudocode, with
`walker_pos` as a 1-based reference position within the current
chromosome:

```
fn run(
    reads: impl Iterator<Item=PreparedRead>,
    fasta: &Fasta,
    tx: &Sender<PileupRecord>,
) -> Result<()> {
    let mut state = WalkerState::new();
    let mut reads = reads.peekable();

    while let Some(peek) = reads.peek() {
        // Chromosome boundary: flush everything from the current chromosome.
        if peek.chrom_id != state.chrom_id {
            state.flush_all(tx, fasta)?;
            state.advance_to_chrom(peek.chrom_id);
            continue;
        }

        // Pull every read that starts at or before the next position we want
        // to process (walker_pos), so the active set is complete for that pos.
        while reads.peek().map_or(false, |r| r.alignment_start <= state.walker_pos) {
            let read = reads.next().unwrap();
            state.admit_read(read, fasta)?;
        }

        // Process events whose ref_pos == walker_pos across the active set,
        // forming or extending open records as needed.
        state.process_position(walker_pos = state.walker_pos, fasta)?;

        // Close records whose footprint is fully behind the walker; push each
        // through `tx`. Send blocks if the channel buffer is full (backpressure).
        state.close_aged_records(tx)?;

        // Expire reads that no longer cover walker_pos + 1 onward.
        state.expire_passed_reads();

        // Advance to the next position. If the active set is empty AND the
        // upstream still has reads, jump to the next read's start (skip uncovered).
        state.walker_pos = state.next_walker_pos(reads.peek())?;
    }

    state.flush_all(tx, fasta)
    // dropping `tx` on return closes the channel; Stage 2's recv loop exits.
}
```

A few subtleties:

- **`next_walker_pos`** returns `walker_pos + 1` if any active
  read still has events at or beyond `walker_pos + 1`, otherwise
  the next pulled read's `alignment_start`. Uncovered positions
  produce no record (per
  [calling_pipeline_architecture.md §"Per-position record layout
  and compression"](calling_pipeline_architecture.md)). Skipping
  empty positions is what keeps off-target regions cheap.

- **`process_position`** is where event folding, open-record
  formation, and mate-overlap resolution happen. It is the
  algorithmic heart of the walker; expanded in §"Open-record
  formation and merging" below.

- **`close_aged_records`** scans the front of `OpenPileupRecordTable`
  for records whose footprint is fully behind the walker — i.e.
  `pos + ref_span ≤ walker_pos` — converts each to a
  `PileupRecord`, and immediately sends it through `tx`. See
  §"Closure rule" below for why this per-record condition is
  sufficient (and why `MAX_RECORD_SPAN` is not the closure
  trigger). The conversion is mechanical: open alleles become
  `AlleleObservation`s; `chain_slots` are sorted and deduped;
  `new_chain_emit` and the slot allocator's `expired_chain_emit`
  populate the record's `new_chains` / `expired_chains`. Each
  `tx.send` may block if Stage 2 has fallen behind and the
  channel buffer is full — that's the intended backpressure.

- **`flush_all`** is the same as `close_aged_records` but
  unconditional — every open record is closed and sent through
  `tx`. Called at chromosome boundaries and at end-of-input.

- **Chromosome boundaries.** Open records never span
  chromosomes (per
  [per_sample_pileup_format.md §"Block layout"](per_sample_pileup_format.md))
  and the `BTreeMap` is reset at each chromosome change. The
  slot allocator is also reset (no chain reaches across
  chromosomes — reads are aligned to one chromosome).

### Closure rule

An open record at position `Q` with current `ref_span = s` is
safe to close as soon as
`walker_pos ≥ Q + s` — i.e., the walker has advanced past the
record's last covered reference position. The reasoning:

1. Reads arrive coordinate-sorted by `alignment_start`, and the
   walker admits every read with `alignment_start ≤ walker_pos`
   before processing the position. Future (not-yet-admitted)
   reads therefore have `alignment_start ≥ walker_pos + 1`.
2. Every event's anchor sits at an `M` position of its read, so
   for any future event the anchor `S ≥ alignment_start ≥
   walker_pos + 1`.
3. A future event with anchor `S` can merge into the record at
   `Q` only if its footprint `[S, S + size)` overlaps the
   record's current footprint `[Q, Q + s)`. That requires
   `S < Q + s`. Combined with `S ≥ walker_pos + 1`, the merge
   is possible only when `walker_pos + 1 < Q + s`.
4. So when `walker_pos ≥ Q + s`, no future event can overlap,
   and the record is safe to close.

Concretely: a SNP record (`ref_span = 1`) at position `Q`
closes the moment the walker advances to `Q + 1`. A 6 bp
deletion record stays open for six walker steps. The walker
never holds records open longer than their own footprint
demands.

#### Why `MAX_RECORD_SPAN` is not the closure trigger

`MAX_RECORD_SPAN` (default 5000 bp) is **the upstream read-span filter**, not
the closure rule. A read whose reference span (CIGAR walk plus
mate gap for paired reads) exceeds `MAX_RECORD_SPAN` is dropped before it
ever reaches the walker (per
[calling_pipeline_architecture.md §"Reads exceeding the
window"](calling_pipeline_architecture.md)). This guarantees
two things downstream:

- Every event's footprint extends at most `MAX_RECORD_SPAN` bases past its
  anchor, so no single event can grow an open record's
  `ref_span` past `MAX_RECORD_SPAN`.
- The number of positions the walker has to remember per open
  record is bounded by `MAX_RECORD_SPAN`, so per-record memory is bounded.

The eager closure rule above plus the `MAX_RECORD_SPAN`
read-filter together bound walker memory to roughly
`O(MAX_RECORD_SPAN × max_concurrent_open_records)`, which is
what the architecture spec's bounded-memory promise delivers
([calling_pipeline_architecture.md §"Memory
bound"](calling_pipeline_architecture.md)).

## Open-record formation and merging

The walker's central operation: at `walker_pos`, fold all events
at that position from all active reads into the right open
record(s), opening new ones and merging existing ones as needed.

### Step 1: collect events at this position

```
struct EventAtPos {
    read_id: u32,
    event:   ReadEvent,
}

let mut hits: Vec<EventAtPos> = Vec::new();
for read in &mut active_reads {
    while let Some(ev) = read.events.get(read.event_cursor as usize) {
        let evpos = match ev {
            Match { ref_pos, .. }                    => *ref_pos,
            Insertion { anchor_ref_pos, .. }         => *anchor_ref_pos,
            Deletion { anchor_ref_pos, .. }          => *anchor_ref_pos,
        };
        if evpos != walker_pos { break; }
        hits.push(EventAtPos { read_id: read.read_id, event: ev.clone() });
        read.event_cursor += 1;
    }
}
```

`hits` is the per-position bundle of contributions from all
active reads. Note that `Insertion` and `Deletion` events fire
at their **anchor** position (one before the variable region),
so they are collected at `walker_pos = anchor_ref_pos`, not at
their interior positions.

### Step 2: resolve mate overlap

For every pair of `(read_id, mate_read_id)` in `hits` whose
events at `walker_pos` are both `Match` (the only event type
where mate-overlap matters — indels at the same anchor on both
mates are handled by deduplication, see §"Mate overlap on
indels"):

1. Compare `bq_baq` of the two `Match` events.
2. Higher BQ wins: keep that event's contribution unchanged.
3. Lower BQ loser: set its `bq_baq` to 0 in the local copy
   used for accumulation. The observation count still
   increments by 1 — the read *did* observe this position —
   but its log-likelihood mass is `ln(1) = 0`.
4. Tie: keep mate 1 (`is_first_mate`).

Disagreement (different bases) follows the same rule: higher
BQ wins, no extra downgrade
([per_sample_caller.md §"Disagreement at an
overlap position"](per_sample_caller.md)).

#### Mate overlap on indels

Two cases:

- **Both mates report the same indel at the same anchor.** Treat
  as one observation, not two. Assign the higher-BQ-proxy event
  to the allele bucket; drop the other. The allele's
  `chain_slots` still gets the (single, shared) `chain_slot_id`
  for the pair. The per-sample caller spec covers per-position
  SNP overlap explicitly but is silent on indel overlap; treating
  the pair as one observation is the consistent extension of
  that rule.

- **Mates disagree on indel presence at the same anchor** (one
  reports the indel, the other reports a clean `Match`). The
  higher-BQ contribution wins; the loser's contribution at this
  position is dropped. If they tie, prefer mate 1.

### Step 3: identify candidate open records per event

Each event has a candidate REF-span footprint:

- `Match` at `p`: footprint is `[p, p+1)`, ref length 1.
- `Insertion` at `p`: footprint is `[p, p+1)`, ref length 1
  (insertions do not consume reference; the inserted bases
  ride on the anchor base).
- `Deletion` at `p`, length `l`: footprint is `[p, p+l+1)`,
  ref length `l + 1` (anchor base + deleted bases).

A new event merges into an existing open record `A` iff the
event's footprint **overlaps** `A`'s current `[pos, pos +
ref_span)`. "Overlaps" here means non-empty interval
intersection (touching intervals — `A.end == event.start` —
are *not* overlapping; they remain separate open records linked
only by phase chain).

Per event, the walker queries `open_records.range(..=event_anchor_pos)`
and checks the last few entries for overlap. In practice only
the immediately preceding open record is a candidate at typical
spacings; a hit further back means we are inside the REF span
of an earlier-opened deletion-extending record, which is rare
but allowed.

If no overlap is found, a fresh open record is created at the
event's anchor position with `ref_span = 1` for
`Match`/`Insertion` or `ref_span = l + 1` for `Deletion`.

### Step 4: fold the event into the chosen open record

If `A` is the chosen open record and `e` is the event:

1. **Widen `A`'s REF span if needed.** If the event's
   footprint extends past `A.pos + A.ref_seq.len()`, refetch
   `A.ref_seq` from the FASTA over the wider span (so
   `A.ref_seq` now covers `[A.pos, event_end)`). The current
   span is always `A.ref_seq.len()`; no separate field tracks
   it.

2. **Rewrite all existing alleles to the new span.** For each
   `OpenAllele` in `A.alleles`:
   - If the allele's previous coverage already ends at
     `A.pos + old_span`, append the reference bases at
     `[A.pos + old_span, A.pos + new_span)` to its `seq`.
   - This is a noop when the span is unchanged.
   - Deletion alleles whose previous representation covered
     the full old span keep their existing length (they're
     still the same biological deletion); the *REF* widens to
     cover more bases, so the deletion's `seq` (anchor base
     only, for example) is shorter than the new REF,
     expressing "more bases deleted relative to the wider
     ref" — exactly freebayes' haplotype-allele convention.

   The rewrite for deletion alleles is the subtle one. Worked
   example: existing open record at 100, REF = "AT", allele
   seqs = ["AT" REF, "A" del-of-101]. Now a new event widens
   REF to "ATC" (covers 100-102). The REF allele becomes
   "ATC". The deletion allele becomes "A" (still anchor base
   only — the deletion now spans 101-102, and the new REF base
   at 102 is also "deleted" relative to the wider span). This
   matches freebayes' equivalent CIGAR-merge construction in
   [Allele.cpp:1454-1470](../../freebayes/src/Allele.cpp#L1454).

3. **Compute `e`'s allele seq under the (possibly widened) `A`.**
   - `Match` at offset `(p - A.pos)` within REF: allele seq is
     `A.ref_seq` with the matching position replaced by the
     read's base. If the base equals the ref base, this is
     REF, not a new allele.
   - `Insertion` of bases `i_seq` at offset `(p - A.pos)`:
     allele seq is `A.ref_seq[..=offset] + i_seq + A.ref_seq[offset+1..]`.
   - `Deletion` of length `l` starting at `p + 1` within
     `A.ref_seq`: allele seq is `A.ref_seq[..=offset]` plus
     `A.ref_seq[offset + 1 + l ..]` (everything outside the
     deleted run).

4. **Locate or create the allele bucket.** Linear scan of
   `A.alleles` for a bucket with matching `seq`; create one if
   none exists. New allele buckets get `scalars =
   FiveScalars::zero()` and an empty `chain_slots`.

5. **Update the allele bucket's scalars.** For the
   contributing read `r`:
   ```
   allele.scalars.num_obs       += 1
   allele.scalars.q_sum         += max(ln_bq_for_read(r, A), r.mq_log_err)
   allele.scalars.fwd           += if r.is_reverse_strand { 0 } else { 1 }
   allele.scalars.placed_left   += if r.alignment_start < A.pos { 1 } else { 0 }
   allele.scalars.placed_start  += if r.alignment_start == A.pos { 1 } else { 0 }
   ```
   `ln_bq_for_read(r, A)` is the read's per-allele BQ
   contribution under this open record in log-error space. It
   is computed in two steps (see §"Compound-event quality"
   below): first the read's per-event BQs over the record's
   span are reduced via `min` to a single per-read BQ; then
   that is converted to `ln(P_err)`. The architecture-mandated
   `max(ln_BQ, ln_MQ)` per-read combination happens here, at
   the point of folding into `q_sum` — not at the
   haplotype-allele construction step. A per-read BQ of 0 maps
   to `ln(1) = 0`.

6. **Add the read's `chain_slot_id` to `allele.chain_slots`** (sorted
   insert + dedup; the slot id may already be present if the
   same read pair contributed an earlier event to this same
   allele in this open record — possible for compound
   haplotypes where both mates support the same combined
   allele).

### Compound-event quality

When a read contributes a compound allele to an open record —
its support comes from multiple events under the record's span
(a SNP and a deletion in the same open record, two adjacent
indels, etc.) — the per-read BQ for that allele is the
**minimum** BQ across the read's events in the open record's
span:

```
ln_bq_for_read(r, A) = ln(P_err(min over events in r ∩ [A.pos, A.pos + A.ref_seq.len()) of bq))
```

For each event the underlying `bq` is `e.bq_baq` for `Match` or
`e.bq_proxy` for indels (the indel BQ proxy already specified).
For `Match` events outside the open record's span, no contribution.

Confirmed against freebayes
[AlleleParser.cpp:3151-3155](../../freebayes/src/AlleleParser.cpp#L3151):
when haplotype alleles are merged, freebayes takes
`min(quality)` and `max(lnquality)` across the constituent
events — equivalent statements (lower BQ = higher
log-error-prob), so freebayes' rule reduces to "weakest base
dominates." This matches our `min` rule.

The architecture-level `max(ln_BQ, ln_MQ)` combination
([calling_pipeline_architecture.md §"The five per-allele
scalars"](calling_pipeline_architecture.md)) is applied
**after** this per-read BQ reduction, at the point of folding
into `q_sum` (Step 5 of "Step 4: fold the event into the chosen
open record"). Freebayes applies the equivalent combination one
layer further out, in
[DataLikelihood.cpp:34](../../freebayes/src/DataLikelihood.cpp#L34),
per-read at likelihood time; we pre-sum into the scalar at
Stage 1, which is what makes the `q_sum` aggregate exactly
reconstruct freebayes' `prodQout` term.

### Step 5: handle compound events from a single read

A read can contribute multiple events that all fall under the
same open record (e.g. SNP at 101 + insertion at 101 in one
read, both anchored at 101). The walker must fold these as **one
combined allele**, not two separate alleles each missing the
other's modification.

Implementation: gather all events for read `r` whose
anchor-or-ref-position fall inside the candidate open record's
(post-widen) REF span, and apply them all to `A.ref_seq`
together to compute `r`'s haplotype string under `A`. Then
locate/create the allele bucket once and update scalars once.

Pseudocode:

```
for read in &active_reads {
    // events_in_window is read's events whose footprints
    // overlap [A.pos, A.pos + A.ref_seq.len()):
    let window_events = collect_events_overlapping(read, A);
    if window_events.is_empty() { continue; }
    let allele_seq = apply_events_to_ref(A.ref_seq, A.pos, &window_events);
    let allele     = A.alleles.find_or_create(allele_seq);
    let ln_bq      = ln_bq_for_read(read, &window_events);   // see §"Compound-event quality"
    update_scalars(allele, read, ln_bq);
    allele.chain_slots.insert_sorted_unique(read.chain_slot_id);
}
```

`ln_bq_for_read` reduces the read's events to a single per-read
BQ via `min` (see §"Compound-event quality" above) before the
scalar update applies `max(ln_BQ, ln_MQ)`.

### Step 6: REF bucket maintenance

After folding all events, the walker checks every `ActiveRead`
that overlaps `A` but contributed no event in this open
record's span (i.e. the read's haplotype string under `A`
equals `A.ref_seq`). Each such read contributes one observation
to the REF bucket.

This is what makes `num_obs` for REF correct: REF is "reads that
match the reference across the open record's full span", not
"reads that just happen to match at the anchor position." A
read with a SNP at 102 inside an open record that spans 100-105
is *not* a REF contributor — it goes into its own
SNP-at-102-rewritten-against-the-wider-span allele bucket.

Computationally this is folded into Step 5 — every active read
gets considered, including those whose `window_events` is
empty: their `allele_seq` is `A.ref_seq` and they fold into the
REF bucket via the same find-or-create path.

## Phase chain slot allocator

The slot allocator and the walker are mostly independent: the
allocator doesn't know what an open record is, only that reads
enter and leave the active set and that each entry/exit
produces lifecycle marks the walker stamps on the next emitted
record.

### State

```
struct SlotAllocator {
    free:                Vec<SlotId>,         // pool of recycled slot ids, sorted descending so pop() yields the lowest
    next_fresh:          SlotId,              // never-used-before slot id, monotonically increasing within the limit
    pending_mates:  AHashMap<Arc<str>, PendingMate>,
    new_marks:           SmallVec<[SlotId; 4]>,
    expired_marks:       SmallVec<[SlotId; 4]>,
}
```

`free` stores recycled slot ids; popping yields the smallest
free id (Vec sorted descending so `pop()` is `O(1)` for the
smallest). `next_fresh` advances when `free` is empty.

### Allocation on read entry

```
fn allocate_for_read(&mut self, read: &PreparedRead) -> SlotId {
    if read.has_mate {
        if let Some(pending) = self.pending_mates.remove(&read.qname) {
            return pending.chain_slot_id;        // second-mate path: reuse the first mate's slot
        }
    }
    let slot = self.free.pop().unwrap_or_else(|| {
        let s = self.next_fresh;
        self.next_fresh += 1;
        if self.next_fresh > MAX_ACTIVE_SLOTS {
            return Err(WalkerError::SlotExhausted);  // hard error, see below
        }
        s
    });
    if read.has_mate {
        self.pending_mates.insert(
            read.qname.clone(),
            PendingMate { chain_slot_id: slot, first_mate_read_id: read.read_id, seen_at: read.alignment_start },
        );
    }
    self.new_marks.push(slot);
    slot
}
```

### Release on read expiry

A slot is released only when **both** mates (or the sole mate of
a solo read) have exited the active set. The walker tracks
per-slot reference counts via a small `slot_refcount: [u8;
MAX_ACTIVE_SLOTS]` array (or a sparse map at higher caps).

```
fn release_slot(&mut self, slot: SlotId) {
    self.slot_refcount[slot as usize] -= 1;
    if self.slot_refcount[slot as usize] == 0 {
        self.free.push_sorted_descending(slot);
        self.expired_marks.push(slot);
    }
}
```

Refcount goes up by 1 on each `allocate_for_read` call (so a
solo read has refcount 1; a pair has refcount 2 after the second
mate is admitted). Refcount goes down by 1 on each read expiry.

### Stamping records

When the walker closes a record at `pos`, it drains the
allocator's `new_marks` and `expired_marks` since the last
emitted record into `PileupRecord.new_chains` and
`expired_chains` respectively. The drain is per-emission (the
walker emits in coordinate order; new/expired marks accumulate
between emissions and reset at each emission).

**[QUESTION]** what if a slot is allocated and released between
two consecutive emitted records (a very short single-read pile
with no surviving open record)? The slot appears in both `new_chains`
and `expired_chains` of the same record. Stage 5's consumer
needs to apply expired *after* new (or just suppress this case
at emit time). Suppress here: if a slot id is present in both
`new_marks` and `expired_marks`, drop it from both — the chain
existed only between records and was never visible.

### Defensive bounds

- `MAX_ACTIVE_SLOTS`: hard cap, default `4096`. Exceeding it is
  a hard error: very high coverage at a repeat region implies
  the run is likely processing pathological input. The cap is
  user-tuneable via `--max-active-chain-slots`.

- `MATE_LOOKUP_WINDOW`: if a pending first-mate entry has
  `seen_at + MATE_LOOKUP_WINDOW < walker_pos`, the entry is
  evicted: the first mate is finalised as solo
  (`pending_mates` cleared, slot refcount stays at 1, the
  slot will release when the first mate's `ActiveRead` exits).
  An eviction counter is bumped in the run summary.

## Mate-overlap mechanics

The walker resolves mate overlap **at folding time** rather than
**at admission time**. That is: when a read enters the active
set, its events are decomposed and stored verbatim — the
walker does not pre-zero any BQ. Only when the walker reaches a
position where both mates contribute does it run the BQ
comparison and zero the loser locally for that position.

Why fold-time rather than admission-time:

- The first-mate reader may not yet know whether the second
  mate exists; without that information, it cannot decide
  whether to zero anything.
- Fold-time resolution keeps the active set's events a
  faithful per-read record, which makes debugging and
  invariant-checking simpler.
- Cost is identical: the comparison happens once per overlap
  position regardless of when it's done.

### Mate lookup at fold time

When `process_position` collects events at `walker_pos`, it
groups them by `read_id`. After grouping, for each read with
both an event at this position **and** an active mate (via
`mate_read_id`), it pairs the events and runs the comparison.

A read's `mate_read_id` is filled when the second mate enters:
the slot allocator finds the first-mate entry in
`pending_mates`, the active-set admission code looks up
the first mate's `read_id`, and both mates' `mate_read_id`
fields are set to point at each other.

### Tie-breaking and disagreement

- BQ tie: prefer mate 1 (flag `0x40`, `is_first_mate` set).
- Different bases at same position, BQ tie: prefer mate 1.
- Different bases at same position, BQ differ: higher BQ wins;
  no extra downgrade.

## Active-set bookkeeping

A read enters the active set when it is admitted from the
upstream stream. It leaves when it can no longer contribute — at
the point where `walker_pos > active_read.read.alignment_end` AND the
walker has processed all events on the read.

Implementation:

- **Entry:** push to `active_reads` with `event_cursor = 0`;
  run `decompose` to fill `events`; allocate slot via the
  allocator; insert `(read_id → new_index)` into
  `active_by_read_id`; if the read's mate is in
  `pending_mates`, link the two reads' `mate_read_id`
  fields (the mate's index is found via `active_by_read_id`,
  not a linear scan).

- **Exit:** when `walker_pos > active_read.read.alignment_end`,
  the read has no remaining events (its event cursor must
  equal `events.len()`; if not, that's an internal-invariant
  violation and the walker panics with read context). Release
  the slot via the allocator. Use `Vec::swap_remove` on
  `active_reads` to drop the read in `O(1)`, then update
  `active_by_read_id`: remove the exiting read's entry, and
  rewrite the entry of the read that was moved into its slot
  (the previous last element of `active_reads`) to point at
  the new index.

The active set is therefore not "every read that overlaps
walker_pos" but "every read that may still contribute." An
over-the-edge read whose alignment ends at walker_pos contributes
its last event at walker_pos and exits at walker_pos + 1 — which
is the right behaviour given our 1-based inclusive `alignment_end`.

## Errors

Per [design_principles.md principle 3](design_principles.md),
every error halts:

| Cause | Variant | Context fields |
|---|---|---|
| Out-of-order read input | `WalkerError::OutOfOrder` | qname, expected min position, actual position |
| Read decoded zero ref bases | `WalkerError::ZeroRefSpan` | qname |
| Slot pool exhausted (`> MAX_ACTIVE_SLOTS`) | `WalkerError::SlotExhausted` | walker_pos, current active read count |
| Open record REF span exceeded `MAX_RECORD_SPAN` (should not happen — upstream rejects oversize reads) | `WalkerError::RecordTooWide` | anchor pos, current span |
| Internal invariant: event cursor not at end on read exit | `WalkerError::Internal` | qname, residual events |
| FASTA fetch failure during open-record widening | `WalkerError::Fasta(io::Error)` | chrom, requested span |

`WalkerError` is a `thiserror`-built enum, per
[design_principles.md principle 6](design_principles.md). The
CLI edge wraps it with `anyhow::Context` carrying the input
sample name and chromosome.

## Run-summary counters

Reported on stderr at end of run, per
[per_sample_caller.md §"Errors"](per_sample_caller.md):

| Counter | Meaning |
|---|---|
| `reads_seen` | reads pulled from upstream |
| `reads_admitted` | reads that entered the active set |
| `mate_lookup_evictions` | first-mate entries timed out at `MATE_LOOKUP_WINDOW` |
| `slot_allocations` | total slot allocations (including reuses across pairs) |
| `slot_high_water` | maximum concurrently-active slots observed |
| `records_emitted` | total `PileupRecord`s sent to Stage 2 |
| `record_widen_events` | times an open record's REF span widened |
| `mate_overlap_positions` | per-position BQ comparisons performed |

## Open decisions and questions

Resolved against freebayes (consultation 2026-05-06,
`AlleleParser.cpp` and `Allele.cpp`):

- **Compound-event quality rule.** Resolved as `min(BQ over
  events in the open record's span)` for the per-read BQ, then
  `max(ln_BQ, ln_MQ)` for the `q_sum` update. Matches freebayes'
  haplotype-allele construction in
  [AlleleParser.cpp:3151-3155](../../freebayes/src/AlleleParser.cpp#L3151).
  The architecture-level MQ combination happens at the scalar-update
  step, not at the haplotype-allele construction step.
- **Indel BQ proxy details.** Window of `l + 2` and `min`
  aggregation confirmed against
  [AlleleParser.cpp:1638](../../freebayes/src/AlleleParser.cpp#L1638)
  / [:1659](../../freebayes/src/AlleleParser.cpp#L1659)
  / [:1714](../../freebayes/src/AlleleParser.cpp#L1714). Insertion
  vs deletion centring differs slightly in freebayes; we follow
  the same shape.
- **First/last CIGAR-op indel rejection.** Drop *only the
  unflanked indel*, not the read. Other observations from the
  same read still contribute. Confirmed against
  [Allele.cpp:1511-1525](../../freebayes/src/Allele.cpp#L1511)
  (`isUnflankedIndel`) and
  [AlleleParser.cpp:3137](../../freebayes/src/AlleleParser.cpp#L3137).
- **Open-record widening: rewrite-of-existing-alleles model.**
  Freebayes maintains alleles as `(alternate_seq, cigar)` pairs
  and merges via CIGAR concatenation
  ([Allele.cpp:1454-1470](../../freebayes/src/Allele.cpp#L1454));
  our literal-string rewrite produces equivalent alleles. The
  worked example in §"Step 4" stands.
- **Zero-observation REF entries.** Kept in the emitted record
  so that `alleles[0]` is unconditionally REF and `ref_span`
  remains derivable as `alleles[0].seq.len()`. The byte cost on
  disk is the Stage 2 encoder's concern: a per-allele
  observed/unobserved flag bit elides the scalar block for
  zero-obs entries. See §"Output interface" for the full
  argument.
- **`SlotId` underlying type.** `u16` with `MAX_ACTIVE_SLOTS =
  4096` (12 bits used, ~16× headroom over the cap). `u8` was
  rejected: the per-element memory saving over `u16` is one
  byte and the entire stored collection of slot ids is small,
  while a 256-slot cap is uncomfortably close to realistic
  pileup depths and would risk silent overflow in repeat
  regions. The cap can still be raised if needed; the type
  choice carries no further constraint.
- **Active-set mate-lookup index.** Maintained from day one as
  a secondary `AHashMap<u32, usize>` (`read_id → vec_index`)
  alongside the primary `Vec<ActiveRead>`. Mate-lookup cost
  scales quadratically with coverage under linear scan — fine
  at 5× but a real fraction of walker time at 30× — and the
  secondary-index discipline (admission inserts both;
  expiry uses `swap_remove` and updates the moved entry) is
  small in a tested codebase. See §"Active read table" for
  the data structure and §"Active-set bookkeeping" for the
  entry/exit pattern.
- **Indel anchored before position 1.** Drop the indel event;
  the rest of the read still contributes normally. Hit by the
  rare `1S2I…` (or similar) shape at `alignment_start = 1`,
  where the leading soft clip leaves `ref_pos = 1` when the
  walker reaches the I, so the anchor would land at
  position 0 (off the chromosome). The CIGAR walk's
  `ref_pos > 1` guard catches this; see §"Edge cases worth
  pinning" for the worked example.
- **Mate overlap on indels.** When both mates of a pair report
  the same indel at the same anchor, treat as a single
  observation: keep the higher-BQ-proxy event, drop the other.
  When they disagree (one reports the indel, the other a
  clean `Match`), the higher-BQ contribution wins. Freebayes
  has no explicit mate-overlap logic, so there's no upstream
  precedent — this is a deliberate project choice extending
  the per-position SNP overlap rule from
  [per_sample_caller.md §"Mate-overlap handling"](per_sample_caller.md)
  to indels. See §"Mate overlap on indels" for the cases.

No questions remain open before implementation.

## Module layout

The pileup walker is one Rust module
(`src/per_sample_caller/pileup_walker.rs`), but internally
splits into a few small files for readability:

```
src/per_sample_caller/
  pileup/
    mod.rs               — public entry: `run(reads, fasta, tx)`
    walker_state.rs      — WalkerState, the loop driver
    active_set.rs        — ActiveReads, admission/exit, mate linking
    decompose.rs         — CIGAR → ReadEvent
    anchor.rs            — OpenPileupRecord, OpenAllele, open-record merging + widening
    five_scalars.rs      — FiveScalars + accumulation
    slot_allocator.rs    — SlotAllocator, lifecycle marks
    mate_overlap.rs      — fold-time BQ comparison
    errors.rs            — WalkerError
```

This is a sketch; the implementation plan
([feature_implementation_plans/](../feature_implementation_plans/))
will commit the final layout. Per [design_principles.md
principle 4](design_principles.md), the public entry takes
already-prepared reads and a channel `Sender<PileupRecord>` —
the file-on-disk dance and the cross-thread orchestration
between the walker thread and the Stage 2 encoder thread live
in `per_sample_caller/mod.rs`.

## Cross-references

- [calling_pipeline_architecture.md](calling_pipeline_architecture.md)
  — Stage 1 architectural commitments: BAQ in-process, five
  scalars, `MAX_RECORD_SPAN` cap, indel BQ proxy, indel anchoring,
  overlapping-event REF extension.
- [per_sample_caller.md](per_sample_caller.md) — Stage 1
  end-to-end shape: input/output, read filters, mate-overlap
  rule, phase-chain slot encoding, parallelism, CLI.
- [phase_chain.md](phase_chain.md) — conceptual model for what
  phase chains represent and what Stage 5 does with them.
- [per_sample_pileup_format.md](per_sample_pileup_format.md) —
  Stage 2 byte layout downstream of this walker.
- [design_principles.md](design_principles.md) — clarity, error
  handling, side-effect placement, naming.
- [freebayes_posterior_gt_probs.md](freebayes_posterior_gt_probs.md)
  — freebayes' formulas the walker's output is sufficient to
  reproduce.
