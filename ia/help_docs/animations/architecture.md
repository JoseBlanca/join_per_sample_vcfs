# Pileup walker — architecture and code-walk

Companion document to [`index.html`](index.html). The HTML page is an
interactive viewer that steps through three scenarios; this document
is the "manual" — what each Rust submodule does, how they hand off
between each other, and where in the source the runtime moments live.

The audience is someone comfortable with Python who wants to make
sense of the Rust code. Where Rust idioms matter (ownership,
`BTreeMap`, `mpsc::SyncSender`), they are flagged inline.

> **Where this fits in the bigger pipeline.** The pileup walker is
> Stage 1 of the per-sample caller. It reads a coordinate-sorted
> stream of `PreparedRead` values (already filtered + BAQ-adjusted
> upstream) and emits a stream of `PileupRecord` values — one per
> covered reference position. The full pipeline is described in
> [`ia/specs/calling_pipeline_architecture.md`](../../specs/calling_pipeline_architecture.md);
> the implementation-ready spec for this stage is
> [`ia/specs/pileup_walker.md`](../../specs/pileup_walker.md).

## Mental model in one paragraph

The walker advances along the reference one position at a time
(`walker_pos`). At every step it (1) admits any new reads whose
`alignment_start ≤ walker_pos`, (2) asks each active read's CIGAR
cursor "do you have any events anchored at `walker_pos`?", (3)
folds those events into a per-position record (opening a new one
or extending an existing one), (4) drops any read the walker has
fully passed, (5) emits any record whose footprint is now strictly
behind the walker, and (6) advances. A "record" can hold open for
several walker steps if a deletion or compound event reaches past
`walker_pos`. That holding-open is the only reason this is more
than a `for ref_pos in 1..=genome_len` loop.

If you come from Python: think of the active set as a list of
generators (one per read), each yielding events tagged by reference
position. The walker is a coordinator that drains "all events at
`walker_pos`" from every generator before stepping forward, and
keeps a small dict of in-flight aggregations (the open records)
that it flushes when they're done.

## The six submodules

[`src/per_sample_caller/pileup/`](../../../src/per_sample_caller/pileup/)
is split into six files plus a `mod.rs` that re-exports the public
types. The split is deliberate — each file owns one concept and
there is a clear handoff between them. Read in this order:

| Module | Lines | Owns |
|---|---|---|
| [`mod.rs`](../../../src/per_sample_caller/pileup/mod.rs) | 263 | Public types: `PreparedRead`, `PileupRecord`, `AlleleObservation`, `FiveScalars`, `SlotId`, `RefBaseFetcher`, `WalkerConfig`. No logic. |
| [`walker.rs`](../../../src/per_sample_caller/pileup/walker.rs) | 628 | The main loop. Owns `WalkerState` and the per-step phase choreography. |
| [`active_set.rs`](../../../src/per_sample_caller/pileup/active_set.rs) | 317 | The set of reads currently in flight (admitted, not yet expired). One entry per read. |
| [`cigar_cursor.rs`](../../../src/per_sample_caller/pileup/cigar_cursor.rs) | 1301 | Per-read state: where each CIGAR op starts on the reference and on the read. Lazily emits events at any queried position. |
| [`open_record.rs`](../../../src/per_sample_caller/pileup/open_record.rs) | 854 | The records currently held open. Keyed by anchor position. Owns the fold logic that turns events into per-allele scalars. |
| [`slot_allocator.rs`](../../../src/per_sample_caller/pileup/slot_allocator.rs) | 595 | Phase-chain slot ids. Pairs mates onto the same slot, recycles slots when reads expire. |
| [`decompose.rs`](../../../src/per_sample_caller/pileup/decompose.rs) | 580 | **Test-only** since the lazy cursor landed. A reference oracle that decomposes a whole CIGAR up-front; used by parity tests against the cursor. |

The "lazy cursor" point in the table is a real surprise worth
calling out: the spec describes eager decomposition into a
`Vec<ReadEvent>` per read, but the production code went lazy —
events are computed on-demand at the position the walker queries.
`decompose.rs` survives only as a parity oracle.

## The walker loop, phase by phase

[`walker.rs:24`](../../../src/per_sample_caller/pileup/walker.rs#L24)
defines the entry point `run`. The interesting code is the
`while` loop at lines 42–89:

```text
while there is more work to do:
  if chromosome boundary:
    flush_chromosome   ← drain everything for current chrom

  while next read's alignment_start ≤ walker_pos:
    admit_read         ← active_set.rs:89

  process_position     ← walker.rs:200ish (folds events into records)
  expire_passed_reads  ← active_set.rs:134
  close_aged_records   ← open_record.rs:142 (drain via BTreeMap)
  advance              ← walker.rs (jumps to next interesting pos)
```

Each phase takes the state from the previous phase and produces a
new state for the next. The viewer's "Step" button advances one
phase at a time so you can watch each transition.

### Phase 1 — admit_read

For every pending read with `alignment_start ≤ walker_pos`:

1. Allocate a fresh `read_id` (monotonic `u32`).
2. Build a `CigarCursor` (precompute the `(ref_pos, read_pos)` pair
   at the start of every CIGAR op, plus a sentinel — see
   [`cigar_cursor.rs`](../../../src/per_sample_caller/pileup/cigar_cursor.rs)).
3. Allocate a `SlotId` from the slot allocator. If the read is the
   second mate of a pair already in flight, **reuse** the first
   mate's slot rather than allocating a new one
   ([`slot_allocator.rs`](../../../src/per_sample_caller/pileup/slot_allocator.rs)).
4. Push to the active set.

Why mate-pair slot sharing matters: at fold time the open-record
code uses the slot id to detect that two contributors at the same
position are the two halves of the same fragment, and applies
samtools' BQ math (sum on agreement, scale-down on disagreement)
instead of double-counting the molecule.

### Phase 2 — process_position

The algorithmic heart. For the current `walker_pos`:

1. **Gather contributors.** For each read in the active set, ask
   the cursor `events_at(walker_pos, &read)`. The cursor returns
   zero or more events anchored exactly at `walker_pos`:
   - A `Match` at every M-op base (one per base, but `events_at`
     returns at most one because at most one M base falls on any
     given reference position).
   - An `Insertion` event if the M→I transition's anchor is here.
   - A `Deletion` event if the M→D transition's anchor is here.

2. **Mate-overlap resolution.** Group contributors by `slot_id`.
   Any slot with ≥ 2 contributors is a mate pair overlapping at
   this position. Apply the deterministic resolution rule (see
   "Mate-pair overlap" section below).

3. **Column depth cap.** If contributor count exceeds
   `max_snp_column_depth` (or `max_indel_column_depth` if any
   contributor has an indel event here), truncate the contributor
   list. This is the only "approximation" in the walker — needed
   for pathological pileups (~PCR-duplicated regions).

4. **Fold.** For each surviving contributor, find or open the
   open-record at `walker_pos`, and fold the contributor's events
   into the right allele bucket.

The "find or open" call is where REF-span widening can happen.
If a contributor brings a `Deletion` event whose footprint
extends past the existing record's REF span, the record's REF
sequence is widened (more bases pulled from the FASTA, every
allele's `seq` extended in lock-step). Scenario 2 in the viewer
shows this exact moment.

### Phase 3 — expire_passed_reads

For every read whose `alignment_end < walker_pos`: drop it from
the active set, release its slot. The slot allocator may put the
slot on a "pending free" list (rather than immediately recycling
it) so that the lifecycle marks `expired_chains` get a chance to
land in the next emitted record before the slot is reused —
otherwise a downstream consumer would see `new[5]` for one
record and `new[5]` again for the next without an intervening
`expired[5]`.

### Phase 4 — close_aged_records

Iterate open records in coordinate order
([`BTreeMap::iter`](https://doc.rust-lang.org/std/collections/struct.BTreeMap.html#method.iter)
over `pos → record`). For each record where `pos + ref_span ≤
walker_pos`: finalise (compute final allele list with phase chain
slot lists), send through the `mpsc::SyncSender<PileupRecord>`,
remove from the map.

The closure rule is provable from the input invariant
(coordinate-sorted reads). The proof is short:

- Reads are sorted by `alignment_start`, so any future read has
  `alignment_start ≥ walker_pos + 1`.
- Every event's anchor sits at an M-op position, so future events
  have `anchor ≥ walker_pos + 1`.
- An event with anchor `S` can merge into a record at `pos` only
  if `S < pos + ref_span` (the record's REF span has to cover the
  event's anchor).
- Combining: a future merge requires `S ≥ walker_pos + 1` AND `S
  < pos + ref_span`. If `walker_pos ≥ pos + ref_span`, no `S`
  satisfies both → safe to close.

That proof is why the walker never needs a fixed look-back
window — the closure rule is per-record, sized to that record's
own footprint.

### Phase 5 — advance

If any active read still has events at or beyond `walker_pos +
1`, set `walker_pos += 1`. Otherwise jump to the next pending
read's `alignment_start` (skip uncovered span). This is what
makes "uncovered position produces no record" naturally fall
out: the walker simply isn't there.

## CIGAR cursor

[`cigar_cursor.rs`](../../../src/per_sample_caller/pileup/cigar_cursor.rs)
is the biggest module by line count. It deserves a closer look.

A `CigarCursor` holds:

- The CIGAR ops (immutable for the lifetime of the read).
- A precomputed `Vec<OpOffset>` — for each op, `(ref_pos at op
  start, read_pos at op start)`, plus a sentinel after the last
  op. This is what makes random-access `events_at(walker_pos)`
  cheap.
- A `mode: CursorMode` — `Linear` (≤16 ops, scan the offset
  table) or `BinarySearch` (>16 ops, binary search). Picked at
  construction.

Two query methods:

- **`events_at(walker_pos, read)`** — returns events whose
  anchor equals `walker_pos`. Used by `process_position`.
- **`events_overlapping(lo, hi, read)`** — returns events whose
  *footprint* intersects `[lo, hi)`. Used at fold time when the
  open-record needs to know "what does this read present as a
  haplotype string under this record's REF window?". A Match
  has a footprint of one base (its anchor). A Deletion has a
  footprint of `[anchor, anchor + deleted_len + 1)` (= the
  anchor base, plus the deleted bases that the deletion's REF
  side covers).

The overlap query is what makes compound REF widening possible:
once a record's REF span grows from 1 to 4, the fold re-queries
each contributor with `events_overlapping(pos, pos + 4)` and
rebuilds that contributor's haplotype string. Every contributor
gets folded exactly once per walker step (enforced by the
`folded_reads` map in `OpenPileupRecord`).

### Adaptor boundary

A subtle filter applied inside the cursor's M-event emit: a base
whose `ref_pos` lies past the mate-pair adaptor boundary is
silently dropped. The boundary is direction-aware:

- **Forward strand:** drop when `ref_pos ≥ adaptor_boundary`.
- **Reverse strand:** drop when `ref_pos ≤ adaptor_boundary`.

The boundary is supplied on `PreparedRead.adaptor_boundary` (or
`None` if it can't be reliably computed — single-end, mate
unmapped, geometry inconsistent). The cursor consults it inline
([`cigar_cursor.rs:69`](../../../src/per_sample_caller/pileup/cigar_cursor.rs#L69))
rather than relying on the read having been pre-trimmed
upstream. This was the G1 finding from
[`ia/reviews/pileup_gatk_comparison_2026-05-08.md`](../../reviews/pileup_gatk_comparison_2026-05-08.md).

### Read-N handling

A read base of `N` is silently skipped on the M emit path
([`cigar_cursor.rs:247`](../../../src/per_sample_caller/pileup/cigar_cursor.rs#L247)).
This is the F5 commit. Reference Ns flow through verbatim; only
read Ns get filtered.

## Open-record store

[`open_record.rs`](../../../src/per_sample_caller/pileup/open_record.rs)
is structured around a `BTreeMap<u32, OpenPileupRecord>` keyed by
anchor position. The choice of `BTreeMap` (rather than `HashMap`)
matters for two reasons:

- Closure scans positions in order via `iter()`, stopping at the
  first record where `pos + ref_span > walker_pos`. With a
  `HashMap` we'd have to either scan all entries or maintain a
  separate sorted index.
- "Find any record overlapping `[lo, hi)`" is a range query, also
  natural on `BTreeMap`.

Each `OpenPileupRecord` holds:

- `pos: u32` — anchor position.
- `alleles: Vec<OpenAllele>` — `alleles[0]` is REF (invariant).
  Every allele's `seq` has the same length (= record's `ref_span`).
- `folded_reads: AHashMap<u32, FoldedReadState>` — for each
  read_id that has contributed to this record, what allele
  index they were assigned and what scalars they added. Used to
  *unfold* a contribution before re-folding it after a REF
  widen.

The widening protocol
([`open_record.rs:281`](../../../src/per_sample_caller/pileup/open_record.rs#L281)):

1. Fetch the new REF bases from the FASTA.
2. Append them to every allele's `seq` (REF gets the literal
   bases, alts get the same bases unless the alt's event
   produces something different in the new region).
3. For every previously folded read, re-evaluate which allele
   bucket the read presents under the wider REF span. The
   read's prior contribution is unfolded from the old bucket;
   the new contribution is folded into the (possibly new)
   bucket.

This is the most work-per-step in the whole walker, but it's
amortised: each widen step adds bases at a per-base cost, and
the number of widens per record is bounded by the longest event
that touches the record (which is bounded by `MAX_RECORD_SPAN
= 5000`).

## Slot allocator

[`slot_allocator.rs`](../../../src/per_sample_caller/pileup/slot_allocator.rs)
mints `SlotId` values (`u16`, 12 bits used, hard cap
`MAX_ACTIVE_SLOTS = 4096`). The id is the durable handle for a
"phase chain" — alleles sharing the same slot id across positions
are linked: they are observed on the same haplotype molecule.

The interesting state:

- `next_fresh: SlotId` — next id to mint if the free list is
  empty.
- `free: Vec<SlotId>` — recycled ids waiting to be reused.
- `pending_free: Vec<SlotId>` — ids whose refcount just hit zero
  but whose `expired[id]` mark hasn't yet been emitted in a
  closed record. Migrated to `free` only after lifecycle marks
  drain.
- `pending_mates: HashMap<Arc<str>, PendingMate>` — qname →
  first mate's slot/read_id, indexed for O(1) second-mate
  lookup.
- `refcounts: HashMap<SlotId, u8>` — how many active reads
  reference each slot. 1 for solo, 2 for an in-flight pair.

The mate-pair refcount trick: when the first mate arrives, the
slot is allocated with refcount **2** (anticipating the second
mate). If the second mate never shows up within
`MATE_LOOKUP_WINDOW = 10_000` bp, the pending entry is evicted
([`slot_allocator.rs:309`](../../../src/per_sample_caller/pileup/slot_allocator.rs#L309))
and the refcount drops back to 1, so the slot releases when the
first mate exits.

## Mate-pair overlap

The trickiest small piece of the walker.
[`walker.rs:394`](../../../src/per_sample_caller/pileup/walker.rs#L394)
(`resolve_mate_overlap_at_pos`).

For each `(slot_id, contributors[])` group with ≥ 2 entries:

- **If either side has an indel event at `walker_pos`:** drop
  the loser entirely (one-observation collapse — we cannot sum
  evidence across mates that disagree on indel placement).
  Tiebreak order: BQ (higher wins) → `is_first_mate` (true
  wins) → smaller `alignment_start` wins.
- **Match-only overlap, bases agree:** keep one contributor with
  `BQ = min(BQ₁ + BQ₂, 200)`; the other contributes nothing
  (its BQ is zeroed at `walker_pos`).
- **Match-only overlap, bases disagree:** keep the higher-BQ
  contributor with `BQ ← 0.8 × BQ`; the other contributes
  nothing.

The "zeroing" is implemented as flags on the contributor
(`bq_zero_in_window`, `bq_override_at_walker_pos`) that the
fold honors when computing per-position BQ. The fold doesn't
*remove* the contributor — it just makes that contributor's BQ
contribution at `walker_pos` zero — so phase-chain attribution
across positions is preserved.

## Five per-allele scalars (`FiveScalars`)

[`mod.rs:170`](../../../src/per_sample_caller/pileup/mod.rs#L170).
The fold maintains, per allele per record, exactly these:

| Field | Updated by |
|---|---|
| `num_obs` | +1 per supporting read |
| `q_sum` (= Σ `max(ln_BQ, ln_MQ)`) | += `max(ln_err_bq, ln_err_mq)` per supporting read |
| `fwd` | +1 if read is forward strand |
| `placed_left` | +1 if read's `alignment_start < record.pos` |
| `placed_start` | +1 if read's `alignment_start == record.pos` |

These are sufficient to reconstruct freebayes' per-genotype
likelihood downstream — that's the durability claim of the
overall pipeline (see
[`ia/specs/calling_pipeline_architecture.md` §"The five
per-allele scalars"](../../specs/calling_pipeline_architecture.md)).

## How to read the source

Recommended reading order, given this map:

1. [`mod.rs`](../../../src/per_sample_caller/pileup/mod.rs) — types only, ~10 minutes.
2. [`walker.rs::run`](../../../src/per_sample_caller/pileup/walker.rs#L24) — the main loop body. Stop at `process_position` calls.
3. [`active_set.rs`](../../../src/per_sample_caller/pileup/active_set.rs) — small file, easy admit/expire.
4. [`slot_allocator.rs`](../../../src/per_sample_caller/pileup/slot_allocator.rs) — read in pieces; the freelist + pending_mates story is the interesting bit.
5. [`cigar_cursor.rs`](../../../src/per_sample_caller/pileup/cigar_cursor.rs) — start with `events_at`, then `events_overlapping`. Skip `BinarySearch` mode on first read.
6. [`open_record.rs::process_position`](../../../src/per_sample_caller/pileup/open_record.rs) — the fold. Read after you understand the cursor.
7. [`walker.rs::resolve_mate_overlap_at_pos`](../../../src/per_sample_caller/pileup/walker.rs#L394) — apply once everything else is in place.

Open [`index.html`](index.html) alongside this and step through
the three scenarios. The viewer names the function being
"executed" at each step; you can pull up the matching Rust file
and follow along.

## Known scope limits of the viewer

The animation is faithful to the algorithm but simplified for
clarity:

- **No BAQ.** Read BQ is the raw value — the production walker
  consumes BAQ-adjusted BQ from upstream. The arithmetic is
  identical; it's just that the input number is different.
- **`q_sum` shown as `Σ BQ`, not `Σ max(ln_BQ, ln_MQ)`.** The
  production scalar combines BQ with MQ in log space; the
  viewer just displays the BQ side because the log arithmetic
  is hard to read at a glance. The combination point is the
  same.
- **Mate-overlap rule shown as "agree → keep one, sum BQ;
  disagree → keep higher, scale 0.8".** Indel-overlap collapse
  is implemented in the production code but not exercised by
  any of the three scenarios.
- **Column depth cap not exercised.** Default `8000` /  `250` —
  none of the scenarios approach it.
- **Phase chain ids shown but not animated across scenarios.**
  Each scenario completes within one or two reads, so the
  chain-id story is only visible inside `mate_overlap`.

These limits are properties of the viewer's pedagogical scope,
not of the production walker.
