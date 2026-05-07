# Pileup walker — lazy CIGAR decomposition

Implementation plan for finding `S2` of
[`ia/reviews/pileup_samtools_comparison_2026-05-07.md`](../reviews/pileup_samtools_comparison_2026-05-07.md):
replace the walker's eager per-read CIGAR decomposition with a
lazy, cursor-driven design so memory and CPU scale with reads, not
with read length × span.

## Motivation

The current walker eagerly turns every read into a
`Vec<ReadEvent>` on admission and stores it on the `ActiveRead`
entry. For Match-heavy reads this materialises one `ReadEvent` per
covered reference position — fine at 150 bp Illumina, fatal for
20–50 kb long reads.

Two concrete hot paths today:

- **Per-read storage.** `ActiveRead.events: Vec<ReadEvent>`
  ([active_set.rs:17–30](../../src/per_sample_caller/pileup/active_set.rs#L17-L30))
  is O(span) per active read.
- **Per-step clones.** Inside `WalkerState::process_position`
  ([walker.rs:201–239](../../src/per_sample_caller/pileup/walker.rs#L201-L239))
  the walker builds one `ReadContribution` per active read per
  walker step, and each one *clones* the read's full event list
  into `full_window_events`. A read of span `L` lives across `L`
  walker steps and is cloned once per step → **O(L²) per read**.
  At 150 bp this is 22 500 clones per read (irrelevant); at 20 kb
  it is 4×10⁸ (fatal).

Replacing both with a CIGAR-cursor design drops per-read storage
to O(num_cigar_ops) and per-step work for clean Match-only reads
to O(1).

### Empirical confirmation

Measured by [`benches/pileup_walker_scaling.rs`](../../benches/pileup_walker_scaling.rs)
on 2026-05-07 (release build, 50 kb window, 30× coverage,
pure-Match reads):

| read length | median time | per-position | scaling vs L=150 |
|---|---|---|---|
|   150 | 1.42 s |    28 µs | 1.0× |
|   500 | 3.02 s |    60 µs | 2.1× |
|  1500 | 8.70 s |   174 µs | 6.1× |
|  5000 | 70.9 s |  1418 µs | 50×  |

Read length grows 33× from L=150 to L=5000; walker time grows
50×. The slope is **super-linear** between L=1500 and L=5000
(8.2× time for 3.3× length), consistent with allocator/cache
pressure on top of the predicted clone work. L=5000 is the
current `MAX_RECORD_SPAN` ceiling; long-read support cannot land
without addressing this.

The benchmark is left in the tree as a regression guard: commit 2
of the migration plan must keep L=150 within ±10% of baseline, and
should produce a step-change improvement at L=1500 and L=5000.

## Scope

In:

- A new `CigarCursor` (or equivalent) type, wrapping a borrow of the
  read's CIGAR plus an offset table, that emits events on demand.
- Replacing `ActiveRead.events` and `event_cursor` with the cursor.
- Replacing `ReadContribution.full_window_events: Vec<ReadEvent>`
  with a window-query API the open-record fold path calls into.
- Removing the per-step clone of the entire event list.
- A long-read benchmark (synthetic, in `benches/`) demonstrating
  the linear-in-`L` runtime and memory.

Out:

- The public `ReadEvent` enum stays as-is — Stage 2 and tests
  consume it indirectly through `AlleleObservation` and the open
  record fold; no public type churn.
- BAQ, mate-overlap, slot allocator semantics — untouched. This is
  a pure refactor with byte-identical pileup output.
- Long-read *input* support end to end (CRAM decoder limits, BAQ
  scaling, MAX_RECORD_SPAN tuning). This plan only removes the
  pileup walker as the bottleneck; downstream stages are out of
  scope.

## Current data flow (the pieces this plan touches)

1. `ActiveSet::admit` calls `decompose(&read)` and stores the
   resulting `Vec<ReadEvent>` plus a cursor index on the
   `ActiveRead`.
2. At each walker step, `WalkerState::process_position`
   ([walker.rs:201–239](../../src/per_sample_caller/pileup/walker.rs#L201-L239)):
   - Pops events whose anchor matches `walker_pos` from
     `active.events[event_cursor..]`.
   - Builds a `ReadContribution { events_at_pos, full_window_events:
     active.events.clone(), ... }` per contributing read.
3. `open_record::process_position`
   ([open_record.rs:486–495](../../src/per_sample_caller/pileup/open_record.rs#L486-L495))
   filters `full_window_events` to events overlapping the open
   record's footprint, then `apply_events_to_ref` consumes the
   filtered slice.

The eager `decompose` in
[decompose.rs](../../src/per_sample_caller/pileup/decompose.rs) is
the function being replaced; its tests stay (re-pointed at the
cursor).

## Target design

### `CigarCursor` — internal type, owned by `ActiveRead`

```rust
/// Lazy walker over one read's CIGAR. Per-active-read state
/// replacing the eager `Vec<ReadEvent>` and `event_cursor`.
///
/// Holds an offset table (one entry per CIGAR op) so a query for
/// "the event whose anchor is `walker_pos`" or "events in window
/// `[lo, hi)`" is O(log num_ops + num_emitted_events) rather than
/// O(span).
pub(super) struct CigarCursor {
    /// Reference to the read's CIGAR. Lifetime tied to the
    /// `ActiveRead` that owns the cursor.
    cigar: Arc<[CigarOp]>,
    /// Per-op cumulative ref/read offsets:
    ///   offsets[i] = (ref_pos_at_op_start, read_pos_at_op_start)
    /// Length is `cigar.len() + 1`; the last entry is end-of-read.
    /// Computed once on cursor construction; ~50 entries even for
    /// long reads.
    offsets: Vec<OpOffset>,
    /// Highest ref position whose events have been *emitted* via
    /// `next_at(walker_pos)`. Replaces `event_cursor`. Used to
    /// guarantee each event is folded exactly once per pair
    /// `(record, read)` (the `folded_reads` invariant in
    /// open_record.rs is unchanged).
    walker_high_water: u32,
}

#[derive(Debug, Clone, Copy)]
struct OpOffset {
    ref_pos: u32,    // 1-based ref position at op start
    read_pos: u32,   // 0-based read offset at op start
}

impl CigarCursor {
    pub(super) fn new(cigar: Arc<[CigarOp]>, alignment_start: u32) -> Self;

    /// Events at exactly `walker_pos`. Match events derived on the
    /// fly from (op, offset). Indels are anchored at the previous
    /// op's last ref position per the existing first/last-op drop
    /// rule. Advances `walker_high_water` to `walker_pos + 1` on
    /// success.
    pub(super) fn next_at(
        &mut self,
        walker_pos: u32,
        seq: &[u8],
        bq_baq: &[u8],
    ) -> SmallVec<[ReadEvent; 2]>;

    /// Events whose anchor falls in `[lo, hi)` — used by the
    /// open-record fold to compute `apply_events_to_ref`. Stateless
    /// (does not advance the persistent cursor): walks the CIGAR
    /// from the op containing `lo` forward to `hi`.
    pub(super) fn events_in_window(
        &self,
        lo: u32,
        hi: u32,
        seq: &[u8],
        bq_baq: &[u8],
    ) -> Vec<ReadEvent>;
}
```

`SmallVec<[ReadEvent; 2]>` is the inline-2 specialisation: a
position can carry at most one Match plus one indel anchored
there, so 2 is the spec maximum. Avoids a heap allocation on the
common path.

The `seq`/`bq_baq` slices stay on `PreparedRead` and are passed in
per call rather than copied into the cursor — the cursor never
needs ownership.

### `ActiveRead` — slimmer

```rust
pub struct ActiveRead {
    pub read_id: u32,
    pub read: PreparedRead,
    pub cursor: CigarCursor,         // replaces `events` + `event_cursor`
    pub chain_slot_id: SlotId,
    pub mate_read_id: Option<u32>,
}
```

No semantic change at the active-set boundary. `admit_read` gains
a CIGAR-arc clone (shared with the cursor) and otherwise looks the
same.

### `ReadContribution` — no full window, no lifetime

```rust
pub struct ReadContribution {
    pub read_id: u32,
    pub chain_slot_id: SlotId,
    pub events_at_pos: SmallVec<[ReadEvent; 2]>,
    pub bq_baq_at_walker_pos: u8,
    pub mq_log_err: f64,
    pub is_reverse_strand: bool,
    pub alignment_start: u32,
    pub is_first_mate: bool,
    /// Set by `resolve_mate_overlap_at_pos` for the overlap loser
    /// in the match-only regime. Consumers of `events_in_window`
    /// zero out the BQ on each emitted event when this is `true`.
    /// Replaces the in-place mutation of the cloned event vector
    /// the eager design did
    /// ([walker.rs:481–491](../../src/per_sample_caller/pileup/walker.rs#L481-L491)).
    pub bq_zero_in_window: bool,
}
```

The contribution carries no borrow. When the open-record fold
needs window events for a contributor, it looks the read up by
`read_id` against an `&ActiveSet` passed alongside the
contributor list:

```rust
pub fn process_position(
    open: &mut OpenPileupRecordTable,
    walker_pos: u32,
    chrom_id: u32,
    contributors: &[ReadContribution],
    active: &ActiveSet,                   // new arg
    fasta: &dyn RefBaseFetcher,
) -> Result<ProcessOutcome, WalkerError> {
    ...
    let active_read = active
        .get_by_read_id(contrib.read_id)
        .expect("contributor's read_id must be live");
    let mut window_events = active_read.cursor.events_in_window(
        rec_pos, rec_end, &active_read.read.seq, &active_read.read.bq_baq,
    );
    if contrib.bq_zero_in_window {
        for ev in &mut window_events { zero_event_bq(ev); }
    }
    ...
}
```

`ActiveSet::get_by_read_id` is already `pub(super)` and goes
through the existing `by_read_id: AHashMap<u32, usize>`
([active_set.rs:40](../../src/per_sample_caller/pileup/active_set.rs#L40)) —
constant-time, no new infrastructure.

**Why not a lifetime?** An earlier draft of this plan put
`active: &'a ActiveRead` directly on the contribution. The
runtime is identical (a Rust lifetime is a compile-time-only
proof; `&'a T` is the same 8 bytes as a raw pointer), and
**memory cost is the same** for both designs:

| design | per-step contribution memory at L=5000, depth 30 |
|---|---|
| eager today | ~6 MB (30 reads × 5000 events × ~40 B, cloned) |
| lazy + `&'a ActiveRead` | ~2 KB (30 × ~64 B) |
| lazy + `read_id` lookup | ~2 KB (30 × ~64 B) |

The reason to prefer the `read_id` form is **complexity, not
memory**: the `<'a>` would ripple into every helper that touches
a `ReadContribution`, and the mate-overlap path mutates the
contributor list (`swap_remove`) while a borrow into `&self.active`
would be live, which the borrow checker will fight. The id form
sidesteps both.

The runtime cost is one hashmap lookup per `events_in_window`
query: ~30 contributors × ~5 records × 50 000 walker steps ≈ 7.5 M
lookups in the L=5000 benchmark workload, comfortably under the
margin of error of the 70.9 s baseline.

## Migration plan — three commits, each shippable

### Commit 1 — introduce `CigarCursor`, prove parity

- Add the `CigarCursor` type (private to the `pileup` module).
- Implement `next_at` and `events_in_window`.
- Keep `decompose` and the eager event list. **Do not change the
  walker yet.**
- Add a property test: for every test read in the existing pileup
  test corpus
  ([pileup/tests.rs](../../src/per_sample_caller/pileup/tests.rs)
  and [pileup/decompose.rs](../../src/per_sample_caller/pileup/decompose.rs)),
  assert that the cursor-emitted event sequence matches
  `decompose(read)` byte-for-byte. Run on every pileup test as a
  parallel oracle.
- No behavioural change. The diff is purely additive.

### Commit 2 — switch the walker to the cursor

- Replace `ActiveRead.events`/`event_cursor` with `cursor`.
- Drop `ReadContribution.full_window_events` and add the
  `bq_zero_in_window` flag. No lifetime parameter on the struct.
- Add the `active: &ActiveSet` argument to
  `open_record::process_position`; switch the fold to call
  `events_in_window` via `active.get_by_read_id(contrib.read_id)`
  instead of filtering a cloned `Vec`.
- Update `resolve_mate_overlap_at_pos`: where it used to zero
  per-event BQ on the loser's cloned event vector, set the
  contributor's `bq_zero_in_window = true` instead.
- Replace `decompose` callers; keep the function itself if it
  serves any tests, otherwise delete.
- Run the full test suite and the existing pileup property tests.
  Output must remain byte-identical (the cursor was proven
  equivalent in commit 1).

### Commit 3 — re-run the existing bench, record the win

The motivating benchmark
[`benches/pileup_walker_scaling.rs`](../../benches/pileup_walker_scaling.rs)
already exists with the eager-decompose baseline numbers shown
above. Commit 3 re-runs it on the cursor-based walker and records
the result.

- Targets at the same `(span = 50 kb, coverage = 30)` workload:
  - **L = 150 (Illumina):** ±10% of the 1.42 s baseline. This is
    the no-regression bar — any slowdown here would be a blocker
    for shipping commit 2 even if long-read numbers improve.
  - **L = 1500:** ≥3× faster than the 8.70 s baseline.
  - **L = 5000:** ≥10× faster than the 70.9 s baseline; per-position
    time should be near flat across L (the clone is gone).
- Also add an indel-heavy variant (e.g. one 10 bp deletion every
  500 bp) so the bench covers more than the pure-Match path.
- Record both runs in
  [`ia/reports/implementations/`](../reports/implementations/)
  with a short note linking back to this plan.

## Test plan

- **Property parity (commit 1).** Cursor and `decompose` emit
  identical event streams on every existing test read.
- **Existing pileup test suite** ([tests.rs](../../src/per_sample_caller/pileup/tests.rs))
  runs unchanged after commit 2.
- **Cursor unit tests** (alongside the new module): Match-only,
  indel-in-middle, indel-at-edges (drop), N-skip, soft-clip,
  hard-clip, padding, consecutive indels, BQ proxy windows. Mirror
  the suite in
  [decompose.rs::tests](../../src/per_sample_caller/pileup/decompose.rs).
- **Edge case for the window query.** A wide deletion record
  whose footprint extends well past the walker_pos: verify
  `events_in_window(rec_pos, rec_end)` returns events from the
  *future* of the cursor's persistent state without advancing it.
- **Mate-overlap loser BQ zeroing.** Reproduce
  `mate_overlap_winner_keeps_bq_loser_zeroed` from existing tests
  with the lazy path; output is byte-identical.

## Risks and mitigations

- **Stale `read_id` in a contribution.** The contribution carries
  only `read_id`; `process_position` looks the read up via
  `ActiveSet::get_by_read_id`. A lookup miss would be a logic bug
  (a contributor outliving its source read) and is asserted, not
  silently handled. The current code has the same class of
  invariant on `event_cursor` ranges, with the same handling
  shape — no new failure mode.

- **`events_in_window` on a wide record.** The function walks the
  CIGAR forward from the op containing `lo`. For a deletion
  spanning thousands of bases the walk is bounded by `MAX_RECORD_SPAN`
  ([mod.rs:38](../../src/per_sample_caller/pileup/mod.rs#L38))
  and by op count (typically <50). No regression vs. today, which
  filters the same cloned vec.

- **Match-event materialisation.** A naive `events_in_window`
  emits one `ReadEvent::Match` per ref position inside an M op,
  which on a wide record reproduces today's allocation. Acceptable
  — the win is in the per-read state and per-step clone, not in
  the per-record fold (which is unavoidable for `apply_events_to_ref`).
  If a future profile says otherwise, the fold itself can be
  rewritten to walk ops directly without materialising Match
  events; out of scope here.

- **Spec note required.** Add a paragraph to
  [`ia/specs/pileup_walker.md`](../specs/pileup_walker.md)
  §"Read decomposition" describing the cursor as the canonical
  representation. Cheap; do this in commit 2.

## Rollback

Each commit is independently revertable. Commit 2 is the only one
with behavioural surface area; if regression is found post-merge,
revert it in isolation while leaving the cursor type in place for
a second attempt.

## Out-of-scope follow-ups

- Walking ops directly inside `apply_events_to_ref` to avoid
  materialising Match events for wide records — only worth doing if
  benchmarks point at it.
- Long-read CRAM input (decoder, BAQ tuning, `MAX_RECORD_SPAN`
  raise). Tracked separately if/when long-read support is decided.
- **Binary-search the offset table.** The `CigarCursor.offsets`
  vector is already sorted by `ref_pos`, so it doubles as an
  index. `events_at(walker_pos)` could `partition_point` for the
  op containing `walker_pos` and inspect at most two ops (the M
  op containing `walker_pos`, and the op possibly anchoring an
  indel at `walker_pos + 1`); `events_overlapping(lo, hi)` could
  start from the first op with `ref_pos >= lo - max_deletion_len`
  (a new field computed once at `new()` to bound left-side
  lookback for deletions whose footprint extends rightward).
  Yields O(log n) lookup vs. the current O(n) op walk. Skipped
  for now because:
  - For Illumina reads (1–5 CIGAR ops) the linear walk plus
    early break already costs essentially nothing, and the
    binary-search overhead would be a wash or slight regression.
  - The benefit shows up at 50+ ops per read (PacBio HiFi and
    longer), which we don't have a benchmark for yet.
  - Adding `max_deletion_len` ships state we'd never measure
    until that benchmark exists.

  Revisit alongside the indel-heavy bench variant called out in
  commit 3, or when long-read CRAM input lands.
