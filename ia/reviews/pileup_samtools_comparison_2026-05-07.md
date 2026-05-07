# Algorithm Review: pileup walker vs. samtools

**Date:** 2026-05-07
**Reviewer:** Claude (algorithm-comparison study)
**Module reviewed:** `src/per_sample_caller/pileup`
**Reference codebase:** `samtools/` (in-tree copy of samtools, not htslib)
**Status:** Advisory — no defects; a small backlog of optional improvements

---

## 1. Scope

This is **not** a defect-finding review. It is a comparative study of
the samtools pileup machinery against our Stage 1 walker, looking for
techniques and corner-case handling we should consider importing.

In-scope project files:

- [src/per_sample_caller/pileup/mod.rs](../../src/per_sample_caller/pileup/mod.rs)
- [src/per_sample_caller/pileup/walker.rs](../../src/per_sample_caller/pileup/walker.rs)
- [src/per_sample_caller/pileup/active_set.rs](../../src/per_sample_caller/pileup/active_set.rs)
- [src/per_sample_caller/pileup/decompose.rs](../../src/per_sample_caller/pileup/decompose.rs)
- [src/per_sample_caller/pileup/slot_allocator.rs](../../src/per_sample_caller/pileup/slot_allocator.rs)
- [src/per_sample_caller/pileup/open_record.rs](../../src/per_sample_caller/pileup/open_record.rs)

Samtools files studied:

- `samtools/bam_plbuf.{c,h}` — thin callback wrapper over htslib's pileup.
- `samtools/bam_lpileup.{c,h}` — "level pileup" used by `tview` to assign
  display rows; structurally analogous to our slot allocator.
- `samtools/bam_plcmd.c` — the `samtools mpileup` driver: filters, BAQ,
  multi-sample coordination, mate-overlap smashing.
- `samtools/consensus_pileup.{c,h}` — alternative pull-based pileup
  used by `samtools consensus`.

Out of scope:

- htslib's core pileup engine (`bam_plp_*`, `overlap_push`,
  `bam_prob_realn`) is not in this repo; behaviour is inferred from
  call sites only.
- Stage 2 encoder, BAQ stage, upstream `cram_input` filters.

## 2. Summary verdict

**No correctness divergences found** between our walker and samtools
where the spec is shared. In fact, on two corner cases we *beat*
samtools (see §3 "Where we already match or beat samtools").

A short list of optional improvements has been distilled, ranked by
priority and effort. They are independent and can be picked off one
at a time.

## 3. Where we already match or beat samtools

These are noted up-front so we don't accidentally regress them while
acting on the findings below.

- **`A1` — Decomposition robustness at read edges.**
  `samtools/consensus_pileup.c:60–65` documents a known bug where
  insertions at the start of a sequence get silently merged with
  soft-clips; Gap5 worked around it by prepending a synthetic `1D` to
  every CIGAR. We handle this cleanly in
  [decompose.rs:94–98](../../src/per_sample_caller/pileup/decompose.rs#L94-L98)
  and [decompose.rs:111–114](../../src/per_sample_caller/pileup/decompose.rs#L111-L114)
  by dropping any indel that is the first/last CIGAR op or whose
  anchor would fall before reference position 1. Regression-tested
  by `first_op_indel_is_dropped_but_rest_of_read_keeps_events` and
  `soft_clip_then_indel_at_alignment_start_one_drops_indel` in
  [decompose.rs](../../src/per_sample_caller/pileup/decompose.rs).

- **`A2` — Mate-overlap evaluated per-position, not per-read.**
  htslib resolves overlap inside `overlap_push`, at admission/push
  time. We resolve at fold time at each `walker_pos` in
  [walker.rs:403–500](../../src/per_sample_caller/pileup/walker.rs#L403-L500).
  Our placement is strictly richer: the mate that loses at one
  position can win at another as BQ varies along the overlap region.
  Don't pull this forward.

- **`A3` — Slot-reuse cooling-off period.**
  `bam_lpileup.c:34` defines `TV_GAP=2`: a freed display level cools
  for two positions before reuse, purely for visualisation
  stability. Our [slot_allocator.rs `pending_free` / `free`
  split](../../src/per_sample_caller/pileup/slot_allocator.rs#L46-L52)
  is the same pattern with a different motivation (preventing
  `expired[S]` and `new[S]` from colliding in one record's
  lifecycle marks). Same idea, deliberately preserved.

## 4. Findings — backlog to act on

Each finding has a stable id (`Si` for "samtools-inspired") so we can
reference them in commits and a one-by-one work plan.

Priority key:
- **High** — clear user-visible benefit, low risk.
- **Medium** — usability or hygiene; nice to have.
- **Low** — only matters under conditions we don't currently target
  (e.g. long reads), or speculative.

### `S1` — Soft warning at slot-allocator high-water mark

- **Priority:** Medium
- **Effort:** Small (one log call, one threshold constant, one test)
- **Status:** Open

**Observation.** `bam_plcmd.c:1234` warns when `max_depth * nfn` is
projected above 1M, before the depth cap is hit. Our slot allocator
errors hard once `MAX_ACTIVE_SLOTS = 4096` is exceeded
([slot_allocator.rs:25](../../src/per_sample_caller/pileup/slot_allocator.rs#L25))
without any earlier signal.

**Proposal.** Emit a one-shot `tracing::warn!` (or equivalent) the
first time `slot_high_water` crosses, say, 75% of `MAX_ACTIVE_SLOTS`
during a run, naming the contig and position so the user can
investigate. Idempotent within a run.

**Rationale.** Lets users diagnose pathological-coverage regions
*before* the run dies, instead of greeting them with an
`OutOfSlots` error and no context.

**Risk.** Effectively none — it's a soft signal that does not change
behaviour on the happy path.

### `S2` — Document and (later) consider lazy CIGAR walking for long reads

- **Priority:** Low today, becomes High if/when long-read support enters scope
- **Effort:** Documentation now; substantial refactor later
- **Status:** Open

**Observation.** `samtools/consensus_pileup.c` (`get_next_base`) walks
CIGAR ops lazily via a single `(cigar_ind, cigar_op, cigar_len)` cursor
inside the per-read state. Our [decompose.rs](../../src/per_sample_caller/pileup/decompose.rs)
eagerly emits one `ReadEvent` per reference position into a
`Vec<ReadEvent>` on admission, stored on each `ActiveRead`
([active_set.rs:17–30](../../src/per_sample_caller/pileup/active_set.rs#L17-L30)).

For Illumina reads (≤300 bp) the eager approach is fine and simpler.
For long reads (PacBio, ONT — tens of kb), each `ActiveRead` would
carry tens of thousands of events; combined with a deep active set,
this is the dominant memory cost. Our `MAX_RECORD_SPAN = 5000` filter
([mod.rs:38](../../src/per_sample_caller/pileup/mod.rs#L38)) bounds the
issue but only because the cap is currently tight.

**Proposal.**

1. **Now (low cost):** add a one-paragraph note to
   [ia/specs/pileup_walker.md](../specs/pileup_walker.md) §"Read
   decomposition" stating that eager decomposition is intentional
   and tied to the short-read assumption; long-read support would
   require switching to a CIGAR cursor.
2. **Later (only if long-read support is on the roadmap):** convert
   `ActiveRead.events: Vec<ReadEvent>` to a stateful CIGAR cursor
   that emits the next event on demand from
   `walker.process_position`. Active-set memory drops from
   O(reads × span) to O(reads).

**Rationale.** Documenting the assumption is cheap insurance against
a silent footgun the next time we widen scope.

**Risk.** Stage 1 (doc only) — none. Stage 2 — touches the hottest
path; would need a dedicated benchmark and a careful diff.

### `S3` — Defensive supplementary/secondary alignment guard at admission

- **Priority:** Low
- **Effort:** Small (requires `is_supplementary`/`is_secondary` to
  flow into `PreparedRead`, plus a `debug_assert!` at admission)
- **Status:** Open

**Observation.** Both samtools and we rely on upstream filters to
exclude supplementary and secondary alignments from the pileup input.
Neither has an in-engine guard. If an upstream filter regression let
one through, we'd silently over-count (the read would be admitted,
take a slot, and contribute to records).

**Proposal.** Extend `PreparedRead`
([mod.rs:69–98](../../src/per_sample_caller/pileup/mod.rs#L69-L98)) with
two booleans (or a single `is_primary_alignment`) and add a
`debug_assert!` at the top of `ActiveSet::admit`
([active_set.rs:96](../../src/per_sample_caller/pileup/active_set.rs#L96)).
Production builds pay nothing; tests trip loudly on a regression.

**Rationale.** Defence-in-depth against an upstream filter bug that
would otherwise be invisible until variant calls disagreed with
ground truth.

**Risk.** Touches `PreparedRead`'s public shape — every test
constructor that hand-builds reads needs the new fields. That's the
main cost; the guard itself is trivial.

### `S4` — Verify (don't necessarily add) a reference-fetch LRU cache

- **Priority:** Medium (depends on benchmark)
- **Effort:** Small to investigate; small to add if needed
- **Status:** Open — investigation first

**Observation.** `bam_plcmd.c:310–357` keeps a 3-slot LRU of contig
sequences to amortise FASTA I/O across positions. Our walker calls
`RefBaseFetcher::fetch(chrom_id, start, length)` per record-open
inside `process_position`
([walker.rs:254–260](../../src/per_sample_caller/pileup/walker.rs#L254-L260)).

If the production wrapper around `noodles_fasta::Repository` already
keeps the active contig in memory, no action is needed. If it
re-reads or re-decompresses on each call, we are paying I/O per
record-open, which on dense regions is one fetch per emitted record.

**Proposal.**

1. Read the production `RefBaseFetcher` impl and confirm whether it
   caches.
2. If it doesn't, wrap it in a 1- or 2-slot contig LRU at the
   walker's call site (or inside the impl, whichever is cleaner).
3. If it does, close this finding with a note in the report.

**Rationale.** Cheap to verify; potentially a real per-sample
runtime win.

**Risk.** None for the investigation. If a cache is added, scope it
narrowly so we don't hide ref-fetch errors behind stale data.

### `S5` — Per-allele observation cap (deferred)

- **Priority:** Low
- **Effort:** Medium
- **Status:** Open — wait for evidence

**Observation.** samtools' mpileup uses a separate, lower cap for
indel columns (`MPLP_MAX_INDEL_DEPTH = 250` vs `MPLP_MAX_DEPTH = 8000`)
because homopolymer-context indels are pathologically deep in real
data. We have a uniform per-walker cap (active slot count) and no
per-allele bound.

**Proposal.** Hold off. Add the cap only if we see a Stage 2 failure
or a runtime spike caused by a single allele's `num_obs` blowing up.
Track this finding here so we don't forget the precedent exists.

**Rationale.** Premature optimisation otherwise; Stage 2 should be
robust to high allele depth in principle.

**Risk.** Adding it now would mean choosing a threshold without
data, which is the worst time to choose one.

## 5. Findings explicitly *rejected*

For the record, so we don't relitigate them:

- **R1 — Pull mate-overlap detection forward to admission.** Rejected
  in favour of the per-position fold (see `A2`).
- **R2 — Switch to a per-position callback model.** Rejected — our
  per-record closure rule is the foundation of `PileupRecord` /
  `AlleleObservation` / phase chains and is committed to in the
  spec.
- **R3 — Recompute BAQ inside the walker.** Rejected — the spec
  ([mod.rs:63](../../src/per_sample_caller/pileup/mod.rs#L63),
  [ia/specs/pileup_walker.md](../specs/pileup_walker.md)) commits to
  the walker treating `bq_baq` as opaque; BAQ is upstream.

## 6. Suggested execution order

If acted on in order, each finding is independent and can be its own
small commit:

1. `S1` — soft high-water warning (one-line user-visible win).
2. `S3` — supplementary/secondary guard (cheap defence-in-depth).
3. `S4` — investigate ref-fetch caching; add if needed.
4. `S2` — long-read assumption note in the spec (doc-only stage).
5. `S5` — leave open; revisit only if observed in production.

`S2` stage-2 (the lazy-CIGAR refactor) is parked until long-read
support is decided.
