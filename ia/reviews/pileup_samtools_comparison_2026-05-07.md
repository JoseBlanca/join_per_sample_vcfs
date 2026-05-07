# Algorithm Review: pileup walker vs. samtools

**Date:** 2026-05-07
**Reviewer:** Claude (algorithm-comparison study)
**Module reviewed:** `src/per_sample_caller/pileup`
**Reference codebase:** `samtools/` (in-tree copy of samtools, not htslib)
**Status:** Advisory — no defects; a small backlog of optional improvements. As of 2026-05-07 the backlog is mostly worked: `S1`, `S3`, `S4`, `S5`, and `S6` are closed (see each finding's Resolution); only `S2` (doc-only stage) remains open.

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
- **Status:** Closed in commit `5c90fbb` (2026-05-07).

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

**Resolution.** Implemented as proposed. Added
`HIGH_WATER_WARN_THRESHOLD = MAX_ACTIVE_SLOTS * 3 / 4` and a one-shot
`high_water_warned: bool` flag to `SlotAllocator` ([slot_allocator.rs](../../src/per_sample_caller/pileup/slot_allocator.rs)).
The warning is `eprintln!` rather than `tracing::warn!` because the
project has no logging crate yet — easy to swap later. The flag is
deliberately preserved across `reset()` so the warning fires at most
once per run, not once per chromosome. Two new tests in the same
file cover the threshold-crossing and reset-preservation properties.

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
- **Status:** Closed in commit `7d875ac` (2026-05-07) — superseded
  by upstream tightening rather than implemented as proposed.

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

**Resolution.** A check of who actually consumes the toggles
(`grep -rn 'drop_secondary\|drop_supplementary' src/`) showed the
two booleans were referenced *only* inside `cram_input.rs` — no
CLI flag, no library caller, no production code. They were public
configuration with no documented use case, so the chosen fix was
strictly stronger than the proposed `debug_assert!`: remove the
toggles from `CramMergedReaderConfig` and unconditionally drop
both flag classes in the cascade. That makes the misconfiguration
*unrepresentable* upstream rather than asserting against it
downstream at the walker boundary, and avoids the public-shape
churn on `PreparedRead`. The drops still surface via
`FilterCounts.secondary` / `.supplementary` so users can audit how
many records were removed. Plan updated at
[../feature_implementation_plans/per_sample_caller_cram_input.md](../feature_implementation_plans/per_sample_caller_cram_input.md)
to keep the original 5-toggle design as a historical reference and
explain the two-step revision (`drop_unmapped` removed in M4 on
2026-05-01, `drop_secondary` / `drop_supplementary` in `7d875ac`).
The producer-agnostic walker-boundary guard remains worth
reconsidering only if a second producer of `PreparedRead` (e.g. a
BAM input path) is added.

### `S4` — Verify (don't necessarily add) a reference-fetch LRU cache

- **Priority:** Medium (depends on benchmark)
- **Effort:** Small to investigate; small to add if needed
- **Status:** Closed (no commit) — no action needed; investigation
  on 2026-05-07 confirmed the dependency already caches.

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

**Resolution.** No code change. Two findings from the
investigation:

1. **No production `RefBaseFetcher` exists yet.** A `grep -rn 'impl
   RefBaseFetcher' src/` returns only `MockFasta`
   ([tests.rs:39](../../src/per_sample_caller/pileup/tests.rs#L39)),
   the in-memory test fetcher. `pileup::run` has no production
   caller; the walker is exercised only from end-to-end tests. The
   CRAM input layer does build a `noodles_fasta::Repository`
   ([cram_files.rs:119–124](../../src/per_sample_caller/cram_files.rs#L119-L124),
   [cram_input.rs:609–616](../../src/per_sample_caller/cram_input.rs#L609-L616)),
   but it is consumed by noodles' CRAM block decoder, not exposed
   to the walker.
2. **`noodles_fasta::Repository` already caches whole contigs
   unboundedly.** Source-grounded against the on-disk crate
   (`~/.cargo/registry/.../noodles-fasta-0.60.0/src/repository.rs`):
   `Repository` wraps `Arc<RwLock<AdapterCache>>` whose `cache`
   field is `HashMap<Vec<u8>, Arc<Sequence>>`. `Repository::get`
   first checks the cache under a read lock; on miss it takes the
   write lock, calls the adapter (`IndexedReader::get` →
   `Reader::query` → seek-and-read the contig from the FASTA via
   the `.fai` index), and inserts the resulting `Arc<Sequence>`.
   There is **no eviction policy** — once a contig is fetched it
   stays in memory until the `Repository` is dropped.

So a future production fetcher that does the obvious thing —
`self.repo.get(name)?` then slice the returned `Arc<Sequence>` —
pays one disk read per contig and then in-memory slicing for
everything afterwards. samtools' 3-slot LRU was solving a problem
we don't have: our walker is sequential single-pass per
chromosome (per [ia/specs/pileup_walker.md](../specs/pileup_walker.md)
§"Chromosome boundaries"), and noodles' "permanent cache per
contig" is strictly stronger than a 3-slot LRU for that pattern.

**Adjacent observation, promoted to its own finding.** Repository's
unbounded cache means a future production fetcher would
accumulate ~3 GB of base data across 24 human chromosomes if it
just shared one `Repository` across the whole run. Initially
captured here as a "wait until we see the memory bill" caveat,
but on closer inspection (a) the design rule is so cheap that
designing for it up front is strictly better, and (b)
`noodles_fasta::Repository::clear()` is in fact public — missed on
the first pass — so the eviction is a one-line operation. Split
out as **S6** below.

### `S6` — Chromosome-boundary eviction in the production `RefBaseFetcher`

- **Priority:** Medium (latent memory blow-up on whole-genome runs)
- **Effort:** Small (one new module, one `Cell<Option<u32>>`,
  one `Repository::clear()` call, a handful of tests)
- **Status:** Closed in commit `453715b` (2026-05-07).

**Observation.** S4's investigation surfaced that
`noodles_fasta::Repository`'s contig cache (`HashMap<Vec<u8>,
Arc<Sequence>>`) is *unbounded* — once a contig is fetched it stays
in memory until the `Repository` drops. The walker is sequential
single-pass per chromosome (per
[ia/specs/pileup_walker.md](../specs/pileup_walker.md)
§"Chromosome boundaries") and never revisits a contig within a run,
so a future production `RefBaseFetcher` that does the obvious
"build one `Repository` at startup, share it across chromosomes"
thing would accumulate ~3 GB of base data across the 24 human
chromosomes for no functional benefit.

**Proposal.** Bake the eviction into the production fetcher rather
than waiting to discover the memory bill on real samples. Once one
chromosome's records are done, the previous contig's
`Arc<Sequence>` is no longer needed; drop it before loading the
next. `noodles_fasta::Repository::clear()` is public (we missed it
on the first pass of S4), so a single Repository can be reused for
the whole run with its cache reset at chromosome boundaries.

Concretely: a `ChromBoundaryRefFetcher` that holds a single
`Repository` plus a `Cell<Option<u32>>` tracking the current
`chrom_id`. On a `fetch` whose `chrom_id` differs from the cell,
call `repository.clear()` and update the cell, then proceed to
`repository.get(name)`. Steady-state cache size is exactly 1
contig.

**Rationale.** Cheap to do correctly the first time; expensive (in
RSS, on whole-genome runs) to do wrong and then have to retrofit.
The hook is natural: the chromosome change is detected by the
fetcher itself from the `chrom_id` argument, so no extra plumbing
into the walker is needed. Performance during a chromosome is
unchanged — the cache holds the active contig in memory between
fetches just as before; only the previous chromosome's bytes are
freed.

**Risk.** Effectively none. Cell-based interior mutability is fine
because the walker is single-threaded per call (`run<I, F>` does
not require `F: Send + Sync` and the walker state is itself not
`Send`/`Sync`). If a future architecture parallelises across
chromosomes within a single fetcher, the contract changes — but
that would require a different design anyway because the
single-Repository-with-clear approach assumes serial chromosome
visits.

**Resolution.** Implemented as proposed in
[src/per_sample_caller/ref_fetcher.rs](../../src/per_sample_caller/ref_fetcher.rs):
`ChromBoundaryRefFetcher` holds one `fasta::Repository` (built
from an `IndexedReader` over the `.fa` + `.fai` files) plus a
`Cell<Option<u32>>` for the current `chrom_id`. On a
chromosome-changing fetch it calls `Repository::clear()` and
updates the cell, then resolves the contig name via the project's
`ContigList` and slices the returned `Arc<Sequence>`. Tests cover
the eviction invariant by inspecting `Repository::len()` (exposed
as `cached_contig_count`): the cache stays at one contig across
chrom changes, regardless of how many chromosomes have been
visited so far.

### `S5` — Per-column depth cap (originally framed as per-allele)

- **Priority:** Low (originally) → Medium (after the per-allele
  vs. per-column nuance was worked out)
- **Effort:** Medium
- **Status:** Closed in commit `8ba0bc0` (2026-05-07).

**Observation.** samtools' mpileup uses a separate, lower cap for
indel columns (`MPLP_MAX_INDEL_DEPTH = 250` vs `MPLP_MAX_DEPTH = 8000`)
because homopolymer-context indels are pathologically deep in real
data. We have a uniform per-walker cap (active slot count) and no
per-allele bound.

**Original proposal.** Hold off. Add the cap only if we see a
Stage 2 failure or a runtime spike caused by a single allele's
`num_obs` blowing up. Track this finding here so we don't forget
the precedent exists.

**Original rationale.** Premature optimisation otherwise; Stage 2
should be robust to high allele depth in principle.

**Original risk.** Adding it now would mean choosing a threshold
without data, which is the worst time to choose one.

**Resolution.** Adopted samtools' caps as defaults rather than
waiting on our own evidence. Reasoning: samtools encodes 15+
years of real-world calling experience; we're not "picking a
threshold without data", we're inheriting a threshold that has
data behind it (just not our data). When we have our own deep
samples we can revisit the numbers; until then samtools' values
are a defensible starting point and a hardcoded "no cap" is
strictly worse than a borrowed-with-attribution cap.

A first-pass design that would have applied the cap *per allele*
was rejected after a more careful look: per-allele clipping
silently biases the allele-frequency estimate at any column
where the dominant allele exceeds the cap (a 99/1 column with
cap=250 would report ~71/29 instead of ~99/1). Per-*column*
truncation is what samtools does for the same reason — it
preserves ratios in expectation. Implementation: truncate the
contributor list to the first N items at fold time;
deterministic-first-N is approximately unbiased because reads in
our active set are not allele-correlated in iteration order, so
no random-sample machinery is needed.

Caps live on a new `WalkerConfig` struct
([pileup/mod.rs](../../src/per_sample_caller/pileup/mod.rs),
defaults `max_snp_column_depth = 8000`,
`max_indel_column_depth = 250`). The cap detection (indel
column? then use the tighter cap) runs at fold time in
[walker.rs:column_depth_cap](../../src/per_sample_caller/pileup/walker.rs).
A new `RunSummary.column_depth_truncations` counter tracks how
many columns hit the cap so QC pipelines can flag pathologically
deep regions. With the current `MAX_ACTIVE_SLOTS = 4096` only the
indel cap can actually fire — the SNP cap is future-proofing for
if/when the slot cap is raised.

This finding also drove the introduction of the per-stage
config-struct convention (mirrors `CramMergedReaderConfig`); the
project's [calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md)
gained a "Configurable parameters" section codifying the rule
that new tunables go in stage configs rather than as bare `pub
const`s.

### `S7` — Adopt samtools' BQ-combining math for match-only mate-overlap

- **Priority:** Medium (statistical correctness improvement; small
  code change)
- **Effort:** Small (one struct field, a helper or two,
  per-position resolver update; no plumbing-shape changes)
- **Status:** Closed in commit `_TBD_` (2026-05-07).

**Observation.** Reading [htslib/sam.c:5803-5942](../../htslib/sam.c#L5803-L5942)
(`tweak_overlap_quality`) shows our walker does not match samtools
on two specific points within a match-only mate-overlap pair at
the current walker position:

1. **Bases agree.** samtools sums the two BQs (capped at 200) and
   stores the sum on one mate, zero on the other
   ([sam.c:5917-5921](../../htslib/sam.c#L5917-L5921)). The
   surviving observation carries *combined* confidence, reflecting
   that two independent calls of the same base genuinely give
   higher confidence than one — Σ-Q-with-cap is the correct
   combination of independent log-error evidence.
2. **Bases disagree.** samtools picks the higher-BQ mate
   per-position and **scales its kept BQ by 0.8**
   ([sam.c:5927, 5931, 5935](../../htslib/sam.c#L5927)) — a
   "we trust this less because the mate disagreed" haircut — and
   zeroes the loser. Our walker keeps the winner's BQ unscaled
   and zeroes the loser.

Our current walker
([walker.rs:393-492](../../src/per_sample_caller/pileup/walker.rs#L393-L492))
treats both cases identically: pick the higher-BQ side, zero the
other, leave the winner's BQ unchanged. That's conservative on
the agree case (we under-credit the combined confidence) and
slightly aggressive on the disagree case (we don't apply the
mismatch-haircut samtools applies). The original review's "our
per-position approach is strictly richer than samtools'" claim
came from a misreading: samtools is also per-position (in the
disagree case) and per-pair-with-BQ-summing (in the agree case);
neither matches the picture of "one global winner per pair".

**Proposal.** Match samtools' math, while keeping our per-walker-pos
resolver shape (we don't pull the decision forward to admission —
that's still the rejected R1):

- **Agree case.** When both mates' Match-event bases agree at
  walker_pos, sum their BQs, cap at 200, store on the keeper
  (`bq_override_at_walker_pos = Some(combined)`); zero the other
  (existing `bq_zero_in_window` mechanism). Keeper choice is
  deterministic; samtools uses a qname hash, but for our model
  the existing `is_first_mate → alignment_start` tie-break is
  equivalent (statistical answer is the same regardless of which
  side carries the combined BQ).
- **Disagree case.** Higher-BQ mate's BQ at walker_pos is scaled
  by 0.8 (samtools-faithful truncation: `(bq as f64 * 0.8) as u8`,
  matching the C `uint8_t` cast); loser zeroed.

A per-position override mechanism (a new
`ReadContribution::bq_override_at_walker_pos: Option<u8>` field) is
needed because the open-record fold pulls window events from the
cursor with original BQs; the walker_pos event needs a post-fetch
BQ rewrite to apply the combined or scaled value.

**Rationale.** samtools' agree-case BQ sum is the right
log-error combination for two independent observations. Our
current behaviour silently underestimates the surviving
observation's confidence in the most common mate-overlap case
(agreement). The disagree-case 0.8 haircut is a smaller effect
but the same direction: samtools' math is principled, ours is
ad-hoc.

**Risk.** Low. Changes downstream `q_sum` values for
mate-overlap pairs (more negative for the agree case, slightly
less negative for the disagree case). `num_obs` is unchanged in
both cases (both mates still count as observations). Existing
tests using BQ values where MQ dominates `q_sum` will pass
unchanged; tests at low MQ_log_err that pin specific `q_sum`
values pin the new math.

**Limitation we accept.** samtools mutates the entire overlap
region's BQ array at admission. Our per-position resolver only
applies the combine/scale at the *current* walker_pos. For
window events at *other* positions in an open record's footprint
(e.g., a deletion record spanning multiple positions), the
keeper/winner's events use original BQ. This means our fold can
under-credit the keeper's window-events confidence relative to
samtools when the open record's footprint extends beyond
walker_pos. The discrepancy is bounded (only relevant for
multi-position records over indels) and was deemed not worth the
complexity of cross-position BQ caching today.

## 5. Findings explicitly *rejected*

For the record, so we don't relitigate them:

- **R1 — Pull mate-overlap detection forward to admission.** Rejected
  in favour of the per-position fold (see `A2`). Note: the original
  rejection text claimed our per-position approach was "strictly
  richer" than samtools'; that claim was based on a misreading of
  htslib (which is also per-position in the disagree case and
  per-pair-with-BQ-summing in the agree case). The corrected
  framing — samtools' admission-time approach trades data
  immutability for cache locality — still doesn't make R1 worth
  pulling forward, but the BQ math itself was worth importing
  (see `S7`).
- **R2 — Switch to a per-position callback model.** Rejected — our
  per-record closure rule is the foundation of `PileupRecord` /
  `AlleleObservation` / phase chains and is committed to in the
  spec.
- **R3 — Recompute BAQ inside the walker.** Rejected — the spec
  ([mod.rs:63](../../src/per_sample_caller/pileup/mod.rs#L63),
  [ia/specs/pileup_walker.md](../specs/pileup_walker.md)) commits to
  the walker treating `bq_baq` as opaque; BAQ is upstream.

## 6. Why we keep the per-record closure model (R2 in depth)

R2 is the largest of the rejected findings — moving to a samtools-style
per-position callback would touch every file in `pileup/`. It is worth
spelling out *why* we don't, because the answer is not "we got there
first and inertia wins"; it is that the two designs serve different
output contracts.

**What phase chains buy us, and where they're consumed.**
[calling_pipeline_architecture.md:416](../specs/calling_pipeline_architecture.md#L416)
makes the per-allele phase chain id a first-class field on every
`AlleleObservation`. It is the mechanism by which the merger
(Stage 5) decides whether two heterozygous calls in the same sample
co-occur on the same haplotype, *without* re-reading the alignments.
That is the entire point: compound-haplotype reasoning at merge time
([architecture §"Compound haplotype evidence"](../specs/calling_pipeline_architecture.md#L1008-L1043)),
chain-evident vs. chain-broken classification per sample, and the
ability for downstream stages to follow phase across long stretches
without ever loading a BAM. Stage 1 is where these chain ids are
assigned, because Stage 1 is the only stage that still has reads in
hand. Strip phasing out of Stage 1 and the entire downstream
contract collapses.

**Why samtools doesn't carry phase.** mpileup's output contract is a
flat per-position record: at each reference column, here are the
stacked bases, qualities, and indel events. There is no slot in that
output for "this base on read X at column 100 is the same haplotype as
this base on read X at column 200" — the link between columns has
already been thrown away. Phasing is left to downstream tools
(WhatsHap, HapCUT, GATK's local reassembly, DeepVariant's per-read
tensors), each of which re-reads the alignments rather than consuming
mpileup output. samtools' callback model is *correct for samtools*
precisely because phasing was never on its contract: a stateless
per-column callback is the cheapest, simplest interface that meets
that contract. The moment phase is on the contract, the callback
needs cross-position state, and the simplicity is gone.

**Why adopting R2 would not save the work.** The performance costs
attributable to our model — the `folded_reads` map, the BTreeMap of
open records, the mate-overlap map rebuilt per position, the chain
slot allocator with `pending_mates`, the slot-cooling buffer — are
overwhelmingly the cost of *tracking read identity across positions*,
not of the record-shaped emit API per se. A samtools-style callback
that also produced phase chain ids would have to maintain the same
state internally. The shape of the output (records vs. positions)
is downstream of the requirement (phase chains across positions),
and changing the shape doesn't lift the requirement. So R2 would
buy us a re-plumbing exercise without lifting the bookkeeping it is
sometimes assumed to lift.

**What R2 *would* trade away.** Records aren't just an emit format;
they are how compound events (an indel and its anchoring match, a
SNP next to an MNP) get folded into a single `PileupRecord` with a
coherent allele list. A per-position callback flattens that — every
column is an independent event, and the per-sample caller would have
to re-stitch compound alleles itself, with less information than the
walker has at fold time (e.g. the per-allele `chain_slots` set). We
would be moving complexity downstream rather than removing it.

**Conclusion.** R2 is the right design *for samtools' output
contract*. It is the wrong design for ours. The per-record closure
model is committed to in
[pileup_walker.md](../specs/pileup_walker.md) and
[calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md)
not as a stylistic preference but because it is the natural shape
for the data the rest of the pipeline consumes. We keep it.

If walker throughput ever becomes a bottleneck, the place to look is
*inside* the closure model — `folded_reads` hashing, per-position
mate-overlap map allocation, BTreeMap range queries — not at the
API shape. Those are local optimisations that don't touch the
data-model commitment.

## 7. Suggested execution order

If acted on in order, each finding is independent and can be its own
small commit:

1. ~~`S1` — soft high-water warning (one-line user-visible win).~~ **Done in `5c90fbb` (2026-05-07).**
2. ~~`S3` — supplementary/secondary guard (cheap defence-in-depth).~~ **Done in `7d875ac` (2026-05-07), via the stronger upstream-tightening route — see S3's Resolution.**
3. ~~`S4` — investigate ref-fetch caching; add if needed.~~ **Closed without code change (2026-05-07) — `noodles_fasta::Repository` already caches; see S4's Resolution.**
4. ~~`S6` — chrom-boundary eviction in the production `RefBaseFetcher`.~~ **Done in `453715b` (2026-05-07) — split out of S4 once `Repository::clear()` was confirmed public; see S6's Resolution.**
5. ~~`S5` — per-column depth cap.~~ **Done in `8ba0bc0` (2026-05-07) — adopted samtools' caps as defaults rather than waiting for our own data; the per-allele framing in the original proposal was rejected as statistically biased; see S5's Resolution.**
6. `S2` — long-read assumption note in the spec (doc-only stage).

`S2` stage-2 (the lazy-CIGAR refactor) is parked until long-read
support is decided.
