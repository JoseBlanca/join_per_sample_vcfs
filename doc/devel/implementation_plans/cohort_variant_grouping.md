# Variant grouping (Stage 4 — cohort-side overlap bundler over per-position pileups)

Proposal date: 2026-05-16.

## Domain intent

The second cohort-side component of the multi-sample SNP caller (see
[doc/devel/specs/calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md)).
It sits directly downstream of the multi-way per-position merger and
bundles overlapping per-position pileups into `OverlappingVarGroup`s:

```
… ─► PerPositionMerger ─► (Stage 3 DUST filter, when it lands) ─► VariantGrouper ─► OverlappingVarGroup stream
                                                                                       │
                                                                                       ▼
                                                                              Stage 5 (per group, parallel)
```

Stage 4's contract, from
[calling_pipeline_architecture.md:988-1051](../specs/calling_pipeline_architecture.md#L988-L1051):

- Walk the upstream per-position stream in genomic order.
- Bundle positions whose alleles share reference coverage in at
  least one sample into a single group.
- Emit a stream of independent groups in genomic order; everything
  downstream (Stage 5 onward) consumes one group at a time.

Stage 4 is sequential by nature — the right boundary of a group can
only be decided after walking past it — but the output stream is
embarrassingly parallel-friendly, which is what Stage 5 cashes in on
via rayon.

## Why now

The multi-way per-position merger has shipped (commit `427d95f`,
[src/cohort/per_position_merger.rs](../../src/cohort/per_position_merger.rs)).
It emits exactly the input Stage 4 wants — `Result<PerPositionPileups,
MergerError>` in strictly increasing `(chrom_id, pos)` order — so the
grouper can be built directly against the merger's output without
touching Stage 3 (DUST) first.

Building Stage 4 next also pins down the `OverlappingVarGroup` type,
which is the input contract for Stage 5 (allele unification + per-sample
likelihood reconstruction). Stage 5 is the largest remaining piece of
the pipeline; unblocking its API surface is high-leverage.

Stage 3 (DUST) is deliberately deferred. It is a pure iterator adaptor
between merger and grouper — constant memory, no API surface for
downstream — and can be slotted in later without rework provided this
plan keeps the grouper generic on its upstream iterator (see "API
shape" below).

## What's already in place

- **Merger** — [src/cohort/per_position_merger.rs](../../src/cohort/per_position_merger.rs).
  Emits `Result<PerPositionPileups, MergerError>` in strict
  `(chrom_id, pos)` order. `PerPositionPileups` carries `chrom_id`,
  `pos`, and `per_sample: Vec<Option<PileupRecord>>` indexed by
  sample order.
- **`PileupRecord`** — [src/per_sample_caller/pileup/mod.rs:403-429](../../src/per_sample_caller/pileup/mod.rs#L403-L429).
  Has `ref_span() -> u32` derived from `alleles[0].seq.len()`
  (REF is always slot 0). This is the per-record reach used to
  decide group extension.
- **Old gVCF-era grouper** — [src/variant_grouping.rs](../../src/variant_grouping.rs)
  and its merging consumer
  [src/genotype_merging.rs](../../src/genotype_merging.rs),
  with the integration test suites in
  [tests/variant_group_test.rs](../../tests/variant_group_test.rs)
  and
  [tests/genotype_merging_test.rs](../../tests/genotype_merging_test.rs).
  Reference material only. The grouper module belongs to the
  gVCF → VCF binary path that this pipeline is replacing, but
  its **integration tests are not disposable**: they encode the
  catalogue of variant-merging edge cases (overlapping
  deletions across samples, compound indels within a sample,
  phase-chain preservation, missing alleles in merged
  haplotypes) that Stage 5 has to reproduce. The old code stays
  in the tree until Stage 5 lands with equivalents — see
  §"Lessons to extract from the old grouper before deleting it"
  for what to carry forward and §"Deferred cleanup: remove the
  gVCF path after Stage 5" for the timing. What carries into
  the new *grouper* directly: the single-pass extension shape,
  the closing rule, and the `skip_group = !is_variable`
  performance idea. What does **not** carry into the new
  grouper: the per-sample VCF iterators, the rayon-parallel
  peek across them (the new grouper sees a single pre-merged
  stream, so there is nothing to peek in parallel), the
  `VarIteratorInfo` shape, and the `Variant` data type.
- **Chromosome id space** — `chrom_id: u32` is the merger's
  contract; chromosome agreement across samples is already enforced
  upstream by `check_chromosome_agreement`. The grouper only ever
  sees one normalised id space.

## Lessons to extract from the old grouper before deleting it

The old gVCF path has two distinct kinds of value: **algorithmic
ideas** that this plan reproduces in the new grouper, and
**test-encoded edge cases** that Stage 5 has to reproduce when it
lands. The grouper code can be deleted as soon as the algorithmic
ideas are in place (this plan's commit), but the integration tests
must stay until Stage 5 has equivalents — see §"Deferred cleanup:
remove the gVCF path after Stage 5" below for the timing.

### Algorithmic ideas reproduced in this plan

Three ideas in [src/variant_grouping.rs](../../src/variant_grouping.rs)
deserve explicit acknowledgement so they survive the deletion:

1. **Drop non-variant positions as early as possible, not at
   close-time.** The old grouper's
   `skip_group = !is_variable` branch
   ([variant_grouping.rs:233-273](../../src/variant_grouping.rs#L233-L273))
   skips positions that are entirely non-variable *before* they
   ever enter a group. Most of a WGS genome is pure REF at every
   covered base; dropping these at the iterator boundary instead
   of after a group is built avoids per-position allocations on
   the hot path. The new grouper inherits this directly — see
   §"Pure-REF positions: dropped at the iterator boundary" and
   the algorithm below.
2. **Single-pass extension with sequential close.** The shape is
   right: one open group, one end-extension per pulled item,
   close-and-emit on the first item that doesn't reach back.
   This is the only correct way to bound the right boundary
   without speculation.
3. **Terminate-on-first-error latching.** The old grouper sets
   `self.done = true` on any error and refuses to yield further
   items
   ([variant_grouping.rs:153-156](../../src/variant_grouping.rs#L153-L156)).
   Same contract here, mirroring the merger.

### Test-encoded edge cases Stage 5 has to reproduce

The integration suite in
[tests/genotype_merging_test.rs](../../tests/genotype_merging_test.rs)
encodes ~25 distinct merging scenarios that Stage 5 will need to
handle on the new pipeline. They are not directly portable — the
old tests build `Variant` records and the new pipeline works on
`PerPositionPileups` / `OverlappingVarGroup` with chain ids — but
the *scenarios* they pin down are exactly the cohort-level
correctness surface Stage 5 has to cover. Categorised:

- **Allele-set unification across samples.** Sample A has a
  deletion at p, sample B has a SNP at p+2 that falls inside
  A's deletion. Stage 5 produces one merged record with both
  alleles in the unified set. Source tests:
  `test_overlapping_deletions`, `test_simple_deletion`,
  `test_deletion_len_2`, `test_simple_insertion`,
  `test_simple_merge`.
- **Within-sample compound alleles.** A single sample carries
  two deletions in the same group; the per-sample genotype
  must reflect both as a compound haplotype, not as
  independent calls. Source tests:
  `test_two_deletions_in_same_sample`,
  `test_two_overlapping_deletions_in_same_sample`,
  `test_het_deletion`.
- **Phase-chain preservation across positions.** Eight tests
  exercise how phase relationships across positions inside a
  group survive merging — when phase is set, when it isn't,
  when a sample lacks the phasing record at one position,
  when a het position interrupts the chain. The new chain-id
  machinery (see Stage 1 §"Chain ids" in the spec) is the
  replacement mechanism; these tests are the truth set it must
  match. Source tests:
  `test_phase_single_het_no_phase_needed`,
  `test_phase_kept_between_hets`,
  `test_phase_false_on_first_position_still_solvable`,
  `test_phase_first_het_not_at_first_position`,
  `test_phase_broken_between_hets_missing_genotype`,
  `test_phase_broken_in_one_sample`,
  `test_phase_lost_in_het`,
  `test_phase_lost_in_het_and_not_recovered`,
  `test_phase_lost_in_het2`,
  `test_phase_maintained_despited_not_phased_in_hom_position`,
  `test_phase_conserved_in_het`.
- **Allele-set hygiene.** What happens when an input record
  carries a non-REF allele that no read supports in this
  cohort, or when a sample's compound haplotype names an
  allele that nothing else in the group emits. Source tests:
  `test_allele_merge_with_non_ref_filtered`,
  `test_missing_allele_in_merged_haplotype`.
- **Group-shape correctness.** What counts as one group vs.
  several, ordering across chromosomes, non-variable groups
  being skipped. Largely covered by this plan's tests already
  but the old suite's `test_grouper_preserves_genomic_order`,
  `test_non_variant_vars_are_removed`,
  `test_group_merging_creates_correct_number_of_vars`,
  `test_binning_with_deletion_spanning_several_snps`, and
  `test_non_variable_groups_are_skipped` are useful
  cross-checks.

**Stage 5's plan must include a porting checklist** that walks
each of the bullets above and pins how the new pipeline
reproduces the scenario — either with a directly equivalent test
on `OverlappingVarGroup` input, or with an explicit waiver where
the new chain-id-based model makes the old test moot. Only when
that checklist is complete does the old test suite's educational
value run out.

## Algorithm: streaming single-pass extension with seed-time filtering

Maintain one open group at a time. The grouper has two phases per
emitted item: **find a variant-bearing seed**, then **extend and
close**.

**Phase A — seed search.** Pull from upstream until we find a
`PerPositionPileups` whose `has_variant_observation` predicate
is true (i.e. at least one `Some` slot has `alleles.len() > 1`).
Items that fail the predicate are dropped immediately — no group
is allocated, no record is cloned, nothing downstream sees them.
Most of the genome traverses this phase as a tight "predicate +
drop" loop.

**Phase B — extend and close.** Once a seed is found, initialise a
group with `start = pos`, `end = pos + max_ref_span - 1` (inclusive,
1-based). Then pull further upstream items:

1. If the item is on the **same chromosome** and `pos <= group.end`:
   add it to the open group (verbatim, regardless of whether it is
   itself variant — REF-only positions covered by a deletion's
   span belong in the group), and extend
   `group.end = max(group.end, pos + max_ref_span - 1)`.
2. Otherwise (different chromosome, or `pos > group.end`): close
   the open group, stash the just-pulled item as the next seed
   candidate (it goes through Phase A on the next call), and emit
   the group.

The closing rule is symmetric with Stage 1's
[open-record closure](../specs/calling_pipeline_architecture.md#L373-L468):
the grouper holds a group open until it can prove no future
position will be drawn into it. The proof: upstream items arrive in
strictly increasing `(chrom_id, pos)`, so any future item has
`(chrom', pos') > (chrom, pos)`; an item can extend `group.end`
only if `pos' <= group.end`, so once `pos' > group.end`, the group
is closed for good.

**Why seed-time filtering, not close-time.** A close-time filter
would allocate the group's `Vec<PerPositionPileups>`, push every
covered pure-REF position into it, then walk the records again at
close-time to decide whether to discard them all. On WGS data the
vast majority of seeds would be discarded, so this would mean
allocating-then-dropping a `Vec<PerPositionPileups>` and N-slot
`Vec<Option<PileupRecord>>`s for ~every covered base. The
seed-time filter pushes the test ahead of any allocation: the
common path is one cheap predicate call and a `continue`, with
zero heap traffic per dropped position. This is the perf idea the
old grouper already encoded with `skip_group = !is_variable` and
it's the same call here.

Memory: `O(widest_emitted_group_in_records × max_alleles_per_position)`,
typically a few hundred bp of reference. Independent of cohort
size (because `PerPositionPileups` already collapses N samples
into one item per position), of genome length, and of how many
pure-REF positions sit between variant-bearing groups (those are
dropped without buffering).

The old grouper used a parallel peek across per-sample VCF iterators
because it had N independent streams to inspect. That parallelism is
already paid in the merger — by the time records reach Stage 4 they
are one stream per position, not N. The grouper therefore stays
single-threaded; rayon comes back in at Stage 5 once independent
groups are out the door.

## What goes in a group

An `OverlappingVarGroup` carries:

- `chrom_id: u32` — chromosome id, in the merger's normalised
  id space.
- `start: u32` — 1-based inclusive position of the group's first
  record.
- `end: u32` — 1-based inclusive position the group's REF spans
  reach to. Derivable from `records` but cached so Stage 5 doesn't
  have to recompute it.
- `records: Vec<PerPositionPileups>` — every per-position pileup
  pulled into this group, in increasing `pos` order. Each one
  carries the full N-slot per-sample vector verbatim.

The per-sample slot vectors come through unchanged: Stage 5 needs to
ask "which samples observed which alleles at which positions in this
group?" and the existing `PerPositionPileups` shape answers that
in `O(1)` per sample lookup. No projection or reshaping happens at
Stage 4.

## Configurable parameters

Per the project's
[per-stage config convention](../specs/calling_pipeline_architecture.md#L81-L122),
Stage 4 ships with a single `GrouperConfig` struct, mirroring the
shape of [`WalkerConfig`](../../src/per_sample_caller/pileup/mod.rs#L105-L146).
The CLI binding is out of scope for this plan (no cohort CLI exists
yet) but the field is part of the library API from day one so test
code, library users, and a future CLI all bind to the same surface.

```rust
pub const DEFAULT_MAX_VAR_GROUP_SPAN: u32 = 10_000;

#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct GrouperConfig {
    /// Hard cap on the reference span of any emitted
    /// `OverlappingVarGroup`. A group whose `end - start + 1` would
    /// exceed this value triggers a [`GrouperError::VarGroupTooWide`]
    /// and the iterator latches `done`. Defaults to
    /// [`DEFAULT_MAX_VAR_GROUP_SPAN`].
    ///
    /// The bound is defensive: in normal operation it is never
    /// reached. Stage 1's `MAX_RECORD_SPAN = 5000` bp caps every
    /// single record's reach, and chain-id reach is bounded by the
    /// same value
    /// ([calling_pipeline_architecture.md:446-457](../specs/calling_pipeline_architecture.md#L446-L457)),
    /// so a single sample cannot contribute a record extending more
    /// than 5000 bp past its anchor. Cross-sample transitive
    /// extension chains can in principle compound, but in practice
    /// the old grouper never observed a real group exceeding a few
    /// hundred bp. The 2× headroom (10000 bp) reflects this: high
    /// enough to never fire on real data, low enough to bound
    /// pathological memory growth on adversarial inputs.
    pub max_var_group_span: u32,
}

impl Default for GrouperConfig {
    fn default() -> Self {
        Self {
            max_var_group_span: DEFAULT_MAX_VAR_GROUP_SPAN,
        }
    }
}
```

**Why hard-error, not log-and-skip.** The walker's analogous cap
([`WalkerConfig::max_record_span`](../../src/per_sample_caller/pileup/mod.rs#L116-L122))
is also a hard error rather than a logged skip. The reasoning carries
over: a group exceeding `max_var_group_span` is either pathological
input (worth investigating, not silently dropping) or an upstream bug
(worth surfacing loudly). A real cohort that hits the cap legitimately
should raise the cap explicitly via the config; that decision belongs
to the user, not to a silent drop. This is the same posture
[the spec describes for walker-internal caps](../specs/calling_pipeline_architecture.md#L128-L133):
"Exceeding it is a hard error … rather than silent memory blow-up on
pathologically deep regions." Stage 1's upstream read-span *filter*
(log-and-skip) is the permissive layer; the cap at the cohort level
is the defensive layer.

**Default value rationale.** `10_000` bp = 2× Stage 1's
`MAX_RECORD_SPAN`. The choice is documented above; revisit if real
data demonstrates either that 10000 is too tight (a legitimate
cohort tripping the cap) or too loose (pathological inputs
consuming excess memory before tripping).

## API shape

New module: `src/cohort/variant_grouping.rs`. (The old root-level
`src/variant_grouping.rs` keeps its name unchanged for the time
being — it will be retired with the rest of the gVCF path, see
"Out-of-scope follow-ups". The two module paths are unambiguous
even though the leaf name is the same.)

```rust
// src/cohort/mod.rs
pub mod per_position_merger;
pub mod variant_grouping;

// src/cohort/variant_grouping.rs
use crate::cohort::per_position_merger::{MergerError, PerPositionPileups};

pub struct VariantGrouper<I>
where
    I: Iterator<Item = Result<PerPositionPileups, MergerError>>,
{ /* fields private */ }

impl<I> VariantGrouper<I>
where
    I: Iterator<Item = Result<PerPositionPileups, MergerError>>,
{
    /// Construct a grouper over `upstream` with default tuning.
    /// Equivalent to `with_config(upstream, GrouperConfig::default())`.
    pub fn new(upstream: I) -> Self;

    /// Construct a grouper over `upstream` with explicit tuning.
    /// Use this when overriding `max_var_group_span` or any future
    /// `GrouperConfig` field.
    pub fn with_config(upstream: I, config: GrouperConfig) -> Self;

    pub fn config(&self) -> &GrouperConfig;
}

impl<I> Iterator for VariantGrouper<I>
where
    I: Iterator<Item = Result<PerPositionPileups, MergerError>>,
{
    type Item = Result<OverlappingVarGroup, GrouperError>;
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct OverlappingVarGroup {
    pub chrom_id: u32,
    /// 1-based inclusive.
    pub start: u32,
    /// 1-based inclusive.
    pub end: u32,
    /// In strictly increasing `pos` order.
    pub records: Vec<PerPositionPileups>,
}

impl OverlappingVarGroup {
    pub fn new(records: Vec<PerPositionPileups>) -> Self; // for tests
    pub fn span(&self) -> (u32, u32, u32) { (self.chrom_id, self.start, self.end) }
}

#[derive(thiserror::Error, Debug)]
#[non_exhaustive]
pub enum GrouperError {
    #[error("upstream: {0}")]
    Upstream(#[from] MergerError),

    /// A group's reference span would exceed
    /// `GrouperConfig::max_var_group_span`. Carries the group's
    /// chromosome and the seed/extension positions so the user can
    /// either investigate the locus or raise the cap.
    #[error(
        "variant group at chrom {chrom_id} starting at {start} would span \
         {attempted_span} bp (>{cap} bp cap); raise GrouperConfig::max_var_group_span \
         or investigate the locus"
    )]
    VarGroupTooWide {
        chrom_id: u32,
        start: u32,
        /// The `pos + ref_span - 1` the next extension would set as
        /// the new `end`, used to compute `attempted_span = (this -
        /// start + 1)`.
        attempted_end: u32,
        attempted_span: u32,
        cap: u32,
    },
}
```

**Two layers of errors.** The grouper has its own error type
(`GrouperError`) that wraps the upstream `MergerError` and adds
grouper-side variants. Today there is exactly one such variant
(`VarGroupTooWide`); future grouper-internal failure modes land here
without disturbing callers. The `#[from] MergerError` conversion
keeps idiomatic `?`-propagation easy.

**Latching.** Once any error is surfaced (upstream or grouper-side),
`done` is set and subsequent `next()` calls return `None`. Same
one-shot contract as the merger
([per_position_merger.rs:215-300](../../src/cohort/per_position_merger.rs#L215-L300)).
The in-progress group is dropped at the moment of error, never
partially emitted.

**Generic on the iterator type.** The merger's
[generic-iterator precedent](multi_way_per_position_iterator.md)
applies: testing is trivial with synthetic
`std::iter`-based fixtures, and a Stage 3 (DUST) adaptor that emits
`Result<PerPositionPileups, MergerError>` will compose without
touching this module. If DUST ever wants a richer error type, we
generalise the grouper's error then, not now.

Typical caller wiring:

```rust
let merger = PerPositionMerger::new(iters, sample_names, chromosomes)?;
let grouper = VariantGrouper::new(merger);
// Or: let grouper = VariantGrouper::with_config(merger, GrouperConfig { max_var_group_span: 50_000, ..Default::default() });
for group in grouper {
    let group = group?;
    // hand to Stage 5
}
```

## Algorithm details

State kept by the grouper:

```rust
struct VariantGrouper<I> {
    upstream: I,
    config: GrouperConfig,
    /// First item of the next group, already pulled from upstream
    /// but not yet folded into a group. `None` before the first
    /// emission and after exhaustion.
    pending_seed: Option<PerPositionPileups>,
    /// Latch: once set, all subsequent `next()` calls return `None`.
    /// Set on natural exhaustion *and* on the first error emitted —
    /// matches the merger's terminate-on-first-error contract.
    done: bool,
}
```

`VariantGrouper::new` does **not** prefetch from upstream. The first
`next()` call pulls the first seed. Rationale: the prefetch in
`PerPositionMerger::new` is justified by the need to populate per-reader
heads for the linear-scan min computation; the grouper has no equivalent
multi-source state to prime. Lazy init keeps construction infallible
and matches `std::iter::Iterator`-implementing adaptor convention.

`Iterator::next()`:

1. If `done`, return `None`.
2. **Phase A — seed search.** Repeatedly:
   - Source the next candidate: take `pending_seed` if `Some`,
     otherwise call `upstream.next()`.
   - On `None` → set `done`, return `None`.
   - On `Some(Err(e))` → set `done`, return `Some(Err(e))`.
   - On `Some(Ok(pp))`:
     - If `has_variant_observation(&pp)` → break out of the
       loop with `pp` as the seed.
     - Else → drop `pp` and continue the loop. (No allocation,
       no clone — the value is consumed and discarded.)
3. Initialise an in-progress group from the seed:
   - `chrom_id = seed.chrom_id`
   - `start = seed.pos`
   - `end = seed.pos + max_ref_span(&seed) - 1`
   - Check the cap up front: if `end - start + 1 >
     config.max_var_group_span`, set `done` and return
     `Some(Err(GrouperError::VarGroupTooWide { chrom_id, start,
     attempted_end: end, attempted_span: end - start + 1, cap:
     config.max_var_group_span }))`. This guards against the
     pathological case where a single seed record's `ref_span`
     alone exceeds the cap — Stage 1's `MAX_RECORD_SPAN` should
     prevent this, but the cap should not rely on that
     invariant holding silently.
   - `records = vec![seed]`
4. **Phase B — extend.** Loop pulling from `upstream`:
   - `None` → set `done`, emit `Ok(finalised group)`.
   - `Some(Err(e))` → set `done`, emit `Some(Err(e))`.
     The in-progress group is dropped — partial groups are
     never emitted because Stage 5 has no way to know the group
     was truncated.
   - `Some(Ok(pp))`:
     - If `pp.chrom_id != chrom_id` or `pp.pos > end`: stash
       `pp` in `pending_seed` (it will be re-evaluated by
       Phase A on the next `next()` call, possibly dropped if
       it is itself pure-REF) and emit `Ok(finalised group)`.
     - Else: compute `attempted_end = max(end, pp.pos + max_ref_span(&pp) - 1)`.
       If `attempted_end - start + 1 > config.max_var_group_span`,
       set `done`, drop the in-progress group, and return
       `Some(Err(GrouperError::VarGroupTooWide { chrom_id, start,
       attempted_end, attempted_span: attempted_end - start + 1,
       cap: config.max_var_group_span }))`.
       Otherwise extend `end = attempted_end` and push `pp` into
       `records`. Note: this branch folds in every item that the
       open group's REF span reaches, **regardless of whether the
       item is itself variant**. A REF-only position covered by a
       deletion's span is legitimate group content — Stage 5 needs
       those records to compute likelihoods for samples whose
       reads spanned the position without a deletion.

Helpers:

```rust
fn has_variant_observation(pp: &PerPositionPileups) -> bool {
    pp.per_sample
        .iter()
        .flatten()
        .any(|r| r.alleles.len() > 1)
}

fn max_ref_span(pp: &PerPositionPileups) -> u32 {
    pp.per_sample
        .iter()
        .flatten()
        .map(|r| r.ref_span())
        .max()
        .unwrap_or(1)
}
```

`has_variant_observation` short-circuits on the first slot with a
non-REF allele, so the hot path (pure-REF position across all
samples) walks every slot but does no allocation. For WGS data
this is the dominant case and the predicate's cost is the only
work the grouper does at most positions.

The `max_ref_span` `unwrap_or(1)` defends against an upstream item
with every slot `None` — which the merger cannot emit by
construction (it only emits at positions where at least one reader
had a record) but is not statically enforced. Falling back to
`ref_span = 1` keeps the grouper from panicking on a pathological
upstream and pins the group's reach to the position itself; such
an item also fails `has_variant_observation` and so is dropped at
seed time anyway.

**Per-item monotonicity is not re-checked here.** The merger already
guarantees strictly increasing `(chrom_id, pos)`
([per_position_merger.rs:241-258](../../src/cohort/per_position_merger.rs#L241-L258)),
and the grouper's only consumer (Stage 5) does not depend on a
second defense. If a future upstream variant relaxes monotonicity
this assumption surfaces as a wrong-grouping bug, not as silent
data loss — caught by the test cases below.

## Pure-REF positions: dropped at the iterator boundary

This is a performance decision more than a policy one. The
architecture spec at
[calling_pipeline_architecture.md:1010-1012](../specs/calling_pipeline_architecture.md#L1010-L1012)
mentions "trivial groups of size one" but does not pin behaviour on
positions that are pure-REF across every sample; the grouper drops
them and the rationale is essentially all about avoided work.

**What is dropped.** A `PerPositionPileups` item where every `Some`
slot has `alleles.len() == 1` (i.e. the slot holds only the REF
allele) and there is no group currently open that the item would
extend into. Such items never seed a group, never get cloned, and
never reach Stage 5.

**What is *not* dropped.** A pure-REF position that falls inside
an open group's REF span is folded into the group verbatim. Stage
5 needs those records to know which samples had read coverage at
the covered reference positions without supporting the spanning
deletion. The drop only happens for pure-REF positions outside any
open group's reach — the ones that would otherwise seed a
brand-new group with nothing to call.

**Why this matters for performance.**

1. **The hot path is short.** WGS samples cover the whole genome at
   every covered base ([Stage 2 record shape](../specs/calling_pipeline_architecture.md#L496-L529));
   the vast majority of positions are pure REF everywhere. With
   seed-time filtering, the per-position work at most positions
   is a single `has_variant_observation` predicate call plus a
   `continue` — no heap allocation, no clone, no group state to
   build up and tear down.
2. **No wasted allocations.** A close-time variant of the filter
   would push every covered pure-REF position into the in-progress
   group's `Vec<PerPositionPileups>` and then walk those records
   again at close-time only to discard the entire group. On WGS
   that is a per-base heap allocation that goes immediately to the
   bin. The seed-time filter pays nothing for the same correctness
   property.
3. **Bounded memory regardless of how sparse the variants are.**
   Pure-REF runs between variant-bearing groups are dropped
   without buffering, so genome-wide memory stays the same whether
   variants are dense or sparse — the grouper only ever holds the
   widest *emitted* group in memory.
4. **Precedent.** The old gVCF grouper already does this with
   `skip_group = !is_variable`
   ([variant_grouping.rs:233-273](../../src/variant_grouping.rs#L233-L273)).
   The seed-time version of the same idea slots into the
   single-pass extension naturally; the old grouper's parallel-peek
   shape made it slightly less direct, but the perf intent is
   identical.

**Opt-out flag** (CLI, not in this plan): a future
`--emit-ref-only-groups` flag can flip the policy for diagnostic
runs that want exhaustive coverage. Implementation is a single
boolean field on the grouper plumbed through from the CLI. Out of
scope for this plan — defer until Stage 5/CLI lands.

## Out of scope

- **Stage 3 (DUST) filter.** Separate iterator adaptor; not built
  here.
- **Cross-group context.** Stage 5 sees one group at a time. No
  cross-group state is maintained by the grouper.
- **Parallelisation across groups.** Stage 5 owns that, via rayon
  on its own input stream.
- **CLI flag for the pure-REF policy.** Deferred; the policy is
  hard-coded to "drop REF-only groups" in this plan and toggled
  later when there is a CLI to bind to.
- **CLI binding for `GrouperConfig`.** The library API ships in
  this plan; the cohort CLI subcommand that binds
  `--max-var-group-span` (and any future `GrouperConfig` fields)
  lands with the cohort entry point. Out of scope here.
- **The Stage 5 record shape.** This plan ships the
  `OverlappingVarGroup` input contract; what Stage 5 emits is a
  separate design.

## Test strategy

All tests are unit-level in `cohort/variant_grouping.rs`'s
`#[cfg(test)]` module. They use synthetic upstream iterators
(`std::iter`-based `Result<PerPositionPileups, MergerError>`
streams) constructed via small fixture builders, mirroring the
merger's existing test style at
[per_position_merger.rs:360-745](../../src/cohort/per_position_merger.rs#L360-L745).

Fixture builders to add:

```rust
fn ref_only_pp(chrom_id: u32, pos: u32, n_samples: usize) -> PerPositionPileups;
fn snp_pp(chrom_id: u32, pos: u32, n_samples: usize, sample_idx: usize) -> PerPositionPileups;
fn del_pp(chrom_id: u32, pos: u32, n_samples: usize, sample_idx: usize, ref_span: u32) -> PerPositionPileups;
```

Each builds a `PerPositionPileups` with the given variant shape at
the given sample slot and `None` in all other slots. `del_pp`
constructs a REF allele with `seq.len() == ref_span` (per the
indel-anchor convention).

Required cases:

- **Empty upstream.** Grouper yields `None` immediately.
- **Single trivial SNP.** One upstream item, one emitted group with
  `start == end == pos`, one record.
- **Two non-overlapping SNPs.** Two single-position groups.
- **Two SNPs sharing a position.** One group at that position
  with `records.len() == 1` and both samples' slots populated.
  (This case really tests the upstream merger's slot-merging more
  than the grouper, but pin it so the contract is explicit.)
- **Deletion drawing in a downstream SNP.** Deletion at p with
  `ref_span = 5` and SNP at p+3 → one group with `start = p`,
  `end = p+4`, two records.
- **Transitive chain.** SNP at p (REF span 1), deletion at p+0
  with span 3 (so end = p+2), SNP at p+1, deletion at p+2 with
  span 3 (extending end to p+4), SNP at p+3 — all bundled into
  one group with `end = p+4` and five records.
- **Deletion at the very end of a chromosome.** Group closes
  cleanly when the next upstream item is on a new chromosome
  even if `pos_next <= group.end_pos` numerically; the
  `chrom_id` change forces a close.
- **Multi-chromosome with no overlap.** Each chromosome's records
  appear in its own group(s); no group spans chromosomes.
- **Pure-REF positions filtered at seed time.** Three consecutive
  REF-only positions followed by a SNP → grouper emits exactly
  one group, for the SNP. (No empty group, no REF-only group.)
  Also pin: a long run of REF-only positions (e.g. 1000 of them)
  followed by a single SNP still emits exactly one group — this
  smoke-tests that the seed-search loop does not accumulate
  per-dropped-position state.
- **REF-only position drawn into a real group.** Deletion at p
  with span 5 covering REF-only positions p+1..p+3, plus a SNP
  at p+4 → one group containing all five records; REF-only
  records pulled in by an open group's extension are **not**
  dropped (the seed-time filter only applies when no group is
  open).
- **Mixed: REF-only run, variant group, REF-only run, variant
  group.** Two emitted groups; no off-by-one on the
  `pending_seed` handoff between Phase B close and Phase A
  re-entry.
- **`pending_seed` re-evaluation drops pure-REF.** A variant
  group at p closes because a pure-REF position at p+10 is past
  `end`; that pure-REF position lands in `pending_seed`, then
  Phase A on the next `next()` call drops it and continues the
  search. The next emitted group is the next variant-bearing
  position, not the stashed pure-REF one.
- **Upstream error mid-stream.** A synthetic iterator yielding
  `Err(MergerError::OutOfOrder { … })` after some records.
  Confirms: the error surfaces from `next()` wrapped in
  `GrouperError::Upstream`; the in-progress group is **not**
  emitted; subsequent `next()` calls return `None`.
- **Upstream error before any group closes.** Same shape, error
  is the first item. Confirms: no group emitted, `done` latches.
- **`VarGroupTooWide` triggered by extension.** Configure
  `GrouperConfig { max_var_group_span: 10, .. }`. Seed at p=100
  with `ref_span = 1`, then extend with records at p=101, p=105,
  and finally p=109 with `ref_span = 5` (so `attempted_end =
  113`, span = 14 > 10). Confirms: the error carries
  `chrom_id`, `start = 100`, `attempted_end = 113`,
  `attempted_span = 14`, `cap = 10`; subsequent `next()` calls
  return `None`; no group is emitted.
- **`VarGroupTooWide` triggered by seed alone.** Configure
  `max_var_group_span: 4`. Seed with `ref_span = 5` at p=100
  (so the seed-time check trips before any extension).
  Confirms: error carries `attempted_end = 104, attempted_span = 5`;
  no group emitted; `done` latches. Pin so the seed-time guard
  is not removed by accident in a future refactor.
- **`max_var_group_span` honoured at the boundary.** Configure
  `max_var_group_span: 5`. Build a group whose final
  `end - start + 1` lands at exactly `5`. Confirms: the group
  emits successfully (the cap is a strict-exceed check, not
  `>=`).
- **`end` correctness across re-extensions.** Group seeded with
  `ref_span = 2` at p, then a record at p+1 with `ref_span = 10`
  extends `end` to p+10; subsequent items at p+5 (`ref_span = 1`)
  do **not** shrink `end`.
- **Single record's `end == pos + ref_span - 1`.** Pin the
  inclusive-1-based convention so a future refactor that
  accidentally flips to exclusive is caught.
- **Empty `per_sample` (all-`None`).** A pathological
  `PerPositionPileups` with every slot `None` must not panic;
  `max_ref_span` falls back to `1` and the item forms a
  single-position group (which then gets filtered as REF-only).

`OverlappingVarGroup::span()` and the basic getters get smoke-test
coverage; they are thin enough that the integration tests above
exercise them implicitly.

## Validation

Inside the dev container (`./scripts/dev.sh`):

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo build --examples`
- `cargo build --benches`

No new benchmark is added in this plan. Benchmarking the grouper in
isolation would over-fit a synthetic generator. The natural place
to bench it is end-to-end once Stage 5 exists and a realistic
`.psp` → group → merged-record pipeline can be driven from real
data.

## Assumptions / silent choices

- **`OverlappingVarGroup` owns `Vec<PerPositionPileups>` directly**,
  not a flattened `Vec<(sample_idx, PileupRecord, pos)>`. Stage 5
  consumes per-position pileups (allele unification works per
  reference position), so preserving the upstream shape avoids a
  conversion. The cost is one `Vec<Option<PileupRecord>>` per
  group record — bounded by N — which is paid one way or another.
- **The grouper is single-threaded.** Cross-sample parallelism is
  already paid at the merger level, and per-group parallelism
  happens downstream at Stage 5. Adding a parallel layer at the
  grouper would compete with both for cores without unlocking new
  work; the algorithm is also inherently sequential at the
  boundary-detection step
  ([calling_pipeline_architecture.md:1028-1039](../specs/calling_pipeline_architecture.md#L1028-L1039)).
- **Errors latch, partial groups are dropped.** A reader error
  mid-group discards the in-progress group rather than emitting it
  truncated. Stage 5 has no way to distinguish a truncated group
  from a complete one — partial output would be a correctness
  hazard. Matches the merger's terminate-on-first-error contract.
- **Pure-REF-only groups are dropped silently.** See §"Pure-REF-only
  groups" above for the rationale. This is the one place the
  grouper makes a policy call rather than a mechanical one; flagged
  explicitly so a future reviewer sees it.
- **`max_ref_span` per position falls back to `1` on an all-`None`
  upstream item.** Defensive; the merger does not emit such items
  in production. Pinned by a test so a future reader knows what
  the fallback means.
- **Grouper does not own the upstream iterator's source.** As with
  the merger, the caller holds the `PspReader`s and the merger
  alongside the grouper. Bundling is deferred until a real
  path-based opener exists (same call as the merger plan's
  "Path-based opener helper" follow-up).

## Risks

- **Stage 3 (DUST) error type compatibility.** The grouper is
  currently hard-coded to `MergerError` as the upstream `E`. When
  DUST lands and it needs a richer error type (e.g. for
  reference-fetch failures), either DUST wraps everything into
  `MergerError` at its boundary, or we generalise the grouper on
  `E` here. Both are minor changes; not pre-built.
- **Group-size pathology under contrived inputs.** Repeated
  re-extensions could in principle grow a group arbitrarily large.
  In practice the upstream `MAX_RECORD_SPAN` bounds any single
  record's reach to 5000 bp, and observed transitive chains in
  the old grouper rarely exceeded a few hundred bp. Addressed in
  this plan via `GrouperConfig::max_var_group_span` (see
  §"Configurable parameters") with a hard-error variant
  (`GrouperError::VarGroupTooWide`) that surfaces the offending
  locus and the cap. Default is 10000 bp; users hitting it on
  real data raise the cap explicitly via the config rather than
  the grouper deciding to drop on their behalf.
- **`PerPositionPileups: Clone` cost in tests.** Same caveat as the
  merger plan — `Clone` is only needed for `PartialEq` assertions
  in tests. Production paths move groups out of the iterator.
  Confirm by inspection that no `clone()` slips into production
  code at code-review time.

## Deferred cleanup: remove the gVCF path after Stage 5

The old code is **not** removed when this plan ships. The new
grouper reproduces the algorithmic ideas from
[src/variant_grouping.rs](../../src/variant_grouping.rs)
immediately, but the test suites in
[tests/genotype_merging_test.rs](../../tests/genotype_merging_test.rs)
and
[tests/variant_group_test.rs](../../tests/variant_group_test.rs)
encode merging edge cases (overlapping deletions across samples,
within-sample compound indels, phase-chain preservation across
~25 specific scenarios — catalogued in §"Test-encoded edge cases
Stage 5 has to reproduce" above) that Stage 5 has to handle
correctly on the new pipeline. Deleting those tests before Stage
5 lands would lose the truth set Stage 5 needs to be validated
against. The tradeoff is dead code lingering in the tree for a
while; that is cheaper than re-deriving the catalogue from real
data later.

**Gating condition.** Removal proceeds only once **all** of the
following are true:

1. Stage 5 (per-group processing) has shipped with its own
   integration tests.
2. The porting checklist demanded at the end of
   §"Test-encoded edge cases" above is complete — every
   merging scenario the old suite pins down is either covered
   by a Stage 5 test on `OverlappingVarGroup` input, or has an
   explicit waiver recorded in the Stage 5 plan.
3. The new cohort CLI (or whatever replaces the gVCF subcommand
   in `src/main.rs`) is in place, so deleting `pipeline.rs`
   does not leave the binary without a non-Stage-1 entry point.

Until all three hold, the old code stays.

**Removal scope** (one commit, once gating conditions are met):

- `src/variant_grouping.rs` — delete the file.
- `src/lib.rs` — remove `pub mod variant_grouping;`.
- `src/pipeline.rs` — the gVCF → VCF entry point uses
  `VarGroupIterator`. The whole pipeline file is part of the
  gVCF path and goes away with the grouper.
- `src/genotype_merging.rs` — consumes `OverlappingVarGroup` and
  `VarIteratorInfo` from the old grouper. Belongs to the same
  gVCF path; delete the file.
- `tests/genotype_merging_test.rs`, `tests/variant_group_test.rs`
  — integration tests for the gVCF path; delete only after the
  porting checklist confirms every scenario is reproduced.
- `src/main.rs` — drop the gVCF CLI subcommand that wires
  `pipeline.rs`. Confirm scope when the commit is prepared.
- Any other artefacts that the gVCF path drags in (e.g.
  `gvcf_parser.rs` if nothing else consumes it). Audit by
  `cargo build` after deleting the modules above and removing
  whatever the compiler points at.

This is **not** a refactor when it finally happens — no
behavioural code from the old files gets ported into the new
tree at cleanup time. The algorithmic lessons are reproduced
now in the new grouper; the test-encoded edge cases are
reproduced as Stage 5 lands; the rest is dead code on a
deprecated execution path. Keeping it around past those
milestones invites someone to read it as authoritative.

## Out-of-scope follow-ups

- **Path-based opener bundling merger + grouper.** Same pattern
  as the merger plan's
  [open-paths helper](multi_way_per_position_iterator.md#L398-L412)
  follow-up. Wait until a real caller (the cohort CLI subcommand)
  has shipped and the ergonomic shape is obvious.
- **`--emit-ref-only-groups` CLI flag.** Toggle the pure-REF
  filtering policy off for diagnostic runs that want exhaustive
  output. Implement when the cohort CLI lands.
- **`max_group_span` cap and the corresponding logged skip
  category.** Add if a real dataset demonstrates pathological
  group growth. Until then, the upstream `MAX_RECORD_SPAN` is the
  effective bound and adding a second cap is speculative.
- **DUST adaptor (Stage 3).** Slot between merger and grouper as a
  separate `Iterator::filter`-shaped struct. Out of scope for this
  plan but explicitly enabled by it: the grouper's generic-iterator
  shape is what lets DUST drop in later without changing this code.
- **Error generalisation `VariantGrouper<I, E>`.** Lift the
  hard-coded `MergerError` to a generic `E` if a Stage 3 adaptor
  needs a richer error type. Don't pre-build.

## File touch list

**This plan's commit** (the new grouper, no removals):

- `src/cohort/mod.rs` — add `pub mod variant_grouping;` alongside
  the existing `per_position_merger` declaration.
- `src/cohort/variant_grouping.rs` — new file: `VariantGrouper`,
  `OverlappingVarGroup`, `GrouperConfig` (with the
  `DEFAULT_MAX_VAR_GROUP_SPAN` const), `GrouperError`, the
  `has_variant_observation` and `max_ref_span` helpers, fixture
  builders, full `#[cfg(test)]` module.
- (No changes to `src/lib.rs` — `pub mod cohort;` is already
  exported.)

No changes to existing types. The grouper consumes
`PerPositionPileups` and `MergerError` from the merger module
unchanged.

**Deferred follow-up commit** (the gVCF cleanup, scope and
gating conditions detailed in §"Deferred cleanup: remove the
gVCF path after Stage 5"): the old grouper, its merging
consumer, and their integration tests come out only after Stage
5 has shipped with equivalent test coverage for the merging
scenarios the old suite encodes. Not done as part of this
plan's commit.
