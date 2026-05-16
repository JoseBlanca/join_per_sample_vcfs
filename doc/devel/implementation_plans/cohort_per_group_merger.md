# Per-group merger (Stage 5 — allele unification + per-sample likelihood reconstruction)

Proposal date: 2026-05-16.

## Domain intent

Stage 5 is the algorithmic heart of the cohort side. It consumes the
`OverlappingVarGroup` stream produced by Stage 4 and emits one
merged multi-sample record per group, ready for Stage 6's posterior
EM:

```
… ─► VariantGrouper ─► PerGroupMerger ─► MergedRecord stream ─► Stage 6 (posterior EM)
                          │
                          └─ rayon-parallel across groups
```

Stage 5's contract from
[calling_pipeline_architecture.md:1052-1814](../specs/calling_pipeline_architecture.md#L1052-L1814):

1. **Unify alleles** across the samples in the group into a single
   merged allele set. Compound haplotype alleles only enter if at
   least one sample has chain-id evidence linking the constituents.
2. **Reconstruct per-sample likelihoods** against the merged allele
   set using freebayes' closed-form formula from the five per-allele
   scalars stored in each `.psp` record.
3. **Reconcile compound alleles** sample-by-sample: chain-evident
   samples evaluate the compound directly; chain-broken samples
   fall back to constituents-independent likelihoods with a
   cohort-prior tag (`CA = 1`) for Stage 6.
4. **Emit a merged record** containing the merged allele set, the
   per-sample scalar table (so Stage 6 can recompute likelihoods
   when contamination correction is active), and the per-sample,
   per-genotype log-likelihood vector.

This is the largest and most subtle module on the cohort side. It
is also the natural home of every correctness obligation that the
old gVCF-era merger's ~25-test catalogue
([tests/genotype_merging_test.rs](../../tests/genotype_merging_test.rs))
already pinned down — these are reproduced as Stage 5 tests on
`OverlappingVarGroup` input (see §"Test strategy").

## Why now

Stage 4 has shipped (commit `6f6cd46`,
[src/cohort/variant_grouping.rs](../../src/cohort/variant_grouping.rs)),
so its output type — `OverlappingVarGroup` carrying
`Vec<PerPositionPileups>` — is fixed and the test fixtures are
proven. Stage 5 is the bottleneck for completing the cohort
pipeline; Stage 6 (posterior EM) consumes Stage 5's output and is
much smaller in scope.

Building Stage 5 also unblocks the **deferred cleanup of the gVCF
path**: the architecture allows deletion of
[src/genotype_merging.rs](../../src/genotype_merging.rs),
[src/pipeline.rs](../../src/pipeline.rs),
[src/variant_grouping.rs](../../src/variant_grouping.rs),
[tests/genotype_merging_test.rs](../../tests/genotype_merging_test.rs)
and
[tests/variant_group_test.rs](../../tests/variant_group_test.rs)
only after Stage 5 reproduces the merging-edge-case catalogue
encoded in those tests
([cohort_variant_grouping.md §"Deferred cleanup"](cohort_variant_grouping.md#L635-L692)).
The porting checklist that gates that cleanup lives inside the
test-strategy section of this plan.

## What's already in place

- **Stage 4 grouper** —
  [src/cohort/variant_grouping.rs](../../src/cohort/variant_grouping.rs).
  Emits `Result<OverlappingVarGroup, GrouperError>` whose
  `records: Vec<PerPositionPileups>` carry the full per-sample
  slot vectors. Pure-REF groups are already filtered upstream by
  the seed-time check; Stage 5 never sees a group whose every
  position is pure-REF across every sample.
- **`PerPositionPileups`** —
  [src/cohort/per_position_merger.rs:36-46](../../src/cohort/per_position_merger.rs#L36-L46).
  `(chrom_id, pos, per_sample: Vec<Option<PileupRecord>>)`.
- **`PileupRecord` / `AlleleObservation` / `AlleleSupportStats`** —
  [src/per_sample_caller/pileup/mod.rs:340-455](../../src/per_sample_caller/pileup/mod.rs#L340-L455).
  Each allele observation already carries:
  - `seq: Vec<u8>` — the allele bytes against this record's REF
    span (anchor-base convention for indels).
  - `support`: the five freebayes scalars
    (`num_obs`, `q_sum = Σ max(ln_BQ, ln_MQ)`, `fwd`,
    `placed_left`, `placed_start`).
  - `chain_ids: Vec<ChainId>` — unique-per-`.psp` `u64`
    identifiers naming the reads that contributed to this allele.
    Sorted, deduplicated.
- **Chain id semantics** —
  [doc/devel/specs/phase_chain.md](../specs/phase_chain.md). Chain
  ids are unique per `.psp` file and **never comparable across
  samples**. Stage 5 only ever intersects chain ids inside a
  single sample's slot vector.
- **Walker compound-extension precedent** —
  [src/per_sample_caller/pileup/open_record.rs:600-787](../../src/per_sample_caller/pileup/open_record.rs#L600-L787).
  The walker already builds compound alleles *within a record's
  footprint*: when a read's events span the record (e.g. SNP at
  101 + deletion at 102–105 on the same read), the read folds
  into a single `AlleleObservation` whose `seq` encodes the
  compound. Stage 5's job is the **cross-record** compound — a
  compound whose constituents live in separate records but whose
  reads carry chain ids that intersect.

## Lessons to extract from the old gVCF merger before its deletion

The old gVCF-era merger lives in
[src/genotype_merging.rs](../../src/genotype_merging.rs) and its
catalogue of edge cases in
[tests/genotype_merging_test.rs](../../tests/genotype_merging_test.rs).
The cleanup of those files is gated by this plan (per
[cohort_variant_grouping.md §"Deferred cleanup"](cohort_variant_grouping.md#L635-L692));
nothing is deleted until Stage 5 reproduces every lesson on the new
input shape. Two distinct categories of value live in the old
code:

### Algorithmic ideas reproduced here

1. **Single-pass per-haplotype allele assembly with deletion-state
   tracking.** The old `create_variant_for_region`
   ([genotype_merging.rs:23-245](../../src/genotype_merging.rs#L23-L245))
   walks the group's variants in order and, for each sample,
   maintains a per-haplotype `positions_left_in_del` counter
   plus per-haplotype allele fragments. When a deletion is active
   in a sample at a position, that haplotype contributes an empty
   string (deleted base) rather than the REF. The Stage 5 allele
   unification step needs an analogous projection: every
   merged-allele's seq is a concatenation of per-position
   contributions across the group's reference span.
2. **Phase-broken-haplotype detection (`first_het_seen` +
   `phase_broken_since_het`).** When phase is broken between
   hets, the old merger emits a missing call (`-1`) for the
   affected haplotype rather than guessing. The new pipeline
   uses chain ids instead of VCF phase bits, but the
   *consequence* is the same: a chain-broken sample at a
   compound allele cannot have its per-haplotype assignment
   reconstructed from per-position observations alone; Stage 5
   falls back to the cohort-prior path (CA = 1) at that compound
   for that sample.
3. **Output-phase rule: `false` if first position unphased or any
   het unphased.** The new pipeline does not emit per-sample
   phase fields in v1 (compound alleles are handled at the
   allele-set level, not as VCF phase). But the implicit
   invariant — "any sample whose phase across the group cannot
   be unambiguously reconstructed must surface that fact" —
   carries forward as the CA flag on the compound, plus
   downstream consumers' ability to filter on CA.
4. **`is_first` initialisation pattern.** Per-sample
   accumulators start in a sentinel state and initialise on the
   first variant the sample participates in
   ([genotype_merging.rs:70-75](../../src/genotype_merging.rs#L70-L75)).
   The new plan keeps the same pattern — a sample with no record
   at any position of the group contributes no scalars and
   appears as a "no-data" entry in the merged record.
5. **Sample offset / multi-iterator indexing.** The old code
   carries `var_iter_idx` because the gVCF path opens N
   independent VCF iterators. The new pipeline merges one stream
   upstream of Stage 4; Stage 5 indexes samples directly by
   their position in `PerPositionPileups::per_sample`. The
   complexity collapses — Stage 5 does not maintain
   per-iterator offsets.
6. **Synthetic-PL fallback for inputs without PL data.** The
   old merger
   ([genotype_merging.rs:289-298](../../src/genotype_merging.rs#L289-L298))
   manufactures PLs from GT calls (99 % confidence on the called
   genotype) when input lacks PLs. The new pipeline **never**
   needs this fallback: Stage 5 reconstructs likelihoods from
   the five per-allele scalars exactly
   ([calling_pipeline_architecture.md:1481-1484](../specs/calling_pipeline_architecture.md#L1481-L1484)),
   which is one of the core wins of the `.psp` format over
   gVCF/PL input.
7. **EM convergence as a separate stage.** The old merger does
   per-variant EM after each merge
   ([genotype_merging.rs:277-363](../../src/genotype_merging.rs#L277-L363)).
   The new pipeline separates that into Stage 6, which runs EM
   cohort-wide across all merged records. Stage 5 stops at the
   per-sample-per-genotype likelihood; Stage 6 owns the
   posterior.

### Test-encoded edge cases that drive the test plan

The old test suite catalogues 26 distinct merging scenarios that
Stage 5 must handle correctly. Categorised:

**A. Basic merging (4 tests).** Simple SNPs, simple insertions,
simple deletions, deletion-of-length-2 — the canonical
single-variant cases that any allele-unification implementation
must get right. Reproduced as Stage 5 tests on synthetic
`OverlappingVarGroup` fixtures.

**B. Allele unification across samples (4 tests).** Overlapping
deletions across samples that must unify into a canonical form,
deletion-SNP overlap, allele-merge-with-non-ref-filtered,
missing-allele-in-merged-haplotype. These are the cohort-level
correctness obligations of allele unification itself.

**C. Within-sample compound alleles (2 tests).**
`test_two_deletions_in_same_sample` and
`test_two_overlapping_deletions_in_same_sample` — exercise the
case where a single sample carries multiple indels in the same
group. The new pipeline's chain ids replace the old VCF phase
bits but pin the same biological scenario: two indels co-occurring
on the same haplotype in one sample must produce a compound
allele, not two independent ones.

**D. Phase-chain preservation (11 tests).** The `test_phase_*`
family — the most subtle slice of the old catalogue. Each test
pins a specific way that VCF phase information across hets in a
group can be preserved, broken, or inferred. The new chain-id
mechanism reproduces the same biological scenarios using a
different mechanism, so the test bodies don't port verbatim —
but every scenario maps to a Stage 5 case (see §"Test strategy
§Phase-equivalent chain-id scenarios" below for the mapping).

**E. Allele-set hygiene (1 test).**
`test_allele_merge_with_non_ref_filtered` — symbolic alleles
like `<NON_REF>` are filtered out before merging. Stage 5 does
not have `<NON_REF>` (the freebayes-scalar reconstruction
generalises it), but the underlying invariant — "the merged
allele set contains only concrete sequence alleles, never
symbolic placeholders" — applies.

**F. Group-shape correctness (3 tests).**
`test_group_merging_creates_correct_number_of_vars`,
`test_non_variant_vars_are_removed`,
`test_analyze_single_group_from_merged_bin`. Mostly Stage 4
concerns (already covered by the new grouper's tests), but Stage
5 must also drop merged records that turn out to have only one
allele after unification (every sample is hom-ref).

### Mapping from old tests to Stage 5 test scenarios

The test-strategy section below has a one-to-one table: for each
old test, the Stage 5 equivalent on `OverlappingVarGroup` input
and what it asserts. That table is the **porting checklist** that
gates the gVCF cleanup commit (see §"Deferred cleanup:
completing the gVCF-path removal").

## Algorithmic alternatives considered

Stage 5's job — turn per-sample observation evidence at overlapping
positions into a single cohort-level record with per-sample
likelihoods — has been solved before by GATK, by freebayes, and by
the old gVCF merger this project ships. Each took a different
approach. Below we describe the four candidates and the rationale
for picking ours.

### Algorithm 1 — GATK GenotypeGVCFs (PL projection with `<NON_REF>`)

**Shape.** Per-sample input is a gVCF whose every site carries an
explicit PL array plus a symbolic `<NON_REF>` allele whose PL
represents "any other allele not enumerated here." At join time,
the cohort merger
([ReferenceConfidenceVariantContextMerger.java](../../gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/ReferenceConfidenceVariantContextMerger.java))
unifies allele sets across samples by remapping each sample's
local alleles to the cohort's union, filling missing-allele slots
with the `<NON_REF>` PL. Genotype posteriors are then computed by
EM
([GenotypeGVCFsEngine.regenotypeVC](../../gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/GenotypeGVCFsEngine.java)).

**Strengths.**

- Per-sample input is self-contained: every site has explicit
  PLs for every allele the sample sees, plus a fallback PL for
  novel alleles. The merger's job is largely bookkeeping.
- Cohort re-callable: adding a new sample doesn't re-process
  existing samples.
- Battle-tested at population scale (gnomAD, UK Biobank).

**Weaknesses for our pipeline.**

- **Per-site PLs are lossy.** The `<NON_REF>` PL is a single
  scalar that has to stand in for any conceivable novel allele,
  including compound haplotypes that span multiple positions.
  When the cohort introduces a compound that one sample
  *should* be evaluated against, GATK cannot reconstruct the
  per-read likelihood for that compound — the per-read evidence
  was already collapsed into per-site PLs upstream.
- **No per-read evidence at join time.** Compounds spanning
  positions are decomposed into per-site records in the gVCF;
  the merger cannot detect "these two alleles co-occur on the
  same haplotype in this sample" because that information was
  thrown away when the gVCF was written. GATK papers over this
  with downstream tools (HaplotypeCaller does some local
  re-assembly, but cross-sample compound calling still requires
  external phasing).
- **Storage overhead grows with allele count.** A diploid site
  with 4 alleles needs 10 PL entries per sample (`C(4+2-1, 2) =
  10`); 6 alleles is 21; 10 alleles is 55. The `<NON_REF>`
  mechanism doesn't shrink this — it still needs PLs against
  every site-local allele plus one against `<NON_REF>`.

**Verdict.** Rejected for our pipeline. The two structural
problems — lossy `<NON_REF>` and absent per-read evidence at
join — together preclude correct compound-haplotype handling at
cohort scale. The architecture spec already commits to a different
shape ([calling_pipeline_architecture.md:476-494](../specs/calling_pipeline_architecture.md#L476-L494)).

### Algorithm 2 — freebayes (combo search over joint BAMs)

**Shape.** Read every sample's BAM into memory at the same
position. Enumerate (or search) joint genotype combos across all
samples; for each combo compute likelihood × prior; marginalise
back to per-sample posteriors. The per-read likelihood uses
`max(ln_BQ, ln_MQ)` plus a multinomial allele-balance term plus a
read-dependence-dampening correction (RDF, default 0.9).

**Strengths.**

- Per-read evidence is available at join time, so haplotype-level
  alleles compose naturally (a SNP + downstream deletion on the
  same read is a single allele end-to-end).
- Implicit integration over `p` (allele frequency) via combo
  marginalisation gives well-calibrated posteriors at small
  cohorts.
- No `<NON_REF>` machinery needed.

**Weaknesses for our pipeline.**

- **Does not scale.** Joint BAMs grow with `N_samples ×
  read_depth × locus_span`. Beyond a few hundred samples it
  doesn't fit in memory; combo search itself grows much faster
  than `N`.
- **No reusable per-sample artefact.** Adding a new sample to a
  cohort means rerunning the joint call across every sample's
  BAM from scratch.
- **Combo enumeration is exponential.** Freebayes prunes the
  search; the pruning has its own failure modes at sites with
  many alleles.

**Verdict.** Rejected at the cohort level. The pipeline borrows
freebayes' per-read likelihood model (`max(ln_BQ, ln_MQ)`, the
multinomial, the five per-allele scalars as sufficient statistics)
but drops the BAM-in-memory and combo-search aspects. The
architecture spec settles this at lines
[1730-1753](../specs/calling_pipeline_architecture.md#L1730-L1753)
("Why EM rather than freebayes' combo search").

### Algorithm 3 — old gVCF merger (string-concat alleles, synthetic PLs)

**Shape.** Per-sample input is a gVCF; the grouper produces an
`OverlappingVarGroup` of `Variant` records;
`create_variant_for_region` walks the group, maintains per-sample
per-haplotype allele fragments + deletion-state counters, and
emits a unified `Variant` with a string-merged allele set. PLs are
either copied from input or synthesised from GT (99 %-confidence
on the called genotype). EM runs per-variant on the resulting PLs
([genotype_merging.rs:277-363](../../src/genotype_merging.rs#L277-L363)).

**Strengths.**

- Cohort-scaleable: linear in sample count, no joint BAMs in
  memory.
- Per-haplotype tracking via VCF phase bits handles the simplest
  compound cases.
- Existing test catalogue pins ~25 merging scenarios.

**Weaknesses for our pipeline.**

- **Input is gVCFs, not `.psp`.** The whole module is built on
  `Variant` records carrying GT/PL/AD — not on the per-allele
  scalar shape the new pipeline emits.
- **VCF phase is coarse.** Phase bits are per-record, not
  per-read; a sample whose reads phase a compound but whose
  per-record phase was lost (typical in gVCF input from a
  caller that doesn't preserve fine-grained phase) cannot have
  its compound reconstructed.
- **Synthetic PLs are an information-discarding shortcut.** A
  GT call collapses the per-read evidence into a single
  "99 %-confidence" PL, which Stage 6 then has to undo.
- **EM-per-variant fragments the cohort calibration.** The old
  merger runs EM separately per group; cross-group information
  (e.g. allele frequency at neighbouring loci) is unavailable
  during that EM.

**Verdict.** Rejected as an input shape but retained as a
*correctness oracle*: the scenarios its tests pin down are the
ones Stage 5 must reproduce, even though the implementation is
completely different. The data shape (`.psp` per-allele scalars +
chain ids) makes Algorithm 4 below structurally cleaner.

### Algorithm 4 (chosen) — scalar reconstruction with chain-id-anchored compounds

**Shape.** Per-sample input is the `.psp` records inside the
`OverlappingVarGroup`. Each `AlleleObservation` carries the five
freebayes-canonical scalars plus a `chain_ids` list naming the
reads that contributed. Stage 5:

1. Builds the merged allele set by projecting every per-sample
   per-position allele onto the group's reference span.
2. Proposes compound haplotype alleles only when at least one
   sample has chain ids linking the constituents across
   records.
3. For each sample, projects its per-position scalars onto the
   merged allele set: scalars for alleles the sample observed
   carry through; alleles the sample never observed get
   `(count = 0, S_a = 0, ...)`.
4. For each sample, computes the closed-form likelihood
   `L(G) = Σ_{a∉G} S_a + multinomial(allele_probs(G), obs_counts_in_G)`
   for every candidate genotype `G` in the merged set.
5. For chain-broken samples at chain-anchored compounds,
   substitutes the constituents-independent likelihood (which
   does not bias the sample away from the compound) and flags
   the sample `CA = 1` so Stage 6 knows to use the
   cohort-derived compound frequency `f_C` as the prior signal
   for that sample at that compound.
6. Emits a `MergedRecord` carrying the merged allele set, the
   per-sample scalar table, and the per-sample, per-genotype
   log-likelihood vector under `c_s = 0`. Stage 6 augments with
   contamination correction if enabled.

**Strengths.**

- **Per-read evidence is preserved through the `.psp` chain
  ids.** A compound proposed by sample A is *anchored* by
  chain-id evidence in A. Other samples (B, C, ...) may evaluate
  the compound differently based on whether their own chain ids
  link the constituents, with the chain-broken fallback handling
  the case where they don't.
- **`<NON_REF>`-equivalent behaviour is automatic.** A sample
  that never observed allele `x` has `S_x = 0` and `count_x =
  0`; those terms simply drop out of the likelihood formula
  ([calling_pipeline_architecture.md:1491-1508](../specs/calling_pipeline_architecture.md#L1491-L1508)).
  No symbolic-allele machinery needed.
- **Cohort-scaleable.** Linear in sample count, no joint BAMs,
  no exponential combo search.
- **Cohort re-callable.** Adding a new sample is one more `.psp`;
  existing samples' `.psp` files are untouched.
- **Compound-handling is principled.** The chain-anchor rule
  prevents the failure mode where a single-sample compound
  proposal silently miscalls every chain-broken carrier; the CA
  flag carries the chain-broken signal forward to Stage 6 so
  the cohort prior carries the burden.
- **Bit-for-bit reproducibility of freebayes' likelihood.** The
  reconstruction from scalars is exact
  ([calling_pipeline_architecture.md:1481-1484](../specs/calling_pipeline_architecture.md#L1481-L1484)),
  not an approximation.

**Weaknesses.**

- **Compound-allele detection requires looking at every chain
  id across positions of the group.** This is `O(records ×
  alleles × chain_ids)` per sample. In practice each
  `AlleleObservation`'s `chain_ids` is short (one or a few
  reads supporting the allele), so the constant is small, but
  the algorithm has to be careful about the asymptotic on
  high-coverage adversarial inputs.
- **Chain-broken samples at chain-anchored compounds need a
  separate code path** (the cohort-prior fallback). This is the
  most subtle piece of the implementation and has the highest
  test-coverage requirement.
- **Allele count can grow with compound proposals.** A group
  with two chain-anchored compounds plus their constituents can
  have 5 alleles where naive merging would have 3. The
  `max_alleles` cap (see §"Configurable parameters") guards
  against this.

**Decision.** Algorithm 4 is the choice. It is what the
architecture spec already commits to
([calling_pipeline_architecture.md:1052-1814](../specs/calling_pipeline_architecture.md#L1052-L1814));
this plan operationalises it. The decision was driven by the
combination of:

1. Per-read evidence preserved in chain ids without the cost of
   per-read storage in the `.psp` (the five scalars suffice for
   likelihood reconstruction;
   [freebayes_posterior_gt_probs.md §Compactness implications](../specs/freebayes_posterior_gt_probs.md)).
2. Cohort scalability matching GATK's gVCF topology but with
   freebayes' algorithmically-correct likelihood model.
3. Compound-handling that avoids both GATK's
   per-site-decomposition failure and freebayes' per-cohort BAM
   re-read.

The combination is novel relative to either upstream caller and
the project's existing design docs spend significant space
defending it
([gatk_vs_freebayes_comparison.md](../specs/gatk_vs_freebayes_comparison.md),
[phase_chain.md](../specs/phase_chain.md)). Stage 5 is where
that combination lands in code.

### Side comparison: the analogous step in our own walker

Our per-sample pileup walker
([src/per_sample_caller/pileup/open_record.rs:600-787](../../src/per_sample_caller/pileup/open_record.rs#L600-L787))
already does a *within-record* version of allele unification: when
a single read's CIGAR events span a record's footprint, the
walker folds the read into one `AlleleObservation` whose seq
encodes the compound. Chain ids are added to that allele's
`chain_ids` list. Stage 5's job is the cross-record version of
the same machinery, with the additional constraint that the
compound must be **chain-anchored** in at least one sample (the
walker's within-record compounds are automatically anchored
because they live in one read).

Patterns from the walker that carry forward:

- **Allele-key = seq bytes.** The walker uses byte-vector equality
  to dedupe alleles within a record; Stage 5 uses the same after
  projecting each per-position allele onto the group's reference
  span.
- **Sort + dedupe of chain ids.** The walker maintains
  `chain_ids` as a sorted, deduplicated `Vec<u64>`
  ([open_record.rs:918-922](../../src/per_sample_caller/pileup/open_record.rs#L918-L922)).
  Stage 5 reuses the same shape for compound-allele chain id
  sets and for the chain-id-intersection check.
- **Subtract-on-rebucket pattern.** The walker subtracts a
  prior contribution before re-adding it when a record widens
  ([open_record.rs:760-762](../../src/per_sample_caller/pileup/open_record.rs#L760-L762)).
  Stage 5 doesn't widen records but uses the same subtract-then-add
  pattern when reassigning observations to compound buckets
  during chain-id resolution.

## Algorithm: allele unification + scalar projection + likelihood reconstruction

### Step 1 — allele unification across the group

The group spans the reference range `[group.start, group.end]`.
Each record in the group sits at position `pos` with REF span
`ref_span` covering `[pos, pos + ref_span - 1]`. Each record's
alleles are encoded against that local REF span using the
walker's anchor convention.

**Project every per-sample allele onto the group's reference
span.** For each `(sample_idx, record, allele_obs)` triple in the
group, build a `MergedAlleleKey` whose bytes encode the haplotype
across the full `[group.start, group.end]` range:

- Positions outside the record's footprint contribute the group's
  reference sequence (fetched from the FASTA via the same
  `RefSeqFetcher` Stage 1 uses).
- Positions inside the record's footprint contribute the allele's
  bytes, anchor-aligned per the walker's convention.

The result is a byte sequence of length `group.end - group.start + 1`
in the simple SNP-only case, or longer for compound or insertion-
bearing alleles (the seq length encodes the indel shape).

**Deduplicate the projected alleles** across samples by byte
equality. The merged allele set is the deduplicated list, with
REF (the unaltered group reference sequence) always at index 0.

**Carry the source-sample-and-allele-observation provenance**
through the projection: each merged allele records which
`(sample_idx, record_idx, local_allele_idx)` triples projected to
it, so Step 3 can sum scalars correctly.

### Step 2 — chain-anchor entry for cross-record compounds

A **compound allele** is one whose support, in some sample, is
distributed across two or more records (i.e. its non-REF bytes
come from more than one `(record_idx, local_allele_idx)` pair).
The architecture spec at
[lines 1062-1081](../specs/calling_pipeline_architecture.md#L1062-L1081)
requires that compound alleles enter the merged set **only if at
least one sample has chain-id evidence** linking all of the
compound's constituents on the same haplotype.

**Chain-id intersection check.** For each candidate compound `C`
whose constituents live at records `r_1, r_2, ..., r_k` for some
sample `s`, look at the chain id list of `s`'s allele observation
at each constituent's record (i.e. the allele observation in `s`
that contributed to this compound's projection). The compound is
**chain-anchored** by `s` if:

```
chain_ids(s, r_1, allele_at_r_1)
  ∩ chain_ids(s, r_2, allele_at_r_2)
  ∩ ...
  ∩ chain_ids(s, r_k, allele_at_r_k)
≠ ∅
```

— i.e. the intersection is non-empty, meaning at least one read
in `s` carries all constituents on the same molecule.

**Compound entry rule.** A candidate compound enters the merged
allele set iff at least one sample chain-anchors it. Compounds
that no sample chain-anchors are rejected; their constituents
appear in the merged allele set as independent per-position
alleles (this is the natural outcome of Step 1 — projection of
each constituent on its own already produced separate merged
alleles).

**Within-record "compounds" are always anchored.** When a single
read's events span a single walker-emitted record (e.g. SNP at
101 + deletion 102–105 on one read, and the walker's record sits
at pos=101 with `ref_span=5`), the resulting `AlleleObservation`
encodes the compound natively as its `seq` bytes. Its `chain_ids`
list is populated by the walker. Stage 5 never has to
intersect-across-records for these — they appear as single
records in the group and project cleanly in Step 1. The
chain-anchor check is only needed for *cross-record* compounds.

### Step 3 — per-sample scalar projection onto the merged allele set

Build a per-sample, per-merged-allele scalar table. For each
sample `s` and each merged allele `M`:

- `count_s(M)` = sum of `support.num_obs` across every
  `(record_idx, local_allele_idx)` that projected to `M` in
  sample `s`.
- `S_s(M)` = sum of `support.q_sum` across the same.
- Similarly for `fwd`, `placed_left`, `placed_start`.

For **single-position alleles** (constituent of no compound or
fully decomposable to one record) the sum is over exactly one
record. For **cross-record compounds** in chain-evident samples:

- `count_s(C) = |chain_ids intersection|` — the number of reads
  in this sample that carry the compound. Each intersecting
  chain id contributes once.
- `S_s(C) = sum over reads in the chain-id intersection of
  max(ln_BQ, ln_MQ)`. **Subtlety:** the `q_sum` scalar in
  `AlleleSupportStats` is a per-allele aggregate across all reads
  supporting that allele in that record — it doesn't split out
  per-read contributions. The architecture spec at
  [lines 1352-1377](../specs/calling_pipeline_architecture.md#L1352-L1377)
  describes the same constraint at the contamination layer and
  resolves it via the **homogeneous-quality approximation**: treat
  every read supporting an allele as sharing the mean quality
  `q_sum / num_obs`. The same approximation applies here for
  cross-record compound `q_sum` reconstruction:

  ```
  S_s(C) ≈ |chain_id ∩| × min_over_constituents( q_sum_at_record_i / num_obs_at_record_i )
  ```

  with `min_over_constituents` because a read's per-position
  quality at the constituent positions is each subject to its
  own BQ/MQ floor — under the indel BQ-proxy convention the
  effective per-read indel quality is the *minimum* in a window,
  and the cross-record compound's effective quality cannot
  exceed any single constituent's
  ([calling_pipeline_architecture.md:305-371](../specs/calling_pipeline_architecture.md#L305-L371)).

- `fwd_s(C)`, `placed_left_s(C)`, `placed_start_s(C)` —
  reconstructed by the same homogeneous approximation: scale the
  per-allele bias counts of any constituent by `|chain_id ∩| /
  num_obs_at_that_constituent`. The bias scalars feed Stage 6
  observation-bias priors, where the approximation error tends to
  be in the noise (the bias priors are themselves order-of-magnitude
  signals, not precision instruments).

**Subtract compound contribution from constituents in chain-evident
samples.** When a sample chain-anchors a compound `C`, the reads
that support `C` also support each constituent's local allele
*individually* in the underlying records (the chain ids point to
those local allele observations). To avoid double-counting, **the
compound's scalars are subtracted from each constituent's
scalars in this sample**: a read that supports both should
contribute to `count_s(C)`, not to `count_s(constituent_1)` and
`count_s(constituent_2)`. After the subtraction, a constituent's
scalars in a chain-evident sample reflect only the reads that
support it *without* supporting the compound (e.g. reads spanning
only the SNP but not the deletion).

**Chain-broken samples don't subtract.** A chain-broken sample
has no compound observations to subtract; the constituent
scalars stand. The chain-broken sample's compound likelihood is
computed via the fallback path in Step 5, which does not consume
`count_s(C)` or `S_s(C)`.

### Step 4 — chain-evident vs chain-broken sample classification

Per sample, per chain-anchored compound: classify the sample as
**chain-evident** (the chain-id intersection across the
compound's constituents in this sample is non-empty) or
**chain-broken** (the intersection is empty — typically because
no read or pair in this sample spans all constituents). Record a
per-sample bit `chain_evident[s][C]` for every compound in the
merged set.

A sample is chain-evident at a non-compound allele trivially —
single-position alleles have no chain-intersection question. The
`chain_evident` table is only meaningful for compound alleles.

### Step 5 — per-sample, per-genotype likelihood

For each sample `s`, enumerate all genotypes `G` over the merged
allele set under the configured ploidy (see §"Configurable
parameters § ploidy"). For each `G`:

- **Standard case** (no chain-broken compound in `G`): apply the
  freebayes-canonical closed-form
  ([calling_pipeline_architecture.md:1186-1193](../specs/calling_pipeline_architecture.md#L1186-L1193)):

  ```
  L_s(G) = Σ_{a ∉ G} S_s(a)
         + log multinomial(allele_probs(G), { count_s(a) : a ∈ G })
  ```

  with `allele_probs(G)` the ploidy-aware expected fractions
  (diploid biallelic `A/A → (1, 0)`, `A/T → (0.5, 0.5)`; higher
  ploidies analogous).

- **Chain-broken-compound case** (`G` contains a compound `C`
  for which sample `s` is chain-broken): the standard formula
  would force every constituent-supporting read into the
  error-cost term and bias the sample away from `C`
  proportionally to depth. Substitute the
  **constituents-independent** likelihood:

  ```
  L_s(G) = Σ_{a ∉ G_decomposed} S_s(a)
         + Σ over constituents c of C of
             log multinomial(p_decomposed_at_pos(c), counts_at_pos(c))
  ```

  where `G_decomposed` is the per-position decomposition of
  `G` (e.g. a genotype `C/REF` where `C` is "SNP at p1 +
  deletion at p2" decomposes to "SNP at p1 het, deletion at p2
  het"), and the multinomials are evaluated per-position
  independently. Conceptually, the chain-broken sample is
  evaluated *as if* the compound's constituents segregated
  independently. Mark `CA = 1` on `s` for `C`. The cohort prior
  `f_C` (estimated by Stage 6's M-step) carries the
  compound-vs-pair distinction for `s`.

  The justification is documented at
  [calling_pipeline_architecture.md:1098-1124](../specs/calling_pipeline_architecture.md#L1098-L1124)
  (the "Compound haplotype consistency" subsection).

### Step 6 — emit the merged record

Wrap the unified allele set, the per-sample scalar table, the
per-sample-per-genotype log-likelihood vectors, the CA flags, and
metadata (chrom_id, group span, ploidy, genotype enumeration
order) into a `MergedRecord` and hand it to Stage 6.

A merged record is **dropped** (Stage 5 emits nothing for the
group) if the merged allele set has only one allele after
unification (every sample is hom-ref by construction). This is
the residual filter after Stage 4's pure-REF drop: Stage 4 sees a
group with variant evidence, Stage 5 unifies and may find every
"variant" was filtered (e.g. by the chain-anchor rule rejecting
the only proposed compound and leaving only REF). Such groups
should not flow to Stage 6.

## What goes in a `MergedRecord` (output shape)

```rust
pub struct MergedRecord {
    pub chrom_id: u32,
    pub start: u32,                          // 1-based inclusive
    pub end: u32,                            // 1-based inclusive
    pub alleles: Vec<MergedAllele>,          // index 0 is always REF
    pub ploidy: u8,
    /// Per-sample scalar table.
    /// `scalars[sample_idx][allele_idx]` is the projected
    /// `AlleleSupportStats` for that sample and merged allele.
    /// A sample with no records in this group has every entry zeroed.
    pub scalars: Vec<Vec<AlleleSupportStats>>,
    /// Per-sample chain-anchor flags.
    /// `chain_anchored[sample_idx][allele_idx]` is `true` for a
    /// (compound) allele where this sample is chain-broken and the
    /// likelihood used the fallback path.
    pub ca_flags: Vec<Vec<bool>>,
    /// Per-sample per-genotype log-likelihood, evaluated under
    /// `c_s = 0` (no contamination correction). Index `[sample_idx]`
    /// is a vector of length `genotype_count(ploidy, alleles.len())`
    /// in the enumeration order defined by `genotype_order()`.
    /// Stage 6 either consumes these directly or re-derives from
    /// `scalars` when contamination correction is active.
    pub log_likelihoods: Vec<Vec<f64>>,
}

pub struct MergedAllele {
    /// Allele bytes projected onto `[start, end]` reference span.
    pub seq: Vec<u8>,
    /// `true` iff this allele is a cross-record compound; the
    /// per-sample CA flag is only meaningful for these.
    pub is_compound: bool,
    /// For `is_compound = true`: the indices of the constituent
    /// per-position records inside the source `OverlappingVarGroup`
    /// plus the local allele index. Used by Stage 6 for `f_C`
    /// estimation and by VCF emission for INFO field annotation.
    pub constituents: Vec<CompoundConstituent>,
}

pub struct CompoundConstituent {
    pub record_idx: usize,
    pub local_allele_idx: usize,
}

/// Canonical genotype enumeration order shared by all samples
/// in the record. For ploidy 2 and 3 alleles: AA, AB, BB, AC, BC,
/// CC (GATK / VCF convention).
pub fn genotype_order(ploidy: u8, n_alleles: usize) -> Vec<Vec<u8>>;
```

`ca_flags` and `scalars` are kept separately rather than rolled
into a per-genotype struct because Stage 6's contamination math
needs raw scalars per (sample, allele) and CA per (sample,
compound allele) — neither is per-genotype.

## Configurable parameters

Per the per-stage config convention. New struct:

```rust
pub const DEFAULT_MAX_ALLELES_PER_RECORD: usize = 6;
pub const DEFAULT_PLOIDY: u8 = 2;

#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct PerGroupMergerConfig {
    /// Cohort-wide ploidy. Plant target population has varied
    /// ploidies; per-sample ploidy is a deferred extension
    /// (see Out-of-scope follow-ups). Defaults to
    /// [`DEFAULT_PLOIDY`].
    pub ploidy: u8,

    /// Hard cap on the number of alleles in a single merged
    /// record. Excess alleles (lowest-evidence first) are
    /// dropped before genotype enumeration, with the dropped
    /// scalars folded into an "<OTHER>" bucket that contributes
    /// only to the error-cost term — equivalent to GATK's
    /// allele subsetting under `maxAlternateAlleles`. Defaults
    /// to [`DEFAULT_MAX_ALLELES_PER_RECORD`].
    pub max_alleles: usize,
}

impl Default for PerGroupMergerConfig {
    fn default() -> Self {
        Self {
            ploidy: DEFAULT_PLOIDY,
            max_alleles: DEFAULT_MAX_ALLELES_PER_RECORD,
        }
    }
}
```

**Default `max_alleles = 6`.** Diploid with 6 alleles gives 21
genotypes per sample; ploidy 4 with 6 alleles gives 126
genotypes. GATK's default is 6 (`maxAlternateAlleles = 6`); we
match it. Real cohorts rarely exceed 4 distinct alleles at a
single site; the cap is defensive against pathological repeat
regions.

**Allele subsetting under the cap.** When the merged allele set
exceeds `max_alleles`, drop the alleles with the lowest
cohort-wide total `count` (sum across samples of `count_s(a)`).
Their per-sample scalars are summed into the kept alleles'
"out-of-`G`" pool so they still contribute to the error-cost
term but are not enumerated as candidate genotypes. REF and any
chain-anchored compound are protected from the cap — a compound
that earned its place in the merged set via the chain-anchor
rule cannot be cap-dropped (dropping it would silently lose the
chain-anchored signal).

**Ploidy is cohort-wide.** A `--ploidy N` flag binds this. Mixed
ploidies (different per-sample ploidies in the same cohort) are
out of scope for v1.

## API shape

```rust
// src/cohort/per_group_merger.rs
use crate::cohort::variant_grouping::{GrouperError, OverlappingVarGroup};
use crate::per_sample_caller::pileup::AlleleSupportStats;

pub struct PerGroupMerger<I>
where
    I: Iterator<Item = Result<OverlappingVarGroup, GrouperError>>,
{ /* fields private */ }

impl<I> PerGroupMerger<I>
where
    I: Iterator<Item = Result<OverlappingVarGroup, GrouperError>>,
{
    pub fn new(upstream: I, ref_fetcher: Arc<dyn RefSeqFetcher>) -> Self;

    pub fn with_config(
        upstream: I,
        ref_fetcher: Arc<dyn RefSeqFetcher>,
        config: PerGroupMergerConfig,
    ) -> Self;

    pub fn config(&self) -> &PerGroupMergerConfig;
}

impl<I> Iterator for PerGroupMerger<I>
where
    I: Iterator<Item = Result<OverlappingVarGroup, GrouperError>>,
{
    type Item = Result<MergedRecord, PerGroupMergerError>;
}

#[derive(thiserror::Error, Debug)]
#[non_exhaustive]
pub enum PerGroupMergerError {
    #[error("upstream: {0}")]
    Upstream(#[from] GrouperError),

    #[error("reference fetch: {0}")]
    RefFetch(#[from] RefFetchError),

    /// Genotype enumeration produced a numerically degenerate
    /// likelihood (NaN, ±inf) for the named sample/genotype.
    /// This signals an internal bug, not a data condition — the
    /// closed-form formula is finite for all valid inputs.
    #[error(
        "degenerate likelihood at chrom {chrom_id} {start}-{end} \
         for sample_idx {sample_idx} genotype_idx {genotype_idx}: {kind}"
    )]
    DegenerateLikelihood {
        chrom_id: u32,
        start: u32,
        end: u32,
        sample_idx: usize,
        genotype_idx: usize,
        kind: DegeneracyKind,
    },
}
```

**Reference fetcher dependency.** Stage 5 needs to project each
per-position allele onto the group's reference span, which
requires fetching the REF bytes for `[group.start, group.end]`.
The same `RefSeqFetcher` trait Stage 1 uses
([src/per_sample_caller/ref_fetcher.rs](../../src/per_sample_caller/ref_fetcher.rs))
applies. `Arc` because Stage 5 runs under rayon and workers share
the fetcher.

**Errors latch.** Per the merger/grouper precedent
([per_position_merger.rs:130-135](../../src/cohort/per_position_merger.rs#L130-L135)):
once any error surfaces, subsequent `next()` calls return `None`.

**Parallelism is internal to `next()`.** The architecture spec
mandates rayon-parallel processing across groups
([calling_pipeline_architecture.md:1755-1768](../specs/calling_pipeline_architecture.md#L1755-L1768))
with an output reorder buffer to preserve genomic order. The
iterator's `next()` consumes a batch of groups from upstream,
processes them in parallel via `rayon::par_iter`, and emits
results in input order. Batch size is internal (configurable via
`PerGroupMergerConfig::batch_size`, default ~16-64; tune against
real data). Memory bound is `O(batch_size × max_alleles ×
N_samples)`.

**Streaming-friendly iteration.** Despite the internal batch,
the iterator's external shape is one `MergedRecord` per call.
Callers can compose it with Stage 6 cleanly.

Typical wiring:

```rust
let merger = PerPositionMerger::new(iters, sample_names, chromosomes)?;
let grouper = VariantGrouper::new(merger);
let stage5 = PerGroupMerger::new(grouper, Arc::clone(&ref_fetcher));
for record in stage5 {
    let record = record?;
    // hand to Stage 6
}
```

## Algorithm details

### Allele unification, line by line

```
fn unify_alleles(
    group: &OverlappingVarGroup,
    ref_seq: &[u8],         // bytes for [group.start, group.end]
    config: &PerGroupMergerConfig,
) -> UnifiedAlleleSet {
    let group_span = (group.end - group.start + 1) as usize;
    let n_samples = group.records[0].per_sample.len();
    let mut allele_map: HashMap<Vec<u8>, AlleleInfo> = HashMap::new();

    // REF is always allele 0.
    allele_map.insert(ref_seq.to_vec(), AlleleInfo::ref_default());

    // For each (sample, record, local allele) project onto the
    // group span and dedupe.
    for (record_idx, pp) in group.records.iter().enumerate() {
        let local_offset = (pp.pos - group.start) as usize;
        let local_span = /* ref_span of this record */;
        for (sample_idx, slot) in pp.per_sample.iter().enumerate() {
            let Some(rec) = slot else { continue };
            for (local_allele_idx, allele) in rec.alleles.iter().enumerate() {
                let projected = project_onto_group(
                    ref_seq,
                    local_offset,
                    local_span,
                    &allele.seq,
                );
                allele_map
                    .entry(projected)
                    .or_insert_with(AlleleInfo::default)
                    .add_source(sample_idx, record_idx, local_allele_idx);
            }
        }
    }

    // Detect cross-record compounds and apply the chain-anchor
    // rule (Step 2).
    let compounds = detect_compound_candidates(&allele_map, group);
    for compound in compounds {
        if !at_least_one_sample_chain_anchors(&compound, group) {
            // Reject; compound's constituents stay as independent
            // alleles, no entry added.
            continue;
        }
        let projected = compound.project(ref_seq);
        allele_map.entry(projected).or_insert_with(...).mark_compound();
    }

    // Apply the max_alleles cap.
    let mut alleles: Vec<MergedAllele> = allele_map.into_alleles();
    if alleles.len() > config.max_alleles {
        alleles = subset_alleles_by_cohort_count(alleles, config.max_alleles);
    }

    UnifiedAlleleSet { alleles, ref_idx: 0 }
}
```

### Compound-candidate detection

A **compound candidate** is a hypothetical allele whose support
in some sample is split across multiple records. To enumerate
candidates without combinatorial explosion:

1. For each sample, build a `chain_id → (record_idx,
   local_allele_idx) set` map by iterating over the sample's
   non-REF allele observations across all records of the group.
   Each chain id will appear in the map for every record where
   that read supported a non-REF allele.
2. Any chain id that appears in the map with two or more
   distinct `(record_idx, local_allele_idx)` entries proposes a
   compound: those entries' alleles co-occurred on this read.
3. Group such proposals by the set of `(record_idx,
   local_allele_idx)` pairs. Each distinct grouping is a
   candidate compound.
4. For each candidate, count how many chain ids across all
   samples support it. If at least one sample's count is ≥ 1
   (i.e. that sample chain-anchors), the candidate enters the
   merged set.

In practice the chain id map is short (a typical sample has tens
to hundreds of chain ids in a group, most spanning only one
record), and the candidate count is small (typically 0–3
distinct cross-record compounds per group). The algorithm is
`O(N_samples × N_chain_ids × N_records_per_chain)`; the constant
is tight.

### Likelihood computation, closed-form

For sample `s` and genotype `G` over the merged allele set:

```
log L_s(G) = error_cost + allele_balance

error_cost = Σ over alleles a ∉ G of S_s(a)
             // alleles include any subset-dropped "<OTHER>"
             // bucket so that pre-cap evidence is preserved

allele_balance = log multinomial(p_G, n_G)
   where:
     p_G[i] = (count of allele a_i in G) / ploidy
     n_G[i] = count_s(a_i) for each a_i ∈ G
     log multinomial(p, n) = log(N!) − Σ log(n_i!) + Σ n_i log p_i
```

Numerically:

- Compute `Σ n_i` once per `(s, G)` for the multinomial normalising
  constant.
- Use `lgamma` (log-gamma) for factorials to avoid overflow.
- A genotype with `p_G[i] = 0` for some `i` with `n_G[i] > 0`
  has `n_G[i] log 0 = −∞`; that genotype's log-likelihood is
  `−∞` (the reads contradict the genotype completely). The
  `f64::NEG_INFINITY` value propagates cleanly through Stage
  6's log-space EM. **Caveat:** `0 × log 0` should evaluate to
  `0`, not `NaN`; the code uses an explicit branch:

  ```rust
  fn xlogy(n: f64, p: f64) -> f64 {
      if n == 0.0 { 0.0 } else { n * p.ln() }
  }
  ```

  This is the convention adopted across the project's existing
  posterior code
  ([genotype_posteriors.rs](../../src/genotype_posteriors.rs))
  and tested there.

### Chain-broken compound fallback, explicit

For sample `s` and genotype `G` containing a compound `C` for which
`s` is chain-broken:

```
L_s(G) = error_cost_decomposed(G, s)
       + Σ over constituent positions p of C of
           log multinomial(p_G_at_pos[p], counts_at_pos[p](s))
```

where:

- `G` is decomposed to its per-position effect: e.g. `C/REF` for
  `C = "SNP at p1 + DEL at p2"` decomposes to "SNP-het at p1,
  DEL-het at p2".
- `error_cost_decomposed` uses per-position scalars rather than
  whole-group scalars: at each constituent position, accumulate
  `S_s(a)` for alleles that the decomposed `G` does not contain
  at that position.

This is structurally the same closed-form as the standard case,
just applied per-position rather than at the compound level.

### Genotype enumeration

Diploid biallelic gives 3 genotypes; diploid triallelic gives 6;
ploidy `P` with `K` alleles gives `C(P + K − 1, K − 1)`. For the
cap `max_alleles = 6` and ploidy ≤ 4, the maximum genotype count
per sample per record is `C(4+5, 5) = 126`. Enumeration uses the
standard combinatorial walk; see
[GenotypeIndexCalculator.indexOfFirstGenotypeWithAllele](../../gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/genotyper/GenotypeIndexCalculator.java)
for the GATK precedent of the indexing math we can reuse.

The genotype enumeration order is **canonical** (lexicographic
on allele indices, e.g. AA, AB, BB, AC, BC, CC for ploidy 2 with
3 alleles) and **shared across samples** in the record. Stage 6
indexes per-sample posteriors by the same order.

## Stage 6 hand-off contract

The `MergedRecord` is everything Stage 6 needs:

- **For `c_s = 0` (default)**: use `log_likelihoods` directly.
  The values are exact (no approximation) and never need
  recomputation across EM iterations.
- **For `c_s > 0` (contamination correction active)**: re-derive
  per-iteration likelihoods from `scalars` using the
  homogeneous-quality approximation
  ([architecture spec §"Scalar-storage approximation"](../specs/calling_pipeline_architecture.md#L1352-L1377)).
  The `log_likelihoods` field becomes a "starting point" for
  iteration 0; subsequent E-steps recompute from `scalars` with
  updated `c_s` and `q_b`.
- **`ca_flags` is consumed by Stage 6's M-step on `f_C`**:
  cohort-derived compound frequencies are estimated from the
  chain-evident samples (CA = 0); chain-broken samples (CA = 1)
  consume those frequencies as their prior.
- **Genotype order** is implicit in the enumeration helper —
  Stage 6 reconstructs the same order via `genotype_order()`.

The merged record is **read-only** from Stage 6's perspective:
Stage 6's EM does not mutate it, only reads scalars and writes
posteriors into its own output table. Each
`MergedRecord` is independent (no cross-record state), which is
what makes Stage 6's EM cohort-wide-not-per-record.

## Out of scope

- **Stage 6 itself.** Posterior engine + EM lives in
  `src/cohort/posterior_em.rs` (or similar), separate plan.
- **Contamination correction inside Stage 5.** The mixture math
  is applied at Stage 6's E-step where `c_s` is iteratively
  estimated. Stage 5 outputs the scalars and the `c_s = 0`
  likelihoods; Stage 6 augments.
- **`f_C` estimation.** Cohort-derived compound frequency is a
  Stage 6 quantity; Stage 5 only emits the CA flags so Stage 6
  knows which samples need it.
- **Observation-bias priors.** The bias scalars (`fwd`,
  `placed_left`, `placed_start`) are carried through verbatim
  in the scalar table. Stage 6 (or a downstream filter) computes
  the bias-correction multipliers.
- **VCF emission.** Stage 5 emits Rust `MergedRecord` structs,
  not VCF lines. A separate cohort `VcfWriter` (modelled on
  [src/vcf_writer.rs](../../src/vcf_writer.rs)) maps
  `MergedRecord` + Stage 6 posteriors → VCF output.
- **Mixed per-sample ploidies.** Single cohort-wide ploidy in
  v1; per-sample ploidy as a deferred extension (see Risks).
- **CLI wiring.** `--ploidy`, `--max-alleles`, and any future
  Stage 5 flags get bound when the cohort CLI subcommand lands.
- **The cohort-only inference of compounds with no chain
  anchor anywhere.** v2 extension scoped at the end of Stage 6
  ([architecture spec §"Future extension: cohort-only inference"](../specs/calling_pipeline_architecture.md#L1910)).

## Test strategy

Stage 5 tests live in
`src/cohort/per_group_merger.rs`'s `#[cfg(test)]` module and use
synthetic `OverlappingVarGroup` fixtures (no `.psp` round-trip).
Fixture builders extend those already in
[cohort/variant_grouping.rs:tests](../../src/cohort/variant_grouping.rs)
to construct `PerPositionPileups` with arbitrary chain-id
configurations.

Required fixture builders (add to a shared test-helpers module):

```rust
fn snp_obs(seq: &[u8], q_sum: f64, num_obs: u32, chain_ids: &[u64]) -> AlleleObservation;
fn record_with_alleles(chrom_id: u32, pos: u32, alleles: Vec<AlleleObservation>) -> PileupRecord;
fn pp_multi_sample(chrom_id: u32, pos: u32, per_sample: Vec<Option<PileupRecord>>) -> PerPositionPileups;
fn group_from(records: Vec<PerPositionPileups>) -> OverlappingVarGroup;
fn ref_fetcher_from_bytes(seq: &[u8], offset: u32) -> Arc<dyn RefSeqFetcher>;
```

### Reproducing the old gVCF merger tests

For every test in
[tests/genotype_merging_test.rs](../../tests/genotype_merging_test.rs),
add the Stage 5 equivalent built on `OverlappingVarGroup` input.
The mapping table — this is the **porting checklist** that gates
the gVCF cleanup commit:

| Old test | Stage 5 equivalent | What's reproduced | Notes |
|---|---|---|---|
| `test_simple_merge` | `stage5_simple_snp_merge_across_two_samples` | merged allele set has REF + 2 alts; both samples evaluated against full set | trivial |
| `test_simple_insertion` | `stage5_simple_insertion_across_samples` | insertion at p=10 in s0 + SNP at p=10 in s1 unify; merged set has 3 alts | uses INS allele bytes |
| `test_simple_deletion` | `stage5_simple_deletion_across_samples` | DEL at p=1 in s0 (homozygous) + SNP at p=2 in s1; merged ref span covers both | covers ref-span projection |
| `test_deletion_len_2` | `stage5_deletion_len_2` | 2 bp deletion spans two positions; merged set has 3 alts | spans cross-record |
| `test_overlapping_deletions` | `stage5_overlapping_deletions_unify` | s0 has DEL "AT→A" + SNP at p=3; s1 has SNP at p=1 + DEL "TT→T"; merged set has 3 alts (AT, AC, GT) | canonical form |
| `test_het_deletion` | `stage5_het_deletion` | both samples het DEL; per-sample likelihood for `DEL/REF` is the maximum | het-allele coverage |
| `test_two_deletions_in_same_sample` | `stage5_two_dels_same_sample_chain_anchored` | s0 has two distinct DEL alleles in same group; uses **chain ids** to anchor as compound iff a read spans both | new chain-id-driven shape |
| `test_two_overlapping_deletions_in_same_sample` | `stage5_two_overlapping_dels_same_sample_chain_anchored` | s0 has two overlapping DELs co-occurring on same haplotype; chain-anchored compound enters merged set | this is the canonical compound case |
| `test_allele_merge_with_non_ref_filtered` | `stage5_pure_ref_only_group_drops` | a group where every sample has only REF allele after unification produces no merged record | Stage 4 already filters most; Stage 5 catches residual |
| `test_missing_allele_in_merged_haplotype` | `stage5_sample_with_no_record_at_some_position` | sample has slots `None` at some positions; scalars zero out cleanly; likelihood degenerates gracefully | tests `<NON_REF>`-equivalent path |
| `test_group_merging_creates_correct_number_of_vars` | covered by emission count test below | n/a | grouper-level, not Stage 5 |
| `test_non_variant_vars_are_removed` | `stage5_residual_pure_ref_after_unification` | a group where unification yields only REF allele produces no record | residual filter |
| `test_analyze_single_group_from_merged_bin` | `stage5_one_record_per_group` | one `MergedRecord` emitted per input group (no record splitting) | structural |

**Phase-equivalent chain-id scenarios.** The `test_phase_*` family
(11 tests) all pin VCF-phase scenarios that the new pipeline
expresses via chain ids. The mapping is **not** test-name-for-test-name
because VCF phase is a per-record bit while chain ids are per-read;
some old scenarios are unrepresentable in the new shape (because the
new shape doesn't have a per-record phase bit to break) and some
new scenarios are unrepresentable in the old shape (because the
old shape doesn't have within-sample compound chain anchoring).
The functional equivalents:

| Old phase scenario | Stage 5 equivalent | Behaviour |
|---|---|---|
| `test_phase_single_het_no_phase_needed` | `stage5_single_het_no_compound_needed` | one het site, no compound proposal | trivial |
| `test_phase_kept_between_hets` | `stage5_compound_chain_anchored_two_hets_one_read` | two hets with a read spanning both → compound enters | chain-anchored |
| `test_phase_false_on_first_position_still_solvable` | `stage5_compound_partial_chain_evidence` | only one sample chain-anchors; that sample is chain-evident, others are chain-broken with CA=1 | tests fallback path |
| `test_phase_first_het_not_at_first_position` | `stage5_compound_offsets_into_group` | compound spans positions p+2, p+5 of a group starting at p | offset arithmetic |
| `test_phase_broken_between_hets_missing_genotype` | `stage5_no_chain_anchor_anywhere_constituents_only` | no sample chain-anchors the compound → constituents emitted independently, no compound in merged set | rejection path |
| `test_phase_broken_in_one_sample` | `stage5_chain_broken_sample_uses_fallback` | one chain-evident sample, one chain-broken sample at the same compound; CA=1 only on broken | core fallback test |
| `test_phase_lost_in_het` | `stage5_chain_broken_at_het_genotype` | het genotype at compound for chain-broken sample uses constituents-independent likelihood | likelihood-shape test |
| `test_phase_lost_in_het_and_not_recovered` | `stage5_chain_broken_persists_across_em` | the chain-broken sample's CA=1 flag carries through to Stage 6 hand-off | hand-off integration |
| `test_phase_lost_in_het2` | `stage5_chain_broken_at_hom_alt` | hom-alt-at-compound for chain-broken sample uses fallback (no error-cost mis-assignment) | likelihood-shape test |
| `test_phase_maintained_despited_not_phased_in_hom_position` | `stage5_hom_position_does_not_affect_chain_anchor` | hom-ref positions do not break chain anchoring | invariance test |
| `test_phase_conserved_in_het` | `stage5_chain_evident_full_path` | a sample that chain-anchors and is itself heterozygous at the compound uses the standard likelihood path | sanity |

### New tests the old suite couldn't express

The old gVCF suite did not have:

- **Multi-sample chain-anchored compounds** with mixed chain-evident
  and chain-broken samples. Add 4–6 tests covering:
  - Both samples chain-evident at the compound (likelihood
    matches standard formula for both).
  - One chain-evident, one chain-broken (different likelihood
    paths, different CA flags).
  - No sample chain-evident, compound rejected (constituents
    only in merged set).
  - Chain-evident sample has a single read supporting the
    compound (boundary: |chain_id ∩| = 1).
- **`max_alleles` cap behaviour.**
  - 7-allele group with `max_alleles = 6`: lowest-evidence
    allele dropped; its scalars roll into the "<OTHER>" bucket;
    remaining 6 alleles enumerated.
  - Chain-anchored compound protected from cap: 6 SNPs + 1
    compound with `max_alleles = 6` keeps the compound, drops
    the lowest-evidence SNP.
- **Ploidy variation.**
  - Triploid sample with 3-allele genotype enumeration.
  - Tetraploid biallelic — 5 genotypes (AAAA, AAAB, AABB, ABBB,
    BBBB); confirm allele-probs vector is `(1, 0)`, `(3/4, 1/4)`,
    `(1/2, 1/2)`, `(1/4, 3/4)`, `(0, 1)`.
- **Edge cases of the closed-form formula.**
  - `0 × log 0` in the multinomial (a genotype excluding an
    allele the sample never observed) — must produce `0.0`,
    not `NaN`.
  - Genotype completely contradicted by reads (`L = −∞`) — must
    propagate cleanly.
  - All samples hom-ref at a chain-anchored compound site —
    the compound is in the merged set (one sample chain-anchored
    it) but no sample's MAP genotype includes it.
- **Reference projection edge cases.**
  - Compound spans a deletion that consumes the anchor base of
    a downstream SNP.
  - Insertion at p=group.start where the projected allele is
    longer than `group.end - group.start + 1`.

### Property tests (proptest)

Two property invariants are worth proptest coverage:

1. **Likelihood under hom-ref input is dominated by the multinomial.**
   For random groups with only REF observations across all
   samples (no variant alleles), every sample's `L(REF/REF)` is
   the maximum across genotypes. (Trivial soundness check on
   the closed-form.)
2. **Compound rejection is symmetric in sample order.** Permuting
   `OverlappingVarGroup.records[i].per_sample` does not change
   whether a compound is chain-anchored or rejected — the
   chain-anchor check is sample-set-symmetric (chain ids are
   per-sample, but the "at least one sample chain-anchors" rule
   is symmetric). proptest randomly permutes samples and asserts
   the same merged allele set up to row reordering.

## Validation

Inside the dev container (`./scripts/dev.sh`):

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib --tests`
- `cargo build --examples --benches`

End-to-end validation against the old gVCF test catalogue is the
porting checklist above. Each row of the table must have a
passing Stage 5 test before the gVCF cleanup commit can land.

No new benchmark in this plan; benchmarking Stage 5 in isolation
over-fits the synthetic generator. The natural benchmark is
end-to-end through Stage 5 + Stage 6 on a realistic cohort, which
comes later.

## Assumptions / silent choices

- **`MergedRecord` is heap-allocated per group**, not pooled.
  Each group is independent and lives on a rayon worker; pool
  reuse would couple workers. Reach for pooling only if profiling
  shows allocator pressure.
- **`scalars[sample_idx]` always has length `alleles.len()`.**
  Samples with no record in the group still get a row of zeros.
  Costs ~40 B per sample per allele but keeps Stage 6's indexing
  uniform.
- **Single cohort-wide ploidy.** Mixed-ploidy cohorts are out of
  scope for v1. Per-sample ploidy would require enumeration
  changes Stage 6 doesn't currently support either.
- **Homogeneous-quality approximation for cross-record compound
  `q_sum`** (Step 3). Matches the contamination-layer approximation
  the architecture spec already commits to.
- **REF allele always at index 0.** Stage 6, VCF emission, and
  any downstream consumer rely on this. The convention matches
  Stage 1 / Stage 2 / the walker.
- **Allele cap drops by cohort-wide total count, not per-sample
  evidence.** A rare allele present in many samples beats a
  same-count allele present in one sample, which is the right
  bias for cohort-level calling.
- **`f64` likelihoods, not `f32`.** Genotype enumeration produces
  log-likelihoods at Phred-equivalent precision; `f32` would lose
  significant bits at deep coverage.
- **`Arc<dyn RefSeqFetcher>` for the reference fetcher.** Workers
  share immutably; the underlying fetcher is internally
  thread-safe (mmap-based, see
  [ref_fetcher.rs](../../src/per_sample_caller/ref_fetcher.rs)).

## Risks

- **Compound-candidate enumeration adversarial inputs.** A region
  where many distinct chain ids each link multiple records could
  in principle produce a combinatorial explosion of candidate
  compounds. In practice chain ids are per-read and reads support
  few alleles, but a defensive cap on candidates (perhaps
  `max_alleles * 2`) is a v2 follow-up if real data trips it.
- **Numerical stability at very deep coverage.** `count` fields
  are `u32`; cumulative `q_sum` is `f64`. At 100000× depth the
  multinomial's `log(N!)` is large enough that subtraction-based
  rearrangements could lose precision; the closed-form
  implementation should use `lgamma` rather than direct
  factorial computation.
- **Per-sample ploidy in v2 is non-trivial.** Genotype enumeration
  becomes per-sample; the per-sample scalar table grows; Stage 6's
  HW prior structure changes. Flag now so the v1 ploidy decision
  is visibly load-bearing.
- **Mixed ploidies vs `f_C` estimation.** Compound frequency
  `f_C` is conceptually a cohort-wide single number; under
  mixed ploidies the meaning of "frequency" shifts and the EM's
  M-step needs adjustment. Out of scope for v1.

## Out-of-scope follow-ups

- **Per-sample ploidy** via `--ploidy-file` mapping
  `sample_id → ploidy`. Defer until a real plant cohort with
  mixed ploidies requires it.
- **Stage 5 batch size tuning.** Default ~16-64 batch; revisit
  after end-to-end profiling.
- **Contamination correction inside Stage 5's `log_likelihoods`
  pre-computation.** Currently Stage 5 always emits `c_s = 0`
  likelihoods and Stage 6 re-derives with contamination active.
  An optimisation would have Stage 5 precompute the
  contamination-augmented likelihood once Stage 6 has converged
  on `c_s`, but Stage 5 is upstream of Stage 6 in the iterator
  chain — would require a second pass.
- **CLI binding for `PerGroupMergerConfig`.** Ships with the
  cohort CLI; not in this plan.
- **VCF emission of `MergedRecord` + Stage 6 posteriors.**
  Separate module.
- **Cohort-only compound inference (no chain anchor)** — v2,
  scoped at end of architecture spec
  [§"Future extension: cohort-only inference"](../specs/calling_pipeline_architecture.md#L1910).

## Deferred cleanup: completing the gVCF-path removal

When this plan ships and every row of the porting-checklist table
in §"Test strategy §Reproducing the old gVCF merger tests" has a
passing Stage 5 test, the old gVCF path becomes deletable. The
cleanup scope, per
[cohort_variant_grouping.md §"Deferred cleanup"](cohort_variant_grouping.md#L635-L692),
is:

- `src/variant_grouping.rs` — old grouper, delete.
- `src/genotype_merging.rs` — old merger, delete.
- `src/pipeline.rs` — gVCF entry point, delete.
- `tests/genotype_merging_test.rs` — old merger tests, delete
  only after the porting checklist confirms every row.
- `tests/variant_group_test.rs` — old grouper tests, delete.
- `src/main.rs` — drop the gVCF CLI subcommand wiring (which
  may temporarily leave the binary without a non-Stage-1 entry
  point; the cohort subcommand lands separately).
- `src/gvcf_parser.rs` and any other gVCF-path-only modules
  the compiler points at after the above.

The cleanup is **one commit**, separate from this plan's commit,
and is the trigger for renaming any cohort-side module that
collides on leaf-name with the old gVCF code (none currently,
but worth pinning).

## File touch list

**This plan's commit** (the new Stage 5 merger, no removals):

- `src/cohort/mod.rs` — add `pub mod per_group_merger;`.
- `src/cohort/per_group_merger.rs` — new file:
  `PerGroupMerger`, `MergedRecord`, `MergedAllele`,
  `CompoundConstituent`, `PerGroupMergerConfig` (with
  `DEFAULT_MAX_ALLELES_PER_RECORD` and `DEFAULT_PLOIDY`
  consts), `PerGroupMergerError`, the
  allele-unification / compound-detection / likelihood
  helpers, the genotype-enumeration helper, fixture builders,
  full `#[cfg(test)]` module covering every row of the
  porting-checklist table plus the new test categories above.
- (No changes to `src/lib.rs` — `pub mod cohort;` is already
  exported.)

No changes to existing types. `PerGroupMerger` consumes
`OverlappingVarGroup` and `GrouperError` from the variant
grouping module unchanged.

**Deferred follow-up commit** (the gVCF cleanup): scope and
gating in §"Deferred cleanup" above. Not done as part of this
plan's commit.
