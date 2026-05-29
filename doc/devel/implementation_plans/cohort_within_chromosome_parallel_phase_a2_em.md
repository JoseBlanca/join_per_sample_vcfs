# Phase A.2 — column-native EM input boundary

Status: planned. Successor to Phase A.1 layer 4 (commit `6e69cbe`).
See the
[Phase A impl report](../reports/implementations/cohort_within_chromosome_parallel_phase_a_2026-05-28.md)
for everything that landed in Phase A.1.

## Where we are after Phase A.1

The cohort var-calling pipeline now reads chunked columnar data
from PSPs, partitions windows into variant groups, and for each
group runs the column-native kernel chain end to end:

```
chunk (columns) ─► layer 1 (unify_alleles_columnar) ─► UnifiedAllelesColumns
                                                              │
                                                              ▼
                ─► layer 2 (project_scalars_columnar) ─► ProjectedScalarsColumns
                                                              │
                                                              ▼
                ─► layer 3 (compute_log_likelihoods_columnar) ─► LogLikelihoodsColumns
                                                              │
                                                              ▼
                ─► build_merged_record_columnar ─► MergedRecord ─► PosteriorEngine ─► PosteriorRecord
                       ↑                                                                      │
                       └── (Phase A.1 adapter; this plan removes it) ────────────────────────►│
                                                                                              ▼
                                                                                         VCF writer
```

The byte-identity gate has been met on the 3-tomato cohort (VCF
body diff vs `main` = 0 lines). Phase A.1's
[`build_merged_record_columnar`](../../../src/var_calling/cohort_block/worker.rs)
is the boundary that still produces a row-shape
[`MergedRecord`](../../../src/var_calling/per_group_merger.rs)
per group — the EM is the last consumer that hasn't been adapted
to the columnar layer outputs.

## Why Phase A.2

The `MergedRecord` boundary is overhead for two reasons:

1. **Clones.** Every per-group call to
   `build_merged_record_columnar` clones four buffers from the
   scratch into the `MergedRecord`:
   - `scalars: Vec<AlleleSupportStats>` — `n_samples × n_kept`.
   - `other_scalars: Vec<AlleleSupportStats>` — `n_samples`.
   - `chain_anchor_flags: Vec<bool>` — `n_samples × n_kept`.
   - `log_likelihoods: Vec<f64>` — `n_samples × n_genotypes`.
   And it builds a `Vec<MergedAllele>` (one `Vec<u8>` per allele)
   from `UnifiedAllelesColumns`. The EM then moves these into the
   emitted `PosteriorRecord`. The clones are not free (small alloc
   + memcpy per group), and the
   [perf review](../../../doc/devel/PROJECT_STATUS.md) flagged the
   per-group allocator pressure as the dominant cost at large N.

2. **The `chain_anchor_flags` table is built twice.** The
   columnar pipeline already encodes the same information in
   `UnifiedAllelesColumns::chain_anchor_counts` — `flag = (count
   == 0)`. The worker's `build_merged_record_columnar` reshapes
   that into a flat `Vec<bool>` so the EM can index it directly;
   the EM only ever reads the flag, never the count. This is
   pure transcoding work — O(`n_alleles × n_samples`) per group —
   that disappears when the EM reads `unified.chain_anchor_count`
   directly.

3. **Cleaner data flow.** Layer 4 is the only remaining place in
   the cohort pipeline that materialises a row-shape per-group
   record. Removing it lets the kernel chain end at the EM's
   inputs rather than at an intermediate adapter.

Phase A.2 ports the EM's input boundary from "consume an owned
`MergedRecord`" to "read from layer 2 + 3 outputs by borrow" while
leaving the EM body (E-step, M-step, mixture, posterior summary,
exact-AF QUAL) unchanged.

## What the EM actually reads

A focused audit of
[`run_em_for_record`](../../../src/var_calling/posterior_engine.rs)
(the row-shape entry point) shows what it consumes from
`MergedRecord`:

| Field | Used by | Access pattern |
|---|---|---|
| `chrom_id` / `start` / `end` | `RecordLocus` for error context | read once |
| `ploidy` | every EM step, genotype enumeration | scalar |
| `n_samples` | every loop bound | scalar |
| `n_genotypes` | every loop bound | scalar |
| `alleles[0].seq.len()` | `classify_allele` → REF span | read once |
| `alleles[a].is_compound` | classify + compound_mask | per allele |
| `alleles[a].seq.len()` | classify → SNP vs indel | per allele |
| `alleles[a].constituents` | not read by EM (forwarded to PosteriorRecord only) | — |
| `scalars[s, a]` | mixture pre-pass (when contamination on) | sample-major flat |
| `other_scalars[s]` | mixture pre-pass | per sample |
| `chain_anchor_flags[s, a]` | EM does not read this directly — but it is forwarded to PosteriorRecord, which the VCF writer uses for the `CA` flag | sample-major flat |
| `log_likelihoods[s, g]` | E-step / mixture pre-pass | sample-major flat |

Everything else on `MergedRecord` (the `Vec<MergedAllele>`'s
seq + constituents) is forwarded straight to
`PosteriorRecord.alleles` for the VCF writer's `ALT` / compound
emit path — the EM body never touches them.

**Implication:** the EM doesn't need `MergedRecord`. It needs
`(locus, ploidy, n_samples, n_genotypes, ref_len_for_classify,
allele_classes_iter, scalars: &[…], other_scalars: &[…],
chain_anchor_flags_iter, log_likelihoods: &[f64])`. Layers 2 + 3
already produce the slice-shaped buffers; layers 1 already encodes
the per-allele flags. Phase A.2 routes them directly.

## What stays unchanged

- The EM math: `e_step` / `e_step_simd` / `m_step_p_hat` /
  `m_step_f_hat_compound` / `run_em_loop`.
- `compute_mixture_log_likelihoods` / `_simd` — they read
  `scalars` + `log_likelihoods` by slice already.
- `summarise_posteriors`, `compute_qual_via_exact_af`.
- `RecordScratch`.
- `PosteriorEngineConfig` (incl. contamination, fixation index
  overrides, etc.).
- `PosteriorRecord` (still row-shape — the VCF writer's input).
  Its `Vec<MergedAllele>` is still built at emit time, but only
  once per group, at the boundary where the writer actually needs
  row-shape allele bytes.

## What changes

### A new EM entry point

`run_em_for_record` keeps its existing row-shape signature for
the test/fuzz surface and gains a thin borrowed-slice sibling
that the worker calls. The shared body lives behind a private
helper.

```rust
/// Borrowed-slice view of one group's EM inputs. Same logical
/// content as a `MergedRecord` minus the per-allele bytes —
/// callers pass an `AllelesView` for the per-allele
/// classification queries.
pub(crate) struct EmInputs<'a> {
    pub locus: RecordLocus,
    pub ploidy: u8,
    pub n_samples: usize,
    pub n_genotypes: usize,
    pub alleles: &'a dyn AllelesView,
    pub scalars: &'a [AlleleSupportStats],
    pub other_scalars: &'a [AlleleSupportStats],
    pub log_likelihoods: &'a [f64],
}

/// Minimal accessor surface for the per-allele queries the EM
/// makes. Two impls are provided:
///
/// - [`MergedAllelesView`] — wraps `&[MergedAllele]`. Used by
///   the row-shape `run_em_for_record` entry (tests / fuzz).
/// - [`ColumnarAllelesView`] — wraps `&UnifiedAllelesColumns`.
///   Used by the new worker path.
///
/// The EM only ever needs `len`, `ref_len`, `is_compound`,
/// `seq_len(idx)`. No constituent walks; no `seq_bytes` access.
pub(crate) trait AllelesView {
    fn len(&self) -> usize;
    fn ref_len(&self) -> usize;
    fn seq_len(&self, allele_idx: usize) -> usize;
    fn is_compound(&self, allele_idx: usize) -> bool;
}

pub(crate) fn run_em_columnar<M: MathBackend>(
    inputs: EmInputs<'_>,
    config: &PosteriorEngineConfig,
    math: &M,
    scratch: &mut RecordScratch,
) -> Result<EmOutputs, PosteriorEngineError>;
```

`EmOutputs` is a small owned bundle the worker turns into a
`PosteriorRecord` together with the per-group allele bytes:

```rust
pub(crate) struct EmOutputs {
    pub allele_frequencies: Vec<f64>,
    pub compound_frequencies: Vec<Option<f64>>,
    pub posteriors: Vec<f64>,
    pub best_genotype: Vec<usize>,
    pub gq_phred: Vec<f64>,
    pub qual_phred: f64,
    pub diagnostics: EmDiagnostics,
}
```

The existing `run_em_for_record(record: MergedRecord, …) ->
PosteriorRecord` becomes a thin shim around `run_em_columnar`
that:
1. Builds a `MergedAllelesView` over `record.alleles`.
2. Calls `run_em_columnar` with the record's slices.
3. Joins `EmOutputs` + `record.alleles` + `record.scalars` +
   `record.other_scalars` + `record.chain_anchor_flags` into
   the `PosteriorRecord` fields the writer consumes.

The trivial-record path (`alleles.len() < 2`) also splits along
the same boundary.

### A new worker emit path

`run_window` stops building `MergedRecord` per group. The
iterator that today yields `Result<MergedRecord, PerGroupMergerError>`
becomes a per-group consumer that:

1. Runs layers 1 + 2 + 3 into the scratch (unchanged).
2. Constructs `EmInputs` borrowing from the scratch + an
   `EmAllelesView` over `UnifiedAllelesColumns`.
3. Calls `run_em_columnar`.
4. Builds a `PosteriorRecord` by:
   - `Vec<MergedAllele>` materialised from the unified columns
     (this is the only `MergedAllele` allocation in the
     pipeline — driven by the VCF writer's row-shape requirement,
     not the EM).
   - `scalars` / `other_scalars` / `chain_anchor_flags` produced
     by:
     - `Vec::clone` from the scratch (simple, same allocs as
       today), OR
     - `std::mem::take` from the scratch (saves the clone, but
       resets capacity — only worth it if the per-record buffer
       allocator pressure dominates).
   - `EmOutputs` fields moved into `PosteriorRecord`.
5. Pushes the `PosteriorRecord` into the driver's output buffer.

`PerGroupMergerError` does not appear anywhere in the new path —
the iterator's error type collapses to
`PosteriorEngineError`. (The current
`ColumnarMergedRecordsIter` only used `PerGroupMergerError` to
satisfy the existing `PosteriorEngine` upstream contract; Phase
A.2 drops that contract.)

### A streaming `PosteriorEngine` shim that drains `EmInputs`

The existing
[`PosteriorEngine`](../../../src/var_calling/posterior_engine.rs)
type is `Iterator<Item = Result<MergedRecord, PerGroupMergerError>>`
based; it owns `RecordScratch` and drives `run_em_for_record`.
After Phase A.2:

- The row-shape `PosteriorEngine` stays as is — backward-compatible
  for the existing test surface and the `var-calling-from-bam`
  path that still uses the streaming row-shape pipeline.
- The cohort worker no longer instantiates `PosteriorEngine`. It
  owns its own `RecordScratch` directly (via
  `ColumnarPipelineScratch`) and calls `run_em_columnar` per group
  in the iterator's `next()` body. The driver's emit loop drains
  `PosteriorRecord`s from `output_buf` exactly as today.

The `PosteriorEngine` type therefore continues to exist as
the row-shape consumer; the cohort pipeline simply stops using it.

## Step-by-step phasing

Each step ships as its own commit + tests; the prior step's
output is the byte-identity oracle.

### Step 1 — split `run_em_for_record` along the row/columnar boundary

Introduce `EmInputs<'a>` + `AllelesView` + `MergedAllelesView` +
`run_em_columnar` (with `EmOutputs`). The existing
`run_em_for_record(record: MergedRecord, …) -> PosteriorRecord`
becomes a thin shim that builds `EmInputs` over `record`'s
slices, calls `run_em_columnar`, joins the result with the
record's allele bytes / scalars / flags into a `PosteriorRecord`.

**Files touched:**
- `posterior_engine.rs` — adds the new types + entry point, keeps
  the row-shape shim.
- `per_group_merger.rs` — no changes (still emits `MergedRecord`).

**Byte-identity oracle:** the existing
`PosteriorEngine::process` on a `MergedRecord` — unchanged
externally, internally now routes through `run_em_columnar`.
All existing `posterior_engine` tests + the cohort CLI
integration tests pass without modification.

This step is a pure refactor of the EM module. No behaviour
change. Diff is heavy in `posterior_engine.rs` (`record.alleles`
→ `inputs.alleles.…`, `record.scalars` → `inputs.scalars`, etc.)
but every hunk is a 1:1 field-access translation.

### Step 2 — `ColumnarAllelesView` + worker uses `run_em_columnar`

`ColumnarAllelesView<'a>` wraps
`&'a UnifiedAllelesColumns` and implements `AllelesView`:

```rust
impl AllelesView for ColumnarAllelesView<'_> {
    fn len(&self) -> usize { self.columns.n_alleles() }
    fn ref_len(&self) -> usize {
        // REF is allele 0; the per-allele seq slice's length is
        // the projected-onto-group-span length, which equals the
        // REF span by construction.
        self.columns.allele_seq(0).len()
    }
    fn seq_len(&self, idx: usize) -> usize {
        self.columns.allele_seq(idx).len()
    }
    fn is_compound(&self, idx: usize) -> bool {
        self.columns.is_compound[idx]
    }
}
```

Worker rewires
[`build_merged_record_columnar`](../../../src/var_calling/cohort_block/worker.rs)
into `build_posterior_record_columnar` that:
- Calls `run_em_columnar` with the scratch slices.
- Allocates `Vec<MergedAllele>` from `UnifiedAllelesColumns` at
  emit time only.
- Builds `chain_anchor_flags` (still row-shape — PosteriorRecord
  shape) by reading `unified.chain_anchor_count(allele, sample)
  == 0` directly into the PosteriorRecord's owned `Vec<bool>`.
  This is the same shape transformation as today's worker, but
  happens at emit time instead of at `MergedRecord` build time —
  net: same allocs, one less buffer threaded through the EM.
- Clones (`Vec::clone`) `scratch.projection.scalars` and
  `scratch.projection.other_scalars` directly into the
  `PosteriorRecord`.

The `ColumnarMergedRecordsIter` is removed (no more
`MergedRecord` upstream); `run_window` becomes a plain `for g in
0..partition.n_groups()` loop that pushes `PosteriorRecord`s into
`output_buf` directly.

**Files touched:**
- `cohort_block/worker.rs` — `MergedRecord` build path removed;
  new `ColumnarAllelesView` + `build_posterior_record_columnar`.
- `cohort_block/driver.rs` — no signature change (output_buf
  shape unchanged).

**Byte-identity oracle:** the layer 4 row-shape path (still
present via the `run_em_for_record` shim). Both paths produce
identical `PosteriorRecord`s; the new worker test compares them
on the existing layer 4 fixtures
([worker.rs:`run_window_columnar_matches_row_shape_pipeline`](../../../src/var_calling/cohort_block/worker.rs)
and the cap+compound regression fixture).

After step 2, the worker no longer produces `MergedRecord`. The
cohort CLI integration tests must still pass byte-identically.
The 3-tomato real-data VCF diff vs the layer 4 baseline must
remain 0 lines.

### Step 3 — drop the unused `MergedRecord` clones (optional)

If profiling after step 2 still shows the per-group clone of
`scratch.projection.scalars` / `other_scalars` on the hot path,
this step replaces the `Vec::clone` with `std::mem::take` (or
`std::mem::replace` against a pre-allocated dummy with the same
capacity as the high-water-mark). The scratch loses capacity
reuse but the per-group allocation drops by 2 Vecs.

**Trigger:** perf review motivates it. Otherwise skip; the clone
cost is small relative to the EM itself.

### Step 4 — perf review

Run
[`benchmarks/tomato1/scripts/perf_scaling_synthetic.py`](../../../benchmarks/tomato1/scripts/perf_scaling_synthetic.py)
on N = 50 / 200 / 1000 against the layer 4 baseline. Measure:
- Wall time delta.
- Peak RSS delta.
- Per-group allocator counts (dhat or `heaptrack`).

Expected gains:
- The `chain_anchor_flags` table build in
  `build_merged_record_columnar` (today's worker) is gone after
  step 2 — saving `n_alleles × n_samples` Vec writes per group.
- The `Vec<MergedAllele>` allocation at EM input time is gone —
  saving one allocation per group + per-allele Vec<u8> seq
  clones. The MergedAllele build still happens at PosteriorRecord
  emit but only once.
- Allocator pressure on the EM's input boundary is removed —
  the EM now reads from the scratch directly, no intermediate
  `MergedRecord` build.

Likely the wins are small relative to the EM body itself (the
mixture pre-pass + the genotype loop dominate); the SIMD review
(Phase D in the original plan) is where the bigger numbers live.
Phase A.2's primary value is data-flow clarity + the
column-native kernel chain reaching all the way to the EM's
inputs without an adapter; the perf delta is a secondary
benefit.

Update
[`doc/devel/reports/implementations/cohort_within_chromosome_parallel_phase_a_2026-05-28.md`](../reports/implementations/cohort_within_chromosome_parallel_phase_a_2026-05-28.md)
and [PROJECT_STATUS](../../../PROJECT_STATUS.md) with the
measured numbers.

## Validation

Each step ships with:

- **Unit tests** against the prior step's EM output on the same
  fixtures as Phase A.1 layer 4
  ([no-compound, multi-position, with-compound, cap-fires,
  cap+compound regression](../../../src/var_calling/cohort_block/worker.rs)).
- **Cohort CLI integration tests**
  ([`tests/cohort_cli_integration.rs`](../../../tests/cohort_cli_integration.rs))
  pass without modification — they exercise `var-calling`
  end-to-end through the new driver, including
  `var_calling_emits_deterministic_vcf_across_runs`.
- **3-tomato real-data byte-identity** vs the layer 4 baseline.
  VCF body diff = 0 lines after every step.

The `posterior_engine.rs` test suite (1000+ tests) must pass
after step 1 without modification — they all consume
`MergedRecord` and exercise the row-shape entry that stays in
place.

## Risks

### `MathBackend` monomorphisation

`run_em_for_record` is generic over `<M: MathBackend>`. The
backend monomorphises four ways
(`ExactMath`, `InterpUnivariateMath`, `InterpUnivariateSimdMath`,
plus test-only backends). The new `run_em_columnar` must preserve
the same generic structure verbatim — the `if M::HAS_LANE_4`
branches that dispatch into `e_step_simd` and
`compute_mixture_log_likelihoods_simd` cannot collapse without
re-monomorphising the SIMD vs scalar paths.

**Mitigation:** the refactor is mechanical — every site that
reads from `record` becomes a read from `inputs`. The branching
structure stays identical.

### `AllelesView` virtual dispatch cost

`&dyn AllelesView` adds a virtual-call per
`is_compound` / `seq_len` read inside the EM. The EM reads
these once per allele per record (in `classify_allele` and the
`compound_mask` fill); the hot loops read from
`scratch.compound_mask: &[bool]` — already a flat slice — not
through the view.

**Mitigation:** keep `compound_mask` exactly as it is — a flat
slice in `RecordScratch`, populated once per record from
`inputs.alleles.is_compound(a)`. The hot loops never see the
view. If profiling later shows virtual-call overhead at the
classify step, switch to a thin generic
`run_em_columnar<M, V: AllelesView>` and inline the view.

### `chain_anchor_flags` table forwarding

`PosteriorRecord` carries `chain_anchor_flags: Vec<bool>` for
the VCF writer's `CA` flag. Today the worker builds this table
in `build_merged_record_columnar` and the EM forwards it
unchanged through `MergedRecord`. In Phase A.2 the table
construction moves to PosteriorRecord emit time and is built
directly from `unified.chain_anchor_count(allele, sample) == 0`.

**Risk:** subtle off-by-one in the iteration order vs the
existing flat-row layout. The existing layer 4 worker code does
this exact transformation (worker.rs lines around the
`chain_anchor_flags` fill); Phase A.2 lifts that code from
"during MergedRecord build" to "during PosteriorRecord emit".

**Mitigation:** byte-identity test against the layer 4 baseline.
The existing fixtures exercise compounds + chain-broken samples
+ cap; if the transformation is wrong by even one element the
EM's output would differ.

### Allele-ref-length semantics

The EM's `classify_allele(allele, ref_len)` uses
`ref_len = alleles[0].seq.len()`. In `UnifiedAllelesColumns` the
per-allele seq is the projection onto the group's reference span;
`allele_seq(0).len()` IS the ref span length by construction
(REF projects to itself). So `ColumnarAllelesView::ref_len()` =
`self.columns.allele_seq(0).len()` — same answer, computed from
the columnar storage.

**Mitigation:** unit-test the column-native classify against the
row-shape classify on the same fixture set as
[`classify_alleles`](../../../src/var_calling/posterior_engine.rs).

### The trivial-record path

`run_em_for_record` has an `n_alleles < 2` early return that
emits `trivial_posterior_record` instead of running the EM.
Phase A.1 layer 4 already routes around this by returning
`Ok(None)` at the worker level (the row-shape merger also
returns `Ok(None)`). Phase A.2 keeps the same skip in the worker
— the EM entry never sees a 1-allele record. If a caller (test
or future code) passes a 1-allele record to `run_em_columnar`,
it should produce the same trivial output.

**Mitigation:** the trivial branch lives at the top of
`run_em_columnar` and produces an `EmOutputs` that the caller
joins with the alleles to form the trivial `PosteriorRecord`.
Unit-tested on the existing trivial-path fixtures.

## Out of scope

- **`PosteriorRecord` columnar refactor.** The original plan's
  `GroupedVariantsBatch` accumulator. The VCF writer still wants
  a row-shape input; refactoring it is a separate workstream.
- **`PileupRecordRef<'a>` borrowed-view threading.** The per-position
  `PileupRecord` allocations live one layer below this plan, in
  the chunk loader.
- **SIMD review (Phase D).** Phase A.2 doesn't touch the SIMD
  lanes inside `e_step_simd` /
  `compute_mixture_log_likelihoods_simd`.
- **`SharedRefFetcher` `'static` cleanup.** Tracked separately in
  the Phase A impl report's deferred-work list.

## Estimated effort

- **Step 1** (split `run_em_for_record`): ~1–2 sessions. Diff is
  heavy in `posterior_engine.rs` but mechanical; the test surface
  is large so each test the refactor breaks signals a missed
  field-access translation.
- **Step 2** (worker uses `run_em_columnar`): ~1 session. Touches
  only `worker.rs`. Byte-identity gate is the existing test set +
  the 3-tomato real-data diff.
- **Step 3** (optional clone elimination): half a session,
  triggered by perf review.
- **Step 4** (perf review): half a session, depending on whether
  re-running the scaling sweep is needed.

Total ~3 sessions, comparable to a single layer of Phase A.1.

## Sequencing

The plan is committed independently of Phase A.1 layer 4. There
is no rush to start Phase A.2 — layer 4 produces byte-identical
output and the column-native kernel chain is complete. Phase A.2
is the optional "clean up the EM input boundary" step that the
original within-chromosome plan called out under "Phase A.1 —
rewrite kernels native columnar."

If the cohort pipeline's perf budget is already met under layer
4, Phase A.2 can ship at the convenience of the maintainer.
Layer 4's correctness gate is the load-bearing one; Phase A.2 is
strictly a cleanup + minor perf step.

## Forward references

- [Phase A impl report](../reports/implementations/cohort_within_chromosome_parallel_phase_a_2026-05-28.md)
  — covers what landed in Phase A.0 + A.1.
- [Original within-chromosome plan](cohort_within_chromosome_parallel.md)
  — §"Phase A — single-threaded chunk loader, columnar from day
  one" anticipated this step under the heading "rewrite kernels
  native columnar"; Phase A.1 + A.2 together discharge that
  obligation.
- [`per_group_merger.rs`](../../../src/var_calling/per_group_merger.rs)
  — the row-shape `MergedRecord` shape that Phase A.2 stops
  building in the cohort pipeline. The type itself stays
  (consumed by the existing `PosteriorEngine` row-shape entry +
  the streaming `var-calling-from-bam` path).
- [`posterior_engine.rs`](../../../src/var_calling/posterior_engine.rs)
  — the module Phase A.2 step 1 refactors.
- [`cohort_block/worker.rs`](../../../src/var_calling/cohort_block/worker.rs)
  — the module Phase A.2 step 2 rewrites.
