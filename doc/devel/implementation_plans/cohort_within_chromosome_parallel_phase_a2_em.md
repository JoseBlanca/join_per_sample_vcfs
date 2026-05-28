# Phase A.2 — native column-native EM

Status: planned. Successor to Phase A.1 layer 4 (commit landing
the column-native pipeline + `MergedRecord` adapter; see the
[Phase A impl report](../reports/implementations/cohort_within_chromosome_parallel_phase_a_2026-05-28.md)).

## Why

Phase A.1 layer 4 made every per-group production-path computation
column-native through to the EM's inputs — `UnifiedAllelesColumns`,
`ProjectedScalarsColumns`, `LogLikelihoodsColumns` — but the EM
itself still consumes a row-shape
[`MergedRecord`](../../../src/var_calling/per_group_merger.rs).
The Phase A.1 worker bridges this by cloning the columnar layer
outputs into a fresh `MergedRecord` per group; the EM then runs
against the clone. This is a faithful port (byte-identical to
`main` on the 3-tomato cohort) but introduces a per-group
`MergedRecord` clone tax — six `Vec` allocations per group
(`alleles`, `scalars`, `other_scalars`, `chain_anchor_flags`,
`log_likelihoods`, `compound_constituents`).

Phase A.2 lifts the EM (and the mixture pre-pass + posterior
summarisation + exact-AF QUAL pass) to consume the columnar layer
outputs by slice/borrow instead of a cloned record. The EM math
itself is already columnar-shaped — only the entry-point signatures
change. Net effect: the per-group `MergedRecord` clone tax goes to
zero, and the EM hot path reads directly from buffers layers 2 + 3
already produced.

## Out-of-scope (deliberately deferred)

- **SIMD review.** Phase D in the original plan. After Phase A.2's
  signature refactor, `compute_log_likelihoods_columnar` (layer 3),
  the mixture pre-pass, and the EM E-step are the candidates;
  measure first.
- **`PileupRecordRef<'a>` deep threading.** Phase A.2 stays at the
  per-group level; the per-position `PileupRecord` allocations
  happen one layer down inside the chunk loader, and are out of
  scope here.
- **`GroupedVariantsBatch` columnar refactor.** The PosteriorRecord
  output shape stays row-shaped for Phase A.2; the VCF writer's
  per-record consumer doesn't yet benefit from a columnar batch.

## Approach — one step at a time

The plan splits the EM port into independently-shippable steps,
each byte-identity-tested against the prior step's EM output on
the same fixture set. Same shape as the Phase A.1 layer
progression.

### Step 1 — `run_em_for_record` boundary refactor

Move `run_em_for_record`'s entry from `record: MergedRecord` to
borrowed slices:

```rust
pub(crate) fn run_em_for_record<M: MathBackend>(
    inputs: EmInputs<'_>,
    config: &PosteriorEngineConfig,
    math: &M,
    scratch: &mut RecordScratch,
) -> Result<PosteriorRecord, PosteriorEngineError>
```

where `EmInputs<'a>` is a borrowed bundle:

```rust
pub(crate) struct EmInputs<'a> {
    pub locus: RecordLocus,
    pub ploidy: u8,
    pub n_samples: usize,
    pub n_genotypes: usize,
    pub alleles: AllelesView<'a>,
    pub scalars: &'a [AlleleSupportStats],
    pub other_scalars: &'a [AlleleSupportStats],
    pub chain_anchor_flags: &'a [bool],
    pub log_likelihoods: &'a [f64],
}
```

`AllelesView<'a>` exposes the per-allele methods the EM actually
calls (`seq.len()` for the REF span, `is_compound` flag access,
constituents iteration for the chain-broken path). A pair of impls
satisfies callers from both sides:

- `MergedRecordAllelesView<'a>` — wraps `&'a [MergedAllele]`. Used
  by the existing row-shape callers (tests, fuzz, the `main`
  branch's posterior-engine tests).
- `ColumnarAllelesView<'a>` — wraps `&'a UnifiedAllelesColumns`.
  Used by the new worker path.

After this step, the EM body is unchanged from `main`; only the
field-access syntax changes from `record.alleles[a]` to
`inputs.alleles.is_compound(a)` etc. Byte-identity is automatic.

**Byte-identity oracle:** the existing
`PosteriorEngine`-consumes-`MergedRecord` path. The worker calls
the new entry with a `MergedRecordAllelesView` built from its
local `MergedRecord` (zero extra copy — just borrows the existing
Vec). Phase A.1 byte-identity tests still pass.

**Files touched:** `posterior_engine.rs`,
`cohort_block/worker.rs` (the row-shape view construction is in
the worker), the existing `PosteriorEngine` consumers (an Arc, a
test fixture). Net diff is heavy in `posterior_engine.rs` but each
hunk is a 1:1 field-access translation.

### Step 2 — Worker switches to `ColumnarAllelesView`

After step 1 is byte-identical, the worker constructs `EmInputs`
from the columnar layer outputs directly:

- `alleles`: `ColumnarAllelesView::new(&scratch.unified)`.
- `scalars`: `&scratch.projection.scalars`.
- `other_scalars`: `&scratch.projection.other_scalars`.
- `chain_anchor_flags`: `&scratch.chain_anchor_flags`.
- `log_likelihoods`: `&scratch.log_likelihoods.log_likelihoods`.

The `MergedRecord` clone path is removed from the worker's hot
loop; the worker keeps no per-group allocations beyond what the
`PosteriorRecord` emit itself does.

**Byte-identity oracle:** the step-1 row-shape path (still
present, exercised by `MergedRecord` consumers in tests).

**Files touched:** `cohort_block/worker.rs` only.

### Step 3 — Wire `PosteriorEngine` for the borrowed iterator

`PosteriorEngine`'s iterator API today is
`Iterator<Item = Result<MergedRecord, PerGroupMergerError>>`. Step
3 generalises to a trait `PosteriorEngineUpstream` with two impls:

- The existing `Iterator<Item = Result<MergedRecord, …>>` for the
  row-shape consumers.
- A new `ColumnarUpstream<'a>` that the worker owns and that
  yields the column-native `EmInputs<'a>` directly (no
  `MergedRecord` ever constructed).

Step 3 is optional — step 2 already removes the per-group clone
tax. Step 3 only matters if profiling later shows the
`MergedRecord` allocation in `next()` (the iterator's own
output stage) remaining on the hot path. Defer until perf review
motivates it.

### Step 4 — Perf review

After step 2 (or 3) lands and byte-identity is confirmed, run
[`benchmarks/tomato1/scripts/perf_scaling_synthetic.py`](../../../benchmarks/tomato1/scripts/perf_scaling_synthetic.py)
on N=50 / 200 / 1000 vs the layer-4 baseline. Expected: per-group
allocation count drops, wall time drops modestly (the EM
dominates the per-group cost; removing 6 small allocs is a small
win unless the allocator was a bottleneck). Update
[PROJECT_STATUS](../PROJECT_STATUS.md) with the final numbers.

## Validation

Each step ships with:

- Unit tests against the prior step's EM output on the same
  fixtures as Phase A.1 layer 4 (no-compound, multi-position,
  with-compound, cap-fires).
- The cohort CLI integration tests
  ([`tests/cohort_cli_integration.rs`](../../../tests/cohort_cli_integration.rs))
  must still pass byte-identically.
- The 3-tomato real-data byte-identity diff vs the layer 4
  baseline must remain 0 lines.

## Risks

- **`PosteriorEngine` is generic over its `MathBackend`.** The
  EM body is monomorphised four ways (`ExactMath`,
  `InterpUnivariateMath`, `InterpUnivariateSimdMath`,
  test-only backends). Step 1's field-access translation must
  preserve this generic structure.
- **`run_em_for_record` is `fn`-private and is called from one
  site (`PosteriorEngine::process`).** No external API breakage,
  but the per_group_merger module's row-shape consumers (tests +
  fuzz harnesses) need to keep building `MergedRecord` for backward
  compatibility. The `MergedRecordAllelesView` shim handles this.
- **Mixture pre-pass and posterior summarisation pass through the
  same EmInputs.** Both consume `scalars` and `chain_anchor_flags`
  by slice; the signature refactor extends to them uniformly.
