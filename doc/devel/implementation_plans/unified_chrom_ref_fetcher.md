# Unified `ChromRefFetcher` trait

Implementation plan to collapse the family of reference-FASTA fetcher
types in [`src/per_sample_pileup/ref_fetcher.rs`](../../src/per_sample_pileup/ref_fetcher.rs)
behind a single narrow trait. Follow-up to the
[reference_fasta_streaming](reference_fasta_streaming.md) plan
(Phases A / B / C, commits `839af2e` / `eac3843` / `b703a8c`).

The earlier phases shrank peak memory in the cohort var-calling path
from 4.86 GB → 3.84 GB on the tomato fixture by replacing the
"cache every contig" fetcher with a per-worker sliding-buffer
streamer. This left the codebase with several fetcher types that
each bake their access-pattern assumptions into the *type* (evicting
vs accumulating; single-chrom vs multi-chrom; Sync vs !Sync) and one
type that's outright dead code.

This plan unifies them behind one contract and uses runtime errors —
not type-level signaling — to surface access-pattern violations.

## Background — today's fetcher menu

Four fetcher types live in
[ref_fetcher.rs](../../src/per_sample_pileup/ref_fetcher.rs) today:

| Type | Lines | Used by | Memory | Threading | Multi-chrom |
|---|---:|---|---|---|---|
| `ChromBoundaryRefFetcher` | [L52-L102](../../src/per_sample_pileup/ref_fetcher.rs#L52) | Stage 1 walker | 1 contig (evicted on chrom change) | !Sync | yes |
| `SyncRefFetcher` | [L131-L150](../../src/per_sample_pileup/ref_fetcher.rs#L131) | Stage 1 BAQ, `var_calling_from_bam` | All visited contigs accumulated | Sync | yes |
| `SingleChromRefFetcher` | [L254-L298](../../src/per_sample_pileup/ref_fetcher.rs#L254) | **dead code** | (was: full contig) | Sync | no |
| `StreamingChromRefFetcher` | [L371-L639](../../src/per_sample_pileup/ref_fetcher.rs#L371) | cohort var-calling `process_one_chromosome` | 1 MB sliding buffer | Sync | no |

All four implement the `RefSeqFetcher` trait from
[`pileup::mod`](../../src/per_sample_pileup/pileup/mod.rs#L510):

```rust
pub trait RefSeqFetcher {
    fn fetch(&self, chrom_id: u32, start_1based: u32, length: u32) -> io::Result<Vec<u8>>;
    fn iter_bases<'a>(
        &'a self,
        chrom_id: u32,
        length: u32,
    ) -> io::Result<Box<dyn Iterator<Item = io::Result<u8>> + 'a>>;
}
```

Two design smells in this menu:

1. **Type-encoded contracts.** The choice of fetcher tells you the
   access shape: single-chrom fetchers panic on a foreign `chrom_id`;
   evicting fetchers panic if you go back to an evicted contig. The
   contract is enforced by *which type you picked at construction*,
   not by anything the consumer code makes visible.
2. **`chrom_id` is half-dead.** Single-chrom fetchers validate the
   `chrom_id` argument and error on mismatch. It's been a "we'll
   check you pass the right one" parameter since Phase B. The
   consumer's chrom_id is implicit in *which fetcher instance* they
   hold; passing it through the API is redundant.

Both observations point to the same refactor: a narrower trait,
fewer impls, runtime errors for contract violations.

## Goal

One trait per access pattern, with all consumers held to the same
shape. Specifically:

- **One `ChromRefFetcher` trait** with a contig-bound API
  (`fetch(start, length)` — no `chrom_id`).
- **One `StreamingChromRefFetcher` impl** for the
  monotonic-forward-with-small-lookback access pattern. Errors with
  `OutOfPattern` on a backward jump past its sliding buffer — the
  signal that some consumer's access pattern doesn't fit.
- **A `WalkerFetcher` wrapper** for multi-chrom consumers (the
  Stage 1 walker, `var_calling_from_bam`). Holds a swappable
  `ChromRefFetcher` inside and rebuilds it when `fetch` sees a new
  contig.
- **No random-access impl yet.** We add one only if a real consumer
  hits the `OutOfPattern` error.

## Spec / supporting documents

- Earlier phases of the broader effort:
  [reference_fasta_streaming.md](reference_fasta_streaming.md)
  (Phase A: streaming MD5 verify; Phase B: per-worker
  `SingleChromRefFetcher`; Phase C: sliding-buffer
  `StreamingChromRefFetcher`).
- Perf review that surfaced the L4 / L5 / L6 fetcher concerns:
  [perf_psp_to_vcf_2026-05-20.md](../reports/reviews/perf_psp_to_vcf_2026-05-20.md)
  §5 L4-L6.
- Pipeline architecture spec:
  [calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md).

## API surface

The single trait, replacing today's `RefSeqFetcher`:

```rust
/// Reference fetcher bound to one contig at construction time.
/// Implementations differ in access-pattern guarantees:
///
/// - [`StreamingChromRefFetcher`] requires monotonic-non-decreasing
///   `start` across `fetch` calls (small look-back within the
///   current sliding buffer is fine). Out-of-pattern access returns
///   [`ChromRefFetchError::OutOfPattern`].
/// - A future random-access impl (built only if a real consumer
///   needs it) would accept any access pattern with higher memory
///   cost.
pub trait ChromRefFetcher {
    /// Total bases in the bound contig.
    fn length(&self) -> u32;

    /// Return `length` uppercased bases starting at 1-based `start`.
    /// Implementations promise SAM-spec canonical bytes
    /// (`A`/`C`/`G`/`T`/`N`), uppercased regardless of how the FASTA
    /// is masked on disk.
    fn fetch(&self, start_1based: u32, length: u32) -> Result<Vec<u8>, ChromRefFetchError>;

    /// Forward sequential iterator over every uppercased base of
    /// the contig, in 1..=length order. Used by DUST's mask
    /// construction (one pass per chrom). Yields per-byte
    /// `Result<u8, ChromRefFetchError>` so refill / I/O errors
    /// can propagate without latching the iterator early.
    fn iter_bases<'a>(
        &'a self,
    ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>;
}

/// Why a `ChromRefFetcher` call failed.
#[derive(Debug, thiserror::Error)]
pub enum ChromRefFetchError {
    /// The requested `[start, start+length)` exceeds the contig's
    /// total length.
    #[error("fetch [{start}, {end}) past contig {contig_name} length {contig_length}")]
    OutOfBounds {
        contig_name: String,
        contig_length: u32,
        start: u32,
        end: u32,
    },
    /// `start_1based` was 0 (1-based coordinate contract).
    #[error("start_1based must be >= 1")]
    InvalidStart,
    /// Streaming impls only: the requested range lies before the
    /// current sliding buffer, i.e. the consumer's access pattern
    /// is not monotonic-non-decreasing. The error message names
    /// the buffer's current origin so the caller can see how far
    /// the streamer has advanced.
    #[error(
        "fetch at base {requested_start} lies before the streamer's current buffer (origin = {buffer_origin}). \
         StreamingChromRefFetcher requires monotonic-non-decreasing access; \
         a random-access impl is needed for this consumer."
    )]
    OutOfPattern {
        requested_start: u32,
        buffer_origin: u32,
    },
    /// Underlying I/O failure (file read, seek, .fai parse, …).
    #[error("I/O failure on contig {contig_name}: {source}")]
    Io {
        contig_name: String,
        #[source]
        source: io::Error,
    },
}
```

Constructor signatures (free functions on each impl, not on the
trait — `Self: Sized` makes `dyn ChromRefFetcher` constructible
otherwise impossible):

```rust
impl StreamingChromRefFetcher {
    /// Open the FASTA at `fasta_path`, look up `contig_name` in the
    /// sibling `<fasta_path>.fai`, and return a fetcher ready to
    /// serve `fetch`. Contig bytes are loaded lazily on the first
    /// `fetch` / `iter_bases` call.
    pub fn new(fasta_path: &Path, contig_name: &str) -> Result<Self, ChromRefFetchError>;

    /// Variant constructor for non-standard `.fai` locations.
    pub fn with_fai_path(
        fasta_path: &Path,
        fai_path: &Path,
        contig_name: &str,
    ) -> Result<Self, ChromRefFetchError>;
}
```

Multi-chrom consumers (walker, from-bam) use:

```rust
/// Adapter that turns a multi-chrom access pattern into a sequence
/// of single-contig [`ChromRefFetcher`] instances. The inner
/// fetcher is swapped when a `fetch_in(chrom_id, ...)` call sees a
/// different `chrom_id` than the currently-bound contig.
///
/// Used by the Stage 1 pileup walker (which visits contigs in
/// order) and by `var_calling_from_bam` (single-sample walker into
/// the cohort pipeline). The walker's access pattern is sequential
/// within a contig plus a strict-forward chrom transition; matches
/// `StreamingChromRefFetcher`'s contract.
pub struct WalkerFetcher {
    fasta_path: PathBuf,
    contigs: ContigList,
    inner: RefCell<Option<(u32, StreamingChromRefFetcher)>>,
}

impl WalkerFetcher {
    pub fn new(fasta_path: PathBuf, contigs: ContigList) -> Self;

    /// Fetch a window on `chrom_id`. Rebuilds the inner fetcher on
    /// chrom transition; subsequent fetches on the same chrom hit
    /// the streamer's sliding buffer.
    pub fn fetch_in(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, ChromRefFetchError>;

    pub fn iter_bases_in<'a>(
        &'a self,
        chrom_id: u32,
    ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>;
}
```

`WalkerFetcher` is **not** itself a `ChromRefFetcher` (the trait is
single-contig). It exposes its own `fetch_in` / `iter_bases_in` with
the chrom_id back in the signature. Walker code keeps a similar
shape to today's calls — just with the new types.

For backwards compat during migration, the old `RefSeqFetcher` trait
stays in place until step 3. New code uses `ChromRefFetcher` or
`WalkerFetcher` directly; we delete the old trait once nothing
imports it.

## Migration order

Discovery-oriented: implement the new trait + impl in step 1, then
migrate consumers one at a time. If a migration surfaces a real
non-monotonic access pattern, we pause and add a random-access impl
before continuing.

```
Step 0 — cleanup
   delete dead SingleChromRefFetcher + fix stale doc references
   (no behaviour change, no test changes)
        │
Step 1 — implement new trait + StreamingChromRefFetcher (new API)
   live alongside today's RefSeqFetcher + StreamingChromRefFetcher
   (old API). No consumer migrated yet. Unit tests for the new
   contract (including a deliberate OutOfPattern firing test).
        │
Step 2(a) — migrate cohort var-calling
   process_one_chromosome switches to the new ChromRefFetcher.
   DUST + PerGroupMerger trait-object usage updates accordingly.
   Mock fetchers in dust_filter tests migrate to the new trait.
   Phase C's "StreamingChromRefFetcher (old API)" becomes the only
   consumer of the old trait that's left to migrate.
        │
Step 2(b) — migrate walker + var_calling_from_bam
   New WalkerFetcher wrapper. Walker code drops
   ChromBoundaryRefFetcher and SyncRefFetcher (for the from-bam
   use). The Stage 1 pipeline_fetchers tuple shrinks accordingly.
        │
Step 2(c) — migrate BAQ
   Per-rayon-worker StreamingChromRefFetcher construction. The
   "shared SyncRefFetcher across rayon workers" pattern in BAQ
   becomes "each worker constructs its own fetcher per batch."
   Removes the last SyncRefFetcher use site.
        │
Step 3 — final cleanup
   Delete the old RefSeqFetcher trait. Delete SyncRefFetcher,
   ChromBoundaryRefFetcher, and the now-unused old-API
   StreamingChromRefFetcher. Module shrinks substantially.
```

Each step is one commit; the test suite stays green at each
intermediate commit. We can stop after any step and either continue
later or revert without entanglement.

## Step-by-step detail

### Step 0 — cleanup (no behaviour change)

**Scope:**

- Delete `SingleChromRefFetcher` struct + impl + 10 tests
  ([ref_fetcher.rs:228-316](../../src/per_sample_pileup/ref_fetcher.rs#L228-L316)
  + the `SingleChromRefFetcher tests` test block at
  [ref_fetcher.rs:943-](../../src/per_sample_pileup/ref_fetcher.rs#L943)).
- Fix the stale `SingleChromRefFetcher` doc references in
  [var_calling.rs:230](../../src/pop_var_caller/var_calling.rs#L230),
  [var_calling.rs:338](../../src/pop_var_caller/var_calling.rs#L338),
  [var_calling.rs:383](../../src/pop_var_caller/var_calling.rs#L383),
  [cohort_driver.rs:436](../../src/pop_var_caller/cohort_driver.rs#L436):
  replace with `StreamingChromRefFetcher` references.
- Fix the `StreamingChromRefFetcher` doc paragraph
  ([ref_fetcher.rs:353,362](../../src/per_sample_pileup/ref_fetcher.rs#L353))
  that compares to `SingleChromRefFetcher` (the comparison is now
  pointless).

**Tests:** existing tests continue to pass. Lib test count drops by
the ~10 tests we delete.

**Acceptance:** 868 lib + 41 integration + 4 doc tests pass; clippy
clean on the touched files.

### Step 1 — new trait + impl, no consumer migration

**Scope:**

- New `ChromRefFetcher` trait in
  [ref_fetcher.rs](../../src/per_sample_pileup/ref_fetcher.rs)
  alongside the old `RefSeqFetcher` (which stays for now).
- New `ChromRefFetchError` enum (uses `thiserror`).
- New `StreamingChromRefFetcher` (or rename today's to
  `LegacyStreamingChromRefFetcher` and let the new one own the
  short name — naming decision in the impl). The new one has the
  new API and the new error type.
- The internal buffer-management code (the `refill`,
  `base_to_file_offset`, the `Source` enum, the `StreamState`
  struct) can be shared between old and new fetchers via free
  functions; both call the same `refill` helper.
- The new fetcher includes the `OutOfPattern` runtime check on
  every fetch: if `start_1based < buf_start_base`, return
  `Err(OutOfPattern { … })`. Forward refills happen transparently
  (no error).

**Design decision on `iter_bases` + `fetch` composition.**

DUST calls `iter_bases()` once at chrom load and walks the whole
contig forward; the sliding buffer ends positioned near the contig's
end. PerGroupMerger then calls `fetch(start, length)` starting near
the first variant group's position, which is typically near the
*start* of the contig. If `OutOfPattern` fired on
`start < buf_start_base`, it would trigger on PerGroupMerger's very
first call after DUST.

**Resolution:** `iter_bases()` resets the sliding buffer on
construction and on Drop. The buffer is then empty when
PerGroupMerger's first `fetch` arrives, which refills from its own
`start` and updates the watermark to that origin. Subsequent
`fetch` calls must be monotonic non-decreasing relative to that
origin; backward jumps beyond the buffer's origin fire
`OutOfPattern`.

Conceptually: `iter_bases()` is a one-shot phase that owns the
buffer for its lifetime; `fetch()` is a separate phase that starts
fresh. The trait doc states this contract explicitly. No extra
memory cost (the buffer is reused, not duplicated). Alternative
designs considered: (a) separate buffer for `iter_bases` (doubles
per-worker memory); (b) drop `OutOfPattern` entirely (loses the
fail-fast contract). Going with the reset-on-Drop approach because
it costs nothing and keeps the contract crisp.

**Tests:**

All new unit tests live in `ref_fetcher.rs::tests`:

1. `chrom_ref_returns_bytes_for_bound_contig` — basic positive path.
2. `chrom_ref_serves_line_wrapped_contig` — 60-col wrap (Ensembl
   default).
3. `chrom_ref_uppercases_soft_masked` — uppercasing contract.
4. `chrom_ref_fetch_past_contig_end_returns_out_of_bounds` — bound
   check.
5. `chrom_ref_fetch_start_zero_is_rejected` — 1-based contract.
6. `chrom_ref_construct_missing_contig_errors` — `.fai` lookup
   failure.
7. `chrom_ref_buffer_refill_on_forward_jump` — large forward jump
   triggers refill, returns correct bytes.
8. `chrom_ref_small_backward_within_buffer` — small look-back
   inside the current 1 MB stays within buffer, no error.
9. **`chrom_ref_backward_beyond_buffer_returns_out_of_pattern`** —
   build a fetcher with a deliberately small buffer (4 KB via a
   test-only constructor), advance past it, do a backward fetch.
   Assert `OutOfPattern` fires with the right field values. Pins
   the contract.
10. `chrom_ref_iter_bases_yields_uppercased_contig_in_order` —
    DUST mask path.
11. `chrom_ref_against_real_file` — smoke test using a real FASTA
    + .fai built via `build_fasta`.
12. `chrom_ref_with_fai_path_alternate_location` — non-standard
    `.fai` location.

**Acceptance:** all existing tests pass; ~12 new tests pass; clippy
clean. The new fetcher is unused at this point.

### Step 2(a) — migrate cohort var-calling

**Scope:**

- [`process_one_chromosome`](../../src/pop_var_caller/cohort_driver.rs#L466)
  constructs the new `StreamingChromRefFetcher` (new API) instead
  of the old one. Drops the `chrom_id` argument from `fetch` calls.
- DUST consumes the new trait. Today `DustFilter` is generic over
  `F: RefSeqFetcher` (old trait); it stays generic over a trait
  but the bound becomes the new one. The `chrom_id` parameter in
  `ensure_mask_for` is no longer threaded through to the fetcher —
  the fetcher knows its bound chrom, and DUST's own
  `chromosomes`-table lookup serves the contig-name and length
  for `chrom_id` resolution. Equivalent change in
  `dust_filter.rs::DustFilter::iter_bases` callsite.
- PerGroupMerger fetches similarly drop `chrom_id`.
- The `SharedRefFetcher` type alias becomes
  `Arc<dyn ChromRefFetcher + Send + Sync>`.
- Mock fetchers in
  [`dust_filter.rs::tests::StubFetcher`](../../src/var_calling/dust_filter.rs#L1060)
  and in `per_group_merger.rs::tests` migrate to the new trait.
  Mostly mechanical — drop the `chrom_id` parameter, return the
  same bytes.
- Old `StreamingChromRefFetcher` (with the chrom_id-aware API)
  becomes unused by `var_calling`. We leave it in place — walker,
  from-bam, and BAQ still need the old trait until later steps.

**Tests:** existing dust_filter unit tests (38), per_group_merger
unit tests, and cohort integration tests (the existing
`cohort_cli_integration.rs` 9 tests) all pass. No new tests in this
step — the contract is unchanged, only the API plumbing changes.

**Acceptance:** all existing tests pass; lib tests stable; cohort
integration tests stable; clippy clean.

### Step 2(b) — migrate walker + var_calling_from_bam

**Scope:**

- New `WalkerFetcher` wrapper in
  [ref_fetcher.rs](../../src/per_sample_pileup/ref_fetcher.rs).
  Internal state: `fasta_path`, `contigs: ContigList`,
  `inner: RefCell<Option<(u32, StreamingChromRefFetcher)>>`.
- `WalkerFetcher::fetch_in(chrom_id, start, length)`:
  - If `inner` is `Some(c, _)` and `c == chrom_id`: delegate to the
    inner fetcher's `fetch(start, length)`.
  - Otherwise: construct a new `StreamingChromRefFetcher` for
    `chrom_id` via `contigs.entries[chrom_id].name`, replace
    `inner`, then delegate.
- `WalkerFetcher::iter_bases_in(chrom_id)`: same shape.
- `stage1_pipeline.rs` swaps the two-fetcher pair
  (`baq_fetcher: SyncRefFetcher` + `walker_fetcher: ChromBoundaryRefFetcher`,
  [stage1_pipeline.rs:111-114](../../src/pop_var_caller/stage1_pipeline.rs#L111-L114))
  for a single `WalkerFetcher`. Wait — BAQ uses the BAQ fetcher
  in a *parallel* context (rayon workers share it). `WalkerFetcher`
  has `RefCell` interior mutability and is `!Sync`. So Stage 1 in
  step 2(b) keeps `SyncRefFetcher` for BAQ and only migrates the
  walker fetcher. BAQ migration is step 2(c).
- `var_calling_from_bam.rs::run_var_calling_from_bam`
  ([L354](../../src/pop_var_caller/var_calling_from_bam.rs#L354))
  constructs a `WalkerFetcher` instead of `SyncRefFetcher`. The
  cohort pipeline downstream of the walker is single-thread in
  this code path (no per-chrom rayon), so `!Sync` is fine.

  Subtlety: from-bam passes the fetcher to `drive_cohort_pipeline`
  as `SharedRefFetcher = Arc<dyn RefSeqFetcher + Send + Sync>`.
  After step 2(a) this becomes `Arc<dyn ChromRefFetcher + Send + Sync>`.
  But `WalkerFetcher` isn't a `ChromRefFetcher` — it's a separate
  type. We need a small adapter that wraps `WalkerFetcher` into
  the trait, hiding the chrom_id swap. Either:
  - **(i)** Make from-bam's path build a wrapper that implements
    `ChromRefFetcher` by routing to a single inner `WalkerFetcher`
    (always passing the same chrom_id, since from-bam is
    single-chrom-at-a-time downstream of the walker — wait, this
    needs verification).
  - **(ii)** Or have `drive_cohort_pipeline` take a multi-chrom
    fetcher, parameterizing it differently.

  I'll defer the exact shape to implementation time after reading
  the from-bam code more carefully. **This is one of the few
  places where the migration could surprise us; flag it during
  implementation.**

- Walker tests: `ChromBoundaryRefFetcher` has its own test set in
  [ref_fetcher.rs::tests](../../src/per_sample_pileup/ref_fetcher.rs#L807)
  (cache eviction invariant, etc.). Those tests become tests
  for `WalkerFetcher` with equivalent behaviour assertions.
- Stage 1 integration tests
  ([pileup_cli_integration.rs](../../tests/pileup_cli_integration.rs))
  pass unchanged.

**Acceptance:** all existing tests pass (lib + Stage 1 integration);
no perf regression on the Stage 1 walker (the inner fetcher's I/O
shape is now sliding-buffer instead of `noodles::Repository`
materialise-once, but the chrom-eviction shape is preserved).

### Step 2(c) — migrate BAQ

**Scope:**

- BAQ today gets the fetcher via the `pipeline_fetchers` tuple from
  [`with_stage1_pipeline`](../../src/pop_var_caller/stage1_pipeline.rs#L89).
  The BAQ engine
  ([per_sample_pileup/baq/engine.rs:107](../../src/per_sample_pileup/baq/engine.rs#L107))
  takes `&dyn RefSeqFetcher` per call. Reads come in batches; rayon
  parallelises the batches. The fetcher is shared (`Sync`).
- New design: BAQ's per-batch processing constructs a fresh
  `StreamingChromRefFetcher` per rayon worker. Each worker holds
  its own fetcher for the duration of the batch; batches are
  per-contig, so the chrom_id is known.
- The change ripples through BAQ's parallelism layer. Need to read
  the BAQ batch driver to see where the rayon split happens and
  how to thread the per-worker fetcher in.
- Open question: what's the lifetime of a "batch"? If a batch
  spans multiple chroms, we'd need a `WalkerFetcher`-style swap
  inside each worker. If batches are always single-chrom, simpler.
  **Confirm at implementation time by reading BAQ's driver.**
- Remove the `SyncRefFetcher` field from `pipeline_fetchers`
  (Stage 1's BAQ fetcher slot). After this step, no production
  code holds a `SyncRefFetcher`.

**Tests:** existing BAQ unit tests + Stage 1 integration tests
pass. The "shared fetcher across rayon workers" semantics change to
"per-worker fetcher" — if a BAQ test currently asserts on shared
state, it'll need to adapt; otherwise mechanical.

**Acceptance:** all tests pass; no perf regression on Stage 1 BAQ.

### Step 3 — final cleanup

**Scope:**

- Delete the old `RefSeqFetcher` trait (`fetch(chrom_id, …)` and
  `iter_bases(chrom_id, …)`) from
  [per_sample_pileup/pileup/mod.rs](../../src/per_sample_pileup/pileup/mod.rs#L505-L555).
- Delete `SyncRefFetcher`
  ([ref_fetcher.rs:131-150](../../src/per_sample_pileup/ref_fetcher.rs#L131-L150))
  + its impl + tests.
- Delete `ChromBoundaryRefFetcher`
  ([ref_fetcher.rs:52-102](../../src/per_sample_pileup/ref_fetcher.rs#L52-L102))
  + its impl + tests, if step 2(b) successfully replaced it.
- Delete the old-API `StreamingChromRefFetcher` (the
  pre-step-1-rename one), keeping only the new-API one with the
  short name.
- Drop the `fetch_from_repository` shared helper if both legacy
  callers are gone.
- The `noodles_fasta::Repository` runtime dependency disappears
  from production code (still used at construction time by the
  streaming impl's `noodles_fasta::io::indexed_reader` path; no
  cache survives).

**Tests:** all existing tests pass.

**Acceptance:** lib test count net change is roughly neutral
(removed legacy fetcher tests are replaced by new fetcher tests
during steps 1, 2b); the module's line count drops substantially.

## Test strategy

Each step is one commit and ships green. Specifically:

- **Step 0** — existing tests, minus the ~10 deleted dead-code
  tests, all pass.
- **Step 1** — adds ~12 unit tests on the new trait + impl. The
  critical one is the `OutOfPattern` test (#9 above), which pins
  the contract.
- **Steps 2(a)/2(b)/2(c)** — each migrates existing tests as needed
  but adds no new ones. The contract is unchanged; only the API
  plumbing changes.
- **Step 3** — net test count drops as legacy tests are deleted;
  coverage of fetcher behaviour stays the same because step 1
  added equivalent coverage on the new impl.

**Real-data validation.** After step 2(a) lands, re-run the
`tmp/measure_rss.sh` benchmark on the tomato fixture. Peak RSS
should be identical to Phase C (3.84 GB ± noise); wall should be
unchanged. The migration is an API change, not an algorithm change.

After step 2(b), the Stage 1 walker's memory usage should drop
by approximately one contig's worth — the streaming buffer replaces
the full-contig cache. On tomato that's ~90 MB; on human chr1 it's
~250 MB.

After step 2(c), Stage 1's BAQ memory drops by `Σ contig.length` (the
accumulating SyncRefFetcher's cache disappears). Per-worker overhead
is one 1 MB streaming buffer instead.

## Risks and discovery items

1. **A consumer turns out to have a non-monotonic access pattern.**
   Likely candidates: BAQ (batches might span boundaries
   unexpectedly), the walker (look-backs > 1 MB on large indels).
   Mitigation: the `OutOfPattern` error tells us exactly which
   consumer + which transition triggered it. Response: either
   constrain the consumer (preferred) or build a random-access
   impl with the same trait (acceptable; same drop-in shape).

2. **`var_calling_from_bam`'s downstream type plumbing.** The
   cohort driver expects a `SharedRefFetcher = Arc<dyn ... + Send + Sync>`.
   `WalkerFetcher` is `!Sync` (RefCell). We need a `Sync`
   adapter wrapping `WalkerFetcher`, or the from-bam path needs a
   different trait alias. **To be resolved at step 2(b)
   implementation time.**

3. **BAQ's rayon worker structure.** The exact granularity of BAQ
   batches matters for whether per-worker fetcher construction
   is cheap (fai parse + file open per batch). If batches are
   small and frequent, we may need a `FetcherFactory` that holds
   the pre-parsed fai and hands out fetchers cheaply.
   **To be resolved at step 2(c) implementation time.**

4. **Test mock churn.** Each consumer's tests have their own
   mock fetcher. Each one needs a 5-10 line update per migration
   step. Tedious but mechanical; not a design risk.

5. **API surface during transition.** Old and new traits coexist
   for steps 1 → 2(c). Documentation in `ref_fetcher.rs` should
   make the lifecycle clear: "old: being removed in step 3; new:
   use for any new code." Avoid having both reachable from the
   same crate root.

## Out of scope (explicitly)

- **Random-access impl.** We only build one if a step 2 migration
  uncovers a real non-monotonic consumer.
- **Compile-time access-pattern marker traits** (e.g. a
  `MonotonicFetcher` marker that DUST requires). Considered;
  rejected because it doubles the type menu we're trying to
  shrink. Runtime errors are fine for a contract that's verifiable
  at integration-test time.
- **Refactoring how test mocks share fixtures across modules.**
  Mocks stay private to their consumers' test blocks.
- **Adding multi-thread access to a single `StreamingChromRefFetcher`
  instance.** Phase B's design is one fetcher per worker; that
  stays.
- **Persisting the streaming buffer across `WalkerFetcher` chrom
  swaps.** When the walker moves chroms, the old streamer is
  dropped wholesale. Reusing the buffer allocation across swaps
  is a possible micro-opt; we'd want a profile to motivate it.

## Estimated effort

- **Step 0**: ~50 lines deleted, no behaviour change. Minutes.
- **Step 1**: ~400 lines (new trait + new impl + ~12 unit tests).
  Half a day.
- **Step 2(a)**: ~250 lines (cohort driver + DUST + per-group
  merger plumbing + mock updates). Half a day, mostly mechanical.
- **Step 2(b)**: ~250 lines (`WalkerFetcher` + walker migration +
  from-bam migration + walker test migration + Sync-adapter
  decision for the from-bam path). Day.
- **Step 2(c)**: ~250 lines (BAQ per-worker construction; needs
  reading BAQ's driver carefully). Day.
- **Step 3**: ~200 lines deleted (old trait + three old impls +
  their tests). Hours.

Total: ~1400 net lines across 6 commits, ~3-4 days of work
spread out. Each commit is independently shippable.
