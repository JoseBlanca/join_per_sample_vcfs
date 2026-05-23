# Reference FASTA streaming — cut memory in `var-calling`

Implementation plan for two related findings in the
[2026-05-20 perf review](../reports/reviews/perf_psp_to_vcf_2026-05-20.md):
**H5** (M5 verify streaming) and **L5** (`SyncRefFetcher` re-shape).
The user has explicitly asked for both: first the startup spike, then
the steady-state memory used by the cohort `var-calling` workers. The
refactor in Phase B is somewhat invasive; the user has signed off on
that scope and there is no schedule pressure.

## Today's memory shape (the diagnosis)

Walking the code for `pop_var_caller var-calling`:

1. [`SyncRefFetcher`](../../src/per_sample_pileup/ref_fetcher.rs#L109-L137)
   wraps `noodles_fasta::Repository`, which holds a
   `HashMap<Vec<u8>, Arc<Sequence>>` cache with **no eviction
   policy**. Its own doc comment ([ref_fetcher.rs:103-108](../../src/per_sample_pileup/ref_fetcher.rs#L103-L108))
   states the steady-state footprint is "every chromosome the run
   visits (~3 GB on a 24-chrom human reference)".
2. The startup step
   [`verify_fasta_matches_psp_chromosomes`](../../src/pop_var_caller/common.rs#L183-L208)
   calls `fetcher.fetch(chrom_id, 1, contig.length)` for **every**
   chromosome to compute its MD5 and compare against the `.psp`
   header. Each call allocates a fresh contig-sized `Vec<u8>`,
   uppercases it byte-by-byte, then feeds it to `Md5::digest`. The
   doc comment ([common.rs:174-182](../../src/pop_var_caller/common.rs#L174-L182))
   is explicit that it also acts as a cache-warmer — by the time the
   `par_iter` over chromosomes begins, **every contig is already
   resident**.
3. Once the parallel work starts, each worker's
   [`DustFilter::ensure_mask_for`](../../src/var_calling/dust_filter.rs#L656-L683)
   and
   [`PerGroupMerger`'s `ref_fetcher.fetch(...)`](../../src/var_calling/per_group_merger.rs#L704)
   read from that already-warm cache. No worker ever drops its
   contig, so the cache stays at all-chroms-resident for the rest of
   the run.

There are therefore two distinct sources of pressure:

- **Startup spike (Phase A target):** the verify loop holds a
  contig-sized buffer + `Repository` cache copy simultaneously for
  each visited contig. On the tomato fixture this is two ~91 MB
  allocations per visited contig in sequence; on a human genome it
  is two ~250 MB allocations for chr1.
- **Steady-state ceiling (Phase B target):** the `Repository` cache
  is never evicted, so peak resident memory during the parallel work
  equals `Σ contig.length` across the whole reference (~3 GB on
  human; ~1.2 GB on the SL4.0 tomato reference).

## Spec / supporting documents

- Perf review (verdict + measurement plan):
  [perf_psp_to_vcf_2026-05-20.md](../reports/reviews/perf_psp_to_vcf_2026-05-20.md)
  §5 H5 (lines 125-131) and §5 L5 (lines 190-196).
- Pipeline architecture spec:
  [calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md).
- The per-chrom parallel plan that exposed L5 as the next ceiling:
  [cohort_per_chromosome_parallel.md](cohort_per_chromosome_parallel.md)
  §"Open work" #2 (the L5 follow-up note).
- Existing single-chrom eviction precedent:
  [`ChromBoundaryRefFetcher`](../../src/per_sample_pileup/ref_fetcher.rs#L30-L83)
  — `!Sync` but proves out the chrom-boundary-clear pattern;
  Phase B is its `Sync` cousin.

## Why this is two phases, not one

The two findings are coupled but separable:

- Phase A is a self-contained ~150-line change. It removes the
  startup spike, parallelises the slowest startup step, and — most
  importantly — **decouples verify from cache state**. Today the
  verify step pre-warms the `Repository`; while that's the case the
  fetcher can't evict at chrom boundaries (worker 0 starts before
  worker 12 has fetched, and worker 0 eviction would have to know
  not to evict 12's chrom-of-interest). After Phase A the cache is
  empty at the start of `par_iter`, and Phase B's eviction logic
  has only one source of cached bytes to manage.
- Phase B changes the fetcher type itself + every construction site.
  It's the larger blast radius. Shipping it after Phase A means the
  Phase B PR has only one moving piece (the fetcher) and a clean
  baseline to compare against.

We could fold them. We won't, because the bisect surface and the
review surface are both smaller this way.

## Scope

### Phase A — Streaming MD5 verify

**In scope:**

- New function `compute_contig_md5_streaming` in
  [src/pop_var_caller/common.rs](../../src/pop_var_caller/common.rs)
  (or a new sibling module) that takes a FASTA path + a
  `ParsedChromosome` (carries name + length + expected MD5) and
  returns the lowercase-hex MD5 of the contig's uppercase bases —
  computed without holding the contig's full byte buffer.
- Replace the body of
  [`verify_fasta_matches_psp_chromosomes`](../../src/pop_var_caller/common.rs#L183-L208)
  with a `chromosomes.par_iter().enumerate()` that calls the new
  function per contig, comparing against the expected MD5.
- The function no longer needs `&impl RefSeqFetcher`; it takes the
  FASTA path directly. The signature changes from
  `verify_fasta_matches_psp_chromosomes(fetcher, chromosomes)` to
  `verify_fasta_matches_psp_chromosomes(fasta_path, chromosomes)`.
  The three callers
  ([var_calling.rs:316](../../src/pop_var_caller/var_calling.rs#L316),
  [estimate_contamination.rs:391](../../src/pop_var_caller/estimate_contamination.rs#L391),
  and the cleanup in [var_calling_from_bam.rs](../../src/pop_var_caller/var_calling_from_bam.rs))
  drop their `let verify_fetcher = SyncRefFetcher::new(...);` line.
- The doc comment at [common.rs:174-182](../../src/pop_var_caller/common.rs#L174-L182)
  is rewritten: the "pre-warms the cache" coupling no longer
  exists; bytes are streamed through MD5 and dropped.
- Unit tests in the same module: golden MD5 of a fixture FASTA
  against pre-computed digests (lowercase + soft-masked + line-wrapped
  fixtures); mismatch detection; missing-contig error path.

**Out of scope (deferred to Phase B):**

- Any change to `SyncRefFetcher` or `RefSeqFetcher`. The fetcher
  still pre-warms via the parallel DUST + per-group merger
  fetches — Phase A just stops doing it at startup.
- Dropping the `noodles_fasta::Repository` runtime dependency.

### Phase B — Per-worker fetcher (steady-state)

> **Design note added during implementation:** the plan originally
> described a single `ChromLeaseSyncRefFetcher` shared across
> workers via `Arc`, with per-`chrom_id` `Mutex<Slot>` and RAII
> `ChromLease`. The user pointed out that since each worker only
> touches its own chrom, the shared-Arc shape was unnecessary
> machinery — per-thread fetchers with no shared state are simpler
> and have zero possibility of cross-thread interference. The
> shipped Phase B uses `SingleChromRefFetcher` (one chrom, owned
> outright by the worker, no Mutex, no Arc-sharing). The plan
> below is preserved as written; the **shipped shape** is documented
> at the top of [src/per_sample_pileup/ref_fetcher.rs](../../src/per_sample_pileup/ref_fetcher.rs).

**In scope (as shipped):**

- New `SingleChromRefFetcher` in
  [src/per_sample_pileup/ref_fetcher.rs](../../src/per_sample_pileup/ref_fetcher.rs)
  alongside the existing two fetchers. Implements
  [`RefSeqFetcher`](../../src/per_sample_pileup/pileup/mod.rs).
  Bound to **one** `chrom_id` at construction; owns the contig's
  uppercased `Vec<u8>` outright; `fetch` calls for any other
  chrom_id are rejected with `InvalidInput`. No shared cache, no
  Mutex, no lease.

**In scope (original plan — superseded; kept for context):**

- New `ChromLeaseSyncRefFetcher` (working name) in
  [src/per_sample_pileup/ref_fetcher.rs](../../src/per_sample_pileup/ref_fetcher.rs)
  alongside the existing two fetchers. Implements
  [`RefSeqFetcher`](../../src/per_sample_pileup/pileup/mod.rs) **and**
  `Sync`. Holds one cache slot per `chrom_id`, lazily populated on
  first fetch by a worker, dropped explicitly via an RAII
  `ChromLease`.
- Per-chrom workers in [`process_one_chromosome`](../../src/pop_var_caller/cohort_driver.rs)
  acquire a `ChromLease` at the start and let it drop at the end.
  Lease drop releases the slot; once the lease count for a slot
  reaches zero, the cached `Arc<Vec<u8>>` is dropped and the
  contig's bytes are freed.
- Drop the runtime use of `noodles_fasta::Repository`. Reads
  go through `noodles_fasta::io::IndexedReader::query` once at
  lease acquisition (streamed into a `Vec<u8>`, uppercased in place),
  then all subsequent `fetch(chrom_id, start, len)` calls slice into
  the cached `Arc<Vec<u8>>`. (noodles' `query` returns a fully-
  buffered `Record`; that's fine here because the lease is exactly
  the "I'm about to read this contig many times" assertion the
  caller is making — one chrom in RAM is the goal.)
- `RefSeqFetcher::fetch` keeps returning `Vec<u8>` (preserves the
  trait's existing shape and downstream call sites; S2 from the
  perf review proposes a borrowed return but explicitly defers it).
- Test coverage:
  - Cache stays bounded: with two leases out across two chroms
    plus a third chrom never leased, `slot_count() == 2`.
  - Lease drop clears the slot: after the lease scope ends,
    `slot_count() == 0`.
  - Concurrent leases on the same chrom share one slot
    (count of populated slots stays at 1 even with two leases
    on the same `chrom_id`).
  - Round-trip equivalence: bytes returned by
    `ChromLeaseSyncRefFetcher` match `SyncRefFetcher` on the same
    fixture (golden test against current behaviour).

**Out of scope (Phase B, deferred):**

- **Window-streaming fetch (don't cache the whole contig).** This
  would drop peak memory further, from `T × max_chrom` to
  `T × buffer_size`, but the per-group merger does many small
  fetches per chrom; without a per-worker buffer cache it would be
  slow. Punted to a v3 if `T × max_chrom` still hurts in practice
  (it doesn't on human at T ≤ 16: 16 × 250 MB = 4 GB worst case,
  realistic average ≪ that because most chroms are small).
- **The S2 perf-review item** — change `fetch` to return
  `&[u8]` / `Arc<Vec<u8>>` to remove the per-call slice allocation.
  Wide ripple through DustFilter / PerGroupMerger / test mocks;
  perf review itself says "defer".
- **`var-calling-from-bam` and `estimate-contamination` migration.**
  Both use `SyncRefFetcher` but the per-chrom parallel pattern only
  exists in `var-calling` today
  ([cohort_per_chromosome_parallel.md "Out of scope"](cohort_per_chromosome_parallel.md#scope)).
  When those two subcommands get their own per-chrom parallel pass,
  they migrate at that time. For Phase B they keep using a
  compatibility `pub type SyncRefFetcher = ChromLeaseSyncRefFetcher`
  alias (or a thin no-lease wrapper) so their call sites don't
  change.

## Design

### Phase A — function shape

The .fai index gives each contig's byte offset in the FASTA, the
number of bases, and the line width (with and without the terminating
`\n`). Streaming MD5 needs all three: seek to the offset, read a
fixed-size buffer, strip newlines, uppercase, feed into `Md5::update`,
stop once the declared number of bases has been consumed.

```rust
/// Streaming MD5 of one contig's uppercase bases. Reads the FASTA
/// from disk in fixed windows; never holds more than `BUF_SIZE` bytes
/// of the contig in memory at once. Returns the lowercase-hex digest
/// that matches `Md5::digest(uppercase_bytes)` on the same contig.
fn compute_contig_md5_streaming(
    fasta_path: &Path,
    contig_name: &str,
) -> Result<String, FastaVerifyError> { ... }
```

Implementation sketch:

1. Open `fasta_path` + sibling `.fai` via
   `noodles_fasta::fai::Reader::read_index`.
2. Look up the `IndexRecord` for `contig_name`; bubble
   `FastaVerifyError::ContigMissing` if absent.
3. Open the FASTA file (own `File`, no shared handle) and
   `seek(SeekFrom::Start(record.offset()))`.
4. Loop with a reused 64 KiB buffer (`[u8; 65_536]` on the stack,
   not heap-alloc'd). For each read:
   - If we see a `>` at column 0 of any line, stop (we've crossed
     into the next contig — guards against fai-inconsistency).
   - Strip `\n` / `\r`.
   - `make_ascii_uppercase` on the surviving slice.
   - `md5.update(&slice)`.
   - Decrement remaining-bases counter.
5. Bail with `FastaVerifyError::TruncatedContig` if EOF arrives
   before the declared length is consumed.
6. Return `format_md5_hex(md5.finalize().into())`.

The 64 KiB buffer is on the worker's stack frame, so under
`par_iter` we pay `T × 64 KiB` of stack — negligible.

**Parallel call site:**

```rust
pub(crate) fn verify_fasta_matches_psp_chromosomes(
    fasta_path: &Path,
    chromosomes: &[ParsedChromosome],
) -> Result<(), FastaVerifyError> {
    chromosomes
        .par_iter()
        .try_for_each(|contig| {
            let fasta_md5 = compute_contig_md5_streaming(fasta_path, &contig.name)?;
            if fasta_md5 != contig.md5 {
                return Err(FastaVerifyError::Md5Mismatch {
                    contig: contig.name.clone(),
                    fasta_md5,
                    psp_md5: contig.md5.clone(),
                });
            }
            Ok(())
        })
}
```

`par_iter` uses the global rayon pool already configured by
`configure_rayon_pool` at the start of `run_var_calling`. No
additional pool setup. Each worker opens its own `File` handle (one
fd per chrom in flight, capped at `min(threads, n_chroms)`).

**Error semantics.** `par_iter().try_for_each` returns on the first
error; in-flight workers run to completion. Matches the current
fail-fast semantics (today the for-loop returns on the first contig
mismatch).

### Phase B — `ChromLeaseSyncRefFetcher`

The trait stays the same; the type changes.

```rust
pub struct ChromLeaseSyncRefFetcher {
    fasta_path: PathBuf,
    contigs: ContigList,
    /// One slot per `chrom_id`. `Mutex` protects (a) lazy
    /// initialisation and (b) the lease counter. Workers reading
    /// the same chrom share the inner `Arc<Vec<u8>>` cheaply.
    slots: Vec<Mutex<Slot>>,
}

struct Slot {
    bytes: Option<Arc<Vec<u8>>>,
    leases: u32,
}

/// RAII handle. While alive, the slot for `chrom_id` is guaranteed
/// non-empty; `fetch` calls for that chrom hit the cache. When the
/// lease drops, `leases` decrements; on zero, `bytes` is dropped and
/// the contig's memory is freed.
pub struct ChromLease<'a> {
    fetcher: &'a ChromLeaseSyncRefFetcher,
    chrom_id: u32,
}

impl ChromLeaseSyncRefFetcher {
    pub fn new(fasta_path: PathBuf, contigs: ContigList) -> io::Result<Self> { ... }

    /// Acquire a lease on `chrom_id`. First lease on an unloaded slot
    /// reads the contig from disk; subsequent leases just increment
    /// the counter.
    pub fn lease(&self, chrom_id: u32) -> io::Result<ChromLease<'_>> { ... }

    /// Test/diagnostic only — number of currently-populated slots.
    #[cfg(test)]
    pub fn populated_slot_count(&self) -> usize { ... }
}

impl<'a> Drop for ChromLease<'a> {
    fn drop(&mut self) {
        let mut slot = self.fetcher.slots[self.chrom_id as usize].lock().unwrap();
        slot.leases -= 1;
        if slot.leases == 0 {
            slot.bytes = None; // Arc drops, contig freed
        }
    }
}

impl RefSeqFetcher for ChromLeaseSyncRefFetcher {
    fn fetch(&self, chrom_id: u32, start_1based: u32, length: u32) -> Result<Vec<u8>, io::Error> {
        let slot = self.slots[chrom_id as usize].lock().unwrap();
        let bytes = slot.bytes.as_ref().ok_or_else(|| {
            io::Error::new(io::ErrorKind::Other,
                "fetch on unleased chrom — caller must hold a ChromLease")
        })?;
        // bounds-check + slice + return Vec<u8> (existing behaviour;
        // contig is already uppercased at lease time so no per-call
        // uppercase here)
        ...
    }
}
```

**Lazy load body.** Inside `lease`, when `slot.bytes.is_none()`:

1. Open `noodles_fasta::io::indexed_reader::Builder::default().build_from_path(&self.fasta_path)?`.
2. `reader.query(...)` for the contig's full extent (1..=length).
3. Extract the bases as `Vec<u8>`; uppercase once with
   `make_ascii_uppercase` (in place, no per-call allocation).
4. Wrap in `Arc::new(...)` and store in the slot.

We pay the uppercase cost once per chrom-lease lifetime; today's
`SyncRefFetcher` uppercases on every fetch (L4 in the perf review).
This is a side-benefit, not the headline win.

**Contention.** Under the per-chrom outer parallelism each worker
only ever leases its own chrom — distinct mutex slots, zero
contention. Within a worker, `DustFilter` and `PerGroupMerger`
serialise their `fetch` calls naturally (they're pull-iterator
chains), so the single per-slot mutex hold-time per fetch is just
the slice operation.

**Why not pre-warm?** The perf review's L5 recommends pre-warming
into `Vec<Arc<Vec<u8>>>` indexed by `chrom_id`. That fixes
contention (the headline of L5) but is the **opposite** of the
memory goal here — pre-warming is exactly the state we're trying
to avoid. We diverge from L5 on the pre-warm question; we keep the
"drop the noodles `Repository` runtime dep" half. The contention
problem L5 was solving disappears anyway because each worker's
slot is disjoint.

### Wiring Phase B into `var_calling.rs`

[`run_var_calling`](../../src/pop_var_caller/var_calling.rs#L230)
currently:

```rust
let fetcher_concrete = SyncRefFetcher::new(...);
verify_fasta_matches_psp_chromosomes(&fetcher_concrete, &chromosomes)?;  // Phase A removes the fetcher arg
let fetcher: SharedRefFetcher = Arc::new(fetcher_concrete);
// ... per-chrom par_iter ...
```

After Phase B:

```rust
let fetcher_concrete = ChromLeaseSyncRefFetcher::new(args.reference.clone(), contigs_from_parsed(&chromosomes))?;
verify_fasta_matches_psp_chromosomes(&args.reference, &chromosomes)?;  // unchanged from Phase A
let fetcher: SharedRefFetcher = Arc::new(fetcher_concrete);
```

Then in
[`process_one_chromosome`](../../src/pop_var_caller/cohort_driver.rs):

```rust
pub(crate) fn process_one_chromosome(...) -> Result<...> {
    // First thing: lease this worker's chrom. Lives until the
    // function returns.
    let _lease = fetcher.lease(chrom_id).map_err(VarCallingCliError::Io)?;

    // ... existing body (DustFilter / PerGroupMerger / PosteriorEngine
    // / writer) — every fetch they make hits the leased slot ...

    // Implicit drop here releases the lease; slot goes to leases=0,
    // bytes drop, contig memory freed before the next batch of
    // workers picks up.
}
```

The `_lease` binding is the only new line in the worker body.

**One subtlety: `SharedRefFetcher` is `Arc<dyn RefSeqFetcher + Send + Sync>`**
([per_group_merger.rs:485](../../src/var_calling/per_group_merger.rs#L485)).
The lease lives on the concrete type, not the trait object — so
the worker needs access to the concrete `ChromLeaseSyncRefFetcher`,
not just the trait. Two options:

- **(a) Extend the trait with `lease(&self, chrom_id) -> Box<dyn Drop>`.**
  Forces every fetcher impl (including the test mocks) to grow a
  no-op lease. Clean but trait-pollution.
- **(b) Pass the concrete fetcher to `process_one_chromosome` alongside the trait-object alias.**
  The worker holds `&ChromLeaseSyncRefFetcher` for the lease + an
  `Arc<dyn RefSeqFetcher + Send + Sync>` for the pipeline stages.
  No trait change, but the worker now has two references to the
  same object.

Decision: **(b)**. Trait pollution for a behaviour that only one
impl ever does is the wrong axis. The plumbing cost is one extra
argument to `process_one_chromosome` and a cheap
`Arc::downcast`-free clone of the concrete `Arc`. The test mocks
keep working unchanged.

A slimmer **(c)** is possible — drop `SharedRefFetcher` entirely
and have the pipeline stages take a generic `F: RefSeqFetcher`
parameter. The trait already supports that, but it would
monomorphise every stage's code over the fetcher type. Wide ripple,
unclear win. Out of scope.

## Test plan

### Phase A tests

Unit tests in
[src/pop_var_caller/common.rs](../../src/pop_var_caller/common.rs)
(or the new sibling module):

1. **`streaming_md5_matches_oneshot_on_short_contig`** — build a
   small FASTA via `build_fasta`, compute MD5 the new way and the
   old way (`Md5::digest(&fetcher.fetch(...))`), assert equality.
2. **`streaming_md5_matches_oneshot_on_line_wrapped_contig`** — same
   but with a FASTA where the contig is wrapped at 60 chars/line
   (the real-world default).
3. **`streaming_md5_uppercases_soft_masked`** — fixture with a
   mixed-case contig (`ACGTacgtNn`); assert digest matches
   `Md5::digest(b"ACGTACGTNN")`.
4. **`streaming_md5_detects_mismatch`** — wire `verify_fasta_matches_psp_chromosomes`
   end-to-end; pass a contig with the wrong expected MD5; assert
   `FastaVerifyError::Md5Mismatch` with the right contig name.
5. **`streaming_md5_missing_contig_in_fasta`** — `.psp` declares
   contig "chr5" but FASTA only has "chr1"; assert the right
   typed error.
6. **`streaming_md5_truncated_contig`** — fai claims 100 bases but
   the FASTA only has 50; assert `FastaVerifyError::TruncatedContig`.
7. **`verify_runs_in_parallel`** — under a 4-thread rayon pool, pass
   100 small contigs; assert correct (smoke test only — no timing
   assertion).

No new integration tests required: the existing
`tests/cohort_cli_integration.rs` already exercises the verify
step end-to-end on real FASTA/PSP pairs.

### Phase B tests

Unit tests in
[src/per_sample_pileup/ref_fetcher.rs](../../src/per_sample_pileup/ref_fetcher.rs):

1. **`lease_loads_chrom_lazily`** — fresh fetcher has
   `populated_slot_count() == 0`; after `lease(0)` returns,
   `populated_slot_count() == 1`.
2. **`lease_drop_frees_chrom`** — `{ let _l = fetcher.lease(0)?; }`;
   after the scope, `populated_slot_count() == 0`.
3. **`fetch_without_lease_errors`** — calling `fetch(0, 1, 4)`
   without holding a lease returns an error (catches the contract
   violation early).
4. **`fetch_within_lease_returns_bytes`** — golden bytes match the
   existing `SyncRefFetcher` on the same fixture.
5. **`concurrent_leases_same_chrom_share_slot`** — two leases on
   chrom 0; `populated_slot_count() == 1`; drop one; still 1;
   drop the other; 0.
6. **`concurrent_leases_different_chroms_independent`** — leases on
   chrom 0 and chrom 1; `populated_slot_count() == 2`; drop chrom
   0's lease; 1; drop chrom 1's; 0.
7. **`lease_then_lease_then_drop_serial`** — lease chrom 0, drop,
   lease chrom 1, drop; `populated_slot_count` is `1, 0, 1, 0`.
8. **`uppercase_normalises_soft_masked`** — bytes returned by
   `fetch` are uppercase even from a mixed-case FASTA (mirrors the
   existing `SyncRefFetcher` test).
9. **`unknown_chrom_id_in_lease_errors`** — out-of-range chrom_id
   surfaces an `io::Error` from `lease`, not a panic.
10. **`lease_concurrent_threads`** — spawn 4 threads, each leases
    its own chrom_id, holds for a short delay, drops; assert no
    panics + final slot count is 0. Sanity check that the Mutex
    composition is sound under real concurrency.

Integration test in
[tests/cohort_cli_integration.rs](../../tests/cohort_cli_integration.rs):

11. **`var_calling_evicts_chrom_after_worker_finishes`** — run the
    existing multi-chrom fixture; instrument the fetcher with a
    `populated_slot_count` poll on a side channel (or check via a
    `chroms_evicted_total` counter exposed on the fetcher). Assert
    that at no point during the run did `populated_slot_count`
    exceed `min(threads, n_chroms)`. Acceptance: the steady-state
    invariant holds.

### Memory validation (manual, recorded in the impl report)

The integration tests prove the eviction logic is correct. The
actual memory win is measured manually on real data:

- Tomato fixture (`tmp/SRR7279727.multichrom.psp`, 13 chroms,
  ~1.2 GB total reference, T=8):
  - **Before:** peak RSS during `var-calling` (`/usr/bin/time -v`).
  - **After Phase A:** peak RSS expected lower by ~max-chrom-size
    (the verify spike's `Vec<u8>` no longer doubles the in-cache
    bytes).
  - **After Phase B:** peak RSS expected to drop to roughly
    `T × max_chrom_size + working_set` ≈ 8 × 95 MB + overhead ≈
    900 MB, down from full ~1.2 GB cached + working set.

If a human-genome PSP is available at impl time, repeat with T=8 on
24 chroms — expected drop from ~3 GB resident to ~2 GB peak
(T × max chrom ≈ 8 × 250 MB ≈ 2 GB).

The acceptance threshold for Phase B is **steady-state peak RSS
strictly less than `Σ contig.length`**. Anything less than that
would mean the eviction isn't actually firing.

## Sequencing

**Phase A — one PR, three commits:**

1. New `compute_contig_md5_streaming` function + its unit tests.
   Standalone; no callers yet.
2. Reshape `verify_fasta_matches_psp_chromosomes` to take the FASTA
   path and use the new function; remove the `verify_fetcher` line
   from all three callers.
3. Update the doc comment at [common.rs:174-182](../../src/pop_var_caller/common.rs#L174-L182)
   (no more pre-warm claim).

Existing integration tests cover end-to-end. No bench regression
expected (verify is a tiny fraction of total wall on real cohort
runs; the perf review's H5 threshold was "≥ 0.5 s wall reduction").

**Phase B — one PR, four commits:**

1. New `ChromLeaseSyncRefFetcher` type + `ChromLease` RAII +
   unit tests (1-10 above). Standalone; no callers yet.
2. `process_one_chromosome` acquires the lease at the top. Build
   the concrete fetcher in `run_var_calling`, pass it through
   `CohortPipelineParams` as a new field alongside the existing
   `SharedRefFetcher` trait object (option (b) from the design).
3. Drop the construction of the old `SyncRefFetcher` in
   `var_calling.rs`; replace with `ChromLeaseSyncRefFetcher`.
   `estimate_contamination.rs` and `var_calling_from_bam.rs` keep
   their `SyncRefFetcher` (the alias / wrapper described under
   "Out of scope") because they don't yet have per-chrom workers.
4. Integration test (11 above) + manual RSS measurement on the
   tomato fixture; numbers go into the impl report.

The Phase B commit-2 + commit-3 split is mostly for review
ergonomics — they touch the same call sites and could collapse to
one if the reviewer prefers.

## Open work / non-goals

1. **Window-streaming fetch (don't cache the whole contig per
   worker).** Phase B caps peak memory at `T × max_chrom`; that's
   ~2 GB on human at T=8, which is good but not minimal. A future
   slice could push it further by having the fetcher serve fetches
   directly from a windowed reader, never materialising the whole
   contig. Cost: many small fetches per chrom (DUST mask construction
   walks the whole contig once; per-group merger reads many small
   windows). A per-worker rolling window (e.g. 1 MB sliding
   buffer) would absorb most of the small-fetch overhead. Defer
   until `T × max_chrom` shows up as a problem.

2. **DUST mask reuse across same-chrom retries.** Today the
   `DustFilter` builds a fresh mask per chrom via `sdust_mask(&seq, …)`
   inside `ensure_mask_for`. Under Phase B the mask still rebuilds
   on every worker invocation (each worker is a fresh `DustFilter`).
   If the same chrom is processed twice (e.g. a retry after a
   transient I/O error), we recompute. Won't fix; the retry case is
   rare and the mask rebuild is bounded by the contig size we've
   already paid for.

3. **`estimate-contamination` migration.** Same upstream shape; same
   `SyncRefFetcher` use. Will migrate to `ChromLeaseSyncRefFetcher`
   when the per-chrom parallelisation slice for that subcommand
   lands (already on the
   [cohort_per_chromosome_parallel.md follow-ups list](cohort_per_chromosome_parallel.md)).

4. **`var-calling-from-bam` migration.** Same — defer with the
   subcommand's own per-chrom slice. The walker uses
   `ChromBoundaryRefFetcher` already (Stage 1 fetcher with
   per-walker eviction), so memory is already bounded for the
   walker path; only the cohort-side post-walker stages would
   migrate.

5. **`RefSeqFetcher` trait borrowed-return (S2).** Out of scope per
   the perf review's own deferral note. Would be a separate slice
   on top of Phase B.

6. **Pre-warm + parallel MD5 fusion.** A future optimisation could
   fold Phase A's streaming MD5 + Phase B's lease-time load into
   one pass: when a worker leases its chrom, also MD5 it on the
   way in. Today verify and lease are temporally separate
   (verify at startup, lease inside the worker), but they read the
   same bytes from disk. Skipped because (a) it couples startup
   correctness checks to runtime work — fail-fast becomes
   fail-when-the-first-worker-touches-its-chrom; (b) the MD5 of
   chr0 should not be paid by the chr0 worker but by all workers in
   parallel at startup, which is what Phase A already arranges.

## Estimated effort

- **Phase A:** ~150 net lines (new function + tests + callsite
  reshape + doc-comment rewrite). Small PR.
- **Phase B:** ~400 net lines (new fetcher type + lease + tests +
  `process_one_chromosome` change + construction-site swap + impl
  report). Medium PR; the biggest risk surface is the
  `noodles::indexed_reader::query` integration (own File handle per
  load, fai parsing), which is exercised by the existing
  `SyncRefFetcher` so the pattern is known.

Total: ~550 net lines across two PRs. No new runtime dependencies;
removes the runtime use of `noodles_fasta::Repository` (kept only
at lease load time via `noodles_fasta::io::IndexedReader`).

---

## Phase C — sliding-buffer streaming fetcher (added after Phase B benchmarks)

Phase B left a ~480 MB resident-cache contribution per run on the
tomato T=8 measurement (each of the 8 in-flight workers held its
whole contig in memory, ~60 MB average). On a human cohort at T=8
the same shape projects to ~2 GB resident in the fetchers alone.
Since the var-calling pipeline reads its contig **sequentially in
monotonically-non-decreasing position order**, a sliding-window
reader with a buffer of a few MB delivers the same data without
holding the contig.

Memory projection on tomato T=8:
- Phase B (measured): 4.54 GB peak, ~480 MB in fetcher caches.
- Phase C (target): ~4.06 GB peak, ~8 MB in fetcher buffers
  (8 workers × 1 MB sliding buffer).
- On human T=8: drops fetcher contribution from ~2 GB to ~8 MB.

### Access pattern analysis

The two consumers of ref bytes in the cohort var-calling path
(verified by grep over `src/var_calling/` and `src/pop_var_caller/`
on 2026-05-23):

1. **DUST mask construction**
   ([dust_filter.rs:667](../../src/var_calling/dust_filter.rs#L667))
   — currently `fetch(chrom_id, 1, contig.length)`. The whole contig
   in one call. `sdust_mask`
   ([dust_filter.rs:389-444](../../src/var_calling/dust_filter.rs#L389-L444))
   is a forward sequential byte-by-byte scan with rolling state; it
   emits low-complexity intervals as it goes and never indexes
   backwards. **Algorithmically streaming-friendly**; the current
   `&[u8]` parameter is a caller-allocates-everything convenience.
2. **PerGroupMerger**
   ([per_group_merger.rs:704](../../src/var_calling/per_group_merger.rs#L704))
   — many small `fetch(chrom_id, start, span)` calls. Spans are
   bounded by `--var-group-max-span` (default 10 KB; see
   [variant_grouping.rs:31](../../src/var_calling/variant_grouping.rs#L31)).
   The grouper guarantees groups are emitted in non-decreasing
   `start` order, so fetches are monotonic.

The Stage 1 (per-sample pileup walker) BAQ stage also fetches but
uses a different fetcher (`SyncRefFetcher` /
`ChromBoundaryRefFetcher`) and a different parallelism shape;
out of scope for Phase C.

`var_calling_from_bam` runs the pipeline single-sample with no
per-chrom outer parallelism and still uses `SyncRefFetcher`; out of
scope for Phase C. Migration tracked as a follow-up.

### Why both DUST and the per-group merger must change in the same slice

If only the per-group merger streams, the DUST mask still calls
`fetch(1, contig.length)`. The fetcher would have to materialise the
full contig to satisfy that single call. `rayon::par_iter` starts
all workers at roughly the same time, so the *peak* across workers
is back to `T × max_chrom_size` during the DUST phase. The
per-group merger savings come after the DUST peak — too late to
lower peak RSS.

So Phase C ships as one PR with two changes: the new streaming
fetcher and the `sdust_mask` signature refactor.

### Scope

**In scope:**

- New `StreamingChromRefFetcher` in
  [src/per_sample_pileup/ref_fetcher.rs](../../src/per_sample_pileup/ref_fetcher.rs)
  alongside the existing fetchers. Bound to one `chrom_id` like
  `SingleChromRefFetcher`. Internal state: an open `File` + a
  reused 1 MB buffer of uppercased non-newline bases + the buffer's
  base-coordinate offset. `fetch(chrom_id, start, length)`:
  - If `[start, start+length)` lies inside the current buffer →
    slice it out (memcpy-only path, no I/O).
  - Otherwise → seek to `start`, refill buffer with up to
    `max(buffer_size, length)` bases, then slice. Under the
    monotonic-forward access pattern this fires once every ~100
    fetches on average (1 MB / 10 KB = 100).
- Refactor `sdust_mask` from `(&[u8], window, threshold) -> SdustIntervals`
  to `(impl Iterator<Item = u8>, length, window, threshold) -> SdustIntervals`.
  The body stays the same; the `for i in 0..=len { let b = if i < len
  { SEQ_NT4[seq[i] as usize] } else { 4 }; ... }` loop becomes
  `for (i, b) in stream.chain(once_eof_sentinel()).enumerate() { ... }`.
  `length` argument is the contig length needed for the EOF sentinel
  trailing iteration.
- New `bases(&self) -> impl Iterator<Item = io::Result<u8>> + '_`
  method on `StreamingChromRefFetcher` for the DUST mask pass. Pulls
  byte-at-a-time over the streamer's buffer (refilling as needed),
  uppercased + newline-stripped. The mask construction goes through
  this rather than `fetch(1, length)`.
- `DustFilter::ensure_mask_for`
  ([dust_filter.rs:656-683](../../src/var_calling/dust_filter.rs#L656-L683))
  swaps `fetcher.fetch(chrom_id, 1, entry.length)` for the
  streaming iterator. The fetcher trait gains a small extension
  trait (or the streaming method lives on `SingleChromRefFetcher` /
  `StreamingChromRefFetcher` directly without going through the
  trait — DUST already knows the concrete fetcher kind via
  `process_one_chromosome`'s wiring).

**Out of scope (deferred):**

- BAQ / Stage 1 migration. Uses `SyncRefFetcher` /
  `ChromBoundaryRefFetcher` already; not in the cohort parallel
  path.
- `var_calling_from_bam` migration. Single-sample walker; different
  shape.
- A streaming fetcher that supports backward seeks gracefully. Today's
  monotonic-forward access pattern means a backward seek is
  *possible* (degrades to "extra refill per backward jump") but
  not exercised. If a future code path needs random access, the
  current shape still works, just slower.

### Design details

**Streaming over fai-indexed FASTA.** noodles' `IndexedReader::query`
buffers the whole record into a `Vec`, so we can't use it for
chunked reads. We roll our own contig reader:

```rust
struct StreamState {
    file: File,
    /// Cached contig metadata: byte offset, line layout. From the
    /// `.fai` we read at fetcher-construction time.
    fai: FaiEntry,
    /// Buffer of uppercased, newline-stripped bases.
    /// Length ≤ BUFFER_SIZE.
    buf: Vec<u8>,
    /// 1-based base coordinate of the first byte in `buf`.
    /// Invariant: bases in `buf` are contig positions
    /// `[buf_start_base, buf_start_base + buf.len())`.
    buf_start_base: u32,
}
```

`refill(target_base_1based, length_hint)`:

1. Compute file offset of `target_base_1based` via `.fai` math:
   `offset = fai.offset + (b-1) / line_bases * line_width + (b-1) % line_bases`
   where `b = target_base_1based`.
2. `file.seek(SeekFrom::Start(offset))`.
3. `buf.clear()`; read raw FASTA bytes 64 KiB at a time
   (`File::read` into a stack scratch); strip newlines and
   uppercase; push into `buf` until `buf.len() >= length_hint`
   capped at `BUFFER_SIZE`.
4. Set `buf_start_base = target_base_1based`.

`fetch(chrom_id, start, length)`:

1. Validate `chrom_id == self.chrom_id` (same as Phase B's check).
2. Lock `inner: Mutex<StreamState>` (uncontended — single-thread
   use, the Mutex only exists to satisfy the `Sync` bound on the
   trait alias).
3. If `start >= buf_start_base && start + length <= buf_start_base + buf.len()`:
   slice `buf[start - buf_start_base ..]`, return.
4. Otherwise: `refill(start, length)`, then slice + return.

**Why `Mutex`, not `RefCell`.** The trait alias is
`Arc<dyn RefSeqFetcher + Send + Sync>`. Phase B's
`SingleChromRefFetcher` is naturally `Sync` because all its fields
are `Sync` (the inner `Vec<u8>` is read-only after construction).
The streaming fetcher needs interior mutability for the buffer
refill, which makes `Sync` non-automatic. `Mutex` is the cheapest
choice that keeps the trait shape; lock contention is zero in
practice (per Phase B's design, each worker owns its fetcher
outright).

**Buffer size constant.** 1 MB (`1024 * 1024`) — 100× the largest
single fetch (10 KB `--var-group-max-span` default). The buffer
sits in the per-worker `StreamingChromRefFetcher`, so peak across
T workers is `T × 1 MB`. Picking the constant:

- 64 KB: too small; some variant-dense regions would trigger
  refills inside a group fetch.
- 256 KB: 25× the max span. Safe but more refills.
- 1 MB: 100× the max span. Roughly one refill per 100 group
  fetches in dense regions; close to free in sparse regions.
- 4 MB: 400× the max span. Wasteful for the small-genome case
  (most contigs are < 100 MB; 4 MB doesn't help and increases
  the per-worker footprint).

1 MB is the sweet spot. Promote to a constant
`STREAMING_REF_BUFFER_BYTES` with a doc comment explaining the
sizing.

**`sdust_mask` refactor.** Current signature:

```rust
fn sdust_mask(seq: &[u8], window: u32, threshold: u32) -> SdustIntervals
```

New signature:

```rust
fn sdust_mask(
    bases: impl Iterator<Item = u8>,
    length: u32,
    window: u32,
    threshold: u32,
) -> SdustIntervals
```

The body's `for i in 0..=len { let b = if i < len { SEQ_NT4[seq[i] as usize] } else { 4 }; ... }`
becomes:

```rust
let len = length as usize;
let mut bytes = bases.fuse();
for i in 0..=len {
    let b: u32 = match bytes.next() {
        Some(byte) => SEQ_NT4[byte as usize] as u32,
        None => 4, // EOF sentinel — flushes any candidates still in `perf`
    };
    // ... rest of the existing loop unchanged
}
```

The iterator API hides whether bases came from a `Vec<u8>` slice
(today's tests) or from a streaming reader (production). Tests
adapt to pass `seq.iter().copied()` and `seq.len() as u32`.

**DUST mask source.** `DustFilter::ensure_mask_for` previously
read the bytes via `self.fetcher.fetch(chrom_id, 1, entry.length)`
on the trait object. The trait `RefSeqFetcher::fetch` stays
unchanged (returns `Vec<u8>` for compatibility with the
PerGroupMerger fetches). For the DUST mask we need a *streaming*
byte iterator on the concrete fetcher type. Options:

- **(a) Extension trait.** New trait `StreamingRefSeqBytes` with
  `fn bases<'a>(&'a self, chrom_id: u32) -> io::Result<Box<dyn Iterator<Item = io::Result<u8>> + 'a>>`.
  Implemented by `StreamingChromRefFetcher`; the trait-object
  alias becomes `Arc<dyn RefSeqFetcher + StreamingRefSeqBytes + Send + Sync>`.
  Cost: nested-trait combinators are ugly, and `var_calling_from_bam`
  uses `SyncRefFetcher` (which doesn't naturally implement streaming
  bytes) — would have to add a fallback impl that re-reads bytes
  buffered.
- **(b) Pass the concrete fetcher to DUST.** `process_one_chromosome`
  constructs the `StreamingChromRefFetcher` and passes it both to
  `DustFilter::new` (which takes a concrete generic-typed fetcher
  or a streaming iterator directly) and to the
  `CohortPipelineParams.fetcher` slot for the per-group merger
  (via `Arc<dyn RefSeqFetcher>`). DUST's `ensure_mask_for` calls
  `streaming_fetcher.bases()` on the concrete type.

(b) is cleaner. DUST already has a `LoadedChrom` cache that
fires once per chrom transition; it can hold a clone of the
streaming fetcher Arc and call the streaming method. We change
`DustFilter::new` to take `Arc<StreamingChromRefFetcher>` directly
instead of `&dyn RefSeqFetcher`, or we keep `RefSeqFetcher` for
the existing tests and add a new constructor `DustFilter::with_stream`.

Decision: take the concrete `Arc<StreamingChromRefFetcher>` in
`DustFilter::new`. The existing tests in `dust_filter.rs` already
build a `MockFetcher` impl of `RefSeqFetcher`; they migrate to
constructing a minimal `StreamingChromRefFetcher` over a tempfile,
or we add a `from_bytes` constructor on `StreamingChromRefFetcher`
purely for tests (in-memory byte buffer instead of a file). I'll
go with the in-memory test constructor — keeps the tests fast and
hermetic, and the production path stays simple.

### Wiring in `process_one_chromosome`

Replaces today's `SingleChromRefFetcher` construction:

```rust
let fetcher_concrete =
    StreamingChromRefFetcher::new(fasta_path, chrom_id, &chrom_entry.name)?;
let fetcher_arc: Arc<StreamingChromRefFetcher> = Arc::new(fetcher_concrete);
// Trait-object view for the per-group merger fetches (small
// windows; uses RefSeqFetcher::fetch through the trait).
let fetcher: SharedRefFetcher = fetcher_arc.clone();
// Concrete view passed to DUST for the streaming mask construction.
let dust_stream = fetcher_arc;
```

Then `DustFilter::new` accepts `dust_stream` directly.

### Test plan

Unit tests on the new fetcher in
[src/per_sample_pileup/ref_fetcher.rs](../../src/per_sample_pileup/ref_fetcher.rs):

1. **`streaming_fetcher_returns_bytes_for_bound_chrom`** — basic
   slice extraction, single-line FASTA fixture.
2. **`streaming_fetcher_serves_line_wrapped_contig`** — 60-col
   wrapped FASTA, multiple fetches at increasing positions, every
   one returns the correct bases.
3. **`streaming_fetcher_buffer_refill_on_forward_jump`** — fetch
   at `start=1`, then `start=BUFFER_SIZE + 100`; assert second
   fetch returns the right bytes (validates the refill path).
4. **`streaming_fetcher_buffer_refill_on_backward_jump`** —
   degenerate case; documents the slower path. Fetch at
   `start=10000`, then `start=100`; assert correct bytes.
5. **`streaming_fetcher_fetch_past_contig_end_errors`** — same as
   `SingleChromRefFetcher`.
6. **`streaming_fetcher_uppercases_soft_masked`** — same as
   `SingleChromRefFetcher`.
7. **`streaming_fetcher_bases_iterator_yields_all_bases`** —
   driver for the DUST mask pass; iterator length = contig length.
8. **`streaming_fetcher_bases_iterator_uppercases`** — same
   uppercase contract on the iterator path.

Unit tests on the refactored `sdust_mask`:

9. Existing tests adapt to the new signature — pass
   `bytes.iter().copied()` and `bytes.len() as u32`. No expected
   output changes; the algorithm is unchanged.

Integration: existing `tests/cohort_cli_integration.rs` covers
end-to-end. We add no new integration test — output VCFs must be
byte-identical to Phase B's (no algorithmic change).

### Acceptance

- **Peak RSS** drops on the tomato fixture by ~400 MB at T=8.
  Acceptance threshold: ≥ 200 MB reduction vs Phase B; ≥ 600 MB
  reduction vs the baseline. Measured via `tmp/measure_rss.sh`
  back-to-back on `eac3843` (Phase B) vs Phase C HEAD.
- **Output equivalence:** VCFs match Phase B's bit-for-bit (same
  test fixture, same record set, same byte content).
- **No wall regression:** Phase C should not be measurably slower
  than Phase B. Streaming fetchers add a few thousand extra
  `read(2)` syscalls per worker (~contig_size / 64 KB), which on a
  warm page cache is microseconds.

### Estimated effort

- `StreamingChromRefFetcher` + tests: ~200 lines.
- `sdust_mask` signature + DUST caller wiring: ~80 lines.
- `process_one_chromosome` constructor swap: ~30 lines.
- Test fixture / `from_bytes` constructor: ~40 lines.

Total: ~350 lines.
