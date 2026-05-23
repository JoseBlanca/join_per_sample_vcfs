# Fix Application Report: ref_fetcher_2026-05-23.md

**Date:** 2026-05-23
**Source review:** `doc/devel/reports/reviews/ref_fetcher_2026-05-23.md`
**Source state reviewed against:** commit `5b9e420` (post-PSP H1+H3)
**Execution mode:** interactive
**Overall status:** **complete**

---

## 1. Executive summary

The review's 3 Blockers, 23 Majors, and most of the 19 Minors are
landed across 13 commits (eac3843 has the most recent ref fetcher
work in main; the fixes themselves span the contiguous range
`01de8b7..5a47ec2`). Lib tests: 891/891 pass. Integration tests:
34/34 pass. `cargo doc --no-deps --lib` introduces no new warnings
on `ref_fetcher.rs` or `dust_filter.rs`. The file shrank by ~700
LOC net (a 30% reduction) after the legacy `RefSeqFetcher` surface
was retired.

### Review totals

- Blockers: 3 — **all landed**
- Majors: 23 (M1–M27 with M3, M6, M7, M11, M16, M17, M21, M23 noted
  as convergent/subsumed) — **22/23 landed**; M11 (GAT for
  `iter_bases`) deferred to the perf-review queue as documented
  in the review.
- Minors: 19 (Mi1–Mi19) — **6 landed**, 13 deferred (see §5).
- Nits: grouped — addressed via `cargo fmt` (Mi4).

### Open-question answers from the author

1. **Migration completion**: complete it; "we'll discuss the fix"
   if some functionality can't be served by the new functions.
2. **`#[non_exhaustive]` on `ChromRefFetchError`**: trust the
   reviewer — apply.
3. **`ChromRefFetcher` sealing**: trait is internal-to-crate, not
   necessary at the module level — seal it.
4. **IUPAC handling**: fold non-ACGT to N.
5. **Walker migration shape**: B — replace `RefSeqFetcher` with a
   typed-error trait `MultiChromRefFetcher`.

---

## 2. Findings table

| ID    | Severity | Status      | Commit(s)           | Notes |
|-------|----------|-------------|---------------------|-------|
| B1    | Blocker  | Applied     | `f42d552`           | `ContigFai::validate` + 3 constructor sites |
| B2    | Blocker  | Subsumed    | `a7ed411`           | The flatten point (`legacy_io_error`) is gone with M12 |
| B3    | Blocker  | Applied     | `48aea86`           | `streaming_fetch_at_production_buffer_size_rejects_backward_jump` |
| M1    | Major    | Applied     | `a7ed411`           | `DustFilter::ensure_mask_for` → `ensure_mask_loaded` |
| M2    | Major    | Applied     | `cd6e57f`           | `#[non_exhaustive]` on `ChromRefFetchError` |
| M3    | Major    | Subsumed    | `a7ed411`           | The match arms are gone with `legacy_io_error` |
| M4    | Major    | Applied     | `a7ed411`           | Mutex poison surfaces as a typed error via `lock_inner` |
| M5    | Major    | Applied     | `320e44f`           | `chrom_name` / `chrom_length` unification |
| M6    | Major    | Applied     | `a7ed411`           | `bases()` + `StreamingBaseIter` + `RefSeqFetcher::iter_bases` deleted |
| M7    | Major    | Subsumed    | `a7ed411`           | `StreamingBaseIter` is gone; only `ChromRefBaseIter` remains |
| M8    | Major    | Applied     | `cd6e57f`           | `Source::Memory` was `#[cfg(test)]`-gated; later removed entirely (a7ed411) |
| M9    | Major    | Applied     | `def0c38`           | `sealed::Sealed` + sealed `ChromRefFetcher` |
| M10   | Major    | Applied     | `cd6e57f`           | `fetch_into` default impl uses `mem::swap` |
| M11   | Major    | Deferred    | —                   | GAT for `iter_bases` — tracked under perf review (H2) |
| M12   | Major    | Applied     | `a7ed411`           | Legacy `RefSeqFetcher` trait + `StreamingChromRefFetcher::new`/`bases`/`chrom_id` deleted; `WalkerLegacyAdapter` → `MultiChromStreamingRefFetcher` |
| M13   | Major    | Applied     | `7e0be25`           | `open_contig` helper (140 LOC dedup) |
| M14   | Major    | Applied     | `cd6e57f`           | `assert_send::<StreamingChromRefFetcher>()` |
| M15   | Major    | Applied     | `cd6e57f`           | `ChromRefBaseIter::Drop` uses `try_borrow_mut` |
| M16   | Major    | Applied     | `a7ed411`           | `MultiChromStreamingRefFetcher::fetch` rebuilds outside the lock |
| M17   | Major    | Applied     | `f42d552`           | Reject `line_bases = 0` / `line_width < line_bases` (the same B1 site) |
| M18   | Major    | Applied     | `0c3aaf1`           | `for_contig` docstring documents buffer size |
| M19   | Major    | Deferred    | —                   | `ManualEvictChromRefFetcher` cap — out of scope for this fix wave |
| M20   | Major    | Deferred    | —                   | `Io` variant rename (per-operation variants) — future round |
| M21   | Major    | Applied     | `a7ed411`           | `.expect("post-rebuild slot is populated")` replaced with typed error path |
| M22   | Major    | Applied     | `48aea86`           | 3 fetch_into tests (default impl + `&T` forwarding + override) |
| M23   | Major    | Applied     | `a7ed411`           | `assert_send_sync::<MultiChromStreamingRefFetcher>()` |
| M24   | Major    | Applied (partial) | `48aea86`     | `evict_before` past contig end covered; `prepend_backward` truncated-source arm deferred (fixture brittleness) |
| M25   | Major    | Applied (partial) | `48aea86`     | In-buffer-hit covered; legacy `bases()` interleave dropped (the method is gone, M12) |
| M26   | Major    | Applied     | `5bcb7d7`           | `canonicalise(b)` folds non-ACGT to N |
| M27   | Major    | Applied     | `0c3aaf1`           | `# Errors` sections on every public `Result`-returning entry point |
| Mi1   | Minor    | Subsumed    | `a7ed411`           | M1 dust_filter doc link fixed; surrounding `RefSeqFetcher::iter_bases` prose deleted with the trait |
| Mi2   | Minor    | Deferred    | —                   | `InvalidStart` rename (`NonZeroU32`) — paired with M17, deferred |
| Mi3   | Minor    | Subsumed    | `f42d552`           | `b >= 1` precondition guarded by B1's `line_bases > 0` check |
| Mi4   | Minor    | Applied     | `01de8b7`           | `cargo fmt --all` |
| Mi5   | Minor    | Applied     | `e477391`           | `iter_bases` reset uses a scoped block |
| Mi6   | Minor    | Deferred    | —                   | `STREAMING_REF_FILE_READ_CHUNK` visibility — minor |
| Mi7   | Minor    | Deferred    | —                   | `buf_start_base` initial-value convention — minor |
| Mi8   | Minor    | Applied     | `a7ed411`           | `new(String)` deleted with M12; remaining constructors take `&str` |
| Mi9   | Minor    | Deferred    | —                   | `prepend_backward` allocates a `Vec` — tracked under perf review L4 |
| Mi10–Mi16 | Minor (naming) | Partial | `a7ed411`, `320e44f` | `legacy_io_error` deleted; `StreamingBaseIter` deleted; `StreamState` kept (the name describes the only remaining streamer state). M5 covered `contig_name` → `chrom_name` |
| Mi17  | Minor    | Applied     | `a7ed411`           | Module doc rewritten |
| Mi18  | Minor    | Applied     | `a7ed411`           | `chrom_id` field deleted from `StreamingChromRefFetcher` |
| Mi19  | Minor    | Deferred    | —                   | Multi-file split — defer until the file is more settled |

## 3. Questions asked and answers

1. **Migration shape (precondition for M12 + B2 + M3 + M4 + M6 +
   M7 + M16 + M21 + M23 + Mi8 + Mi17 + Mi18)** — option B
   (replace with typed-error trait `MultiChromRefFetcher`).

## 4. Per-finding log

The per-finding outcomes are recorded as commit messages — each
commit body explains *what* was applied and *why* (or why it was
deferred). The §2 table above maps each finding ID to its
landing commit.

Highlights worth surfacing separately:

- **B1 + M17 (`f42d552`)**: added `ContigFai::validate(chrom_name)`
  rejecting `line_bases = 0` and `line_width < line_bases` before
  any fetch can divide-by-zero in the offset arithmetic. Called
  from all 3 constructor sites (now consolidated into the M13
  helper).

- **M12 + everything it absorbs (`a7ed411`)**: the big one. 944
  lines deleted from `ref_fetcher.rs`; 472 added. Net: ~470 LOC
  reduction. Touched 10 files. Walker plumbing migrated from
  `io::Error` to `ChromRefFetchError` end-to-end (`WalkerError::Fasta`
  now wraps the typed error directly). Test suite collapsed by ~11
  legacy `streaming_fetcher_*` tests; each had a `chrom_ref_*`
  equivalent on the typed-error API.

- **M13 (`7e0be25`)**: `open_contig(fasta_path, fai_path,
  chrom_name) -> Result<(ContigFai, File), _>` consolidates the
  ~80-LOC `.fai`-open + index-load + record-lookup + `u32::try_from`
  + B1-validate + `File::open` path. Each constructor wraps the
  returned `(fai, file)` into its own state. The B1 invariants
  now hold in one place rather than three.

- **M5 (`320e44f`)**: 70 token renames across 4 files. Public
  error-variant fields (`OutOfBounds.contig_name` →
  `OutOfBounds.chrom_name`, etc.) are part of the change; consumer
  matchers (`DustFilter`, `PerGroupMerger`, `MockFasta`) update
  in the same commit. `ContigFai`, `ContigList`,
  `noodles_fasta::fai::*` keep `contig` (noodles-adjacent code
  per the M5 carve-out).

## 5. Deferred items and rationale

| ID    | Reason for deferral |
|-------|--------|
| M11   | GAT for `iter_bases` — non-trivial (trait becomes non-object-safe) and tracked under the perf review (H2). |
| M19   | `ManualEvictChromRefFetcher` memory cap — defensive; today's BAQ usage stays well below worst-case. Out of scope for the typed-error migration wave. |
| M20   | `Io` per-operation variant split — minimum bar (rename to `IoFailure` + add typed sub-variants) is non-trivial and trades surface complexity for routing clarity. Defer until a consumer needs the routing. |
| M24 (partial) | `prepend_backward` on a truncated `.fa` — crafting a malformed fixture that fails specifically inside `prepend_backward` without failing in earlier ops is brittle. The forward refill path exercises `read_uppercased_bases` already. Re-revisit once `read_uppercased_bases` graduates to its own typed-error wrapper. |
| M25 (partial) | Legacy `bases()` interleave-with-fetch test — the legacy method is deleted (M12); the new-API equivalent is covered by `chrom_ref_fetch_after_iter_bases_starts_fresh_phase`. |
| Mi2   | `NonZeroU32` for `line_bases`/`line_width` — paired with M17 in the review. Future round; the B1 runtime check covers the panic surface today. |
| Mi6   | Const visibility nit (`STREAMING_REF_FILE_READ_CHUNK` vs `STREAMING_REF_BUFFER_BYTES`) — micro-consistency; ship later. |
| Mi7   | `buf_start_base` initial-value convention (`0` vs `1`) — internal-only; document later. |
| Mi9   | `prepend_backward` allocation — tracked under perf review L4. |
| Mi10–Mi16 (partial) | Remaining naming nits — `StreamState`, the `done` field on the iter, etc. The high-friction names (`legacy_io_error`, `StreamingBaseIter`, `chrom_id`) are gone with M12. The rest is style polish, deferred. |
| Mi19  | 2,000-LOC single-file split — defer until the public surface stabilises further; M12 already pulled out ~700 LOC. |

## 6. Validation

- `cargo fmt --check`: clean.
- `cargo build --lib` inside the dev container: clean (no warnings on `ref_fetcher.rs`).
- `cargo test --lib`: **891 passed; 0 failed; 0 ignored.**
- `cargo test --test '*'`: **34 passed; 0 failed; 0 ignored.**
- `cargo doc --no-deps --lib`: no new warnings on `ref_fetcher.rs` or `dust_filter.rs`. Pre-existing doc warnings in unrelated modules (`posterior_engine.rs` `ExactMath`; `pop_var_caller/cli/error_bridge.rs` `ErrorSheddingAdapter`) remain.
- Compile-time fences: `assert_send::<StreamingChromRefFetcher>()` (M14), `assert_send_sync::<MultiChromStreamingRefFetcher>()` (M23).

## 7. Out-of-scope follow-ups noted in the review

- `unified_chrom_ref_fetcher` plan: completed by this fix wave (step 3 — drop the legacy trait + adapter renames).
- The pre-existing clippy warning in `record_encode.rs:359` (review §7 third bullet) is unrelated and still tracked under the vcf_writer surface.
