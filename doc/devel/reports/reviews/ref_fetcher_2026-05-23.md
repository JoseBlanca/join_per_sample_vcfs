# Code Review: ref_fetcher
**Date:** 2026-05-23
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** `src/per_sample_pileup/ref_fetcher.rs` (the reference-FASTA fetcher module, 2,224 LOC; ~770 LOC tests)
**Status:** **Request-changes** — 3 Blockers, 23 Major, 19 Minor, plus Nits. Production code paths exercised by callers are correct, but the public-trait contracts have material gaps (untrusted-input panics, lossy legacy error wrap, missing `#[non_exhaustive]`, contract-untested invariants) that should be fixed before this module ships as a stable API surface.

---

## 1. Scope

- **What was reviewed:** the single file [src/per_sample_pileup/ref_fetcher.rs](../../../../src/per_sample_pileup/ref_fetcher.rs) in full.
- **Reviewed against:** branch `main` at commit `5b9e420` (post-PSP H1+H3).
- **In-scope files:** [src/per_sample_pileup/ref_fetcher.rs](../../../../src/per_sample_pileup/ref_fetcher.rs).
- **Deliberately out of scope:**
  - Direct callers ([cohort_driver.rs](../../../../src/pop_var_caller/cohort_driver.rs), [var_calling_from_bam.rs](../../../../src/pop_var_caller/var_calling_from_bam.rs), [dust_filter.rs](../../../../src/var_calling/dust_filter.rs), [per_group_merger.rs](../../../../src/var_calling/per_group_merger.rs), [baq/engine.rs](../../../../src/per_sample_pileup/baq/engine.rs), [stage1_pipeline.rs](../../../../src/pop_var_caller/stage1_pipeline.rs), [pileup_to_psp.rs](../../../../src/per_sample_pileup/pileup_to_psp.rs)) — they have their own review surfaces.
  - The legacy `RefSeqFetcher` trait definition in [pileup/mod.rs](../../../../src/per_sample_pileup/pileup/mod.rs) — only its impls in this file are in scope.
- **Categories dispatched (all 10):**
  - `reliability` — always.
  - `errors` — always.
  - `naming` — always.
  - `defaults` — public trait + constructor surface.
  - `idiomatic` — always.
  - `refactor_safety` — always.
  - `unsafe_concurrency` — uses `Mutex`, `RefCell`, `Arc`.
  - `smells` — always.
  - `tooling` — crate has Cargo.toml + the file shows fmt drift + a deny-level rustdoc lint actually firing.
  - `extras` — accepts untrusted input (`--reference` path + contig-name strings + 1-based position arithmetic on attacker-influenced `.fai`); on the hot path; defines a public trait.

## 2. Verdict

**Request-changes.**

Three Blockers must be fixed (one is a reachable panic on attacker-influenced input; one drops contig identity from every legacy diagnostic; one means the load-bearing `OutOfPattern` contract has no production-buffer-size regression test). Many Majors are convergent across 2–6 categories — the stale doc link, the missing `#[non_exhaustive]`, the constructor duplication, and the dual legacy/new surface on `StreamingChromRefFetcher` are all named by multiple sub-agents independently.

## 3. Execution status

| Command | Status | Result |
|---|---|---|
| `cargo fmt --check` | ran | exit 1 — ~50 fmt diffs in this file |
| `cargo clippy --lib --all-features` | ran | exit 0 — **0** ref_fetcher-specific warnings |
| `cargo clippy --lib --all-features -- -D warnings` | ran | exit ≠ 0 — fails ONLY on pre-existing `clippy::let_and_return` in [record_encode.rs:359](../../../../src/var_calling/vcf_writer/record_encode.rs#L359), out of scope |
| `cargo test --lib -- per_sample_pileup::ref_fetcher` | ran | **34/34 pass** |
| `cargo doc --no-deps --lib` | ran | exit ≠ 0 — **the deny gate is actually firing**: [Cargo.toml](../../../../Cargo.toml#L18-L21) sets `[lints.rustdoc] broken_intra_doc_links = "deny"`, but rustdoc emits a warning for the stale link at [ref_fetcher.rs:282](../../../../src/per_sample_pileup/ref_fetcher.rs#L282). |
| `cargo audit` | **NOT RUN** | the `cargo-audit` subcommand is not installed in the dev container |

Findings labeled "Needs verification": **0** — every cited location was read.

## 4. Open questions and assumptions

1. **Migration completion timeline for `WalkerLegacyAdapter` + `legacy_io_error`.** Affects M3, M4, M5, M12, M20, and several smells findings. The module's own prose (lines 8–39) describes the unified-fetcher migration as in-flight, the walker adapter as transitional, and `legacy_io_error` as a flatten-to-`io::Error` shim that goes away in step 3. If step 3 is being scheduled, several Major findings collapse (the duplicate iter, the legacy-trait impl, the lossy error wrap). If it's deferred indefinitely, the Majors against the legacy surface need fixing in place. *Author needs to decide.*

2. **`ChromRefFetchError` evolution policy.** Affects M2 and the smells finding about over-broad `Io` variant. The error is a public type but the crate isn't published. Is the intent that internal callers freely match exhaustively today and accept breakage on variant additions, or is the type expected to evolve (in which case `#[non_exhaustive]` should land now)? See M2.

3. **Public-trait stability bar.** Affects M9 (`ChromRefFetcher` not sealed) and M22 (untested `fetch_into` default impl + `&T` forwarding). If `ChromRefFetcher` is meant to be impl'd outside this module, sealing it is wrong but the test coverage of the trait contract (drop-reset, default `fetch_into`, blanket `&T` impl) needs to be much higher. If it's meant to be internal only, sealing is correct.

4. **Acceptable IUPAC handling.** Affects M26 (stable-output contract violation). The trait doc at lines 562–564 promises SAM-canonical `{A, C, G, T, N}` but the code does `b.to_ascii_uppercase()` which passes `R/Y/S/W/K/M/B/D/H/V` and any other ASCII byte through unchanged. Is the contract correct (fold non-`ACGT` to `N`), or should the docs be tightened to "uppercase-canonicalised SAM bases" (no `N` fold)? Downstream code (DUST sdust mask, BAQ alignment) hits this.

## 5. Top 3 priorities

1. **B1** — `ContigFai::base_to_file_offset` panics on a `.fai` with `line_bases = 0` (which `noodles_fasta::fai::Record` will parse without error). The `--reference` path is attacker-influenced. **Fixed value class:** any malformed `.fai`.
2. **B2** — `legacy_io_error` drops `contig_name` from the `Io` variant when flattening to `io::Error`. Every walker-side diagnostic loses contig identity (`"file ended early"` instead of `"chr5: file ended early"`). One mechanical edit.
3. **B3** — the load-bearing `OutOfPattern` contract has **only** a 4 KiB test-buffer regression. A future refactor that silently swallows the violation (falling through to a refill instead of erroring) would compile and ship without the test catching it.

## 6. Findings

### Blocker

#### B1: [src/per_sample_pileup/ref_fetcher.rs:141-146](../../../../src/per_sample_pileup/ref_fetcher.rs#L141-L146) — `ContigFai::base_to_file_offset` panics on `line_bases = 0` (untrusted-input parser surface)
- **Categories:** extras (Blocker), idiomatic (Minor — same site)
- **Confidence:** High
- **Why it matters:** `noodles_fasta::fai::Record` parses a `.fai` with `line_bases = 0` without error (verified at `~/.cargo/registry/.../noodles-fasta-0.61.0/src/fai/record.rs:184`). The `--reference` path is attacker-influenced. `ContigFai::base_to_file_offset` does `zero_based / self.line_bases as u64` and `zero_based % self.line_bases as u64` — both panic on division-by-zero when `line_bases == 0`. `line_width = 0` and `line_width < line_bases` produce wrong offsets (no panic but the read targets garbage bytes).
- **Fix:** add a typed `ChromRefFetchError::InvalidIndex { contig: String, problem: &'static str }` variant (this is the same edit that fixes M2's `#[non_exhaustive]` issue). Guard at construction in all three constructor sites:
  ```rust
  // for_contig_internal / new / for_contig_with_fai_path, after pulling the fai Record:
  if line_bases == 0 {
      return Err(io::Error::new(io::ErrorKind::InvalidData,
          format!(".fai line_bases is 0 for contig {contig_name}")));
  }
  if line_width == 0 || line_width < line_bases {
      return Err(io::Error::new(io::ErrorKind::InvalidData,
          format!(".fai line_width ({line_width}) is invalid (< line_bases {line_bases})")));
  }
  ```
  (or the equivalent typed variant for the new-API constructors). Add a unit test that builds a corrupted `.fai` and asserts the typed error.

#### B2: [src/per_sample_pileup/ref_fetcher.rs:1049-1065](../../../../src/per_sample_pileup/ref_fetcher.rs#L1049-L1065) — `legacy_io_error` drops `contig_name` from `Io` variant
- **Categories:** errors (Blocker), smells (cross-cat)
- **Confidence:** High
- **Why it matters:** `legacy_io_error` is the single flatten-to-`io::Error` site for every consumer of the legacy `RefSeqFetcher` surface — that's the entire Stage 1 walker + BAQ. The `Io` arm at line 1051 is `ChromRefFetchError::Io { source, .. } => source` — it discards the `contig_name` field and returns only the inner `io::Error`. Every diagnostic from the walker loses the contig that was being fetched. (The `OutOfBounds` / `OutOfPattern` arms use `e.to_string()` for the message, which calls only the parent `Display` — they don't preserve `source()` either, but they at least format the contig name into the message string.)
- **Fix:**
  ```rust
  fn legacy_io_error(e: ChromRefFetchError) -> io::Error {
      match e {
          ChromRefFetchError::Io { source, contig_name } => io::Error::new(
              source.kind(),
              format!("contig {contig_name}: {source}"),
          ),
          // ... other arms unchanged (they already format the contig into the message)
      }
  }
  ```
  Add a regression test: build a fetcher for `chrX`, drive a fetch that triggers an `Io` error (e.g. truncated `.fa`), assert the resulting `io::Error::to_string()` contains `"chrX"`.

#### B3: [src/per_sample_pileup/ref_fetcher.rs (test module)](../../../../src/per_sample_pileup/ref_fetcher.rs) — `OutOfPattern` contract has no production-buffer-size regression test
- **Categories:** reliability (Blocker)
- **Confidence:** High
- **Why it matters:** The streaming-vs-random-access distinction is the load-bearing API contract of this module (documented at lines 540–557). The only test exercising it (around line 1890; verify locally) uses a custom-sized 4 KiB buffer. The production path uses `STREAMING_REF_BUFFER_BYTES = 1 MiB`. A refactor that silently swallows the OutOfPattern violation and falls through to a refill (e.g. by reordering the `start_1based < state.buf_start_base` check) would compile and ship — the tests would still pass on the 4 KiB buffer, and the cohort var-calling path would silently degrade to repeated whole-contig refills instead of erroring.
- **Fix:** add a test
  ```rust
  #[test]
  fn streaming_fetch_at_production_buffer_size_rejects_backward_jump() {
      // Build a synthetic FASTA with one 4 MiB contig so the production
      // buffer size (1 MiB) is exercised, not the test-only 4 KiB.
      let (_dir, path) = build_fasta_with_contig_length("chr0", 4 * 1024 * 1024);
      let f = StreamingChromRefFetcher::for_contig(&path, "chr0").unwrap();
      // Forward fetch, then backward jump past the buffer's origin.
      f.fetch(2 * 1024 * 1024, 16).unwrap();
      let err = f.fetch(1, 16).unwrap_err();
      assert!(matches!(err, ChromRefFetchError::OutOfPattern { .. }));
  }
  ```
  Reliability sub-agent's full challenge-test list is in [§8](#8-missing-tests-to-add-now).

### Major

#### M1: [src/per_sample_pileup/ref_fetcher.rs:282](../../../../src/per_sample_pileup/ref_fetcher.rs#L282) — Stale intra-doc link `DustFilter::ensure_mask_for` (renamed to `ensure_mask_loaded`)
- **Categories:** tooling, refactor_safety, idiomatic, reliability (cross-cat), errors (Nit), extras (cross-cat) — **6 categories converge**
- **Confidence:** High
- **Why it matters:** [Cargo.toml:18-21](../../../../Cargo.toml#L18-L21) sets `[lints.rustdoc] broken_intra_doc_links = "deny"`. `cargo doc --no-deps --lib` emits a warning for this link, so the deny gate either isn't run in CI or is being bypassed. The target was renamed in commit `8ad10f6` (the cohort `ChromRefFetcher` migration); the current method is at [dust_filter.rs:704](../../../../src/var_calling/dust_filter.rs#L704).
- **Fix:**
  ```diff
  -    /// construction in
  -    /// [`crate::var_calling::dust_filter::DustFilter::ensure_mask_for`].
  +    /// construction in
  +    /// [`crate::var_calling::dust_filter::DustFilter::ensure_mask_loaded`].
  ```
  Then audit why the deny gate isn't firing in CI — `cargo doc --no-deps -- -D warnings` should be in the precommit / CI script.

#### M2: [src/per_sample_pileup/ref_fetcher.rs:499](../../../../src/per_sample_pileup/ref_fetcher.rs#L499) — `ChromRefFetchError` missing `#[non_exhaustive]`
- **Categories:** errors, refactor_safety
- **Confidence:** High
- **Why it matters:** The module-doc describes this error as evolving; the in-scope-now `legacy_io_error` matches it exhaustively. Without `#[non_exhaustive]`, adding any variant (e.g. `InvalidIndex` from B1, or a new `Truncated` from a future bytes-split tightening) is a breaking change for every external matcher and a silent miss for any future wildcard match elsewhere in the tree.
- **Fix:**
  ```diff
   /// Why a [`ChromRefFetcher`] call failed.
  +#[non_exhaustive]
   #[derive(Debug, thiserror::Error)]
   pub enum ChromRefFetchError {
  ```

#### M3: [src/per_sample_pileup/ref_fetcher.rs:1050-1064](../../../../src/per_sample_pileup/ref_fetcher.rs#L1050-L1064) — `legacy_io_error` arms use `..` rest pattern, hiding future fields
- **Categories:** refactor_safety
- **Confidence:** High
- **Why it matters:** The `OutOfBounds { .. }` and `OutOfPattern { .. }` arms absorb every field. A future field on either variant compiles silently here even though the new field might need to change the `io::ErrorKind` mapping. This is the same single flatten point as B2 — silent miss routes errors wrong.
- **Fix:** name every field explicitly with `_`:
  ```rust
  ChromRefFetchError::OutOfBounds { contig_name: _, contig_length: _, start: _, end: _ } => ...
  ChromRefFetchError::OutOfPattern { requested_start: _, buffer_origin: _ } => ...
  ```
  (Lands together with B1's new variant addition; the compiler will then force a touch.)

#### M4: [src/per_sample_pileup/ref_fetcher.rs:1001](../../../../src/per_sample_pileup/ref_fetcher.rs#L1001) — `WalkerLegacyAdapter::fetch` silently swallows `Mutex` poison
- **Categories:** errors, unsafe_concurrency (cross-cat)
- **Confidence:** High
- **Why it matters:** `self.inner.lock().unwrap_or_else(|e| e.into_inner())` accepts a poisoned mutex without surfacing the fact. The walker can then continue with potentially-inconsistent inner state (post-panic mid-rebuild, the inner is `None` and the next fetch will rebuild it; in other poison scenarios the state is unknown). The project memory `feedback_no_logs_use_errors` says to surface invariant violations as typed errors, not logs.
- **Fix:** either propagate poison as an `io::Error` (preferred for a transitional adapter), or document the safe-by-luck argument inline:
  ```rust
  let mut slot = match self.inner.lock() {
      Ok(g) => g,
      Err(_poison) => {
          // The walker is single-threaded inside a chrom worker;
          // the only way to poison this mutex is for the inner
          // StreamingChromRefFetcher's rebuild path to panic, which
          // panic-frees the rebuild critical section. The post-poison
          // state is the inner = None invariant from PoisonError::into_inner;
          // we restart the rebuild on the next fetch.
          return Err(io::Error::new(
              io::ErrorKind::Other,
              "WalkerLegacyAdapter mutex poisoned (prior panic during rebuild)",
          ));
      }
  };
  ```

#### M5: Module-wide — `chrom` vs `contig` naming inconsistency for the same concept
- **Categories:** naming
- **Confidence:** High
- **Why it matters:** The file uses `chrom_id: u32` (field, errors), `contig_name: String` (field, errors), `for_contig` (constructor on the new trait), `OutOfBounds { contig_name, contig_length }` (error fields). The two words refer to the same thing inside this codebase but are mixed within a single error variant ([ref_fetcher.rs:506-511](../../../../src/per_sample_pileup/ref_fetcher.rs#L506-L511)). Project memory `reference_genomics_abbreviations` settles on `chrom`; the noodles `fai` crate uses `contig`. The mix raises the reader's burden every time.
- **Fix:** unify on `chrom` for the project's surface (matches `ChromRefFetcher` trait, `ChromRefBaseIter`, the perf reviews); keep `contig` only for noodles-adjacent code (`ContigFai`, the noodles `.fai` parser). Concretely:
  ```diff
  - OutOfBounds { contig_name: String, contig_length: u32, ... }
  + OutOfBounds { chrom_name: String, chrom_length: u32, ... }
  ```
  And rename `contig_name: String` field on `StreamingChromRefFetcher` to `chrom_name`. Touches every error-construction site (mechanical; one PR).

#### M6: [src/per_sample_pileup/ref_fetcher.rs:283-289, 293-319, 651-657, ...](../../../../src/per_sample_pileup/ref_fetcher.rs#L283-L289) — `StreamingChromRefFetcher::bases()` iter is exposed THREE ways
- **Categories:** naming
- **Confidence:** High
- **Why it matters:** The same per-base iterator is reachable as `StreamingChromRefFetcher::bases()` (inherent, returns `StreamingBaseIter`), as `RefSeqFetcher::iter_bases` (returns the same), and as `ChromRefFetcher::iter_bases` (returns the new `ChromRefBaseIter`). Two of these have no production caller (see M12). The third is the right one.
- **Fix:** depends on M12 — remove the legacy `bases()` + `StreamingBaseIter` + `RefSeqFetcher::iter_bases` impl on this type.

#### M7: [src/per_sample_pileup/ref_fetcher.rs:376, 782](../../../../src/per_sample_pileup/ref_fetcher.rs#L376) — `StreamingBaseIter` vs `ChromRefBaseIter` names don't communicate the legacy/new distinction
- **Categories:** naming
- **Confidence:** Medium
- **Why it matters:** Both names are passable in isolation. Together they describe iteration-style differences ("streaming" vs "chrom-ref") when the actual distinguisher is "this serves the legacy `RefSeqFetcher` trait" vs "this serves the new `ChromRefFetcher` trait". The convention `WalkerLegacyAdapter` already uses for this is in the codebase.
- **Fix:** subsumed by M12 (delete `StreamingBaseIter`). If kept, rename to `LegacyChromRefBaseIter`.

#### M8: [src/per_sample_pileup/ref_fetcher.rs:149-156](../../../../src/per_sample_pileup/ref_fetcher.rs#L149-L156) — `Source::Memory` variant lives in production binaries
- **Categories:** idiomatic (Major), smells (Minor), refactor_safety (verified)
- **Confidence:** Medium
- **Why it matters:** The variant is gated with `#[cfg_attr(not(test), allow(dead_code))]` rather than `#[cfg(test)]`. The discriminant is present in production builds, the match arms at lines 159–176 expand to two-variant dispatches in release code, and `Source` carries an extra ~32 bytes per `StreamState` for a variant that can never be constructed outside tests. Refactor_safety confirms the match exhaustiveness is preserved either way.
- **Fix:** `#[cfg(test)] Memory(io::Cursor<Vec<u8>>),` plus matching `#[cfg(test)]` arms in `Source::seek_to` and `Source::read_chunk`. The compiler will then drop the variant from production builds.

#### M9: [src/per_sample_pileup/ref_fetcher.rs:558](../../../../src/per_sample_pileup/ref_fetcher.rs#L558) — `ChromRefFetcher` trait is not sealed
- **Categories:** idiomatic
- **Confidence:** Medium
- **Why it matters:** The trait has a non-trivial cross-call contract (the `iter_bases` Drop-reset, the monotonic-forward `fetch` invariant) that only `StreamingChromRefFetcher` honours and only one consumer (`DustFilter`) relies on. A downstream impl could silently violate the contract; adding a new method (e.g. `fetch_into` was added in commit `503b7cc`) is a breaking change. Open question 3 asks the author whether the trait is intended for external impls.
- **Fix:** if the trait is internal only, seal it:
  ```rust
  mod sealed { pub trait Sealed {} }
  pub trait ChromRefFetcher: sealed::Sealed { ... }
  impl sealed::Sealed for StreamingChromRefFetcher {}
  impl<T: ChromRefFetcher + ?Sized> sealed::Sealed for &T {}
  ```

#### M10: [src/per_sample_pileup/ref_fetcher.rs:595-612](../../../../src/per_sample_pileup/ref_fetcher.rs#L595-L612) — `fetch_into` default impl allocates instead of swapping
- **Categories:** idiomatic
- **Confidence:** High
- **Why it matters:** The docstring promises "allocation-free in the hot path when `dst` is reused"; the default impl is `let v = self.fetch(start, length)?; dst.clear(); dst.extend_from_slice(&v); Ok(())` — that's exactly the alloc-then-extend pattern the method exists to avoid. For non-`StreamingChromRefFetcher` impls (or any future overrideless impl), the contract silently doesn't hold.
- **Fix:**
  ```diff
   fn fetch_into(&self, start_1based: u32, length: u32, dst: &mut Vec<u8>)
       -> Result<(), ChromRefFetchError>
   {
  -    let v = self.fetch(start_1based, length)?;
  -    dst.clear();
  -    dst.extend_from_slice(&v);
  -    Ok(())
  +    let mut v = self.fetch(start_1based, length)?;
  +    std::mem::swap(dst, &mut v);
  +    Ok(())
   }
  ```
  Or document the divergence and add `#[must_override]`-style comment.

#### M11: [src/per_sample_pileup/ref_fetcher.rs:582-587](../../../../src/per_sample_pileup/ref_fetcher.rs#L582-L587) — `iter_bases` returns `Box<dyn Iterator>` on the DUST hot path
- **Categories:** idiomatic, extras (cross-cat — performance guard)
- **Confidence:** Medium
- **Why it matters:** Already raised in the 2026-05-23 perf review as H2 (deferred). Re-flagged here because it's the canonical example of "static dispatch on a generic surface" — `DustFilter` is generic over `F: ChromRefFetcher`, but the trait method itself re-boxes the inner iter, defeating monomorphisation.
- **Fix:** GAT — `type BaseIter<'a>: Iterator<Item = Result<u8, ChromRefFetchError>> where Self: 'a;`. Trade-off: trait becomes non-object-safe. Open question 3 affects whether this is acceptable. Already tracked under the perf review's open items.

#### M12: [src/per_sample_pileup/ref_fetcher.rs:194-237, 283-289, 292-374, 376-414](../../../../src/per_sample_pileup/ref_fetcher.rs#L194-L237) — Vestigial legacy surface on `StreamingChromRefFetcher`
- **Categories:** smells (Major), idiomatic (Minor)
- **Confidence:** High (`grep -rn StreamingBaseIter\\|RefSeqFetcher.*StreamingChromRefFetcher` across `src/`)
- **Why it matters:** Three pieces on this type have no production caller after the cohort migration:
  - `StreamingChromRefFetcher::new` (line 194): legacy multi-chrom constructor with `chrom_id` parameter. Replaced by `for_contig` (line 656).
  - `bases()` + `StreamingBaseIter` (lines 283, 376): legacy bases iter. Replaced by `iter_bases() -> ChromRefBaseIter` (line 870-ish).
  - `RefSeqFetcher` impl on `StreamingChromRefFetcher` (line 292): the walker now goes through `WalkerLegacyAdapter`, which delegates to `RefSeqFetcher::fetch` on its **inner** `StreamingChromRefFetcher`. The inner's `iter_bases` impl is still reachable, but `WalkerLegacyAdapter` doesn't expose it — the inner is fetch-only.
- **Fix:** verify with `grep` (the smells sub-agent did), then delete. If `WalkerLegacyAdapter` still calls `RefSeqFetcher::fetch` directly on the inner, keep only that one method; drop the `iter_bases` impl. Save ~80–100 LOC and remove M6's and M7's underlying duplication.

#### M13: [src/per_sample_pileup/ref_fetcher.rs:690-775, 1131-1211](../../../../src/per_sample_pileup/ref_fetcher.rs#L690-L775) — ~80 lines of `.fai`-open + index-load + error-wrap duplicated across three constructors
- **Categories:** smells (Major), idiomatic (Minor), extras (Minor), reliability (cross-cat)
- **Confidence:** High
- **Why it matters:** `StreamingChromRefFetcher::for_contig_internal` and `ManualEvictChromRefFetcher::for_contig_with_fai_path` repeat the same sequence: derive `<fasta>.fai` path → `fai::fs::read` → linear-search by contig name → `u32::try_from` overflow checks on `length`/`line_bases`/`line_width` → open `<fasta>`. A third partial copy lives in the legacy `new`. The B1 fix needs to land in all three sites unless they're unified.
- **Fix:** extract
  ```rust
  fn open_contig(
      fasta_path: &Path,
      contig_name: &str,
  ) -> io::Result<(ContigFai, File)> {
      // The ~80 LOC body, returning (ContigFai, File). Validates
      // line_bases > 0 + line_width >= line_bases (the B1 fix).
  }
  ```
  Then each constructor wraps `(fai, file)` into its own state struct. Tested once, B1's invariants enforced once.

#### M14: [src/per_sample_pileup/ref_fetcher.rs:104-121](../../../../src/per_sample_pileup/ref_fetcher.rs#L104-L121) — `StreamingChromRefFetcher`'s `Send + !Sync` is held by accident of field choice
- **Categories:** unsafe_concurrency
- **Confidence:** Medium
- **Why it matters:** The cohort concurrency model documented at lines 113–119 depends on the type being `Send` (the `Arc<dyn ChromRefFetcher + Send>` alias) and `!Sync` (interior mutability via `RefCell`). Today this holds because `RefCell<StreamState>` auto-derives both bounds correctly. A future field change (e.g. adding an `Arc<Whatever>` or wrapping in `Mutex` — both ironic) could silently flip `Sync` and let the cohort driver share the fetcher across threads via accident, defeating the H1 refactor.
- **Fix:** add an explicit marker + compile-time fence:
  ```rust
  pub struct StreamingChromRefFetcher {
      // ... fields ...
      inner: RefCell<StreamState>,
      _not_sync: PhantomData<Cell<()>>,  // explicit: !Sync
  }
  // ... in a const block near the type ...
  const _: () = {
      const fn assert_send<T: Send>() {}
      assert_send::<StreamingChromRefFetcher>();
      // No Sync assertion — !Sync is the contract.
  };
  ```

#### M15: [src/per_sample_pileup/ref_fetcher.rs:828-840](../../../../src/per_sample_pileup/ref_fetcher.rs#L828-L840) — `ChromRefBaseIter::drop` can double-panic
- **Categories:** unsafe_concurrency
- **Confidence:** Medium
- **Why it matters:** `Drop::drop` does `let mut state = self.fetcher.inner.borrow_mut();`. If `next()` panicked while holding its own `borrow_mut` (which can happen on a future bug — out-of-bounds buffer index after a refill regression), `Drop` re-enters and the `borrow_mut` panics again → double-panic → abort.
- **Fix:**
  ```diff
   impl Drop for ChromRefBaseIter<'_> {
       fn drop(&mut self) {
  -        let mut state = self.fetcher.inner.borrow_mut();
  -        state.buf.clear();
  -        state.buf_start_base = 0;
  +        if let Ok(mut state) = self.fetcher.inner.try_borrow_mut() {
  +            state.buf.clear();
  +            state.buf_start_base = 0;
  +        }
  +        // If the borrow fails, an outer next() panicked mid-borrow;
  +        // skip the reset rather than double-panic. The next iter_bases
  +        // call resets the buffer anyway (see line ~870).
       }
   }
  ```

#### M16: [src/per_sample_pileup/ref_fetcher.rs:1001](../../../../src/per_sample_pileup/ref_fetcher.rs#L1001) — `WalkerLegacyAdapter::fetch` holds its `Mutex` across the inner `for_contig` rebuild (file open + `.fai` parse + I/O)
- **Categories:** unsafe_concurrency
- **Confidence:** Medium
- **Why it matters:** The lock is acquired at line 1001 and held across the full chrom-swap (file open, `.fai` parse, the actual delegated fetch). When multiple Stage 1 reader threads share this adapter (it's the only `+ Send + Sync` fetcher), all of them serialise through this lock — the rebuild can take milliseconds (`.fai` parse is ~26 KB of text). The lock-holding window also makes the `Mutex` poison swallow (M4) more dangerous: a panic mid-rebuild poisons the lock and leaves `inner = None`, which is the safe-by-luck recovery state — but only because of M4's documented invariant.
- **Fix:** build outside the lock, then swap in:
  ```rust
  // pseudocode
  let needs_rebuild = {
      let slot = self.inner.lock().unwrap_or_else(...);
      match &*slot { Some((id, _)) if *id == chrom_id => false, _ => true }
  };
  if needs_rebuild {
      let new_inner = StreamingChromRefFetcher::for_contig(...)?;  // outside lock
      let mut slot = self.inner.lock().unwrap_or_else(...);
      *slot = Some((chrom_id, new_inner));
  }
  let slot = self.inner.lock()...;
  // delegated fetch
  ```
  This shortens the critical section from "rebuild + fetch" to "swap + fetch".

#### M17: [src/per_sample_pileup/ref_fetcher.rs:141-146](../../../../src/per_sample_pileup/ref_fetcher.rs#L141-L146) — `ContigFai::base_to_file_offset` uses `as` casts for u32→u64
- **Categories:** idiomatic
- **Confidence:** Medium
- **Why it matters:** `(b - 1) as u64` and the `line_bases as u64` casts are infallible u32→u64 widening so they don't lose data, but the `b - 1` subtraction wraps on `b == 0`. There's no precondition check. B1 already covers the panic surface (`line_bases = 0`); this is the same site, requesting tighter typing.
- **Fix:** introduce `NonZeroU32` on `line_bases`/`line_width`:
  ```rust
  struct ContigFai {
      seq_offset: u64,
      length: u32,
      line_bases: NonZeroU32,
      line_width: NonZeroU32,
  }
  fn base_to_file_offset(&self, b: NonZeroU32) -> u64 { ... }
  ```
  The `NonZeroU32` parse falls into the constructor where B1's check already lives.

#### M18: [src/per_sample_pileup/ref_fetcher.rs:656](../../../../src/per_sample_pileup/ref_fetcher.rs#L656) — `StreamingChromRefFetcher::for_contig` hard-codes `STREAMING_REF_BUFFER_BYTES`; not surfaced
- **Categories:** defaults
- **Confidence:** High
- **Why it matters:** Callers see `StreamingChromRefFetcher::for_contig(&path, "chr0")` and have no way to know the resident memory footprint without reading the source. A larger contig (human chr1) at the production buffer size is fine; an experimental caller picking the wrong access pattern (random-access at a tighter budget) needs to either compose a different fetcher or accept the silent 1 MB cost.
- **Fix:** either:
  (a) document the buffer size on the constructor's docstring (one line — minimum bar);
  (b) split into `for_contig` (default 1 MB) and `for_contig_with_buffer_size(path, contig, buf_bytes)` (explicit override). Production callers stay on `for_contig`.

#### M19: [src/per_sample_pileup/ref_fetcher.rs:1122](../../../../src/per_sample_pileup/ref_fetcher.rs#L1122) — `ManualEvictChromRefFetcher` has no memory cap
- **Categories:** defaults
- **Confidence:** High
- **Why it matters:** The fetcher's whole point is "caller-managed eviction" (lines 1067-1106), and the worst case if the caller forgets to call `evict_before` is one whole contig resident — ~250 MB on human chr1. The constructor doesn't take a cap or even warn. A reviewer reading just the constructor signature has no idea this could happen.
- **Fix:** add an explicit `with_max_buffer_bytes(max: usize)` builder method (or a constructor parameter) that errors on `fetch` calls whose resident range would exceed the cap. Default cap = `STREAMING_REF_BUFFER_BYTES * 4` (= 4 MB) — chosen so BAQ's typical access patterns work but a pathological one trips. Documented at the constructor.

#### M20: [src/per_sample_pileup/ref_fetcher.rs:530-537](../../../../src/per_sample_pileup/ref_fetcher.rs#L530-L537) — `ChromRefFetchError::Io` is a mechanism-named over-broad catch-all
- **Categories:** errors
- **Confidence:** Medium
- **Why it matters:** The `Io` variant currently covers `.fai` parse, contig-not-found, `try_from` overflow, file open, seek, refill, and arithmetic-overflow input validation. The variant name describes the mechanism (`io::Error`), not the failed operation. The errors checklist rule: "variant names describe the failed operation, not the mechanism."
- **Fix:** split into per-operation variants over a future migration round. Minimum: rename to `IoFailure` to make the catch-all nature explicit, and add typed variants (`FaiParseError`, `ContigNotFound`, `InvalidIndex` — the B1 fix is one of them). Each carries `contig_name` so the legacy flatten in B2 doesn't have to special-case.

#### M21: [src/per_sample_pileup/ref_fetcher.rs:1024](../../../../src/per_sample_pileup/ref_fetcher.rs#L1024) — `.expect("post-rebuild slot is populated")` without `// PANIC-FREE:` annotation
- **Categories:** errors
- **Confidence:** Medium
- **Why it matters:** Project convention (see existing `// PANIC-FREE:` comments at e.g. [reader.rs:1201](../../../../src/per_sample_pileup/psp/reader.rs#L1201)) is that every surviving `.expect()` is annotated with a load-bearing invariant comment. This one isn't. The invariant is correct (the rebuild branch sets `Some(...)` immediately above), but the reader has to reconstruct that. Spread across the codebase, this is a silent-drift risk.
- **Fix:** prepend the standard annotation:
  ```rust
  // PANIC-FREE: the immediately-prior `if needs_rebuild { *slot = Some(...) }`
  // branch sets the slot when this code runs; the `match slot { None => ..., Some(...) => }`
  // would have returned earlier otherwise.
  let (_, inner) = slot.as_ref().expect("post-rebuild slot is populated");
  ```
  Or restructure to avoid the `.expect()` entirely.

#### M22: [src/per_sample_pileup/ref_fetcher.rs:600, 619, 856](../../../../src/per_sample_pileup/ref_fetcher.rs#L600) — `ChromRefFetcher::fetch_into` default impl + `&T` forwarding + `dst.clear()` invariant untested
- **Categories:** reliability
- **Confidence:** High
- **Why it matters:** The default `fetch_into` is the contract that downstream impls inherit if they don't override. The `&T` blanket impl ([line 619](../../../../src/per_sample_pileup/ref_fetcher.rs#L619)) is the seam the cohort driver crosses with `&*fetcher`. The `dst.clear()` behaviour is documented at line 595–612 ("`dst` is cleared at the start") — a regression that prepended stale bytes would silently corrupt the per-group ref window.
- **Fix:** add three tests:
  ```rust
  #[test] fn fetch_into_default_impl_clears_dst() { /* feed a stub that only overrides `fetch`, prime dst with "GARBAGE", call fetch_into, assert dst.starts_with(canonical_bytes) and len == length */ }
  #[test] fn fetch_into_through_ref_forwarding_works() { /* same shape, called via `&fetcher.fetch_into(...)` to exercise the blanket impl */ }
  #[test] fn streaming_fetch_into_returns_exact_length_from_slab() { /* the optimized override; assert dst.len() == length after a fetch_into spanning a refill */ }
  ```

#### M23: [src/per_sample_pileup/ref_fetcher.rs:976](../../../../src/per_sample_pileup/ref_fetcher.rs#L976) — `WalkerLegacyAdapter` `Send + Sync` claim untested at compile time, no concurrent-fetch regression
- **Categories:** reliability
- **Confidence:** High
- **Why it matters:** The adapter is the only `Send + Sync` fetcher in the module. Both M14 (auto-derived bounds held by accident) and M16 (lock-holding window) apply to it. A regression that flipped `Sync` (e.g. by adding a non-`Sync` field) would silently break the Stage 1 walker; today only the integration tests would catch it, indirectly.
- **Fix:**
  ```rust
  const _: () = {
      const fn assert_send_sync<T: Send + Sync>() {}
      assert_send_sync::<WalkerLegacyAdapter>();
  };

  #[test]
  fn walker_adapter_serves_two_threads_alternating_chroms() {
      // Build two-chrom FASTA, share Arc<adapter> across two threads,
      // each fetching alternately from chr0 and chr1. Run for N
      // iterations; assert no deadlock and bytes match the per-chrom
      // expected payload.
  }
  ```

#### M24: [src/per_sample_pileup/ref_fetcher.rs:~1284, ~1366](../../../../src/per_sample_pileup/ref_fetcher.rs#L1284) — `ManualEvictChromRefFetcher::evict_before` past contig end + `prepend_backward` truncated-source untested
- **Categories:** reliability
- **Confidence:** Medium (line numbers approximate — verify on author's read)
- **Why it matters:** Both edge cases are documented in code comments but lack tests. `evict_before(end_pos+1)` is a legitimate caller intent (clear the buffer); a refactor that mis-handles it would corrupt subsequent fetches. `prepend_backward` on a `.fa` truncated mid-line is a corruption case that should surface a typed error.
- **Fix:** add the two regression tests.

#### M25: [src/per_sample_pileup/ref_fetcher.rs:~1222, 283](../../../../src/per_sample_pileup/ref_fetcher.rs#L1222) — `ManualEvictChromRefFetcher::fetch` in-buffer hit + `StreamingChromRefFetcher::bases()` interleave with `fetch` untested
- **Categories:** reliability
- **Confidence:** Medium
- **Why it matters:** The in-buffer-hit path (line ~1222) is the warm-cache happy path for BAQ; the interleave-then-fetch path (line 283 calling pattern) is the legacy contract used by the walker. Both are documented but untested.
- **Fix:** add the two tests.

#### M26: [src/per_sample_pileup/ref_fetcher.rs:562-564](../../../../src/per_sample_pileup/ref_fetcher.rs#L562-L564) — Stable-output contract violated: `to_ascii_uppercase` passes IUPAC ambiguity codes through unchanged
- **Categories:** extras
- **Confidence:** Medium
- **Why it matters:** The trait docs promise "SAM-spec canonical (`A`/`C`/`G`/`T`/`N`) regardless of how the FASTA is masked on disk". `refill` and `read_uppercased_bases` call `b.to_ascii_uppercase()`, which only folds case — `R/Y/S/W/K/M/B/D/H/V` (IUPAC ambiguity), `U` (RNA), `-` (gap), and any other ASCII byte pass through. Downstream code (DUST's sdust mask, BAQ alignment, the per-group merger's ref slice) all see the raw byte. Whether this is a bug or a docstring overclaim is Open question 4.
- **Fix:** if the docs are correct, fold non-`ACGT` to `N`:
  ```rust
  fn canonicalise(b: u8) -> u8 {
      match b.to_ascii_uppercase() {
          b'A' | b'C' | b'G' | b'T' => b.to_ascii_uppercase(),
          _ => b'N',
      }
  }
  ```
  applied per byte in `refill`'s output loop. If the docs are wrong (the contract is "uppercase, no folding"), retire the `ACGTN` promise and document the actual semantics.

#### M27: [src/per_sample_pileup/ref_fetcher.rs:282](../../../../src/per_sample_pileup/ref_fetcher.rs#L282) — and four other public `Result`-returning entry points lack `# Errors` sections
- **Categories:** extras
- **Confidence:** High
- **Why it matters:** `StreamingChromRefFetcher::for_contig`, `::new`, `ManualEvictChromRefFetcher::for_contig`, the trait's `fetch` / `iter_bases` / `fetch_into` — none document `# Errors`. Rust API guidelines (C-FAILURE) require it on every `pub fn ... -> Result<_, _>`.
- **Fix:** add `# Errors` sections enumerating the variants each constructor / method can return. Mechanical.

### Minor

#### Mi1: [src/per_sample_pileup/ref_fetcher.rs:282](../../../../src/per_sample_pileup/ref_fetcher.rs#L282) and surrounding docstring rotation
- **Categories:** errors (Nit), reliability (cross-cat)
- Subsumed by M1; listed separately because the broader docstring around the stale link also references `RefSeqFetcher::iter_bases` consumers that no longer exist after M12.

#### Mi2: [src/per_sample_pileup/ref_fetcher.rs:512-514](../../../../src/per_sample_pileup/ref_fetcher.rs#L512-L514) — `InvalidStart` name/message describes the rule rather than the operation
- **Categories:** errors
- **Fix:** rename to `StartCoordZero` or `start: u32` becomes `NonZeroU32` and the variant goes away entirely (paired with M17).

#### Mi3: [src/per_sample_pileup/ref_fetcher.rs:139](../../../../src/per_sample_pileup/ref_fetcher.rs#L139) — `ContigFai::base_to_file_offset` precondition `b >= 1` is undocumented and unchecked
- **Categories:** errors
- **Fix:** subsumed by M17 (`NonZeroU32` lifts to the type system).

#### Mi4: Module-wide — ~50 `cargo fmt` diffs in this file
- **Categories:** tooling (Minor), naming (Nit), idiomatic (Nit)
- **Fix:** run `cargo fmt`. The file is not auto-formatted (imports unsorted, some `format!` lines exceed 100 cols).

#### Mi5: [src/per_sample_pileup/ref_fetcher.rs:1018](../../../../src/per_sample_pileup/ref_fetcher.rs#L1018) — `drop(state)` in `iter_bases` is load-bearing but uncommented
- **Categories:** unsafe_concurrency
- **Fix:** either wrap the borrow in a scoped block (cleaner) or comment the load-bearing drop:
  ```rust
  {
      let mut state = self.inner.borrow_mut();
      state.buf.clear();
      state.buf_start_base = 0;
  } // state dropped here; ChromRefBaseIter borrow_mut below would deadlock otherwise.
  ```

#### Mi6: [src/per_sample_pileup/ref_fetcher.rs:68](../../../../src/per_sample_pileup/ref_fetcher.rs#L68) — `STREAMING_REF_FILE_READ_CHUNK` is private but `STREAMING_REF_BUFFER_BYTES` is public
- **Categories:** defaults
- **Fix:** match the visibility — either both `pub(super)` or both `pub`.

#### Mi7: [src/per_sample_pileup/ref_fetcher.rs:182](../../../../src/per_sample_pileup/ref_fetcher.rs#L182) and `ManualEvictChromRefFetcher::for_contig` — inconsistent initial `buf_start_base` convention (`0` for streaming, `1` for manual-evict)
- **Categories:** defaults (Nit), reliability (cross-cat), smells (cross-cat)
- **Fix:** pick one — `0` as "no buffer yet" or `1` as "ready for the first 1-based base". Document the convention on the field. Touches both fetcher constructors.

#### Mi8: [src/per_sample_pileup/ref_fetcher.rs:189-237 and 651-657](../../../../src/per_sample_pileup/ref_fetcher.rs#L189-L237) — `new` takes `contig_name: String`, `for_contig` takes `&str`
- **Categories:** idiomatic
- **Fix:** subsumed by M12 (delete `new`). If kept: harmonise to `&str`; `to_string` at the call site is one line.

#### Mi9: [src/per_sample_pileup/ref_fetcher.rs:1346-1366](../../../../src/per_sample_pileup/ref_fetcher.rs#L1346-L1366) — `ManualEvictChromRefFetcher::prepend_backward` allocates a `Vec` scratch then `splice`s it in
- **Categories:** idiomatic (Minor), allocations (raised by perf review L4)
- Already tracked under the perf review. Listed here so the apply-fixes pass doesn't miss it.

#### Mi10–Mi16: Naming refinements
- **Categories:** naming (Minor)
- `legacy_io_error` reads as a noun; the `inner` field is generic across two types; `buf` field across two types; `StreamState` "Stream" is vague; missing `_BYTES` suffix on the file-read chunk const; `done` field is a bare adjective; the local `slot` in `WalkerLegacyAdapter::fetch` is generic. See [naming.md](../../../../tmp/review_2026-05-23_ref_fetcher/naming.md) for the full per-rename suggestions.

#### Mi17: [src/per_sample_pileup/ref_fetcher.rs:7-17](../../../../src/per_sample_pileup/ref_fetcher.rs#L7-L17) — Module doc says "Both adapters" but only one adapter exists
- **Categories:** smells
- **Fix:** delete the second adapter from the prose (the `SingleChromLegacyAdapter` was removed in the Step-2 migration).

#### Mi18: [src/per_sample_pileup/ref_fetcher.rs:106](../../../../src/per_sample_pileup/ref_fetcher.rs#L106) — `chrom_id` field doc says "unused under the new API"; field itself is set to `0` sentinel
- **Categories:** smells, defaults (cross-cat)
- **Fix:** delete the field entirely once M12 lands (the legacy `RefSeqFetcher` impl is what reads it). If kept, replace the prose with the actual contract: "set to 0 by `for_contig`; only read by the legacy `RefSeqFetcher` impl below".

#### Mi19: Module-wide — 2,224 LOC single file holds three fetcher types
- **Categories:** smells
- **Fix:** split into `streaming.rs` + `manual_evict.rs` + `walker_adapter.rs` + `chrom_ref_fetcher_trait.rs` + a sibling `tests/` directory. Net positive after M12 lands (the legacy surface goes away).

### Nits

Grouped, not enumerated:
- 50 `cargo fmt` diffs — run `cargo fmt`.
- Missing `#[must_use]` on constructors that return owned fetchers.
- Std-import grouping at lines 41–46 doesn't follow the project's "std → external → internal" convention.
- Several test-only `pub(super) fn` accessors on `ContigFai`/`StreamState` could be `#[cfg(test)] pub(super)`.
- A few `// XXX:` / `// TODO:` style comments without dates or issue references.
- `cargo audit` is not installed in the dev container — recommend `cargo install cargo-audit` in the Containerfile.

## 7. Out of scope observations

- The `unified_chrom_ref_fetcher` migration plan referenced in the module-doc (lines 24–39) is still in flight (Step 3 of the plan: drop the legacy `RefSeqFetcher` trait + the walker adapter). When that ships, M12 / M3 / M4 / M6 / M7 / M16 collapse.
- The `legacy_io_error`-and-walker pair shares a lifecycle with [stage1_pipeline.rs:26-123](../../../../src/pop_var_caller/stage1_pipeline.rs#L26-L123) (the only consumer of the legacy `RefSeqFetcher` surface today). When the walker migrates to `ChromRefFetcher`, this whole surface disappears.
- The pre-existing clippy warning in [record_encode.rs:359](../../../../src/var_calling/vcf_writer/record_encode.rs#L359) blocks `cargo clippy --lib --all-features -- -D warnings` from succeeding. Out of scope here; tracked under the vcf_writer review surface.

## 8. Missing tests to add now

Grouped by function under test. Each test name follows `function_returns_expected_on_condition`. Code stubs at the relevant per-category file (`tmp/review_2026-05-23_ref_fetcher/reliability.md`); below is the synthesized list.

### `StreamingChromRefFetcher::fetch`

- `streaming_fetch_at_production_buffer_size_rejects_backward_jump` — B3.
- `streaming_fetch_returns_invalid_index_error_on_zero_line_bases_fai` — B1.
- `streaming_fetch_returns_invalid_index_error_on_zero_line_width_fai` — B1.
- `streaming_fetch_returns_invalid_index_error_on_line_width_less_than_line_bases` — B1.
- `streaming_fetch_after_iter_bases_drop_starts_fresh_phase` — partial coverage exists; extend to assert the buffer reset is visible to a subsequent backward `fetch` (the iter-then-fetch interleave that M25 flags).
- `streaming_fetch_returns_io_error_with_contig_name_on_truncated_fasta` — B2 regression.

### `StreamingChromRefFetcher::fetch_into`

- `fetch_into_default_impl_clears_dst_before_writing` — M22.
- `fetch_into_returns_exact_length_from_streaming_slab_across_refill` — M22.
- `fetch_into_via_ref_forwarding_impl_works` — M22.

### `StreamingChromRefFetcher::iter_bases`

- `iter_bases_resets_buffer_on_drop_after_partial_walk` — Mi-tier reliability gap.
- `iter_bases_returns_chrom_ref_fetch_error_io_on_truncated_fasta` — error-path.
- `iter_bases_drop_does_not_double_panic_when_next_panicked_mid_borrow` — M15 (this is a stress test; can use `std::panic::catch_unwind`).

### `ManualEvictChromRefFetcher::fetch`

- `manual_evict_fetch_in_buffer_hit_is_alloc_free` — M25; uses `dhat-heap` or asserts buffer pointer equality.
- `manual_evict_fetch_after_append_forward_returns_correct_bytes` — happy path.
- `manual_evict_fetch_after_prepend_backward_returns_correct_bytes` — happy path.
- `manual_evict_evict_before_past_contig_end_clears_buffer` — M24.
- `manual_evict_prepend_backward_returns_io_error_on_truncated_source` — M24.

### `WalkerLegacyAdapter::fetch`

- `walker_adapter_serves_two_threads_alternating_chroms` — M23.
- `walker_adapter_returns_io_error_on_mutex_poison` — M4.

### `ContigFai::base_to_file_offset`

- `contig_fai_base_to_file_offset_panics_on_b_zero` — Mi3 (documents the precondition; subsumed by M17's `NonZeroU32`).

### `legacy_io_error`

- `legacy_io_error_preserves_contig_name_in_io_variant` — B2.
- `legacy_io_error_preserves_io_kind_for_out_of_bounds` — variant-mapping regression.

## 9. What's good

- **Per-worker ownership invariant is documented at the type definition** ([ref_fetcher.rs:113-119](../../../../src/per_sample_pileup/ref_fetcher.rs#L113-L119)) — the H1 perf-review fix rests on this invariant being load-bearing, and it's stated explicitly where future readers will look first.
- **`ChromRefFetchError` variants describe the failure mode at the error site, not just the mechanism** ([ref_fetcher.rs:498-538](../../../../src/per_sample_pileup/ref_fetcher.rs#L498-L538)) — `OutOfBounds`/`OutOfPattern`/`InvalidStart` carry the relevant coordinates and contig name in their fields; only `Io` is the catch-all (M20 flags that).
- **No `unsafe`, no `async`, no thread spawning, no channels** — the concurrency surface is intentionally small (`RefCell` for single-thread interior mutability + one `Mutex` on the walker adapter). The unsafe_concurrency review found three concrete polish items but no unsoundness.
- **`#![forbid(unsafe_code)]` at the crate root** ([src/lib.rs:16](../../../../src/lib.rs#L16)) — propagates to this module; refactor mistakes can't sneak in `unsafe { }`.
- **Test fixtures use real on-disk FASTA via `tempfile`** rather than mocks — the bench fixture pattern from [psp_reader_perf.rs](../../../../benches/psp_reader_perf.rs) is replicated correctly here for the streaming/manual-evict tests.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check` — expected to fail until Mi4 is applied (run `cargo fmt`).
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings` — currently fails on a pre-existing unrelated warning; fixing it is out of this review's scope but should be a CI gate.
- `./scripts/dev.sh cargo doc --no-deps --lib -- -D warnings` — currently fails on M1's stale link; CI is not gating this. Fix M1, then add `cargo doc --no-deps -- -D warnings` to the precommit script.
- `./scripts/dev.sh cargo test --lib -- per_sample_pileup::ref_fetcher` — currently 34/34 pass; will grow to ~55 after §8.
- `cargo install cargo-audit && ./scripts/dev.sh cargo audit` — requires installing the subcommand in the dev container first.

### Author response convention

Address each finding by its identifier (`B1`, `M5`, `Mi10`, etc.) with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the four open questions in §4 first — several Majors collapse depending on the answers.
