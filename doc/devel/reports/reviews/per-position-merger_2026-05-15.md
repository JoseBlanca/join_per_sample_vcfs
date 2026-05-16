# Code Review: per-position-merger
**Date:** 2026-05-15
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** New cohort-side multi-way per-position merger (commit `427d95f` on `main`)
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** PR-shaped commit adding a new public module to the crate. First cohort-side stage of the multi-sample SNP caller (linear-scan k-way merge over per-sample `.psp` iterators + companion chromosome-agreement helper).
- **Reviewed against:** commit `427d95f`, branch `main`.
- **In-scope files:**
  - [src/cohort/mod.rs](../src/cohort/mod.rs)
  - [src/cohort/per_position_merger.rs](../src/cohort/per_position_merger.rs)
  - [src/lib.rs](../src/lib.rs#L19) (one-line addition only — `pub mod cohort;`)
- **Deliberately out of scope:**
  - [doc/devel/implementation_plans/multi_way_per_position_iterator.md](../doc/devel/implementation_plans/multi_way_per_position_iterator.md) and [doc/devel/reports/implementations/multi_way_per_position_iterator_2026-05-15.md](../doc/devel/reports/implementations/multi_way_per_position_iterator_2026-05-15.md) — docs, not code.
  - Upstream `PspReader` / `PileupRecord` / `PspReadError` / `ParsedChromosome` — not changed by this commit.
  - Pre-existing rustfmt drift in [src/per_sample_caller/pileup/walker.rs](../src/per_sample_caller/pileup/walker.rs) and [src/per_sample_caller/psp/reader.rs](../src/per_sample_caller/psp/reader.rs) — surfaced under §7.
  - Pre-existing `cargo doc` failure in [src/per_sample_caller/pileup_to_psp.rs](../src/per_sample_caller/pileup_to_psp.rs) — surfaced under §7.
- **Categories dispatched** (one sub-agent per category, in parallel):
  - `reliability` — always applies.
  - `errors` — always applies; the module defines a new error type.
  - `naming` — always applies.
  - `defaults` — new public API surface.
  - `idiomatic` — always applies.
  - `refactor_safety` — always applies; contains a manual `Debug` impl.
  - `smells` — always applies.
  - `tooling` — crate-level change with `Cargo.toml`.
  - `extras` — public API on a per-position hot path; PR-shaped diff (intent-vs-code check).
  - **Skipped:** `unsafe_concurrency` — no `unsafe`, no `Arc`/`Mutex`/atomics/channels/`async`/thread spawning.

## 2. Verdict

**Approve-with-changes.** No Blockers — clippy + 516 tests + integration tests + example/bench builds all green. The Major findings cluster into three categories: (a) error-type design (the helper and the merger share one enum, and the `Reader` variant double-prints its cause), (b) one unjustified `.expect()` in production code, (c) test-coverage gaps for three correctness-load-bearing branches. (a) and (c) are harder to fix once downstream stages depend on the surface; recommend addressing before merging additional cohort-side work.

## 3. Execution status

Commands run inside the dev container (`./scripts/dev.sh`); output is quoted verbatim from the implementation phase where re-running would produce identical output:

| Command | Exit | Result |
|---|---|---|
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | `Finished \`dev\` profile [unoptimized + debuginfo] target(s) in 1.64s` — clean. |
| `cargo test --all-features --lib` | 0 | `test result: ok. 516 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.06s`. 18 new tests in `cohort::per_position_merger::tests`. |
| `cargo test --all-features --tests` | 0 | All integration tests pass. |
| `cargo build --examples` | 0 | Clean. |
| `cargo build --benches` | 0 | Clean. |
| `cargo doc --no-deps` | non-zero | Pre-existing failure in [src/per_sample_caller/pileup_to_psp.rs:5-6](../src/per_sample_caller/pileup_to_psp.rs#L5-L6) (`no item named \`psp\` in scope`). Unrelated to this commit — surfaced under §7. |
| `cargo fmt --check` | non-zero | In-scope files clean. Pre-existing drift in [src/per_sample_caller/pileup/walker.rs](../src/per_sample_caller/pileup/walker.rs) and [src/per_sample_caller/psp/reader.rs](../src/per_sample_caller/psp/reader.rs) — surfaced under §7. |

Commands not run:
- `cargo audit` — no new dependencies; tooling sub-agent confirmed `Cargo.toml` already lists `thiserror = "2.0.18"`.

Findings labeled "Needs verification": 0.

## 4. Open questions and assumptions

1. **Should `check_chromosome_agreement` errors live in their own enum?** The helper today returns `MergerError`, but only one variant (`ChromosomeMismatch`) is reachable from it. Splitting requires deciding whether the path-based opener helper (deferred per the plan) will wrap both via `#[from]` or expose them separately. — affects **M1**.
2. **Should the `chrom_id` field on `ChromosomeMismatch` be `Option<u32>` or stay `u32` with documented `0` sentinel?** Today the count-mismatch arm uses `chrom_id: 0` as a placeholder; downstream consumers may already be matching on this field. — affects **Mi1**.
3. **Is the `Locus` type worth promoting to a shared module now (e.g. `src/per_sample_caller/pileup/locus.rs` → re-export to cohort) or should the merger ship a private `struct Locus` and unify later?** The walker has its own `struct Locus`; promoting needs one upstream touch. — affects **M4**.
4. **Are downstream cohort stages expected to construct `PerPositionPileups` directly?** If only the merger ever produces these, locking construction with a `_private: ()` field is cheap insurance against future drift. If downstream stages will synthesise them (e.g. for fixtures), keep the fields open. — affects **Mi12**.

## 5. Top 3 priorities

1. **M1: Split `MergerError` into two error types.** `check_chromosome_agreement`'s signature advertises three impossible variants. Once callers start matching on the helper's `Result`, every fix becomes a breaking API change. → [src/cohort/per_position_merger.rs:54](../src/cohort/per_position_merger.rs#L54)
2. **M2: Drop `: {source}` from the `Reader` variant's `Display`.** The `#[source]` link already renders the cause in `anyhow`-style chains; today every log line carries the inner error twice. → [src/cohort/per_position_merger.rs:75](../src/cohort/per_position_merger.rs#L75)
3. **M6: Add the test for "refill error in a tied advance".** The most complex branch in `next()` (partial mutation of `heads`, then early-return Err, then `done` latch) is unverified. → [src/cohort/per_position_merger.rs:282-289](../src/cohort/per_position_merger.rs#L282-L289)

## 6. Findings

### Major

#### M1: [src/cohort/per_position_merger.rs:54](../src/cohort/per_position_merger.rs#L54) — Single `MergerError` enum funnels two unrelated operations
- **Confidence:** High
- **Categories:** errors (filed); naming, smells, idiomatic, defaults (cross-category notes)
- **Problem:** `MergerError` mixes the merger's own failures (`SampleCountMismatch`, `Reader`, `OutOfOrder`) with the standalone helper's only failure (`ChromosomeMismatch`). `pub fn check_chromosome_agreement<R: Read + Seek>(readers: &[PspReader<R>]) -> Result<Vec<ParsedChromosome>, MergerError>` therefore advertises three impossible variants. The errors checklist rule is "one error type per fallible operation, not per crate or module".
- **Why it matters:** Callers that match exhaustively (and `MergerError` is `#[non_exhaustive]`, so `_` arms are forced) write dead arms for three irrelevant variants. Splitting later becomes a breaking API change once downstream cohort stages start using the helper.
- **Suggested fix:**
  ```rust
  #[non_exhaustive]
  #[derive(Error, Debug)]
  pub enum ChromosomeAgreementError {
      #[error("sample {sample_idx} ({sample_name}) chromosome {chrom_id} disagrees with sample 0: {detail}")]
      Mismatch {
          sample_idx: usize,
          sample_name: String,
          chrom_id: u32,
          detail: String,
      },
  }
  ```
  Drop `ChromosomeMismatch` from `MergerError`. A future path-based opener that wants a single bubble-up type can wrap both via `#[from]` at its own boundary.

#### M2: [src/cohort/per_position_merger.rs:75](../src/cohort/per_position_merger.rs#L75) — `Reader` variant's `Display` duplicates the source's message
- **Confidence:** High
- **Categories:** errors
- **Problem:** `#[error("reader {sample_idx} ({sample_name}) failed: {source}")]` interpolates the inner `PspReadError` directly into the parent message. The `#[source]` attribute already exposes it through `Error::source()`.
- **Why it matters:** With this `Display`, an `anyhow`-style chain renderer prints the inner error twice — once as the suffix of the parent line, once as the chained cause. The clean parent-vs-cause separation that `#[source]` exists to provide is lost.
- **Suggested fix:**
  ```rust
  #[error("reader {sample_idx} ({sample_name}) failed")]
  Reader {
      sample_idx: usize,
      sample_name: String,
      #[source]
      source: Box<PspReadError>,
  },
  ```
  The variant still carries enough context to identify the failed operation on its own line; the cause renders via the source chain.

#### M3: [src/cohort/per_position_merger.rs:250](../src/cohort/per_position_merger.rs#L250) — Non-test `.expect()` without a structural fallback
- **Confidence:** High
- **Categories:** errors
- **Problem:** Inside `Iterator::next` the merger re-scans `heads` to find which reader's head matches `min_key`, with a panic safety net:
  ```rust
  let offender = self
      .heads
      .iter()
      .position(|h| h.as_ref().is_some_and(|r| (r.chrom_id, r.pos) == min_key))
      .expect("min_key was derived from a non-None head");
  ```
  The invariant holds today (the immediately preceding `iter().flatten()` loop set `min_key` from a non-`None` head and nothing `take`s between then and now). But the project's convention is either (a) a `// PANIC-FREE:` comment naming the invariant, or (b) restructuring to remove the panic site.
- **Why it matters:** A future refactor that reorders the scans, or introduces a mutation between them, would turn this into a production panic instead of a `Result`. The structural fix removes the risk entirely.
- **Suggested fix (preferred — structural):** Fold the offender lookup into the initial scan; no second pass, no `.expect()`.
  ```rust
  let mut min: Option<((u32, u32), usize)> = None;
  for (idx, head) in self.heads.iter().enumerate() {
      if let Some(record) = head {
          let key = (record.chrom_id, record.pos);
          if min.is_none_or(|(m, _)| key < m) {
              min = Some((key, idx));
          }
      }
  }
  let Some((min_key, offending_sample_idx)) = min else {
      self.done = true;
      return None;
  };
  // …and use `offending_sample_idx` in the OutOfOrder branch.
  ```
  **Alternative (minimum):** Add the marker comment naming the invariant.

#### M4: [src/cohort/per_position_merger.rs:130](../src/cohort/per_position_merger.rs#L130) — Raw `(u32, u32)` tuple where the crate has a `Locus` domain type
- **Confidence:** High
- **Categories:** naming
- **Problem:** Five sites (`last_emitted: Option<(u32, u32)>` at line 130, `min_key: Option<(u32, u32)>` at line 227, the `else`-bound rebinding at line 234, the monotonicity compare at line 241, and the `last_emitted = Some(min_key)` write at line 293) all use a bare `(u32, u32)` tuple to mean "genomic locus". The crate already has the concept named — `struct Locus { chrom_id: u32, pos: u32 }` in [src/per_sample_caller/pileup/walker.rs:259](../src/per_sample_caller/pileup/walker.rs#L259) and `type Locus = (usize, u64)` in [src/per_sample_caller/cram_input.rs:607](../src/per_sample_caller/cram_input.rs#L607). The project's naming standard (CLAUDE.md / feature-implementation skill) explicitly cites this pattern.
- **Why it matters:** Reading `if min_key <= last` requires mentally re-deriving that `.0` is chromosome and `.1` is position. The merger is the cohort-side anchor for this concept — leaving the tuple raw cements a third spelling for what should be the crate's shared vocabulary.
- **Suggested fix:** Either promote the walker's `Locus` struct to a shared `pub(crate)` module and import it, or define a module-local `struct Locus { chrom_id: u32, pos: u32 }` with `#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]`. See **Open Question 3** for the placement choice. Pairs naturally with **M3**'s structural fix and with **Mi6** (`min_key` → `min_locus`).

#### M5: [src/cohort/per_position_merger.rs:240-258](../src/cohort/per_position_merger.rs#L240-L258) — Missing test: equality branch of the strict-monotonicity check
- **Confidence:** High
- **Categories:** reliability
- **Problem:** The guard is `if min_key <= last`. Only the strict-less-than case is covered (`out_of_order_record_is_detected` uses positions 5 then 3). The equality case — a reader re-emitting the same key — is the canonical regression for "strictly increasing" vs "non-decreasing" and is not asserted. A future change that loosened `<=` to `<` would still pass the entire suite.
- **Why it matters:** Strict monotonicity is the merger's load-bearing invariant for every downstream stage that keys on `(chrom_id, pos)`. Silently allowing duplicates would corrupt downstream aggregation.
- **Suggested fix:** Add (full body in §8):
  ```rust
  #[test]
  fn duplicate_key_from_same_reader_is_out_of_order() { /* see §8 */ }
  ```

#### M6: [src/cohort/per_position_merger.rs:282-289](../src/cohort/per_position_merger.rs#L282-L289) — Missing test: refill error during a tied advance
- **Confidence:** High
- **Categories:** reliability
- **Problem:** The advance loop visits tied readers in `sample_idx` order. If reader *j* errors during refill while reader *k > j* is also at the min, the code partially mutates `heads` (`heads[j].take()`) and then early-returns `Err`. The existing `reader_error_mid_stream_surfaces_and_latches_done` uses two readers at *different* positions, so the "tied-set with error in the middle" path is unverified. This is the most complex branch in `next()`.
- **Why it matters:** A bug here would either resurrect a swallowed record or silently lose one. Reliability rules explicitly require coverage of "every error variant from every code path".
- **Suggested fix:** Add (full body in §8):
  ```rust
  #[test]
  fn refill_error_in_tied_advance_latches_done() { /* see §8 */ }
  ```

#### M7: [src/cohort/per_position_merger.rs:309-311](../src/cohort/per_position_merger.rs#L309-L311) — Missing test: documented "reference-string mismatch is non-fatal" invariant
- **Confidence:** High
- **Categories:** reliability
- **Problem:** `check_chromosome_agreement`'s doc commits: "*mismatched reference strings that happen to agree on every per-chromosome field are not fatal here.*" None of the tests constructs two headers whose `WriterHeader.reference` differ while every `ChromosomeEntry` agrees and asserts `Ok`.
- **Why it matters:** Callers feeding files produced against differently-named FASTA paths but identical contig MD5s rely on this. A silent regression that expanded the check to include `reference` would reject valid cohorts.
- **Suggested fix:** Add (full body in §8):
  ```rust
  #[test]
  fn chromosome_agreement_differing_reference_is_not_fatal() { /* see §8 */ }
  ```

### Minor

#### Mi1: [src/cohort/per_position_merger.rs:326](../src/cohort/per_position_merger.rs#L326) — `chrom_id: 0` is a sentinel placeholder when the chromosome count differs
- **Confidence:** High
- **Categories:** defaults (filed); errors, idiomatic, smells, naming, reliability (convergent cross-category notes — six sub-agents surfaced this)
- **Problem:** When chromosome counts disagree, the helper constructs `ChromosomeMismatch { chrom_id: 0, ... }`. The literal `0` is not a real chromosome id — it's a fallback because the variant's struct shape requires one. A caller matching `ChromosomeMismatch { chrom_id, .. }` will read `0` and conclude the divergence is on chromosome 0.
- **Why it matters:** Hidden default; downstream consumers either misattribute or have to parse `detail` to disambiguate, which the variant's own docstring tells them not to do. Multiple categories independently surfaced this as a design smell.
- **Suggested fix (preferred):** Change `chrom_id` to `Option<u32>`; the count-mismatch case carries `None`.
  ```rust
  ChromosomeMismatch {
      sample_idx: usize,
      sample_name: String,
      /// `None` when the divergence is a count mismatch
      /// (no specific chromosome to blame).
      chrom_id: Option<u32>,
      detail: String,
  },
  ```
  See **Open Question 2** for the type-change decision.

#### Mi2: [src/cohort/per_position_merger.rs:75](../src/cohort/per_position_merger.rs#L75) — `Reader` variant name is mechanism-flavored and ambiguous
- **Confidence:** Medium
- **Categories:** errors
- **Problem:** Named after the source object, not the failed operation. Emitted from two sites — initial prefetch in `new()` and refill during `next()` — with different recovery stories (construction-time failure vs mid-stream data/IO).
- **Suggested fix:** Split into two variants:
  ```rust
  PrefetchRecord { sample_idx, sample_name, source: Box<PspReadError> },
  RefillRecord  { sample_idx, sample_name, source: Box<PspReadError> },
  ```
  Or, minimum, rename to `AdvanceRecord` (operation-named) and document the prefetch/refill duality.

#### Mi3: [src/cohort/per_position_merger.rs:134](../src/cohort/per_position_merger.rs#L134) — `done: bool` should read as a predicate
- **Confidence:** High
- **Categories:** naming
- **Problem:** Bare adjective; the naming rule asks for `is_done` / `has_terminated` so use sites (`if self.done`) read as predicates.
- **Suggested fix:** Rename to `is_done` (or, given the doc highlights latch semantics, `is_terminated`). Propagate to the three use sites.

#### Mi4: [src/cohort/per_position_merger.rs:45](../src/cohort/per_position_merger.rs#L45) — `per_sample` field is a bare adjective phrase
- **Confidence:** Medium
- **Categories:** naming
- **Problem:** `per_sample` qualifies a missing noun ("per-sample what?"). The naming rule lists bare adjective phrases as forbidden.
- **Suggested fix:** Rename to `per_sample_records` (or `records_by_sample`). Propagate to the construction sites at 267, 294, 297 and to the test assertions.

#### Mi5: [src/cohort/per_position_merger.rs:124](../src/cohort/per_position_merger.rs#L124) — `heads` field is a bare common noun
- **Confidence:** Medium
- **Categories:** naming
- **Problem:** `heads: Vec<Option<PileupRecord>>` is documented as "peeked next record per reader". The plural alone doesn't carry that meaning at use sites like `self.heads[sample_idx].take()`.
- **Suggested fix:** Rename to `peeked_heads` (closest to the existing doc) or `peeked_records`.

#### Mi6: [src/cohort/per_position_merger.rs:227](../src/cohort/per_position_merger.rs#L227) — `min_key` names the structural role, not the domain concept
- **Confidence:** Medium
- **Categories:** naming
- **Problem:** The value is the minimum genomic locus across peeked heads; `min_key` reads as a generic min-finding placeholder.
- **Suggested fix:** Rename to `min_locus` (pairs with **M4**).

#### Mi7: [src/cohort/per_position_merger.rs:246](../src/cohort/per_position_merger.rs#L246) — `offender` is a bare participle
- **Confidence:** Medium
- **Categories:** naming
- **Problem:** The value is the sample index of the regressing reader, which is exactly what gets put into `MergerError::OutOfOrder.sample_idx` two lines later.
- **Suggested fix:** Rename to `offending_sample_idx`. (Goes away entirely if **M3**'s structural fix is taken.)

#### Mi8: [src/cohort/per_position_merger.rs:330](../src/cohort/per_position_merger.rs#L330) — `base` / `this` in the chromosome-agreement zip loop
- **Confidence:** Medium
- **Categories:** naming
- **Problem:** `for (chrom_id, (base, this)) in baseline.iter().zip(other.iter()).enumerate()` — `this` is a generic placeholder; the naming rule lists `this`/`thing`/`item` as forbidden.
- **Suggested fix:**
  ```rust
  for (chrom_id, (baseline_chrom, other_chrom)) in
      baseline.iter().zip(other.iter()).enumerate()
  { ... }
  ```

#### Mi9: [src/cohort/per_position_merger.rs:372](../src/cohort/per_position_merger.rs#L372) — Test helpers `rec` and `names` are abbreviated / generic
- **Confidence:** Medium
- **Categories:** naming
- **Problem:** `fn rec(...)` and `fn names(...)` are called from every test in the module. `rec` is non-Rust-standard abbreviation; `names` is a bare plural (names of what?). Sibling test helpers in `psp/reader.rs:1563` and `psp/writer.rs:897` use `record`.
- **Suggested fix:** Rename to `record(chrom_id, pos)` and `sample_names(n)`; propagate mechanically.

#### Mi10: [src/cohort/per_position_merger.rs:369-370](../src/cohort/per_position_merger.rs#L369-L370) — Test type aliases shadow `Iterator::Item`
- **Confidence:** Low
- **Categories:** naming
- **Problem:** `type Item = Result<PileupRecord, PspReadError>;` re-uses the iterator-protocol name `Item` for a different concept inside the test module.
- **Suggested fix:** Rename to `ReadResult` and `TestReaderIter`, or leave as-is given the narrow test-module scope.

#### Mi11: [src/cohort/per_position_merger.rs:202-212](../src/cohort/per_position_merger.rs#L202-L212) — Public accessors lack `///` doc comments
- **Confidence:** High
- **Categories:** idiomatic (filed); extras (convergent — same finding); defaults, reliability (cross-category notes)
- **Problem:** `sample_names`, `chromosomes`, `n_samples` are `pub` on a re-exported `pub struct`, but none carries a `///` summary. The rule asks for a one-line summary on every `pub` item.
- **Suggested fix:**
  ```rust
  /// Sample names in reader order, as supplied to [`Self::new`].
  pub fn sample_names(&self) -> &[String] { ... }

  /// Agreed-upon chromosome list (validated via
  /// [`check_chromosome_agreement`]). Metadata only; not
  /// length-checked against the readers.
  pub fn chromosomes(&self) -> &[ParsedChromosome] { ... }

  /// Number of samples; equals `per_sample.len()` on every emission.
  pub fn n_samples(&self) -> usize { ... }
  ```

#### Mi12: [src/cohort/per_position_merger.rs:38](../src/cohort/per_position_merger.rs#L38) — `PerPositionPileups` has open `pub` fields with a length invariant
- **Confidence:** Medium
- **Categories:** idiomatic (filed); smells (cross-category)
- **Problem:** Three `pub` fields plus a doc-stated invariant `per_sample.len() == n_samples()`. Today the merger is the only producer, but the type is `pub` — nothing structurally prevents external struct-literal construction with a wrong-sized vector. The rule asks for "control construction at the type level" via a private marker.
- **Why it matters:** Future cohort stages that index `per_sample[sample_idx]` against `PerPositionMerger::n_samples()` will turn a handwritten short vector into an out-of-bounds panic rather than a construction-time error.
- **Suggested fix (option A):**
  ```rust
  pub struct PerPositionPileups {
      pub chrom_id: u32,
      pub pos: u32,
      pub per_sample: Vec<Option<PileupRecord>>,
      _private: (),
  }
  impl PerPositionPileups {
      pub(crate) fn new(...) -> Self { ... }
  }
  ```
  See **Open Question 4**.

#### Mi13: [src/cohort/per_position_merger.rs:39](../src/cohort/per_position_merger.rs#L39) — `PerPositionPileups::chrom_id` is undocumented while siblings are
- **Confidence:** High
- **Categories:** extras
- **Problem:** Inside the same struct, `pos` (line 41) and `per_sample` (lines 43-44) have rustdoc; `chrom_id` (line 39) does not. Inconsistency reads as an oversight.
- **Suggested fix:**
  ```rust
  /// Index into [`PerPositionMerger::chromosomes`]. Matches
  /// [`PileupRecord::chrom_id`].
  pub chrom_id: u32,
  ```

#### Mi14: [src/cohort/per_position_merger.rs:141](../src/cohort/per_position_merger.rs#L141) — Manual `Debug` impl picks fields without exhaustive destructure
- **Confidence:** High
- **Categories:** refactor_safety
- **Problem:** The manual `Debug` impl reaches into `self` field-by-field. The struct has six fields; the impl prints four. A new field will silently fail to appear in `Debug` output because there's nothing to make the compiler complain.
- **Why it matters:** `Debug` drifts silently because no caller asserts on its output. Future debugging on a freshly added field will be misleading.
- **Suggested fix:** Destructure `Self` exhaustively (no `..`) so adding a field forces an edit here.
  ```rust
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
      let Self {
          readers: _,
          heads: _,
          sample_names,
          chromosomes,
          last_emitted,
          done,
      } = self;
      f.debug_struct("PerPositionMerger")
          .field("sample_names", sample_names)
          .field("n_chromosomes", &chromosomes.len())
          .field("last_emitted", last_emitted)
          .field("done", done)
          .finish()
  }
  ```

#### Mi15: [src/cohort/per_position_merger.rs:330](../src/cohort/per_position_merger.rs#L330) — Duplicated per-field comparison in `check_chromosome_agreement`
- **Confidence:** High
- **Categories:** smells
- **Problem:** The three field comparisons (`name`, `length`, `md5`) build a `ChromosomeMismatch` with the same `sample_idx`, `sample_name`, `chrom_id`, differing only in the `detail` string. Three near-identical occurrences — the smells rule's threshold for extraction.
- **Why it matters:** Each new field on `ParsedChromosome` requires another paste-and-edit; easy to drift one branch (wrong format string, swapped `this`/`base`).
- **Suggested fix:**
  ```rust
  let mismatch = |detail: String| MergerError::ChromosomeMismatch {
      sample_idx,
      sample_name: sample_name.clone(),
      chrom_id: chrom_id as u32,
      detail,
  };
  if base.name != this.name {
      return Err(mismatch(format!("name {:?} vs {:?}", this.name, base.name)));
  }
  // …
  ```

#### Mi16: [src/cohort/per_position_merger.rs:179](../src/cohort/per_position_merger.rs#L179) — Repeated `Some(Ok) / None / Some(Err)` shape and `Reader` construction
- **Confidence:** High
- **Categories:** smells
- **Problem:** Both `new()` (lines 179-191) and `Iterator::next` (lines 275-290) drain a reader via the same three-arm `match`, building the same `MergerError::Reader { sample_idx, sample_name: ..., source: Box::new(source) }` on the error arm. Two occurrences — borderline, but the error-construction half is mechanical and benefits most from extraction.
- **Suggested fix:** Extract only the error-construction half:
  ```rust
  fn reader_error(
      sample_idx: usize,
      sample_names: &[String],
      source: PspReadError,
  ) -> MergerError {
      MergerError::Reader {
          sample_idx,
          sample_name: sample_names[sample_idx].clone(),
          source: Box::new(source),
      }
  }
  ```

#### Mi17: [src/cohort/per_position_merger.rs:312-358](../src/cohort/per_position_merger.rs#L312-L358) — `check_chromosome_agreement` missing single-reader and non-zero-`chrom_id` tests
- **Confidence:** High
- **Categories:** reliability
- **Problem:** Boundary-class coverage rule asks for 0, 1, MAX, single-element cases. Empty (`readers.len() == 0`) is covered; `readers.len() == 1` (short-circuits the skip-1 loop) is not. Separately, every mismatch test puts the divergence on `chromosomes[0]`, so the `chrom_id` field never carries a non-zero value — a regression hard-coding `chrom_id: 0` instead of using the enumeration index would pass.
- **Suggested fix:** Add `chromosome_agreement_single_reader_returns_its_list` and `chromosome_agreement_reports_divergence_on_correct_chrom_id` (full bodies in §8).

#### Mi18: [src/cohort/per_position_merger.rs:202-212](../src/cohort/per_position_merger.rs#L202-L212) — Public accessors `sample_names()` / `chromosomes()` have no direct tests
- **Confidence:** High
- **Categories:** reliability
- **Problem:** Tests reach the merger via `n_samples()` and emitted `PerPositionPileups` shape. `sample_names()` and `chromosomes()` returning the slices handed at construction is not asserted; a refactor that swapped storage to a smaller projection would not break any test.
- **Suggested fix:** Add `accessors_round_trip_constructor_args` (full body in §8).

#### Mi19: [src/cohort/per_position_merger.rs:215-300](../src/cohort/per_position_merger.rs#L215-L300) — No property-based test for the merger's algebraic invariant
- **Confidence:** Medium
- **Categories:** reliability
- **Problem:** The merger is a pure function from `Vec<Vec<Result<PileupRecord, _>>>` (all-Ok, all-monotone) to `Vec<PerPositionPileups>`; the reliability rules ask for property-based tests on any function with an algebraic law. Hand-rolled examples cover up to 3 readers; they won't catch e.g. an off-by-one bug when N≥4 with varying tie groups.
- **Suggested fix:** Add a `proptest` target (`dev-dependencies: proptest = "1"`); full body in §8.

#### Mi20: [src/cohort/per_position_merger.rs:159-160](../src/cohort/per_position_merger.rs#L159-L160) — `chromosomes` constructor parameter's role is ambiguous in the doc
- **Confidence:** Medium
- **Categories:** convergent cross-category notes from errors, defaults, smells, naming
- **Problem:** The doc says `chromosomes` "is metadata only and is not length-checked against the readers". But the field is exposed via `chromosomes()` (Mi11), used by `check_chromosome_agreement` callers as the authoritative shared list, and named identically to the helper's return. The "metadata only" framing reads as if a per-sample length check were expected.
- **Why it matters:** Multiple sub-agents surfaced this as a docs/API tension. A caller reading the doc may wonder what relationship the merger expects between `chromosomes` and `readers`.
- **Suggested fix:** Rewrite the constructor doc to drop the "metadata only" framing and say what the parameter *is*: "The agreed-upon chromosome list — typically obtained from [`check_chromosome_agreement`]. The merger does not validate this against the readers' headers; that's the caller's responsibility." Pairs with **Mi11**'s accessor doc.

### Nits

Grouped (4 nits, none enumerated):

- [src/cohort/per_position_merger.rs:265](../src/cohort/per_position_merger.rs#L265) — `(0..self.n_samples()).map(|_| None).collect()` is idiomatically `vec![None; self.n_samples()]` (`Option<PileupRecord>` is `Clone` via `PileupRecord`'s derive). The inline comment explains why the *outer* loop iterates over a local, but the `Vec` construction immediately above doesn't need a range. Surfaced by **idiomatic**, **smells**, **refactor_safety**, **reliability**, **extras** — convergent style nit.
- [src/cohort/per_position_merger.rs:335, 343, 351](../src/cohort/per_position_merger.rs#L335) — `chrom_id as u32` casts a `usize` loop index without `try_from`. Harmless in practice (no header has 2³² chromosomes) but the project's defaults pass usually flags `as` casts; consider `u32::try_from(chrom_id).expect(...)` to document the bound.
- [src/cohort/per_position_merger.rs:248-250](../src/cohort/per_position_merger.rs#L248-L250) — `.expect("min_key was derived from a non-None head")` reads better as `unreachable!("...")` since the message is a logic invariant. (Subsumed by **M3** if the structural fix is taken.)
- [src/cohort/per_position_merger.rs:660-744](../src/cohort/per_position_merger.rs#L660-L744) — `check_chromosome_agreement` tests build both readers from `writer_header(_)`, which hardcodes `sample: "sample"`. The `sample_name` field of the emitted `ChromosomeMismatch` is `"sample"` for both readers, so the assertion can't distinguish. Single-line fix: `h1.sample = "B".to_string();` then `assert_eq!(sample_name, "B")`.

## 7. Out of scope observations

These pre-existing issues were surfaced by sub-agents but are not from this commit. Follow up separately.

- **[src/per_sample_caller/pileup_to_psp.rs:5-6](../src/per_sample_caller/pileup_to_psp.rs#L5-L6) — pre-existing `cargo doc` failure ("no item named `psp` in scope").** `.github/workflows/ci.yml:43-46` runs `cargo doc --no-deps --lib --all-features` with `RUSTDOCFLAGS: -D warnings`; this should already be breaking CI on `main`. Last touched by an unrelated commit. **Follow-up:** fix the intra-doc link in a separate PR.
- **Pre-existing `rustfmt` drift in [src/per_sample_caller/pileup/walker.rs](../src/per_sample_caller/pileup/walker.rs) and [src/per_sample_caller/psp/reader.rs](../src/per_sample_caller/psp/reader.rs).** Project toolchain pinned at 1.95; those files were last touched under an earlier rustfmt and have format diffs that `cargo fmt --check` flags. Unrelated to this commit. **Follow-up:** run `cargo fmt --` workspace-wide in a dedicated commit.
- **[Cargo.toml:26-32](../Cargo.toml#L26-L32) — narrow `[lints.clippy]` set.** Project intentionally restricts to `fallible_impl_from = "warn"` and `fn_params_excessive_bools = "warn"`; the comment says "re-enable when those modules get their own clippy pass". Pre-existing policy; flagging only as a future broadening candidate.
- **Hot-path `Vec<Option<...>>` allocation per emission** ([src/cohort/per_position_merger.rs:265-266](../src/cohort/per_position_merger.rs#L265-L266)). On the per-position hot path (~10⁹ emissions × N≈10³ samples). The implementation plan explicitly defers benchmarking — "benching it in isolation would over-fit a synthetic generator" — so this is a follow-up for the performance category to revisit once Stages 3-5 drive realistic load. **Possible fix:** pre-allocate an internal scratch buffer and swap it out into the emitted struct, eliminating the per-emission alloc/free.

## 8. Missing tests to add now

All tests live in [src/cohort/per_position_merger.rs](../src/cohort/per_position_merger.rs)'s existing `#[cfg(test)] mod tests` block. Grouped by function under test.

### `PerPositionMerger::next`

#### `duplicate_key_from_same_reader_is_out_of_order` (covers M5)
Pins the equality branch of the `min_key <= last` monotonicity check.
```rust
#[test]
fn duplicate_key_from_same_reader_is_out_of_order() {
    let a = iter_from(vec![Ok(rec(0, 5)), Ok(rec(0, 5))]);
    let mut merger = PerPositionMerger::new(vec![a], names(1), Vec::new()).unwrap();
    let first = merger.next().unwrap().unwrap();
    assert_eq!((first.chrom_id, first.pos), (0, 5));
    let err = merger.next().unwrap().unwrap_err();
    assert!(
        matches!(err, MergerError::OutOfOrder { chrom_id: 0, pos: 5, .. }),
        "got {err:?}"
    );
    assert!(merger.next().is_none(), "done must latch after OutOfOrder");
}
```

#### `refill_error_in_tied_advance_latches_done` (covers M6)
Pins the partial-mutation branch in the tied-readers advance loop.
```rust
#[test]
fn refill_error_in_tied_advance_latches_done() {
    // Readers 0, 1, 2 all tied at (0, 7). Reader 1 errors on refill.
    let a = iter_from(vec![Ok(rec(0, 7))]);
    let b = iter_from(vec![Ok(rec(0, 7)), Err(fake_err())]);
    let c = iter_from(vec![Ok(rec(0, 7))]);
    let mut merger =
        PerPositionMerger::new(vec![a, b, c], names(3), Vec::new()).unwrap();
    let err = merger.next().unwrap().unwrap_err();
    let MergerError::Reader { sample_idx, .. } = err else {
        panic!("expected Reader error, got {err:?}");
    };
    assert_eq!(sample_idx, 1);
    assert!(merger.next().is_none());
    assert!(merger.next().is_none());
}
```

#### `merger_next_handles_u32_max_position` (challenge test — boundary)
```rust
#[test]
fn merger_next_handles_u32_max_position() {
    let a = iter_from(vec![Ok(rec(0, u32::MAX - 1)), Ok(rec(0, u32::MAX))]);
    let b = iter_from(vec![Ok(rec(0, u32::MAX))]);
    let mut merger =
        PerPositionMerger::new(vec![a, b], names(2), Vec::new()).unwrap();
    let p0 = merger.next().unwrap().unwrap();
    assert_eq!(p0.pos, u32::MAX - 1);
    assert!(p0.per_sample[0].is_some() && p0.per_sample[1].is_none());
    let p1 = merger.next().unwrap().unwrap();
    assert_eq!(p1.pos, u32::MAX);
    assert!(p1.per_sample[0].is_some() && p1.per_sample[1].is_some());
    assert!(merger.next().is_none());
}
```

#### `merger_next_handles_one_empty_one_populated` (challenge test — exhausted-reader-not-requeried)
```rust
#[test]
fn merger_next_handles_one_empty_one_populated() {
    let empty: TestIter = iter_from(Vec::new());
    let b = iter_from(vec![Ok(rec(0, 1)), Ok(rec(0, 2))]);
    let mut merger =
        PerPositionMerger::new(vec![empty, b], names(2), Vec::new()).unwrap();
    for pos in 1..=2 {
        let p = merger.next().unwrap().unwrap();
        assert_eq!((p.chrom_id, p.pos), (0, pos));
        assert!(p.per_sample[0].is_none());
        assert!(p.per_sample[1].is_some());
    }
    assert!(merger.next().is_none());
}
```

### `PerPositionMerger::new`

#### `merger_new_accepts_empty_iterator_per_reader` (challenge test — all readers exhausted at prefetch)
```rust
#[test]
fn merger_new_accepts_empty_iterator_per_reader() {
    let a: TestIter = iter_from(Vec::new());
    let b: TestIter = iter_from(Vec::new());
    let mut merger =
        PerPositionMerger::new(vec![a, b], names(2), Vec::new()).unwrap();
    assert_eq!(merger.n_samples(), 2);
    assert!(merger.next().is_none());
    assert!(merger.next().is_none());
}
```

### `PerPositionMerger` accessors

#### `accessors_round_trip_constructor_args` (covers Mi18)
```rust
#[test]
fn accessors_round_trip_constructor_args() {
    let chroms = vec![ParsedChromosome {
        name: "chr1".to_string(),
        length: 100,
        md5: "0".repeat(32),
    }];
    let merger = PerPositionMerger::<TestIter>::new(
        Vec::new(),
        Vec::new(),
        chroms.clone(),
    )
    .unwrap();
    assert!(merger.sample_names().is_empty());
    assert_eq!(merger.chromosomes(), chroms.as_slice());
    assert_eq!(merger.n_samples(), 0);
}
```

### `check_chromosome_agreement`

#### `chromosome_agreement_differing_reference_is_not_fatal` (covers M7)
```rust
#[test]
fn chromosome_agreement_differing_reference_is_not_fatal() {
    let mut h0 = writer_header(1);
    let mut h1 = writer_header(1);
    h0.reference = "GRCh38.fa".to_string();
    h1.reference = "GRCh38_other_path.fa".to_string();
    let r0 = psp_reader_with_header(h0);
    let r1 = psp_reader_with_header(h1);
    let chroms = check_chromosome_agreement(&[r0, r1])
        .expect("reference-string drift must not be fatal per the doc-comment");
    assert_eq!(chroms.len(), 1);
}
```

#### `chromosome_agreement_single_reader_returns_its_list` (covers Mi17 part 1)
```rust
#[test]
fn chromosome_agreement_single_reader_returns_its_list() {
    let r0 = psp_reader_with_header(writer_header(2));
    let chroms = check_chromosome_agreement(std::slice::from_ref(&r0)).unwrap();
    assert_eq!(chroms.len(), 2);
}
```

#### `chromosome_agreement_reports_divergence_on_correct_chrom_id` (covers Mi17 part 2)
```rust
#[test]
fn chromosome_agreement_reports_divergence_on_correct_chrom_id() {
    let h0 = writer_header(3);
    let mut h1 = writer_header(3);
    h1.chromosomes[2].length = 42;
    let r0 = psp_reader_with_header(h0);
    let r1 = psp_reader_with_header(h1);
    let err = check_chromosome_agreement(&[r0, r1]).unwrap_err();
    let MergerError::ChromosomeMismatch { chrom_id, .. } = err else {
        panic!("expected ChromosomeMismatch, got {err:?}");
    };
    assert_eq!(chrom_id, 2);
}
```

#### `check_chromosome_agreement_count_mismatch_chrom_id_is_sentinel` (challenge test — pins Mi1's sentinel until Mi1 is resolved)
```rust
#[test]
fn check_chromosome_agreement_count_mismatch_chrom_id_is_sentinel() {
    let h0 = writer_header(1);
    let h1 = writer_header(2);
    let r0 = psp_reader_with_header(h0);
    let r1 = psp_reader_with_header(h1);
    let err = check_chromosome_agreement(&[r0, r1]).unwrap_err();
    let MergerError::ChromosomeMismatch { chrom_id, detail, .. } = err else {
        panic!("expected ChromosomeMismatch, got {err:?}");
    };
    assert_eq!(chrom_id, 0, "count-mismatch arm uses chrom_id=0 sentinel");
    assert!(detail.contains("count"));
}
```

### Cross-cutting

#### `merger_emits_strict_sorted_union` (covers Mi19 — property)
Requires `proptest = "1"` as a `dev-dependency`.
```rust
proptest! {
    #[test]
    fn merger_emits_strict_sorted_union(
        inputs in proptest::collection::vec(
            proptest::collection::vec(0u32..200, 0..30)
                .prop_map(|mut v| { v.sort_unstable(); v.dedup(); v }),
            0..5,
        ),
    ) {
        let readers: Vec<TestIter> = inputs
            .iter()
            .map(|positions| {
                iter_from(positions.iter().map(|&p| Ok(rec(0, p))).collect())
            })
            .collect();
        let n = readers.len();
        let merger =
            PerPositionMerger::new(readers, names(n), Vec::new()).unwrap();
        let emitted: Vec<u32> = merger.map(|r| r.unwrap().pos).collect();
        let mut expected: Vec<u32> = inputs.iter().flatten().copied().collect();
        expected.sort_unstable();
        expected.dedup();
        prop_assert_eq!(emitted, expected);
    }
}
```

## 9. What's good

- **Inline rationale on `Box<PspReadError>`** at [src/cohort/per_position_merger.rs:67-72](../src/cohort/per_position_merger.rs#L67-L72) — the hot-path Result-size concern is documented at the variant where the choice was made, so the next reader doesn't have to reverse-engineer it.
- **`#[non_exhaustive]` on `MergerError`** ([src/cohort/per_position_merger.rs:52](../src/cohort/per_position_merger.rs#L52)) — forward-compatible with the plan's deferred heap-merger and path-opener variants without breaking downstream matches.
- **Synthetic `std::iter`-based test fixtures** ([src/cohort/per_position_merger.rs:380-407](../src/cohort/per_position_merger.rs#L380-L407)) — the merger's generic-iterator API lets tests drive it without on-disk `.psp` round-trips; only `check_chromosome_agreement` needs real readers. Clean separation of concerns.
- **Latch-on-error contract** ([src/cohort/per_position_merger.rs:288](../src/cohort/per_position_merger.rs#L288)) consistent with the walker's terminate-on-first-error design — one mental model across the pipeline.
- **Implementation report** at [doc/devel/reports/implementations/multi_way_per_position_iterator_2026-05-15.md](../doc/devel/reports/implementations/multi_way_per_position_iterator_2026-05-15.md) documents the silent choices (defensive length check, boxed error variant, `Option`-less `ChromosomeMismatch.detail`) alongside the plan, so the rationale survives even after the implementation report ages.

## 10. Commands to re-verify

The orchestrator ran these — re-run before merging the fixes:

- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --all-features --lib`
- `./scripts/dev.sh cargo test --all-features --tests`
- `./scripts/dev.sh cargo build --examples`
- `./scripts/dev.sh cargo build --benches`

New invocations the review introduces:

- `./scripts/dev.sh cargo test --all-features --lib cohort::per_position_merger::tests::duplicate_key_from_same_reader_is_out_of_order` (and the other new test names above) — run individually as you add them.
- If **Mi19** is taken: `./scripts/dev.sh cargo test --all-features --lib cohort::per_position_merger::tests::merger_emits_strict_sorted_union` after adding `proptest` to `dev-dependencies`.

---

### Author response convention

Address each finding by its identifier (e.g. "M1 fixed in <hash>", "Mi8 deferred to <issue>"). Answer open questions in §4 first.
