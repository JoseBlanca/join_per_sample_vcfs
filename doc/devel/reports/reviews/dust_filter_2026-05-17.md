# Code Review: dust_filter
**Date:** 2026-05-17
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** PR diff `a00148a` — Stage 3 sdust low-complexity filter v1 (new module `src/var_calling/dust_filter.rs`).
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** PR diff `a00148a` adding `src/var_calling/dust_filter.rs` (Stage 3 sdust low-complexity filter v1).
- **Reviewed against:** commit `a00148a` on `main`.
- **In-scope files:**
  - [src/var_calling/dust_filter.rs](../../../src/var_calling/dust_filter.rs)
  - [src/var_calling/mod.rs](../../../src/var_calling/mod.rs) (1-line `pub mod dust_filter;`)
- **Deliberately out of scope:**
  - `PROJECT_STATUS.md`, the implementation report, and the plan — navigation/docs, not code.
  - The cloned `sdust/` reference C source — vendored read-only for porting reference, not project tooling.
  - All other modules — pre-existing, untouched.
- **Categories dispatched** (9, one sub-agent each; `unsafe_concurrency` skipped because the module has no `unsafe` / `Arc` / `Mutex` / atomics / channels / `async` / thread spawning):
  - `reliability` — always.
  - `errors` — always.
  - `naming` — always.
  - `defaults` — new public config type (`DustFilterConfig`) with `Default` impl.
  - `idiomatic` — always.
  - `refactor_safety` — always.
  - `smells` — always.
  - `tooling` — adds new module to a crate.
  - `extras` — parser/scanner, hot path, accepts untrusted reference bytes, stable output, "diff matches stated intent" against the plan and impl report.

## 2. Verdict

**Request-changes.** No Blockers, but the `pileups.pos - 1` underflow at line 620 (**M1**) is a silent wrong-result mode on a documented contract violation, and the error-Display chain doubling (**M6**) plus missing `// PANIC-FREE:` comment (**M7**) are project-convention violations the previous Stage 6 reviews flagged in similar code. The test coverage gaps on `sdust_mask` (no property test — **M2**, no boundary test on the half-open `is_masked` semantics — **M3**, no out-of-order regression test — **M4**) leave the load-bearing structural invariants of the algorithmic core unprotected against future refactors. **M5** is an honest mis-spec in the impl report: the test claimed to pin the strict-`>` threshold boundary actually pins nothing.

The shipped code passes `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, and `cargo test --all-targets --all-features` cleanly; the algorithmic port from `lh3/sdust` was verified byte-for-byte against the vendored C source by the extras agent and no port-fidelity issues were found. The findings below are about correctness contracts, error discipline, and test coverage — not the algorithm itself.

## 3. Execution status

| Command | Exit | Result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean, no warnings |
| `cargo test --all-targets --all-features` | 0 | 693 lib + 109 integration tests pass (per binary: 693, 3, 5, 25, 26, 17, 2, 7, 4, 8, 12; three empty test binaries) |
| `cargo doc --no-deps` | not run | flagged at **Mi21** — `Cargo.toml` declares `broken_intra_doc_links = "deny"` but doc-build was not exercised on a PR with 10+ intra-doc links. |
| `cargo audit` | not run | project does not use cargo-audit. |

Findings labeled "Needs verification": 0.

## 4. Open questions and assumptions

1. **Tracing policy.** The defaults checklist requires a `tracing::debug!` event when a `Default`-derived config is used; the user-memory entry "No logs — promote to typed errors" prefers typed errors over logs and notes that "logs don't get read." Does the project want `tracing::debug!` on filter construction (the rule), or does the broader no-logs preference override? Affects **Mi6** (was Major in the defaults agent's report; downgraded here pending the answer).
2. **`window` upper bound.** Is the `4096` cap on `window` defended by anything concrete (cohort scale, FASTA-line wrap, allocator budget), or is it intuition? Affects **Mi7**.
3. **`sdust_mask` panic vs `Result` symmetry.** `DustFilterConfig::new` rejects `window` outside `3..=4096` via `Result`; `sdust_mask` asserts `window >= 3` and panics. Two enforcement paths for the same invariant. Affects **Mi24**. (Compatible options: (a) keep `sdust_mask` as a "validated precondition" function with the assert and document the precondition by reference to `DustFilterConfig::new`; (b) make `sdust_mask` return `Result`.)
4. **Out-of-order upstream.** `PerPositionMerger` documents strict monotonic per-chromosome emission. If the contract were violated, the filter would silently produce wrong pass/skip decisions because the sweep advances monotonically. Should the filter assume the contract and document it (no defensive check) or actively reject violations? Affects **M4**.

## 5. Top 3 priorities

1. **M1 — `pileups.pos - 1` underflow on `pos == 0`.** Silent wrong filtering in release, panic in debug. Add `InvalidPos` variant + `checked_sub`. Loadbearing for the project's "no silent failures" convention. ([dust_filter.rs:620](../../../src/var_calling/dust_filter.rs#L620))
2. **M2 — No property test on `sdust_mask`.** Pure function over a structured input domain; six golden snippets are not enough to protect the "sorted, non-overlapping, in-bounds" structural invariant the iterator layer trusts. Add a proptest (or seeded-random invariant loop). ([dust_filter.rs:332](../../../src/var_calling/dust_filter.rs#L332))
3. **M6 — Error `Display` doubles the source chain.** `#[error("upstream merger failed: {0}")]` and `#[error("reference fetch failed for chrom {chrom_id}: {source}")]` interpolate the `#[source]` cause into the parent's message, breaking anyhow-style chain rendering. Drop the `{0}`/`{source}` interpolation. ([dust_filter.rs:463](../../../src/var_calling/dust_filter.rs#L463))

## 6. Findings

### Blocker

_None._

### Major

#### M1: src/var_calling/dust_filter.rs:620 — `pileups.pos - 1` underflows silently on `pos == 0`

**Categories:** errors, reliability (convergent)
**Confidence:** High
**Problem:** `let pos0 = pileups.pos - 1;` at the upstream-to-sdust coordinate boundary uses unchecked `u32` subtraction. `PileupRecord::pos` is documented 1-based but the field is `pub` and the type system does not enforce non-zero. If an upstream item arrives with `pos == 0` (malformed mock, future merger bug, hand-crafted PSP), the subtraction wraps to `u32::MAX` in release mode (panics in debug). `is_masked` then checks `start <= u32::MAX < finish`, which is `false` for any realistic mask — the record silently passes through and masking is effectively disabled for that record. Per the user-memory note "No logs — promote to typed errors," contract violations on upstream input are exactly the case a typed error variant exists to surface.
**Suggested fix:**
```rust
// New variant on DustFilterError:
#[error("upstream emitted invalid 1-based pos 0 on chrom_id {chrom_id}")]
InvalidPos { chrom_id: u32 },

// In DustFilter::next, replace `let pos0 = pileups.pos - 1;`:
let pos0 = match pileups.pos.checked_sub(1) {
    Some(p) => p,
    None => {
        self.done = true;
        return Some(Err(DustFilterError::InvalidPos {
            chrom_id: pileups.chrom_id,
        }));
    }
};
```
Plus a regression test `filter_surfaces_invalid_pos_zero_and_latches` matching the existing latch-test pattern.

#### M2: src/var_calling/dust_filter.rs:332 — `sdust_mask` lacks property/fuzz coverage

**Category:** reliability
**Confidence:** High
**Problem:** `sdust_mask` is a pure function over a structured input domain (ACGT/N byte slice + two `u32` knobs) and is the entire algorithmic guts of Stage 3 — every Stage 4-bound record traces through it. The reliability checklist explicitly flags this class as requiring property-based or fuzz tests. Today's coverage is six golden snippets + eleven shape assertions + one panic check; nothing holds for *all* inputs. A future refactor of `shift_window` / `save_masked_regions` / `find_perfect` could violate the "sorted, non-overlapping, in-bounds" structural invariant that `is_masked` trusts (lines 582-588), and the existing tests would not catch it on inputs outside the golden set.
**Suggested fix:**
```rust
#[test]
fn sdust_mask_invariants_hold_on_random_inputs() {
    use proptest::prelude::*;
    proptest!(|(
        seq in proptest::collection::vec(
            prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'),
                        Just(b'T'), Just(b'N'),
                        Just(b'a'), Just(b'c'), Just(b'g'), Just(b't')],
            0usize..512),
        window in 3u32..=128,
        threshold in 1u32..=200,
    )| {
        let out = sdust_mask(&seq, window, threshold);
        // Each interval is non-empty and in-bounds.
        for &(s, f) in &out {
            prop_assert!(s < f);
            prop_assert!(f as usize <= seq.len());
        }
        // Sorted by start, non-overlapping.
        for w in out.windows(2) {
            prop_assert!(w[0].1 <= w[1].0);
        }
        // Perfect intervals never span an N.
        for &(s, f) in &out {
            for b in &seq[s as usize..f as usize] {
                prop_assert!(SEQ_NT4[*b as usize] < 4);
            }
        }
    });
}
```
If `proptest` is not yet a dev-dep, the equivalent loop over ~256 inputs using the existing `xorshift` generator in tests is acceptable.

#### M3: src/var_calling/dust_filter.rs:581 — Half-open boundary semantics of `is_masked` not pinned by any test

**Category:** reliability
**Confidence:** High
**Problem:** `is_masked` treats intervals as half-open: `start <= pos0 && pos0 < finish` (line 587) plus a sweep-advance `finish <= pos0` (line 583). Both off-by-ones are risky in different directions. The module's own doc-comment (lines 44-49) calls out the 1-based→0-based conversion as the load-bearing invariant. Existing tests cover whole-interval cases (every position inside or every position outside); none probe `pos = start` or `pos = finish`.
**Suggested fix:** Use the `homopolymer_with_flanks` golden snippet (mask `(50, 100)`) and probe 1-based positions 50, 51, 100, 101:
```rust
#[test]
fn filter_half_open_boundary_pos_at_start_masked_pos_at_finish_passes() {
    let bases = GOLDEN_SNIPPETS[0].bases.to_vec();
    let len = bases.len() as u32;
    let chromosomes = vec![chrom("chr1", len)];
    let fetcher = StubFetcher::new(vec![bases]);
    // pos = 50 (0-based 49, just before mask) -> pass
    // pos = 51 (0-based 50, first masked base) -> drop
    // pos = 100 (0-based 99, last masked base) -> drop
    // pos = 101 (0-based 100, just past mask) -> pass
    let upstream = positions_iter(vec![(0, 50), (0, 51), (0, 100), (0, 101)]);
    let filter = DustFilter::new(
        upstream, fetcher, chromosomes, DustFilterConfig::default(),
    );
    let kept: Vec<u32> = filter.map(|r| r.unwrap().pos).collect();
    assert_eq!(kept, vec![50, 101]);
}
```

#### M4: src/var_calling/dust_filter.rs:599 — `DustFilter::next` is not regression-tested against out-of-order upstream positions

**Category:** reliability
**Confidence:** High
**Problem:** `is_masked` advances `self.sweep` monotonically (lines 582-584). The `PerPositionMerger` contract is "strict monotonic per chromosome," but this assumption is not documented in `DustFilter` and not tested. An upstream that regresses (test mock, future bug, contract drift) produces silently wrong pass/skip decisions — the worst failure mode. The other upstream-contract violations (`OutOfOrder`, unknown chrom_id, ref-fetch failure) all have latching tests; this one is the gap. See also **Open question 4**.
**Suggested fix:** Either (preferred) add a defensive check at the top of `next` comparing `pileups.pos` to a remembered `last_pos` on the same `chrom_id`, surfacing a typed `DustFilterError::OutOfOrder` and latching; or at minimum add a regression test that pins the current behavior so the next maintainer does not accidentally rely on tolerance:
```rust
#[test]
fn filter_is_undefined_on_out_of_order_within_chrom() {
    // Documents current behavior: filter consumes out-of-order
    // positions without error, but its mask decision for the
    // regressing position may be wrong. If you add the
    // defensive check, expect Err here instead.
    let seq = vec![b'A'; 128];
    let chromosomes = vec![chrom("chr1", 128)];
    let fetcher = StubFetcher::new(vec![seq]);
    let upstream = positions_iter(vec![(0, 64), (0, 5)]);
    let filter = DustFilter::new(
        upstream, fetcher, chromosomes, DustFilterConfig::default(),
    );
    let out: Vec<_> = filter.collect();
    assert_eq!(out.len(), 0); // current: both positions dropped
}
```

#### M5: src/var_calling/dust_filter.rs:759-765 — `sdust_high_threshold_disables_masking` does not pin the strict-`>` boundary

**Category:** extras
**Confidence:** High
**Problem:** Both the plan and the impl report claim a dedicated test pins the `>` vs `≥` threshold semantics. The shipped test runs a 128-bp poly-A at `T = 5000`, where the density bar `T/10 = 500` is mathematically unreachable by any window — the test would pass under `>`, `≥`, `<`, or `<=` flipped one way. It is a sensitivity smoke ("does the knob exist?"), not a boundary test. The strict-`>` invariant is currently protected only by the golden vectors, which have no design-time guarantee of including a boundary-density case.
**Suggested fix:** Either (a) construct a snippet whose top candidate's density equals exactly `T/10` and assert it is **not** masked, paired with a `T/10 + ε` snippet asserted **is**; or (b) retract the plan/report claim and add an explicit comment to `sdust_high_threshold_disables_masking` clarifying that strict-`>` is only protected by the golden vectors. (a) is more defensive; (b) is honest about what the test does.

#### M6: src/var_calling/dust_filter.rs:463-484 — `DustFilterError` Display strings double the source chain

**Category:** errors
**Confidence:** High
**Problem:**
```rust
#[error("upstream merger failed: {0}")]
Upstream(#[source] Box<PerPositionMergerError>),

#[error("reference fetch failed for chrom {chrom_id}: {source}")]
RefFetch { chrom_id: u32, #[source] source: std::io::Error },
```
`{0}` and `{source}` interpolate the cause's `Display` into the parent message. A consumer rendering the chain (or `anyhow`) will print the cause twice — once inline in the parent string and once as the next chain link. The errors-checklist convention is: let `source()` carry the cause; the parent's `Display` describes only what *this* layer was trying to do.
**Suggested fix:**
```rust
#[error("upstream merger failed")]
Upstream(#[source] Box<PerPositionMergerError>),

#[error("reference fetch failed for chrom {chrom_id}")]
RefFetch {
    chrom_id: u32,
    #[source]
    source: std::io::Error,
},
```

#### M7: src/var_calling/dust_filter.rs:182 — `expect("deque non-empty at cap")` lacks the `// PANIC-FREE:` comment the project requires

**Category:** errors
**Confidence:** High
**Problem:** The only `expect` in non-test code. The invariant is real (we just tested `len() >= cap` and `cap >= 1`), but the project convention — surfaced in the contamination_estimation and posterior_engine reviews — pairs every `expect` / `unwrap` with an inline `// PANIC-FREE:` comment naming the invariant. Mechanical fix; same shape as elsewhere in the codebase.
**Suggested fix:**
```rust
// PANIC-FREE: window.len() >= cap (just tested) and cap >= 1
// (window_cap >= SD_WLEN = 3 validated by DustFilterConfig::new).
let s = self.window.pop_front().expect("deque non-empty at cap") as usize;
```
Or, preferred, restructure to `if let Some(s) = self.window.pop_front()` with `else` unreachable and a comment — same shape, no panic in the type.

#### M8: src/var_calling/dust_filter.rs:184, 188, 202, 204-206, 271 — Unchecked counter subtractions rely on undocumented invariants

**Categories:** errors, smells (convergent on `window.len() - big_l`)
**Confidence:** Medium
**Problem:** `shift_window`, `save_masked_regions`, and `find_perfect` decrement `u32`s with plain `-=` and use `usize` subtractions (`self.window.len() - self.big_l as usize` at 202, `win_len - self.big_l - 1` at 271) without `checked_sub` or `// INVARIANT:` comments. The algorithm is correct as written; the next reader who edits a branch in `shift_window` has nothing in the code that flags a broken invariant. A wrap in release mode would not panic — it would silently distort scores and produce wrong masks, the worst failure mode for a silent filter. The trim invariant `big_l <= window.len()` is documented in `find_perfect`'s comment block (lines 263-266) but not at the dependent sites in `shift_window`/`save_masked_regions`.
**Suggested fix:** Add `// INVARIANT:` comments at each decrement site naming the matching increment site that guarantees non-zero, or switch to `debug_assert!(self.cw[s] >= 1, ...)` paired with plain ops. Apply uniformly across lines 184, 188, 202, 204-206, 271. For example at line 184:
```rust
// INVARIANT: `cw[s]` was incremented when triplet `s` entered the
// window via shift_window's push; it has not been decremented
// since because each triplet is popped at most once on eviction.
debug_assert!(self.cw[s] >= 1);
self.cw[s] -= 1;
```
The `debug_assert!` form provides loud failure under tests without adding release-time cost; the comment form alone is also acceptable per the errors checklist.

### Minor

#### Mi1: src/var_calling/dust_filter.rs:428,482 — Magic upper bound `4096` repeated; lower bound `3` restates `SD_WLEN`

**Category:** naming
**Confidence:** High
**Problem:** `DustFilterConfig::new` validates `3..=4096`; the error variant prints `"must be in 3..=4096"`. The `4096` upper bound has no named constant; the `3` lower bound restates `SD_WLEN`. Two sites must agree by hand.
**Suggested fix:** Introduce `const MAX_DUST_WINDOW: u32 = 4096;` with a doc-comment justifying the cap, then use `(SD_WLEN..=MAX_DUST_WINDOW)` in `new` and a parameterised `#[error(...)]` template:
```rust
#[error("invalid sdust window {window}: must be in {}..={}", SD_WLEN, MAX_DUST_WINDOW)]
InvalidWindow { window: u32 },
```

#### Mi2: src/var_calling/dust_filter.rs:146 — Field `res` is the forbidden generic noun `result`

**Category:** naming
**Confidence:** High
**Problem:** `SdustState.res: SdustIntervals` uses `res` (abbreviation of `result`) as a field name. The naming rule explicitly bans the generic noun. Unlike `cw`/`cv`/`rw`/`rv`/`big_l`, `res` is not part of the algorithm's formula — just a buffer name, no C-traceability claim.
**Suggested fix:** Rename to `masked_intervals: SdustIntervals` and update the two use sites.

#### Mi3: src/var_calling/dust_filter.rs:518 — Boolean field `done: bool` is a bare adjective

**Category:** naming
**Confidence:** High
**Problem:** Naming rule requires predicate form on booleans (`is_finished`, `has_terminated`). The field's doc clarifies it is specifically the "iterator-exhausted-or-errored" latch; readers seeing `self.done = true` at call sites get no semantic anchor.
**Suggested fix:** Rename `done` → `is_finished` throughout.

#### Mi4: src/var_calling/dust_filter.rs:641 — Test helper `high_complexity` is a bare adjective

**Category:** naming
**Confidence:** High
**Problem:** `fn high_complexity(n: usize) -> Vec<u8>` — modifier without noun. Used at ~9 call sites in tests.
**Suggested fix:** `fn high_complexity_bases(n: usize) -> Vec<u8>`.

#### Mi5: src/var_calling/dust_filter.rs:952 — Test helper `chrom` is a noun used as a constructor

**Category:** naming
**Confidence:** Medium
**Problem:** Function names should be verbs. `chrom(name, length) -> ParsedChromosome` reads like an accessor.
**Suggested fix:** `fn make_chrom(name: &str, length: u32) -> ParsedChromosome`.

#### Mi6: src/var_calling/dust_filter.rs:443-454 — `Default` applies behaviorally-significant values without a `tracing` event

**Category:** defaults
**Confidence:** High
**Problem:** The defaults checklist requires a `tracing::debug!` event when a `Default`-derived config is used, identifying field, value, and provenance. The module emits no tracing events anywhere. This finding is **downgraded** from Major because the user-memory entry "No logs — promote to typed errors. Logs don't get read." pushes back on the same axis. See **Open question 1**.
**Suggested fix:** Resolve the tracing policy at the project level. If logs are wanted, add a `ConfigSource::{Default, Explicit}` marker to `DustFilterConfig` and emit a `tracing::debug!` from `DustFilter::new`. If logs are not wanted, document the trade-off in the `Default` impl's doc comment ("the resolved config can be inspected at runtime via `DustFilter::config()`; we do not log it").

#### Mi7: src/var_calling/dust_filter.rs:421-432 — Upper bound `4096` on `window` not tied to a concrete source

**Category:** defaults
**Confidence:** Medium
**Problem:** Doc says "guard against accidental misuse — the spec's largest reasonable choice is well under this." Unlike the default value (`64`, cited to `lh3/sdust`/minimap2) and the lower bound (`3 = SD_WLEN`, algorithmic requirement), the upper bound has no citation. A future contributor hitting the bound has nothing concrete to anchor a decision on. See **Open question 2**.
**Suggested fix:** Either tighten the doc to cite a concrete reason (e.g. "1024 is already 16× the minimap2 default; nothing in the SNP-calling pipeline has been tuned beyond that"), or drop the upper bound entirely and rely on `u32` correctness (the algorithm is correct for any `window >= 3` modulo the `u64`-widening guards already in place).

#### Mi8: src/var_calling/dust_filter.rs:415-418 — `DustFilterConfig` fields `window` and `threshold` lack `///` doc comments

**Category:** defaults
**Confidence:** High
**Problem:** Public type, private fields, no field-level docs. Documentation lives only on the constants and `::new`; rustdoc on the struct shows two undescribed `u32` fields.
**Suggested fix:**
```rust
pub struct DustFilterConfig {
    /// sdust `W` — maximum subinterval length the algorithm
    /// considers. Defaults to [`DEFAULT_DUST_WINDOW`]; validated to
    /// `3..=4096` in [`Self::new`].
    window: u32,
    /// sdust `T` — score-density threshold. Masking condition is
    /// `10·s > T·L`; the density bar is `T/10` triplet-pairs per
    /// triplet. Defaults to [`DEFAULT_DUST_THRESHOLD`]. Unbounded;
    /// a very large value effectively disables masking.
    threshold: u32,
}
```

#### Mi9: src/var_calling/dust_filter.rs:434,438,548 — Public accessors lack `///` doc comments

**Categories:** idiomatic, extras (convergent)
**Confidence:** High
**Problem:** `DustFilterConfig::window`, `DustFilterConfig::threshold`, and `DustFilter::config` are `pub` with no `///`. `DustFilterConfig::new` has an inline paragraph but no `# Errors` block.
**Suggested fix:**
```rust
impl DustFilterConfig {
    /// Build a validated config. `window` must lie in `3..=4096`.
    ///
    /// # Errors
    /// Returns [`DustFilterError::InvalidWindow`] if `window` is
    /// outside that range.
    pub fn new(window: u32, threshold: u32) -> Result<Self, DustFilterError> { ... }

    /// Validated sdust window size (`W`).
    pub fn window(&self) -> u32 { self.window }

    /// Validated sdust score-density threshold (`T`).
    pub fn threshold(&self) -> u32 { self.threshold }
}

impl<I, F> DustFilter<I, F> ... {
    /// Returns the validated configuration in effect.
    pub fn config(&self) -> DustFilterConfig { self.config }
}
```

#### Mi10: src/var_calling/dust_filter.rs:230-240 — `let mut merged` bool collapses to `if let … else`

**Category:** idiomatic
**Confidence:** High
**Problem:** Mutable bool used purely to communicate between two adjacent `if` blocks that are mutually exclusive.
**Suggested fix:**
```rust
if let Some(last) = self.res.last_mut()
    && candidate.start <= last.1
{
    last.1 = last.1.max(candidate.finish);
} else {
    self.res.push((candidate.start, candidate.finish));
}
```

#### Mi11: src/var_calling/dust_filter.rs:244-250 — `while let` + `if/else break` collapses to `while … is_some_and(…)`

**Category:** idiomatic
**Confidence:** High
**Suggested fix:**
```rust
while self.perf.last().is_some_and(|p| p.start < start_threshold) {
    self.perf.pop();
}
```

#### Mi12: src/var_calling/dust_filter.rs:220-229 — Manual `len() + index` where `Vec::last` reads more directly

**Category:** idiomatic
**Confidence:** High
**Suggested fix:**
```rust
let Some(&candidate) = self.perf.last() else { return; };
if candidate.start >= start_threshold {
    return;
}
```

#### Mi13: src/var_calling/dust_filter.rs:257-316 — `i64` reverse cursor in `find_perfect` should be `(0..hi).rev()` over `u32`

**Categories:** idiomatic, smells (convergent; naming and refactor_safety also flagged as cross-cat)
**Confidence:** Medium
**Problem:** The downward sweep uses an `i64` cursor purely so the loop guard can compare against zero. Every body read casts back to `u32` or `usize`. The `if self.big_l + 1 > win_len { return; }` early-guard already handles the empty-range case, so the negative-index termination is not load-bearing. A `(0..prefix_len).rev()` loop keeps the index in `u32`, removes both casts inside the body, and is structurally identical to the C `int i` downward iteration (Rust just doesn't need a sentinel).
**Suggested fix:**
```rust
let win_len = self.window.len() as u32;
let prefix_len = win_len.saturating_sub(self.big_l);
let mut c = self.cv;
let mut r = self.rv;
let mut max_r: u32 = 0;
let mut max_l: u32 = 0;
for i in (0..prefix_len).rev() {
    let t = self.window[i as usize] as usize;
    r += c[t];
    c[t] += 1;
    let new_r = r;
    let new_l = win_len - i - 1;
    // ... rest of body, `i` already u32, casts removed ...
}
```

#### Mi14: src/var_calling/dust_filter.rs:499-519 — Co-dependent fields `current_chrom` / `current_mask` / `sweep`

**Category:** idiomatic
**Confidence:** Medium
**Problem:** Three private fields whose meaningful truth table has fewer rows than the cross product. The invariant is upheld only by convention in `ensure_mask_for`.
**Suggested fix:**
```rust
struct LoadedChrom {
    chrom_id: u32,
    mask: SdustIntervals,
    sweep: usize,
}

pub struct DustFilter<I, F> /* ... */ {
    // ...
    loaded: Option<LoadedChrom>,
    is_finished: bool,
}
```
`ensure_mask_for` assigns `self.loaded = Some(LoadedChrom { … })` atomically; same byte footprint, one fewer way to leave the iterator inconsistent.

#### Mi15: src/var_calling/dust_filter.rs:1083 — Test match pattern uses bare `..` on `RefFetch`

**Category:** refactor_safety
**Confidence:** High
**Problem:** `Some(Err(DustFilterError::RefFetch { chrom_id: 0, .. })) => {}` absorbs any future field added to the variant. The outer enum is `#[non_exhaustive]` but struct-style-variant inner fields are not.
**Suggested fix:**
```rust
Some(Err(DustFilterError::RefFetch { chrom_id: 0, source: _ })) => {}
```

#### Mi16: src/var_calling/dust_filter.rs:106 — Public `SdustIntervals = Vec<(u32, u32)>` is primitive obsession

**Category:** smells
**Confidence:** Medium
**Problem:** Tuple-as-interval semantics conveyed by convention. Every consumer must remember `.0 = start`, `.1 = finish`. A misread silently produces a coordinate bug. The naming is also worse than it could be: `last.1`, `current_mask[sweep].1`, etc.
**Suggested fix:** Introduce a `BedInterval { start: u32, end: u32 }` newtype with `contains(pos)`, update `save_masked_regions`/`is_masked`/`find_perfect`'s insertion site and the golden-snippet literals. A small `bed(s, e)` test helper absorbs the test-side noise.

#### Mi17: src/var_calling/dust_filter.rs:200,280,370 — `10·X > T·Y` density predicate open-coded in three places

**Category:** smells
**Confidence:** High
**Problem:** The score-density predicate with `u64` widening recurs at three sites, each with its own widening boilerplate. Smells rule "Duplicated logic" triggers at 3+.
**Suggested fix:**
```rust
/// `true` iff `10·score > threshold · length`, with u64 widening
/// to avoid overflow at large `threshold`. Mirrors the score-density
/// gate from `sdust_core` (sdust/sdust.c:78, 111, 148).
#[inline]
fn dust_density_exceeds(score: u32, length: u32, threshold: u32) -> bool {
    (score as u64) * 10 > (threshold as u64) * (length as u64)
}
```

#### Mi18: src/var_calling/dust_filter.rs:559-573 — `self.done = true` mutated inside `ok_or_else` / `map_err` closures

**Categories:** errors, naming, smells, idiomatic (convergent across four agents)
**Confidence:** High
**Problem:** Latching the iterator is part of the documented contract, but the latch-set assignments are hidden inside closure bodies. Functional but easy to drop on a refactor that replaces `ok_or_else` with `?` and a custom error constructor.
**Suggested fix:**
```rust
let entry = match self.chromosomes.get(chrom_id as usize) {
    Some(e) => e,
    None => {
        self.done = true;
        return Err(DustFilterError::UnknownChromId { chrom_id });
    }
};
let seq = match self.fetcher.fetch(chrom_id, 1, entry.length) {
    Ok(s) => s,
    Err(source) => {
        self.done = true;
        return Err(DustFilterError::RefFetch { chrom_id, source });
    }
};
```

#### Mi19: src/var_calling/dust_filter.rs:257 — `find_perfect`'s `big_l <= window.len()` invariant asserted only by comment

**Category:** reliability
**Confidence:** Medium
**Problem:** The trim invariant is real but documented only in prose. A future refactor of `shift_window` that broke it would silently produce wrong (empty-ish) output via the early return at line 268-270.
**Suggested fix:**
```rust
fn find_perfect(&mut self, start: u32) {
    let win_len = self.window.len() as u32;
    debug_assert!(self.big_l <= win_len,
        "shift_window must keep big_l <= window.len(); big_l={} win_len={}",
        self.big_l, win_len);
    if self.big_l + 1 > win_len { return; }
    // ...
}
```

#### Mi20: src/var_calling/dust_filter.rs:843 — Golden-vector test doesn't record `lh3/sdust` commit SHA

**Category:** reliability
**Confidence:** Medium
**Problem:** The comment block at lines 788-802 explains the `expected` values were generated from `lh3/sdust` at impl time but does not pin the upstream commit. If a future bug-fix to the C source changes its output, our golden vectors silently drift from "what `lh3/sdust` produces."
**Suggested fix:** Documentary only — add the upstream commit SHA to the comment block:
```rust
// generated once at implementation time (2026-05-17) by running
// `lh3/sdust -w 64 -t 20` at commit <SHA> against each snippet's
// bases and copying its BED-style output verbatim.
```

#### Mi21: tooling — `cargo doc --no-deps` not run; PR enables `broken_intra_doc_links = "deny"` with 10+ intra-doc links

**Category:** tooling
**Confidence:** High
**Problem:** `Cargo.toml` declares `[lints.rustdoc] broken_intra_doc_links = "deny"`, but the verification recap shows `cargo doc --no-deps: not run`. The new module ships with intra-doc links to `DustFilter`, `DustFilter::next`, `sdust_mask`, `DustFilterConfig::new`, `DustFilterConfig::default`, `PileupRecord::pos`, `PerPositionPileups`, `PerPositionMergerError`, `RefSeqFetcher`, `Self::new`. The lint only fires when rustdoc runs.
**Suggested fix:** Run `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features` inside the container; if it passes, note that in the fixes-applied report. If any link fails to resolve, qualify it the same way `PileupRecord::pos` already is on line 46.

#### Mi22: src/var_calling/dust_filter.rs:788-802 — Module-level doc and spec disagree on memory model

**Category:** extras
**Confidence:** High
**Problem:** `calling_pipeline_architecture.md` (spec, line 981) says "Memory is bounded by the window: the scorer maintains O(w) state … does not need a whole-genome mask in memory." The implementation loads the whole chromosome's reference plus the full per-chromosome masked-interval list. The plan and impl report document this trade openly; the module-level doc does not acknowledge the spec departure.
**Suggested fix:** Either update the spec to reflect the per-chromosome-batch design, or add a one-line breadcrumb to the module doc:
> This trades the spec's stated O(w) memory bound (calling_pipeline_architecture.md §"Stage 3") for simplicity and a direct one-to-one port of `sdust_core`; see [`implementation_plans/dust_filter.md`](../../doc/devel/implementation_plans/dust_filter.md) §"Streaming model" for the rationale.

#### Mi23: src/var_calling/dust_filter.rs:347 — `i as u32 + 1 - l` can silently wrap on `seq.len() >= 2^32`

**Category:** extras
**Confidence:** Medium
**Problem:** `sdust_mask` is `pub` and accepts arbitrary `&[u8]`. Inside the loop, `i as u32 + 1` wraps if `i >= 2^32 - 1`. No realistic eukaryotic chromosome is that large (human chr1 is ~2.5×10^8), but the "untrusted input" framing in the review scope keeps this in play. The output type `Vec<(u32, u32)>` cannot represent honest intervals beyond `2^32 - 1` regardless.
**Suggested fix:** Either return `Result<SdustIntervals, DustFilterError>` from `sdust_mask` with a `ReferenceTooLong` variant, or assert at the top of `sdust_mask`:
```rust
let len_u32 = u32::try_from(seq.len())
    .expect("sdust_mask: seq longer than u32::MAX is unsupported");
```
The non-panicking variant is preferred under graceful-degradation framing.

#### Mi24: src/var_calling/dust_filter.rs:333 — Dual validation paths for `window`: `sdust_mask` asserts, `DustFilterConfig::new` returns `Result`

**Categories:** errors, smells (convergent cross-cat)
**Confidence:** High
**Problem:** Two enforcement paths for the same invariant with subtly different ranges (`sdust_mask` asserts `>= 3`; `DustFilterConfig::new` validates `3..=4096`). If a caller bypasses `DustFilterConfig` and calls `sdust_mask` directly with `window = 2`, the panic is hit; if they pass `4097`, the panic is *not* hit but no validation has happened either. See **Open question 3**.
**Suggested fix:** Either (a) keep `sdust_mask` as a "validated precondition" function and document the precondition by reference to `DustFilterConfig::new` (the `# Panics` section already names the lower bound; extend it to mention the project's convention is to go through `DustFilterConfig::new`); or (b) make `sdust_mask` take `DustFilterConfig` instead of two raw `u32`s, eliminating the dual path entirely.

### Nits

Grouped, not numbered:

- **dust_filter.rs:271, 285, 358, 378** — `as u32` / `as i64` casts in the hot loop without `// safe:` comments where the safety reasoning is non-trivial. A one-word `// safe: big_l < win_len` next to each would save a reviewer the re-derivation.
- **dust_filter.rs:182** — `expect("deque non-empty at cap")` message slightly imprecise — the deque is at *or above* capacity when this fires. `"deque non-empty above eviction threshold"` more accurate.
- **dust_filter.rs:392** — Comment "we match that behaviour exactly so the golden-vector test against the cloned binary holds" — after the rework, the test is purely against the committed const; tighten to "so the committed golden vectors hold."
- **dust_filter.rs:499** — Consider `#[must_use = "DustFilter is lazy; iterate it to filter positions"]` on the struct.
- **dust_filter.rs:730-774** — Hoist `use std::iter::repeat_n;` to the `mod tests` use-block instead of fully qualifying four times.
- **dust_filter.rs:118** — `PerfInterval::r` and `PerfInterval::l` single-letter field names are documented at the struct level but read as opaque at use sites (e.g. `p.r as u64 * max_l as u64`). Acceptable per the project's stated C-traceability policy; the cost is real.
- **dust_filter.rs:144** — `PerfInterval`/`perf` abbreviation acceptable per C-traceability; `CandidateInterval`/`candidates` would read more naturally if breaking C-name parity is ever on the table.
- **dust_filter.rs:200,280,291,299,370** — Literal `10` (the sdust score-density multiplier) recurs in five places; a named `SDUST_DENSITY_FACTOR: u32 = 10` would tighten the link, at the cost of one layer of indirection.
- **dust_filter.rs:906** — `StubFetcher::fail_for` reads as a predicate name with no predicate noun. `fail_for_chrom_id` would be explicit. Private test helper, low cost to leave.

## 7. Out of scope observations

- **`DustFilter::new` takes `chromosomes: Vec<ParsedChromosome>` by value with no type-level link to the upstream merger.** The caller is expected to pass `merger.chromosomes().to_vec()`; the type system does not enforce this is the same table the merger holds. API-design follow-up — would slot into the cohort CLI wiring PR. Convergent: errors (cross-cat) and smells (cross-cat).
- **`self.chromosomes.get(chrom_id as usize)` assumes `chrom_id` is a positional index.** Upstream-contract follow-up. If `chrom_id` is an opaque identifier (which the merger's `PerPositionMergerError::OutOfOrder` text implies it might be — it carries a raw u32), the lookup is semantically wrong. Worth pinning the assumption in a comment or via a `HashMap<u32, ParsedChromosome>`.
- **`ensure_mask_for` loads the whole chromosome's bases into memory.** ~250 MB transient for human chr1. This is a documented design tradeoff in the plan ("Streaming model" §) and the impl report ("Tradeoffs and follow-ups" §). Not a new finding; tracked in the plan as a follow-up to revisit if real cohort runs show memory pressure.
- **`find_perfect`'s `self.perf.insert(j, …)` is `O(W²)` per window in the worst case.** The C source `lh3/sdust` has the same shape (`memmove` inside a `for` loop at `sdust.c:120`). Intentional port fidelity; performance follow-up if profiling shows it matters. Worth keeping a comment near the insertion site noting that `Vec::insert` is `O(n)` and the chosen structure mirrors the C reference.
- **Fetcher length validation.** `ensure_mask_for` does not verify the fetcher returned `entry.length` bytes — `sdust_mask` runs on whatever the fetcher hands back, and the resulting interval coordinates are in those bytes' index space. Probably a property of `ChromBoundaryRefFetcher` to enforce upstream, but worth noting at the boundary.

## 8. Missing tests to add now

Per the reliability checklist's challenge-tests pass:

### `sdust_mask`

- **`sdust_mask_invariants_hold_on_random_inputs`** (see **M2** above) — property-based / seeded-random invariant test covering sorted, non-overlapping, in-bounds, and "no masked base is non-ACGT."
- **`sdust_mask_at_minimum_window_does_not_panic`** — input class: `window = SD_WLEN = 3`. Catches off-by-one in `triplet_cap = window_cap - SD_WLEN + 1` at the `SdustState::new` construction site (line 158).
- **`sdust_threshold_zero_masks_aggressively`** — input class: `threshold = 0`. Catches any future short-circuit that treats `T=0` as "disable."
- **`sdust_threshold_max_does_not_overflow`** — input class: `threshold = u32::MAX`. Catches removal of the `u64` widening at the density gates.
- **`sdust_returns_empty_on_seq_shorter_than_wlen`** — input class: `seq.len() < 3` (i.e. `b""`, `b"A"`, `b"AC"`). Catches panics or wrong output on very short references.
- **`sdust_handles_seq_of_exactly_wlen_bases`** — input class: `seq.len() == 3`. Catches the loop-body never firing.
- **`sdust_returns_empty_on_all_n`** — input class: all-N input. Catches panic / wrong output on long N runs.
- **`sdust_handles_long_n_run_after_low_complexity`** — input class: low-complexity prefix then long N gap. Pins the flush-perf-on-N branch with a non-empty `perf`.
- **`sdust_intervals_finish_le_seq_len`** — structural assertion across edge cases. Catches finish > seq.len().

### `DustFilter`

- **`filter_half_open_boundary_pos_at_start_masked_pos_at_finish_passes`** (see **M3**) — pins half-open semantics at `pos = start` and `pos = finish`.
- **`filter_is_undefined_on_out_of_order_within_chrom`** (see **M4**) — pins current behavior or adds the typed-error fix.
- **`filter_surfaces_invalid_pos_zero_and_latches`** (see **M1**) — `pos == 0` from upstream surfaces `InvalidPos` and latches.
- **`filter_latches_after_upstream_exhaustion`** — second `next()` call after upstream exhaustion stays `None`.
- **`filter_reloads_mask_when_chrom_revisited`** — pin the current "reload on every chrom_id change" semantic so a future merger contract change can't make this silently wasteful.
- **`filter_handles_zero_length_chromosome`** — `ParsedChromosome { length: 0, .. }`.

### `DustFilterConfig::new`

- **`config_new_accepts_threshold_zero_and_max`** — the validation intentionally accepts any `u32`; no test pins that.

## 9. What's good

1. **No `..Default::default()` anywhere** in production or tests; every struct literal names every field — adding a field to `DustFilterConfig`, `SdustState`, or any test fixture will fail compilation at every call site. Mirrors the contamination_estimation review's Mi13 direction.
2. **Two-layer split** (`sdust_mask` pure function + `DustFilter` iterator adaptor) — the algorithmic core is testable in isolation against the cloned binary; the iterator wrapper is just plumbing. Easy to reason about each layer separately.
3. **`u64` widening on every density comparison** (lines 200, 280, 290-291, 298-299, 370) preempts the integer-overflow class of bugs the unchecked decrements (**M8**) still leave open elsewhere.
4. **`#[non_exhaustive]` on `DustFilterError`** with operation-named variants and no foreign types in public variants — matches the project convention and leaves room for future variants without breaking match sites that use catch-all arms.
5. **Module-level doc comment** lays out the algorithm, streaming model, and coordinate-conversion rule explicitly before the code, in plain English — readers do not have to reverse-engineer the contract from the implementation.

## 10. Commands to re-verify

Re-run after applying fixes:

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --all-targets --all-features`
- `./scripts/dev.sh env RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features` *(introduced by this review; see Mi21)*
- New test invocations for the additions listed in §8: `cargo test --lib var_calling::dust_filter`.

### Author response convention

Address each finding by its identifier (M1 … Mi24) with `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Open questions 1-4 in section 4 are upstream of several findings — resolve those first.
