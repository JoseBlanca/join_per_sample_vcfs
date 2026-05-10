# Hot loops & codegen checklist

**Purpose.** Surface inner-loop patterns that block autovectorization, force redundant bounds checks, or otherwise leave compiler optimizations on the table.

**Triggers.** Tight numeric loops, byte / bit processing, slice indexing, generic code in a hot path, dynamic dispatch where static would do, regex / parse loops.

**Skip when.** Code is not numerical or per-record processing in performance-relevant call paths. If dispatched anyway, write `No findings.`

## Rules

- **Iterators are zero-cost most of the time.** Idiomatic `iter().map(...).filter(...).sum()` compiles to the same code as a hand-written `for`. Do not file findings rewriting iterator chains as index loops without a measurement showing the chain hurts. Known exception: `Iterator::fold` has historically optimized worse than the equivalent `for` in some shapes — file at Likely only with a benchmark.
- **`copied()` over `cloned()` for `Copy` types.** `iter().copied()` lets LLVM emit register copies; `cloned()` goes through the `Clone` trait machinery, which is the same for `Copy` types but blocks some optimizations.
- **Bounds checks: avoid by structure, not by `unsafe`.**
  - `chunks_exact(N)` over `chunks(N)` removes a per-iteration tail check and unlocks autovectorization.
  - Iterators (no indexing) eliminate bounds checks the compiler cannot prove redundant.
  - A leading `assert!(idx < slice.len())` (or `assert_eq!(slice.len(), N)`) hoists the bounds check out of the loop body — counterintuitively, more `assert!`s can be faster.
  - `unsafe { *slice.get_unchecked(i) }` is the last resort, requires a `// SAFETY:` comment, and must be paired with a benchmark that proves the win.
- **Autovectorization needs the compiler's confidence.** It vectorizes when the loop body is straight-line, the trip count is known or bounded by `chunks_exact`, the data is contiguous, and there are no data-dependent branches. Filed: a tight numeric loop with an `if` inside that could be a `select` (e.g., `if x > 0 { acc + x } else { acc }` → `acc + x.max(0)`); non-`#[inline]` helpers preventing inlining of the loop body. Cross-reference: `methodology.md` on `target-cpu` — autovectorization to AVX2/AVX-512 needs the right target.
- **`#[inline]` is a hint, not a free win.** Apply when (a) the function is in another crate and called from a hot loop (LLVM cannot cross-crate inline without LTO) or (b) it is small and called many times. Mark single-call helpers `#[inline(always)]` only with a benchmark; over-inlining bloats the icache.
- **Cold paths get `#[cold]`.** Error formatters, panic helpers, slow paths inside a fast-path function — `#[cold]` improves branch prediction layout and keeps the cold body out of icache. The `#[cold]` attribute is the closest stable equivalent to `likely` / `unlikely` (the intrinsics are nightly-only).
- **Static dispatch in hot loops.** `fn f<T: Trait>(x: T)` (monomorphizes) over `fn f(x: &dyn Trait)` (vtable indirection) when the type set is closed and the call is in an inner loop. Cost: more codegen, possible icache pressure. File only when a vtable call is on the per-record path.
- **Default `HashMap` hasher is slow for non-adversarial workloads.** `std::collections::HashMap` uses SipHash for HashDoS resistance; for internal maps with trusted keys, `rustc_hash::FxHashMap` (or `ahash::AHashMap`) is faster. Use the default for any map keyed on data from outside the process; switch to FxHash for internal maps where the keys are integers or short owned types and the gain is benchmarked.
- **`HashMap::entry` over double lookup.** `if let Some(v) = m.get(k) { ... } else { m.insert(k, ...) }` does two hashes; `m.entry(k).or_insert_with(...)` does one.
- **Avoid recomputing in tight loops.** Hoist invariants. `Regex::new` inside a loop recompiles per iteration; `s.chars().count()` is `O(n)`; allocations in the closure passed to `map` show up at scale.
- **`format!`, `write!`, and `Display::fmt` are not equivalent.** Inside a hot serializer, write directly to a reused `String` or `Vec<u8>` with `write!(buf, ...)` (and do not `unwrap` the result if the writer is `String` or `Vec<u8>` — both are infallible). Cross-reference: `allocations.md` on `format!` allocating a fresh `String`.
- **Match on bytes for ASCII parsers.** `match b { b'A' => ..., b'C' => ..., ... }` over `match c { 'A' => ..., 'C' => ... }` for ASCII-only protocols (VCF, BED, FASTA): byte matches lower the trailing UTF-8 decoding step out of the loop.
