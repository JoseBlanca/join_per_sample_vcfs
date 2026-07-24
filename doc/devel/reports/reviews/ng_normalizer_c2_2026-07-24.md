# Code Review: ng normalizer — C2 (algorithm 1c, `FixpointLeftAligner`)

**Date:** 2026-07-24
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** working-tree — the 1c additions to `src/ng/alignment/left_align_structured.rs`
**Status:** Approve-with-changes

---

### 1. Scope
- Algorithm 1c: `FixpointLeftAligner` + the generic `drive_to_fixpoint` loop + `FIXPOINT_MAX_ITERATIONS` (plan C2; spec §6; arch §Module home). Lives in `left_align_structured.rs` alongside 1a, per the arch.
- In-scope: the 1c additions only (1a reviewed under B1).
- Categories dispatched (7, parallel): reliability, errors, naming, idiomatic, smells, refactor_safety, extras.

### 2. Verdict
Approve-with-changes.

### 3. Execution status
- `cargo fmt --check` → 0, clean · `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean · `cargo test --lib` → 0, 2328 passed / 4 ignored (post-fix); `left_align_structured` 18 → 19.

### 4. Open questions and assumptions
1. **Fixpoint/panic logic confirmed correct** — already-leftmost converges in 1 iteration; a k-shift in exactly 2 (1a is a one-pass idempotent fixpoint); a non-converging inner exhausts the cap and panics.
2. **Both `#[should_panic]` tests confirmed genuine** — cannot pass for the wrong reason: the stub never slices (no OOB path), the real-1a cap-of-one case is in-bounds and read-consistent (no `left_align_indels` panic), and the `expected` substring matches only the fixpoint panic (not the `debug_assert`/bounds-check messages).
3. **Intent confirmed** — 1a-to-fixpoint; cap fails loudly (panics) not swallowed; thin wrapper doing no shifting; A2-verified + the deliberately-capped test.

### 5. Findings (all Minor/Nit — no Blocker/Major)

- **Mi1 (naming): `to_fixpoint` borrowed the `to_*` by-value-conversion idiom but mutates in place and panics.** **Applied**: renamed `drive_to_fixpoint`.
- **Mi2 (reliability): `FIXPOINT_MAX_ITERATIONS` had no regression guard pinning it above the minimum.** **Applied**: a compile-time `const _: () = assert!(FIXPOINT_MAX_ITERATIONS >= 2, …)`.
- **Mi3 (reliability): no already-leftmost 1c test (the one-iteration path).** **Applied**: `an_already_leftmost_input_reaches_the_fixpoint_without_shifting`.
- **Mi4 (errors): panic message omitted the alignment/read/reference (not reproducible from a log).** **Applied**: appended them, matching the property oracle's panic.
- **Mi5 (smells): the agreement test used an ad-hoc positional `(u64, Vec, &[u8], &[u8])` tuple (two adjacent `&[u8]` invite a swap).** **Applied**: switched to the file's shared `Case` struct.
- **Nits**: stub `NeverStabilizes` (a verb phrase) → **applied** rename to `NonConvergingNormalizer` (types are nouns). `FIXPOINT_MAX_ITERATIONS` `pub` vs `pub(crate)` → **won't fix** (symmetry with 1b's `pub MAX_PASSES`). Iteration/pass vocabulary drift vs 1b → **won't fix** (a 1c iteration and a 1b pass are genuinely different concepts).

### 6. Out of scope observations
None.

### 7. Missing tests to add now
Applied (18 → 19): the already-leftmost one-iteration path. The two cap-panic paths, the two-iteration convergence, and the 1a-agreement were already present and confirmed genuine.

### 8. What's good
- `drive_to_fixpoint` is generic over the inner normalizer, so the cap-exhaustion panic is driven honestly by a `NonConvergingNormalizer` stub — the real inner never reaches it, and the test proves the panic fires on genuine non-convergence, not a slice accident.
- The whole-`Alignment` fixpoint comparison (derived `PartialEq` over both fields) means a `reference_offset`-only change also counts as movement — refactor-safe.

### 9. Commands to re-verify
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::alignment::left_align_structured`
