# Error handling — errors never pass silently

**Purpose.** No silent failures, no swallowed `Result`s, no panics in production paths, every error carries enough context to be debugged from a log line alone.

**Triggers.** Search the in-scope code for: `unwrap`, `expect`, `panic!`, `unreachable!`, `todo!`, `unimplemented!`, `let _ =`, `.ok()` (without recovery), `if let Ok(_)`, `Box<dyn Error>`, plain integer arithmetic on untrusted input, `Drop` impls, every `?`. Read every error type and `From`/`Into` impl.

**Skip when.** Never skipped.

## Rules

- **No `unwrap()` or `expect()`** in non-test code unless paired with a `// PANIC-FREE:` comment naming the invariant that guarantees the call cannot fire. Prefer restructuring to eliminate the call entirely. Each violation without the comment is at least Major; on a Blocker-class path (data loss, security), upgrade.
- **No `panic!()` or `unreachable!()`** in non-test code unless the branch is genuinely unreachable given the type system, with a `// UNREACHABLE:` comment explaining why. Prefer making the state unrepresentable.
- **No `todo!()` or `unimplemented!()`** in code submitted for review.
- **No discarded `Result` or `#[must_use]` value** (`let _ = …`, `.ok()` without recovery, `if let Ok(_) = …` ignoring the error) without a comment naming the specific failure mode being discarded and why it is safe to ignore. Silent discard of a meaningful `Result` is Blocker.
- **Errors are typed** (`thiserror` for libs, `anyhow` for bins). No `Box<dyn Error>` in public library APIs. Error enums are scoped to the module or operation that produces them, not crate-wide. Catch-all `Other(String)` variants require justification.
- **Errors carry context.** Each `?` across a meaningful boundary adds context identifying *what was being attempted* and *which input caused it*. Root causes are preserved via `#[source]` or `#[from]`, never flattened into a string.
- **`From`/`Into` conversions are total and lossless.** Lossy or fallible conversions use `TryFrom`/`TryInto`. (Pairs with `clippy::fallible_impl_from = "deny"`.)
- **Integer arithmetic on values derived from untrusted input** uses `checked_*`, `saturating_*`, or `wrapping_*` explicitly. Plain operators are only exceptionally acceptable on values whose range is provably bounded by the type system; in release builds they wrap silently — that is data corruption, not a panic.
- **Recovery points emit logs.** Every retry, fallback, default branch, or caught-and-ignored error emits a `tracing` event at `warn` or higher with structured fields (error chain, causing input, recovery action). Silent recovery is silent failure.
- **`Drop` impls do not panic.** Failures during cleanup are logged via `tracing` and swallowed.
