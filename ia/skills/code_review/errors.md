# Error handling — errors never pass silently

**Purpose.** No silent failures, no swallowed `Result`s, no panics in production paths, every error carries enough context to be debugged from a log line alone.

**Triggers.** Search the in-scope code for: `unwrap`, `expect`, `panic!`, `unreachable!`, `todo!`, `unimplemented!`, `let _ =`, `.ok()` (without recovery), `if let Ok(_)`, `Box<dyn Error>`, plain integer arithmetic on untrusted input, `Drop` impls, every `?`. Read every error type and `From`/`Into` impl, paying attention to: third-party error types appearing in `pub` variants (`sqlx::`, `reqwest::`, foreign `Error` types), missing `#[non_exhaustive]` on public error types, redundant `Error` suffix or mechanism-named variants (`Io`, `Parse`, `Other`) instead of operation-named ones.

**Skip when.** Never skipped.

## Rules

- **No `unwrap()` or `expect()`** in non-test code unless paired with a `// PANIC-FREE:` comment naming the invariant that guarantees the call cannot fire. Prefer restructuring to eliminate the call entirely. Each violation without the comment is at least Major; on a Blocker-class path (data loss, security), upgrade.
- **No `panic!()` or `unreachable!()`** in non-test code unless the branch is genuinely unreachable given the type system, with a `// UNREACHABLE:` comment explaining why. Prefer making the state unrepresentable.
- **No `todo!()` or `unimplemented!()`** in code submitted for review.
- **No discarded `Result` or `#[must_use]` value** (`let _ = …`, `.ok()` without recovery, `if let Ok(_) = …` ignoring the error) without a comment naming the specific failure mode being discarded and why it is safe to ignore. Silent discard of a meaningful `Result` is Blocker.
- **Errors are typed** (`thiserror` for libs, `anyhow` for bins). Catch-all `Other(String)` variants require justification.
- **One error type per fallible operation, not per crate or module.** Test: if two operations have different failure modes *or* should display different messages, split them. A crate-wide `Error` enum that everything funnels through is wrong — it produces poor error chains, forces unrelated callers to match on irrelevant variants, and prevents extracting modules into separate crates.
- **No third-party error types in public variants.** Public variants must not hold `sqlx::Error`, `reqwest::Error`, `Box<dyn Error>`, or other foreign error types directly — that bakes the dependency into your public API and forces every downstream onto your dependency version. Wrap in a private newtype (`pub struct DbError(sqlx::Error)` with the inner field private) or hide behind an opaque inner repr, as `std::io::Error` does. Internal-only error types may expose foreign types freely.
- **Public error types are `#[non_exhaustive]`.** Both the outer struct (when using the struct-plus-`Kind` pattern below) and any `Kind` enum, so variants and fields can be added without a semver break. Internal-only error types may omit it.
- **Variant names describe the failed operation, not the mechanism.** `ReadConfig`, `ParseHeader`, `OpenIndex` — not `Io`, `Parse`, `Other`. A reader matching on `Io` learns nothing; matching on `ReadConfig` learns *what* failed. No redundant `Error` suffix on variants (the type already says `Error`).
- **Errors carry context.** Each `?` across a meaningful boundary adds context identifying *what was being attempted* and *which input caused it*. With `thiserror`, context lives in fields on the error (struct fields or variant fields), and is attached at the `?` site with `.map_err`, e.g. `.map_err(|source| ReadHeaderError { path: path.clone(), kind: source.into() })?` for the struct-plus-`Kind` shape, or `.map_err(|source| Error::ReadHeader { path: path.clone(), source })?` for a flat enum. `#[source]` marks the cause field so it appears in the error chain; `#[from]` only generates a `From` impl and must not be used when the inner type could originate from multiple sites (see the `From` rule). Root causes are preserved through `#[source]`/`#[from]`, never flattened into the parent's `Display` string.
- **`Display` messages follow stdlib convention.** Concise, lowercase, no trailing punctuation. Describe *what was attempted* and *which input* (e.g. `"failed to read VCF header from {path}"`); let `source()` carry the cause. Do not concatenate the source's message into the parent's `Display` — that breaks `anyhow`-style chain rendering.
- **`From`/`Into` conversions are total, lossless, and semantically faithful.** Lossy or fallible conversions use `TryFrom`/`TryInto`. Do not impl `From<Inner>` when `Inner` errors could originate from more than one site — that erases provenance and lets every `?` collapse into the same variant regardless of where it fired. Implement `From` only when one origin exists; otherwise convert explicitly at the `?` site with context. (Pairs with `clippy::fallible_impl_from = "deny"`.)
- **Integer arithmetic on values derived from untrusted input** uses `checked_*`, `saturating_*`, or `wrapping_*` explicitly. Plain operators are only exceptionally acceptable on values whose range is provably bounded by the type system; in release builds they wrap silently — that is data corruption, not a panic.
- **Recovery points emit logs.** Every retry, fallback, default branch, or caught-and-ignored error emits a `tracing` event at `warn` or higher with structured fields (error chain, causing input, recovery action). Silent recovery is silent failure.
- **`Drop` impls do not panic.** Failures during cleanup are logged via `tracing` and swallowed.

## Public library error shape

For an operation whose failures carry contextual data (file path, record number, input span), prefer the struct-plus-`Kind` pattern over inlining context into every variant:

```rust
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
#[error("failed to read VCF header from {path}")]
pub struct ReadHeaderError {
    pub path: PathBuf,
    #[source]
    pub kind: ReadHeaderErrorKind,
}

#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum ReadHeaderErrorKind {
    #[error("I/O error")]
    Io(#[source] std::io::Error),
    #[error("invalid header line {line}")]
    InvalidLine { line: usize, #[source] cause: ParseError },
}
```

The struct carries the operation's context once; the `Kind` enum stays small and matchable. This is how `std::io::Error` and most well-designed library errors are shaped, and it makes the existing "errors carry context" rule mechanical to satisfy: the context lives on the struct, not stringified into each variant's message.

Attach the context at the `?` site:

```rust
let header = parse_header(&bytes).map_err(|kind| ReadHeaderError {
    path: path.clone(),
    kind,
})?;
```

If the inner operation returns a foreign error type (e.g. `io::Error`), convert into the `Kind` first — either via a `From<io::Error> for ReadHeaderErrorKind` impl when `io::Error` has exactly one origin in this operation, or by mapping to the specific `Kind` variant explicitly when it doesn't (per the `From` rule).
