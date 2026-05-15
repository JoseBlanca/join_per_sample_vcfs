# Defaults must be announced

**Purpose.** A default a caller cannot see at the call site is not announced, no matter how well-documented it is elsewhere. The API itself is the primary surface for default visibility; docs and logs are reinforcements, not substitutes.

**Triggers.** Read every `Default` impl, every `new()`/`builder()` constructor, every parameter typed `Option<T>`, every numeric or string literal acting as a fallback, every config struct, every CLI/env-driven setting.

**Skip when.** The scope contains no public API surface and no configuration values — no defaults to announce. The orchestrator should not dispatch this category for trivial pure-data snippets.

## Rules

- **Visible at the call site.** Functions whose behavior depends on implicit values take an explicit configuration parameter. Never bury a behaviorally-significant default behind a no-argument constructor: `Client::new(url)` hides the timeout; `Client::builder(url).build()` or `Client::new(url, ClientConfig::default())` exposes it.
- **`Default` is discouraged.** A type may implement `Default` only if (a) the defaults are behaviorally inert (zero, empty, identity) or (b) an equivalent explicit constructor exists whose docs enumerate every default chosen. `#[derive(Default)]` on a config struct whose fields drive runtime behavior is a finding.
- **No hidden `Option` fallbacks.** A parameter typed `Option<T>` that silently falls back to a non-`None` default is a hidden default. Either take `T` directly, take a documented config struct, or name the parameter to reveal the fallback and document the default in a `# Defaults` section of the function's `///` docs.
- **Inspectable at runtime.** Every effective default is recoverable from a running instance via a public accessor or structured dump (`Debug`, `Serialize`, startup `tracing` event listing all resolved values). "Read the source" is not an acceptable answer to "what timeout is this using?".
- **Logged when applied, not only at startup.** When a default is *used* — including lazy, per-request, and fallback defaults — emit a `tracing` event at `debug` or higher identifying the field, the value, and why no explicit value was provided.
- **No magic numbers (defaults specifically).** Every behavioral default lives in a named `pub const` with a doc comment explaining the choice, the units, and the source (spec, measurement, vendor recommendation). Public API docs reference the constant by name (`/// Defaults to [`DEFAULT_TIMEOUT`]`) so the doc and the value cannot drift.
- **Builder visibility.** Configuration types use the builder pattern or explicit struct literals. `.build()` methods document every field filled from defaults; the builder's `Debug` impl distinguishes explicitly-set fields from defaulted ones.
- **Environment-dependent defaults.** Document precedence (CLI > env > config file > built-in), test it (one test per layer plus one per pairwise override), and emit a startup `tracing` event recording each effective value and the source it came from. An operator must be able to answer "where did this value come from?" from logs alone.
