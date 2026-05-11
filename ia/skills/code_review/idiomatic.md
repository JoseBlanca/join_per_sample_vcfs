# Idiomatic Rust checklist

**Purpose.** Code reads as Rust the language and ecosystem expect — iterators, exhaustive matches, type-level invariants, minimal lifetimes, no needless allocation.

**Triggers.** Read every match expression, every function signature, every struct/enum definition, every use of `clone()` / `to_string()`, every `Box`/`Vec`/`String` parameter, every `impl Copy`, every visibility modifier.

**Skip when.** Never skipped.

## Rules

- **Iterators over index loops.** `?` over `match` for propagation. `if let` / `let ... else` over nested matches.
- **Pattern matches on internal enums are exhaustive** — no `_` arm that would silently absorb a new variant. For `#[non_exhaustive]` external enums or primitives where exhaustive matching is impossible, the catch-all includes a `// REVIEW ON UPGRADE:` comment. (Pairs with `clippy::wildcard_enum_match_arm = "deny"`.)
- **Make illegal states unrepresentable.** Triggers a finding when:
  - (a) a function takes a primitive with domain meaning (use a newtype: `Cents`, `UserId`);
  - (b) two `bool` parameters at the same call site (use an enum);
  - (c) a struct's fields are only valid in certain combinations (split into an enum);
  - (d) a return type allows an "empty" value where empty is invalid (use `NonEmpty<T>` or a constructor-validated newtype);
  - (e) a numeric field has a documented range (use a newtype with a fallible constructor).
- **Control construction at the type level.** When a type's invariants are established by a constructor, prevent construction by other paths rather than relying on convention:
  - **Private unit field** (`_private: ()`) on a struct with otherwise-`pub` fields blocks struct-literal construction outside the defining module while keeping field reads public. Use when a `pub` struct must be constructed only via a validating constructor.
  - **Sealed traits** for `pub trait`s that must not be implementable downstream: nest a private supertrait in a `mod private { pub trait Sealed {} }` and require it as a bound. Use when adding a method must remain non-breaking, or when correctness depends on knowing every implementor.
  - `#[non_exhaustive]` is the cross-crate equivalent for structs and enums.
- **`#[must_use]`** is required on builders, status-check functions, RAII guards, transaction handles, and any type whose construction implies an obligation. (Pairs with `clippy::must_use_candidate = "warn"`.)
- **Visibility is minimized.** Private by default; `pub(crate)` for cross-module use; `pub` only for the external API. A `pub` item with no caller outside its module, or with no `///` doc comment, is a finding.
- **Minimize mutable scope.** Two forms reduce the surface where partially-initialized or stale state is observable:
  - **Temporary mutability via shadowing**: `let mut buf = …; populate(&mut buf); let buf = buf;` rebinds the name to immutable so the rest of the function cannot mutate. Trigger a finding when a `mut` binding is mutated only during initialization but the rest of the function still sees `&mut`.
  - **Scope-based initialization**: `let buf = { let mut v = Vec::new(); …; v };` confines `mut` to the inner scope and exposes only the finished value to surrounding code.
- **Lifetimes are correct and minimal.** No `'static` added to silence the borrow checker (use `Arc`, owned data, or a different ownership boundary); no elidable annotations written explicitly; no lifetime parameter appearing only once in a signature.
- **Trading `'a` for `Arc` is a deliberate choice, not a default.** Replacing a lifetime parameter with `Arc<T>` is the right move when ownership is genuinely shared (multiple threads, multiple long-lived owners, graph-shaped data) or when threading `'a` through three or more types is making the API viral. It is the wrong move when (a) the data has a single clear owner and a borrow suffices — refcount bumps and heap allocation have a real runtime cost; (b) the inner type still borrows (`Arc<Cow<'a, str>>`, `Arc<&'a T>`) — the `'a` is still there and the `Arc` only added indirection; (c) the change is cosmetic, hiding an unclear ownership tree rather than resolving it. When `Arc` *is* the right call and the data is immutable after construction, prefer `Arc<str>` over `Arc<String>` and `Arc<[T]>` over `Arc<Vec<T>>` — same ownership semantics, one fewer indirection.
- **Function parameters use the most general bound the implementation requires.** `impl AsRef<Path>`, `impl IntoIterator<Item = T>`, `impl Into<String>`. Over-specific signatures force needless conversions on callers.
- **No needless allocation.** `.clone()` where a borrow suffices, `.to_string()` / `.to_owned()` where `&str` works, `String`/`Vec<T>` parameters where `&str`/`&[T]` work, `format!` where `Display` or `write!` would do, `.collect::<Vec<_>>()` followed by another iteration, `Box<T>` where `T` fits by value.
- **`Copy` discipline.** Implemented only on small (≤ ~16 bytes), value-like types where bitwise copy is semantically meaningful. Each `.clone()` at a call site is justified by ownership or replaced with a borrow.
- **Imports are grouped** (`std`, external, internal) and individual `use` statements are preferred over globs (except preludes and test modules). Re-exports are intentional.
