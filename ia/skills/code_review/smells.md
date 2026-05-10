# Code smells to flag

**Purpose.** Heuristics that prompt a closer look — each one is a candidate finding, not an automatic one. When a smell triggers, propose a concrete refactor or accept the existing code with a one-line justification.

**Triggers.** Skim every function for length and branching, every struct for cohesion, every parameter list for primitive obsession and boolean parameters, every `Deref` impl on a non-pointer type, every `#[allow(...)]`, every `TODO`/`FIXME`, every commented-out block.

**Skip when.** Never skipped.

## Rules

- **Long or branchy functions.** Longer than ~75 lines, or with more than ~4 distinct decision points (`if`, `match` arm, `?`, loop condition). Use `cargo clippy -- -W clippy::cognitive_complexity` to surface the worst. A long flat state machine may be clearer than its decomposed alternative — justify and move on.
- **Deeply nested control flow** (>3 levels). Extract, invert with early `return`, or use `let ... else`.
- **Duplicated logic** in three or more places. Extract to a function, trait, or macro. Two occurrences may be left with a follow-up note.
- **"God" structs.** Triggers when (a) different methods operate on disjoint subsets of fields, (b) more than ~7 fields with no clear grouping, or (c) the struct requires `Arc<Mutex<Self>>` because too many things mutate too many fields. Split along the field-disjointness boundary.
- **Primitive obsession.** `String` for IDs/emails/paths/URLs, `u64` for money, tuples for structured data. Each occurrence is a finding pointing to a newtype.
- **Boolean parameters at call sites.** `do_thing(true, false)` is unreadable. Replace with an enum or named struct fields. (Pairs with `clippy::fn_params_excessive_bools = "deny"`.)
- **Stringly-typed APIs.** `String`/`&str` parameters or returns representing values from a finite set. Use an enum.
- **Long `Result`/`Option` `match` chains** where `?`, `.map`, `.and_then`, `.ok_or`, or `let ... else` would read better. Conversely, deeply nested `.and_then` chains where sequential `?` would be clearer.
- **Macros where a generic function would work.** `macro_rules!` is justified for syntax functions cannot express, not for avoiding type parameters.
- **Comments explain *why*, not *what*.** The *what* belongs in code — better names, smaller functions, self-documenting types. Inverse smell: missing `///` docs on `pub` items, or invariants that live only in the author's head.
- **Dead code, commented-out code, or `#[allow(...)]` without justification.** `#[allow(...)]` requires a comment naming (a) the specific lint, (b) why suppression is correct here, (c) the condition for removal. `#[allow(dead_code)]` on genuinely unused items is a finding — delete them.
- **`TODO`/`FIXME`/`XXX`** without an issue reference, owner, and removal condition.
- **`Deref` polymorphism.** Implementing `Deref<Target = Inner>` on a wrapper to "inherit" `Inner`'s methods confuses method resolution, breaks `&Wrapper` / `&Inner` substitutability in generic code, and silently re-exports internal API on every change to `Inner`. `Deref` is for genuine smart pointers (`Box`, `Rc`, `Arc`, `MutexGuard`) where the target *is* the conceptually pointed-to data, not for wrapper structs that happen to contain another type. Fix: composition with explicit forwarding methods, or a newtype with the inner value exposed via a named accessor.
- **`#![deny(warnings)]` in source.** Couples the crate's build to whatever the toolchain decides to lint next; a routine compiler upgrade can break the build for reasons unrelated to the code. Enforce `-D warnings` in CI instead, where the toolchain is pinned.
