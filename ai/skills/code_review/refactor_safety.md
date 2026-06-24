# Refactor safety — let the compiler flag the change

**Purpose.** When a struct gains a field, an enum gains a variant, or a slice's expected length changes, the compiler should flag every site that assumed the old shape. Code that still compiles but silently does the wrong thing after a refactor is the dominant defensive-Rust failure mode and is almost always avoidable at the language level.

**Triggers.** Read every struct literal, every match on a slice, every `Default` use in initialization, every manual `PartialEq`/`Eq`/`Hash`/`Clone`/`Debug` impl, every partial destructure with `..`.

**Skip when.** Never skipped.

## Rules

- **Match on slices, not length-check then index.** `match s { [a, b] => …, [a, b, c] => …, _ => … }` instead of `if s.len() == 2 { let a = s[0]; let b = s[1]; … }`. Each expected length adds a compile error at every site that assumed a different one; index-based code panics at runtime. Trigger a finding when reachable code indexes a slice with a literal after an `if .len() == N` check, or any time a fixed-shape destructure would replace a length-conditional index. (Pairs with `clippy::indexing_slicing = "deny"`.)
- **No `..Default::default()` in struct literals.** Spell every field explicitly. When a new field is added the compiler flags every literal that needs to think about it; `..Default::default()` silently fills it with whatever the `Default` impl picks, which is rarely what the call site intended for the new field. The destructuring-shorthand alternative — `let StructName { a, b, c, .. } = StructName::default()` and rebuild explicitly — is acceptable when many fields would be repeated, *provided* the destructure is itself exhaustive (no `..` on the destructure pattern).
- **Destructure all fields in trait impls that depend on the field set** — `PartialEq`, `Eq`, `Hash`, manual `Clone`, `Debug` impls that pick which fields to print. `let Self { a, b, c } = self;` (no `..`) forces a new field to surface in the impl. A `..` catch-all when destructuring an internal type in such an impl is a finding.
- **Named field patterns over bare `_`.** `Foo { has_fuel: _, .. }` rather than relying on `..` alone when ignoring a field deliberately. The named form documents *which* field is being ignored and forces a compile error if the field is renamed; `..` silently absorbs structural change.
