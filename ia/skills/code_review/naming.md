# Naming & readability checklist

**Purpose.** Names that communicate domain meaning, consistent across the crate, with constants for non-obvious literals.

**Triggers.** Read every function signature, type definition, field name, variable name, and module name in scope. Search for numeric and string literals.

**Skip when.** Never skipped.

## Rules

- **Function names are verbs or verb phrases** (`parse_header`, `compute_checksum`). Pure accessors returning a field-like value are an exception and use the noun (`len`, `name`, `as_bytes`) per Rust convention. Avoid placeholders (`header`, `data2`, `do_stuff`).
- **Types, fields, and variables name what they *are* in the domain**, not their generic shape. Forbidden when used alone: `data`, `info`, `manager`, `handler`, `helper`, `util`/`utils`, `tmp`, `foo`/`bar`, `thing`, `object`, `item`, `value`, `result`. Acceptable as suffixes when the prefix carries the meaning (`parse_result`).
- **Names describe domain meaning, not structural shape.** Don't encode the number of fields, their primitive types, or their container in a type or field name (`FiveScalars`, `IntPair`, `ThreeBytes`, `U64Map`, `ConfigStruct`). Such names go stale the instant the shape changes — a sixth field makes `FiveScalars` a lie — and force the reader to read every field before guessing what the type *represents*. Pick the domain term (`AlleleEvidence`, `Coordinate`, `ByteRange`). Exception: when the shape *is* the meaning and is fixed by an external contract (`Sha256Hash`, `Rgb`, `Ipv4Addr`).
- **Acronyms follow Rust convention.** `HttpClient`, not `HTTPClient`; `parse_url`, not `parse_URL`.
- **No abbreviations** except (a) Rust standard conventions (`len`, `ctx`, `buf`, `cfg`), (b) acronyms universally known in the crate's domain (`HTTP`, `AST`), or (c) variables whose lifetime is bounded to ~3 lines (`i` in a loop). Domain-specific abbreviations require a glossary entry in the crate-level docs.
- **Boolean variables, fields, and `-> bool` functions read as predicates** (`is_valid`, `has_children`, `contains_key`, `should_retry`). No noun-form booleans (`validity`, `children`, `retry`).
- **The same concept uses the same name everywhere in the crate.** `user_id` in one signature is `user_id` in all signatures — not `uid`, `user`, or `id`.
- **Module names reflect the concept they own** (`parser`, `connection`, `tariff`), not the layer or pattern (`models`, `services`, `helpers`, `core`, `common`). A reader scanning `src/` should be able to predict each module's contents from its name alone.
- **Names are as long as they need to be to be unambiguous in their scope, and no longer.**
- **No magic numbers.** Numeric and non-trivial string literals whose meaning is not obvious from surrounding code live in a named `const` with a short doc comment naming units and source (spec section, measurement, vendor recommendation). `pub const` for cross-module / public-API defaults, private `const` for module-local intent, function-scoped `const` for purely local values. Triggers a finding when:
  - the literal would need to be kept in sync with anything else (another code site, a spec, a config file, an external tool's default);
  - a reader has to look beyond the immediate function to know what the value represents;
  - the same literal recurs in three or more places without a shared name.

  Trivial literals are exempt: `0`, `1`, `±1` in arithmetic and indexing, and fixed buffer sizes whose only meaning is the literal itself (`[0u8; 16]` for an MD5). The rule is about readability, not numerology — the question is "would a reader meeting this literal know what it means?", not "is this a number?".
