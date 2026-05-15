---
name: rust-clone-audit
description: Use this skill when auditing Rust code for unnecessary `.clone()`, `.to_string()`, `.to_owned()`, `.to_vec()`, `Arc::clone`, `Rc::clone`, `String::from(&...)`, `format!("{x}")` followed by re-borrow, or `iter().map(|x| x.clone()).collect()` calls. Trigger on phrases like "audit my clones", "are these clones necessary", "review for needless allocations", "spot the cheap copies", "is this clone justified", or when reviewing Rust code that introduces new clones in a hot path. Pair with the broader `rust-code-review` skill when the audit is part of a full review; use standalone when the user is specifically asking about clones or allocation patterns.
---

# Rust Clone Audit

You are auditing Rust code for unnecessary clones — copies that allocate when a borrow would do, or that pay deep-copy cost when a refcount bump or no copy at all suffices. This skill differs from clippy's lints by **weighing complexity-of-fix against likely speedup** and explicitly recommending "leave alone" or "benchmark first" for borderline cases. Surfacing every clone is noise; surfacing the wrong ones drives over-engineering.

## Core principles

- **Measure cost in two axes:** (a) likely speedup if removed, (b) complexity of the fix. A high speedup × cheap fix is a slam dunk; a small speedup × big fix is not worth proposing.
- **Default to silence on Tier D.** A clone that cannot plausibly affect performance is not a finding — reporting it adds noise and trains the reader to ignore future findings.
- **Recommend benchmark, don't simulate it.** When the fix's value is uncertain, instruct the author to run criterion. Never invent or estimate perf numbers.
- **`Arc::clone` / `Rc::clone` are refcount bumps, not deep copies.** The skill must explicitly down-rank these; lumping them with `Vec::clone` is a foot-gun and the most common false-positive in naïve clone audits.
- **Test files, error paths, and one-shot startup code are Tier D.** Optimising clones there is unhelpful churn.

## Tier rubric

Every finding falls into one of four tiers. **Only Tier A and B are surfaced by default**; Tier C is surfaced with a "measure first" prefix; Tier D is never reported.

### Tier A — Always fix
*High-confidence speedup, trivial fix (≤ 5 lines, no API change).*

Examples:
- Clone of an already-moving value: `let x = ...; let y = x.clone(); consume(y);` where `x` is unused after.
- `.to_string()` immediately re-borrowed as `&str`: `func(&s.to_string())` where `func` accepts `&str`.
- `.clone()` on a `Copy` type (clippy's `clone_on_copy` catches some, not all).
- `iter().map(|x| x.clone()).collect()` → `.cloned()` (or just `.iter()` if the caller can borrow).
- Cloning a key just before it goes out of scope to look it up in a map.

Surface immediately with a one-line diff and a one-sentence justification.

### Tier B — Fix if cheap
*Likely small or context-dependent speedup, fix involves a localised signature relaxation or borrow-split.*

Examples:
- Owned `String`/`Vec<T>` parameter where `&str`/`&[T]` would do, callers all have owned values to pass.
- `String::from(&s)` to convert `&String` → `String` when the API can be relaxed to `&str`.
- Cloning to lengthen a borrow that could be split-borrowed via destructuring (see *Pattern: field-level destructure* below).
- `format!("{}", x)` where `x: Display` and the result is immediately re-borrowed.

Surface with the relaxation diff. If the API change cascades to multiple call sites or crosses a `pub` boundary, **downgrade to Tier C**.

### Tier C — Needs measurement
*Speedup uncertain, fix complexity is non-trivial (lifetime annotations across module boundaries, `Cow`, `Rc`/`Arc`, ownership restructuring, or `Send`-affecting changes in async code).*

Surface with explicit "benchmark first" framing and a sketch of the criterion harness. Do **not** apply the fix without measurements; the refactor cost may exceed the speedup, and the new shape may regress readability.

### Tier D — Leave alone
*Cost is structurally bounded or the fix is high-risk for unmeasured benefit.*

Examples:
- `Arc::clone` / `Rc::clone` — refcount bump (atomic increment, no allocation). Tier D unless inside a documented hot inner loop.
- Clones in `Drop`, error paths, panic handlers, or test fixtures.
- One-shot startup code (config loading, CLI parsing).
- Clones whose removal would require `unsafe`, lifetime parameters on public APIs, or breaking changes for unmeasured benefit.

**Do not report Tier D findings.** If you skipped many, mention the count and category in the review introduction: "Skipped 12 `Arc::clone` refcount bumps and 4 test-fixture clones as Tier D."

## Footgun catalog

Five patterns clippy misses or mis-handles. Call them out by name when you see them.

### 1. Clones on already-moving values across split control flow
```rust
let x = build();
match cond {
    A => use_a(x.clone()),
    B => use_b(x.clone()),
}
// `x` unused after — both clones redundant.
```
If `x` is unused after the `match`, both clones are redundant; the value can move into whichever arm fires:
```rust
let x = build();
match cond {
    A => use_a(x),
    B => use_b(x),
}
```
Clippy's `redundant_clone` misses this when the move-out points are split across arms.

### 2. `.to_string()` immediately re-borrowed
```rust
fn caller(s: &str) {
    inner(&s.to_string());      // allocates a String, then re-borrows as &str
}
fn inner(s: &str) { ... }
```
If `inner` already accepts `&str`, drop the `to_string`:
```rust
fn caller(s: &str) {
    inner(s);
}
```
Same pattern with `format!("{}", x)` where `x: Display` and the result is immediately read as `&str` — replace with a direct `Display` call, or pass `x` through.

### 3. `iter().map(|x| x.clone()).collect()` instead of `.cloned()`
```rust
let owned: Vec<String> = strings.iter().map(|s| s.clone()).collect();
```
Replace with `.cloned()` when the output really needs to be owned. Better: return the iterator and let the caller decide whether to borrow or own:
```rust
fn each_string(&self) -> impl Iterator<Item = &str> + '_ { ... }
```
Use `.copied()` for `Copy` types — same shape.

### 4. `Arc::clone` / `Rc::clone` mis-classified as a real clone
These are refcount bumps: an atomic increment, a decrement on drop, no allocation. Cost is nanoseconds. **Do not flag as a finding** unless they sit inside a documented hot inner loop where the surrounding work is also nanosecond-scale. Even then, prefer to share a `&Arc<T>` borrow over taking another reference count when the API permits.

### 5. Clones used to dodge a borrow-checker error
```rust
let data = obj.field.clone();   // workaround for &/&mut conflict
do_thing(&data);
obj.other_method();
```
The clean fix often involves restructuring: split-borrow via destructuring, hoist the read out of the conflict's scope, or reshape ownership. **This is Tier C — recommend benchmarking first.** If the surrounding work is non-trivial, the clone may already be cheap relative to it, and the refactor is risky.

Worked example: this project's pileup `process_position` (Mi6 fix in `ia/reviews/pileup_2026-05-09.md`) replaced `rec.ref_seq.clone()` per affected record with field-level destructure of `&mut OpenPileupRecord`, saving up to a 5 KB clone per record per walker step. The fix added two lines and changed one helper signature — Tier B because the speedup was bounded, the fix was small, and a regression test already covered the affected path.

## Pattern: field-level destructure

When a clone exists because one borrow blocks another on the same struct, the borrow-checker accepts independent borrows of *disjoint fields*. Destructure through `&mut` to expose them:

```rust
let rec = map.get_mut(&key).unwrap();
let MyStruct { read_field, mut_field, .. } = rec;
// read_field: &mut FieldA, mut_field: &mut FieldB — independent borrows.
// Can pass `read_field` to a function (coerced to &FieldA) while mutating
// `mut_field` in the same scope.
```

Cost: one extra line, no API change. Benefit: drops a `Vec`/`String` clone per call site. **Tier B** by default; **Tier A** if the call site is inside a hot loop.

## Benchmark idiom (Tier C only)

When recommending a Tier C fix, instruct the author to verify with criterion. The standard pattern:

```bash
# In the project's container/dev shell:
cargo bench --bench <name> -- --save-baseline before
# Apply the candidate fix.
cargo bench --bench <name> -- --baseline before
```

Decision rule (criterion's defaults):
- "Performance has improved" with p < 0.05 and change beyond the 2% noise threshold → accept the fix.
- "No change detected" → revert; downgrade to Tier D. The clone wasn't load-bearing.
- "Performance has regressed" anywhere → revert. The fix introduced overhead elsewhere (extra indirection, worse cache behaviour).
- Sub-microsecond benchmarks: require ≥ 5% change and a realistic input size guarded by `black_box(...)` to avoid loop-invariant code motion or dead-code elimination silently invalidating the result.

If no relevant bench exists, the **first** recommendation is "add a criterion bench for this hot path" — not "remove the clone."

## Reporting policy

### For each Tier A or B finding
- `path/to/file.rs:LINE` — **Tier A** / **Tier B**
- One-line problem statement (what the clone is, why unnecessary).
- One-line fix (diff or replacement snippet).

### For Tier C
- `path/to/file.rs:LINE` — **Tier C — measure first**
- One-line problem statement.
- Sketch of the criterion harness or pointer to an existing bench.
- Trailing line: "Apply the fix only if criterion reports a real improvement."

### Style nits
Group small style observations (`.cloned()` vs `.copied()`, formatting consistency) into a single "Style nits" section. Do not enumerate as separate findings.

### Tier D
Skip entirely. Optionally: a one-line summary in the review intro naming the count and category.

## Stop-and-report conditions

If you encounter any of the following, **do not silently work around them** — flag explicitly:

- A clone on an `unsafe` boundary (FFI handle, raw pointer container). Removing it might violate a soundness invariant.
- A clone on a type whose `Clone` impl has documented side effects (rare but real — some database handles do this).
- A clone whose removal would require changing a `pub` API: surface the trade-off explicitly; do not silently rewrite the public surface.
- More than ~10 clones in one function: the function likely needs a structural rethink, not a clone-by-clone audit. Recommend a separate refactor and stop the audit on this function.

## What this skill is not

- **Not a clippy replacement.** Run `cargo clippy --all-targets --all-features -- -W clippy::redundant_clone -W clippy::clone_on_copy -W clippy::needless_pass_by_value` first; clippy catches the obvious cases. The skill's value is in the cases clippy doesn't see and the cases clippy flags but the right answer is "leave it."
- **Not a benchmarking framework.** The skill recommends measurement; the author runs the measurements.
- **Not a Cow / SmallVec / Rc tutorial.** If a fix would benefit from those types, mention them by name and link to their docs (`std::borrow::Cow`, the `smallvec` crate, `std::rc::Rc`); the skill body stays focused on judgement, not type primer.

## Reusable invocation prompt

> Audit the following Rust code for unnecessary clones per the **rust-clone-audit** skill. Apply the four-tier rubric: surface Tier A and B as findings; surface Tier C with explicit "measure first" framing; suppress Tier D entirely. For each finding give a file:line, a one-line problem statement, and a concrete fix. Do not report `Arc::clone` / `Rc::clone` unless they are inside a documented hot loop. Run `cargo clippy -- -W clippy::redundant_clone -W clippy::clone_on_copy` first to filter out the obvious cases. If you find yourself wanting to apply a Tier C fix, stop and recommend a criterion bench instead.

---

## Sources

- [agentskills.io specification](https://agentskills.io/specification) — frontmatter rules, single-file 500-line guidance.
- [anthropics/skills — skill-creator SKILL.md](https://github.com/anthropics/skills/blob/main/skills/skill-creator/SKILL.md) — pushy-description pattern, progressive disclosure.
- [anthropics/claude-code — pr-review-toolkit/code-reviewer.md](https://github.com/anthropics/claude-code/blob/main/plugins/pr-review-toolkit/agents/code-reviewer.md) — confidence-threshold reporting rubric (adapted to the four-tier model here).
- [`std::clone::Clone` documentation](https://doc.rust-lang.org/std/clone/trait.Clone.html) — the load-bearing line "always explicit and may or may not be expensive" justifies the measurement tier.
- This project's pileup Mi6 fix in [`reviews/pileup_2026-05-09.md`](../reviews/pileup_2026-05-09.md) — worked example of the field-level destructure pattern (Tier B, ~5 KB clone per call dropped).
