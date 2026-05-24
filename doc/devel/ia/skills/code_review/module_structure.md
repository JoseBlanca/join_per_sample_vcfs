# Module structure and dependency direction

**Purpose.** Module placement and import direction encode the project's architecture. As a crate grows, back-references accumulate where a "shared" module reaches into a pipeline-stage module — silent coupling that the type system does not flag and `cargo check` will never object to. This category surfaces those.

**Triggers.** Multi-file or crate-level scope. For every top-level peer module in scope, sweep its `use` lines for back-references into pipeline-stage modules. For every type defined inside a pipeline stage, sweep its consumers across the crate.

**Skip when.** Pure single-file snippets with no peer-module context.

## Vocabulary

- **Peer module.** A top-level entry under `src/` that names what it *is* — a data model (`pileup_record`), an algorithm (`baq`), a file format (`fasta`, `psp`, `vcf`), a utility (`iter_ext`).
- **Pipeline-stage module.** A top-level entry that names what the pipeline *does* at that point (`pileup`, `var_calling`, `pop_var_caller`). Stage modules compose peer modules; the inverse should not happen.
- **Back-reference.** A peer module importing from a pipeline-stage module. The dependency direction contradicts the framing.

## Rules

- **Peer modules don't import from pipeline-stage modules.** Sweep each peer's `use crate::*` lines. Any import from a pipeline-stage path is a finding. Fix: (a) lift the depended-on type to a peer too, (b) move the dependent module inside the stage, or (c) decouple the consumer via a trait (see "trait-based decoupling" below). The trap repeats — "just move the directory" without addressing back-references gives the peer module a directory location that misrepresents its dependencies.

- **Data interchange types belong at top level.** A type defined inside one stage but consumed by three or more sites outside it — especially across stages — is the pipeline's interchange, not the stage's private type. Lift to its own peer (`src/<type>.rs` for leaf modules, `src/<type>/` when submodules will follow). Heuristic: `grep -rl <TypeName> src/ tests/ benches/ examples/`; if the consumers cross stage boundaries, the type is misplaced. This rule also covers traits, errors, and constants — anything that flows between stages.

- **`pub mod` ⇒ honest dependency direction.** A `pub` module's `use` lines surface in rustdoc and become part of the crate's public surface. A `pub` peer module importing from a pipeline-stage module exposes implementation through the API. Either reduce visibility (`pub(crate)`), decouple via a trait, or accept the honest framing "this is THIS project's implementation, not a generic library" — but make that choice explicitly, not by accident.

- **Module name encodes role.** Peer modules name what they are (algorithm / format / data model / utility). Pipeline-phase modules name what the pipeline does at that point. Nesting a data-model module inside a pipeline-phase module (e.g., a hypothetical `pileup/record/`) blurs the taxonomy — downstream consumers that never touch the walker would import a record type from a path that implies a relationship to the walker. Prefer flat top-level peers.

- **Misplaced traits.** A trait defined in module A but with all implementations in module B is in the wrong place. The trait belongs alongside its implementations (or at the boundary the implementations target). Symmetric inverse: a trait that has exactly one implementor and was *not* added to decouple a `pub` API is over-abstracted unless a second implementor is genuinely imminent.

- **`foo/mod.rs` alone in its directory.** A directory containing only `mod.rs` (no sibling submodules, no near-term plans for any) signals "expecting siblings" without delivery. Either add siblings or flatten to `foo.rs`. Same import path, less filesystem noise, more honest about the module's leaf-ness. Use `foo/mod.rs` only when the module already has or imminently will have submodules.

- **`super::super::*` references.** Cross-module references that walk multiple levels up are fragile under restructuring and obscure the dependency graph — a reader at the use-site has to trace the file's location to interpret the path. Prefer crate-absolute paths (`crate::foo::bar`) for anything beyond the immediate parent. `super::*` (one level) is fine; `super::super::*` (two levels or more) is suspicious unless the file genuinely belongs deep in its hierarchy.

- **Module inception (`foo/foo.rs`).** Triggers `clippy::module_inception` and is awkward at the use-site (`crate::foo::foo::Item`). Rename the inner file to something specific (`driver.rs`, `core.rs`, `impl.rs`). Common artefact of promoting an existing flat file into its own directory without renaming.

- **Trait-based decoupling for `pub` consumers on hot paths.** When a `pub` consumer (writer, encoder, exporter) takes a concrete input type from a pipeline stage but reads the data through a stable, finite set of accessors, propose a trait with read-only methods. Costs zero runtime under monomorphisation — the trait dispatch inlines back to a direct field access — so it avoids the per-record allocation tax of a view struct. Worth proposing when: (a) the consumer is `pub` and (b) the back-reference is the kind of coupling external consumers shouldn't see. Not worth the trait churn when the consumer is internal-only and the input type is stable.

## Cross-category interactions

- **Smells overlap.** "God struct" (`smells.md`) often co-occurs with "data interchange type misplaced" — the god struct accumulates fields because the upstream consumers each grow a need. Splitting the type can be the right move even before lifting it to a peer.
- **Refactor_safety overlap.** `super::super::*` references are also a refactor-safety hazard (they don't fail at compile time when the file moves; they fail at runtime if the path now resolves to a different item). File the finding here; cite refactor_safety as convergent.
- **Naming overlap.** When promoting a type to a peer, the module name and the type's role are scrutinised together. If a misplaced-type finding is also a naming-precision finding, cite both.
