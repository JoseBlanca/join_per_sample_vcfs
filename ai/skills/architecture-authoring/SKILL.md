---
name: architecture-authoring
description: Use this skill when writing or revising an architecture doc — the code-facing "types & interfaces" companion to a spec, under doc/devel/**/arch/. It defines the structure of an arch doc: the module home, the types and trait signatures as they will appear in code, crisp decision records, and a reconciliation table that converges on existing types instead of inventing new ones. Trigger on "write the architecture", "draft the arch doc / interfaces for X", "turn this spec into interfaces", or after a spec settles and before the implementation plan. Rests on a spec (spec-authoring); feeds the plan (implementation-plan-authoring). Borrows naming/defaults/errors/module-structure discipline from rust-code-review.
---

# Architecture-doc authoring

An **architecture doc is the code-facing distillation of a spec.** Where the spec answers
*what and why*, the arch doc answers ***how it appears in code***: the module home, the types
and trait signatures implementations plug into, the decisions in crisp record form, and the
map from these names to existing code. It points back to the spec for every *why* — it does
not re-argue them.

Exemplars: [`doc/devel/ng/arch/ng_step_interfaces.md`](../../../doc/devel/ng/arch/ng_step_interfaces.md)
(shared vocabulary + step traits), [`doc/devel/ng/arch/module_layout.md`](../../../doc/devel/ng/arch/module_layout.md)
(the module tree), and the per-step
[`doc/devel/ng/arch/read_filtering.md`](../../../doc/devel/ng/arch/read_filtering.md). Read
the closest match before writing.

## The governing stance

> **Signatures are illustrative; the contract is the deliverable.**

An arch doc commits to *shapes and contracts*, not final code. State each type and trait as it
will appear, with a doc comment that says what it is and the invariant it holds — then let the
implementation adjust the exact signature to the real API. What must be right is the contract
(what goes in, what comes out, what holds), the names, and the module boundaries.

Corollary: **don't over-specify what code review will settle.** Verbose copied-in schemas, full
method bodies, and tuning constants belong to the code, not the arch doc — pulling them forward
is a common RFC/design-doc failure mode (detail better suited to code review than doc review).
Sketch the interface; leave the body to the implementer.

## Structure (adapt; keep the spine)

1. **Header / status block.** "architecture draft (date), companion to `../spec/<name>.md`",
   plus a one-line pointer to the shared arch docs and the naming convention. If the shared
   step-interfaces doc defers this step's "real interface", say so — this doc *is* it.
2. **Module home.** Where it lives in the `src/` tree and why (file vs folder, which module it
   shares). Defer the tree-wide rules to the module-layout doc; state only this unit's placement.
3. **The types.** The domain newtypes it seeds/uses and the local types, as `rust` blocks with
   doc comments. Explanation above each block (progressive disclosure still applies).
4. **The interfaces.** The traits and the public surface, with doc-commented signatures and a
   one-paragraph **contract** (invariants, error behaviour, laziness, ownership).
5. **Design decisions — decided.** Crisp records, one line of *why* each, cross-referenced to
   the spec section that argued it. Open items marked `OPEN:`.
6. **Reconciliation with existing code.** The convergence table (below).
7. **Open items** — impl-time confirmations and any genuinely-unresolved boundary.
8. **Test/bench shape** (brief) — where tests live and what the regression anchor is.

## Code-shape discipline (borrowed from `rust-code-review` and the codebase)

The arch doc is where the code's conventions are locked in before code exists. Apply the
naming, defaults, errors, and module-structure rules from
[`rust-code-review`](../rust-code-review/SKILL.md) (its `code_review/naming.md`,
`defaults.md`, `errors.md`, `module_structure.md` checklists) at the *design* stage:

- **Newtypes for domain scalars; nouns for types, verbs for functions.** Every same-primitive
  quantity that could be transposed gets its own type so the compiler refuses the mix
  (`ContigId`, `MapQual`, `Bp`). Constrained scalars hide their field behind a checked
  constructor; unconstrained ones keep a `pub` field and `.get()`.
- **Modules named for the concept they own**, never a layer or pattern (`allele_candidates`,
  `likelihood` — never `utils`, `common`, `services`, `models`).
- **Defaults are visible:** named `pub const`s with a doc comment giving units and source; no
  magic numbers. `Option<T>` means "absent", not a sentinel `0`.
- **Errors are explicit and future-proof:** `#[non_exhaustive]` `thiserror` enums with a
  per-variant doc comment stating when it fires (see `src/ng/ref_seq.rs::RefSeqError`). Decide
  the error surface here — fallible constructor vs. fatal, `Result` vs. panic.
- **One swappable trait per step, impls side by side** (the bake-off shape) — but a step with
  no competing implementations is a single file with no trait ceremony. State which this is.
- **Reuse over rewrite.** New modules stand on existing code; the reconciliation table names
  what they build on. Prefer a single home for a shared contract so impls "cannot drift"
  (`ref_seq.rs::validate_window`).

## Decision records — the `— decided` form

Mirror `ng_step_interfaces.md` §5: **one decision per record**, each a bolded claim + one or two
sentences of *why*, ending with the spec cross-ref. A good record states **what was rejected and
why**, not only what was chosen (the ADR discipline — the rejected alternative is what makes the
record worth keeping). Open questions use `OPEN:` inline so a reader can grep the unresolved set.
Distinguish a genuine open *design* item from an impl-time confirmation (a concrete type to pin
when coding) — the latter is not a decision.

## Revision pass (before calling an arch doc done)

1. Does every trait/type carry a **contract**, not just a signature?
2. Is there a **reconciliation row for every new name**, and is the parity/reuse target named?
   Did you check you are not minting a duplicate of an existing type?
3. Do the **decision records** each cross-ref the spec section that argued the *why* (so this
   doc stays distillation, not re-argument)?
4. Are **naming, defaults, errors, module placement** all pinned per the code-shape discipline?
5. Are **open items** split into real open design questions (`OPEN:`) vs impl-time confirmations?
6. Run the **clear-technical-writing** pass (explanation before each code block; names say what
   the value is).
7. Is the doc **cross-linked** both ways with its spec, and discoverable from any index/README?

The bar: an implementer can build from this doc and the spec without re-deciding anything, and
the implementation-plan author can turn it straight into ordered steps.

## Sources & influences

- **Architecture Decision Records** (Nygard; MADR, adr.github.io) — the one-decision-per-record
  form, status, and *considered options with trade-offs*; the arch doc's decision records are
  ADRs embedded in the interface doc.
- **Google "Design Docs"** (industrialempathy.com) — sketch APIs and data storage without
  verbose copied schemas; keep detail out of the doc when it belongs in code review.
- The project's [`rust-code-review`](../rust-code-review/SKILL.md) category checklists
  (`naming`, `defaults`, `errors`, `module_structure`) and
  [`clear-technical-writing`](../clear-technical-writing/SKILL.md) — applied at design time so
  the conventions are locked in before code exists.
- Grounding examples: the repo's own `ng_step_interfaces.md`, `module_layout.md`,
  `arch/read_filtering.md`, and the code they describe (`src/ng/`).
