---
name: spec-authoring
description: Use this skill when writing or revising a design spec — a "what to build and why" document under doc/devel/**/spec/ (a proposal, a design spec for a module/step/feature, or a decision record). It defines the structure and discipline of a spec: draw the boundaries, decide with recorded rationale, defer nothing silently, and track open questions honestly. Trigger on "write a spec", "draft the design for X", "spec out this feature/module/step", or when settling a design before any code. Prose mechanics defer to the clear-technical-writing skill; the next doc downstream is architecture-authoring.
---

# Design-spec authoring

A **spec answers *what* to build and *why*.** It is the home of rationale — the place a
decision and its reasoning live so no one has to reconstruct them later. It settles the
design; it is *not* code-facing (that is the architecture doc, `architecture-authoring`) and
*not* a build order (that is the implementation plan, `implementation-plan-authoring`).

The canonical exemplars in this repo are
[`doc/devel/ng/spec/read_filtering.md`](../../../doc/devel/ng/spec/read_filtering.md) and
[`doc/devel/ng/spec/ref_seq.md`](../../../doc/devel/ng/spec/ref_seq.md). Read one before
writing a new spec — match its shape.

## Prose comes first — defer to `clear-technical-writing`

A spec is read-facing prose for a domain expert. The
[`clear-technical-writing`](../clear-technical-writing/SKILL.md) skill governs *how* to write
every sentence, name, and term; **invoke it, don't restate it.** The five rules that bite
hardest in specs:

- **Explain before you formalize.** Lead each section with a plain-English "what this does /
  why", *then* the type, formula, or signature. A symbol-only summary line is the signal a
  plain lead is missing.
- **A name says what the value *is*.** In the types a spec introduces, prefer the domain term;
  no single-letter names outside a worked derivation.
- **Jargon once, never stacked.** Define each acronym/term on first use; a "Jargon, once:"
  aside or a short glossary earns its keep when a doc carries several.
- **Relocate precision, don't delete it.** Plain lead in the body, exact detail in the next
  sentence or a code block. Both survive.
- **The reader test:** would a domain expert who has never seen this understand it on the first
  pass, without looking anything else up? If no, it is not done.

## What a spec is *for* (the discipline)

- **Goals and non-goals — state both.** List what the design must achieve, and — just as
  valuable — the **non-goals**: things that could reasonably be goals but are deliberately
  excluded. Pair them with the scope-fixing properties and an explicit **"does not"** list.
  Boundaries drawn on paper stop the design creeping in code. (read_filtering §1: three
  properties + a "does not" list.)
- **Decide, and record why — with the alternatives.** Every real choice is a decision *plus its
  reasoning*, and for choices that had a genuine fork, the **options considered and why the
  rejected ones lost**. A record without its rationale loses its value the moment circumstances
  change and no one can tell whether it still applies. (read_filtering §7 records each decision
  against the alternative it beat — e.g. `RecordSource` vs. an iterator of records.)
- **Show the trade-offs.** The design section is where trade-offs live: given the goals and
  non-goals, *why this shape best satisfies them*. A spec that asserts a design without weighing
  the tension it resolves is thin.
- **Defer nothing silently.** Out-of-scope work gets a **"Deferred, with a recommended home"**
  entry — where it *should* live and why not here. Nothing vanishes without a forwarding address.
- **Address the cross-cutting concerns, briefly.** Name how the design meets the ones that
  apply — for this repo most often performance/memory, the error model, and concurrency — a
  sentence or two each, not a treatise.
- **Track open questions honestly.** Resolved items keep the reasoning that closed them; open
  ones carry a leaning and "confirm before code". Closing a question means writing down the
  decision, not deleting the question.
- **Reuse over rewrite.** If the design ports or builds on existing code, map it: a table of
  *what → existing code → how it is reused*. Name the parity oracle.

## When a spec is warranted

Write one when the design has real trade-offs, cross-cutting impact, or decisions worth
recording for later. **Skip it when the solution is obvious and low-trade-off** — a spec that
reads as an implementation manual with no alternatives weighed is a sign the work should have
gone straight to an arch doc or to code. Match length to stakes: roughly a page for incremental
work, more for a foundational module.

## Structure (adapt to the subject; keep the spine)

1. **Header / status block.** Date(s), what changed and when, cross-refs to the companion arch
   and plan docs, and a one-line state (e.g. "No code yet — this settles the design"). Update
   the status line when the design shifts materially — and when a decision *reverses*, record the
   supersession there rather than silently rewriting the old text (an append-only decision trail,
   ADR-style).
2. **What it is — goals, non-goals, and what it is not.** The must-achieve list, the
   deliberately-excluded non-goals, the scope-fixing properties, and the explicit "does not" list.
3. **Foundations / conventions this doc establishes** (if it is the first in a series): the
   newtype rules, error pattern, config conventions, module-placement rule it pins down.
4. **The design, section by section.** Each section: plain-English lead → the decision → the
   precise detail (a `rust` block for types/traits, a table for a set). Decisions carry their
   rationale — and the alternative they beat — inline.
5. **The types** (reference blocks): only the types this doc actually touches; explanation
   above each block.
6. **Reuse map** to existing code (`| what | existing code | reuse |`).
7. **Deferred, with a recommended home.**
8. **Open questions** — resolved (with reasoning) and open (with leaning). Retitle to
   "Resolved decisions & open questions" once most are closed.

Not every spec needs every section; keep the ones that carry weight for the subject.

## Conventions the repo's specs follow

- **Decision records read `— decided` / `resolved: X` / `OPEN`.** A resolved item states the
  choice *and* what was rejected and why (read_filtering §7). An open item states the options,
  a leaning, and "confirm before code".
- **`Option<T>` is "absent", not a sentinel; defaults are named `pub const`s with their source;
  types are validated newtypes when constrained.** The spec names these conventions; the arch
  doc and code inherit them.
- **Boy-scout notes** flag choices deliberately left to the implementer to improve for
  readability — distinct from open *design* questions.
- **Cross-reference precisely.** Link the companion arch/plan docs and any sibling specs by
  relative path; when a fact is settled elsewhere, point rather than restate.

## Revision pass (before calling a spec done)

1. Run the **clear-technical-writing revision pass** over the whole doc (names, jargon, section
   openings, mechanics, the reader test).
2. Are **goals and non-goals** both stated, and is there a "does not" list where scope could creep?
3. Does every **decision** carry its *why* and the **alternative it beat**? Could a reader
   reconstruct the reasoning and the trade-off?
4. Is anything out-of-scope dropped **without a recommended home**?
5. Are the applicable **cross-cutting concerns** (performance/memory, errors, concurrency)
   addressed, at least briefly?
6. Are **open questions** honest — resolved ones show their reasoning, open ones a leaning?
7. Is the **reuse map** present and is the parity oracle named, if this is a port?
8. Is the **status line** current, with any reversed decision recorded (not silently rewritten)?

The bar: a spec the owner never has to send back with "*why* did we decide that?" — and that
the architecture doc can be distilled from without re-deciding anything.

## Sources & influences

- **Google "Design Docs"** (industrialempathy.com/posts/design-docs-at-google) — context &
  scope, **goals and non-goals**, **alternatives considered**, cross-cutting concerns, and the
  rule to *skip the doc when the solution is obvious*.
- **Architecture Decision Records** (Nygard's "Documenting Architecture Decisions"; MADR,
  adr.github.io) — status + context + decision + consequences, **considered options with their
  trade-offs**, and the discipline that a decision without its rationale rots.
- **Diátaxis** and the project's [`clear-technical-writing`](../clear-technical-writing/SKILL.md)
  skill — explanation vs. reference, one term per concept, the reader test.
- Decisive, though, are the repo's own exemplars (`read_filtering.md`, `ref_seq.md`); the
  external frameworks corroborate them, they don't replace them.
