---
name: implementation-plan-authoring
description: Use this skill when writing or revising an implementation plan — the ordered build order for a settled design, under doc/devel/**/impl_plan/. It defines the structure of a plan: scope in/out, the principles that fixed the order, preconditions, milestones of small committable steps (each with a checkbox, dependencies, and a source cross-ref), review checkpoints, and a verification table. Trigger on "write the implementation plan", "turn this design into a build order", "plan out how to build X", or after a spec and architecture are settled. A plan turns design into order; it invents no design. Its output is executed step-by-step by plan-driven-implementation.
---

# Implementation-plan authoring

An **implementation plan turns a *settled* design into build order** — the sequence of small,
committable steps that construct the thing, each pointing back at the spec/arch section it
realizes. It is **not a place for new design.** If writing the plan surfaces an unresolved
design question, stop and take it back to the spec (`spec-authoring`) or architecture
(`architecture-authoring`) — do not decide it in the plan.

The plan's consumer is the [`plan-driven-implementation`](../plan-driven-implementation/SKILL.md)
skill, which executes it step by step (implement → review → apply-fixes → commit). Write steps
as the **smallest committable units** that skill can run.

Exemplars: [`doc/devel/ng/impl_plan/foundations.md`](../../../doc/devel/ng/impl_plan/foundations.md)
and [`doc/devel/ng/impl_plan/read_filtering.md`](../../../doc/devel/ng/impl_plan/read_filtering.md).
Read one before writing a new plan — match its shape.

## Preconditions for writing a plan

- **The design is settled** — a spec and an architecture doc resolve the *what* and *how*, and
  the plan can cite them. If they are not final, the plan is premature.
- **The plan cites, it does not re-derive.** Every step names its **Source** (the spec/arch
  section it builds) and its **Depends** (earlier steps it needs). No rationale is re-argued;
  link it.

## Structure (keep the spine)

1. **Header / status block.** "draft, date"; one line stating this turns the settled design
   (link spec + arch) into build order and is **not** a place for new design. If it follows an
   earlier plan, link it.
2. **Scope — In / Out.** What this plan builds, and what it explicitly does **not** (each Out
   item pointing to the later plan that owns it — nothing dropped silently).
3. **Principles (how the order was chosen).** Name *why* the steps are in this order (below).
4. **Preconditions (already in place).** What must exist before step 1 — the reuse targets, the
   parity oracle, prior milestones done. The executor verifies these hold before starting.
5. **The steps.** Milestones (A, B, C …) of numbered sub-steps (A1, A2 …). Each sub-step:
   - a leading **`☐`** checkbox (flipped to **`✅`** as it completes),
   - one or two lines: what to build,
   - ***Depends:*** earlier steps, and ***Source:*** the spec/arch section.
6. **Checkpoints.** After each milestone, a blockquote **`> Checkpoint X: … Pause for review.`**
   — a hard human-review pause the executor honors.
7. **Verification summary.** A `| milestone | proven by |` table — how each milestone is proven
   (unit tests, byte-parity vs an oracle, integration fixture).
8. **Out of scope (next plans).** The follow-on work and where it will live.

## Principles that fix the order (state the ones that apply)

These recur across the repo's plans; name the ones that shaped *this* order:

- **Types first, then implementation** — within every milestone (project rule).
- **Simplest impl first, as the test oracle for the next.** Build the synthetic/in-memory
  version first; each earlier impl is a parity oracle for the next (foundations: InMemory →
  Resident → Windowed).
- **The algorithmic heart before the plumbing.** Build and test the decision/algorithm in
  isolation before the I/O, seams, and wiring that feed it (read_filtering: the cascade before
  the record source and iterator).
- **Reuse over rewrite.** A port calls the existing predicates/constants as-is and supplies
  only its own driver; the plan names what is reused and re-derives nothing.
- **Verify against ground truth.** The north-star test is parity with a production oracle (or
  an earlier impl), not self-consistency — a port is correct when it matches what it ports.
- **Incremental, with pauses.** One milestone, then stop for review. Steps are committable
  units, not a big-bang.
- **Container builds / ungated** (repo specifics): `cargo` via `./scripts/dev.sh`, a native host
  build at completion; gate behind a feature only if compile time bites.

## Step-writing conventions

- **A step is one committable unit** — small enough for one implement→review→commit loop, large
  enough to be worth a commit. Pure-scaffold steps may pair with an adjacent step, but say so.
- **Ordered by dependency, cheapest-firing first where order is free.** Make `Depends` explicit
  so the executor (and a reader) can see the critical path.
- **Every step cites a `Source`.** A step with no spec/arch section behind it is a design gap —
  resolve it upstream, not in the plan.
- **Checkboxes are the live progress signal.** `plan-driven-implementation` flips `☐ → ✅` as it
  commits each step; write them so that flip is unambiguous (one checkbox per committable step).
- **Name the verification for each milestone**, and prefer an external oracle over
  self-consistency.

## Revision pass (before calling a plan done)

1. Is the **design fully settled** upstream — does any step smuggle in a decision the spec/arch
   didn't make? (If so, kick it back upstream.)
2. Does **every step cite a Source** and its **Depends**?
3. Is the **order justified** by the stated principles (types-first, oracle-first,
   heart-before-plumbing)?
4. Does each **milestone have a checkpoint** and a **verification** row, with a real oracle?
5. Is every **Out-of-scope** item handed to a named later plan (nothing dropped)?
6. Are **Preconditions** concrete and checkable (the executor confirms them before step 1)?
7. Are steps sized as **committable units** with unambiguous checkboxes?
8. Run the **clear-technical-writing** pass over the prose.

The bar: `plan-driven-implementation` can pick up this plan and run it end to end — every step
buildable, testable against an oracle, and traceable to the design — without ever needing a
design decision mid-run.

## Sources & influences

- The project's [`plan-driven-implementation`](../plan-driven-implementation/SKILL.md) skill —
  the consumer of this doc; write steps as the committable units it executes (implement →
  review → apply-fixes → commit), with checkpoints as its hard pauses.
- **Anthropic skill-authoring best practices** (platform.claude.com) — the *plan → validate →
  execute* pattern and verifiable intermediate outputs, mirrored here as ordered steps proven
  against an oracle at each milestone checkpoint.
- **RFCs / delivery planning** (pragmaticengineer.com) — a plan's job is de-risked delivery, not
  new design; design belongs in the spec/arch it cites.
- Grounding exemplars: the repo's own `impl_plan/foundations.md` and `impl_plan/read_filtering.md`.
