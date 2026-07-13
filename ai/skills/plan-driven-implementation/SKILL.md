---
name: plan-driven-implementation
description: Use this skill to execute an existing written implementation plan step by step, running a rigorous per-step loop that composes the other skills — implement (rust-feature-implementation) → review (rust-code-review) → apply fixes (apply-code-review-fixes) → commit — and pausing for the user whenever reality departs from the spec, architecture, or plan. Trigger on "implement the plan", "follow the implementation plan", "run the plan step by step", "build <feature> following its impl plan", or being handed a plan doc to execute.
---

# Plan-driven implementation

Drive a feature from a **written implementation plan** to committed, reviewed code, one step
at a time, composing three existing skills per step.

A plan is a plan, and an architecture is a **proposal** — not a spec to follow to the letter.
Writing the code is where their gaps and rough edges show, so the implementer keeps normal
latitude: **minor deviations that preserve the intent are expected** — made and recorded as
you go, not escalated. What this skill will not do is *redesign*: when a divergence is large
enough to change the design or ripple beyond the current step, it **stops and asks the user**
rather than quietly reshaping the plan. The discipline is not "never deviate" — it is
**"deviate small and say so; escalate big."**

It orchestrates, in order, for every step:

1. **Implement** — [`rust-feature-implementation`](../rust-feature-implementation/SKILL.md) (types first).
2. **Review** — [`rust-code-review`](../rust-code-review/SKILL.md) on the step's diff.
3. **Apply fixes** — [`apply-code-review-fixes`](../apply-code-review-fixes/SKILL.md) on that review.
4. **Commit** — one commit for the step, once validation is green.

Then the next step, until the plan is done.

## Preconditions (stop if not met)

- **A written implementation plan exists** with ordered steps and milestones (e.g. under
  `doc/devel/**/impl_plan/`). If there is no plan, stop and say so — writing the plan is a
  separate, prior task, not this skill's job.
- **The design is settled** — the plan cites a spec and architecture that resolve the *what*
  and *how*. This skill fills in *code*, not open design questions.
- **The container is the autonomy tool — work there freely.** The dev container
  (`./scripts/dev.sh`, per `CLAUDE.md`) exists to grant full permissions to build, test, and
  run tooling **without prompting the user at every step**. Default rule: if a step needs
  something the host won't permit, do it in the container rather than stopping to ask —
  permission and tooling friction is *never* a reason to interrupt. The user is interrupted
  only by a **plan/design divergence** (stop-and-ask). Asking for permission is allowed, but only
  when the permissions granted in the container are not enough.

At the start, read `PROJECT_STATUS.md` and the plan in full, and confirm the plan's
preconditions ("already in place") actually hold in the current code. Note the plan's
milestone checkpoints — they are human-review pause points.

## Core principles

- **The plan guides; it doesn't bind to the letter.** Follow its step order, granularity,
  and intent by default. Small local improvements an implementer would normally make — 
  a helper split differently, a better-fitting data structure, adapting to a
  reused API's real shape — are fine: make them and **record them** (in the step's impl
  report and commit). Reordering, merging, or dropping steps, or inventing scope, is a bigger
  move — reserve it for a recorded reason, and if it touches the design, stop and ask.
- **One step, one loop, one commit.** Each step goes fully through implement → review →
  apply → commit before the next begins. Do not batch steps through a single stage.
- **Absorb small divergences; escalate big ones.** The judgment call at the heart of this
  skill: a *minor* departure that keeps the intent is the implementer's to make and note; a
  *significant* one — something the design didn't anticipate, a contradiction between the plan
  and the spec or the actual code, a change that ripples beyond this step — is a
  **stop-and-ask**. When unsure which side of the line you're on, treat it as significant and
  ask. Either way, **no departure is silent** — small ones are recorded, big ones are raised.
- **Delegate, don't duplicate.** The three composed skills own their own quality bars,
  reports, and `PROJECT_STATUS.md` updates. This skill sequences them and centralizes the
  decision to pause; it does not re-implement their checklists.
- **Real validation only.** A step is not committed until its build is actually green in the
  container. Never fabricate command output or wave a step through.

## The per-step loop

For each step, in plan order:

### 1. Implement
Invoke `rust-feature-implementation` for **this step only**. Hand it:
- the step's contract from the plan, and the exact spec/architecture sections the step cites;
- the reuse targets the plan names (existing types/functions to build on).

**Suppress its interactive planning gate.** The plan *is* the approved design, so do not have
the sub-skill re-present a plan and re-ask for comments — that gate is already satisfied. Its
**"no silent assumptions" rule still holds**, though: if implementing the step forces a choice
the plan left unspecified *and that choice changes direction*, that is a **stop-and-ask**
(below), not a silent default.

Types first, tests included — both are already this skill's discipline. Output: the step's
source + tests, an implementation report, `PROJECT_STATUS` at `implemented`.

### 2. Review
Invoke `rust-code-review` scoped to the step's **working-tree diff** (the uncommitted changes
from stage 1). It fans out its per-category sub-agents, writes a review report, and moves
`PROJECT_STATUS` to `reviewed`. Do not pre-filter its findings.

### 3. Apply fixes
Invoke `apply-code-review-fixes` on that review report. It applies the `Apply`-class findings
conservatively, writes a fix-application report, and moves `PROJECT_STATUS` to `fixes-applied`
(or `shipped`).
- If it raises an **`Ask`** finding (a real user decision), or a **Blocker that reveals a
  design flaw** rather than a coding slip → **stop-and-ask**. Do not guess.
- If the review found nothing actionable, this stage is a recorded no-op — note it and proceed.

### 4. Commit
One commit for the step, **only once validation is green in the container**
(`cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`,
`cargo test --all-targets --all-features`). Commit together:
- the step's source + tests,
- the three reports (impl / review / fixes-applied),
- the `PROJECT_STATUS.md` updates.

Message: `<type>(<scope>): <step-id> — <step title>` — e.g.
`feat(ng): A3 — InMemoryRefSeq + tests`. End with the repo's `Co-Authored-By` trailer.

If validation cannot be made green without a design decision → **stop-and-ask** (do not commit
red, and do not weaken the checks to pass).

### Then
Mark the step **done** in the plan doc (flip its `☐` to `✅`), and move to the next step.

## Step granularity and milestone checkpoints

- **A "step" is the plan's smallest committable unit.** Follow the plan's granularity. If the
  plan's steps are so fine that a stage would be near-empty (e.g. a pure-scaffold step),
  tightly-coupled adjacent steps *may* share one loop iteration — but only explicitly, named
  in the commit, never silently.
- **Milestone checkpoints are hard pauses.** When the plan marks a checkpoint ("pause for
  review"), stop after committing that milestone's last step and hand back to the user before
  starting the next milestone. Summarize what landed and what's next.

## On plan completion

When the plan's last step is committed:

- **Also verify on the host.** Every step was already validated in the container; as a final
  check, run the build and test suite **natively on the host** too, to confirm the code works
  outside the container environment. `cargo` is host-allow-listed, so this needs no prompt; if
  the host toolchain is unavailable, note that and skip rather than failing the run.
- **Summarize the run** — steps landed, commits, reports produced, and any deferred
  follow-ups — and hand back to the user.

## Latitude vs. stop-and-ask

The plan and architecture are proposals; implementing them is where judgment lives. Choose a
response to any departure by its **size and blast radius**, not by whether it's a departure at
all.

**Absorb and note** — the implementer's latitude. Make the change, record it in the step's
impl report and commit message, and keep going:

- a clearer name, a helper factored differently, a better-fitting data structure or signature;
- adapting to a reused API's *actual* shape (a slightly different signature or return than the
  plan or an architecture sketch assumed) when the intent still holds;
- a local ordering tweak between trivial sub-steps;
- anything the code review would treat as a normal implementation choice.

**Stop and ask** — the divergence changes the design or reaches beyond the step. State
precisely what diverged, the options, and a recommendation, then wait:

- **A real gap:** the step needs a decision the spec, architecture, *and* plan all leave open,
  and the choice would shape the design.
- **A true contradiction:** the plan disagrees with the spec, or the actual code contradicts
  the plan in a way that **invalidates the step's approach** — not merely a signature to adapt to.
- **Blast radius beyond the step:** the change would alter a shared type or interface, another
  component, or a later step.
- **An `Ask` from apply-fixes**, or a **review Blocker that is really a design problem**.
- **Validation unfixable without a design decision.**
- **A milestone checkpoint** (a planned pause).

When unsure which side you're on, treat it as stop-and-ask. When you pause, do not partially
commit a broken step — leave the tree in the last green state (or clearly describe the
in-progress diff) so the user can decide cleanly.

## What this skill does *not* do

- It does **not** make *design-level* decisions — those go back to the spec/architecture via
  a stop-and-ask. (Minor implementation choices are the coder's normal latitude — see above.)
- It does **not** skip the review or apply stages to move faster, and does **not** relax
  validation to get a step committed.
- It does **not** edit the *design* in the spec, architecture, or plan — only the plan's step
  checkboxes (progress) as steps complete. A needed design change is a stop-and-ask.

## Progress tracking

- Keep a live task list (one entry per plan step) so the user can see position in the run.
- Flip the plan's step `☐ → ✅` after each committed step.
- Lifecycle status in `PROJECT_STATUS.md` is handled by the three sub-skills per step; this
  skill does not update it directly.

## Reusable prompt template

> Execute this implementation plan using the **plan-driven-implementation** skill.
>
> - **Plan:** `<path to the impl-plan doc>`
> - **Scope:** `<which milestone(s)/steps to run this session, or "through the next checkpoint">`
> - **Design authority:** `<the spec + architecture docs the plan rests on>`
>
> Requirements:
> 1. Run every step through implement → review → apply-fixes → commit, in order, one step at a time.
> 2. Follow the plan's step order and granularity; do not reorder, merge, or skip silently.
> 3. Suppress the sub-skills' own planning gates (the plan is approved); centralize pausing here.
> 4. Commit a step only when its build is green in the container; never fabricate results.
> 5. Absorb minor deviations that keep the intent (and record them); stop and ask me only when a divergence is significant — it changes the design or reaches beyond the step — and at every milestone checkpoint.

---

*Draft, 2026-07-13. Name and details open to revision — see the handoff notes accompanying this draft.*
