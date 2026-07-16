---
name: spec-authoring
description: Use this skill when writing or revising a design spec — a "what to build and why" document under doc/devel/**/spec/ (a proposal, a design spec for a module/step/feature, or a decision record). A spec has two readers, the coder who will build it and the manager checking the right thing gets built, and answers their five questions and nothing else: what are we building, how does it relate to the rest of the project, what will bite the coder, what is decided vs open, and how do we know it works. The discipline: record facts instead of manufacturing rationale, ground every claim about code in the code, defer nothing silently, and treat length as a correctness property — an unread spec has failed. Trigger on "write a spec", "draft the design for X", "spec out this feature/module/step", or when settling a design before any code. Prose mechanics defer to the clear-technical-writing skill; the next doc downstream is architecture-authoring.
---

# Design-spec authoring

A **spec answers *what* to build and *why*.** It is the home of rationale — the place a
decision and its reasoning live so no one has to reconstruct them later. It settles the
design; it is *not* code-facing (that is the architecture doc, `architecture-authoring`) and
*not* a build order (that is the implementation plan, `implementation-plan-authoring`).

## Who you are writing for — get this wrong and everything else follows

**Two readers, both real, neither of them a reviewer:**

- **the coder** who will write the software. They know the domain better than you. They need to
  build the thing without losing a day to something you knew and didn't say.
- **the manager** organising the project, making sure the *right* features get built. They need to
  know what this costs, what it depends on, what else it moves, and what is still unknown.

Between them they ask **five questions, and a spec answers those and nothing else:**

| | question | who asks |
|---|---|---|
| 1 | **What software are we going to build?** | both |
| 2 | **How does this piece relate to the rest of the project** — what it stands on, and what else must change? | manager: scheduling · coder: dependencies |
| 3 | **What will bite me?** — the traps that cost a day if nobody says them | coder |
| 4 | **What is decided, and what is still open?** | manager: future work + risk · coder: don't build on it |
| 5 | **How do we know it works?** | coder: the tests · manager: the definition of done |

**Anything serving neither reader is there to make the author look rigorous.** That is the cutting
test, and it is sharper than "is this a fact or an argument?" because it also catches *facts nobody
needs*. Run it paragraph by paragraph: **which of the two is this for?** No answer → delete it.

### Why the second reader is load-bearing, not decoration

A coder-only spec cannot catch a whole class of error — the kind where the instruction is clear and
the *purpose* is gone. Real case from this repo: a spec for a module inside a research lab wrote
**"Goal: reuse the catalog's admission policy."** To the coder that reads fine — it is precise and
actionable. To the manager it is a disaster: the module exists to *question* that policy, and the
sentence quietly abandons the experiment the whole thing was funded for. **The coder had no way to
see it.** The manager would have caught it in a sentence.

So when you write a goal, ask the manager's question too: *does this still describe the thing we are
trying to find out?*

### The bar — one per reader

- **Coder: does this save them a day at the keyboard?**
- **Manager: can they tell what it costs, what it blocks, and whether it is still the right thing to
  build?**

**Not "could this survive review."** That bar is what produces defensive specs: it puts an adversary
in the room, and against an adversary "we inherited this and never measured it" feels like a
weakness to dress up — so you dress it up, and now the doc is long, defensive, and in places
invented. To *either* real reader that same sentence is exactly what they needed. A spec that
survives review and goes unread has failed completely; a scrappy list of traps that gets read has
succeeded.


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

  **But this applies to choices that were actually made, and alternatives that were actually
  live.** It is *not* a licence to attach an argument to every sentence. Two ways it goes wrong,
  both worse than saying nothing — see *Facts, not arguments* below: manufacturing a rationale for
  something that never had one, and inventing a strawman alternative in order to defeat it.
- **Show the trade-offs — where there are any.** Where the design resolves a real tension, say what
  it traded. Where it doesn't, don't stage one. A spec that asserts a design without weighing a
  genuine tension is thin; a spec that manufactures tension to look rigorous is worse, because now
  the reader has to work out which of your trade-offs are real.
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

### Facts, not arguments

**The single most common failure in this repo's specs, and the one that makes them long *and*
wrong.** A spec is a record, not a case. Write down what is true and what was decided; attach a
reason only where one exists.

- **Never manufacture a rationale.** When a value is inherited, guessed, or copied from a tool that
  never explained itself, **say that**. *"Period 2–6: the catalog's setting, inherited, never
  measured"* is a **complete entry** — and a far more useful one than an invented justification,
  because it marks the decision as **soft**. In a lab, "this number is a guess" is the most
  actionable line you can write: it says what is safe to move. An invented reason freezes it.
- **Watch for the tell.** If you are writing a *because* and you cannot point at a measurement, a
  constraint, or a line of code — stop. You are about to invent one. Real examples of the failure:
  "confirming an existing decision with new reasoning" (there was no reasoning — it was inherited);
  "`bundle_threshold` exists to serve `flank_bp`" (they were two unrelated histories that collided
  at the same number).
- **Do not promote an expedient into a goal.** "Match the catalog so the parity test works" is a
  scaffolding decision. Written up as *"Goal: reuse the catalog's admission policy"* it inverted
  the purpose of the entire module, which exists to *question* that policy.
- **Softness is content.** Mark what is inherited, unmeasured, guessed, or a starting value. That
  is the map of where the design can move — often the most valuable thing in the document.

### Ground every claim in the code

**Do not write what you have not read.** In practice the split is stark: the lines that earn a
spec's keep are the ones found by opening files — a `contig_seq.len()` that silently means
something else inside a window, a function that is `cfg(test)`-gated so production cannot call it,
a `u32::try_from(..).unwrap_or(u32::MAX)` that clamps instead of failing. The lines that embarrass
it are the ones reasoned out from memory.

- Cite `file:line` for claims about existing code. If you did not look, either look or cut the
  claim.
- **Reasoning is cheaper than verifying, which is exactly why it is dangerous.** When you are
  writing to justify (above), you need reasons faster than you can check them — the two failures
  compound.
- Trace the code before asserting what it does. Claims like "this rule is over-aggressive" or "the
  consumer doesn't care how this arrives" are checkable, and were both wrong when checked.

## When a spec is warranted

Write one when the design has real trade-offs, cross-cutting impact, or decisions worth
recording for later. **Skip it when the solution is obvious and low-trade-off** — a spec that
reads as an implementation manual with no alternatives weighed is a sign the work should have
gone straight to an arch doc or to code.

### Length is a correctness property, not a matter of taste

An unread spec has failed regardless of what is in it. This is not a style preference — it is the
difference between the document doing its job and not.

**The anchor is the exemplars: `ref_seq.md` is 342 lines, `read_filtering.md` is 594.** Those are
foundational modules. If you are past ~600 for a step that mostly orchestrates existing code, you
are not being thorough — **you are arguing**, and every manufactured justification costs five to ten
lines (see *Facts, not arguments*). One draft in this repo reached 2,043 lines and the owner's
verdict was *"I haven't read it. I'd rather write the code myself."* That is the failure mode, and
it is total.

Two traps on the way there:

- **Thoroughness is not quality.** Alternatives-considered, supersession records, decision tables —
  each defensible alone, collectively a document nobody opens. Do not optimise the proxy.
- **Length is a symptom; treat the cause.** When told a doc is too long, the reflex is to trim a
  section. The cause is almost always that you are justifying rather than recording — cut *that*
  and the length follows.

## Structure (adapt to the subject; keep the spine)

1. **Header / status block.** Date(s), cross-refs to the companion arch and plan docs, and a
   one-line state (e.g. "No code yet — this settles the design"). Keep it to a few lines.

   **Record a supersession only once the reversed decision was *acted on*** — code shipped against
   it, or another doc written to it — because then someone needs the archaeology and git will not
   tell them *why*. **While the design is still being argued out, git is the history: fix the text
   and move on.** The append-only ADR trail is for decisions that escaped into the world, not for a
   draft. Applied to a doc revised five times in one afternoon before a line of code exists, it
   produces a status block longer than some sections and a §-per-reversal of pure noise. If in
   doubt: has anything downstream been built on the old decision? No → just change it.
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

**Cut first.** Every other check below adds; if you only ever run those, the doc only ever grows.

1. **Cut what serves neither reader.** Go paragraph by paragraph and ask **which of the two is this
   for — the coder or the manager?** No answer, delete it. This catches more than "am I defending a
   choice?" does, because it also catches facts nobody needs. Watch specifically for: a *because* you
   cannot ground; an alternative nobody proposed; a trade-off you staged; a heading that announces
   what you are about to explain instead of explaining it; a paragraph proving that a definition is a
   definition; a supersession record for a decision nothing was built on.
2. **Check the length against the exemplars** (342 / 594). Over it? You are arguing — go back to 1.
3. **Check every claim about code.** `file:line`, or you looked. If neither, cut it.
4. Run the **clear-technical-writing revision pass** (names, jargon, section openings, mechanics,
   the reader test).
5. Are **goals and non-goals** both stated, and is there a "does not" list where scope could creep?
   Is anything in here **someone else's decision** — what a downstream consumer does with your
   output is not yours to make.
6. **Read the goals as the manager.** Do they still describe the thing we are trying to find out, or
   have they quietly become "do what the existing code does"? This is the failure a coder cannot
   see (§*Who you are writing for*).
7. Does every **decision that had a real fork** carry its *why* and the alternative it beat? And is
   every value that *didn't* — inherited, guessed, copied — **marked as soft** rather than dressed
   up (*Facts, not arguments*)?
8. Is anything out-of-scope dropped **without a recommended home**?
9. Are the applicable **cross-cutting concerns** (performance/memory, errors, concurrency)
   addressed, at least briefly?
10. Are **open questions** honest — open ones with a leaning and the **measurement that would settle
   them**?
11. Is the **reuse map** present and is the parity oracle named, if this is a port?

**The bar, one per reader:** does this save the *coder* a day at the keyboard, and can the
*manager* tell what it costs, what it blocks, and whether it is still the right thing to build?

Not "could this survive review" — that bar is what produces defensive, unread specs, because it puts
an adversary in the room instead of two colleagues. The trap list earns its keep. The paragraph
defending a name never did.

A good self-check: **which lines here would I have found only by opening a file?** If the answer is
"few", the spec is mostly opinion, and its length is mostly noise.

## When a spec fails, fix this skill

Every concrete rule above — the manufactured rationale, the append-only trail on an unbuilt draft,
the 2,043-line doc, the goal that abandoned its own experiment — is here because a real spec in this
repo failed that way and the owner said so. That is why they bite: they are not imported best
practice, they are this project's own scar tissue.

So the rule is standing: **when a spec gets sent back, the spec is only half the fix.** Ask what in
this file permitted it — usually there is something, and usually it is a good instruction applied
without judgement. Being wrong is a net gain if the process improves; being wrong twice the same way
means it didn't.

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
