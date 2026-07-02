---
name: clear-technical-writing
description: Use this skill whenever writing or editing prose meant for a human reader — spec/architecture/design docs, READMEs, doc comments (///), commit messages, PR and review reports, and design explanations in chat. It encodes the project owner's standing requirement: clarity and readability first, precision served without sacrificing them. Trigger whenever producing or revising natural-language explanation, naming a type/field/function, or introducing a term, formula, or acronym in reader-facing text.
---

# Clear Technical Writing

The reader is a **domain expert** — for this
project, a geneticist who knows the biology but is not a statistician
or a Rust internals specialist. Write so that reader understands it on the
first pass, without decoding jargon or looking things up elsewhere.

**The governing rule:** clarity and readability come first. Precision is
required too, but **never at the expense of clarity** — when they seem to
conflict, restructure so you get both (plain explanation first, precise
detail right after), don't drop the clarity. The technical-writing
literature calls that layering *progressive disclosure*: reveal the idea,
then the detail, then the rigour — never all at once.

**Scope.** This is about internal design/spec/architecture docs, doc
comments, commit messages, and review/PR prose — explanation and reference
written for fellow engineers. It deliberately ignores the published-web-doc
concerns (SEO, page accessibility/WCAG, readability scores, analytics): the
only metric here is whether a colleague understands the writing on first
read.

This skill exists because the same failure recurs: precision reached for
*through* shorthand, jargon, and cryptic names, leaving prose that is exact
but unreadable. The rules below name that failure's specific forms, each
with a real before/after from this codebase. They are consistent with the
Google and Microsoft developer-writing style guides, the Diátaxis
documentation framework, and plain-language guidance (sources at the end);
the project-specific examples are what make them stick.

## The reader test (run on every sentence, name, and paragraph)

> Would a reader who knows the domain but has never seen this code
> understand this **without looking anything else up**?

If no, it is not done. This is the single check that catches all the
failures below.

## Know which mode you are writing — and don't blur them

The Diátaxis framework separates documentation into four modes; design and
spec docs live in two of them, and mixing them is a common source of
confusion:

- **Explanation** — *understanding-oriented*: what the thing is, why it
  works this way, how it connects to the rest. Plain prose, framed in the
  reader's domain. ("A model that answers: at this GC content, what depth
  would one copy produce?")
- **Reference** — *information-oriented*: the precise, complete facts — the
  exact fields, types, defaults, formulas. Terse and exhaustive. (The
  `struct` block with each field's units.)

Keep the two **adjacent but distinct**: lead a section with explanation,
then drop into reference. The recurring mistake is a reference-style line
(`gc_rel = depth/(scale·curve(GC))`) standing in for the explanation that
should precede it. That is what Rule 1 forbids. (How-to and tutorial are
the other two modes — rare in this repo's design docs, but the same
no-blurring rule applies.)

## Rule 1 — Explain before you formalize

Lead with what the thing *does*, in plain English framed in the reader's
domain. Only then give the formula, type, or signature. A formula or a
shorthand phrase is **not** an explanation — it is the thing being
explained.

- ❌ "`SingleCopyCoverageModel` — fit depth∼GC curve + single-copy scale so
  `gc_rel = depth/(scale·curve(GC))`."
- ✅ "A per-sample model that answers *at a window of this GC content, what
  read depth would one copy produce in this sample?* It holds the sample's
  typical one-copy depth level and a GC-bias multiplier; dividing observed
  depth by their product gives how many copies' worth of reads are there."

The same applies to a section: open with a "what this does, in words"
paragraph **before** any equation. If you find yourself writing a one-line
summary made of symbols, that is the signal a plain-English lead is missing.

## Rule 2 — A name must say what the value *is*, not its topic

A reader hitting the identifier cold should know what the value is, not
merely what it relates to. Length is fine; obscurity is not. Prefer the
established domain term (the reader already knows it).

- ❌ `gc_rel` (tells the reader "something about GC", not what it is)
- ✅ `relative_copy_number` (the actual quantity; standard CNV vocabulary)
- ❌ `CoverageModel` → ✅ `SingleCopyCoverageModel` (says what it models)
- ❌ struct/section called `calibration` when it actually estimates a
  **prior probability** → name it for what it is: `ParalogPrior`,
  `prior_probability`. A fancy-but-wrong label is worse than a plain
  accurate one; check the label matches the actual content.

**No single-letter names in code.** `k → alt_reads`, `n → total_reads`,
`F → inbreeding_coefficient`. Single letters "explain nothing" (owner's
words). Verbs for functions, nouns for types (see
`rust-feature-implementation`'s naming standards — this skill is its prose
counterpart).

**Math symbols are the one exception, and only inside a derivation.** In a
worked formula, conventional symbols for the quantities being manipulated
(`p` allele frequency, `q` carrier frequency, `T`/`m` copy counts, `σ`) are
clearer than long names — `2pq(1−F)` beats
`2·allele_freq·(1−allele_freq)·(1−inbreeding_coefficient)`. But: define each
symbol inline on first use, and keep the durable *interface* names
(`alt_reads`, …) in the types and the API. Symbols live in the maths;
names live in the code.

## Rule 3 — No unexplained jargon or acronyms; never stack them

Define every acronym/term on first use, and in the glossary if the doc has
one. Do not pile several into one sentence.

- ❌ "EM for π over the genome-wide LRs, EB posterior σ(LR + logit π), and
  FDR q-values." (six undefined terms in one sentence — unreadable)
- ✅ A paragraph that says: the score tells us which explanation fits
  better but not how *likely* a locus is to be a paralog; for that we need
  how common paralogs are genome-wide; we estimate that rate from the data
  by an iterative loop, combine it with each locus's score to get a
  probability, and rank loci to assign a false-discovery rate. *Then*, if
  needed, the precise notation, with every symbol defined.

**Know which jargon to explain.** Match the reader: for a geneticist,
explain the *statistics and software* terms, keep the *biology* brief. Do
not explain what the reader already knows, and do not leave unexplained
what they don't.

## Rule 4 — Don't borrow opaque software jargon for a concept that has a plain description

Name the concept in plain words; don't gesture at it with a term of art.

- ❌ "a set of pure kernels"
- ✅ "a module that does the statistics, exposed through pure functions
  (same inputs → same output, no side effects)"

"Kernel" is terser but is also overloaded (stats kernel, OS kernel, matrix
kernel) and decodes to nothing for the reader. If a plain phrase says it,
use the plain phrase.

## Rule 5 — No internal or private labels in reader-facing text

Step codes, plan-internal names, and prior-conversation shorthand mean
nothing to a fresh reader. Describe the thing, not its label.

- ❌ "mirrors `dump_sample_summary.rs`, the D2 tool"
- ✅ "mirrors the existing `dump_sample_summary.rs` we already use to
  cross-check the coverage summaries"

If a label is load-bearing (a real cross-reference), spell out what it is
the first time, then it may be used.

## Rule 6 — Relocate precision, don't delete it

When precision threatens clarity, **move** it, don't drop it: plain lead in
the body, exact detail in a following sentence, a code block, or a
footnoted formula. The reader gets the idea first and the rigour second.
Both survive. (This is progressive disclosure applied at paragraph scale.)

## Rule 7 — Sentence-level mechanics

The standard plain-language mechanics, because they compound with the rules
above:

- **One idea per sentence; keep sentences short.** Split a sentence that
  needs two breaths to read.
- **Active voice, present tense.** "The walker emits one record per covered
  position", not "one record per covered position is emitted by the
  walker."
- **Plain words over fancy ones.** "stops running" not "is terminated";
  "use" not "utilise"; "enough" not "sufficient". Reserve the precise
  technical term for when it is the precise technical term.
- **One term per concept — and one concept per term.** Pick a single name
  for each thing and never switch (don't drift between "window", "tile",
  and "bin" for the same object). Inconsistent terminology makes a reader
  wonder whether two names mean two things.
- **Use a list or table for a set of items.** Three or more parallel
  things (fields, options, cases) read better enumerated than buried in a
  sentence.
- **One idea per paragraph.** Start it with the point; let the rest support
  that point.

## A glossary earns its keep

If a doc carries more than a few domain terms a reader might not know, give
it a short **Vocabulary/Glossary** section, grouped so the reader can skip
what they know (e.g. genetics terms brief, statistics terms fully
explained). Two sample entries set the depth:

> **q-value.** The false-discovery-rate analogue of a p-value: flag
> everything with q ≤ 0.01 and ~1% of what you removed was real. It answers
> "of the calls I throw out, what fraction am I throwing out by mistake?"

> **Marginalising over a parameter.** Averaging the likelihood over all
> plausible values of an unknown knob (weighted by a prior) instead of
> picking the single best-fitting value — which charges a model for its
> flexibility (the Occam penalty).

## Revision pass (do this before considering writing done)

1. Read each **name**: does it say what the value *is*? (Rule 2)
2. Read each **acronym/term**: defined on first use? Not stacked? (Rule 3)
3. Read each **section opening**: plain-English "what it does" before any
   formula? (Rule 1)
4. Scan for **internal labels** (step codes, "the X tool"), **borrowed
   jargon** ("kernel"), and **symbol-only summaries**. (Rules 4, 5)
5. Where precision crowded out clarity, **restructure** rather than cut.
   (Rule 6)
6. Check **mechanics**: short active-voice sentences, plain words, one term
   per concept, lists for sets, explanation kept distinct from reference.
   (Rule 7, modes)
7. Final pass: apply the **reader test** to the whole piece.

The goal is writing the owner never has to send back with "what the hell is
that?" — clear on the first read, precise on the second.

## Sources & influences

- **Google developer documentation style guide** and **Microsoft Writing
  Style Guide** — sentence-level clarity, active voice, plain words,
  consistent terminology.
- **Diátaxis** (diataxis.fr) — the explanation/reference/how-to/tutorial
  modes and the rule not to blur them.
- **plainlanguage.gov** / general plain-language guidance — word choice,
  one idea per sentence.
- The decisive material, though, is the **project's own before/after
  examples** above; the external guides corroborate, they don't replace
  them.
