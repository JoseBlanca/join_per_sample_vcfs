## Design principles

A few commitments shape both this specification and the code that will
implement it. They take precedence when a concrete design decision is
ambiguous.

1. **Clarity and readability are paramount.** They apply equally to
   this specification and to the code that implements it — not one more
   than the other. A reader new to the project should be able to
   understand both without needing to consult additional sources.
   Cleverness that hurts readability is never a good trade.
   - *In the spec:* where a choice is not obvious, the reason is spelled
     out; where a trade-off exists, both sides are named.
   - *In the code:* identifiers name what they are; control flow reads
     top to bottom; abstractions are introduced only when they earn
     their complexity; decisions that can't be understood from the
     names alone get a short comment explaining *why* (never *what* —
     the code says that already).
   - *No magic numbers, as a concrete example of readability.* A bare
     `20` or `30` mid-function tells a reader nothing; the name of a
     `const` is where the intent lives. Every numeric (or
     non-trivial string) literal whose meaning is not obvious from
     immediate context lives in a named `const` (`pub const` if part
     of the public API, private otherwise) with a short doc comment
     identifying units and source — spec section, measurement,
     vendor recommendation. Trivial literals are exempt: `0`, `1`,
     `±1` in arithmetic and indexing, and fixed buffer sizes whose
     only meaning is the literal itself. The rule is a readability
     rule, not a numerology rule: the question to ask is "would a
     reader meeting this literal know what it means?", not "is this
     a number?".

2. **No guessing at user intent — every default and decision is
   explicit.** Every parameter default, every fallback behaviour, every
   assumption about the input is named in this document. Hidden
   heuristics and "sensible" values that only exist in the source
   accumulate into systems nobody can debug later. If a behaviour isn't
   written here, the implementation should not silently invent it.

3. **Errors must not pass silently.** Corrupt inputs, unexpected
   formats, out-of-range values, missing files — each produces a clear,
   actionable error message and halts. The default posture of every
   component is "reject and report," never "guess and proceed." It is
   far better to fail loudly at Stage 1 than to emit subtly wrong
   genotypes at the final stage.