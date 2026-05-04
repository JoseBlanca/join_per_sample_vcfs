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

4. **Push side effects to the edges of the program — when it costs
   nothing to do so.** The pure parts of the code are easier to
   reason about, easier to test, and easier to reuse: given the same
   inputs they produce the same outputs, no filesystem state, no
   network, no clock, no hidden globals. When a function can be
   structured so that I/O happens at its caller and the function
   itself only transforms data, prefer that arrangement. The most
   common shape this takes in this project is splitting a "do the
   work" function from a "set up the inputs" function: the work
   function takes already-prepared streams or buffers; a thin
   convenience wrapper opens the files and calls into it.
   - *Example:* the per-sample caller's CRAM reader exposes
     `new(paths, fasta, config)` for the convenient case, and
     `pub(crate) from_open_crams(streams, contigs, sample, config)`
     for the testable core. The merge / order / dedup / filter logic
     lives in the latter; tests inject pre-built record streams and
     skip the file-on-disk dance entirely. The public API stays
     simple; the test suite gets fast, deterministic, mock-free
     coverage of the logic that matters.
   - *Not an absolute rule.* When pushing the side effect outward
     would force a generic parameter through every caller, leak a
     dependency type into the public API, or just add ceremony for
     no real testing or reuse benefit, leave it where it is. The
     value is in the parts where it pays off without friction, not
     in mechanically applying the principle. The question to ask is
     "does the pure half become independently useful?", not "is
     there I/O happening here?".

5. **Don't leak implementation details or third-party API shapes
   into our own public APIs.** Every type from a dependency that
   appears in our public signature is a tie that's hard to undo
   later: changing the dependency or upgrading it across a breaking
   release becomes a public-API change for us, even if the behaviour
   we expose is identical. Wrap dependency types at the boundary —
   our own owned struct, our own enum, our own error — and let the
   wrapping happen once, at the place where the dependency is
   actually used. Internal modules can talk to each other in
   dependency-native types if it's cheaper that way; the rule is
   about what we expose outward.
   - *Example:* the CRAM reader yields our own `MappedRead`, not
     noodles' `RecordBuf`. The conversion happens in one place, in
     the per-record decoding helper. If we ever swap noodles for
     another CRAM library — or noodles changes its record type
     across a major release — only the conversion helper changes;
     downstream stages and external callers do not.
   - *Not an absolute rule either.* For types that are genuinely
     standard across the ecosystem (`std::path::Path`,
     `std::io::Read`, integer types) the wrapping is pointless. The
     test is "is this type effectively part of the dependency's
     identity, or is it a vocabulary term shared by everyone in the
     domain?". Wrap the former; pass the latter through.

6. **Typed errors at module boundaries; anyhow context at the program
   edge.** Principle 3 says errors must not pass silently. This one is
   about *what shape* errors take. Two layers, each with a different
   job:
   - *Inside a module / library boundary:* errors are **typed enums**
     built with `thiserror`, scoped to the operation that produces
     them (`CramInputError`, `VcfParseError`, `BaqError`). Callers
     that need to react to a specific failure can match on the
     variant; tests can assert on the exact variant rather than
     stringly-comparing messages. Each variant carries enough context
     in its fields to point at the offending input — file path,
     QNAME, record position — so the failure is debuggable from the
     error alone. Catch-all `Other(String)` variants are a bad sign:
     they collapse failure modes that callers should be able to
     distinguish.
   - *At the orchestrator / CLI edge:* errors flow through
     `anyhow::Result` with `.with_context(|| format!(...))` added at
     each meaningful boundary. The context names *what was being
     attempted* and *which input* — `opening CRAM {path}`, `validating
     header against {fasta}`, `merging position {chrom}:{pos}`. The
     final error a user sees is a chain that reads top-down like a
     story: "failed to open CRAM /a.cram: contig list mismatch with
     /b.cram: name disagreement at index 3 (chrA vs chrB)". `anyhow`
     handles the chain plumbing; we provide the prose.
   - *Connecting them.* Typed errors implement `std::error::Error` and
     bubble up into `anyhow::Result` cleanly via `?`. Source chains
     are preserved with `#[source]` / `#[from]` on the typed enum, so
     the underlying cause survives all the way to the top.
   - *Not an absolute rule.* For one-shot binaries with no internal
     library shape, `anyhow` everywhere is fine. For a tiny module
     where one error variant is plausible, a single `thiserror` enum
     with one variant is overkill — `io::Error` directly is fine.
     The split is a guideline for places where both layers genuinely
     exist.

7. **Type names describe data shapes; variable and field names
   might describe roles.** A type alias or struct should name *what the
   value is*, not what it's used for. Roles could live on the binding
   (variable, field, argument), where surrounding context already
   constrains their meaning. The same data shape often plays
   different roles in different parts of the code; a shape-named
   type composes naturally with `Option` / `Vec` / `HashMap` to
   express presence, multiplicity, or keying at the storage site,
   while a role-named type stops being usable the moment a second
   use site appears.
   - *Example:* the CRAM merge's `type Locus = (usize, u64)` —
     contig-list index plus 1-based position — appears in three
     different roles: `prev_per_file: Vec<Option<Locus>>` (per-file
     order tracker), `window_anchor: Option<Locus>` (duplicate-
     window anchor), and as the inner half of `argmin_head`'s
     scratch `Option<(usize, Locus)>` (best stream index plus its
     locus). The same shape gets a single name; each role is named
     where the value is bound, and the surrounding `Option` /
     `Vec` make storage semantics explicit. An earlier iteration
     used `type PerFileOrder = Option<(usize, u64)>`, which baked
     the role into the type name — and as a result the same shape
     went unnamed in `window_anchor` and `argmin_head`, where the
     role was different.
   - *Newtype wrappers are a different tool.* `struct UserId(u64)`
     exists precisely to *create* a type-level role distinction
     between values that share a shape, and that is the whole point
     of reaching for one. The principle here is about type
     *aliases* and *owned domain types*; for those, naming the
     shape keeps the door open for reuse where naming the role
     closes it.
   - *Not an absolute rule.* When a data shape has exactly one
     plausible role across the whole module and you can't picture
     it appearing elsewhere, a role-named alias is fine — there is
     no reuse cost to pay. The test is "would this shape have an
     independent meaning if it appeared somewhere else?", not "is
     this a tuple?".