# Indel normalization (left-alignment)

Implementation plan for the Stage 1 indel-normalization slice of the
calling pipeline. Specified in
[calling_pipeline_architecture.md §"Indel normalization (left-alignment)"](../specs/calling_pipeline_architecture.md#L224)
(commit `8d7dd17`), with supporting updates to the Stage 1 conceptual
list, the Stage 2 indel-anchoring convention, and the Stage 3 DUST
rationale.

**Problem.** A short-read aligner places an indel inside a tandem
repeat or homopolymer at an arbitrary offset among equally-scoring
positions. Two reads carrying the *same* biological indel present it at
different CIGAR offsets, so recorded verbatim they become different
`(anchor, REF, ALT)` triples. The Stage 5 merge unifies alleles by
exact sequence match, so these near-duplicates bucket as distinct
alleles, support fragments, and the genotyper drops the site. This is
the dominant cause of low indel recall on the HG002 benchmark, and DUST
does not fix it (every indel in a region that *passes* the complexity
threshold still needs its evidence consolidated).

**Fix.** Left-align each indel against the reference in Stage 1 before
it is recorded as an allele observation — the canonical
`bcftools norm` / `vt normalize` / GATK form, matching GIAB truth.
Identical indels then consolidate onto one `(anchor, REF, ALT)` and
support sums.

## What we learned from freebayes and GATK

Both reference callers in the tree solve exactly this, and both do it
the **same way**: they left-align by **rewriting the read's whole CIGAR**
once, not by shifting indels in isolation. This plan adopts their
approach rather than the per-op scheme an earlier draft proposed.

- **freebayes — left-aligns during calling, per read, on by default.**
  [`AlleleParser.cpp:2055-2058`](../../freebayes/src/AlleleParser.cpp#L2055)
  calls `stablyLeftAlign(currentAlignment, refWindow)` right after the
  MAPQ filter and **before** parsing alleles, where `refWindow` is
  exactly the read's own footprint
  `currentSequence.substr(start, end - start + 1)`. The core
  ([`LeftAlign.cpp`](../../freebayes/src/LeftAlign.cpp)) walks the CIGAR,
  records every indel, then for each indel does two shifts:
  (1) a **repeat-unit shift** that only tries shift distances dividing
  the indel length (`do { ++i } while (length % i != 0)` — handles
  multi-base SSR motifs); (2) a **single-base exchange/rotation** loop
  for homopolymers (`seq = seq.last + seq[..last]` — the right-rotation).
  Crucially the shift condition requires the moved bases to match
  **both the reference and the read's own bases**
  ([LeftAlign.cpp:126-133](../../freebayes/src/LeftAlign.cpp#L126),
  [:161-169](../../freebayes/src/LeftAlign.cpp#L161)), and stops at the
  previous indel (collision → merge). It then merges adjacent same-class
  indels and **iterates `leftAlign` until the CIGAR is stable**
  (`stablyLeftAlign`, `maxiterations = 50`).
  `bamleftalign` is the same logic as a standalone BAM-rewrite tool.
- **GATK — same shape, single right-to-left pass.**
  [`AlignmentUtils.leftAlignIndels`](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/read/AlignmentUtils.java#L731)
  early-returns when there is no indel, then traverses the CIGAR
  **right to left**, accumulating each indel's ref/read ranges and
  left-aligning it into the preceding **alignment block only** (a
  `maxShift` bounded by that block prevents colliding into the previous
  indel — collisions merge). The shift core
  [`normalizeAlleles`](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/read/AlignmentUtils.java#L818)
  **trims first** (shared right base, then shared left base for
  parsimony), then shifts left while `nextBaseOnLeftIsSame` **across all
  sequences (ref + read)** and `lastBaseOnRightIsSame`. The read-level
  driver
  [`LeftAlignIndels.apply`](../../gatk/src/main/java/org/broadinstitute/hellbender/tools/LeftAlignIndels.java#L74)
  skips unmapped reads and reads with `≤ 1` CIGAR element, and — the one
  edge case worth importing — when a deletion left-aligns all the way to
  the read start it is **removed and the read's start position is bumped
  right** (`leadingDeletionBasesRemoved`).

Distilled, the learnings that change the plan:

1. **Rewrite the CIGAR, don't shift per-op.** Producing a left-aligned
   CIGAR and rebuilding the cursor from it means the cursor's emission,
   offsets, and binary-search pruning need **no change at all** — the
   earlier draft's per-op `left_shift` field, insertion-rotation logic,
   and (most importantly) the binary-search-pruning hazard all
   **disappear**, because the rewritten CIGAR's `ref_pos` values already
   sit at the normalized positions and stay monotonic.
2. **Check the read's bases, not just the reference.** Both callers
   require the shifted bases to agree on ref *and* read. This is
   strictly more conservative than pure reference-only `bcftools norm`:
   it won't shift an indel past a position where the read disagrees with
   the reference. In the consolidation case (reads matching the
   reference through the repeat) it converges to the same canonical
   `(anchor, REF, ALT)` as bcftools — so it still matches truth — while
   refusing to over-shift through a read-side SNP.
3. **Trim for parsimony first**, then roll left. Matches the VCF
   convention and GATK's `normalizeAlleles`.
4. **Bound the shift by the preceding alignment block / read footprint**
   and merge indels that collide. This is the spec's "overlapping
   events" handling, grounded in both impls; the freebayes refWindow
   confirms the spec's "roll stops at the read's footprint" bound.
5. **Iterate to stability** (freebayes) or do GATK's single right-to-left
   pass; both reach the same fixed point. We port GATK's right-to-left
   pass — it is the more elegant model for a Rust port and is O(cigar).
6. **Skip trivial reads** (no indel, or `≤ 1` CIGAR op) on a fast path.
7. **Handle a leading deletion** by dropping it and bumping
   `alignment_start` (GATK), reconciled with our existing first/last-op
   indel rejection.

## Decision: ref+read dual check (Option A)

The shift core requires the moved base to match **both the reference and
the read** (what freebayes and GATK do), rather than reference-only
(`bcftools norm`). The two agree whenever the read matches the reference
through the repeat — i.e. on every clean repeat-region indel, which is
the consolidation case — and diverge only when a read carries a mismatch
inside the repeat. There, Option A stays faithful to the read (never
records an indel the read's own bases contradict, upholding the `.psp`
faithfulness principle) and refuses to roll over a sequencing error;
reference-only would not. Revisiting reference-only is a **data-driven
follow-up gated on the HG002 numbers** (only relevant if reads with
frequent in-repeat sequencing errors are shown to lose recall), not a
now decision.

## Scope

In:

- A pure `left_align_cigar` routine — a Rust port of GATK's
  `leftAlignIndels` + `normalizeAlleles` (cross-checked against
  freebayes' `LeftAlign.cpp`) — taking the CIGAR, the read footprint's
  reference bases, the read bases, and the read start, returning a
  normalized `Vec<CigarOp>` plus any `leading_deletion_bases_removed`.
- Wiring that, **at read-admit time** (where the `MultiChromRefFetcher`
  is reachable), fetches the read's footprint reference, rewrites
  `read.cigar` in place, and bumps `alignment_start` if a leading
  deletion was removed — then builds the `CigarCursor` from the
  normalized CIGAR (cursor code unchanged).
- A fast path that skips reads with no indel op or `≤ 1` CIGAR element.
- A debug-only invariant check that left-alignment preserved the read's
  mismatch count against the reference (freebayes' `countMismatches`
  guard) — catches a buggy port immediately.
- Unit tests for the pure routine against canonical `bcftools`/`vt`
  cases: homopolymer run, dinucleotide SSR, multi-base motif,
  deletion vs. insertion, already-canonical no-op, multi-indel
  merge/collision, leading-deletion removal, roll-stops-at-footprint,
  and a read-side-SNP-blocks-shift case (the ref+read distinction).
- Cursor-vs-`decompose`-oracle parity preserved automatically (both
  consume the same rewritten CIGAR); add fixtures with normalized input.
- A walker-level test on a constructed repeat fixture showing two reads
  carrying the same biological indel at different CIGAR offsets collapse
  to one record with summed support.

Out:

- Re-deriving the algorithm. We port the established trim+left-roll from
  the two in-tree references; we do not invent a variant.
- Right-alignment / `--keep-left`-style alternatives. Left is the
  GIAB/truth convention.
- Re-tuning the DUST threshold (the spec notes it becomes less
  load-bearing; revisiting its default is a separate follow-up —
  Stage 3 §"Parameters and opt-out").
- Any change to the Stage 5 merge. The `.psp` already carries canonical
  coordinates, so the existing exact-match unification consolidates them
  unchanged.
- Multi-allelic decomposition / MNP splitting.

## Why in read-prep, before the cursor is built

Left-alignment moves an indel's anchor **leftward**. The walker collects
an event only when `walker_pos == anchor`
([`process_position`](../../src/pileup/walker/driver.rs#L374) →
[`events_at`](../../src/pileup/walker/cigar_cursor.rs#L478)). A
post-collection rewrite would point at a position the walker already
passed — so normalization must happen **before the cursor is built**.
Rewriting the `PreparedRead`'s `cigar` in `BaqEngine::process` (the
per-read prep stage, upstream of the walker) achieves that and keeps the
walker free of any normalization concern. This mirrors freebayes doing
it before allele parsing.

### Timing is safe (no event is lost)

The walker still sees every normalized anchor because:

1. [`advance`](../../src/pileup/walker/driver.rs#L499) steps strictly
   `+1` whenever the active set is non-empty; it only jumps over gaps
   when *no* read is active
   ([driver.rs:519](../../src/pileup/walker/driver.rs#L519)). So while a
   read is alive the walker visits **every** position it spans.
2. Left-alignment never moves an indel left of the read's own footprint
   (the shift is bounded by the preceding alignment block), and
   `remove_deletions_at_ends = false` leaves `alignment_start` fixed, so
   `alignment_start ≤ normalized_anchor ≤ alignment_end` and the read
   stream stays coordinate-sorted.

The read is therefore active across the whole interval containing the
normalized anchor, and the walker steps one-by-one across it — it is
guaranteed to land on the normalized anchor with the read active.

## The algorithm (pure routine)

Port of GATK `leftAlignIndels` (`AlignmentUtils.java`), cross-checked
against freebayes `LeftAlign.cpp`:

1. **Early return** if the CIGAR has no indel op.
2. **Right-to-left pass.** Carry each indel's ref-range and read-range;
   when the preceding op is an alignment block, left-align the indel
   into it via the shift core, bounded by `maxShift = block length`
   (so it cannot cross the previous indel — a collision merges them).
3. **Shift core** (`normalizeAlleles`): trim the shared right base, then
   the shared left base (parsimony); then while
   `start_shift < maxShift` and the next base on the left is equal
   **across reference and read** and the last base on the right is equal,
   shift left by one.
4. **Rebuild** the CIGAR from the shifted ranges, emitting the new
   match lengths the shifts imply, and merging adjacent same-class
   indels. A deletion that reaches the read start is dropped and
   reported via `leading_deletion_bases_removed`.

Because the read bases never move — only CIGAR op boundaries do — the
cursor extracting an insertion's `seq` from `read.seq[read_off..]`
against the **new** CIGAR automatically yields the correctly-rotated
inserted sequence. No separate rotation step is needed.

## Wiring

- Rewrite happens in
  [`BaqEngine::process`](../../src/pileup/per_sample/baq_engine.rs) — the
  per-read transform that `baq_stream` already runs in parallel — via the
  shared [`indel_norm::left_align_prepared`](../../src/pileup/walker/indel_norm.rs)
  helper, after the HMM and `mapped_to_prepared`. The cursor is later
  built from the already-normalized CIGAR in the walker, unchanged.
- Reference window: **reused from BAQ**, which already fetches a window
  `[xb, xe]` wider than the read footprint. No extra fetch; the read's
  first aligned base is at index `pos_0 - xb`. A raw-ASCII copy of the
  window (`BaqEngine::raw_ref`) is kept only for indel-bearing reads.
- Fast path: skip the raw-window copy + rewrite entirely when the CIGAR
  has no indel op — the overwhelming majority of reads.

## Interactions to confirm

- **BAQ ordering.** BAQ is computed upstream (per-sample pileup engine)
  and produces a per-read-base `bq_baq` array. Because left-alignment
  moves only CIGAR boundaries and never the read bases, that array stays
  index-valid after the rewrite. The indel BQ-proxy window is in *read*
  coordinates ([`indel_bq_proxy_*`](../../src/pileup/walker/decompose.rs#L181)),
  so it is unaffected and is automatically "centred on the normalized
  anchor" as the spec requires. Plan keeps BAQ before normalization;
  flag for review that we are not re-running BAQ on the normalized CIGAR.
- **First/last-op indel rejection** and the **`anchor < 1` guard** in
  the cursor ([cigar_cursor.rs:422](../../src/pileup/walker/cigar_cursor.rs#L422),
  [:441](../../src/pileup/walker/cigar_cursor.rs#L441)) still apply to
  the normalized CIGAR. Reconcile with GATK's leading-deletion removal
  so the two rules agree on indels that reach a read end.
- **Phase-chain ids** are untouched (membership, not coordinates).
- **Overlapping-event records** (Stage 2 §"Overlapping events extend the
  anchor REF") fall out of the CIGAR rewrite + the existing record
  formation, since collisions merge during normalization.

## Phasing (incremental)

1. **Pure routine + tests.** ✅ Done (commit `db4acc9`). Ported
   `left_align_cigar` / `normalize_alleles` + a CIGAR builder into
   `src/pileup/walker/indel_norm.rs`; 10 unit tests over the canonical
   cases. No walker wiring.
2. **Read-prep wiring (in the parallel BAQ stage).** ✅ Done.
   Normalization runs in `BaqEngine::process`
   ([baq_engine.rs](../../src/pileup/per_sample/baq_engine.rs)) — the
   existing per-read transform that `baq_stream` already fans out over a
   rayon `par_iter`. It **reuses the reference window BAQ fetches** (no
   second fetch; the read's first aligned base sits at `pos_0 - xb`
   within it) and rewrites the `PreparedRead`'s CIGAR via the shared
   `indel_norm::left_align_prepared` helper, with the debug
   mismatch-invariant. `remove_deletions_at_ends = false` keeps
   `alignment_start` fixed so the chunk stays coordinate-sorted; the
   cursor rejects any resulting first/last-op deletion.

   **Why here, not the walker.** BAQ is placement-invariant (its HMM
   realigns the read against the window regardless of input indel
   placement), so order doesn't matter for correctness — which freed us
   to put normalization in the *parallel* prep stage rather than the
   serial walker `admit_read`. Net wins: free parallelism, zero extra
   fetch, pure walker. (An earlier draft wired it into `admit_read`; that
   was reverted.)

   **Normalization is mandatory, including under `--no-baq`.**
   `BaqEngine::process` takes an `apply_baq: bool` that gates **only** the
   BAQ HMM: every read is left-aligned either way (`bq_baq` is the raw
   qual when BAQ is off). Both the `pileup`/`var-calling` path
   (`stage1_pipeline`) and `var-calling-from-bam` now route *both* modes
   through `BaqStream` with `apply_baq = !no_baq`; the separate
   `prepare_passthrough` no-baq map (which skipped normalization) is
   deleted. `baq_skip_counts` is still reported as `None` ("baq:
   disabled") under `--no-baq` to preserve that signal. So `--no-baq`
   gives exactly no-BAQ — nothing more, nothing less — and still
   normalizes.

   Tests: engine-level normalization tests in `baq_tests.rs` (deletion +
   insertion left-align to the leftmost; two reads at different offsets
   normalize to one anchor — the consolidation proof at the prep layer)
   and a non-stripping unit test in `indel_norm`. One pre-existing
   mixed-CIGAR BAQ fixture had a latent off-by-one (CIGAR consumed 15
   read bases, seq was 16) that the normalizer's invariant correctly
   surfaced — fixed to be consistent.
3. **Integration + benchmark.** Review the deliberate output diffs on
   real pileup/cohort data (normalization *changes* output where repeats
   exist — not a byte-identical gate; the synthetic integration fixtures
   showed no drift); HG002 indel-recall re-measurement is the downstream
   acceptance signal.

Pause between phases per the project's incremental-steps convention.

## Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features
  -D warnings`, `cargo test`.
- New unit tests (pure routine), cursor/oracle parity, walker collapse
  test.
- HG002 benchmark indel-recall re-measurement is the real acceptance
  signal but lives downstream of this slice (benchmarks harness); noted
  as the headline follow-up rather than gated here.
