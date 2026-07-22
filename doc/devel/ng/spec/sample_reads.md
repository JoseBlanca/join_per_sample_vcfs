# ng — the sample's reads: several files, one ordered stream

*Status: design spec, 2026-07-20. **All open questions resolved** (owner, 2026-07-20): O2 duplicate =
error, O3-merge → rebuild ng-owned, O4 sample-name agreement, O5 → `src/ng/read/input/`, arms unified
by an enum. Ready for its arch doc. Split out of the former `read_ingestion.md` (owner, 2026-07-20)
together with [`alignment_file.md`](alignment_file.md): that document owns everything whose subject is
**one alignment file**, this one owns everything whose subject is **the sample**. It is the layer the
region analysers (the pileup on the generic path, the tract extractor on the STR path) pull from: they
ask for a contig interval and get every read of the sample overlapping it, coordinate-sorted, whichever
file each read lives in. Under [`ng_proposal.md`](ng_proposal.md) §1 and
[`../arch/module_layout.md`](../arch/module_layout.md). Naming: **STR** in prose, `ssr` in code. The
Code-facing companion: [`../arch/sample_reads.md`](../arch/sample_reads.md) (types & interfaces).*

---

## 1. What this is — scope and non-goals

A sample is **usually several files, sometimes one**. The dominant reason for several is that the same
sample was **sequenced in several experiments** (separate runs, libraries, often separate projects);
per-lane and per-chromosome splits are the less common causes.

That distinction drives most of this document, because per-experiment files behave differently from
per-chromosome ones:

- They **span the same coordinate range**, so the merge interleaves at essentially every position
  rather than concatenating disjoint stretches. Ties are routine, not incidental (§3.2).
- They carry **different read groups, chemistries, and possibly different aligner versions**, so "are
  these really all the same sample?" is a live question (§3.1) and "how did each file behave?" is worth
  reporting separately (§3.3).

This module does exactly three things: check what can only be checked *across* files (§3.1), merge k
ordered streams into one (§3.2), and present a single entry point whose one-file and k-file arms are
indistinguishable to the caller (§3.4). It contains **no reader logic** — everything about opening,
validating, seeking and filtering a file is [`alignment_file.md`](alignment_file.md).

### The design was chosen (owner, 2026-07-20)

- **Iterator layers, and no merge at all when there is one file.** The merge is a stage that may be
  absent, not a special case with k=1 inside it (§3.2, §3.4). The two arms are unified by an **enum,
  never a `Box<dyn Iterator>`** — dynamic dispatch would block inlining on the hottest loop in the
  module (§3.4).
- **Filtering happens below, never here.** Merging a read only to drop it is wasted work; the per-file
  chain already filters, and it does so *before* a `MappedRead` is ever built — the reader→filter seam
  passes an undecoded reused buffer, so dropped reads cost no allocation (`alignment_file.md` §3.2).
  Every read reaching this module is therefore one that survived.
- **Counts stay per input file**, never pre-summed (§3.3).
- **The merge is deliberately concrete, not generic** (§5).

### Non-goals (deliberately excluded)

- **Everything per-file** → [`alignment_file.md`](alignment_file.md): the open gate, the `@SQ`
  reconciliation, index handling, region seeking, `ReadFilter` composition, within-file order
  verification.
- **Cross-*sample* merge.** This merges the files *of one sample*. A cohort of samples is a set of
  these streams, combined by whatever drives the cohort — not this module (§6, and the note in §5).
- **Filtering policy** — step 1's, and applied one layer down (`read_filtering.md`).

---

## 2. Where it sits, and what it assumes

```
file 1 ─▶ AlignmentFile ─▶ region source ─▶ ReadFilter ─▶ order-verify ─┐
file 2 ─▶ AlignmentFile ─▶ region source ─▶ ReadFilter ─▶ order-verify ─┤
   ⋮        (alignment_file.md — one validated handle each)             ▼
                                              §3.2 argmin k-way merge  ─▶ ONE ordered
                                              (+ same-file-twice check)   MappedRead stream
                                                                          per region
                                        ┌──────────────────────────────────────┤
                                        ▼                                      ▼
                               generic: pileup walker              STR: tract extractor
                               (read_preparation_generic)          (read_preparation_ssr)
```

**Two preconditions this module relies on and does not re-check** — both established by
`alignment_file.md` §3.1/§3.2, stated here because the merge is only correct because of them:

1. **`ref_id == ContigId` for every read of every file.** The per-file open gate proves each file's
   `@SQ` list *is* the reference's contig table, in order. This is what makes it sound to compare
   `(ref_id, pos)` *across* files: a contig index means the same thing everywhere. Without it the merge
   would silently interleave reads from different contigs — the permutation hole reference_info.md §1
   named.
2. **Each input stream is already coordinate-monotonic.** The per-file order-verify adapter is the
   single authority for read order in ng; the merge trusts it rather than re-checking. (Production
   checks in the merge because its readers do not; ng flips the check to where the file is read. A
   `debug_assert` may keep production's guard as cheap belt-and-braces.)

The output contract is a stream of `Result<MappedRead, IngestError>` in coordinate order — the same
`MappedRead` step 1 already yields and step 2 already consumes, so nothing downstream changes shape.

`MappedRead.source_file_index` is preserved through the merge, and it is **downstream-meaningful, not
just provenance**: because the usual reason for several files is several experiments, the file index is
the sample's *batch* label — the grouping a per-batch error model would key on (the STR parameter
pre-pass freezes ε per sample-group for exactly this reason). This module guarantees the tag survives;
what downstream does with it is not its concern.

---

## 3. The concerns

### 3.1 The cross-file check: one sample name

Each file's own `@RG SM` is validated at its open (`alignment_file.md` §3.1, check 4). What cannot live
there, because it is not a property of a single file, is **agreement**: the k files must name the *same*
sample. That check runs in `SampleReads::open` immediately after the k per-file opens — still before any
read flows. Disagreement → `SampleNameMismatch { files, names }`.

This guard earns its keep precisely because of the several-experiments reality (§1). When the files are
per-lane splits of one run, pulling in a foreign file is unlikely. When they are separate experiments
gathered by hand from different projects, **grabbing the wrong file is the realistic failure mode**, and
the `SM` agreement is the only thing that catches it. (Resolves the pre-split open question O4: yes.)

### 3.2 The merge

Merge the sample's per-file streams into one coordinate-ordered stream. This is production's
`SegmentMergedReads` (`segment_merge.rs`), re-expressed for ng:

- **Argmin k-way merge.** One head per file; each step emits the head with the smallest `(ref_id,
  pos)`, ties broken to the lowest file index (deterministic). A linear argmin over a handful of files,
  not a heap: for the k values that occur here (a few files) a linear scan of a small contiguous array
  beats a binary heap's `O(log k)` on constants and locality, and production reached the same
  conclusion. Because the files usually come from different experiments and so cover the same
  coordinate range (§1), **ties are routine, not incidental** — the tie-break is what makes output
  order reproducible, and is tested as such (T6).
- **The same-file-twice check runs only on a tie, which is what makes it free.** Two reads at different
  positions cannot be the same read, so the check has nothing to do unless the argmin already found
  ≥2 heads sharing `(ref_id, pos)` — and it found that out as a by-product of the comparison it had to
  make anyway. Within a tie, compare cheapest-first: `flag` (a `u16`) before `qname` (a short string),
  so the string compare happens only for reads that agree on position *and* flags.
  Deduplication proper is a *within-file* concern — the `0x400` flag set by the aligner/markdup is a
  per-read predicate, so it is a filter and runs **below** the merge, in step 1. Across files there is
  no such thing as a legitimate duplicate: reads from different experiments are different reads. So
  this comparison is not deduplication at all, it is an **input-sanity check**. A match on
  `(qname, flag, ref_id, pos)` means the caller passed **the same file twice** (or two files with
  overlapping content) → **hard error** (`DuplicateReadAcrossFiles`), never a silent drop. (Resolves
  the pre-split open question O2: error. Matches production.)
  *And it fires almost immediately.* The fault it detects is duplicated **input**, not a duplicated
  read, so every read collides — the very first overlapping position trips it. The check therefore does
  not need to be exhaustive across the stream to be reliable, which is the second reason it can be this
  cheap.
- **No drops of its own.** The merge filters nothing; its rejections are *errors*, not filtered reads.

The merge is fused: the first `Err` is yielded once (`Some(Err(_))`) and then `None`, so
`for read in merged { let read = read?; … }` surfaces it and cannot mistake it for clean EOF (the
`ReadFilter`/`SegmentMergedReads` convention).

**Every read of every sample passes through here, so the per-read budget is what matters.** The
comparisons are not the thing to worry about — k tuple compares per read, on small integers. Two other
costs are, and the implementation is constrained to avoid both:

- **Compare keys, not reads.** Hold the k head keys as a compact `[(ContigId, Position); k]`-shaped
  array *beside* the head slots, so the argmin scans a few contiguous integers in L1 rather than
  chasing a pointer into each `MappedRead` to reach its coordinates. The keys are copied once when a
  head is refilled, not re-read per comparison.
- **Move each `MappedRead` exactly once, and never clone it.** A `MappedRead` owns its sequence,
  qualities and CIGAR; it is the big object in this pipeline. The head slots are `Option<MappedRead>`,
  emitted with `Option::take` and refilled from the source — one move per read, no copy, no allocation.
  A `Peekable`-based design that hands out `&Result<MappedRead, _>` and then moves out of the peek slot
  is easy to write and quietly does the same work twice; and any `clone()` in this loop is a bug, not
  a slow path.

The target the implementation should be held to: **the merge adds O(k) integer comparisons and one
move per read, and allocates nothing per read.** Anything beyond that is a defect at this volume. A
benchmark over a synthetic multi-file stream belongs with the implementation (§7).

### 3.3 Counts stay per input file

The `ReadFilterCounts` (`read_filtering.md` §4) are reported as a vector indexed in parallel with
`source_file_index`, **not summed**. Drop rates are a per-experiment property — one bad run shows up as
a single file with an anomalous MAPQ or mismatch drop rate, and summing erases exactly that signal. A
caller that wants a total can add them up; this module refuses to decide that for them.

### 3.4 The entry point, and the missing merge

One type per sample, built once, then queried per region:

```
SampleReads::open(files, reference_info, filter_config, index_policy) -> Result<SampleReads, IngestError>
SampleReads::reads_in_region(&self, contig, start, end) -> Result<SampleRegionReads<'_>, IngestError>
     where SampleRegionReads: Iterator<Item = Result<MappedRead, IngestError>>   // ordered
SampleReads::counts() -> &[ReadFilterCounts]     // per file, §3.3
SampleReads::sq_md5s() -> &[&[Option<Md5>]]      // per file, forwarded for the deferred check below
```

`open` opens each file through `AlignmentFile::open` (which validates it), then checks the k files agree
on one sample name (§3.1). `reads_in_region` asks each `AlignmentFile` for its region chain and — **only
when there is more than one file** — wraps them in the merge.

**The deferred assembly check is forwarded, not performed here.** Each file's `@SQ M5` tags are
captured at its open and cross-checked against the reference's true per-contig digests once the caller
joins `reference_info`'s background verification, before output is committed
([`alignment_file.md`](alignment_file.md) §3.1). This module only makes the per-file tags reachable so
the caller can run that check for all k files at one point; it neither performs the comparison nor
merges the tags, since a disagreement must name *which file* is aligned to the wrong assembly.

**Both arms must present one type to the caller, and that type is an enum — decided (owner,
2026-07-20), not a `Box<dyn Iterator>`.** `SampleRegionReads` is `Single(..) | Merged(..)`, with
`Iterator::next` matching on the variant. The cost of the alternative is not mainly the indirect call
itself but that it is **opaque to the optimiser**: a `Box<dyn Iterator>` cannot be inlined into the
consumer's loop, so the whole per-read chain stops being optimised across the boundary — on a hot path
carrying millions of reads through ~10⁶ region queries. The enum's `match` is a predictable branch that
inlines away, and it is no more code to write. It lives with the entry point
rather than inside the merge module: the single-file arm is the *absence* of a merge, not a variant of
one. T11 asserts the caller cannot tell the arms apart.

---

## 4. Error model

The `read_filtering.md`/`reference_info.md` house pattern, same as one layer down:

- **A read is filtered** — tallied per file (§3.3), never an error, never silent.
- **A run-level failure is fatal** — an `Err`, never swallowed into a drop. *At `open`*: anything
  `AlignmentFile::open` rejects for any of the k files, plus this module's own
  `SampleNameMismatch` — all before any read flows. *Mid-stream*: anything an input stream yields, plus
  `DuplicateReadAcrossFiles`; yielded once, then the iterator fuses.

This module's own variants are `SampleNameMismatch { files, names }` and `DuplicateReadAcrossFiles`;
everything else arrives wrapped from `AlignmentFileError`. Whether `IngestError` wraps that enum or
flattens it → arch.

---

## 5. The merge is deliberately concrete (and when to revisit)

A generic k-way merger — "merge any k sorted iterators" — is tempting and was considered. **Not built,
on purpose.** The reasoning, recorded so it is not re-litigated from scratch:

- **The one plausible second caller may not want this operation at all.** The cohort layer (§6) also
  combines k sorted streams, but it likely wants to **group** — at position p, hand me all N samples'
  reads at p — not to **interleave** one read at a time. Those are different functions with different
  signatures, and a generic interleave does not give you a group. Building the abstraction now means
  guessing which of the two is the shared thing, and guessing wrong leaves an abstraction one caller
  uses and the other works around.
- **The same-file-twice check does not generalise.** It compares `(qname, flag, ref_id, pos)`, which
  only means something for reads. A generic merger would need it as a hook, and that hook makes the
  generic version *harder to read* than the concrete one for the single caller that exists today.

**If it is generalised later, two constraints:**

- **Do not put a `Sortable`-style trait on the item.** The right shape is `Iterator<Item = T>` plus a
  key extractor `Fn(&T) -> K where K: Ord`. A trait forces `MappedRead` to have one canonical ordering,
  which is wrong — it is ordered by `(ref_id, pos)` here, but another caller might want `(ref_id, end)`
  or `(qname)`. A key function costs nothing and keeps the ordering decision at the call site.
- **Keep the concrete implementation shaped so the refactor is mechanical.** Key extraction goes in one
  named place, and the same-file-twice check stays visibly separate from the argmin. Rust lets a type
  parameter be added to a working struct without redesigning it, so the generalisation is cheap *at the
  moment it is justified* — which is the second real use, not the first anticipated one.

**Where the question gets reopened:** the cohort layer (§6). If it turns out to want `merge` rather than
`group_by_position`, generalise then, with two callers in hand and this module's tests already written.

---

## 6. Module layout, and the reuse question

Home (owner, 2026-07-20): the `read/input/` submodule (shared with `alignment_file.md` §6, which owns
`open_bam.rs` and `region_query.rs`). This spec owns:

```
src/ng/read/input/
├── mod.rs    – SampleReads::{open, reads_in_region, counts, sq_md5s}; SampleRegionReads
│               (Single|Merged); the cross-file sample-name check; the IngestError enum; re-exports
└── merge.rs  – the argmin k-way merge + the same-file-twice check (§3.2)
```

**Decided (owner, 2026-07-20): rebuild ng-owned.** Production's `segment_merge::SegmentMergedReads` is
`pub(crate)` (ng-reachable) and already does argmin + cross-file duplicate detection over the same
`MappedRead`, so reuse was a live option. It loses on four counts, three of them behavioural:

- **It re-checks per-file order**, which is redundant here — ng moved that check down to where the file
  is read (§2, precondition 2), so production's guard would duplicate work on every read.
- **It has no merge-free single-file arm** (§3.4); a one-file sample would pay merge machinery for
  nothing.
- **It reports counts in production's shape**, not the per-file vector this module requires (§3.3).
- **Its head handling is unaudited against the per-read budget** (§3.2). Adopting it would mean
  measuring it anyway — and if it fell short, patching production, which the ng freeze forbids.

Against that, the merge is small: an argmin over a few keys, a tie-break, and one comparison. Rebuilding
is cheaper than adapting.

Together with `alignment_file.md`'s O3-file (also rebuild), **ng now owns the whole read path over
noodles, reusing only pure helpers** — the `read_filtering` precedent applied end to end. The two were
decided separately on their own reasoning, which is what splitting the specs was for; the shared
outcome is a result, not a package deal.

**Deferred — cross-sample cohort assembly.** Whatever drives the cohort holds N `SampleReads`; not here.
That layer is also where §5's generic-merger question gets reopened.

---

## 7. Test obligations (T-list — the acceptance bar)

- **T6 — two-file merge** interleaves in coordinate order with correct `source_file_index` tags; **reads
  at the same position in both files break to the lower file index**, and the whole output is
  byte-identical across repeated runs (the tie-break is load-bearing because same-experiment-range files
  tie constantly, §3.2); `ReadFilterCounts` are reported **per file**, each equal to that file's own
  drop tally.
- **T7 — the same file passed twice → `DuplicateReadAcrossFiles`** (fused), raised at the **first**
  colliding position rather than after draining the region (§3.2); **distinct reads at one locus are
  not duplicates** (both kept), including reads that share `(ref_id, pos)` and differ only in `qname`,
  which is the case the cheap-first comparison order must not get wrong.
- **T14 — the merge's per-read cost.** A benchmark over a synthetic k-file stream, asserting the
  budget §3.2 sets: no per-read allocation, and no `MappedRead` clone in the merge loop. Cheap to
  check with a counting allocator or `dhat`; worth having because the failure is silent — a `clone()`
  or a redundant move costs throughput without changing a single output byte.
- **T11 — the single-file arm has no merge and matches the merged one.** A one-file sample yields the
  same reads as the same file passed as one of two inputs (the other empty over the region); asserted
  through the same `SampleRegionReads` type, so the caller cannot tell the arms apart.
- **T12b — two files naming *different* samples → `SampleNameMismatch` at `open`**, before any read
  (§3.1). (The single-file half is `alignment_file.md` T12a.)

T-numbers are kept from the pre-split spec so review notes stay traceable; T1–T5, T8–T10, T12a and T13
live in `alignment_file.md`.

Fixtures: tiny in-memory BAM/CRAM built with noodles writers + an in-crate CSI/CRAI build (production's
`segment_merge` test module is the template), over a small multi-contig reference indexed by
`reference_info`.

---

## 8. Deferred, with a home

- **Cross-sample cohort assembly** → whatever drives the cohort holds N `SampleReads` (§6); also where
  the generic-merger question reopens (§5).
- **Mixed BAM+CRAM in one sample** → production rejects it (`preflight` mixed-format error). ng inherits
  that unless a reason to allow it appears.
- **A generic k-way merger** → §5, on the second real use.

---

## 9. Open questions

### Still open

*None. Both specs are fully decided; what remains is arch-level shape (exact types, error nesting,
whether `IngestError` wraps or flattens `AlignmentFileError`).*

### Resolved (owner, 2026-07-20)

- **O3-merge — rebuild vs reuse the merge** → **rebuild ng-owned** (§6): production's merge re-checks
  order redundantly, has no single-file arm, reports counts in the wrong shape, and is unaudited
  against the per-read budget.
- **O5 — module home** → `src/ng/read/input/`, shared with `alignment_file.md` (§6).
- **The two arms are unified by an enum, not `Box<dyn Iterator>`** (§3.4).
- **O2 — cross-file identical read** → hard error, not a dedup; checked only on an argmin tie (§3.2).
- **O4 — sample-name check at `open`** → yes, extended to cross-file agreement (§3.1).
