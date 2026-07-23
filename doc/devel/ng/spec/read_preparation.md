# ng — read preparation

*Status: design spec, rewritten 2026-07-23, replacing the 2026-07-14 draft of this file. Defines
what read preparation does, how much work it gives each read, and how it composes with the
alignment module ([`alignment.md`](alignment.md)). The two path specs —
[`read_preparation_generic.md`](read_preparation_generic.md) and
[`read_preparation_ssr.md`](read_preparation_ssr.md) — still hold for their path-specific detail,
except where §8 records a change. **No code yet.** Naming: **STR** in prose, `ssr` in code.*

---

## 1. What read preparation is

Read filtering decided *which* reads survive. Read preparation decides **what each surviving read
looks like as evidence at its locus**, and produces it: a per-read transform from a filtered read
into a locus-ready observation the evidence-gatherer downstream can use.

Its distinctive job is **deciding how much work a read needs**. Most reads need almost none — the
mapper placed them well and there is nothing to fix. Some need their indels rewritten into a
canonical form. A few need their alignment redone from scratch because the mapper's answer is not
trustworthy. Preparation makes that choice per read and carries it out.

**It does not compute alignments itself.** Lining a read up against the reference belongs to the
alignment module, which knows nothing about caller steps. Preparation *composes* that module: it
picks an algorithm, supplies the inputs, and interprets the result. Keeping the boundary there is
what lets an alignment algorithm be swapped and measured without touching preparation, and vice
versa.

Preparation is **marker-aware**, and that is its one structural split: a locus is either generic
(SNP/indel) or a microsatellite, and the two paths share this contract but almost nothing else —
different work, different outputs, different consumers.

**Non-goals.** Preparation never *drops* a read for a whole-read property — that was filtering. It
does not decompose a read into per-position events (the pileup walker's job, and it needs
locus-column context that preparation deliberately does not have). It does not generate candidate
alleles or compute read likelihoods. And it does not reassemble: **local haplotype reassembly is out
of scope for ng, not deferred** — the production caller already calls generic loci better than GATK
without it, so it buys nothing here, and it would break the per-read independence every mode below
relies on (§6).

---

## 2. How much work a read needs — the generic path's three modes

On the generic path, preparation picks one of three modes per read. All three produce the **same
output type**, which is what makes them interchangeable.

| mode | what it does | when | cost |
|---|---|---|---|
| **pass through** | nothing to the placement | the read shows no insertions or deletions | none |
| **canonicalize** | rewrite the read's indels into their leftmost equivalent spelling; optionally cap base qualities by alignment confidence | the read has indels and its placement is trusted | cheap |
| **re-align** | discard the mapper's line-up and compute a fresh one from the read's bases | the read's placement is not trusted (§4) | a full alignment per read |

**Pass-through is a fast path, not a different answer.** Left-alignment shifts indels; a read with
no indels has nothing to shift, so canonicalizing it is provably a no-op. Recognizing that from the
read's own alignment record and skipping the work changes nothing about the result. (Whether it also
skips the base-quality capping is a separate decision — that step reads alignment *confidence*, which
is not identically neutral just because there are no indels. §9.)

**Canonicalize is about spelling, not quality.** The same insertion or deletion can be written at
several equivalent reference positions when it sits in or near a repeat — the gap slides without
changing a single base of the result. Left-alignment picks the leftmost of those spellings so that
equivalent variants get an identical one. This matters because differently-spelled equivalents look
like different variants, and the reads supporting them scatter across several weak candidates
instead of pooling into one strong one. The operation itself lives in the alignment module
([`alignment.md`](alignment.md) §6); preparation decides when to apply it.

**Re-align is the only mode that questions the mapper.** The other two accept the read's placement
and work within it; this one throws it away and computes a new line-up with a best-path alignment
algorithm ([`alignment.md`](alignment.md) §4.1). It is the expensive mode and the rare one — and it
is the only route by which a mis-placed read can be rescued, so its trigger (§4) determines how much
of that class the caller ever recovers.

---

## 3. The STR path always aligns

The microsatellite path has no equivalent of pass-through or canonicalize: **measuring the read's
repeat is the preparation**, and measuring it requires an alignment. Every read that reaches the STR
preparer is aligned against its locus with the repeat-aware best-path aligner
([`alignment.md`](alignment.md) §4.2), and the repeat is read off the two flank boundaries.

Canonicalizing instead would not help. In a repeat, shifting the mapper's indel to its leftmost
spelling does not recover how many units the read carries — which is the only quantity the STR path
wants.

**This was decided against the alternative, not by default.** Production considered a two-tier
scheme that trusted the mapper's alignment for clean-looking reads and aligned only the doubtful
ones, and dropped it in favour of aligning every spanning read. A fast path that counts units
directly from a clean read remains parked as something to measure, not as a decision already taken —
so if the cost of aligning everything ever justifies a fast tier, the comparison is a measurement
away, and §9 keeps it open.

---

## 4. Choosing the mode — the part that is not settled

Pass-through and canonicalize are chosen from the read's own alignment record: does it carry indels
or not. That decision needs nothing but the read.

**Re-align is different: it needs a judgement the read cannot make about itself.** "The mapper's
answer here is not trustworthy" is a property of the *place*, not of one read — a region where reads
disagree with each other, pile up mismatches in the same column, or clip at the same offset. Nothing
in the current ng step map produces that judgement:

- **Region typing** classifies the reference (microsatellite, repeat cluster, satellite, generic).
  It says what the reference *is*, not how well reads mapped to it — it never looks at reads.
- **The evidence-gatherer** does see the reads and could discover it, but it runs *after*
  preparation, so a verdict it produces arrives too late for the read it should have changed.

So the trigger needs either a new producer or a deliberate two-pass arrangement, and picking one is
open (§9). Recording it here rather than assuming a default matters because **the trigger, not the
algorithm, sets how much this mode is worth**: an aligner that never fires rescues nothing.

---

## 5. What preparation produces

Each path owns its output type; the two are never held interchangeably, because the router picks the
path before preparation runs and the two go to different consumers. Both are specified in their path
specs and only summarized here:

- **Generic** — a still-decomposable read: the placement (canonical or re-aligned), the bases, the
  qualities, and the per-read values the pileup walker needs. The walker turns it into per-position
  evidence. Reused from production unchanged.
- **STR** — the repeat the read shows: either the measured repeat between both flanks, or, when the
  read ran off its own end mid-repeat, the part it did prove plus which flank it was anchored to (a
  lower bound on the length, kept rather than discarded).

A read that the transform itself defeats yields **no observation** — a normal, tallied per-read
result, never a silent drop and never a masked run failure (§7).

---

## 6. The interface — per read, statically dispatched, with reused buffers

Three properties every implementation upholds:

- **Per read, and independent of every other read.** A read is prepared in isolation: no mate, no
  neighbouring read, no locus-column context. This is what makes preparation parallel with
  deterministic output, and it is why the genuinely pairwise work (reconciling overlapping mates'
  qualities) sits downstream in the gatherer instead. It is also why reassembly cannot be a mode
  here: assembling haplotypes needs every read in a region at once.
- **It re-places what the read already says; it does not invent sequence.** Canonicalizing an indel
  that is already there, redoing a line-up, measuring a repeat that was already sequenced — all
  re-express the read's own content.
- **No usable observation is a result, not an error** (§7).

**Dispatch is resolved at compile time.** Preparation runs on every read of every sample — billions
of calls — so which implementation runs is fixed by a generic type parameter the compiler
specializes, never by a trait object (`Box<dyn …>`): a virtual call and its indirection per read is
a cost this path cannot carry. The per-read *mode* (§2), which genuinely varies read to read, is a
matched enum rather than a second dispatch mechanism.

**Buffers are caller-owned and reused.** The alignment algorithms preparation calls need matrices,
and allocating them per read is the other cost this path cannot carry. So preparation threads a
reusable scratch value, exactly as the existing read-likelihood models do.

```rust
pub trait ReadPreparer {
    /// What this implementation needs to know about the locus — path-owned: nothing on the
    /// generic path, the microsatellite locus (motif, boundaries, flanks) on the STR path.
    type Locus;
    /// The per-read observation this implementation emits — path-owned (§5).
    type Prepared;
    /// Reused buffers, including those the alignment algorithms need — allocated once per
    /// worker, never per read.
    type Scratch: Default;

    /// Prepare one filtered read. `None` means this read yields no usable observation here —
    /// a tallied per-read outcome, not a run error (§7).
    fn prepare_read(&self, read: &MappedRead, locus: &Self::Locus,
                    scratch: &mut Self::Scratch) -> Option<Self::Prepared>;
}
```

The implementation **holds its own reference accessors** as fields and fetches what it needs around
each read's span; there is no reference-window argument. The generic transform needs two views of
the reference at once (raw bytes for left-alignment, canonical bytes for the quality-capping model),
which a single materialised window cannot carry.

---

## 7. Error model

Two outcomes, and they must not be confused:

- **A read produces no usable observation** — `prepare_read` returns `None` and the reason is
  tallied. Normal. The reasons differ by path (the generic path's quality-capping step can decline a
  read; on the STR path a read may anchor no flank at all, or fail the base-quality gate) and each
  path spec enumerates its own.
- **A reference fetch fails** — a contig mismatch, a window past a contig end. This is a broken run,
  **fatal**, and surfaced as such. It is never folded into a per-read `None`, so that `None` always
  means "this read, unusable here" and never hides a broken reference.

Every `None` is counted by reason. A run that silently prepares nothing must be indistinguishable
from a run that prepared everything only by reading the counts.

---

## 8. What changed from the 2026-07-14 draft

Recorded so the path specs are not silently contradicted:

- **Alignment moved out.** The earlier draft had each preparer owning its alignment machinery. The
  alignment algorithms now live in their own module ([`alignment.md`](alignment.md)), which knows
  nothing about caller steps; preparation composes them. The STR path spec's `ViterbiScratch`-holding
  preparer sketch becomes "holds the repeat-aware best-path aligner and its scratch".
- **Re-align is a new generic mode.** The generic path spec describes itself as the
  *trust-the-mapper* implementation that "does not realign". It gains a third mode (§2) for reads
  whose placement is not trusted. Its two existing modes are unaffected.
- **Reassembly is out of scope, not deferred.** Both the shared draft and the generic path spec list
  local reassembly as a future sibling to bench. It is not: production already beats GATK on generic
  loci without it (§1). Those entries should be struck.
- **The trait gained a reused-scratch parameter** (§6). The earlier sketch was
  `prepare_read(&self, read, locus)`; the alignment algorithms it now calls need matrices that
  cannot be allocated per read.
- **Dispatch is explicitly static.** The earlier draft did not say; at this call volume a trait
  object is not an option (§6).

---

## 9. Open questions

- **What flags a region as "not to be trusted", and when?** (§4) Region typing never looks at reads;
  the gatherer looks at reads but runs too late. Either a new producer or a two-pass arrangement.
  This decides whether the re-align mode is worth its cost, so it should be settled before that mode
  is built rather than after.
- **Does pass-through skip the base-quality capping too, or only the left-alignment?** (§2)
  Left-alignment is provably neutral on an indel-free read; the capping step is not obviously so,
  since it measures alignment confidence rather than indel placement. Cheap to settle by comparing
  output on indel-free reads with the step on and off.
- **Is a fast path for clean microsatellite reads worth it?** (§3) Aligning every spanning read is
  the decision in force, with a direct unit-counting fast path parked. Whether the saved time is
  worth a second code path — and whether the two agree on clean reads — is a measurement nobody has
  run.
- **Which reads reach the STR preparer at all.** The lower-bound observation class only exists if
  the read selection upstream admits partially-covering reads; production's spanning gate excludes
  exactly those. Owned by the STR path spec, repeated here because a change to selection is easy to
  miss and makes the class unreachable.
