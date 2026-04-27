# Stage 1 — per-sample caller (CRAM → per-position records)

**Status:** Draft, 2026-04-27. Detailed specification of Stage 1 of the
multi-sample calling pipeline. Slots underneath
[calling_pipeline_architecture.md](calling_pipeline_architecture.md),
which already commits to BAQ-in-process, five scalars per allele,
VCF-style indel anchoring, kept sub-threshold observations, no records
on uncovered positions, and per-allele phase chain identifiers. This
document specifies how those outputs are produced from the input
CRAM(s); the byte-level layout of the artefact written to disk is
Stage 2's contract and is not redefined here.

Open design questions are tagged **[DECISION]** (a choice we have
to make) or **[QUESTION]** (a sub-question or detail to confirm).

## Purpose and scope

One Stage 1 invocation processes **one sample** end to end:

```
sample.cram ───┐
sample.cram ───┼──► Stage 1 per-sample caller ──► sample.psf
sample.cram ───┘
                  reference.fasta (+ .fai)
```

Inputs are one or more **coordinate-sorted CRAM files belonging to the
same biological sample** plus the reference FASTA the CRAMs were
aligned against. Output is a single `.psf` artefact (Stage 2 format)
that downstream stages can read as a stream of per-position records
without ever touching CRAM again.

A "sample" here means the SM tag of the read groups in the input
CRAMs. Read groups can be many; the SM tag must be one and the same
across every input file (see §"Multi-CRAM ingestion" below).

What Stage 1 does *not* do, by design:
- It does not call variants (no minimum read counts, no minimum allele
  ratio, no thresholding of any kind beyond per-read filters listed
  below).
- It does not consult any other sample's data.
- It does not build a `.psf` index. Records are written in coordinate
  order and the file is consumed sequentially by Stage 3.

## Inputs

### CRAM files

- One or more files. Each must be coordinate-sorted (`SO:coordinate`
  in `@HD`).
- All files must share the same `@SQ` list, in the same order, with
  matching lengths (and `M5` checksums when present). A mismatch is a
  hard error.
- All `@RG` SM tags across all files must be identical. A read-group
  SM that disagrees with the others is a hard error — Stage 1 refuses
  to silently mix samples.
- Indexes (`.crai`) are **not required** for Stage 1's main streaming
  path; Stage 1 reads each file from start to end. They become
  required only if `--region` is used (see §CLI).

### Reference FASTA

- Indexed FASTA (`.fai`). Required to decode CRAM in the first place.
- The FASTA's contigs and lengths must agree with the CRAM `@SQ`
  list. The CRAM `@SQ M5` checksums (when present) must match the
  reference; mismatches are a hard error, never a warning.
- Stage 1 also uses the reference directly for two of its own
  purposes:
  - REF-string lookup when emitting allele strings.
  - BAQ HMM input (see §BAQ).

So the reference is required by both the CRAM decoder and the caller
itself; we hand the same FASTA to both.

## Output

Stage 1's product is a stream of **per-position records** consumed by
Stage 2's encoder. Each emitted record has the fields fixed by the
architecture doc:

- `delta_pos` (distance from the previous emitted position, encoded
  by Stage 2; Stage 1 just emits absolute positions and lets Stage 2
  delta-encode).
- chromosome id (Stage 2 carries it on chromosome change).
- Phase-chain lifecycle markers: `new_chains` (slot ids that started
  since the previous emitted position) and `expired_chains` (slot ids
  that ended). See §"Phase chain identifiers" for the slot model.
- `n_alleles` and the per-allele payload: the literal allele sequence,
  the five scalars below, and the list of chain slot ids that
  contributed this observation.

**Per-allele scalars** (committed in the architecture doc, repeated
here as the work product Stage 1 must produce per allele observed at
each position):

| Scalar | Source |
|---|---|
| observation count | number of supporting reads |
| `Σ max(ln_BQ, ln_MQ)` | sum over supporting reads of `max(ln_BQ_BAQ, ln_MQ)`, i.e. with BAQ-adjusted BQ |
| forward-strand count | reads on the forward strand |
| placed-left count | reads whose mapped 5′ end is to the left of the current position (freebayes' `placedLeft`) |
| placed-start count | reads whose mapped 5′ end *is* the current position (freebayes' `placedStart`) |

`ln_BQ` here is `ln(P(error))` derived from `min(BQ, BAQ)` after BAQ
capping; `ln_MQ` is `ln(P(misalignment))` derived from the read's MAPQ.
Indel allele BQ uses the min-over-window proxy already specified in
the architecture doc.

The output stream is sequential, monotonically non-decreasing in
`(chrom_id, pos)`; Stage 1 never seeks backwards.

## Library choice for CRAM I/O

CRAM decoding is non-trivial (slice-based, reference-driven, multiple
codec versions, optional embedded reference, rANS/gzip/bzip2/lzma
external-block codecs). We are not implementing it ourselves. We
will use the `noodles` library.

`noodles` (`noodles-cram`, `noodles-sam`, `noodles-fasta`,
`noodles-bgzf`).** Pure Rust, actively maintained, no C dependency.
Streaming API; we can drive the slice/container level directly if we
want fine-grained parallelism.

Target CRAM 3.x only. At header time, detect the version byte and emit a clear error if it's 4.x

## Multi-CRAM ingestion

A sample's reads may be split across several coordinate-sorted CRAM
files (typical when a sample has been re-sequenced, when lanes were
mapped separately, or when a bioinformatics pipeline emitted one CRAM
per chromosome / per shard). Stage 1 has to walk all of them as if
they were a single coordinate-sorted stream.

### Pre-flight validation

Before any reads are pulled, Stage 1 reads the header of every input
CRAM and checks:

1. `@HD SO:coordinate` on every file.
2. The `@SQ` list — names, lengths, and `M5` when present — is
   identical across files (and matches the FASTA). Order matters:
   coordinate sort is meaningful only relative to a `@SQ` order, so
   inputs must agree on it.
3. Every `@RG` carries an SM tag. Every SM across every file is the
   same string. The set of `@RG` IDs is otherwise unconstrained — a
   sample with N read groups across M files just lists them all.

A failure in any of these halts Stage 1 with a specific error
message naming the offending file and field. This is consistent with
the project's "errors must not pass silently" principle.

Should the CLI also accept an explicit `--sample-name`
argument used to *select* which SM to keep when a CRAM (rarely)
contains multiple? We don't need it right now, so let's follow the
YAGNI advice and not implement it until the case comes.

### Peek-and-scan coordinate merge

Once validated, Stage 1 wraps each CRAM decoder's record iterator in
a `BufferedPeekable` (see [buffered_peekable.md](buffered_peekable.md))
and merges them by `(ref_id, pos)` using the same peek-and-scan idiom
that `variant_grouping.rs` already uses to merge per-sample variant
streams. Each peekable exposes the head record without consuming it;
on every step we scan the heads of all peekers, pick the smallest,
and consume from that one:

```
peekers = [BufferedPeekable::with_buffer_size(decoder.records(), 1)
           for each input CRAM]

while any peeker has a head:
    let smallest_idx = argmin over peekers of peeker.peek()?.key()
        // key = (ref_id, pos, file_index); file_index breaks ties
        // so the merge is deterministic
    let read = peekers[smallest_idx].next().unwrap()?
    feed read into the pileup walker
```

`file_index` is the tiebreaker on equal coordinates so the merge is
deterministic. The scan cost per pulled read is `O(k)` in the number
of input files; with `k` typically ≤ a handful (and rarely more than
~16 even for per-lane shards), this is dwarfed by CRAM decode and
BAQ. No parallel peek is needed at this scale — Stage 1's per-sample
peek runs sequentially.

Stage 1 constructs the peekers with `buffer_size = 1` because the
CRAM decoder underneath already does its own slice-level batching;
an extra outer buffer would just delay records without saving work.
See [buffered_peekable.md](buffered_peekable.md) §"Why it exists in
the project" for the rationale on per-caller buffer sizes.

The merge loop also validates that each iterator's head record is
not before its predecessor (the file's own coordinate sort
invariant); a violation halts Stage 1, naming the file and offending
read.

#### Why peek-and-scan rather than a min-heap k-way merge

A min-heap merge would be `O(log k)` per pull instead of `O(k)`, but
at `k ≤ ~16` this is single-digit operations either way and the
choice is dominated by code clarity. Two reasons we prefer
peek-and-scan here:

- **Project consistency.** `variant_grouping.rs` already merges
  per-sample streams by peek-and-scan
  (`compute_next_span_seed_and_skip_non_variable`,
  [variant_grouping.rs:172](../../src/variant_grouping.rs#L172)).
  Using the same idiom in Stage 1 means a reader does not have to
  switch mental models between the per-sample and cohort sides of
  the pipeline. Both sides will share the `BufferedPeekable`
  utility once `VarIterator` is retrofitted to use it (see
  [buffered_peekable.md](buffered_peekable.md) §"Follow-up:
  retrofit `VarIterator`").
- **Validation in the peek.** The merge loop can naturally
  validate order, skip records, or filter on the fly — the same
  shape `variant_grouping.rs` exploits. A heap entry is just an
  ordering key and is awkward to combine with side checks.

This is the only place Stage 1 cares about the multi-file shape;
everything downstream of the merge sees one coordinate-sorted stream
of `(read, alignment-position, decoded-sequence-and-quals)` tuples.

### Duplicate-read detection across CRAMs

If two input CRAMs both contain the same read — same `QNAME`, same
flag, same alignment coordinate — that is an upstream bug. Two
input files should never legitimately emit the same physical
fragment twice. Stage 1 detects this and aborts with a hard error
naming both source files and the offending `QNAME`; no warn-and-continue.

**No global hash table.** Detection happens at the position
boundary, not by hashing every read across the whole file. A true
duplicate necessarily produces two reads at the same `(ref_id, pos)`,
which the peek-and-scan merge surfaces in close proximity. We
therefore only need to compare each pulled read against the small
window of reads already accepted at the current start coordinate
(bounded by per-position depth — typically tens of reads at 30×
coverage, far less at the 2–10× the pipeline targets). The check
is `O(depth)` per read at start time and adds no memory beyond the
already-active read set. This catches the bug cheaply without ever
materialising a per-file `QNAME` hash.

### Why merge at the read level rather than at the per-position output

A natural alternative is to run the full Stage 1 pipeline once per
CRAM and merge the per-position records at the end. We do not do
this:

- The per-position records are *per-allele scalar accumulations*.
  Merging them would require summing scalars from multiple files
  per allele per position — possible, but it doubles the
  bookkeeping and makes the phase chain id story messier (chain
  ids must be unique across the merged stream).
- Reads spanning the boundary between two files (extremely common
  for per-chromosome shards if anyone ever produces such a thing,
  rare for per-lane shards but still possible) would split their
  evidence between two intermediate files, and the per-allele
  scalars would have to be reconciled across them.
- Read-level merge is dirt cheap and gives us a single,
  conceptually simple downstream pipeline. We pay for one peek per
  input file per pulled read, with the number of input files
  typically a small handful — negligible against CRAM decode and
  BAQ.

## Read filters (per read, before any per-position accumulation)

Every read coming out of the merged stream is run through a small
filter cascade. Reads that fail any check are dropped *before* any
allele evidence is accumulated. Filters are:

| Filter | Default | Rationale |
|---|---|---|
| Unmapped (`flag & 0x4`) | reject | no alignment, no allele |
| Secondary (`flag & 0x100`) | reject | the primary alignment of the same read carries the canonical placement |
| Supplementary (`flag & 0x800`) | reject | same reason; chimeric segments are reported separately as primary |
| QC fail (`flag & 0x200`) | reject | aligner-marked failure |
| Duplicate (`flag & 0x400`) | reject | duplicate-marked reads inflate coverage without adding evidence |
| MAPQ < `--min-mapq` | reject | low-confidence placement |
| Read length 0 / no SEQ | reject | nothing to score |
| First or last CIGAR op is `I` or `D` | reject for that indel only | already specified in architecture doc — placement is untrustworthy |

Defaults:

- `--min-mapq 20` (CLI-tuneable) — matches bcftools' default and
  is the freebayes ballpark too. Conservative; the CLI exposes
  it as `--min-mapq INT` so users can override per run.
- No `--min-bq` filter at the read level. BAQ replaces "low-BQ
  bases are bad" by capping effective BQ; bases keep their evidence
  and contribute proportionally to their quality. A hard threshold
  here would drop information that the joint stage may later want.
- max-coverage cap by default. Adding `--max-depth N` is cheap
  and handy for runaway repeat regions; recommended even though
  Stage 3's DUST filter handles most of the same cases.

### Mate-overlap handling

For paired-end reads where the two mates' alignments overlap (a
common short-fragment situation), both mates report the same physical
bases of the fragment as two "observations." Treating them as
independent inflates evidence and biases per-allele scalars.
Stage 1 collapses each overlap so every reference position is
counted once.

**Per-position BQ comparison, not arbitrary clipping.** At each
reference position covered by both mates of the same fragment,
keep the observation from the mate with the higher BAQ-adjusted
base quality (`min(BQ, BAQ)`); zero out the other mate's
contribution at that position. This is strictly better than
bcftools' "always pick mate 1" rule at low coverage — every BQ
unit preserved matters at the 2–10× target — and costs almost
nothing on top of the pair-aware bookkeeping the per-pair phase
chain id already needs (see §"Phase chain identifiers").

**MAPQ is not used as a tiebreaker.** Both mates of a pair almost
always carry the same MAPQ, so MAPQ does not discriminate at the
per-base level. The decision is per-position, on per-position data.

**Disagreement at an overlap position.** If the two mates call
different bases at the same overlap position (one says A, the
other C — usually a sequencing error in one of them), the
higher-BQ call wins as-is, no further downgrade. (a) is simple
and consistent with the rest of the spec; (b) downgrade-on-disagree
and (c) drop-from-both are plausible alternatives but add a
special case for what is a small fraction of overlap positions.

**Implementation hook.** Coordinate-sorted reads put mate 1 in the
active set before mate 2 arrives, so when mate 2 enters, the
walker can look up mate 1's BAQ-adjusted BQ at each overlapping
ref position from mate 1's active-set entry. Per overlapping base:
one comparison plus one keep-or-zero choice. Bounded by
per-fragment overlap (typically tens of bases). The "zero" is
applied by setting the losing mate's BAQ-adjusted BQ to 0 at that
position so it contributes zero log-likelihood mass; bases are not
deleted (which would complicate CIGAR walking).

**Tie-breaking on equal BQ.** Prefer mate 1 (read flag `0x40`)
deterministically.

### N CIGAR ops (skipped, RNA-seq data accepted)

Reads with `N` CIGAR ops (RNA-seq splices) are processed normally:
the `N` span produces no allele evidence (the read does not observe
those reference positions, so they get no contribution from this
read), and the read continues on the other side of the splice with
its bases contributing to allele evidence at their respective
reference positions. The read's phase chain identifier covers the
whole physical molecule — both sides of the splice — even though
the chain's coordinate footprint has a gap in the middle.

This means **RNA-seq CRAMs are valid input**, not a misuse to be
warned about. BAQ runs on each contiguous alignment segment
between `N` ops independently; the segment boundaries fall
naturally on the splice boundaries. For DNA workflows the rule is
a no-op (no `N` ops appear).

### Soft-clipped and hard-clipped bases

Soft-clipped bases (`S` CIGAR op) **do not contribute to allele
evidence**. The aligner has marked them as not aligned to the
reference, and the project follows the same rule as bcftools and
freebayes: ignore them entirely. Hard-clipped bases (`H`) are
already absent from `SEQ` by construction and therefore contribute
nothing automatically.

## BAQ application

BAQ runs *after* the per-read filter and *before* any allele
accumulation, on each surviving read. Heng Li's 2011 algorithm:

1. For each read, take a small reference window covering the read's
   alignment span, with a few dozen bases of flank on each side.
2. Run the BAQ HMM (gap-open / gap-extend transitions, base-quality
   driven match/mismatch emissions) to compute, per read base, the
   posterior probability the base is correctly aligned at its
   reported position.
3. Convert that posterior into a BAQ score and **cap** the read's
   BQ at `min(BQ, BAQ)` per base.

### Filter + BAQ as a single per-read stage

Both the read filter and BAQ are purely per-read: they only need the
read itself plus the reference FASTA, no shared state, no cross-read
dependency. Stage 1 fuses them into **one** processing stage rather
than chaining two with intermediate buffering. The merge produces a
stream of reads in coordinate order; we accumulate a small batch
(default ~1000 reads, tuneable), run filter-then-BAQ in parallel via
rayon's `par_iter_mut()` within the batch, and hand the surviving,
BAQ-adjusted reads to the pileup walker in coordinate order. Filter
runs first inside each worker so reads that get rejected never pay
the BAQ HMM cost.

Coordinate order is preserved between batches; within a batch the
work is fully parallel. The walker only ever sees BAQ-adjusted
qualities, so all five scalars downstream are computed with
`min(BQ, BAQ)`.

**Mate-overlap handling is not part of this stage.** It requires
looking up the partner mate from the active set, which is sequential
walker state, so it lives inside the pileup walker (see
§"Mate-overlap handling"). Keeping it out of the parallel stage
avoids cross-worker synchronisation on the active-set lookup.

BAQ implementation source.

- Port `samtools/bam_md.c`'s BAQ routine to Rust directly. ~300
  lines, well-tested logic, no surprises.
- `noodles-sam`'s does not have a BAQ


Extension model: extended-BAQ
(`samtools calmd -E`). `-E` is more aggressive at indel sites and
is bcftools' default.

## Pileup walking

After filtering and BAQ, each read enters the pileup walker. The
walker iterates the reads in coordinate order, maintains an
**active read set** (reads whose alignment overlaps the current
reference position), and emits one per-position record at every
position covered by at least one active read.

### Per-read CIGAR consumption

For each read, walk the CIGAR ops alongside the reference position:

- `M`/`=`/`X`: at each consumed reference position, the read
  contributes one **base call** (the read's base at that ref
  position) with quality = BAQ-adjusted BQ.
- `I` (insertion): the inserted bases are anchored at the
  *previous* reference position (VCF anchoring). The anchor
  position's allele carries the literal `anchor_base + inserted_bases`
  and the indel BQ proxy (min of `l + 2` BAQ-adjusted BQs centred on
  the indel — already specified in the architecture doc).
- `D` (deletion): the deleted span produces one allele entry at
  the *position before* the deletion (anchor), with allele =
  anchor base alone and `ref_span = l + 1`. Interior positions
  (the actually-deleted bases) get no contribution from this read.
- `N`: skip the span (no evidence at those positions; the read
  resumes contributing on the other side of the splice).
- `S`: soft-clip — ignored.
- `H`: hard-clip — not in SEQ, no action.
- `P`: padding — no action.

A read whose first or last CIGAR op is `I` or `D` does *not*
contribute to that indel allele (per architecture doc). It can
still contribute to surrounding non-indel positions.

### Position emission rule

Stage 1 emits a record at position `p` exactly when:

1. There is at least one active read covering `p` *with a base or
   indel observation at p*. (A read whose CIGAR carries an `N` over
   `p` is in the active set but does not produce evidence at `p`,
   and on its own is not enough to emit a record.)
2. After accumulating evidence at `p`, the record has at least one
   allele.

In particular: a position covered only by reads whose CIGAR skips
it (all-`N`) produces no record. A position covered by reads but
all of them are mate-overlap-clipped to BQ 0 still produces a
record — the *observation* exists, the *quality* just happens to
be zero. (The five scalars handle this faithfully: the obs count
goes up, the quality scalar contribution is `ln(0.5)` or whatever
zero-BQ resolves to.)

### Allele extraction at a position

At each emitted position the walker collects, across all active
reads contributing at that position:

- For `M`/`=`/`X`: a SNP-shaped allele = the literal base from the
  read.
- For an `I` anchored at this position: an INS-shaped allele =
  `anchor_base + inserted_bases`, with `ref_span = 1`.
- For a `D` anchored at this position: a DEL-shaped allele =
  anchor base alone, with `ref_span = l + 1`.
- For combinations (a SNP and a deletion both anchored at the same
  position, etc.): the architecture doc's "extend the anchor REF"
  rule kicks in — REF is widened to span the longest event, and
  every co-occurring allele at this position is rewritten against
  that extended span.

Allele strings are uppercase `{A,C,G,T,N}` literals. No `<*>`,
no `<NON_REF>`, no `*`. (Reads with `N` bases contribute an `N`
literal; the joint stage handles `N` consistently.)

### Five-scalar accumulation
I 
For each (position, allele) `allele_evidence`, the walker maintains
the five running scalars listed in §"Output." Updates are O(1) per
contributing read:

```
allele_evidence.num_obs       += 1
allele_evidence.q_sum         += max(ln_BQ_BAQ, ln_MQ)     // already in log space
allele_evidence.fwd           += (read on forward strand) ? 1 : 0
allele_evidence.placed_left   += (read 5′ < pos) ? 1 : 0
allele_evidence.placed_start  += (read 5′ == pos) ? 1 : 0
```

`ln_BQ_BAQ` for indel alleles is the min-over-window proxy already
specified.

### Active-set bookkeeping

The active set is a small structure keyed by read pointer; reads
enter when their first non-soft-clip aligned position ≤ current
position, and leave when their last aligned position < current
position. With coordinate-sorted input the entry/exit operations
are amortised O(1) per read. Memory is bounded by the maximum
read pile depth × read length, which on 2–10× coverage is trivial.

## Phase chain identifiers

Phase chain identifiers let Stage 5 reason about compound haplotypes
(see architecture doc §"Compound haplotype alleles"). Stage 1's job
is to assign them and emit them per-record in a form Stage 2 can
write verbatim — no hashing, no global index, no delta-encoding work
for Stage 2.

### What a phase chain is

A **phase chain** is a maximal set of allele observations that all
come from the same physical read or read pair. If sample s's read r
supports allele A at position p1 *and* allele B at position p2, then
A's observation at p1 and B's observation at p2 belong to the same
chain — which lets Stage 5 know that, in this sample, A and B
co-occur on the same haplotype.

### Per-record encoding: slot ids + lifecycle markers

Chains are short-lived: a read or pair contributes for at most ~1 kb
(PE fragment length), so the number of *concurrently active* chains
at any reference position is bounded by per-position depth —
typically 5–15 at the project's 2–10× coverage target, rarely above
~50 even at 30×. That makes a tiny slot-based encoding much cheaper
than full chain identifiers.

Each per-position record carries:

- `new_chains`: list of slot ids that started since the previous
  emitted position (a chain begins when its read enters the active
  set).
- `expired_chains`: list of slot ids that ended since the previous
  emitted position (a chain ends when its read or pair has fully
  left the active set).
- For each allele: `chain_slots`, the small list of slot ids of the
  chains that contributed this observation (encoded as a bitmap over
  currently-active slots, or a short list, at Stage 2's
  discretion — but the data Stage 2 sees is just slot ids).

Slot ids are small integers (`u8` is enough for the project's
coverage target; `u16` is the defensive choice and still cheap).
They are *recycled* once a chain has fully expired: the next new
chain can occupy a freed slot. The lifecycle markers ensure no
consumer ever confuses two chains that successively occupied the
same slot — see §"Slot identity invariant" below.

This is intentionally analogous to VCF's `|` / `/` phase markers:
within an open phase block the slot id is the linkage; the markers
declare where blocks open and close. The invariant the format
preserves is the same one Stage 5 needs: "two observations carrying
the same active slot id come from the same physical molecule."

### Slot allocation

Stage 1 maintains an internal table of currently active chain slots.
When a new read pair (or solo read) enters the active set:

1. Allocate the lowest free slot id from the table.
2. Tag every allele observation that read pair contributes to with
   that slot id.
3. Add the slot id to the next emitted record's `new_chains` list.

When a read pair fully leaves the active set:

1. Release its slot id (the slot becomes available for reuse).
2. Add the slot id to the next emitted record's `expired_chains`
   list.

Within a single record, multiple alleles may reference the same
slot only in the rare case of mate-overlap with disagreement. Across
records, a slot id refers to the same chain for the duration of its
active lifetime, which is bounded by read/pair span.

Phase chain slots are not shared across samples. Stage 1 only emits
within-sample chains; cross-sample reasoning is the joint stage's job.

### Read pairs share one slot

Both mates of a paired-end read share a single slot id. An allele
observed on mate 1 at position p1 and an allele observed on mate 2
at position p2 come from the same physical DNA molecule, so they
are evidence for the same compound haplotype in this sample.
Issuing two distinct slots would force Stage 5 to reconstruct the
linkage from QNAMEs (or, more likely, miss it entirely); doing it
once at Stage 1 keeps the linkage explicit at no extra cost.

**Assignment mechanics.** Stage 1 maintains a small QNAME → slot id
map of pairs whose first mate has been seen but whose second mate
has not yet been processed. When a read arrives at the per-read
stage:

1. Look up its QNAME in the map.
2. If found: reuse that slot id, then drop the entry — once both
   mates have been tagged, the QNAME → slot binding is no longer
   needed.
3. If not found: allocate a fresh slot id, use it for this read,
   and register `(QNAME → slot)` if the read flags say a mate
   exists.

**Solo reads.** Single-end runs, mate-filtered-out cases, and
malformed inputs where only one mate appears all converge to the
same shape: the surviving mate keeps the slot id it was issued on
first-seen. Stage 5 cannot tell a chain-of-one from a chain-of-two,
and does not need to — the slot semantics are identical: "alleles
tagged with this slot were observed on the same physical molecule."

**Map size.** Bounded by reads whose mate has not yet been seen at
the current merge position. With coordinate-sorted CRAMs and
typical short-read fragment sizes (a few hundred bp to ~1 kb), this
is on the order of tens to a few hundred entries even at 30×
coverage and shrinks linearly with coverage. We add a defensive
upper bound (e.g. 1 Mbp window past the first-mate position) so
malformed inputs cannot grow the map without limit; entries that
exceed the window are released and their first mate is treated as
solo. A whole-file QNAME hash table is therefore unnecessary.

### Slot identity invariant

At any reference position, every slot id present in `chain_slots`
of any allele refers to a *single, unique chain* — the chain that
currently occupies that slot. Slot ids are reused across the file,
but never within an active lifetime: a slot enters `expired_chains`
in the record at which its chain ends, and only after that
expiration can the same slot id reappear in `new_chains` of a later
record. The lifecycle markers thus preserve the conceptual
invariant "different chains are never confused" without ever
materialising global chain identifiers.

This is sufficient for everything Stage 5 needs. Walking records
sequentially, Stage 5 maintains its own active-slot view (apply
`expired_chains`, then apply `new_chains`, then read each allele's
`chain_slots`); the compound-haplotype check between two records
becomes a small intersection of slot-id sets. No hashes, no
per-file index, and no delta-encoding decisions for Stage 2 — it
just writes the integers Stage 1 hands it.

## Output emission

Stage 1's output side is a thin layer that:

1. Receives a per-position record from the pileup walker (held in
   memory: `chrom_id`, `pos`, the `new_chains` and `expired_chains`
   slot lifecycle markers, and the list of `allele_evidence`
   entries with their five scalars and per-observation slot lists).
2. Hands it to the Stage 2 encoder (compression, framing,
   variable-length encoding — all Stage 2's responsibility).

Stage 1 buffers a small fixed number of recently-emitted positions
internally (so the active-set walker can backfill scalars when a
read with a large insertion arrives slightly out of left-edge
order — rare but possible at insertion boundaries), and flushes
positions to Stage 2 as soon as the active set guarantees no more
contributions can arrive at them.

## Errors

Per the project principle, every problem is a hard error with a
specific message:

- CRAM header mismatch (SQ, SO, SM disagreement across files).
- FASTA mismatch (`@SQ` MD5 vs. FASTA), missing FASTA, missing
  `.fai`.
- Out-of-order reads within a single CRAM (would mean the file is
  not actually coordinate-sorted).
- BAQ HMM numerical failure (should not happen; if it does, abort
  with read context).
- I/O errors (CRAM read, FASTA read, output write) propagate as
  errors, not warnings.

Warnings (logged to stderr, do not halt):

- Reads filtered out — counts only, not per-read messages.

## Parallelism

- **Across samples:** Stage 1 is run independently per sample. The
  outer driver (CLI / shell pipeline) parallelises across samples
  by spawning one Stage 1 process per sample, or by running them
  serially. We do not bake cross-sample parallelism into one
  process.
- **Within a sample:** the work is pipelined.
  1. CRAM decoders run in their own threads (noodles' container
     iterator is naturally streaming; one decoder thread per input
     CRAM).
  2. The peek-and-scan merge runs on a coordinator thread.
  3. Filter + BAQ run as a single rayon-parallel per-read stage
     over batches of merged reads (purely per-read; no cross-read
     dependency). Coordinate order is preserved between batches;
     within a batch the work is fully parallel.
  4. The pileup walker runs single-threaded — it is inherently
     sequential per sample (active-set bookkeeping, mate-overlap
     resolution, allele accumulation).
  5. The Stage 2 encoder runs on its own thread, reading from a
     bounded channel fed by the walker.

  The pipeline is bounded by either CRAM decode or filter+BAQ
  depending on coverage. Both are amenable to thread-count tuning
  via a single `--threads` flag.

## CLI

Sketch:

```
caller per-sample \
    --reference ref.fa \
    --output sample.psf \
    [--threads N] \
    [--min-mapq 20] \
    [--max-depth INT] \
    [--region chrom[:start-end]] \
    sample_part1.cram sample_part2.cram ...
```

- Reference is mandatory.
- Output path is mandatory; `--output -` writes to stdout.
- One or more CRAMs as positional args.
- `--region` requires `.crai` for every input CRAM. Without
  `--region`, indexes are not consulted.
- Threading defaults to all logical cores.
- The sample name does *not* appear on the CLI; it is read from
  the CRAMs' `@RG SM` tag (and validated to be unique across files,
  per §"Multi-CRAM ingestion"). **[DECISION]** confirm.

## Module layout (preview)

The Rust source layout this spec implies, for a future
implementation plan:

```
src/
  per_sample_caller/
    mod.rs              — CLI entry, orchestration
    cram_input.rs       — header validation, multi-file peek-and-scan merge
    read_filter.rs      — per-read filters (flags, MAPQ, ...)
    baq.rs              — BAQ HMM (port of samtools' implementation)
    pileup_walker.rs    — active-set, CIGAR consumption, allele extraction
    five_scalars.rs     — per-allele scalar accumulation
    phase_chain.rs      — slot allocator/recycler + per-read tagging
                          + new/expired marker emission
    psf_writer.rs       — Stage 2 encoder (separate spec)
```

The `psf_writer.rs` module is sized by Stage 2's spec, not this one.
