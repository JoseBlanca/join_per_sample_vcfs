# ng — the alignment file: a validated handle and its ordered region reads

*Status: design spec, 2026-07-20. **All open questions resolved** (owner, 2026-07-20): O1 MD5 handling,
O3-file → rebuild ng-owned, O5 → `src/ng/read/input/`, missing `M5` → not an error, reader pooling →
built now. Ready for its arch doc. Split out of the former `read_ingestion.md` (owner, 2026-07-20)
together with [`sample_reads.md`](sample_reads.md): this document owns everything whose subject is
**one alignment file**, that one owns everything whose subject is **the sample**. It opens and
validates a single indexed BAM/CRAM against the reference and yields its reads for a region in
**coordinate order**, filtered by step 1. It closes the two holes
[`reference_info.md`](reference_info.md) §1/§7 explicitly deferred — the `@SQ`↔reference
**permutation hole** and the **coordinate-sort-order check**. Under [`ng_proposal.md`](ng_proposal.md)
§1 and [`../arch/module_layout.md`](../arch/module_layout.md). Naming: **STR** in prose, `ssr` in
code. Code-facing companion: [`../arch/alignment_file.md`](../arch/alignment_file.md) (types &
interfaces).*

---

## 1. What this is — scope and non-goals

This module turns *one alignment file on disk* into *a validated handle that answers region queries
with an ordered stream of reads*. Two things, in this order:

1. **Open and validate.** Open an indexed BAM/CRAM, prove it is coordinate-sorted, prove its `@SQ`
   list is the reference's contig table, parse its index, and extract its sample name (§3.1).
2. **Serve regions.** On a `(contig, start, end)` query, index-seek and stream every overlapping read
   in coordinate order, filtered by step 1's `ReadFilter`, verifying as it goes that the file really
   is sorted (§3.2).

Nothing here knows that a sample can have several files. Assembling k of these into one stream is
[`sample_reads.md`](sample_reads.md)'s job, and it is a strict layer above: it consumes what this
module produces and adds no reader logic of its own.

### The design was chosen (owner, 2026-07-20)

- **Indexed region queries, not a whole-file linear scan.** The analysers pull reads region by region
  (an STR locus is tiny; the pileup walks a stretch), so this is index-seek + bounded scan per region,
  mirroring production's `segment_reader`. A `.bai`/`.csi`/`.crai` next to the input is **required**
  (build-if-missing is a caller policy, §8). The whole-genome linear scan is the deferred alternative
  (§8).
- **Validate inside the open.** Every check runs in the function that opens the file, so an
  unvalidated handle never exists and no read can flow from a file that has not passed (§3.1).
- **Order is verified while streaming, and a violation is fatal.** The whole design assumes
  coordinate-sorted input. That assumption is *checked*, not trusted: every read yielded is compared
  against the previous one, and a regression is a hard error (§3.2). A file that claims
  `SO:coordinate` but is not sorted fails loudly rather than producing quietly wrong calls.
- **Two seams, not one.** The reader→filter boundary passes an *undecoded, reused buffer*, so a dropped
  read is never allocated as a `MappedRead`; from the filter onward every stage is an iterator over
  `Result<MappedRead, _>`, which is what lets the sample layer add or omit a merge above them without
  this module knowing (§3.2).
- **Open once, seek per query.** The reader and the parsed index live for the whole run; a query is a
  seek, never a file open or an index parse. Baseline, not optimisation (§3.3).
- **Wrong-assembly detection is deferred, not skipped.** `@SQ M5` tags are captured at open and
  compared against the reference's true digests when the caller joins `reference_info`'s background
  verification — so startup never blocks on a genome read, and the fault is still caught before any
  output is committed (§3.1).

### Non-goals (deliberately excluded)

- **Anything cross-file.** The merge, the same-file-twice check, the cross-file sample-name agreement,
  per-file count aggregation → [`sample_reads.md`](sample_reads.md).
- **Filtering policy.** Which reads survive is step 1's job (`ReadFilter`,
  [`read_filtering.md`](read_filtering.md)). This module *composes* the filter in (§3.2) but owns no
  thresholds of its own beyond the cheap region-overlap drop that is a *reader* concern, not a filter
  (chunk-edge slop the index over-returns — dropped, uncounted).
- **Per-base evidence work** — adaptor masking, mate-overlap reconciliation, CIGAR decomposition. Per
  `read_filtering.md` §6 these are downstream (`pileup/` / STR prep), not here.
- **Building the reference table.** The canonical contig table comes from
  [`reference_info`](reference_info.md); this module *consumes* a `ReferenceInfo` and reconciles
  against it. It never parses a `.fai` or a FASTA itself.
- **Random-access re-entry / seeking backwards within a served region.** A region query is a forward
  scan; the consumer that wants an earlier region issues another query.

---

## 2. Where it sits — the read path

```
reference.fa ─▶ reference_info ─▶ ReferenceInfo  (canonical ContigId order, per-contig MD5)
                                        │
                                        ▼  §3.1 open gate: SO=coordinate + @SQ == reference
one file  ──▶ index preflight ──▶  AlignmentFile (validated handle: reader + parsed index + SM)
                                        │
                                        ▼  §3.2 reads_in_region(contig, start, end)
                       region source ─▶ ReadFilter ─▶ order-verify ─▶ ordered MappedRead stream
                       raw in-region     (step 1)     (hard error
                        RecordBuf                    on regression)
                                        │
                                        ▼
                          sample_reads.md  (k of these → one stream)
```

The output contract is a stream of `Result<MappedRead, _>` in coordinate order — the same `MappedRead`
(`bam/alignment_input.rs`) step 1 already yields and step 2 already consumes, so nothing downstream
changes shape.

`MappedRead.source_file_index` (already a field) records which file a read came from. This module sets
it and never reads it; it is **downstream-meaningful, not just provenance** — see
[`sample_reads.md`](sample_reads.md) §2, where the several-files-means-several-experiments reasoning
lives.

---

## 3. The concerns

### 3.1 Open and validate (the permutation-hole fix)

`ContigId(i)` means "the i-th contig of the reference, in `.fai`/`ReferenceInfo` order"
([`types.rs`](../../../../src/ng/types.rs)). A read's `ref_id` is the i-th `@SQ` of *its* file. Every
reference fetch (`ReadFilter` #8) — and, one layer up, the merge — assumes these two `i`s coincide.
Today `ReadFilter::new` checks only that each `@SQ` index *resolves* against the reference, so **a file
whose `@SQ` list is a permutation of the reference passes, then fetches the wrong contig for every
read** (reference_info.md §1).

The fix is a fail-fast gate. **Every check lives inside the function that opens the file** — a file is
either opened *and* validated, or it is an `Err`; there is no window in which an unvalidated handle
exists. The checks, in order:

1. **`@HD SO` must be `coordinate`.** Missing or anything else → `NotCoordinateSorted`. (Production's
   `extract_header`, `alignment_input.rs:297`.)
2. **The file's `@SQ` list must equal the reference's contig table — same contigs, same order, same
   lengths** (and same MD5 where both carry one). Built from the SAM header exactly as production's
   `extract_header` does (`alignment_input.rs:306`), compared with
   [`ContigList::first_disagreement`](../../../../src/fasta/mod.rs) (order-aware; `pub(crate)`, ng may
   call it) against `ReferenceInfo::contig_list()`. Any disagreement → `ContigReconcile { file,
   detail }`, naming the first differing field and index. A permutation disagrees on `name` at the
   first transposed position; a re-labelled contig on `length`.
3. **The index must load** — parsed *here, once*, and held for the life of the handle (§3.3).
4. **The `@RG` records must name exactly one sample** (`SM`), extracted here and exposed on the handle.
   Production's `extract_single_sample_name` (`alignment_input.rs:378`) does the same but is
   module-private, so ng reads the header itself (arch §5). *Agreement
   across several files is not checked here — it cannot be, it is not a property of one file;
   `sample_reads.md` §3.1 owns it.*

**The invariant this establishes, which the layer above depends on:** once the gate passes,
`ref_id == ContigId` **by construction** for every read of the file. That is what makes it sound for
`sample_reads.md` to compare `(ref_id, pos)` *across* files without re-checking anything, and it is
also what makes `ReadFilter`'s own resolve-probe redundant (kept or dropped — an arch call, §6).

**The M5 decision (flagged by reference_info.md §7).** A `.fai`-only `ReferenceInfo` carries `md5:
None`, so the MD5 comparison is a no-op (a wildcard, matching `ContigEntry::eq`) and reconciliation is
name+length+order only — the same guarantee production ships. A `ReferenceInfo` read through the
**`Fasta` arm** carries the true per-contig MD5, so reconciliation additionally catches a file aligned
to *a different assembly with the same contig names and lengths* — a check production cannot make (it
trusts the `@SQ M5`).

**The MD5s cost almost nothing, but they are not available when we open the files.** Both halves of
that matter, and they point at a design rather than a trade-off.

*Almost nothing*, because the digests ride along in a pass that happens anyway: they are folded into
the single streaming pass that reconstructs the geometry — one `Md5` per contig and one for the whole
reference, fed from the same uppercase window buffer (`reference_info.rs:600-612`, finalized
`:659-686`). That pass is paid to verify the `.fai` when one exists (reference_info mandates
verification — there is no trust-the-index fast path) or to *build* the `.fai` when one does not, and
`ReferenceInfoCache` is single-flight (`:836-897`) so it happens once per run however many files
reconcile against it. The digests' true marginal cost is CPU over bytes already in cache.

*Not available at open*, because of the deliberate non-blocking startup: with a `.fai` present,
`reference_info` hands back the cheap `.fai` table immediately and runs the verifying genome pass on a
**dedicated spawned thread** (`reference_info.rs:963-989`, an OS thread rather than a rayon worker so a
seconds-long read cannot starve the pool). Those are **two different `ReferenceInfo` objects** —
distinct cache keys, `Fai` vs `Fasta{Some}` — and the one this module receives while the pass is still
running carries `md5: None`. So on the common path there is nothing to compare at open, whatever policy
we write.

**The design: capture the tags at open, compare after the join.** Wrong-assembly detection does not
need the digests *at open* — it needs them *before anything is published*, which is exactly when the
background verification is joined already. So:

1. **At open, capture each `@SQ` line's `M5` tag** into a per-contig `Option<Md5>`, indexed by
   `ContigId` (the §3.1 gate has just proved that indexing is the reference's). This is free — the tags
   are already in the header being parsed for check 2 — and it is *all* this module does with them.
2. **Foreground reconciliation stays name + length + order**, plus the digest comparison in the case
   where the `ReferenceInfo` in hand does carry MD5s (the `.fai`-absent arm, which is synchronous, so
   it does). Reads flow immediately; nothing waits on the genome.
3. **The caller cross-checks after joining the verification handle**, at the same point and for the
   same reason it already joins for the `.fai` check: before committing output. This module exposes the
   captured tags and a pure comparison for that (§3.4); it does not own the sequencing.

This buys a check production cannot make — production trusts the `@SQ M5` — at no startup cost, in
exchange for detecting the fault *late*: a wrong-assembly run wastes its work before aborting. That is
the right trade for something catastrophic and rare, and it is the same fail-late-but-before-commit
shape `reference_info` already uses for the `.fai` itself.

**No `M5` is not a problem — the check simply does not apply.** A contig is compared only when *both*
sides carry a digest. Many BAM/CRAMs have no `M5` in their `@SQ` lines at all, and a `.fai`-only
reference has none either. **This is never an error and never a warning** (owner, 2026-07-20): a file
without `M5` is ordinary input, and refusing or nagging would punish the common case for missing an
opportunistic extra check. The run proceeds exactly as it would have without the feature.

`AssemblyCheck { compared, total }` is still returned, as *information* rather than a complaint — it
costs nothing, and it lets a caller that wants to say "assembly verified for 18 of 24 contigs" do so.
A caller that ignores it is behaving correctly.

**This is the one true inversion from production.** Production builds the *canonical* contig table
from the **first file's `@SQ`** and validates the FASTA and the other files against *that*
(`load_pileup_inputs`, `validate_fasta_agreement`). ng builds the canonical table from the
**reference** (`reference_info`, the sole builder — reference_info.md §1 mandate) and validates every
file against it. The reference is the authority; a file that disagrees is the thing that is wrong.

### 3.2 The region query

Given a validated handle and a `(contig, start, end)` (1-based inclusive), yield every read whose
reference footprint overlaps the interval, **in coordinate order**. Mechanically this is production's
`segment_reader.rs` design, re-expressed for ng:

- **BAM:** resolve `contig`→`ref_id`, query the index for candidate chunks, seek each,
  `read_record_buf` in order, drop records that miss the interval (chunk-edge slop — uncounted, a
  reader concern not a filter drop), and **early-stop** once a record on the target contig starts past
  `end` (the sort guarantees nothing later overlaps).
- **CRAM:** walk the `.crai` for containers on the target contig, seek + decode each container's
  slices to `RecordBuf`s, drop non-overlapping, and early-stop once a container *starts* past `end`.

**The chain, and its two different seams.** The per-file chain has three stages but only *two* of the
joins are iterator joins, and the distinction is load-bearing rather than pedantic:

```
region source ══(fills one reused RawRecord buffer)══▶ ReadFilter ──▶ order-verify ──▶ Result<MappedRead, _>
      ↑                                                     ↑                              ↑
  no MappedRead exists yet          decode() runs here, for survivors only     uniform item type from here on
```

- **The reader→filter seam passes an undecoded buffer, not a product.** `RecordSource::read_next` fills
  a caller-owned `RawRecord` in place, reusing its allocations
  ([`filtering.rs`](../../../../src/ng/read/filtering.rs) `RecordSource`); `RawRecord::flag`/`mapq` are
  cheap field reads on that undecoded record, and `RawRecord::decode` — the call that actually builds a
  `MappedRead` — is invoked by the filter **only for reads that survive the cheap cascade**. So the
  whole pass allocates one record buffer, and a dropped read is never materialised. This is also why
  `RecordSource` is a fill-a-buffer trait rather than an `Iterator`: an iterator must yield something
  owned, which would force either the decode or a borrow of the reused buffer that the type system
  rejects (the lending-iterator problem, `read_filtering.md` §5).
- **From the filter onward the item type is uniform** — `Result<MappedRead, _>` through order-verify
  and, one layer up, the merge. *This* is the composability the sample layer relies on.

This chain is deliberately the *complete* product of this module: the sample layer either uses it as-is
(one file) or merges k of them (`sample_reads.md` §3.2), and this module is indifferent to which.

**How it composes with ng.** The region source produces **raw, in-region records** and hands them to
ng's existing `RecordSource`/`RawRecord` seam so ng's **`ReadFilter`** (unchanged, all nine filters
including the reference-dependent #8) does the filtering and decode. Concretely: two new region-query
`RecordSource` implementations — the index-driven siblings of the whole-file `BamRecordSource` /
`CramRecordSource` that already live in [`read/filtering.rs`](../../../../src/ng/read/filtering.rs) —
each doing only the seek + overlap + early-stop, then `ReadFilter` wraps one and yields
`Result<MappedRead, _>` in region order. So the region reader is *reader logic only*; the *policy*
(what to drop) stays in step 1, one place, as `read_filtering.md` §2.5 requires.

**Filter here, not above.** Filtering is part of this chain rather than something the sample layer
applies to merged output. Merging a read only to drop it is wasted work, and doing it here means a
one-file sample and a k-file sample behave identically by construction.

**And filter here, not *inside* the reader — the seam already makes that unnecessary.** The natural
worry is that keeping `ReadFilter` a separate stage means building a `MappedRead` and then throwing it
away, and that pushing the predicates down into the BAM/CRAM read loop would avoid it. It would not,
because the reader is not an upstream producer of finished reads: **the filter owns the loop and calls
into the reader**, and decoding happens inside the filter:

```
loop {
    source.read_next(&mut buf)?;    // refills ONE reused buffer — no MappedRead
    …cheap flag / mapq checks on the undecoded buf…   // a drop here is counted and never decoded
    let read = buf.decode()?;       // ← the MappedRead is built HERE, survivors only
    …filters #7–#9, which genuinely need sequence/CIGAR/reference…
    yield read
}
```

So both costs are already avoided — no per-read allocation, and no `MappedRead` for a dropped read —
*without* moving policy into the reader, which would break the one-place-for-thresholds rule
(`read_filtering.md` §2.5) and split the drop tallies away from where they are kept.

**The residual, for a measured pass (§8).** `NoodlesRawRecord` wraps a noodles `RecordBuf`, and BAM's
`read_record_buf` fully decodes each record — name, CIGAR, sequence, qualities into the buffer's
`Vec`s — *before* `flag()`/`mapq()` are consulted. The allocations are reused, so this is decode work,
not allocation churn, but it is spent on records that will be dropped. noodles' lazy `bam::Record` (a
byte-slice view with cheap `flags()`/`mapping_quality()`) could filter before materialising. Two
reasons it is not done now: it is **BAM-only** (CRAM decodes at container granularity regardless), and
its value scales with the drop rate — negligible if most reads pass, real on heavily
duplicate-marked data. Measure before building, and if built, keep the thresholds in `ReadFilter` and
the counting where it is.

**Order verification — the guarantee, not an assumption.** Everything downstream of this module is
built on reads arriving in genome order, so the module **proves it while streaming** rather than
trusting the header. A thin adapter over the filtered stream keeps the last emitted `(ref_id, pos)`
and, for each read it is about to yield, hard-errors `OutOfOrderRead { file, previous, current }` if
the new key is less than the last. There is no tolerance, no warning-and-continue, and no silent
re-sort: an unsorted file is a fatal input error the user must see.

Four points that fix the exact semantics:

- **The key is `(ref_id, pos)`, so the check covers both axes.** A position regression within a contig
  and a contig-order regression (a read on contig 5 after a read on contig 7) are the same violation.
  Because the open gate has already proved `ref_id == ContigId` (§3.1), "contig order" here means the
  reference's order, which is the order everything else in ng means.
- **Equal keys are legal.** Several reads may start at the same position; the check rejects only a
  strict decrease.
- **It sits on the *filtered* stream.** Dropped reads (unmapped, low-MAPQ) cannot break the
  monotonicity of what survives, and checking after the filter means the guarantee is about exactly the
  reads the caller sees.
- **The check is per region stream, and deliberately does not span queries.** The last-key state lives
  in the iterator, not on the handle. A caller is entitled to query region B and then region A — that
  is a new forward scan, not a regression (§1, non-goals). Carrying the state on the handle would turn
  legitimate random access into a spurious error.

`@HD SO` (§3.1 check 1) and this check are complementary, not redundant: the header check is cheap and
fails at `open` before any work, while this one catches the file that *claims* to be sorted and is not
— the case the header check structurally cannot see.

**This is the single authority for read order in ng.** The merge one layer up trusts it rather than
re-checking (production checks in the merge because its readers do not; ng flips the check to where
the file is read — see [`sample_reads.md`](sample_reads.md) §2, precondition 2). Within a single
indexed region the scan is usually monotonic already; the check earns its keep across chunk and
container boundaries, and as the guard that a mis-sorted file fails loudly rather than feeding a
silently-wrong merge.

**Counts.** The chain reports this file's `ReadFilterCounts` (`read_filtering.md` §4) — one tally, for
one file. Aggregation across files is refused on purpose; see `sample_reads.md` §3.3 for why.

### 3.3 Cost: open once, seek per query

**This is the baseline, not an optimisation.** The STR path issues on the order of 10⁶ region queries
(one per locus), and a sample multiplies that by its file count, so anything paid *per query* is paid
~10⁶·k times. The rule that makes it affordable: **the opened reader and the parsed index are owned by
the handle and live for the whole run**; a query is then an in-memory index lookup plus one seek, and
never a file open or an index parse.

Production establishes exactly this and says so — `load_alignment_index` runs once per input at startup
(`index_preflight.rs:120`, from `alignment_input.rs:703`), the `Arc`-wrapped index is moved into
`BamFile`/`CramFile` (`segment_reader.rs:457`, `:679`), and the module doc (`segment_reader.rs:9-17`)
frames the whole design as removing the per-call re-open + header re-parse an earlier version paid. ng
does the same; the index load is check 3 of the open gate (§3.1), and T13 asserts it.

**Drive the index and the reader separately, as production does.** Production does not use noodles'
`IndexedReader`/`Reader::query` at all. It calls `BinningIndex::query(target, interval)` on the
already-parsed index to get chunks (`segment_reader.rs:485`), then seeks its own open reader per chunk
and linear-scans with the sorted early-stop (`:594`, `:622`). ng should follow this rather than the
convenience type, because it keeps the index-once / seek-per-query split explicit in the code instead
of hidden inside noodles.

**One production inefficiency not to inherit.** `fetch_mapped_reads` rescans the `.crai` from entry 0
on every call (`segment_reader.rs:1043`) — `continue` past off-contig entries, `break` only once past
the segment end. For a late contig in a many-contig `.crai` over 10⁶ loci that is a repeated O(n)
prefix scan per locus. The `.crai` is sorted, so ng's CRAM region source **binary-searches or carries a
cursor** across queries (production's own `CramSegmentReads::next_index_record` already carries one;
the caching path does not).

**Reader pooling is built now, not deferred (owner, 2026-07-20)** — because it is an **API-shape
decision, not an optimisation.** ng will be parallelised; the expensive thing to retrofit is not the
pool but the borrowing shape around it. If `reads_in_region` takes `&mut self`, concurrent use is
impossible without redesigning every caller; if it takes `&self` and borrows a handle from an
internally-shared pool, parallelism drops in later with no signature churn. Since `&self` requires
interior mutability anyway, **the pool is the mechanism that makes the right signature possible** —
building it later would mean writing the wrong signature first.

Production's shape is the model: `Mutex<Vec<Handle>>` on the file struct (`segment_reader.rs:463`,
`:687`), `borrow_handle` pops or opens fresh (`:508-514`), return happens in the iterator's `Drop`
(`:665-669`), the mutex guards only the `Vec` pop/push and is never held during iteration
(`:521-530`). Crucially the **index and header are not pooled** — they stay `Arc`-shared on the file
struct, so no pooled reader ever re-parses them (§3.3's guarantee survives pooling unchanged).

**Still deferred** (§8, measured pass): the CRAM decoded-container cache (`CachingCramReader` — FIFO
over whole decoded containers, capacity 3, on the reasoning that a locus window overlaps at most two,
`segment_reader.rs:946`/`:996`). Note it **interacts with pooling** and does not simply layer on:
production gives each worker its own non-pooled `CachingCramReader` (`segment_reader.rs:1121-1123`,
"One per worker (not pooled)") precisely because a FIFO container cache shared across workers would
thrash. So the CRAM path's eventual shape is per-worker caches over a shared index, not a pooled
handle — worth knowing now so the pool's design does not assume otherwise.

**A sorted batch of loci could later be served by one forward sweep** instead of N seeks — production's
SNP path already works this way (one index query per contig, one streaming scan, `cli.rs:463-476`),
while its STR path does not (one query per locus window, `fetch_reads.rs:176-185`). The layered-stream
design makes a sweep just another region source, so the interface must not preclude it; but it is **not
built now** — the container cache recovers most of the CRAM cost and BAM seeks are cheap (§8).

### 3.4 The entry point

```
AlignmentFile::open(path, reference_info, filter_config, index_policy) -> Result<AlignmentFile, AlignmentFileError>
AlignmentFile::sample_name() -> &str                       // from @RG SM, §3.1 check 4
AlignmentFile::sq_md5s() -> &[Option<Md5>]                 // @SQ M5 tags, by ContigId, §3.1
AlignmentFile::reads_in_region(&self, contig, start, end) -> Result<RegionReads<'_>, AlignmentFileError>
     where RegionReads: Iterator<Item = Result<MappedRead, AlignmentFileError>>   // ordered, filtered

// the deferred assembly check — a pure comparison the caller runs after joining
// reference_info's VerificationHandle, before committing output (§3.1)
check_assembly(observed: &[Option<Md5>], verified: &ReferenceInfo) -> Result<AssemblyCheck, AssemblyMismatch>
     where AssemblyCheck { compared: usize, total: usize }   // report, never silently "fine"
```

`open` runs the whole gate (§3.1) and leaves a handle owning the reader and the parsed index;
`reads_in_region` builds the region source → filter → order-verify chain (§3.2). **It takes `&self`,
not `&mut self`** — the load-bearing signature choice (§3.3): a handle comes from the internal pool,
and the returned iterator returns it on `Drop`. That is what will let N threads query one
`AlignmentFile` concurrently without touching a single call site. `check_assembly` is a
free function, not a method: it must be callable once the file handles are long gone, and keeping it
pure makes it trivially testable. The `ReadFilter`
needs a reference accessor for #8 — the same `RawRefSeq` the rest of ng uses. (Exact ownership — one
`RawRefSeq` shared across files vs one per file, and how it derives from the `ReferenceInfo` — is an
arch concern, and it is the one place this module's shape is influenced by there being several files.)

---

## 4. Error model

Two outcomes, the `read_filtering.md`/`reference_info.md` house pattern:

- **A read is filtered** — tallied in `ReadFilterCounts`, never an error, never silent. Step 1's job.
- **A run-level failure is fatal** — surfaced as an `Err`, never swallowed into a drop. Two moments:
  *at `open`* (index missing or unparseable, `SO != coordinate`, contig reconciliation fails, the file
  names two samples, the file cannot be opened) — returned before any read flows; and *mid-stream* (a
  truncated file, a decode failure, an out-of-order read, a #8 reference-fetch failure) — yielded once
  in the item stream, then the iterator **fuses**: the first `Err` is yielded once (`Some(Err(_))`) and
  then `None`, so `for read in reads { let read = read?; … }` surfaces it and cannot mistake it for
  clean EOF (the `ReadFilter` convention).

A sketch of the error enum (`#[non_exhaustive]` thiserror, each variant carrying its path): `Index`
(preflight/load), `NotCoordinateSorted`, `ContigReconcile { file, detail }`, `MultipleSampleNames`,
`Open`, `OutOfOrderRead { file, previous, current }` (§3.2 — carries both keys so the message can name
where the file breaks), `Filter(ReadFilterError)` (wrapping step 1's source/decode/reference arms),
and `Region` (an invalid `(start,end)` or a contig absent from the reference). `AssemblyMismatch
{ file, contig, expected, observed }` is raised by `check_assembly` (§3.1), which runs long after the
stream is done — so it is a *separate* error type returned from that function, not a variant that can
appear in the item stream. Exact shape → arch.
`sample_reads.md` adds its own cross-file variants on top rather than extending this enum (arch call).

---

## 5. Reuse map — what ng calls, what ng owns

The ng rule established by `read_filtering.md`: **reuse the pure helpers, own the driver.** `src/bam/`
is read-side production code the ng freeze permits reusing (the freeze is about not *editing*
production and not depending on trf-mod — read_filtering already reuses `bam/alignment_input` heavily).

| concern | reuse (call as-is) | ng owns (new, modelled on production) |
|---|---|---|
| index policy | `index_preflight::{preflight_alignment_indexes, load_alignment_index, AlignmentIndex}` (`pub`) | — |
| contig compare | `ContigList::first_disagreement` (`pub(crate)`) | the reconcile call + the `NotCoordinateSorted`/`@SQ`-extraction glue |
| canonical table | `reference_info::ReferenceInfo` / `contig_list()` | consuming it as the authority (§3.1 inversion) |
| CRAM reference | `fasta::Repository` (build via noodles, as ng's `CramRecordSource` already does) | — |
| sample name | ~~`alignment_input::extract_single_sample_name`~~ — **module-private** (`alignment_input.rs:378`), so not ng-reachable; `extract_header` (`:292`) likewise | ng reads `@HD SO` / `@SQ` / `@RG SM` off the noodles `sam::Header` itself (arch §5) |
| region read | noodles `BinningIndex::query` on the pre-parsed index + own seek/scan (**not** `IndexedReader`/`Reader::query`, §3.3) + `decode_cram_container` shape | the region-query `RecordSource` impls, incl. the `.crai` cursor |
| filtering | ng `ReadFilter` / `MappedRead` (already ng-owned) | composing it into the per-file chain |

**Decided (owner, 2026-07-20): rebuild ng-owned, modelled on production — position (B).** Production's
`segment_reader::AlignmentFile` is `pub(crate)` (ng-reachable) and *already* does region-query + order
guard yielding the same `MappedRead`, and already gets the index-once/seek-per-query structure right
(§3.3), so full reuse (A) or a hybrid raw-record adapter (C) were both live options. (B) wins for the
same reasons `read_filtering.rs` and the existing ng `Bam/CramRecordSource` were built that way:

- **The full step-1 policy runs in one place.** Production's readers apply only the cheap filter subset,
  so (A) would need a #8 wrapper bolted on and ng's filter policy would live in two places — the thing
  `read_filtering.md` §2.5 exists to prevent. Under (B) `ReadFilter` wraps the region source and all
  nine filters run naturally (T9).
- **ng keeps its own contig authority.** (A) inherits production's `@SQ`-derived contig table, so the
  §3.1 gate would have to front a component that internally believes something else. (B) has one
  canonical table, the reference's, end to end.
- **It extends a seam ng already owns.** The region readers become index-driven siblings of the
  existing whole-file `Bam/CramRecordSource`, reusing the `RawRecord` buffer protocol (§3.2) rather
  than adapting a foreign iterator into it.

**What this costs:** the most code of the three options, and ng re-derives the seek/early-stop logic
production already has. Mitigated by reusing every *pure* helper (the table above) and by production's
`segment_reader` remaining the model to copy — including its index-once structure (§3.3) and excluding
its `.crai` prefix-rescan (§3.3).

This decision is independent of the *merge* reuse decision (`sample_reads.md` O3-merge, also rebuild);
the specs were split partly so the two could be answered separately.

---

## 6. Module layout (to confirm in the arch doc)

Home (owner, 2026-07-20): a `read/input/` submodule beside `read/filtering.rs`, reusing that file's
`RecordSource`/`RawRecord` seam. This spec owns:

```
src/ng/read/input/
├── open_bam.rs     – the open-and-validate gate (SO, @SQ, index load, @RG SM) (§3.1)
└── region_query.rs – the BAM/CRAM region-query RecordSource impls (§3.2), incl. the .crai cursor
```

`sample_reads.md` §6 owns the rest of the directory. `read/` already houses steps 1+2
(`module_layout.md` principle 1, note b); this is step 1's input edge, so it belongs here rather than
as a top-level sibling.

*Naming note, for the arch doc to settle or overrule:* `open_bam.rs` also opens CRAM. The owner chose
it as the clearer name, on the everyday reading where "BAM" stands for the alignment file. The
alternative, if the mismatch grates later, is `open_alignment_file.rs`, which matches the type
(`AlignmentFile`) and this spec's own name — at the cost of length inside a directory already called
`input`.

---

## 7. Test obligations (T-list — the acceptance bar)

- **T1 — permutation is caught.** A BAM whose `@SQ` is a permutation of the reference → `open` errors
  `ContigReconcile` naming the first transposed contig. (The hole reference_info.md §1 named;
  mutation-verify against a *resolves-only* check, which passes the permutation.)
- **T2 — length/name/count mismatch caught** at `open`, each naming the right field/index; the digest
  comparison at `open` fires **only when the `ReferenceInfo` carries MD5s** (Fasta arm), a no-op under
  a `.fai`-only table.
- **T2b — the deferred assembly check** (§3.1). A BAM whose `@SQ M5` disagrees with the verified
  reference's per-contig MD5 → `check_assembly` errors `AssemblyMismatch` naming the contig and both
  digests, *after* reads have flowed (proving the check is not at open and does not block startup). A
  BAM with **no `M5` tags** → no error, and `AssemblyCheck { compared: 0, total: n }` so the caller can
  report that nothing was verified. Mixed (some contigs tagged) → `compared` counts exactly the tagged
  ones.
- **T3 — `SO != coordinate` (and missing) → `NotCoordinateSorted` at open**, before any read.
- **T4 — a genuinely out-of-order file errors mid-stream**, in all four of its aspects (§3.2):
  - **T4a** — `SO:coordinate` in the header but a planted *position* regression in the data →
    `OutOfOrderRead` naming both keys, fused after. (The header check cannot catch this; mutation-verify
    that removing the streaming check lets the file through.)
  - **T4b** — a planted *contig-order* regression (a read on an earlier `ref_id` after a later one) →
    the same error.
  - **T4c** — several reads sharing one start position are **not** an error (equal keys are legal).
  - **T4d** — querying region B and then region A on the same handle is **not** an error (the check is
    per region stream, not per handle).
- **T5 — region query** returns exactly the overlapping reads, coordinate-ordered, with chunk-edge slop
  dropped and *uncounted*; the BAM early-stop and CRAM container-stop exercised (a read well past `end`
  is not scanned).
- **T8 — BAM/CRAM parity.** The same reads written as BAM and as CRAM (same reference) produce the same
  ordered `MappedRead` sequence through `reads_in_region`.
- **T9 — the full step-1 filter runs.** A #8 high-mismatch read is dropped by the served stream, proving
  `ReadFilter` is composed in rather than just the cheap subset — the property that decided O3-file in
  favour of the rebuild (§5), so this test is now unconditional.
- **T10 — fatal-is-fatal, fused.** A truncated file mid-region → one `Err`, then `None`; no partial
  silent EOF.
- **T12a — a file whose `@RG` records name two samples → error at `open`** (the cross-file half is
  `sample_reads.md` T12b).
- **T13 — the index is parsed once.** Many region queries against one opened handle do not re-open the
  file or re-parse the index (assert via an instrumented reader or an open-count probe) — the guarantee
  §3.3 rests on, and the one that decides whether 10⁶ STR queries are affordable.

T-numbers are kept from the pre-split spec so review notes stay traceable; T6/T7/T11/T12b live in
`sample_reads.md`.

Fixtures: tiny in-memory BAM/CRAM built with noodles writers + an in-crate CSI/CRAI build (production's
`segment_reader` test module and `index_preflight` tests are the templates), over a small multi-contig
reference indexed by `reference_info`.

---

## 8. Deferred, with a home

- **CRAM decode-once container cache** → a measured performance pass; `CachingCramReader`
  (`segment_reader.rs:996`) is the reuse target (§3.3). Correctness first, then speed. *Note what is
  **not** deferred: opening the reader once and parsing the index once are baseline, not optimisation
  (§3.3, T13).*
- ~~Reader pooling~~ → **no longer deferred**; built now because it fixes the `&self` signature that
  parallelism will need (§3.3).
- **Pre-`RecordBuf` filtering for BAM** → applying the cheap flag/MAPQ predicates to noodles' lazy
  `bam::Record` before materialising a `RecordBuf`, saving decode work on doomed records (§3.2,
  "the residual"). BAM-only; value scales with the drop rate; measure first, and keep thresholds and
  counting in `ReadFilter` if it is built.
- **Batched-locus forward sweep** → serving a sorted batch of loci with one scan instead of N seeks,
  the way production's SNP path already works (`cli.rs:463-476`) and its STR path does not
  (`fetch_reads.rs:176-185`). Just another region source; build it only if measurement asks for it.
- **Build-if-missing index policy** → a caller flag threaded into `open` (`preflight_alignment_indexes`
  already takes `build_if_missing`); ng's default (error vs build) is a CLI decision, not this module's.
- **Whole-genome linear scan** (the un-indexed alternative the owner set aside) → a second
  `reads_in_region`-free entry point if a use case wants it; the existing whole-file
  `Bam/CramRecordSource` already implement the linear source.

---

## 9. Open questions

### Still open

*None. Every question this spec raised has been answered; remaining choices are arch-level (exact type
and error shapes, `open_bam.rs` vs `open_alignment_file.rs` — §6).*

### Resolved (owner, 2026-07-20)

- **O1 — MD5 in reconciliation** → foreground reconciliation is name+length+order; `@SQ M5` tags are
  captured at open and the digest comparison is deferred to the caller's join of `reference_info`'s
  background verification — non-blocking startup *and* wrong-assembly detection, at the price of
  catching it before output rather than before reads. Contigs lacking a digest on either side are
  skipped and counted, never silently passed (§3.1).
- **O3-file — rebuild vs reuse the region reader** → **(B) rebuild ng-owned** (§5), for
  one-place-for-policy, ng's own contig authority, and extending the `RecordSource` seam ng already
  owns. Independent of `sample_reads.md`'s O3-merge.
- **O5 — module home** → `src/ng/read/input/`, with `open_bam.rs` and `region_query.rs` (§6). Same
  answer applies to `sample_reads.md`.
- **Missing `@SQ M5` → not an error, not a warning** (§3.1). The check does not apply; the run
  proceeds normally.
- **Reader pooling → built now, not deferred** (§3.3). It is what makes `reads_in_region(&self)`
  possible, and that signature is what parallelism will need.
