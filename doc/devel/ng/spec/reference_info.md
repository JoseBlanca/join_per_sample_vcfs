# ng — reading a reference's info: contig table, digest, and index (`reference_info`)

*Status: design spec (2026-07-17, revised 2026-07-18). **No code yet.** The missing direction of
the reference stack: every contig table in the tree is built from an alignment file's `@SQ`
headers and the `.fai` is only ever validated *against* one. This module builds the table **from
the reference itself** — a FASTA (optionally checked against a `.fai`) or a `.fai` alone — and
returns a `ReferenceInfo`: the contig table, the content digest, and the reconstructed index.
Companion to [`ref_seq.md`](ref_seq.md), whose impls all take a `ContigList` they cannot build.
**Code-facing companion (types & interfaces):** [`../arch/reference_info.md`](../arch/reference_info.md).*

*Naming: **STR** in prose, `ssr` in code. This module is marker-agnostic, so `ssr` does not
appear.*

---

## 1. What it is

**A module that answers one question: what contigs does this reference have, in what order,
how long are they, where do they sit in the file, and — if we were shown the bases — what are
their MD5s?** In: one to two paths. Out: a `ReferenceInfo` — a `Vec<ContigInfo>` (name, length,
geometry, md5 per contig) plus, when a FASTA was read, a whole-reference digest (§2). It projects
to a production `ContigList` ([`fasta/mod.rs:59-62`](../../../../src/fasta/mod.rs)) for the seams
that need one.

Nothing in the tree does this. `extract_header` builds a `ContigList` from a `sam::Header`'s
`@SQ` records ([`alignment_input.rs:306-313`](../../../../src/bam/alignment_input.rs));
`validate_fasta_agreement` reads the `.fai` only to *compare* it against that table
([`alignment_input.rs:426-484`](../../../../src/bam/alignment_input.rs)) and throws it away.
So a consumer holding a reference and no alignment file has no table, and ng has two:

1. **`typed-regions`** has no alignment file at all
   ([`typed_regions_cli.md`](typed_regions_cli.md) T1) — `WindowedRefSeq::new` demands a
   `ContigList` ([`ref_seq.rs:524`](../../../../src/ng/ref_seq.rs)) and there is nowhere to
   get one.
2. **Read ingestion**, where the hole has teeth. `ContigId` is documented as an index into
   the table *in `.fai` order* ([`types.rs:6-8`](../../../../src/ng/types.rs)), and
   `ReadFilter::new` ([`filtering.rs:602-614`](../../../../src/ng/read/filtering.rs)) assumes
   `ContigId(i)` is the i-th `@SQ` — but it checks only that each index *resolves*, with a
   zero-length fetch. **A BAM whose `@SQ` list is a permutation of the `.fai` passes, then
   fetches the wrong contig for every read in filter #8.** Only a missing contig is caught
   (`read_filter_new_rejects_a_contig_missing_from_the_reference`). Closing that needs two
   tables compared; this module builds the one that does not exist yet. The comparison itself
   is **not** here — §7.

**Once this exists, it is the *only* way ng builds a `ContigList` from a reference file**
(owner, 2026-07-17). Every ng consumer that has a reference and needs its contig table —
`typed-regions`, read ingestion, and `WindowedRefSeq::new` / `ResidentRefSeq::new`, which take
a `ContigList` they cannot build ([`ref_seq.rs:367`, `:524`](../../../../src/ng/ref_seq.rs)) —
calls this module. No consumer hand-rolls `ContigEntry`s from `.fai` records or `@SQ` headers
again; a second builder is a second place for the order contract (T1), the duplicate check
(T2) and the M5 rule (§3.4) to drift. The rule has **one exception, and it is not one**:
`InMemoryRefSeq::from_named_contigs` ([`ref_seq.rs:271`](../../../../src/ng/ref_seq.rs)) builds
its table from bytes it already holds in memory, not from a reference *file*, so it is a
different operation, not a bypass. The two current hand-rolled sites are **test fixtures**
(`ref_seq.rs`'s `build_fasta`, [`anchor.rs:79`](../../../../src/ng/region_typing/anchor.rs));
migrating them onto this module is optional cleanup that would buy them real coverage, but the
mandate binds *non-test* consumers, of which the first is not yet built.

**And the mandate is wider than the table: this is the *only* place ng parses a `.fai`** (owner,
2026-07-18). The reference stack's fai-driven readers — ng's `RawChromReader`, and through it
`WindowedRefSeq` — need per-contig geometry (`offset`, `line_bases`, `line_width`) to seek, and
that geometry is exactly what `ContigInfo` carries (§2). Today they parse the `.fai` themselves,
redundantly and *repeatedly*: `WindowedRefSeq::new` is handed a contig table (one parse), then
`RawChromReader::for_contig` re-reads the whole `.fai` on **every contig transition**
([`raw_chrom_reader.rs:131`](../../../../src/ng/raw_chrom_reader.rs)). The consolidation: this
module parses the `.fai` **once** (cached §3.7, validated §3.8, all the geometry guards in one
place), producing the full `Vec<ContigInfo>`, and the readers take a `&ContigInfo` for geometry
instead of re-parsing. The split stays clean — **this module owns `.fai` *parsing*;
the readers own base *reading*** over the geometry it hands them (reading bases is a non-goal
here, `RefSeq`'s job) — so the dependency is one-way (readers → `reference_info`) with no cycle,
since §4's own pass hand-rolls its loop and calls no reader. Threads share one immutable
`Arc<ReferenceInfo>`; each worker keeps its own windowed buffer over it.

**This consolidation is ng-internal, and that is the frozen-production rule, not a shortfall.**
Production's `crate::fasta::StreamingChromRefFetcher` / `ManualEvictChromRefFetcher` **cannot**
call this module — production must not depend on ng (the B2 guard, [`typed_regions.md`](typed_regions.md) §9) — so
they keep their own `.fai` reading, frozen. "The only place that parses a `.fai`" is therefore a
statement about *ng*, which is the only scope that matters here. The concrete first target is
`RawChromReader` taking a `&ContigInfo` rather than calling `open_contig`/`fai::fs::read`; it is a
follow-up refactor of existing ng code, not part of this module's build, and it is recorded in §7.

**Goals.** Build a `ReferenceInfo` (contig table + digest + index, §2) from a reference — from a
FASTA (optionally with a `.fai`) or from a `.fai` alone (§3.1). When the FASTA is read, compute each contig's MD5 and the whole
reference's, in **one streaming pass**, holding one buffer rather than one contig (§4). When a
`.fai` is supplied with the FASTA, **prove they describe the same genome** — completely enough
that a stale `.fai` cannot survive (§3.3). Carry the reconstructed index the FASTA pass computes
anyway, and **write a `.fai`** from it, so ng can index a reference that lacks one (§3.9). Offer
a **caller-held, thread-safe cache** so a cohort's parallel readers pay the genome read once and
wait on it, not N times (§3.7). Offer an opt-in **verify-in-the-background** entry point so a
caller gets the `.fai` info now and the genome-read verification runs off-thread, its error
surfacing on join (§3.10). And offer **`read_reference_verifying_or_creating_fai`**, a batteries-included entry point that
takes a FASTA and does the right thing — background-verify if the `.fai` is there, scan-and-index
if it is not (§3.11). Stay a leaf: this module knows `crate::fasta` and noodles, and nothing else
about ng (§3.6).

**Non-goals.** **Reading bases for anyone.** That is `RefSeq`'s job ([`ref_seq.md`](ref_seq.md));
this module reads the FASTA once, to describe it, and hands back no sequence bytes (it does hand
back the reconstructed index — §3.9 — but never the bases themselves). **Reconciling against an alignment file's
`@SQ`**, and the coordinate-sort check that stands on it (§7). **Narrowing anything to `u32`**
(T5). **Parallelising a single read** — the streaming pass is serial by nature (§3.3 rejects
the offset-seeking that would let it fan out); the cache makes *many* reads cheap, it does not
split *one* (§7's parallelism note).

---

## 2. The interface

The whole surface, and deliberately small — a caller names its files and gets a table:

```rust
/// What the caller hands in. Two ways to describe a reference (§3.1); the variant names the
/// cost — `Fai` is a small file read, `Fasta` reads the whole reference (§3.2). Owns its
/// paths (`PathBuf`) to stay lifetime-free through the reader and the cache; a `PathBuf` clone
/// next to a genome read does not register.
pub enum ReferenceSource {
    /// The index alone. Names, order and lengths; no MD5s.
    Fai(PathBuf),
    /// The FASTA, **optionally** cross-checked against a `.fai`. `fai: None` reads the FASTA
    /// alone (names, order, lengths, every MD5, the reference digest — §3.4). `fai: Some(p)`
    /// reads the same **and** verifies the FASTA against that index (§3.3); a disagreement is
    /// an error. Both cost a whole-genome read — the index only adds a check (§4).
    Fasta { fasta: PathBuf, fai: Option<PathBuf> },
}

/// **Everything about one contig, in one place.** Its name and length (what a BAM `@SQ`
/// carries), the `.fai` geometry that locates it in the file (§3.8), and — once the bases are
/// read — its MD5 (the `@SQ M5`, §3.4). Named for what it *is*, not where a field came from:
/// some columns are read from the `.fai`, but that is circumstantial; the contig's info is the
/// point. This is ng's own type; it is a *superset* of any production type (`ContigEntry` has
/// name/length/md5 but no geometry; a `.fai` row has geometry but no md5 — `ContigInfo` unifies
/// what production splits).
pub struct ContigInfo {
    pub name: String,
    pub length: u64,
    /// Byte offset of the first base (`faidx.5` OFFSET).
    pub offset: u64,
    /// Bases per line (`faidx.5` LINEBASES).
    pub line_bases: u64,
    /// Bytes per line incl. newline (`faidx.5` LINEWIDTH).
    pub line_width: u64,
    /// The `@SQ M5` (§3.4). `None` until the FASTA is read — a `.fai` has no MD5 column.
    pub md5: Option<[u8; 16]>,
}

/// What reading a reference yields: the whole-assembly digest and every contig's info.
pub struct ReferenceInfo {
    /// The **whole-reference** MD5 — every contig's uppercased bases *concatenated* in file
    /// order (§3.4; the `ssr-catalog` convention, and **not** a SAM field — SAM has no
    /// reference-wide digest). Distinct from the per-contig `@SQ M5`, which lives on each
    /// `ContigInfo`. `None` from a `.fai`-only read.
    pub md5: Option<[u8; 16]>,
    /// Every contig, **in file order** — the order that defines `ContigId` (T1). One source of
    /// truth: name, length, geometry and md5 all here, so there is nothing to keep in step
    /// (this is the merge that retired the old `ContigList` + `fai_records` split, §3.9).
    pub contigs: Vec<ContigInfo>,
}

impl ReferenceInfo {
    /// **Transitional** bridge to the production `ContigList` (one `ContigEntry { name, length,
    /// md5 }` per contig) for ng's own consumers that still speak it today — the `RefSeq` impls,
    /// `ContigTable`, `GenomeRegions` validation. Those migrate to `ContigInfo` (§3.6); when they
    /// have, this goes away. A fresh allocation; the geometry has no place in a `ContigList` and
    /// is dropped. The *only* place `ContigList` surfaces from this module.
    pub fn contig_list(&self) -> ContigList;
}

/// Read a reference. **Pure** — no shared state, no memoisation, no filesystem probing beyond
/// opening the paths it was handed; every call does the full work the source names (§3.1). A
/// single-threaded, single-call consumer calls this directly. A consumer that reads the same
/// reference more than once, or from several threads at once, holds a [`ReferenceInfoCache`]
/// instead (§3.7).
pub fn read_reference_info(source: ReferenceSource)
    -> Result<ReferenceInfo, ReferenceInfoError>;

/// Write a `.fai` from contig info — the five-column `faidx.5` form (§3.8), `md5` ignored (a
/// `.fai` has no MD5 column), byte-identical to `samtools faidx`. Atomic (`.tmp` + rename),
/// because a truncated `.fai` is silently valid and every later reader would trust it (§3.9).
pub fn write_fai(contigs: &[ContigInfo], out: &Path) -> io::Result<()>;

/// Read the `.fai` **now** and verify the FASTA against it on a **background thread** (§3.10),
/// with **both reads going through the cache**. The returned `ReferenceInfo` is the `.fai` read
/// (`get_or_read(Fai)`) — names, lengths, order, geometry, `md5: None` — available immediately,
/// so the caller need not block on the genome read. The initial `.fai` read is synchronous, so a
/// missing or malformed `.fai` errors here, before any thread is spawned. The background thread
/// runs `get_or_read(Fasta { fai: Some })` — which **verifies** (§3.3), computes the MD5s (§3.4),
/// **and populates the cache**: a later `get_or_read` of that key is a hit, and a concurrent one
/// *while the verify runs* waits on the same single-flight slot instead of reading again (§3.10).
/// Takes `&Arc<ReferenceInfoCache>` because the background thread outlives the call and holds a
/// clone. The `VerificationHandle` carries the background outcome.
pub fn read_fai_verify_in_background(
    cache: &Arc<ReferenceInfoCache>,
    fasta: PathBuf,
    fai: PathBuf,
) -> Result<(Arc<ReferenceInfo>, VerificationHandle), ReferenceInfoError>;

/// The pending result of the background FASTA verification. **`#[must_use]`, and it must be
/// joined before the caller commits any output** — an unverified reference is a possibly-wrong
/// run, and this handle is the *only* place that error surfaces (§3.10).
#[must_use]
pub struct VerificationHandle { /* a crossbeam result channel + the worker's join handle */ }

impl VerificationHandle {
    /// Non-blocking: `true` once verification has finished (so a caller can defer `join` until
    /// it will not block). Does not consume the result.
    pub fn is_finished(&self) -> bool;
    /// Block until verification finishes and take its result — the **verified** `ReferenceInfo`
    /// (now carrying the MD5s the `.fai` could not) on success, the geometry or read error on
    /// failure.
    pub fn join(self) -> Result<Arc<ReferenceInfo>, ReferenceInfoError>;
}
// Drop without join() abandons the run's only verification: it logs a warning to stderr (the
// house convention) rather than swallow a possibly-pending failure silently. It does not block.

/// **The batteries-included entry point (§3.11): read a reference the fast, convenient way,
/// verifying its `.fai` if present or creating one if not.** Derives the sibling `<fasta>.fai`
/// and does the right thing:
/// - **`.fai` present** → read it now and verify the FASTA in the background
///   (`read_fai_verify_in_background`). Returns immediately, before the scan finishes, with
///   `md5: None` info and `Some(handle)` the caller must `join` before committing output.
/// - **`.fai` absent** → scan the FASTA now (verified, MD5s), **write** the `.fai` beside it,
///   cache, and return `(verified info, None)` — nothing to join.
///
/// This is the *one* place the module adopts the `<fasta>.fai` naming convention (an exception
/// to §3.6) and the *one* place a read can write a file. A `.fai`-write failure is **fatal**
/// (owner's call, §3.11): a caller on a read-only reference dir uses the primitives instead
/// (`cache.get_or_read(Fasta { fai: None })`, no write).
pub fn read_reference_verifying_or_creating_fai(cache: &Arc<ReferenceInfoCache>, fasta: PathBuf)
    -> Result<(Arc<ReferenceInfo>, Option<VerificationHandle>), ReferenceInfoError>;
```

`ReferenceSource` is a noun, per the house rule.

The module's **second** public item is the optional cache — a caller-held object, not a
global, so the reader above stays a pure function and the cache's lifetime is a thing the
caller can see (§3.7):

```rust
/// A thread-safe, compute-once cache over `read_reference_info`. The caller constructs one
/// (normally one per run), shares `&cache` across its worker threads, and every thread that
/// asks for a reference already being read **waits** for that read rather than starting its
/// own (§3.7). Sharable by `&`: `Send + Sync`.
pub struct ReferenceInfoCache { /* private */ }

impl ReferenceInfoCache {
    pub fn new() -> Self;

    /// The result is `Arc`-shared, so a hit is a pointer clone, not a table copy. Keyed on
    /// the source **and** the file stats (§3.7) — a `Fai` and a `Fasta` of the same file are
    /// distinct entries (different tables), and so are `Fasta { fai: None }` and
    /// `Fasta { fai: Some }` (only the second verified).
    pub fn get_or_read(&self, source: ReferenceSource)
        -> Result<Arc<ReferenceInfo>, ReferenceInfoError>;
}
```

**No enum constructors, no filesystem probing — the module reads the paths it is given.**
This is the modularity line (§3.6): **checking a FASTA against a `.fai` is contig-list work and
lives here; discovering *which* `.fai` — e.g. the sibling `ref.fa.fai` by naming convention —
is a caller's policy and does not.** So the caller builds a `ReferenceSource` by picking the
variant and passing the paths; there is nothing to construct through, and the reader never
stats a file it was not handed.

The one place the two meet is the sibling naming convention, and it is offered as a **pure path
function** the caller may use — never as behaviour the module performs on its own:

```rust
/// `<fasta>` + ".fai" — the reference naming convention, as a path, no I/O. A caller that
/// wants sibling behaviour composes it; the module does not apply it. (ng copies the five
/// lines of `bam/alignment_input.rs`'s private `with_fai_extension` — the standing rule,
/// `ng/mod.rs:11-15`.)
pub fn sibling_fai_path(fasta: &Path) -> PathBuf;
```

A caller that wants "read the FASTA, and check the sibling `.fai` when one is there" writes the
policy in the open, where its filesystem probe is visible:

```rust
let fai = sibling_fai_path(&fasta);
let source = ReferenceSource::Fasta {
    fasta,
    fai: fai.exists().then_some(fai),   // present → checked (§3.3); absent → read alone
};
```

The check still fires for free whenever a `.fai` is passed (the genome read is already
committed, so the index only adds a comparison), and a **stale** sibling still turns a would-be
success into an error — that value is intact. What moved out is the *decision to go looking*,
which is a convention about file layout, not a fact about the reference's contigs.

*Rejected: a `from_fasta_checking_sibling_fai` constructor that probes for the sibling itself.*
It reads nicely at a call site, but it puts a filesystem-naming policy inside a module whose job
is to describe a reference — the exact coupling the owner drew the line against (2026-07-17).
The three lines above are the whole of it, and they belong to the caller that knows its own file
layout.

`ReferenceInfoError` is a `#[non_exhaustive]` `thiserror` enum, the house pattern: the `.fai`
read, the FASTA read, a malformed FASTA (T3), a duplicate contig name (T2), the
fasta-vs-fai disagreement (§3.3), a FASTQ index where a FASTA one was expected (§3.8), a
compressed reference (§3.8, deferred), and a `.fai`-write failure (from `read_reference_verifying_or_creating_fai`'s write
step, §3.11), each carrying the path and enough detail to act on. Nothing
in a supplied reference may panic — its bytes are attacker-influenced input, which is the reason
`RawChromReader` validates its own `.fai` fields rather than trusting them
([`raw_chrom_reader.rs:38-41`](../../../../src/ng/raw_chrom_reader.rs)).

---

## 3. Decisions

### 3.1 Two sources, and the source names the cost

**Owner's call, 2026-07-17: the module builds a contig list from a FASTA (optionally checked
against a `.fai`) or from a `.fai` alone — nothing else.** The two are not redundant; they are
the two things a caller can be holding, and they answer a different amount at a different cost:

| source | reads | gives | costs |
|---|---|---|---|
| `Fai(path)` | the index | names, order, lengths | one small file |
| `Fasta { fai: None }` | the whole reference | + every contig MD5 + the reference digest | a genome read |
| `Fasta { fai: Some }` | the whole reference + the index | the same, **and** a proof the two agree (§3.3) | a genome read |

The last two are one variant, because they are **one code path**: reading the FASTA. An index,
when supplied, only adds a comparison at the end (§4) — it changes nothing the reader does up to
that point, which is why it rides as an `Option` on the `Fasta` variant rather than a variant of
its own. And the FASTA read is not a redundant twin of the `.fai` read: it is the only source of
the MD5s (§3.4), and the only one that can answer at all when no `.fai` exists.

### 3.2 Decided: a supplied `.fai` is always checked — no trust-the-index fast path

**Owner's call, 2026-07-17.** When a caller passes a `.fai` alongside a FASTA
(`Fasta { fai: Some }`), the FASTA is read and the two are fully verified (§3.3). There is no
"trust the index, skip the read" mode and no `verify: bool`.

*Rejected: trust the `.fai` by default, verify on request.* It is cheaper, and it is what the
tree does today (`validate_fasta_agreement` never opens the FASTA — it compares two
descriptions of it and trusts both). But an option that defaults to off is an option nobody
turns on, and the failure it guards is silent wrong bases, not a crash. A check that runs is
worth more than a check that exists. So the check is not optional *given the inputs*: passing a
`.fai` is what asks for it, and passing one always means it runs.

**The cost is legible in the source the caller chose (§3.1).** A caller that cannot afford a
genome read passes `Fai` and gets names/order/lengths with `md5: None` — which is harmless
downstream, because an absent MD5 is a wildcard
([`fasta/mod.rs:43-55`](../../../../src/fasta/mod.rs)). A caller that wants the MD5s and the
verification passes `Fasta`. Neither is a hidden default: the expense is written at the call
site, in the variant name.

### 3.3 Decided: "they match" means the whole index, not just names and lengths

The comparison reconstructs the **full `.fai` record** for each contig from the FASTA —
name, length, offset, `line_bases`, `line_width` (`noodles_fasta::fai::Record::new` takes
exactly those five; they are what the format carries) — and compares it against the on-disk
index field for field.
`fai::Record` derives `PartialEq`, so the comparison is `==` and the work is in reporting
*which* field of *which* contig disagreed.

#### What the standard tools check: nothing

Asked what criterion the standard tools use for "the FASTA and the `.fai` match", the answer
from the source is that **the question does not exist in htslib. The `.fai` is trusted
absolutely.** Three paths, all of vendored `htslib/faidx.c`:

- **Loading** (`fai_load3_core`, `faidx.c:587-710`) — if the `.fai` is there, read it and open
  the FASTA. No comparison, no mtime check, nothing. The only branch that looks at the FASTA's
  content is the one that *builds* an index because none exists.
- **Parsing** (`fai_read`, `faidx.c:399-460`) — checks that four numbers `sscanf` per line.
- **Using** (`fai_retrieve`, `faidx.c:736-818`) — seeks by pure arithmetic on the index's own
  numbers (`offset + beg/line_blen*line_len + beg%line_blen`) and strips line terminators by
  pointer arithmetic on `line_blen`/`line_len` rather than by looking for `\n`. Its only
  guards are `line_blen <= 0` (`faidx.c:748` — the divide-by-zero `RawChromReader` copied) and
  a short read at EOF.

So on a re-wrapped FASTA with a stale index, htslib hands back a "sequence" with newlines
embedded in it, or the wrong bases, and reports success. **Owner's call, 2026-07-17: we do not
follow this.** htslib is a library whose contract is that the index is an input — garbage in,
garbage out — and it has nowhere to put such a check. We are the layer that *hands* it the
index, which is exactly the place the check belongs, and §3.2 already decided to pay for it.

But the source answers the real question better than a convention would. `fai_retrieve` shows
precisely which fields a fetch's correctness rests on: **`offset`, `length`, `line_bases` and
`line_width`** — all four, with the geometry doing the load-bearing work and nothing but
`line_blen > 0` standing behind it. That is the criterion, derived from what breaks rather
than from what anyone documents, and it is the list §3.3 checks.

#### Names and lengths alone are not enough, and the counterexample is a one-liner

Run
`fold -w 60 ref.fa > rewrapped.fa` and keep the old index: every name matches, every length
matches, and every offset and line width is wrong. On a single-contig reference even the
offset matches — the header line is unchanged — and only `line_bases` gives it away. So a
names-and-lengths check passes a reference whose every fetch will return the wrong bases,
because **the fetchers do arithmetic on the geometry, not on the lengths**
([`raw_chrom_reader.rs`](../../../../src/ng/raw_chrom_reader.rs) turns
`(offset, line_bases, line_width)` into a file position for every read). Checking what the
consumers actually depend on is the only version of this check worth its genome read.

*Rejected: `ContigList::first_disagreement`* ([`fasta/mod.rs:69-100`](../../../../src/fasta/mod.rs)).
It is `pub(crate)`, so ng may call it, and the brief expected it here — but it compares
name/length/MD5 and **has no notion of offset or line geometry**, which is precisely the half
that matters (above). It cannot be extended, because `src/fasta/` is frozen. It stays the
right tool for the job it was written for — comparing two *tables* — which is the `@SQ`
reconciliation deferred in §7, not this. This module writes its own record-wise comparison
over `fai::Record`s.

*Rejected: check the offsets by seeking, without a genome read.* Seek to each declared offset,
confirm a header sits where it should. Cheap, and catches a truncated or shuffled index — but
it cannot see a wrong length or wrong bases, and it misses the re-wrap case entirely. It buys
a fraction of the check for a fraction of the cost, and §3.2 already decided the cost is paid.

*Rejected: verify by hashing through the `.fai` (production's `compute_contig_md5_streaming`
shape — seek to `record.offset()`, hash `length` bases,
[`common.rs:250-299`](../../../../src/pop_var_caller/common.rs)).* **This is circular and the
trap is worth naming**: it reads the bases *the index points at*, so it can only ever confirm
that the index agrees with itself. Every fai-driven reader in the tree — `RawChromReader`,
`StreamingChromRefFetcher`, `ManualEvictChromRefFetcher` — has this property. **A `.fai`
cannot be verified by anything that uses it to read.** The pass must be sequential, from byte
zero, parsing the FASTA on its own terms. That is why §4 is not a call to an existing fetcher.

### 3.4 Decided: the digest is the SAM spec's `@SQ M5`, computed the way samtools computes it

**Owner's call, 2026-07-17: the SAM spec defines a hash of a contig's sequence, and that is
the hash we match.** Both digests are MD5 over **bases** — not over the file's bytes. They are
computed **because the bases are already
going past** (§3.2 decided the genome read; the digest is one extra `update` per window on
bytes we are holding anyway), and they are returned rather than checked here because their
consumers are elsewhere:

- **The per-contig MD5 is what an alignment file's `@SQ M5` gets compared against** — the
  reconciliation deferred in §7, and the reason it will have teeth at all (T6). That check is
  the concrete answer to "who needs this": it is the difference between *this BAM names the
  same contigs* and *this BAM was aligned to these bases*.
- **The reference digest** identifies the assembly as a whole — comparable to a `.cat`
  header's `reference_md5` (below) — and, being a concatenation in file order, is the one
  value that would catch a reordering (T1).

This module computes and hands them over. It does not decide what disagreement means; a
caller that has something to compare against does.

- **Per contig** → `ContigEntry::md5`. It has no choice: that field means `@SQ M5`
  ([`alignment_input.rs:331-354`](../../../../src/bam/alignment_input.rs) decodes it from the
  header), and M5 is defined on uppercased bases. A file-byte digest there would be a
  different number wearing the same name, and `ContigEntry`'s `PartialEq` would compare it
  against a real M5 and report a mismatch that isn't one.
- **Whole reference** → `ReferenceInfo::reference_md5`: the uppercased bases of every
  contig, concatenated in file order, which is exactly `ssr-catalog`'s `reference_md5`
  ([`catalog/mod.rs:196-211`](../../../../src/ssr/catalog/mod.rs)). So an ng table and a `.cat`
  header can be checked against each other, and re-wrapping the FASTA does not change it.
  Being a concatenation in order, it also pins the **order** — two references with the same
  contigs shuffled have different digests.

*Rejected: MD5 of the file's raw bytes.* It answers "is this the same file", which nothing
asks; two identical assemblies wrapped differently would disagree, and it is comparable to no
M5 anywhere in this repo or outside it. *(Rejected: both, named separately — free to compute,
but a second digest with no consumer is a field someone will compare by mistake.)*

#### The rule, and where it comes from

The SAM spec's `@SQ M5` is *"MD5 checksum of the sequence in the uppercase, excluding spaces
but including pads (as `*`)"*. The prose is not the thing to implement, though — the thing to
implement is **what writes the M5 values we will be compared against**, and that is
`samtools dict` (`samtools/dict.c:78-86`):

```c
for (i = k = 0; i < seq->seq.l; ++i)
    if (seq->seq.s[i] >= '!' && seq->seq.s[i] <= '~')
        seq->seq.s[k++] = toupper_c(seq->seq.s[i]);
hts_md5_update(md5, (unsigned char*)seq->seq.s, k);
```

**Keep the bytes in `[0x21, 0x7E]` — printable, non-space — uppercase them, hash them.** That
is the spec's sentence made exact: every whitespace byte (space, tab, `\n`, `\r`) is out
because it is below `0x21`; a pad `*` (`0x2A`) and every IUPAC code are in, because they are
not.

*(The other producer of real M5s is Picard's `CreateSequenceDictionary`, and BAMs in this lab
come from both. Verified from source, 2026-07-18: they **agree on every well-formed FASTA and
diverge on exactly one pathology.** Picard upcases then hashes whatever its FASTA reader
returned (`picard/util/SequenceDictionaryUtils.java`, `makeSequenceRecord` + `md5Hash`), and
that reader strips only line terminators and trailing per-line whitespace, not spaces or tabs
in the *middle* of a sequence line (`htsjdk/.../FastLineReader.java` `atEoln` stops only at
`\n`/`\r`). samtools' `>= '!' && <= '~'` filter drops those bytes; Picard keeps and hashes
them. So on a FASTA with an embedded mid-line space or tab the two producers disagree, and the
M5 is genuinely undefined there — and note it is **Picard** that departs from the spec's
"excluding spaces", not samtools. We follow samtools, which is the spec read literally; no real
FASTA carries mid-line whitespace, so the divergence is theoretical. Recorded because if it
ever bites, it bites at the `@SQ` reconciliation (§7) — a BAM's M5 vs ours — not here, and the
reader on the far side should know the number can be ill-defined rather than assume a bug in
one side.)*

**One rule, three fields.** htslib defines the same predicate in `faidx.c:48-50` —
`isgraph_(c) = c > ' ' && c <= '~'` — and `fai_build_core` counts the index's numbers with it
(`faidx.c:268-275`): `line_bases` and `length` count only `isgraph` bytes, while `line_width`
counts raw bytes including the terminator. So the *same* predicate that defines the M5 also
defines the `.fai` lengths our §3.3 comparison reconstructs. Implement it once. Two things
then fall out rather than needing their own handling: a `\r\n` file is correct for free (`\r`
is not a base, and is in `line_width`), and so is any FASTA with trailing spaces.

*Rejected: production's rule* — skip `\n` and `\r`, uppercase everything else
([`common.rs:285-290`](../../../../src/pop_var_caller/common.rs)). It hashes spaces and tabs,
which the spec and both producers exclude, so on a FASTA carrying them our M5 would disagree
with the `@SQ M5` of every BAM aligned to it — a spurious mismatch on somebody's real data
(T6), and one we would have to explain rather than fix, since `src/` is frozen. **This is a
deliberate divergence from production, and the first place ng's copy differs from the original
by being *more* correct rather than differently shaped.** Production is not wrong for its own
job: `verify_fasta_matches_psp_chromosomes` compares its number against a `.psp`'s, and both
sides are its own. It is wrong as a rule for M5.

**Uppercasing is not `canonicalise`** ([`fetcher.rs`](../../../../src/fasta/fetcher.rs)):
canonicalisation folds IUPAC codes to `N`, which would change the digest. `R` hashes as `R`.

### 3.5 Decided: home is a new `src/ng/reference_info.rs`

Named for its primary output, `ReferenceInfo` — the `ref_seq.rs` → `RefSeq` convention, a module
named for the type it exists to produce. The earlier name `contig_list.rs` was retired for the
same reason `ReferenceContigs` was (§3.9 note): the module produces the contig table *and* the
digest *and* the reconstructed index, so a contig-only name undersold it. Not `contig_table.rs`
either, and that near-miss is still worth stating: `ContigTable` is already a trait in
`ref_seq.rs` ([`ref_seq.rs:212-216`](../../../../src/ng/ref_seq.rs)), and it is a **consumer
capability** — "this reference can tell you how long its contigs are", a bound the typed-region
walk demands. This module is the **producer**; giving the file that name would put two different
things under one word.

*Rejected: fold it into `ref_seq.rs`, beside the trait and the impls that consume its output.*
That is the strongest argument for the other side — `ResidentRefSeq::new` and
`WindowedRefSeq::new` both take a `ContigList` they cannot build
([`ref_seq.rs:367`, `:524`](../../../../src/ng/ref_seq.rs)), so the caller currently has to
hand-roll one, which is why every construction site in the tree is a test fixture
([`anchor.rs:80-90`](../../../../src/ng/region_typing/anchor.rs) builds `ContigEntry`s by
hand). But `ref_seq.rs` is 1280 lines and its own spec says it splits into `ref_seq/` when it
grows ([`ref_seq.md`](ref_seq.md) §5); and this module's dependencies point nowhere near it
(§3.6). A file that shares no types with its neighbours does not belong in their file.

### 3.6 Decided: a leaf — it depends on `crate::fasta` and noodles, and nothing else

**Owner's requirement, 2026-07-17: a simple interface, no other assumptions about ng.** Its
output types — `ReferenceInfo`, `ContigInfo` — are **ng's own plain data**, and that is
deliberate: **ng codes what is clearest for ng and does not carry production types for
compatibility's sake** (owner, 2026-07-18 — port-back is a later concern, not a design
constraint now). So `ContigInfo` is the everything-about-a-contig type production never had, not
a wrapper around `ContigEntry`. The module imports `crate::fasta::{ContigEntry, ContigList}` only
because `contig_list()` builds them for ng's not-yet-migrated consumers (§2, a transitional
bridge) and `first_disagreement` is reused later (§7). It does **not** import `ContigId`,
`RefSeq`, `ContigTable`, or any ng *step* type — a consumer needs no ng-step vocabulary to use it.

Read-only reuse across the fence is fine and costs production nothing
([`ng/mod.rs:11-15`](../../../../src/ng/mod.rs)). Two specific reaches, both `pub(crate)` and
both read-only: `crate::fasta::{ContigEntry, ContigList}`, and — only if a caller wants hex —
`pop_var_caller::common::format_md5_hex` ([`common.rs:60`](../../../../src/pop_var_caller/common.rs),
reachable in-crate, as `ssr/catalog/mod.rs:211` already demonstrates). The B2 guard
(`rg 'use crate::ng' src/ --glob '!src/ng/**'`) is unaffected: the dependency points from ng
into production, which is the safe direction.

**`ContigId` is absent on purpose**, though this module is what defines it: `ContigId(i)` is
the i-th entry of the list this returns, in file order (T1). Naming the type would buy nothing
— a `Vec` is already indexed — and would tie a leaf to ng's vocabulary.

**The same line runs through the interface (§2), and it is the owner's, 2026-07-17.** The module
does one thing — describe a reference's contigs — and the test for whether a behaviour belongs
here is whether it is about *the reference's contigs* or about *how a caller found its files*.
**Checking a FASTA against a `.fai` is the former: it lives here.** **Discovering which `.fai` to
use — the sibling naming convention, a CLI's `--reference` flag, anything about file layout — is
the latter: it does not.** So the *pure* layer takes explicit paths and never probes the
filesystem for a file it was not handed; a caller wanting sibling behaviour composes the pure
`sibling_fai_path` with its own existence check (§2). This is why there is no
`from_fasta_checking_sibling_fai` constructor and no reasoning about any particular consumer in
the pure layer: the moment the *reader* knows what a `--reference` flag is, the leaf has a second
root.

**The one sanctioned exception is `read_reference_verifying_or_creating_fai` (§3.11)** — a *named* convenience that
explicitly does adopt the sibling convention (and can write a `.fai`). It does not blur the line;
it stands on the far side of it on purpose, so the pure primitives stay pure and the caller picks
the layer. A single named orchestrator that bundles the file-layout policy is not the same as the
reader silently absorbing it.

### 3.7 Decided: a caller-held, thread-safe, compute-once cache

**Owner's call, 2026-07-17.** The whole-genome MD5 pass (§3.2) is the expensive arm — seconds
to tens of seconds on a real reference — and the consumer that makes it matter is **parallel
per-sample readers**: a cohort validates each sample's alignment file against the *same*
reference, so N worker threads would otherwise each read the whole genome to build the same
table. The cache turns that into one read the others wait on. This mirrors production's cohort
parallelism, so it is a real shape, not a hypothetical one.

**A caller-held object, not a global static.** `ReferenceInfoCache` is a value the caller
constructs (normally once per run) and shares by reference. Rejected: a process-wide `static`
that `read_reference_info` transparently memoises. It is more convenient at the call site
and it is what "call the function again, get the same answer for free" would suggest — but it
is hidden global mutable state in a module sold as a simple pure function (§3.6), it leaks
between tests unless made clearable, and a stale-hit (T8) becomes a process-wide fact instead
of a bounded one. The pure `read_reference_info` stays; the cache is a thin wrapper over
it, so the single-threaded caller pays for none of this.

**The interface** is §2's `get_or_read`. Two properties fix its shape:

- **The result is `Arc`-shared.** `ReferenceInfo` is immutable once built, so a hit hands
  back an `Arc` clone (a pointer bump), never a copy of the table. `Arc<ReferenceInfo>` is
  `Send + Sync`, so every waiting thread gets the same shared value.
- **`&self`, not `&mut self`.** The cache mutates through interior locks, so `&cache` is what
  the workers hold — a single `ReferenceInfoCache` shared across threads, no per-thread
  copy, no outer `Mutex` the caller has to manage.

**The key** is `(which source, the file stats)`, and both halves are load-bearing:

```
Fai(p)                         -> key = Fai( stat(p) )
Fasta { fasta, fai: None }     -> key = FastaAlone( stat(fasta) )
Fasta { fasta, fai: Some(f) }  -> key = FastaWithFai( stat(fasta), stat(f) )

stat(p) = ( path as given, file length, mtime )
```

- **The source discriminant is in the key** because the three shapes return *different* results
  for the same bytes — `Fai` yields `md5: None`; `Fasta { None }` yields the table with MD5s but
  *unverified*; `Fasta { Some }` yields it *verified against that index* (§3.1). A stats-only key
  would let a cheap `Fai` call poison the entry a later `Fasta` call needs (handing back an
  MD5-less table), or let an unverified `Fasta { None }` satisfy a `Fasta { Some(stale) }` that
  ought to error. The `fai` stat is in the key of the last shape for the same reason: two
  different indices checked against one FASTA are two different questions. (Owner confirmed,
  2026-07-17.)
- **The path is in the key** so two different references never collide; **size and mtime** so
  the *same* path, re-stat'd after the file changed, is a new key and misses rather than
  returning stale bytes — within the limits T8 records.
- **Paths are used as given, not canonicalised.** Two spellings of one file (`ref.fa` vs
  `./ref.fa`, a symlink) are two entries and compute twice — harmless (same answer), and
  cheaper than the extra `realpath` syscall and its failure mode. If that redundancy ever
  matters, canonicalise here; it changes nothing observable but the hit rate.

**The single-flight algorithm**, explicitly — this is the part with the subtle correctness:

1. **Stat the file(s)** to build the key. A stat failure is a `ReferenceInfoError`, returned and
   **not** cached (there is no key to cache it under).
2. **Lock the map**, a `Mutex<HashMap<Key, Arc<Mutex<Option<Arc<ReferenceInfo>>>>>>`.
   Get-or-insert the key's **slot** (`Arc<Mutex<Option<…>>>`). **Clone the slot `Arc` and
   unlock the map at once** — the map lock is held only for this lookup, never across the read,
   so a thread working on one reference never blocks a thread asking for another.
3. **Lock the slot.** If it holds `Some(arc)`, clone the `Arc` and return — a hit.
4. If it holds `None`, call the pure `read_reference_info` **with the slot lock still
   held**. On `Ok`, store `Some(Arc::new(result))` in the slot and return the clone. On `Err`,
   **leave the slot `None`** and return the error.

Holding the slot lock across the read is exactly the single-flight contract: a second thread
on the *same* key blocks at step 3 until the first finishes, then falls through to a hit — it
waits, it does not start a second genome read. A thread on a *different* key holds a different
slot lock and runs fully in parallel.

**Errors are not cached** (step 4). A malformed FASTA would re-fail identically, so caching it
looks safe — but a *transient* I/O failure (an NFS stall, a disk timeout) would then be cached
permanently for that `(source, stat)`, turning a blip into a dead reference for the rest of the
run. Leaving the slot `None` costs a retry on the rare error and is correct on the common one.
This is also why the slot is a hand-rolled `Mutex<Option<…>>` and not
`OnceLock::get_or_init`: the latter has no stable fallible form and would cache the failure.

**No eviction, and that is not a leak.** Entries live for the cache's lifetime. The number of
*distinct* references a run touches is ~1 (a cohort shares one reference), so the map holds a
handful of `Arc`s at most; an LRU here would be machinery guarding against a size that does not
occur. If a future consumer ever streams thousands of references through one long-lived cache,
that is the point to add bounds — with a `log` line when it drops one, per the house "no silent
caps" rule. Recorded so the absence is a decision.

**When the cache is not usable, bypass it — never fail because of it.** If a platform or
filesystem cannot report mtime (`Metadata::modified()` errors), `get_or_read` computes directly
via `read_reference_info` and returns the answer uncached, rather than erroring. The cache
is an optimisation; it must never be the reason a correct read does not happen.

### 3.8 The `.fai` format we accept — and there *is* a spec

The `.fai` is not folklore: it has a written spec, **`htslib/faidx.5`** (`man 5 faidx`, shipped
in htslib and the samtools distribution), backed by the reference implementation this doc keeps
citing. It is not an hts-specs standard the way SAM/VCF/CRAM are, but it is a documented format,
and knowing its exact contents is what lets §3.3 compare field-for-field and §4 reproduce it
byte-for-byte. A `.fai` line is TAB-delimited: **five columns for FASTA, six for FASTQ.**

| # | column | meaning |
|---|---|---|
| 1 | `NAME` | first word of the `>` header (to first whitespace); the same string as `@SQ SN` |
| 2 | `LENGTH` | total bases; the same number as `@SQ LN` |
| 3 | `OFFSET` | byte offset (from 0) of the sequence's first base |
| 4 | `LINEBASES` | bases per line |
| 5 | `LINEWIDTH` | **bytes** per line, *including* the newline |
| 6 | `QUALOFFSET` | *FASTQ only* — byte offset of the first quality character |

These five are exactly `noodles_fasta::fai::Record`'s fields, exactly what §3.3 reconstructs and
compares, and exactly what §4 emits. Two constraints the spec states, both of which §4 already
leans on: **within a sequence every line has the same width except a possibly-shorter last one**
(htslib's *"Different line length in sequence"*, `faidx.c:291`), and **`LINEWIDTH − LINEBASES`
is the terminator width** — so CR-LF is handled by that difference, not by a special case. The
spec's own worked example is the proof and doubles as a test vector (§6): one FASTA, LF then
CR-LF, `LINEBASES` fixed at 30 while `LINEWIDTH` goes 31→32 and every `OFFSET` shifts by the
extra CR. §4's predicate reproduces both columns, which is why our reconstruction equals
`samtools`' `.fai` on either line ending.

**Two things the format admits that we do not accept, each a decision, not an accident:**

- **A FASTQ index (six columns).** A reference is FASTA; a `.fq.fai` is out of scope. Detect the
  sixth column and **reject** with a clear error rather than misparse — htslib itself
  distinguishes them (*"Possibly this is a FASTA index, try using faidx"*, `faidx.c:439`), so a
  `QUALOFFSET`-bearing index is a recognisable, nameable wrong input, not a corruption.
- **A bgzip-compressed reference** (`ref.fa.gz`, whose `.fai` offsets are into the *uncompressed*
  stream and which needs a companion `.gzi` to seek). **Deferred** (owner, 2026-07-18) — see §7.
  Today the module reads a plain FASTA: §4 parses raw file bytes from offset zero, and a `.gz`
  would parse as compressed garbage, so a compressed reference is an error now, not a silent
  misread. That matches ng's existing reference stack, which is uncompressed-only
  (`RawChromReader` copies production's plain-`File` offset arithmetic, `raw_chrom_reader.rs`).

### 3.9 Decided: `ContigInfo` carries everything, and ships a `.fai` writer

**Owner's call, 2026-07-18.** The FASTA pass (§4) computes every `faidx.5` column — `offset`,
`line_bases`, `line_width`, alongside the MD5s (§3.8). Rather than keep that geometry in a second
list beside the contig table, **one type carries all of it**: `ContigInfo` (§2) holds name,
length, geometry and md5 together, `ReferenceInfo.contigs` is a `Vec<ContigInfo>`, and `write_fai`
turns those into a `.fai`. This closes what §7 deferred and reverses what §8 leaned — soundly,
because the reason §8 said no has been removed.

**Why §8's "no" no longer holds.** §8 objected that returning the index would put
`noodles_fasta::fai::Index` — a noodles type — in ng's public output, breaking the leaf rule
(§3.6). The objection was about the *type*, not the *data*. So the geometry rides on **`ContigInfo`,
ng's own plain struct**, and the `.fai` is built from it. No noodles in the signature; the
objection is spent, not violated.

**One list, not a table plus its shadow.** An earlier draft carried both a `ContigList` and a
parallel `Vec<FaiRecord>`, with name+length duplicated across them and a projection to keep in
step. `ContigInfo` retires that: name, length, geometry and md5 are fields of one struct, so there
is nothing to keep consistent — the redundancy is gone at the root, not managed. Where ng's
own not-yet-migrated consumers still want a `ContigList` (the `RefSeq`/`ContigTable` seam,
`GenomeRegions` validation), `ReferenceInfo::contig_list()` projects one on demand (§2, §3.6) — a
transitional bridge, removed when they take `ContigInfo`.

**The writer's real use is the FASTA-without-index case.** `contigs` is populated from every
source — parsed from the `.fai` (`Fai`, `Fasta { Some }`) or computed (`Fasta { None }`) — but
writing a `.fai` for a reference that *already has one* is pointless. The case that matters is a
FASTA whose index does not exist: `read_reference_info` reads it once, and
`write_fai(&result.contigs, &sibling_fai_path(&fasta))` drops the index beside it. That is ng
indexing a reference itself, which the reference convention today makes a `samtools faidx`
prerequisite (`cli.rs:123-125`) — a prerequisite this makes optional.

**Two correctness pins.** `write_fai`'s output is **byte-identical to `samtools faidx`** — same
five TAB columns, same `\n`, no trailing blank line — so a `.fai` we write is indistinguishable
from one samtools wrote, and §6 tests exactly that against the oracle. And it writes
**atomically** (`.tmp` + rename), because a half-written `.fai` is silently valid — it has no
header, no checksum, nothing that fails a reader — so every future run would trust a truncated
index (the same hazard `typed_regions_cli.md` §6 records for its own output).

### 3.10 Decided: an opt-in "verify in the background" entry point

**Owner's call, 2026-07-18.** The genome read is the slow part and it is usually the *first*
thing a run does, so the rest of the program waits on it. This entry point hides that latency:
read the `.fai` now (fast — names, lengths, order, geometry), hand the caller that info
immediately, and run the FASTA verification (§3.3 geometry + §3.4 MD5s) on a background thread.
If everything matches, nothing happens; if the `.fai` is stale or the FASTA unreadable, the
error surfaces when the caller **joins** the handle. Opt-in: the default stays the blocking
`read_reference_info`, because the discipline below is a real cost a caller takes on knowingly.

**It reuses the two arms *through the cache*, it is not new verification.**
`read_fai_verify_in_background` is a thin wrapper: the foreground is
`cache.get_or_read(Fai(fai))`; the background thread is
`cache.get_or_read(Fasta { fasta, fai: Some(fai) })`. Both go through the cache on purpose (next
paragraph). No new check, no second implementation of anything — only the threading.

**crossbeam, not async; and rayon has a place, but not this one.** async would import a runtime
for a single background task in a tree that has none; a dedicated thread delivering its result
over a **crossbeam** channel is the right size (owner's instinct, confirmed). The task does **not**
go on a rayon worker: a seconds-long blocking genome read would occupy a pool thread meant for
short, stealable, CPU-bound work and starve it. rayon's place is *inside* the verification — the
deferred two-phase parallel hash (§7): reconstruct geometry serially from byte zero (phase 1, and
the §3.3 check), then hash the contigs in parallel by **phase 1's own** offsets (not the `.fai`'s,
so not circular). So "crossbeam and rayon" maps exactly: crossbeam owns the one background thread
and its result; rayon, when phase 2 lands, fans the hashing within it.

**The whole value is the error surfacing, so the handle refuses to lose it.** Proceeding on the
`.fai`'s unverified info is *optimistic*: if verification fails, the work done meanwhile was on
possibly-wrong data. That is tolerable only if the caller aborts on the late error **before
committing output** — and a background handle makes forgetting easy (drop it → detached thread →
the failure vanishes → a silent wrong-reference run). So `VerificationHandle` is `#[must_use]`,
`join` is mandatory before output-commit, and `Drop` without `join` **logs** rather than swallows.
A caller wanting the discipline enforced structurally can run inside a `crossbeam::thread::scope`,
where the verification thread must join at scope end and its `Result` is checked there — forgetting
becomes a compile-shaped impossibility rather than a convention.

**The md5s arrive by *return*, not by in-place mutation — and that is deliberate.** The tempting
shape is to hand the caller `ContigInfo`s with `md5: None` and have the background thread fill each
`md5` in place once the scan finishes. It reads as "even better for the caller" — one object that
completes itself — but it fights Rust and loses. To mutate data the caller is holding, the data
must be shared-mutable across threads: `md5` stops being a plain `Option<[u8; 16]>` and becomes a
`Mutex`/atomic, every read has to lock, and — worse — the caller cannot tell *when* a field is
filled without joining anyway, so it reads `None`-or-`Some` nondeterministically until it does. So
`ContigInfo` stays a plain immutable struct and the background thread **returns** the verified
`ReferenceInfo` (md5s filled) through `join`. The "upgrade" is then explicit and deterministic:
before `join`, `md5: None`; after, you use the returned value with `md5: Some`. No locks, no race,
and `ContigInfo`'s plain-struct cleanliness (the point of §2's merge) survives. Immutable data plus
message-passing is the idiom; the returned value *is* the message.

**Two facts that make the optimism safer than it looks.** A stale `.fai` from a re-wrap keeps
correct *names and lengths* — only the geometry breaks (§3.3) — so the contig table the caller
proceeds on is usually right even when the `.fai` is stale; what a stale `.fai` corrupts is the
*geometry the readers seek by* (§1's wider mandate), so this mostly protects `RawChromReader`, and
a reader that wants to can simply `join` before its first fetch. And `join` is not merely a verdict:
on success it returns the **verified** `ReferenceInfo`, MD5s included, so the caller's knowledge
legitimately *upgrades* from unverified-without-digests to verified-with — the digests the `Fai`
read could not give (T6) arrive with the all-clear.

**The cache is not an afterthought — it is why this is worth doing.** A background verify whose
result nobody could reach would have paid for a genome read and thrown it away. So both reads go
**through** the cache (§2), and the cache's single-flight machinery (§3.7) makes that do three
things at once, not just "insert":

- **The verified result is cached** under the `Fasta { Some }` key. Any later `get_or_read` of it
  — from a worker that was not handed the `join` result — is a **hit**, not a second genome read.
- **A concurrent reader of the same key coordinates, it does not duplicate.** Because the
  background thread holds that key's single-flight slot for the whole read (§3.7), a worker that
  calls `get_or_read(Fasta { Some })` *while the verify is still running* **waits on the slot** and
  gets the one result — it does not start its own read. The background verify *is* the single-flight
  owner of that key; `join` and a concurrent `get_or_read` observe the same `Arc`.
- **The `Fai` foreground caches too**, so a worker wanting only the cheap table hits as well.

On failure the slot is left empty (errors are not cached, §3.7), so a stale `.fai` does not poison
the key — the error reaches the caller through `join`, and any retry re-reads honestly. **The cache
must be `Arc`-shared** for this: the detached background thread holds an `Arc<ReferenceInfoCache>`
clone, which is why the entry point takes `&Arc<ReferenceInfoCache>` rather than `&self` (§2).

### 3.11 Decided: `read_reference_verifying_or_creating_fai`, the batteries-included entry point

**Owner's call, 2026-07-18.** The most convenient way to get a reference's info is to hand over a
FASTA and let the module do the right thing. `read_reference_verifying_or_creating_fai(cache, fasta)` (§2) is that: it
derives the sibling `<fasta>.fai` and branches —

- **`.fai` present** → `read_fai_verify_in_background` (§3.10): read the index now, verify the
  FASTA in the background, return immediately with `(info, Some(handle))`.
- **`.fai` absent** → scan the FASTA now (`get_or_read(Fasta { None })` — verified, MD5s), **write**
  the `.fai` beside it, and return `(info, None)`.

It is pure composition of §2's primitives; it adds no new reading or verification. What it adds is
two things the pure layer deliberately withholds, and each is a conscious exception:

**It is the one place the module adopts the `<fasta>.fai` convention.** §3.6 keeps sibling
discovery out of the core reader on purpose — one call must not behave differently by what sits
next to the FASTA. `read_reference_verifying_or_creating_fai` is the opposite by design: a *named* orchestrator whose whole
job is "do the conventional thing." The exception is contained — the pure primitives
(`read_reference_info`, `get_or_read`, `write_fai`, `sibling_fai_path`) still probe nothing — so a
caller that wants control drops to them, and one that wants convenience takes this. Both exist on
purpose.

**It is the one place a read can write a file, and a write failure is fatal.** The no-`.fai`
branch creates the index on disk. If that write fails — a read-only reference dir, a full disk —
`read_reference_verifying_or_creating_fai` returns `Err` (a `.fai`-write `ReferenceInfoError`, §2), even though the scan
produced the info. That is the owner's call: `read_reference_verifying_or_creating_fai` promises *fully set up or nothing*.
The escape hatch is explicit and documented — a caller that cannot or will not write the index
skips `read_reference_verifying_or_creating_fai` and calls `cache.get_or_read(Fasta { fai: None })` directly, which reads and
caches the info and writes nothing. So the fatal contract costs a read-only-dir caller one known
line, not a lost scan.

**The return is asymmetric, and honestly so.** `(info, Some(handle))` on the fai-present path
(unverified now, `join` later); `(info, None)` on the scan path (already verified, nothing pending).
`None` is not a missing handle — it is "there is nothing to wait for." A caller wanting uniformly
verified info writes `if let Some(h) = handle { h.join()?; }`; after that line the info is verified
on both paths. The alternative — always return a pre-completed handle — buys uniformity at the cost
of pretending there is work to await when there is not, so `Option` is the honest shape.

---

## 4. The streaming pass

**One pass, one buffer, never one contig.** The `Fasta` variant runs this loop whether or not a
`.fai` was supplied; when one was, its result is compared to the index afterwards (§3.3).

**The one predicate.** A byte is a **base** iff it is in `[0x21, 0x7E]` — printable,
non-space. That is htslib's `isgraph_` (`faidx.c:48-50`) and `samtools dict`'s M5 filter
(`samtools/dict.c:78-80`), and §3.4 explains why it is the only rule in play: it defines the
M5, `line_bases` and `length` all three. Everything else in the loop counts *bytes*.

```
open the FASTA
for each record:
    read the definition line          -> name (up to the first whitespace)
    note the byte offset now reached  -> the record's `offset`
    read the sequence, line by line:
        line_width = every byte of the line, terminator included  (\r\n counts as two)
        line_bases = the bytes of the line that are 0x21..=0x7E
        the first line fixes both; every later line must match, except the last,
            which may be shorter                                  (`faidx.c:277-292`)
        each base: to_ascii_uppercase, feed the contig MD5 and the reference MD5
        length += line_bases
    emit (name, length, offset, line_bases, line_width, contig md5)
```

Resident memory is one read buffer — production's `compute_contig_md5_streaming` uses 64 KiB
([`common.rs:35`, `:263-295`](../../../../src/pop_var_caller/common.rs)) and there is no reason
to differ. The two MD5 states are fed the same bytes, so the reference digest costs one extra
`update` per window, not a second pass.

Reconstructing `line_width`/`line_bases` this way is what makes the §3.3 comparison a
comparison of like with like: the numbers are built by htslib's own rule, so a disagreement is
the `.fai` being wrong about the file, never our arithmetic being wrong about the `.fai`.

**Why ng parses the lines itself, rather than calling a reader — including our own.** Four
readers could plausibly do this, and each is *almost* it. The one to address first is the one
§3.6 invites — *"we already have a memory-efficient streaming reader, why not reuse it?"*:

- **`crate::fasta`'s streaming fetchers** (`StreamingChromRefFetcher` /
  `ManualEvictChromRefFetcher`, and ng's own `RawChromReader`) — the reference stack's
  low-memory readers, built for exactly the "don't hold a whole contig" goal this pass shares.
  They are nonetheless the **wrong tool**, because they are **`.fai`-driven**: each seeks to the
  byte offset the `.fai` declares and does its line arithmetic from `line_bases`/`line_width`
  ([`raw_chrom_reader.rs`](../../../../src/ng/raw_chrom_reader.rs) turns
  `(offset, line_bases, line_width)` into a file position). So they cannot read a FASTA that has
  **no** `.fai` (the `Fasta { None }` case — nothing to drive them), and using them to verify a
  `.fai` is the **circular** check §3.3 forbids: they read the bytes the index points at, so they
  can only confirm the index agrees with itself. They read bases *at positions a trusted index
  gives*; this pass must do the opposite — read from byte zero to *discover* the geometry and
  hold the index to account. The memory win is not lost by declining them: the hand-rolled loop
  below is equally one-buffer streaming, and being from-byte-zero it *can* build and verify.
- **`noodles_fasta::io::Reader::records()`** — what `ssr-catalog` uses
  ([`catalog/mod.rs:195-210`](../../../../src/ssr/catalog/mod.rs)) and what
  `anchor.rs:112-120` uses. It materialises **each contig whole** (`rec.sequence()` is the
  entire sequence), which for a 250 Mb chromosome is exactly the allocation the owner's
  streaming requirement forbids, and exactly the shape [`ref_seq.md`](ref_seq.md)'s 14.6 GB
  lesson is about. The catalog can afford it because it hands whole contigs to `trf-mod`
  anyway. We cannot.
- **`Reader::read_definition` + `Reader::sequence_reader()`** — genuinely streaming: a
  `BufRead` over raw bases with the newlines already stripped, stopping at the next
  definition. It would give name + length + MD5 in bounded memory. But stripping the newlines
  destroys `line_bases` / `line_width` — so no geometry and, with it, **no line-consistency
  check**: a FASTA with a ragged interior line reads clean through it. §3.3 needs the geometry
  and §3.8 settled that we keep it, so this reader is out.
- **`noodles_fasta::io::Indexer::index_record`** / `fasta::fs::index` — reconstructs the index
  properly (it is what `samtools faidx` does, including the uniform-geometry errors this pass
  needs, so it *would* do the line-consistency check for us), streaming — but it consumes the
  bases without handing them over, so the MD5s would cost a **second** genome read.

The pattern is clear: every existing reader drops exactly one of the three things this pass needs
at once — **from-byte-zero parsing, the line geometry, and the bases** — because no other
consumer in either tree ever needed all three together. `sequence_reader` keeps bases, drops
geometry; `Indexer` keeps geometry, hides bases; the `.fai`-driven fetchers need an index to
exist before they will read at all.

So the loop above is Indexer's logic with a digest folded into it; it is ~40 lines
(`indexer.rs`'s `index_record`, plus a byte counter), and copying it is in character —
`RawChromReader` is the same move, a faithful copy of production's fetcher with one behaviour
changed, kept diffable against the original on purpose
([`raw_chrom_reader.rs:36-42`](../../../../src/ng/raw_chrom_reader.rs)). Its geometry
errors come with it: non-uniform `line_bases`, non-uniform `line_width`, an empty sequence.

**`Fai` reads the index and stops**: `noodles_fasta::fai::fs::read`, the same call
`validate_fasta_agreement` makes ([`alignment_input.rs:437-442`](../../../../src/bam/alignment_input.rs)),
then one `ContigEntry { name, length, md5: None }` per record in order. It never opens the FASTA,
so there is no line structure to check — the line-consistency guard above is a `Fasta`-arm
concern only; here the `.fai`'s own numeric validity (`line_bases > 0`, `line_width ≥ line_bases`)
is what a parse guards, the same fields htslib's `fai_insert_index` checks (`faidx.c:107-116`).

---

## 5. What will bite you

Every item here came from opening the file it cites.

**T1 — the order *is* the contract, and it is the whole reason this module exists.**
`ContigId(i)` means "the i-th entry" ([`types.rs:6-8`](../../../../src/ng/types.rs)), so the
`Vec` order this returns silently defines every `ContigId` downstream. Both readers preserve
file order (a `.fai` is written in FASTA order; the pass walks the file), so this is free —
but it is free *by construction*, not by check, and a future "sort the contigs for tidiness"
would silently remap every read in the caller. The `reference_md5` is order-dependent (§3.4),
which is the one thing that would catch it. Say so where the order is established.

**T2 — duplicate contig names resolve to one contig, silently.** `ContigList` does not check
for them, and `ResidentRefSeq` resolves a fetch by `ContigId → name → Repository::get`
([`ref_seq.rs:360-435`](../../../../src/ng/ref_seq.rs)) — so if a reference carries `chr1`
twice, two distinct `ContigId`s fetch the same bases and no error is ever raised. Nothing in
the tree rejects this today, because nothing in the tree builds a table from a reference. This
module is the first place that *can*, and it is the only place with the whole list in hand.
**Reject duplicates** — and know that this is **stricter than htslib**, deliberately.

`fai_insert_index` **warns and drops**: *"Ignoring duplicate sequence"*, keep the first, return
success (`htslib/faidx.c:123-127`). That is on both htslib paths —
building an index and reading one — so `samtools faidx dup.fa` succeeds, and **the `.fai` it
writes has fewer records than the FASTA has contigs**. Two consequences not to be surprised
by: a duplicate-carrying reference *does* have a `.fai` in the wild, and comparing our FASTA
pass against it (§3.3) would fire a count disagreement on a file samtools was happy with. We
reject before reaching that, and the reason we do not merely warn is `ContigId`: dropping the
second `chr1` silently **renumbers every contig after it** (T1). A warning is a defensible
answer for a library that resolves contigs by name; it is not one for a caller that resolves
them by position.

**T3 — the FASTA is user input, so every malformed shape is an error, not a panic.** The
pass's own arithmetic gives it the same exposure `RawChromReader` documents — a `line_bases`
of 0 divides by zero downstream ([`raw_chrom_reader.rs:38-41`](../../../../src/ng/raw_chrom_reader.rs)).
Bases before the first `>`; a definition with an empty name; a contig with no sequence lines;
non-uniform line widths in the middle of a contig (htslib errors here too — *"Different line
length in sequence"*, `faidx.c:291`). A `\r\n` file needs no special case: §4's one predicate
excludes `\r` from the bases and `line_width` counts it, which is what htslib does and
therefore what `samtools`' own `.fai` says. Each is a `ReferenceInfoError` variant with the contig name and the byte
offset, never an `assert!` — the sweep runs in `--release` (`ng_proposal.md` §10), where a
`debug_assert` would hand back a wrong answer instead of a message.

**T4 — a path handed in that does not open is an error, never a silent fall-back to something
more expensive.** The module reads exactly the paths it is given (§2), so this is short: if
`Fai(p)` names a file that isn't there, error — do **not** "help" by looking for a FASTA to
read instead, because `Fai` and `Fasta` cost three orders of magnitude apart (§3.1) and that
help would turn a typo into a silent ten-minute genome read. Likewise a `Fasta { fai: Some(f) }`
whose `f` is missing errors on open; it does not quietly degrade to `fai: None` and skip the
check the caller asked for. Error the way `validate_fasta_agreement` does
([`alignment_input.rs:432-436`](../../../../src/bam/alignment_input.rs)). The graceful
"check the index if it happens to be there, read alone if not" belongs to the caller, who
expresses it by choosing `fai: Some`/`None` — the module never guesses (§2).

**T5 — the `u64 → u32` narrowing is not this module's problem, and taking it on would be a
bug.** `ContigEntry::length` is `u64` ([`fasta/mod.rs:37-43`](../../../../src/fasta/mod.rs))
and `regions::ContigBounds::length` is `u32` ([`regions.rs:57-63`](../../../../src/regions.rs)),
and the brief flags the gap — but it is the **consumer's** gap. This module returns the table
`ContigList` describes, at `ContigList`'s width; a >4 Gb contig is a fact about the reference,
and refusing to *describe* one would refuse a reference for a caller that never builds a
`RegionSet`. The rule when a consumer does narrow is settled and unchanged — **fail rather
than fold** ([`typed_regions.md`](typed_regions.md) §4), `u32::try_from` and an error, the way
`narrow_for_fetcher` already does it ([`ref_seq.rs:103-130`](../../../../src/ng/ref_seq.rs)).
The one precedent that casts with `as` is a test fixture
([`anchor.rs:98`](../../../../src/ng/region_typing/anchor.rs)) and is not the model.

**T6 — an absent MD5 is a wildcard, so `md5: None` is safe and `md5: Some` is a claim.**
`ContigEntry`'s `PartialEq` treats `None` as matching anything
([`fasta/mod.rs:43-55`](../../../../src/fasta/mod.rs)). That is what makes `Fai`'s `md5: None`
harmless. It also means the `Fasta` variant hands a downstream comparison **real teeth for the
first time** — a table with MD5s compared against an `@SQ` list with M5s will now actually
compare them, and any deviation in §3.4's folding rule surfaces there as a spurious mismatch
on somebody else's data. That is the reason §3.4 copies production's byte handling rather
than improving it.

**T7 — the name is the first word, in both files, and they must agree.** A `.fai` record's
name stops at the first whitespace, so `>chr1 AC:CM000663.2` indexes as `chr1`. The pass must
take `noodles_fasta::record::Definition::name()` and not the whole line (`description()` is
the rest, and is separate), or every `Fasta { fai: Some }` comparison fails on every real
reference — the description field is where assemblies put their accessions, and both GRCh38 and
the tomato reference carry them.

**T8 — the cache key is `(path, size, mtime)`, and that is not a content hash.** mtime
resolution is 1–2 s on some filesystems, and mtime moves *backwards* on a restore-from-backup,
a `cp -p`, a `touch -d`. So a file rewritten within the resolution window, to the same size,
gives a **stale cache hit** — the old table returned for new bytes, silently (§3.7). The
strong key is the content MD5 we compute anyway, but
it cannot be a *pre-read* key — you would have to read the file to know it. **We accept the
hole, on one assumption we state rather than hide: a reference genome is a write-once,
read-only input, not edited mid-run.** Rewriting a reference under a running analysis, to an
identical size, within a two-second window, is not a workflow — it is operator error that no
stat key survives. So the rule for the cache is exactly: **do not reuse a
`ReferenceInfoCache` across files that change.** For references, they don't. (Owner
confirmed the trade, 2026-07-17.) A brief TOCTOU rides along — the file could change between
the stat-for-key and the read — and is the same class, accepted on the same ground.

---

## 6. Tests

The fixture already exists: `ref_seq.rs`'s `build_fasta` writes a multi-contig FASTA and a
matching `.fai` to a tempdir and returns the `ContigList`
([`ref_seq.rs:883-915`](../../../../src/ng/ref_seq.rs)) — but it writes **one line per
contig**, so it cannot exercise line geometry at all. This module needs a builder that wraps
at a chosen width; that is the first thing to write, and it is what T3 and §3.3's tests stand
on. Home: the module's own `#[cfg(test)]`, per the in-crate convention.

- **The sources agree on what they can both see.** `Fai(p)`, `Fasta { fai: None }` and
  `Fasta { fai: Some }` over one reference return the same names, order, lengths and geometry;
  only the MD5s differ (`None` from `Fai`, `Some` from the FASTA reads). This pins §3.1's claim
  that they are one description reached two ways.
- **`contig_list()` projects faithfully** (§2, §3.9). Assert the projected `ContigList` has one
  `ContigEntry` per `ContigInfo`, same name/length/md5 in the same order — the seam the production
  consumers rely on, and the proof the merge did not lose the compatibility the old stored
  `ContigList` gave.
- **The MD5s are the ecosystem's, and `samtools` is the oracle.** The strongest available
  check, and it needs no `samtools` at test time: run `samtools dict` **once** on the committed
  golden reference (`tests/data/tandem_repeat/synthetic_ref.fa`, which `anchor.rs:105-120`
  already reads), commit its `M5` values as golden constants, and assert we reproduce them.
  That pins §3.4 against the actual producer rather than against our reading of it — the whole
  point of the M5 being a *shared* number. Same for `LN`, which `samtools dict` prints from the
  same predicate (`dict.c:86`), and for the `.fai` `samtools faidx` writes for that file (§3.3
  compares against exactly those fields).
  Then the cheap companions: per-contig against `Md5::digest(filtered_uppercase_bases)`
  one-shot (production's own streaming-matches-one-shot shape,
  [`common.rs:452-489`](../../../../src/pop_var_caller/common.rs), including the soft-mask
  case), and the reference digest against the golden `.cat` header built from that same FASTA.
- **`write_fai` is byte-identical to `samtools faidx`** (§3.9). Read the golden FASTA, call
  `write_fai(&result.contigs, out)`, and assert the bytes equal the committed
  `samtools faidx` `.fai` — the same oracle as above, now over the writer. A round-trip rider:
  feed that written `.fai` back as `Fasta { fai: Some(written) }` and assert it verifies clean,
  which proves `read`→`write`→`read` is a fixpoint. And that `write_fai` is **atomic** — on a
  simulated mid-write failure, no partial `.fai` is left at the destination.
- **The predicate is the spec's, at its edges.** A sequence line carrying a space or a tab
  hashes as if it were not there, and `line_bases` does not count it — the case that
  distinguishes §3.4's rule from production's, and the reason we diverged. A `*` pad and an
  IUPAC `R` **do** hash, and are not folded to `N` (that would be `canonicalise`).
- **The re-wrap is caught** (§3.3's counterexample, and the reason the check is what it is).
  Write a FASTA, index it, re-wrap the FASTA at a different width, and assert
  `Fasta { fai: Some(stale) }` errors — naming `line_bases`, not just "mismatch".
  **Mutation-verify it**: weaken the check to names-and-lengths and this test must fail. It is
  the only test that distinguishes this spec from the cheaper one.
- **A single-contig re-wrap is caught too.** The offset is unchanged there (T-argument in
  §3.3), so this is the case an offsets-only check would pass.
- **Reordering is caught, and by the digest.** A FASTA whose contigs are permuted against its
  `.fai` errors; and its `reference_md5` differs from the original's (T1).
- **Duplicate names are rejected** (T2), from both `Fai` and `Fasta`.
- **A supplied index is optional, a missing supplied index is not** (§2, T4). `Fasta { fai: None }`
  reads cleanly with no index — assert it does not fail for want of one. `Fasta { fai: Some(f) }`
  with `f` absent errors on open — it does **not** silently degrade to reading the FASTA alone.
  And `sibling_fai_path` is a pure path function: assert it appends `.fai` and touches no
  filesystem (the module does no probing — that is the caller's).
- **Malformed input errors rather than panics** (T3): `line_bases: 0`, bases before the first
  `>`, an empty contig, a truncated last line, `\r\n`. The `\r\n` case must additionally
  assert the reconstructed offsets equal the `.fai` `samtools` would write.
- **The `faidx.5` worked example, verbatim, as a test vector** (§3.8). The spec's own two-line
  FASTA, indexed under LF (`one 66 5 30 31` / `two 28 98 14 15`) and CR-LF
  (`one 66 6 30 32` / `two 28 103 14 16`), pins §4's reconstruction against the format's authors,
  not just against local `samtools`. It is the cheapest possible proof that our five columns are
  *their* five columns.
- **A FASTQ index is rejected, not misparsed** (§3.8). Feed a six-column `.fq.fai` and assert a
  named error (a FASTA reader was handed a FASTQ index), never a silent wrong `LINEWIDTH`.
- **A compressed reference is a clean error** (§3.8, deferred). Feed a `ref.fa.gz` and assert the
  read errors rather than parsing compressed bytes as a FASTA.
- **Streaming, not resident.** The one property with no natural assertion — a test cannot
  watch RSS. What it *can* do is assert the pass never holds a contig: read a FASTA whose
  single contig is larger than any buffer the code allocates and assert it completes; or
  express it structurally, by keeping the pass's only owned byte buffer a fixed-size array
  (which `compute_contig_md5_streaming` does — `[u8; FASTA_MD5_BUFFER_SIZE]` on the stack,
  [`common.rs:263-264`](../../../../src/pop_var_caller/common.rs)) so a contig-sized
  allocation cannot compile. Prefer the second; note the first as what a reviewer should look
  for.

The cache (§3.7) is tested on its own, against the pure reader it wraps:

- **A hit computes once.** Wrap `read_reference_info` behind a counter (a test-only source
  that increments on each real read), call `get_or_read` twice on the same source, assert the
  underlying read ran **once** and the two `Arc`s point at the same allocation
  (`Arc::ptr_eq`).
- **Single-flight: the second thread waits, it does not recompute.** Spawn many threads on one
  key against a source whose read blocks on a barrier the test releases; assert the real read
  ran exactly once and every thread got the same `Arc`. This is the property the whole section
  exists for, and it is the one a single-threaded test cannot see.
- **Different keys run in parallel and do not serialise.** Two distinct references, two
  threads; neither blocks on the other's slot. Structurally: assert the map lock is not held
  across a read (a read that never completes must not wedge a `get_or_read` on a different
  key).
- **The source is in the key** (§3.7): `Fai(p)` then `Fasta(p)` on the same file are two
  entries — the second returns a table *with* MD5s, proving the first did not poison it (T6-adjacent).
- **An error is not cached** (§3.7 step 4): a source that fails once then succeeds (a
  test-only source flipping after its first call) must succeed on the retry — the failure left
  the slot empty. Mutation-verify: make the slot cache errors and this test must fail.
- **mtime-unavailable falls back, it does not fail.** Harder to force portably; at least assert
  the bypass path returns the same value as `read_reference_info` for a source whose stat
  the test stubs as unsupported.

The background-verify entry point (§3.10) is tested for both outcomes, since the failure path is
its entire reason to exist:

- **Clean reference: the info is available before the join, and the join upgrades it.** Assert the
  foreground `Arc<ReferenceInfo>` has `md5: None` and correct names/lengths immediately, and that
  `join` returns `Ok` with `md5: Some` — the digest arriving with the all-clear.
- **Stale `.fai`: the error surfaces on `join`, not before** (§3.10). Build a FASTA, index it,
  re-wrap the FASTA, and assert the foreground read *succeeds* (a re-wrap keeps names/lengths) but
  `join` returns the geometry `Err`. This is the test the whole feature exists for — the late error
  must actually arrive. Mutation-verify: make the background thread swallow its error and this test
  must fail.
- **`is_finished` never blocks, `join` always resolves.** Poll before completion → `false` without
  blocking; after → `true`; `join` then returns without hanging.
- **Abandoning the handle warns, does not swallow.** Drop a handle whose verification failed
  without joining, and assert the warning reaches stderr (not a silent pass). Pairs with the
  `#[must_use]` lint as the two guards against a lost error.
- **The background verify populates the cache** (§3.10). After a clean `join`, call
  `cache.get_or_read(Fasta { fai: Some })` for the same reference behind the read counter and
  assert it is a **hit** — the genome was read once, not twice.
- **A concurrent reader coordinates, it does not double-read** (§3.10). Start the background verify
  against a source whose read blocks on a barrier; on another thread call
  `cache.get_or_read(Fasta { fai: Some })` for the same key; assert it blocks until the barrier
  releases and then returns the *same* `Arc` the `join` yields — real read count exactly one. This
  is the single-flight-through-the-cache property, and the reason the entry point takes the cache.
- **A failed verify does not poison the key** (§3.10, §3.7). After a stale-`.fai` `join` returns
  `Err`, a later `get_or_read(Fasta { Some })` re-reads (and re-fails) rather than returning a
  cached error.

`read_reference_verifying_or_creating_fai` (§3.11) is tested on both branches, since it *is* the branch:

- **`.fai` present → background verify, no write.** Assert the returned info is immediate
  (`md5: None`), a handle is present, the sibling `.fai` on disk is **unchanged** (nothing
  written), and `join` returns verified info.
- **`.fai` absent → scan, write, no handle.** Assert the returned info is verified (`md5: Some`)
  with no handle, and that a `.fai` now exists beside the FASTA, **byte-identical to
  `samtools faidx`** (reuses §3.9's oracle — the write went through `write_fai`).
- **A `.fai`-write failure is fatal** (§3.11). Point `read_reference_verifying_or_creating_fai` at a FASTA whose sibling
  `.fai` cannot be written (read-only dir), no index present, and assert it returns `Err` — not a
  silent success that dropped the write. Then assert the escape hatch: `get_or_read(Fasta { None })`
  on the same reference succeeds and writes nothing.

- **`@SQ` ↔ `.fai` reconciliation, and with it the `ReadFilter::new` permutation hole** (§1)
  → **read ingestion's spec**, where the alignment file lives. It is `ContigList::first_disagreement`
  ([`fasta/mod.rs:69-100`](../../../../src/fasta/mod.rs) — `pub(crate)`, so ng may call it) over
  this module's output and the `@SQ` table, and the reason it is not here is that this module
  must not know what a BAM is (§3.6). Two notes for whoever picks it up: the check must
  compare **order**, not just membership, since a permutation is the bug; and the `Fai` arm's
  `md5: None` makes the M5 comparison a no-op (T6), so an M5 check needs the `Fasta` arm and
  therefore a genome read — a real decision, not a detail.
- **The coordinate-sort-order check** that stands on that reconciliation → same home.
  `extract_header` already rejects a non-coordinate-sorted alignment file
  ([`alignment_input.rs:297-303`](../../../../src/bam/alignment_input.rs)); what is missing is
  the reference-order half.
- **Consolidating ng's fai-driven readers onto `ContigInfo`** (§1's wider mandate) → a follow-up
  refactor of existing ng code. `RawChromReader` stops calling `open_contig` / `fai::fs::read`
  and takes a `&ContigInfo` (geometry) instead; `WindowedRefSeq::new`, which already receives a
  contig table from this module, receives the full `ReferenceInfo` (geometry included), so the
  `.fai` is parsed once per run rather than once per contig transition
  ([`raw_chrom_reader.rs:131`](../../../../src/ng/raw_chrom_reader.rs)). This refactor also drops
  `WindowedRefSeq`/`RawChromReader` off `ContigList` onto `ContigInfo`, retiring the transitional
  `contig_list()` bridge for them (§2, §3.6). ng-only (production's frozen readers keep their own
  parse). Deferred because it touches built ng code and is not needed to stand this module up — but
  it is the payoff that makes this module *the*
  `.fai` reader, not merely *a* contig-table builder.
- **A bgzip-compressed reference** (`ref.fa.gz` + `.fai` + `.gzi`) → whoever first needs one
  (owner deferred, 2026-07-18; §3.8). The shape is known: §4's from-byte-zero pass must run over
  a bgzf reader instead of raw bytes, the `.fai` `OFFSET`s address the *uncompressed* stream, and
  seeking needs the companion `.gzi` — htslib's `fai_load3` already carries the whole pattern
  (`faidx.c:587-710`, the `bgzf`/`.gzi` branch). Until then a compressed reference is a clean
  error, not a misread. This also unblocks ng's uncompressed-only reference stack as a whole, so
  it is likely done there, not here alone.
- ~~**Writing a `.fai`.**~~ **Built (§3.9)** — no longer deferred. `read_reference_info`
  returns `Vec<ContigInfo>` and `write_fai` writes them, so a reference that lacks an index can be
  indexed by ng rather than by a `samtools faidx` prerequisite
  ([`cli.rs:123-125`](../../../../src/pop_var_caller/cli.rs)).
- **Parallelising a *single* read.** The cache (§3.7) makes *many* reads cheap; it does not
  split *one* — a single genome read is serial. `verify_fasta_matches_psp_chromosomes` fans its
  per-contig MD5s across rayon ([`common.rs:208-238`](../../../../src/pop_var_caller/common.rs))
  and gets `max_chrom / threads` instead of the serial sum — but it can only do that
  **because** it seeks by `.fai` offset, which §3.3 rejects as circular. A from-byte-zero pass
  is inherently serial. → whoever measures one read and finds it hurts *even cached* (i.e. the
  first, uncached read is itself too slow); the answer would be a two-phase shape (index first,
  then hash in parallel), i.e. two reads, and it should be a measurement that buys it. The
  cache is the cheaper win and is why this is unlikely to be reached.

---

## 8. Open questions

*(Settled, 2026-07-17: **the check covers the line geometry**, so the pass writes its own
parse loop rather than using noodles' `sequence_reader()` (§4). The question was whether to
follow the standard tools' criterion; the source says they have none (§3.3), and the owner's
ruling is explicit: **duplicate contig names and FASTA↔`.fai` disagreement are refused, no
matter what any other tool does.** A check that `fold -w 60` defeats is not what §3.2 decided
to pay a genome read for.)*
- **Does the offset counting survive `noodles`' buffering?** (§4). The pass needs the byte
  offset reached after each definition line. Wrapping the file in a counting `BufRead` gives it
  *if* the reader consumes exactly what it parses and never reads ahead past a record boundary.
  That is plausible from `indexer.rs`'s own structure — it tracks `self.offset` the same way,
  by counting what each read returns — but it is a fact about a dependency, so **check it at
  implementation time**; the fallback is to parse the definition line with `BufRead::read_until`
  directly and not use `noodles`' reader at all, which is what the Indexer effectively does.
*(Settled, 2026-07-18: **`ReferenceInfo` carries the reconstructed index** — as `ContigInfo`
fields, not a separate list — and a `.fai` writer ships (§3.9). This reverses the earlier "no",
which rested on returning a *noodles* type in the public output; the geometry rides on
`ContigInfo`, ng's own plain struct, so the objection is gone and the data is already computed.
Owner's call.)*
*(Settled, 2026-07-17: **the digests ride on `ReferenceInfo`.** They are near-free once the
bases are being read, and they have named consumers — see §3.4. A separate
`reference_digest(&Path)` would read better and cost a second genome read, which decides it.
Owner's call; recorded here because the question is a natural one to re-ask.)*
