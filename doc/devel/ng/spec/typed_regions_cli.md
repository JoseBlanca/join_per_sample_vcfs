# `pop_var_caller_exp` — the experiment binary, and its typed-region subcommand

*Status: design spec (2026-07-17). **No code yet.** Companion to
[`typed_regions.md`](typed_regions.md), which settles the walk this command drives; nothing about the
walk's behaviour is re-decided here. Two things at once, because the first needs the second to be
worth building: **a second binary, `pop_var_caller_exp`, that is ng's command surface** (§2), and
**its first subcommand, `type-regions`** (§2.1 on). It is the first code outside `src/ng/` that ng's
step-3 plan permits (§5, T7).*

*Naming: **STR** in prose, `ssr` in code.*

---

## 1. What it is

**A subcommand that runs step 3's walk and writes the answer to a file** — in a binary of its own, so
the production CLI does not grow an experiment's knobs. In: a reference FASTA, an optional BED of
regions, and the walk's knobs. Out: one text file whose rows are the typed regions — contig, span,
kind, and for a microsatellite the motif and how many times it repeats.

Today the walk exists only as a library iterator that nothing outside `src/ng/` calls
([`region_typing/mod.rs:892`](../../../../src/ng/region_typing/mod.rs), and its only driver is a
`#[cfg(test)]` anchor, [`anchor.rs:125-138`](../../../../src/ng/region_typing/anchor.rs)). Its output
therefore cannot be looked at, diffed between settings, or handed to anything. **This command is what
makes step 3's output an artefact**.

**Goals.** Give ng a command surface (§2). Drive the shipping
stack (FASTA + `.fai` → `WindowedRefSeq` → `TypedRegionIterator`) end to end, streaming. Expose the
STR-detection knobs of `TypedRegionConfig` at their defaults, so a sweep needs no recompile (spec §5's
"every rule must be a parameter" reaching the command line). Emit a file whose **rows** carry what the
reference alone can say about each region — for a microsatellite that means motif, period, repeat count
and purity. It also carries the run parameters that made it: the reference, and every config value, so
two files are comparable and a file is reproducible (§3.4). Report the counts the walk already tallies,
honestly (T9).

**Non-goals.** **Deciding anything about a region** — this command types, it does not route or mask;
what a consumer does with a `Satellite` row is theirs (spec §1). **A locus catalog** — `ssr-catalog`
is production's; this is a *partition of the genome*, a
different artefact with a different reader. **Indexing** (§8). **Reading the file back into ng** —
nothing consumes it (§8). **A parallel walk** — the walk is single-threaded by decision (spec §7); the
one other thread is the reference verification `reference_info` runs in the background (T1), which the
walk never touches.

**It does not:** open an alignment file; sort (the walk emits in genomic order,
[`mod.rs:575`](../../../../src/ng/region_typing/mod.rs)); filter by kind (that is `grep`, §3.3);
collect the partition in memory (§6). It **may write a `.fai`**: if the reference has no sibling
index, it creates one so the next run starts faster (T1) — the only file it writes besides the output.

---

## 2. A second binary — `pop_var_caller_exp`

**Owner's call, 2026-07-17: ng's commands do not go in `pop_var_caller`.** The production binary has
seven subcommands ([`cli.rs:63-90`](../../../../src/pop_var_caller/cli.rs)), and every one is
something the lab runs on real data. ng's are a different kind of thing: they exist to answer
questions (spec §10), their knobs *are* the experiment.


[`Cargo.toml:129-131`](../../../../Cargo.toml):

```toml
[[bin]]
name = "pop_var_caller_exp"
path = "src/main_exp.rs"
```

`src/main_exp.rs` stays thin, the way `src/main.rs` is thin ([`main.rs:1-4`](../../../../src/main.rs));
the CLI lives in a new `src/pop_var_caller_exp/` mirroring `src/pop_var_caller/` — `cli.rs` for the
top-level `Parser` plus the `Subcommand` enum, then one module per subcommand owning its `Args`
struct, its `run_*`, and its `#[non_exhaustive]` error enum
([`mod.rs:6-31`](../../../../src/pop_var_caller/mod.rs) is the layout to copy). Plus `pub mod
pop_var_caller_exp;` in [`lib.rs`](../../../../src/lib.rs).

**It is a second binary in the same crate, and what that separates is the command surface, not the
code** — the spec should not oversell it. The `pop_var_caller` binary still links a library containing
ng, and both binaries link the same one. What splits is what a user can *invoke*, which is what the
complaint was about.

Staying in-crate is also what keeps the exp binary cheap: `pop_var_caller::common` is `pub(crate)`
([`mod.rs:8-9`](../../../../src/pop_var_caller/mod.rs)), so it can reuse
`DEFAULT_BUFFERED_IO_CAPACITY`, `basename` and `current_command_line` without a production edit.

**One thing not to do later: do not move `src/ng/` into a crate of its own.**

**The binary is named for experiments, not for ng**.

### 2.1 The `type-regions` subcommand

The shape copies `ssr-catalog`
([`ssr_catalog.rs:14-54`](../../../../src/pop_var_caller/ssr_catalog.rs)): `--reference` /
`--output` as bare `PathBuf`s, `--regions` optional, every knob `default_value_t` a named `pub const`
under `help_heading = "Advanced"`.

```
pop_var_caller_exp type-regions --reference ref.fa --output regions.tsv [--regions r.bed] [knobs]
```

Every knob's default is a named `pub const`, cited so nobody re-derives it. **The four rows tagged
§2.3 carry the short-read defaults the consts do not yet hold** (`--min-period`, `--max-str-len`,
`--bundle-threshold`, `--min-copies`): there the citation points at the const *to change*, and §2.3 says what
to — the consts still ship the catalog values until the CLI slice applies them (kept as record,
2026-07-19). The other rows already ship the value shown.

| flag | default | source |
|---|---|---|
| `--min-period` / `--max-period` | 1 / 6 | `DEFAULT_MIN_PERIOD`, `DEFAULT_MAX_PERIOD` ([`segment_criteria.rs:355-358`](../../../../src/ng/region_typing/segment_criteria.rs)); the STR-detection period range — one flag, §2.2 |
| `--max-str-len` | 100 | `DEFAULT_MAX_STR_LEN` ([`mod.rs:240`](../../../../src/ng/region_typing/mod.rs)); §2.3 |
| `--window-bp` | 100 000 | `DEFAULT_WINDOW_BP` ([`mod.rs:244`](../../../../src/ng/region_typing/mod.rs)) |
| `--bundle-threshold` | 30 | `DEFAULT_BUNDLE_THRESHOLD` ([`segment_criteria.rs:348`](../../../../src/ng/region_typing/segment_criteria.rs)); §2.3 |
| `--min-purity` | 0.8 | `DEFAULT_MIN_PURITY` ([`segment_criteria.rs:339`](../../../../src/ng/region_typing/segment_criteria.rs)) |
| `--min-score` | 0 | `DEFAULT_MIN_SCORE` ([`segment_criteria.rs:344`](../../../../src/ng/region_typing/segment_criteria.rs)) |
| `--min-copies` | `6,4,4,3,3,3` | `MinCopies::default()` ([`segment_criteria.rs:449-458`](../../../../src/ng/region_typing/segment_criteria.rs)); §2.3 |
| `--scan-match-reward` / `--scan-mismatch-penalty` / `--scan-min-copies` | 2 / 7 / 2 | [`tandem_repeat.rs:105-113`](../../../../src/ng/tandem_repeat.rs) |

**`--min-copies` is a table, not a number**, and no existing subcommand takes one. **Decided (owner,
2026-07-19): exactly six comma-separated values, one per period 1..6 — `6,4,4,3,3,3` — and a list of
any other length is a hard error** (a clap parse failure: usage, exit 2, before `run` is reached).
**No catch-all, no colon.** `MinCopies`'s `for_wider_periods` slot is structurally unreachable here —
the walk asserts the period ceiling ≤ `MAX_MOTIF_LEN` (6)
([`segment_criteria.rs:940-945`](../../../../src/ng/region_typing/segment_criteria.rs)), so no
period-7+ interval can reach the floor — so the parser fills it with a fixed inert `3` rather than
making the user type a value that decides nothing. The `value_parser` lives in the exp binary's own
`cli/parsers.rs` — production's
([`cli/parsers.rs`](../../../../src/pop_var_caller/cli/parsers.rs)) is the pattern to copy, not a
place to add to — and builds a `MinCopies` via `MinCopies::new`
([`segment_criteria.rs:419-430`](../../../../src/ng/region_typing/segment_criteria.rs)).

**`--max-str-len` must not be smaller than `--bundle-threshold`**, and the walk rejects the pair with an
error before it does any work — see T3. The CLI does not need its own check.

### 2.2 One period range

**`--min-period`/`--max-period` is one flag, default 1..=6, and it means exactly the range asked
for.** It sets [`SsrSegmentCriteria::periods`](../../../../src/ng/region_typing/segment_criteria.rs) —
the single range the walk both detects and classifies; there is no second, wider "scan range" behind
it. So `--min-period 3` types period-3 and longer, and nothing shorter leaks through. *(Why there is
only one range is step 3's decision, settled in [`typed_regions.md`](typed_regions.md) §2.2 — the CLI
exposes the result, it does not re-argue the mechanism.)*

**These are *detection* defaults.** The same parameters will later drive a separate routing /
empirical-stutter study, which sets its own (§2.3, §9); this command only describes the genome's
architecture — what the sequence is where. Keeping the two apart is the owner's rule (2026-07-18), and
it is why the period range needed a second look at all.

### 2.3 Defaults tuned for short-read data — changed in the struct, not just the CLI

**Owner, 2026-07-18: the library `DEFAULT_*` consts move, and the CLI inherits them.** These are ng's
own constants (`src/ng/`), so nothing in production changes. The reason for editing the struct rather
than the flag is the CLI's own principle — *`Default` = what the lab runs* — which is a lie if the
struct still says 1000 and only the flag says 100: every other caller of the library (tests, later ng
steps) would get the wrong number. So the source of truth moves, and this table follows it.

**`--max-str-len` 1000 → 100.** This one field is the **satellite cap** and the **scan margin**
at once (`TypedRegionConfig`, [`mod.rs:207-238`](../../../../src/ng/region_typing/mod.rs)), so 100
does two things:

- *A tract longer than 100 bp is a `Satellite`, not an STR.* With 150 bp reads, a read spans a tract
  plus an anchor each side only up to ~`read_len − 2·bundle_threshold` ≈ `150 − 60 = 90` bp. A tract longer
  than that cannot be genotyped from a single read whatever we label it, so past ~100 bp the STR route
  has nothing to offer and `Satellite` (*"don't look for loci here"*, spec §2.1) is the honest type.
  100 is a round number at that read-length limit, not a measured one — **soft, and the point of the
  knob is to sweep it** (spec §10).
- *The window fetches core ± 100 bp instead of ± 1000.* Cheaper, and still whole: no non-satellite
  repeat can exceed the cap, so the margin always covers one.

**`--bundle-threshold` 50 → 30.** 30 bp is more than enough unique sequence to anchor a short read to a
locus; 50 was inherited from GangSTR and unmeasured (spec §10 already flags it as the open question
*"how much flank does the analysis actually need?"*). It stays the bundle radius too (spec §2.4), so
lowering it also loosens what counts as a bundle — fewer neighbours are "too close". Soft, and swept.

**ng's `Default` is no longer the catalog's settings, and that costs nothing that matters** (owner,
2026-07-18: *"the golden standard is trf-mod; we can always compare with it"*). The trf-mod comparison
never ran at ng's `Default` — `the_walk_reproduces_the_golden_catalog_through_the_shipping_stack`
([`anchor.rs:254-265`](../../../../src/ng/region_typing/anchor.rs)) reads trf-mod's own parameters
from the `.cat` header and walks at *those*. So the golden-standard check is untouched; it just needs
to pin the last two knobs it currently inherits from `Default` (`min_copies`, `max_str_len`)
explicitly from the header, which is a two-line fix and exactly what `typed_regions.md` §8.1 already
asked for.

The only casualty is `default_matches_the_frozen_catalog_params`
([`segment_criteria.rs`](../../../../src/ng/region_typing/segment_criteria.rs)), which asserts
`SsrSegmentCriteria::default() == CatalogParams::default()`. That coupling is what these changes drop
on purpose — **delete the test**, it was bookkeeping, not the oracle. When this lands,
`typed_regions.md` §5.2 (*"v1 starts at the catalog's settings"*) gets the same one-line correction:
ng's `Default` is the lab's short-read values, and the trf-mod comparison pins trf-mod's explicitly.
Not done here — the defaults are not all final and there is no code yet.

**`--min-copies` `[10,5,4,3,3,3]` → `[6,4,4,3,3,3]`** (period 1..6; CLI syntax `6,4,4,3,3,3`, §2.1). The floor is the copy number at which a repeat starts to **stutter**, because a
repeat that does not stutter is one the generic SNP/indel caller handles fine — only a stuttering one
needs the STR route. The change is at the mononucleotide floor: **10 → 6**, so a homopolymer of ~6+
routes to STR handling, which is where Illumina's indel error rises (the germline-slippage threshold
is higher, ~9 units, but the read *artifact* onset is what a caller must react to). The rest tracks
the shorter-motif-stutters-more literature; periods 4–6 stay at 3. Every number is a **starting value,
soft and swept** (spec §10) — the shape is literature-backed, the exact floors are not measured. The
reasoning and sources are in [§9](#9-open-questions).

**`--min-period` 2 → 1** (`DEFAULT_MIN_PERIOD`, [`segment_criteria.rs`](../../../../src/ng/region_typing/segment_criteria.rs)).
The same stutter research decided this: mononucleotides stutter most, so they belong on the STR path,
and the `--min-copies` floor of 6 for period 1 is what a `--min-period 1` default exists to *use* — at
`--min-period 2` that floor is dead. So the default now classifies homopolymers of ≥ 6 bp as period-1 STR
loci. **Code consequence to flag for the implementer:** at the old default a homopolymer was `Generic`
(nothing classified period 1), and a test pins exactly that — `a_homopolymer_is_generic_and_does_not_take_its_neighbour_with_it`
([`mod.rs`](../../../../src/ng/region_typing/mod.rs)). Under the new default a qualifying homopolymer
is a period-1 `SsrSegment`, so that test must be re-pinned to `--min-period 2` (to keep testing the
"period 1 not classified" path) *and* a new one added for the default's "homopolymer is a locus"
behaviour. Genomes are dense in homopolymers, so this is the change that most alters what the partition
looks like — expect many more `SsrSegment` segments, which is the intended architecture description.

---

## 3. The output format

### 3.1 Decided: typed columns, with `members` as a JSON cell

BED's first three columns say *where*. This file must also say *what*, and for a microsatellite the
"what" is structured: motif, period, repeat count, purity. That is the whole reason for a format
section — the coordinates are the easy part.

**Decided (owner, 2026-07-19): a fixed set of typed columns, `.` for absent — and the one genuinely
nested field, a bundle's members, carried as a JSON array inside its cell.** The first four columns
(`chrom`, `start`, `end`, `kind`) are always filled, so the file is a valid BED to any tool that reads
the first three; the rest are filled only by the kinds that have them.

```
#chrom  start  end   kind        motif  period  copies  purity  members
chr1    0      999   generic     .      .       .       .       .
chr1    999    1029  ssr_locus   CAG    3       10.0    1.0     .
chr1    1029   1990  generic     .      .       .       .       .
chr1    1990   2040  ssr_bundle  .      .       .       .       [{"start":1990,"end":2010,"period":2},{"start":2020,"end":2040,"period":5}]
chr1    2040   4240  satellite   .      .       .       .       .
```

**Why typed columns.** Every column but `members` is a scalar, so the format is naturally tabular:
`awk '$4=="ssr_locus" && $6==3'` is the whole query language, the sweep §10 wants is a
`sort | uniq -c`, and pandas reads it (`read_csv(sep='\t', comment='#')`) as a flat, already-typed
frame. It stays a BED our own `--regions` accepts (§3.2), and it matches the one output precedent in
the tree — `ssr-catalog`'s `#`-header TSV ([`io.rs:32`](../../../../src/ssr/catalog/io.rs)).

**Why JSON in the `members` cell, and only there.** A bundle's members are the one field that is a
list rather than a scalar; everywhere else JSON would only wrap a value the column already holds.
Encoding *that* cell as a JSON array — rather than a bespoke `1990-2010:2` mini-language — means
`json.loads(row["members"])` gives the structure where it helps, with no mini-language for a consumer
to reimplement and without pushing every other column through `jq`. It is what GFF does with its
attributes column: tabular where the data is tabular, structured in the one cell that is not. The
array must serialise deterministically — fixed key order, no incidental whitespace — for §6's
byte-identity.

*Rejected — a pure JSON row per line (`chrom start end kind {json}`).* Structured detail with no
padding, but the whole file pays for one cell's benefit: not greppable (every question needs `jq`), no
`--regions` round-trip (§3.2), no reader of that shape in the repo, and pandas reads it *less* readily
than the TSV (a nested column to normalise). The Python-analysis case for it is a wash — a typed TSV
loads as cleanly — so the members cell alone did not justify reshaping every row.

### 3.2 Decided: plain text, not bgzip

`ssr-catalog` writes bgzip ([`io.rs:238-271`](../../../../src/ssr/catalog/io.rs)) and this does not,
because the two artefacts are read differently. The catalog is **random-access**: the STR caller asks
it for one locus, which is what tabix and a bgzf frame are for. This file is read **end to end**, by
an experimenter comparing two configs. Plain text keeps `awk`, `diff` and `head` working with no
pipe, and it keeps one property that is otherwise lost: **the file is a BED our own `--regions`
accepts**, because `RegionSet::from_bed_path` opens a plain file
([`regions.rs:150-156`](../../../../src/regions.rs)) and ignores `#` comments and every column past
the third ([`regions.rs:99-115`](../../../../src/regions.rs)). So
`grep -P '\tsatellite\t' out.tsv > mask.bed` is a mask, and feeding the file back to `--regions`
re-walks exactly what it describes. Size is not an argument either way: a whole-genome partition is
roughly two rows per finding, so tomato is tens of megabytes.

*Rejected: bgzip for consistency with `ssr-catalog`.* It buys an index nothing asks for (§8) and
costs the round-trip above.

### 3.3 Decided: every base gets a row; no `ref_seq` column

**`Generic` rows are emitted.** The file *is* the partition — spec §2.3's spine is that
concatenating the regions reconstructs what was asked for — and that property is checkable in the
file only if the generic runs are in it. A `--kinds ssr_locus` filter would be `grep` with extra
steps, and it would silently destroy the invariant the file's own test asserts (§7).

**No `ref_seq` column**, unlike the catalog ([`io.rs:32`](../../../../src/ssr/catalog/io.rs)). The
catalog stores tract+flank bases because its consumer genotypes loci without a FASTA in hand. This
file's consumer has the FASTA — it is a required flag — and every base is derivable from the
coordinates. Purity is the one field that is not derivable cheaply, so purity stays.

### 3.4 The header block

The `ssr-catalog` shape ([`io.rs:56-69`](../../../../src/ssr/catalog/io.rs)): `## key: value` lines,
then one `#`-prefixed column header, then rows. It carries the tool version, the reference path, and
**every knob's value** — the config that produced the file, so two files are comparable and any one is
reproducible. That is the header's whole job, and it is a goal, not a nicety (§1): a sweep compares
files, and a file that cannot say what config made it is noise.

```
## tool: pop_var_caller_exp type-regions
## version: <crate version>
## reference: /path/to/ref.fa
## min_period: 1
## max_period: 6
## max_str_len: 100
## window_bp: 100000
## bundle_threshold: 30
## min_purity: 0.8
## min_score: 0
## min_copies: 6,4,4,3,3,3
## scan_match_reward: 2
## scan_mismatch_penalty: 7
## scan_min_copies: 2
#chrom	start	end	kind	...
```

**Every config value appears, resolved** — not just the flags the user typed, or a file run with
defaults would record nothing and be incomparable to one that set them explicitly. The
`TypedRegionConfig` the run built is the source of truth; serialise *it*, so a knob added later cannot
be forgotten here — the one `--min-period`/`--max-period` range (§2.2) among them.

**`window_bp` is in the block but is not a comparison key**: it
is a memory knob that must not change the output (spec §2.3), so two files that differ only in
`window_bp` should be identical below the header. Recording it keeps the file reproducible without
implying it is part of the experiment.

**No `reference_md5`**, and this is a departure worth its line. The catalog carries one
([`catalog/mod.rs:205-207`](../../../../src/ssr/catalog/mod.rs)) because it streams the whole FASTA
anyway. Here the digest *is* computed — the background verification (T1) reads the whole FASTA — but it
lands only when `VerificationHandle::join` returns, which is *after* the walk, whereas the header is
written *before* it. Putting the MD5 in the header would force the header, and so the walk, to **block
on that background pass** — defeating the parallelism T1 exists for. So it stays out: the reference
path and its knobs identify the run, the header stays reproducible from path + config alone, and a
caller wanting the MD5 can `md5sum` the FASTA.

---

## 4. Records, per kind

`TypedRegion` is `{ region: GenomeRegion, kind: RegionKind }`
([`mod.rs:145-149`](../../../../src/ng/region_typing/mod.rs)). The mapping to columns:

| kind | `chrom`/`start`/`end` | the rest |
|---|---|---|
| `Generic` | the region | all `.` — the region *is* the whole claim |
| `Satellite` | the region | all `.` — likewise |
| `SsrSegment(SsrSegment)` | the region, which **is** the tract ([`mod.rs:552-561`](../../../../src/ng/region_typing/mod.rs)) | `motif` = `SsrSegment::motif()` bytes, `period` = `motif.period()`, `copies` = derived (T5), `purity` = `SsrSegment::purity_fraction()` |
| `SsrBundle { tracts }` | the hull | `members` = a JSON array, one object per tract: `{"start":…,"end":…,"period":…}` (T4 on the coordinates) |

`SsrSegment`'s accessors are [`segment_criteria.rs:284-330`](../../../../src/ng/region_typing/segment_criteria.rs);
`RepeatInterval`'s fields are
[`tandem_repeat.rs:207-217`](../../../../src/ng/tandem_repeat.rs). `SsrSegment::chrom()` already carries
the contig **name**, so no row needs the `ContigId` → name lookup that `GenomeRegion` would
otherwise force — except that `Generic`, `Satellite` and `SsrBundle` carry no `SsrSegment`, so the command
needs the `ContigList` to name them anyway. Keep one table for the whole run (T8).

---

## 5. What will bite you

Every item here came from opening the file it cites.

**T1 — read the reference through `reference_info::read_reference_verifying_or_creating_fai`, and
join the verification handle before you rename the output.** This one call
(`(cache: &Arc<ReferenceInfoCache>, fasta: PathBuf) -> (Arc<ReferenceInfo>, Option<VerificationHandle>)`,
ng's `reference_info`, landing from `ng-contig-table`; owner, 2026-07-19) does exactly what this
command wants, in two cases (owner, 2026-07-19):

- ***`.fai` present*** → it reads the index in the foreground (names, order, lengths — cheap, so the
  walk can start at once) **and verifies the FASTA against that index on a background thread**,
  returning `Some(handle)`. The walk runs while the check runs; a stale `.fai` (one whose FASTA no
  longer matches) is caught there, not by producing a wrong partition.
- ***`.fai` absent*** → it reads the FASTA once, **writes the sibling `.fai`** (so the next run takes
  the fast path), and returns `None`. A `.fai`-write failure is fatal (`ReferenceInfoError::FaiWrite`).

`.contig_list()` gives the `ContigList` for the walk. **The trap is the handle**: `VerificationHandle`
is `#[must_use]`, and its verification is only worth running if it can stop a bad run — so
`handle.join()?` **before** the atomic rename (§6), never after. Skip it and a partition built against
a stale `.fai`'s contig table is renamed into place, silently wrong. Joining at the end costs almost
nothing: the background pass overlapped the whole walk. Do **not** reach for the bare
`ReferenceSource::Fai`/`Fasta` arms — `Fai` alone skips the verification, `Fasta` blocks the walk on
the whole-genome read.

**T2 — the contig length narrows, and `as` is not the answer.** `ContigEntry::length` is `u64`
([`fasta/mod.rs:37-43`](../../../../src/fasta/mod.rs)); `ContigBounds::length` is `u32`
([`regions.rs:57-63`](../../../../src/regions.rs)). The only precedent casts —
`length: e.length as u32` ([`anchor.rs:98`](../../../../src/ng/region_typing/anchor.rs)) — and it is
a test. B2 deleted exactly this shape from ng and recorded the rule: **ng fails rather than folds**
(`typed_regions.md` §4). Use `u32::try_from` and surface an error. A >4 Gb contig is unrepresentable
to `RegionSet` whatever we do, so failing is the honest answer, not a limitation we are choosing.

**T3 — do not re-validate `--max-str-len` against `--bundle-threshold` in the CLI; just propagate.**
`TypedRegionIterator::over_regions` refuses the pair itself, with
`TypedRegionError::MarginNarrowerThanFlank` and before any work
([`mod.rs`](../../../../src/ng/region_typing/mod.rs)). It used to `assert!` — right when the knobs
were library-only, wrong the moment they became flags, since user input must not panic. `ssr-catalog`
checks its own equivalent pair in the CLI
([`ssr_catalog.rs:78-83`](../../../../src/pop_var_caller/ssr_catalog.rs)) and is **not** the precedent
to copy here: a second check in the CLI is a second place for the rule to drift, and this one already
carries both numbers in its message.

**T4 — one row, two coordinate systems, in opposite directions.** The hull is a `GenomeRegion`:
**1-based inclusive** ([`types.rs:55-59`](../../../../src/ng/types.rs)). Its member tracts are
`RepeatInterval`s: **0-based half-open**, and already re-based to contig coordinates by the walk
([`mod.rs:1458-1465`](../../../../src/ng/region_typing/mod.rs);
[`tandem_repeat.rs:198-217`](../../../../src/ng/tandem_repeat.rs) explains why they stay 0-based).
The output is BED — 0-based half-open. So in the same row the region's start needs `-1` and the
members need **nothing**. Write the conversion once, in one helper, or it will be written twice and
one of them will be wrong.

**T5 — "the number of repeats" is not stored anywhere.** `SsrSegment` has no copies accessor; the
catalog's rule is that everything but the seven stored columns is *derived*
([`io.rs:11-13`](../../../../src/ssr/catalog/io.rs)). It is `tract_len() / period()`
([`segment_criteria.rs:312-324`](../../../../src/ng/region_typing/segment_criteria.rs)), and it is **fractional in
general** — a tract is rarely a whole number of motifs, which is part of what purity measures.
trf-mod's own BED calls the fractional value `copyNum`
([`trf.rs:10-17`](../../../../src/ssr/catalog/trf.rs)), so a fractional column is the field's
convention, not a wart. Say in the header which it is.

**T6 — `Motif` has no `Display`.** Only a hand-written `Debug` that prints `Motif("CAG")`, quotes
included ([`segment_criteria.rs:184-193`](../../../../src/ng/region_typing/segment_criteria.rs)). Use the catalog's
idiom ([`io.rs:163-165`](../../../../src/ssr/catalog/io.rs)): `motif()` returns by value (it is
`Copy`), so bind it before borrowing, then `str::from_utf8`.

**T7 — the B2 guard grep will fail on this commit, and it must be re-aimed, not deleted.** The
step-3 plan's standing rule was "nothing outside `src/ng/`", and B2 hardened it into
`rg 'use crate::ng' src/ --glob '!src/ng/**'` matching nothing (`typed_regions.md` §9). But that grep
is a proxy for the rule that actually matters — **production must not depend on ng** — and it was
exact only while ng had no driver. `src/pop_var_caller_exp/` is a driver, and its dependency points
the safe way: the exp binary depends on ng, ng gains nothing, production still depends on nothing in
ng. So the guard becomes `rg 'use crate::ng' src/ --glob '!src/ng/**' --glob
'!src/pop_var_caller_exp/**'`, and **the excluded path is now load-bearing**: it is what says the
`pop_var_caller` tree is the guarded one. Widen it once more and the guard is a comment.

Mechanically: `region_typing` is not re-exported from `ng/mod.rs`
([`ng/mod.rs:27-31`](../../../../src/ng/mod.rs) re-exports only `ref_seq` and `types`), so either add
the re-export or name `crate::ng::region_typing::…` in full.

**T7a — `format_error_chain` is private to `src/main.rs`** ([`main.rs:30-44`](../../../../src/main.rs)),
and `main_exp.rs` needs the same rendering or the two binaries report the same error differently. It
is 15 lines with its own `#[cfg(test)]` module ([`main.rs:76+`](../../../../src/main.rs)).
**Leaning: hoist it into the library** (a pure move, no behaviour change, and both binaries then
share one rule) rather than copy it — two copies of an error-rendering rule that must agree is how
they stop agreeing. It is the one production file this work touches, and it touches it by moving a
function, not by changing what the production binary does. Rejected: copy the 15 lines (cheap now,
and the copies diverge the first time someone fixes the substring heuristic in one of them).

**T8 — build one contig table and use it twice, or the walk's error is unreadable.**
`over_regions` is fallible because `GenomeRegions` was validated against *a* contig table and nothing
ties it to the reference's ([`mod.rs:892-905`](../../../../src/ng/region_typing/mod.rs)). If they
disagree the failure is `RefSeqError::UnknownContig(ContigId)`, which renders as
`unknown ContigId(7)` — an index, not a name ([`ref_seq.rs:52-53`](../../../../src/ng/ref_seq.rs)) —
arriving through `TypedRegionError::Reference`
([`mod.rs:1640-1641`](../../../../src/ng/region_typing/mod.rs)). Feed the *same* `ContigList` to
`WindowedRefSeq::new` ([`ref_seq.rs:524`](../../../../src/ng/ref_seq.rs)) and to `GenomeRegions` and
that error is unreachable — but only if you do.

**T9 — four of the five rejection counters are structurally zero.** `RejectionCounts` has five fields
([`segment_criteria.rs`](../../../../src/ng/region_typing/segment_criteria.rs)), and the walk reaches exactly **one**
of the gates behind them: `FlankClamped` (pinned by
`segment_criteria::the_walk_reaches_only_one_of_classifications_five_gates`). Print all five as a flat list and
`purity: 0` reads as *"no impure tracts in this genome"* — a wrong answer wearing a measurement's
clothes. **None of these zeroes is a property of classification**: each is caused by a stage *upstream* of
the gate — the scanner's maximal-scoring segments for `Purity` and `NoCleanTrim`, the pre-filter for
`CopyFloor` and `Compound` — so which are zero moves whenever those stages do. Do not hard-code the
list into the summary's prose. Two more facts it must not misstate: `rejected_by_reason` does **not**
partition `repeat_bp_with_no_locus` ([`mod.rs:290-298`](../../../../src/ng/region_typing/mod.rs)), and
`repeat_bp_with_no_locus` is exact only once the iterator is exhausted
([`mod.rs:286-289`](../../../../src/ng/region_typing/mod.rs)) — which this command does, so it is exact
here.

**T10 — `--regions` costs time, not correctness, and the cost is unbounded by the BED.** The scan set
is **whole contigs**: a BED chooses what is emitted, never what is scanned
([`typed_regions.md`](typed_regions.md) §2.5). So every finding comes back the same object the
whole-genome run reports — but a 10 kb `--regions` span on a 90 Mb chromosome **scans all 90 Mb**.
Peak memory is unchanged (same window); wall time is not, and it is a function of the contig, not of
what was asked for. That is a surprising bill for a small BED and the help text should say so, rather
than let a user infer that a narrow `--regions` is a quick run.

---

## 6. Cross-cutting

**Errors.** The house pattern, inherited whole: a `#[non_exhaustive]` `thiserror` enum per subcommand
and `run_typed_regions(&args) -> Result<(), TypedRegionsCliError>`
([`ssr_catalog.rs:56-74`](../../../../src/pop_var_caller/ssr_catalog.rs)); `main_exp.rs` walks the
source chain and exits 1, the way `main.rs` does ([`main.rs:30-73`](../../../../src/main.rs)) — see
T7a on sharing that walk rather than copying it. No anyhow in the CLI path. Variants this one needs,
each a `#[from]` where it can be: `ReferenceInfoError` from reading, verifying, or writing the reference index (T1), the contig-length narrowing (T2), `BedError`
([`regions.rs:315-379`](../../../../src/regions.rs)), `TypedRegionError`, and output I/O. **Not** the
`--max-str-len` / `--bundle-threshold` pair — `TypedRegionError` already carries it (T3).

**Nothing the user can type may panic.** Every knob is a flag now, so each of the walk's config
guards is reachable from a typo, and a typo deserves a message and exit 1 — not a backtrace. That is
why T3's guard is an error. When the remaining knobs grow guards, they follow the same rule; a
`debug_assert` is doubly wrong here, since a sweep runs in `--release` (spec §10) and would get the
wrong *answer* rather than the panic.

**Memory — stream, do not collect.** `ssr-catalog` collects every locus before writing
([`catalog/mod.rs:226-228`](../../../../src/ssr/catalog/mod.rs)) because it sorts and MD5s. This
command must not: the walk's whole design is that it holds ~102 kb and three coordinates
(`typed_regions.md` §6), and `for r in iter { writer.write(r?)?; }` keeps that promise. Collecting a
whole-genome partition puts back exactly the memory the walk exists to avoid. It also needs no sort —
the walk emits in genomic order ([`mod.rs:575`](../../../../src/ng/region_typing/mod.rs)), the same
contract `CatalogWriter` states ([`io.rs:236-237`](../../../../src/ssr/catalog/io.rs)).

**A truncated output file is silently valid, so write atomically.** Two precedents exist: `pileup`
writes `.tmp` then renames ([`cli.rs:621-643`](../../../../src/pop_var_caller/cli.rs)); `ssr-catalog`
writes in place ([`catalog/mod.rs:221-230`](../../../../src/ssr/catalog/mod.rs)). Take `pileup`'s,
because of a property peculiar to this artefact: a partition cut off halfway through chromosome 7 is
indistinguishable from a complete partition of a smaller genome — nothing in the file says how far it
was supposed to go. A half-written catalog at least fails its own reader's header checks.

**The reference verification joins into that same barrier.** The rename is the commit point, so
`VerificationHandle::join()?` goes **just before it** (T1): the walk streams to `.tmp` while the FASTA
is checked in the background, and only a clean join renames `.tmp` into place. A verification failure —
a stale `.fai` — aborts the run with the partition still under its temp name, never published. This is
the one place the two threads meet, and it meets at the commit, which is exactly right.

**Determinism** is the regression anchor. The walk is a pure function of reference, regions and
config (`typed_regions.md` §7), so the same three inputs must give a byte-identical file — including
the `##` block, which means the key order is fixed
([`io.rs:56-69`](../../../../src/ssr/catalog/io.rs) already does this deliberately) and no timestamp
goes in it. *(`ssr-catalog`'s header has a `date` field ([`io.rs:60-68`](../../../../src/ssr/catalog/io.rs));
ours should not, for this reason.)*

**Progress.** There is no logging framework anywhere in `src/` — the convention is `eprintln!` to
stderr, an announce line at the start and a `key: k=v` summary at the end
([`cli.rs:369-376`, `cli.rs:822-859`](../../../../src/pop_var_caller/cli.rs)). Follow it; the summary
is `TypedRegionCounts`, subject to T9.

---

## 7. Tests

`region_typing/anchor.rs` already drives the shipping stack from a real multi-contig FASTA with a
`.fai` on disk, through `WindowedRefSeq` and the public iterator
([`anchor.rs:125-138`](../../../../src/ng/region_typing/anchor.rs)). This command's tests stand on
that fixture and test **the file**, not the walk again:

- **Round-trip.** Parse the written rows back and assert they are the iterator's own output, field
  for field. This is what pins T4's two-directions conversion; a fixture with a bundle is required,
  or the member coordinates are never exercised.
- **The partition invariant, in the file.** Concatenate the rows' spans and assert they reconstruct
  the requested regions exactly (`typed_regions.md` §2.3, now over the artefact). This is the one
  test that would catch a writer dropping generic rows, a header line miscounted as a row, or an
  off-by-one at T4.
- **Determinism.** Two runs, byte-identical (§6).
- **`--regions` read-back.** Feed the command its own output as `--regions` and assert it walks the
  same spans (§3.2's claimed property, which is otherwise an argument).
- **The flag-pair guard is a CLI error, not a panic** (T3). Mutation-verify it: delete the check and
  the test must fail with a panic, not pass.

Home: `src/pop_var_caller_exp/typed_regions.rs`'s own `#[cfg(test)]` module, per the per-subcommand
convention. Not `tests/` — the fixture builders it wants are in-crate, which is the same reason E3's
anchor is in-crate (§2).

---

## 8. Deferred, with a home

- **A tabix/CSI index**, and with it bgzip (§3.2) → whoever first needs random access into this file;
  the `.cat` reader/writer is the shape to copy ([`io.rs:238-333`](../../../../src/ssr/catalog/io.rs)).
- **A reader for this format** → the first ng step that consumes typed regions from disk rather than
  from the iterator. Today nothing does, and a reader written now would be a guess at its shape.
- **ng's other drivers** — a `pileup` subcommand, the gatherer, the sweep runners → their own specs.
  They land in this binary; that is what it is for (§2).
- **Deriving the file from something other than a FASTA** (e.g. typing regions from an existing
  catalog) → not a thing anyone has asked for; noted only so its absence is a decision.

---

## 9. Open questions

- **Which output format? — resolved 2026-07-19 (§3.1):** typed columns, with a bundle's `members`
  carried as a JSON cell. A pure-JSON row was rejected — it forfeits the terminal/`awk` workflow and
  the `--regions` round-trip for one cell's benefit, and the Python-analysis case for it is a wash
  against a typed TSV that pandas reads directly.
- **What happens to the subcommand names when a second experiment arrives?** Today
  `pop_var_caller_exp type-regions` is unambiguous because ng is the only inhabitant (§2.1). A second
  lab wanting its own `type-regions` would force either a prefix (`ng-type-regions`) or nesting
  (`ng type-regions`). **Leaning: cross that bridge then** — renaming one subcommand in an
  experiment binary costs nothing, and prefixing now taxes every ng command for a second lab that may
  never exist. Recorded so the rename is a decision rather than a surprise.
- **Should `pop_var_caller_exp` be excluded from release artefacts? — resolved 2026-07-19: leave it
  in.** It builds by default with `cargo build`; the tree ships no packaging today, so there is
  nothing to exclude it from, and the §2 separation is about the *command surface*, not the shipped
  file list. Whoever first cuts a release can revisit.
- **Should the run summary print the four structurally-zero counters at all?** (T9). **Leaning:
  print them, labelled** — the alternative is hiding a column that will become live the moment §10
  moves a knob, and a hidden zero is worse than an explained one. What would settle it: the first
  reader who misreads one.
- **`--min-copies` default — resolved 2026-07-18, and still soft.** The knob is the copy-number floor
  per period: below it, a tract is not classified as an STR. What it *should represent* (owner) is **the
  copy number at which a repeat starts to stutter**, because a repeat that does not stutter is one the
  generic SNP/indel caller handles fine, and only a stuttering one needs the STR route. **Decided:
  `[6,4,4,3,3,3]`, catch-all 3** — the sequencing-stutter reading (§2.3). It stays a starting value,
  swept against GIAB in spec §10's period × length frontier.

  **The literature is why this is not simply "the owner's guess, implemented", and the reasoning is
  worth keeping** because the two candidates pointed opposite ways for the mononucleotide floor:

  - *Germline slippage* (Lai & Sun 2003, MBE — the direct *rate vs number-of-repeat-units* study)
    puts the threshold at which slippage rises above background at **9 units for mononucleotides, 4
    for di–hexa**. So by *polymorphism* onset a mononucleotide needs the *most* copies before it
    varies — the opposite of the owner's "smaller for shorter motifs", and close to the old
    `[10,5,4,3,3,3]`.
  - *PCR / sequencing stutter* is a different quantity: the artifact in the reads, which is what makes
    the generic caller choke, and for Illumina it is dominated by homopolymers well below 9 units (a
    poly-A of 6–8 already errors). By *read-artifact* onset the mononucleotide floor should be lower —
    which recovers the owner's intuition, by a different road.

  A caller routes on "will the generic path mis-handle this?", which is the read artifact, so the
  mononucleotide floor drops to **6**. The germline number (9) is recorded because it is the cleaner
  citation and because the sweep may show it is the better routing threshold after all.

  *Sources (2026-07-18):* Lai & Sun 2003, *Relationship Between Microsatellite Slippage Mutation Rate
  and the Number of Repeat Units*, Mol. Biol. Evol. 20(12):2123 —
  [oup.com](https://academic.oup.com/mbe/article/20/12/2123/978562) (threshold 9 units for mono, 4
  for di–hexa; slippage rate rises exponentially with repeat number); repeat-unit-length vs stutter
  strength (shorter motifs stutter more; tetra N−4 ≈ 6–7 %, penta N−5 < 1 %) —
  [PMC6412005](https://pmc.ncbi.nlm.nih.gov/articles/PMC6412005/), a direct in-vitro stutter-noise
  model. The germline-vs-sequencing-stutter distinction is ours, not either paper's.
