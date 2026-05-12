# BAQ — per-read Base Alignment Quality

Implementation plan for the BAQ stage of Stage 1 of the calling
pipeline. Specified in
[calling_pipeline_architecture.md §"Per-read likelihood quality"](../specs/calling_pipeline_architecture.md)
and §"Why BAQ in-process rather than requiring `samtools calmd -r`".

BAQ is the local-realignment HMM from
[Heng Li 2011](https://academic.oup.com/bioinformatics/article/27/8/1157/227268)
that computes, per base, the posterior probability that the base is
correctly aligned, and caps effective base quality at
`min(BQ, BAQ)`. The capped values are the BQ that feeds the
per-read likelihood (`max(ln_BQ, ln_MQ)`) downstream.

The walker is already wired to consume BAQ-capped qualities through
[`PreparedRead.bq_baq`](../../src/per_sample_caller/pileup/mod.rs#L232);
today that field carries raw BAM BQ as a placeholder (see the
docstring at
[mod.rs:229-231](../../src/per_sample_caller/pileup/mod.rs#L229-L231)).
This slice fills in the capping.

## Scope

In:

- A `baq` module under
  [`src/per_sample_caller/`](../../src/per_sample_caller/) containing
  a pure HMM-Glocal port and a per-read driver that turns a
  `MappedRead` plus a reference slice into a `PreparedRead` with
  `bq_baq = min(BQ, BAQ)`.
- A `BaqConfig` (parameters: gap-open probability, gap-extension
  probability, band half-width, minimum input BQ floor) and a
  `BaqEngine` that owns the scratch buffers so a per-read call does
  not allocate on the hot path.
- Integration into the existing CRAM → pileup pipeline: where today
  `MappedRead` flows directly into the walker (or its placeholder
  bridge), it will flow through `BaqEngine::process` first.
- A new `FilterBucket::BaqRejected` skip counter, consistent with
  the other per-sample skip categories ("low MAPQ, unmapped,
  first-or-last-CIGAR indel" — see
  [architecture spec §410](../specs/calling_pipeline_architecture.md#L410)).
- Parity tests against htslib's `probaln_glocal` driven from
  pre-generated fixtures (option B in §"Cross-checking the port"
  below).
- Rayon-parallel BAQ across a batch of reads at the pipeline seam.

Out:

- Re-deriving the HMM math. We port and cross-check, we do not
  rewrite the algorithm. The paper plus the two reference
  implementations are the spec.
- Tuning parameters against truth sets. The defaults are inherited
  (see §"Parameters: htslib vs GATK"). A future calibration study
  could revisit them; out of scope here.
- Supporting `BAQ_EXTEND` (extended BAQ — htslib's `-E` mode). The
  architecture spec commits to `min(BQ, BAQ)` capping, not extended
  BAQ. Implementing extended BAQ would be a follow-up, not a default.
- Reading or writing the BAM `BQ` / `ZQ` tags. We are in-process —
  we compute BAQ fresh on every run and feed it into `bq_baq`
  directly. The tags exist for `calmd`'s round-trip-to-disk model,
  which we explicitly avoid.
- Long-read tuning beyond what falls out of using htslib's
  short-vs-long auto-select (see §"Algorithm port" below). The
  walker already gates on `MAX_RECORD_SPAN` (default 5000 bp)
  upstream.

## Current data flow

```
CRAM decoder → MappedRead ───────────────────────────→ PreparedRead → walker
                ▲                                            ▲
                │                                            │
                cram_input.rs                              today: identity-with-raw-BQ
                                                          (mod.rs:229-231 placeholder)
```

After this slice:

```
CRAM decoder → MappedRead → BaqEngine::process ─────→ PreparedRead → walker
                              │  needs:                      ▲
                              │  - ref slice for read span   │
                              │  - rayon for parallelism     bq_baq = min(BQ, BAQ)
                              └─ writes:
                                 PreparedRead.bq_baq[i] = min(qual[i], q[i])
                                 (q[i] from probaln_glocal output)
```

The reference slice already exists upstream: the per-sample caller
holds a FASTA via
[`per_sample_caller/ref_fetcher.rs`](../../src/per_sample_caller/ref_fetcher.rs)
which the pileup walker uses through
[`RefSeqFetcher`](../../src/per_sample_caller/pileup/mod.rs).
The BAQ stage shares that fetcher; no new IO machinery is needed.

The `MappedRead` → `PreparedRead` conversion currently lives at the
pileup walker's input boundary. After this slice it moves into the
BAQ stage and the conversion stops being trivial: it carries the
capped quality plus the existing fields.

## Decisions to make

### Parameters: htslib vs GATK

The architecture spec asserts that samtools/htslib and GATK
"converge on the same parameters"
([architecture spec §289-295](../specs/calling_pipeline_architecture.md#L289-L295))
and quotes `DEFAULT_GOP = 40` → P(gap) = 1e-4 from
[`BAQ.java`](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/baq/BAQ.java).
**They do not actually agree on the gap-open probability.** Reading
htslib's
[`realn.c:115`](../../htslib/realn.c#L115):

```c
probaln_par_t conf = { 0.001, 0.1, 10 }; // Illumina
```

— `d = 0.001` (P(gap) = 1e-3, GOP-equivalent 30 in Phred), not
`1e-4`. The `10` is overwritten to `bw = 7` (or wider, for
indel-induced length-diff overflow) at
[realn.c:195-197](../../htslib/realn.c#L195-L197):

```c
bw = 7;
if (abs((xe - xb) - (ye - yb)) > bw)
    bw = abs((xe - xb) - (ye - yb)) + 3;
```

So the agreement is real for `e = 0.1` and `bw = 7`, and the
disagreement is real for `d`: samtools-default uses **1e-3**, GATK
uses **1e-4**.

Recommendation: **pick samtools' 1e-3**. Reasons, in order:

1. samtools' `calmd -r` is what bcftools mpileup / bcftools call
   pipe their reads through by default. The architecture spec
   already commits to "calibration parity with bcftools" elsewhere
   (§"Why min indel quality"); the same logic applies here.
2. htslib is the canonical C implementation; GATK explicitly
   describes itself as a port of it
   ([BAQ.java:174-175](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/baq/BAQ.java#L174-L175)).
   When a port and its original disagree on a default, the original
   wins by default unless there is a reason to prefer the port.
3. The disagreement should be noted in a follow-up review of the
   architecture spec; that paragraph's claim of agreement is wrong
   and should either be corrected to "agree on `e`, `bw`,
   `minBaseQual`; disagree on gap-open probability" or replaced
   with the rationale for whichever value we pick.

Concrete values, locked by this plan:

| parameter | value | source |
|---|---|---|
| `d` (gap-open prob) | `0.001` | htslib realn.c:115 |
| `e` (gap-extension prob) | `0.1` | both |
| `bw` (band half-width) | `7`, expanded to `abs((xe-xb)-(ye-yb)) + 3` if the indel-induced length difference exceeds 7 | htslib realn.c:195-197 |
| `min_input_bq` (floor on input BQ before HMM) | `4` | GATK BAQ.java:76-78; htslib does not floor explicitly, but in practice Q4 is the lowest BQ Illumina emits |
| long-read auto-select | `d = 1e-7, e = 0.1` for `l_qseq > 1000` | htslib realn.c:117-130 |

The long-read branch matters if/when long-read CRAMs become a
target; until then it is dead code path that costs nothing.

### Where the BAQ pass runs

Two viable seams:

- **(a) Between CRAM decode and the walker.** Per-batch BAQ pass
  parallelised across reads, then the walker consumes
  already-capped `PreparedRead`s as it does today.
- **(b) Inside the walker, on read admission.** BAQ runs serially
  as part of `admit_read`, with the walker holding the rayon pool
  open for batched calls.

**Recommendation: (a).** The walker is already a sequential state
machine; mixing the embarrassingly-parallel BAQ pass into it would
force batch-then-drain semantics that the walker is not built for.
The architecture spec also frames BAQ as "applied inline and
parallelised across reads"
([§135-137](../specs/calling_pipeline_architecture.md#L135-L137))
and the cram-input plan already calls BAQ out as a separate slice
between filter and walker
([per_sample_caller_cram_input.md:46-49](per_sample_caller_cram_input.md#L46-L49)).
(a) matches both.

The walker keeps treating `bq_baq` as opaque (architecture spec
§"Why BAQ in-process" — the walker does not re-apply BAQ).

### Re-compute vs. trust existing `BQ` tag

htslib has three modes
([htslib/sam.h:2140-2196](../../htslib/sam.h#L2140-L2196)):
`BAQ_APPLY` (use cached `BQ` if present, otherwise compute),
`BAQ_REDO` (recompute even if `BQ` is present), and combinations
with `BAQ_EXTEND`.

**Recommendation: always recompute.** We do not see the BAM `BQ`
tag in `MappedRead` at all today, and the architecture spec frames
BAQ as a fresh per-run computation rather than a cached tag. This
sidesteps the `bq` ⇄ `zq` tag rewriting in
[`realn.c:138-177`](../../htslib/realn.c#L138-L177), which is
purely a calmd round-trip concern.

If a user wants to trust a pre-computed BAQ, they can run
`samtools calmd` upstream and we will silently ignore the `BQ`
tag they wrote. That is a non-feature, not a regression: their
upstream BAQ is at most as good as ours, and ours runs in-process
on every read anyway.

## Algorithm port

The core algorithm is the forward-backward HMM in
[`probaln.c:probaln_glocal`](../../htslib/probaln.c#L77). The
function takes a reference slice, a query, a per-base quality
array, and parameters; it writes an alignment-state vector and a
per-base posterior-error array.

The port has two layers:

### Layer 1 — pure HMM (`baq::probaln_glocal`)

A direct Rust translation of the htslib function. Signature:

```rust
pub(super) fn probaln_glocal(
    ref_seq: &[u8],      // 0/1/2/3/4 encoding (A,C,G,T,N)
    query:   &[u8],      // 0/1/2/3/4 encoding
    iqual:   &[u8],      // Phred quality, floored to `min_input_bq` by caller
    cfg:     &BaqConfig, // d, e, bw
    state:   &mut [i32], // out: alignment state per query base, length = query.len()
    q:       &mut [u8],  // out: Phred posterior error per query base, length = query.len()
) -> Result<(), BaqOverflow>;
```

The implementation mirrors htslib byte-for-byte, modulo idiomatic
shape:

- Scratch buffers (`f`, `b`, `s`, `qual2prob` lookup) live on a
  `BaqEngine` (Layer 2) and are passed in via a borrow so a
  per-read call does not allocate.
- `qual2prob` is the lookup table htslib builds at
  [`probaln.c:46`](../../htslib/probaln.c#L46) (`10^(-q/10)`).
  Build it once in `BaqConfig::new` and reuse.
- The overflow checks in
  [`probaln.c:96-104`](../../htslib/probaln.c#L96-L104) become a
  `BaqOverflow` error variant. htslib returns `INT_MIN`; we return
  a structured error so the driver can count it.
- Floating-point: htslib uses `double` for the DP tables (`f`,
  `b`, `s`, `m`) and `float` for the per-base quality lookup
  (`qual[]`, `g_qual2prob[]`) and the input parameters
  (`probaln_par_t.d`, `.e`). The port mirrors that exactly — `f64`
  for DP, `f32` for the lookup and config — to keep the byte-for-byte
  parity test cheap.

### Layer 2 — per-read driver (`baq::BaqEngine::process`)

The driver replaces `realn.c:sam_prob_realn`. It owns the scratch
buffers (HMM tables, translated reference, translated query, state
and quality arrays) so a batch of reads reuses them. Signature:

```rust
pub struct BaqEngine {
    cfg: BaqConfig,
    // scratch — reused across reads
    f: Vec<f32>,
    b: Vec<f32>,
    s: Vec<f32>,
    state: Vec<i32>,
    q: Vec<u8>,
    tref: Vec<u8>,
    tquery: Vec<u8>,
}

pub enum BaqOutcome {
    Capped(PreparedRead),         // success: bq_baq = min(BQ, BAQ)
    Skipped(BaqSkipReason),       // count and drop; matches FilterBucket::BaqRejected
}

pub enum BaqSkipReason {
    Unmapped,                     // BAM_FUNMAP — should not reach BAQ given upstream filter, but guard
    EmptyQuery,
    QualAbsent,                   // qual[0] == 0xff (BAM "no quality" sentinel)
    NoMatchInCigar,               // realn.c:192 — entirely soft-clipped / indels-only
    ContainsRefSkip,              // BAM_CREF_SKIP — RNA-seq, out of scope
    HmmOverflow,                  // probaln_glocal returned INT_MIN
    RefWindowPastChromEnd,        // ref slice ran out before xe — clamp + downgrade rather than crash
}

impl BaqEngine {
    pub fn new(cfg: BaqConfig) -> Self;
    pub fn process(
        &mut self,
        read: &MappedRead,
        ref_fetcher: &dyn RefSeqFetcher,
    ) -> BaqOutcome;
}
```

The driver's responsibilities, in order:

1. **Walk the CIGAR** to find `(xb, xe, yb, ye)` — first and last
   reference/query offsets covered by `M/=/X` ops. Reject reads
   with `BAM_CREF_SKIP` (RNA-seq) and reads with no match positions
   at all. Mirrors
   [`realn.c:178-193`](../../htslib/realn.c#L178-L193).
2. **Pick the band width.** Default `bw = 7`; widen to
   `abs((xe-xb)-(ye-yb)) + 3` if a single indel forces it
   (`realn.c:195-197`).
3. **Extend the reference window** by `bw/2` beyond the read's
   mapped span on each side, plus the soft-clipped length on the
   relevant side
   ([realn.c:200-203](../../htslib/realn.c#L200-L203)). Clamp at
   reference boundaries; if the right edge would run past chrom
   end, accept the clipped window — htslib does the same at
   [realn.c:231-234](../../htslib/realn.c#L231-L234).
4. **Translate** ACGT(N) → 0/1/2/3/4 for both ref and query into
   the scratch buffers.
5. **Floor** input qualities to `min_input_bq = 4`. htslib does
   not do this explicitly (Illumina rarely emits below Q4); GATK
   does ([BAQ.java:76-78](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/baq/BAQ.java#L76-L78)).
   Flooring costs ~one branch per base and prevents a Q0/Q1
   numerical edge in the HMM. Keep the floor.
6. **Call `probaln_glocal`.** On overflow, return
   `Skipped(HmmOverflow)`.
7. **Cap.** Walk the CIGAR again. For each `M/=/X` base at query
   index `i` covering reference offset `x`:
   - if `state[i] & 3 != 0` (HMM says this base is an
     insertion, not a match) or `state[i] >> 2 != x - xb + (i - y)`
     (HMM aligned this base to a different reference position),
     set `bq_baq[i] = 0`;
   - otherwise `bq_baq[i] = min(qual[i], q[i])`.
   Soft-clipped and insertion-emitted positions keep their original
   BQ (htslib's `apply_baq` only touches match positions —
   [realn.c:244-265](../../htslib/realn.c#L244-L265)). Deletion
   positions carry no read base.

Output: a `PreparedRead` with `seq`, `cigar`, `mapq`, `flag`,
`adaptor_boundary`, `mate_*` carried through from `MappedRead`,
and `bq_baq` from step 7.

### Why a port and not an FFI binding

We have htslib in the tree
([`/home/jose/devel/join_per_sample_vcfs/htslib/`](../../htslib/));
binding to `probaln_glocal` via `bindgen` would be ~50 lines and
would avoid re-implementing the HMM. We do not do this:

- The whole project is Rust-without-htslib (noodles + our own
  walker). Pulling in htslib for one function adds a C
  dependency, breaks the static-build story, and forces every
  CI image to carry zlib/libdeflate.
- The function is small (~150 lines of dense C). The risk of a
  port bug is bounded by the parity-test corpus (see below).
- The parity test gives us a permanent, run-on-CI guarantee that
  the port matches htslib byte-for-byte. That is stronger than
  the looser "we depend on the upstream version we shipped"
  guarantee an FFI binding gives.

## Cross-checking the port

**No FFI.** Pure-Rust tests, fed by pre-generated fixtures lifted
directly from htslib's own realn test suite.

htslib ships three realn fixture triples under
[`htslib/test/`](../../htslib/test/):

- `realn01.sam` + `realn01.fa` + `realn01_exp.sam`
- `realn02.sam` + `realn02.fa` + `realn02_exp.sam` (plus `-a` and
  `-e` variants for the `BAQ_APPLY` and `BAQ_EXTEND` modes)
- `realn03.sam` + `realn03.fa` + `realn03_exp.sam`

The `_exp.sam` files are what `sam_prob_realn(..., BAQ_APPLY)`
writes: the input record plus a `BQ:Z:` tag where each character
is `qual[i] - min(qual[i], q[i]) + 64`
([htslib/realn.c:266](../../htslib/realn.c#L266)). Decoding it back
to the capped quality is one subtract per base:

```
bq_baq[i] = qual[i] - (BQ[i] - 64)
```

These three fixtures already exercise the edges we care about
(per the `@CO` comments in each input SAM): reads that overhang
the reference, soft-clipped ends, all-insert / all-soft-clip
degenerate reads, and reads with deletions. That covers most of
the categories the previous draft enumerated.

Fixture handling, in this repo:

- Copy the relevant `realn0N.{sam,fa}` and `realn0N_exp.sam`
  files into [`tests/data/baq/`](../../tests/data/baq/).
- Add a `README.md` next to them recording the upstream path
  (`htslib/test/realn0N.*`) and the htslib commit hash they were
  taken from, so a future engineer can refresh them.
- The test parses the input SAM with `noodles-sam` and the FASTA
  with `noodles-fasta` (both already on the manifest), builds
  `MappedRead`s, calls `BaqEngine::process`, and asserts
  `bq_baq[i] == qual[i] - (BQ[i] - 64)` for every match-position
  base in every fixture read. The parity bar is exact byte
  equality on every match-position base.

If a category turns out to be under-covered by the three upstream
fixtures (e.g. a long homopolymer is something we want explicit
coverage of and is not in `realn01-03`), we can extend the corpus
by running `samtools calmd -Br` on a small synthetic CRAM
ourselves and committing the result alongside, following the same
SAM-plus-FASTA pattern. The README documents the regeneration
recipe so the addition is reproducible.

Reasons for fixtures-not-FFI:

1. No C build in the Rust test suite — `cargo test` works on
   anyone's machine without htslib build deps.
2. The fixtures are the actual ground truth users' pipelines have
   been calibrated against (`samtools calmd -Br` output). An FFI
   test would only prove byte-parity with htslib's *current*
   version, not with that historical artefact.
3. They already exist upstream, so commit 1 is mostly file copies
   plus a parser, not corpus generation.

## Migration plan — four commits, each shippable

### Commit 1 — `probaln_glocal` port + parity fixture

- Add the `baq` module under
  [`src/per_sample_caller/`](../../src/per_sample_caller/):
  `probaln.rs` (Layer 1), `mod.rs`, `errors.rs` (`BaqOverflow`,
  `BaqSkipReason`).
- Copy htslib's `realn01-03.{sam,fa}` + `realn0N_exp.sam` into
  [`tests/data/baq/`](../../tests/data/baq/) and add a `README.md`
  there recording the upstream path and commit hash (per
  §"Cross-checking the port" above).
- Add the parity test in `baq/tests.rs` (alongside the module).
- **No integration.** The port is a private island; the walker is
  unaffected.
- Pass criterion: every fixture read's per-base `q` matches
  samtools' `BQ` tag exactly.

### Commit 2 — `BaqEngine::process` + filter-bucket counter

- Add Layer 2 (`BaqEngine`, `BaqOutcome`, `BaqSkipReason`).
- Add `FilterBucket::BaqRejected` to
  [`cram_input.rs`'s filter-count machinery](../../src/per_sample_caller/cram_input.rs)
  (search for `FilterBucket` and add the variant + counter slot).
- Unit-test the driver on synthetic `MappedRead`s for each
  `BaqSkipReason` plus a happy path.
- **Still no integration.** The walker is unaffected.

### Commit 3 — pipeline integration

- At the seam where `MappedRead` becomes `PreparedRead` today,
  insert a rayon-parallel BAQ pass over a batch of reads. The
  batch size is bounded by `MAX_RECORD_SPAN` × concurrency
  (already memory-bounded by the walker's own filter); a starting
  point is `chunk_size = 1024` reads per rayon task, tuned by
  the bench below.
- Each thread owns its own `BaqEngine` so scratch buffers do not
  contend.
- Wire `Skipped(_)` outcomes into the per-sample run summary
  alongside the existing skip categories.
- Run the full pileup test suite. The walker's tests are
  byte-identical because they construct `PreparedRead`s directly
  with explicit `bq_baq` values — they do not go through BAQ.
- The end-to-end test (CRAM → `.psf`) will produce *different*
  output than today, because today's pipeline passes raw BQ
  through. The committed `.psf` golden file for the existing
  end-to-end test needs to be regenerated; commit the new file
  in the same commit as the integration.

### Commit 4 — benchmark + record

- Add a `benches/baq.rs` Criterion bench: per-read BAQ throughput
  and BAQ-as-fraction-of-end-to-end-pipeline time.
- Record numbers under
  [`ia/reports/implementations/`](../reports/implementations/)
  with a short note linking back to this plan.
- Target: BAQ should not dominate the per-sample wall clock. If
  it does (>30% of total CPU time after rayon), tune the
  per-task chunk size; do not pre-optimise.

## Test plan

- **Parity (commit 1).** Every fixture read's per-base posterior
  matches samtools `BQ`. Exact byte equality.
- **Edge cases (commit 2).** One unit test per `BaqSkipReason`
  variant. Plus a happy-path test: a synthetic read with one
  mismatch in a clean reference context — BAQ should not lower
  its BQ much; another with a 3 bp deletion mid-read — BAQ
  should lower BQ on the flanks of the deletion. Assert
  qualitative behaviour, not specific values (those are pinned
  by the parity fixture).
- **Walker invariance (commit 3).** The existing pileup test
  suite at
  [`pileup/tests.rs`](../../src/per_sample_caller/pileup/tests.rs)
  must pass unchanged. It does not exercise BAQ — it constructs
  `PreparedRead`s with explicit `bq_baq` — so any test breakage
  is a sign of a regression in the walker, not BAQ.
- **End-to-end (commit 3).** Re-record the existing CRAM → `.psf`
  golden output with BAQ on. Diff against the pre-BAQ output and
  spot-check: BAQ should drop BQ around indels and keep BQ
  elsewhere. Commit the new golden file.
- **Counter (commit 3).** The per-sample run summary's
  `BaqRejected` count is reported and is nonzero on a constructed
  CRAM containing a read that triggers `HmmOverflow` (e.g. by
  having an absurd CIGAR — easier to construct than to find
  organically).

## Risks and mitigations

- **Numerical determinism across architectures.** The HMM uses
  `f32`. Different CPUs / compilers can produce slightly different
  results from the same float sequence (FMA, reassociation). htslib
  has historically been stable here; if the parity test starts
  failing on a non-x86 CI runner, the mitigation is to add a
  `[cfg(target_arch)]`-gated tolerance to the parity test (off by
  default) rather than switch to `f64` (which would lose byte parity
  with samtools entirely).

- **Reference fetcher granularity.** The ref window extension can
  reach past `MAX_RECORD_SPAN` by `bw/2 + soft_clip_len`. The
  reference fetcher must serve those bytes. The existing fetcher
  ([ref_fetcher.rs](../../src/per_sample_caller/ref_fetcher.rs))
  already serves arbitrary ranges via the FAI index; no new
  machinery is needed, but the BAQ driver must clamp at chrom
  bounds (step 3 above) rather than asking the fetcher for bytes
  past the end.

- **Long-read regression.** The htslib auto-select at `l_qseq >
  1000` swaps `d` from `1e-3` to `1e-7` — a very different gap
  prior. The architecture spec does not commit to long-read
  behaviour and the walker filters reads over `MAX_RECORD_SPAN`
  (default 5000) upstream, so this branch only fires on reads
  the walker will accept. Worth flagging in the bench results
  if long-read CRAMs are ever a target — the long-read branch
  is not exercised by the short-read fixture and so is not
  covered by the parity test.

- **`MappedRead`-to-`PreparedRead` field carrying.** Every field
  the walker reads off `PreparedRead` must be populated by the
  driver, not left as `Default`. The field list at
  [mod.rs:215-274](../../src/per_sample_caller/pileup/mod.rs#L215-L274)
  is the source of truth. Risk is a field added to `PreparedRead`
  later and not wired through the driver; mitigation is the
  walker's existing
  [`PreparedRead::length`](../../src/per_sample_caller/pileup/mod.rs#L288)
  invariant check (`seq.len() == bq_baq.len() == ...`) — it would
  fire immediately on a populate-by-default field that breaks the
  invariant.

- **Architecture-spec correction.** The spec's claim that the two
  reference implementations agree on parameters
  ([architecture spec §289-295](../specs/calling_pipeline_architecture.md#L289-L295))
  is incorrect for the gap-open probability — htslib uses 1e-3,
  GATK uses 1e-4. The plan picks 1e-3 (see "Decisions to make"
  above) and the spec paragraph should be corrected in a small
  follow-up edit, ideally in commit 1 of this plan so the spec and
  the code agree at every revision.

## Out-of-scope follow-ups

- **Extended BAQ** (`BAQ_EXTEND` in htslib's nomenclature). Mode
  where BAQ values are not capped at `min(BQ, BAQ)` but instead
  exported as a separate quality with the BAQ-derived "spread"
  applied. Only worth implementing if downstream consumers
  explicitly want it; the freebayes-style likelihood Stage 1
  feeds wants capped BQ.

- **Long-read parameter tuning.** htslib's `d = 1e-7` for
  `l_qseq > 1000` is a reasonable starting point but is not
  parity-tested by the short-read fixture. When/if long-read
  CRAMs become a real target, build a long-read fixture and
  cross-check.

- **Trusting an upstream `BQ` tag.** If a future workflow runs
  `samtools calmd` upstream and wants us to honour the cached
  result, add a `BaqConfig::trust_existing_bq_tag = true` path
  that reads `MappedRead.bq_tag` (which would need adding) and
  short-circuits the HMM. Not motivated today.

- **Cross-validation against GATK BAQ.java.** A nice-to-have
  third-party check that we do not statistically diverge from
  GATK at parameter `d = 1e-3` (since GATK's published default
  is `d = 1e-4`, we expect a small systematic difference — but
  it is worth quantifying once on a real-data corpus).
