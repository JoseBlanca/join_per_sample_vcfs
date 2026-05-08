# Algorithm Review: pileup walker vs. freebayes

**Date:** 2026-05-08
**Reviewer:** Claude (algorithm-comparison study)
**Module reviewed:** `src/per_sample_caller/pileup`
**Reference codebase:** `freebayes/` (in-tree copy; the
`AlleleParser` + `Allele` + `LeftAlign` layer only — the Bayesian
model layer is downstream of the pileup-equivalent stage and out
of scope)
**Status:** Advisory — open backlog to triage. Companion to
[pileup_samtools_comparison_2026-05-07.md](pileup_samtools_comparison_2026-05-07.md);
freebayes is the *closer* relative because the calling-pipeline
spec already inherits its scalar set, indel-anchor convention,
indel BQ-proxy window, and first/last-CIGAR-op indel-rejection
rule. This study looks for further lessons in the
read-registration / decomposition layer specifically.

---

## 1. Scope

This is **not** a defect-finding review. It is a comparative
study of freebayes' read-registration / pileup-equivalent layer
against our Stage 1 walker, looking for techniques and corner-case
handling worth importing or — equally usefully — for places where
our model is already cleaner than freebayes' and we should not
accidentally regress.

In-scope project files:

- [src/per_sample_caller/pileup/mod.rs](../../src/per_sample_caller/pileup/mod.rs)
- [src/per_sample_caller/pileup/walker.rs](../../src/per_sample_caller/pileup/walker.rs)
- [src/per_sample_caller/pileup/active_set.rs](../../src/per_sample_caller/pileup/active_set.rs)
- [src/per_sample_caller/pileup/decompose.rs](../../src/per_sample_caller/pileup/decompose.rs)
- [src/per_sample_caller/pileup/cigar_cursor.rs](../../src/per_sample_caller/pileup/cigar_cursor.rs)
- [src/per_sample_caller/pileup/slot_allocator.rs](../../src/per_sample_caller/pileup/slot_allocator.rs)
- [src/per_sample_caller/pileup/open_record.rs](../../src/per_sample_caller/pileup/open_record.rs)
- [src/per_sample_caller/cram_input.rs](../../src/per_sample_caller/cram_input.rs)

Freebayes files studied:

- `freebayes/src/AlleleParser.{h,cpp}` — the read-registration
  driver and per-position iterator.
- `freebayes/src/Allele.{h,cpp}` — allele types, complex-merge
  helpers (`mergeAllele`, `subtractFromEnd/Start`,
  `addToStart/End`).
- `freebayes/src/LeftAlign.{h,cpp}` — `stablyLeftAlign`, applied
  during BAM ingest before allele registration.
- `freebayes/src/Parameters.{h,cpp}` — the user-visible knobs
  this layer exposes (mismatch budgets, coverage cap, BQ floor,
  `--use-min-indel-quality`).
- `freebayes/src/IndelAllele.{h,cpp}`, `Sample.{h,cpp}` — typed
  helpers used by the registrar.

Out of scope:

- The freebayes Bayesian model layer (`Genotype.cpp`,
  `DataLikelihood.cpp`, `Marginals.cpp`, `Result*.cpp`,
  `Bias.cpp`, `Dirichlet.cpp`, `Ewens.cpp`, `Multinomial.cpp`,
  `Contamination.cpp`). The `.psf` scalar set was already
  derived from this layer in
  [freebayes_posterior_gt_probs.md](../specs/freebayes_posterior_gt_probs.md);
  no new evidence is sought here.
- Stage 5 / Stage 6 — the merger and posterior engine.
- Haplotype-window construction (`buildHaplotypeAlleles`,
  `getCompleteObservationsOfHaplotype`,
  `getPartialObservationsOfHaplotype`). That is freebayes'
  *evaluation*-time machinery; we have no analogue and don't
  need one (Stage 4's grouping plus Stage 5's per-group
  reconstruction does the same job from `.psf` data).

## 2. Summary verdict

**Freebayes is the algorithmic ancestor of our Stage 1 design,
not the rival.** Our spec has already adopted the freebayes
scalar set, the `l + 2` indel-BQ window, the first/last-CIGAR-op
indel-rejection rule, and the haplotype-allele
"REF-extension on overlap" convention. The five findings below
are mostly small refinements: places where freebayes implements
something our spec is silent on (left-alignment, per-read
mismatch budgets, BQ floor for mismatch counting) and places
where the spec should explicitly *contrast* with freebayes so
that future readers understand the deliberate divergence
(useMinIndelQuality choice, soft-clip and N-base handling).

Two genuine wins for our design — places we beat freebayes and
should not regress — are noted up front in §3.

## 3. Where we already match or beat freebayes

These are noted up-front so we don't accidentally regress them
while acting on the findings below.

- **`A1` — Mate-overlap dedup at all.**
  Freebayes performs **no** mate-overlap deduplication. A grep
  through `freebayes/src/AlleleParser.cpp` for "overlap" yields
  only VCF header text and haplotype-window terminology — no
  overlap-resolution code. `RegisteredAlignment` carries
  `isPaired` / `isMateMapped` for tagging, but both mates of an
  overlapping fragment contribute their bases independently to
  the per-position counts. Combined with freebayes' use of MQ
  for the reference allele's quality, this means a double-counted
  ref/alt evidence vote at every overlap position. We do
  per-position mate-overlap with samtools-style summed-BQ math
  for agree pairs and 0.8-haircut math for disagree pairs (see
  [pileup_samtools_comparison_2026-05-07.md `S7`](pileup_samtools_comparison_2026-05-07.md#s7--adopt-samtools-bq-combining-math-for-match-only-mate-overlap)).
  Strictly better than freebayes here; do not pull this forward.

- **`A2` — Decomposition rejects indels at *both* ends of the read.**
  Freebayes rejects **deletions** at the first or last CIGAR op
  (`AlleleParser.cpp:1681-1682`) but does *not* apply the same
  guard to insertions (`AlleleParser.cpp:1706-1770` has no
  equivalent check). We reject *both* insertions and deletions
  at the first or last CIGAR op
  ([decompose.rs:106](../../src/per_sample_caller/pileup/decompose.rs#L106),
  [decompose.rs:120](../../src/per_sample_caller/pileup/decompose.rs#L120)),
  on the same flanking-evidence rationale that freebayes spells
  out in its deletion comment (`AlleleParser.cpp:1677-1679`).
  Insertions at read tails are just as untrustworthy — they
  almost always indicate adapter readthrough or local
  misalignment — so the symmetric rule is the right one.
  Regression-tested by
  `first_op_indel_is_dropped_but_rest_of_read_keeps_events` and
  `soft_clip_then_indel_at_alignment_start_one_drops_indel` in
  [decompose.rs](../../src/per_sample_caller/pileup/decompose.rs).

- **`A3` — Streaming record emission with an explicit close boundary.**
  Freebayes keeps `registeredAlignments` as `map<long unsigned
  int, deque<RegisteredAlignment>>` keyed by alignment *end*
  position (`AlleleParser.h:196`,
  `AlleleParser.cpp:2095`). Alignments are removed lazily when
  the position passes them
  (`AlleleParser.cpp:2127-2163` `removeRegisteredAlignmentsOverlappingPosition`,
  invoked only as a side-effect of the coverage-cap path; the
  steady-state cleanup is `updateRegisteredAlleles` at line
  2172, which only nulls *alleles* whose
  `position + referenceLength < currentPosition`). There is no
  per-record "this record is now safe to close and emit" step —
  freebayes' output unit is a per-position genotype call, not a
  per-record observation summary. Our spec
  ([calling_pipeline_architecture.md §"Open-record closure and
  the read-span filter"](../specs/calling_pipeline_architecture.md))
  commits to the explicit per-record closure rule and the
  `MAX_RECORD_SPAN` upstream filter that bounds it. Architectural
  divergence we want to keep — see §6.

  *Empirical addendum (2026-05-08).* The Scope C microbenchmark
  ([pileup_freebayes_bench_c_2026-05-08.md](../reports/implementations/pileup_freebayes_bench_c_2026-05-08.md))
  measures the bookkeeping cost of the two data structures in
  isolation (50 kb walker span, no fold work, no I/O, parity-checked
  event stream). At an Illumina-shaped distribution (90% `span=1`,
  10% `span ∈ [2, 10]`) our `BTreeMap<u32, OpenSlot>` with
  merge-on-overlap is **2.05–2.20× faster** than a freebayes-shape
  flat allele pool with linear-sweep ageing across the full
  10–1000× coverage sweep. At a deletion-heavy stress
  distribution (50% `span ∈ [10, 100]`) the gap reaches **17×**
  because each long-footprint event keeps a pool entry alive for
  10–100 walker steps, turning the linear retain into `O(N²)`,
  while the merge-on-overlap rule keeps the BTreeMap small. The
  architectural argument in §6 was correct on output-contract
  grounds; the bench shows the design choice is also
  performance-positive, not just performance-neutral.

- **`A4` — Per-column depth handling preserves allele-frequency expectation.**
  Freebayes' `--limit-coverage` cap at `AlleleParser.cpp:2075-2092`
  is a **full-position purge**: when the per-position counter
  exceeds `skipCoverage`, every alignment overlapping the
  position is dropped (`removeRegisteredAlignmentsOverlappingPosition`)
  and the position is permanently marked in
  `coverageSkippedPositions`, so that *future* reads touching it
  also get skipped. This is much more aggressive than samtools'
  per-column cap and biases allele frequencies in mid-coverage
  regions where the cap fires intermittently. Our `S5` work
  ([pileup_samtools_comparison_2026-05-07.md `S5`](pileup_samtools_comparison_2026-05-07.md#s5--per-column-depth-cap-originally-framed-as-per-allele))
  picked the samtools-style "truncate the contributor list to
  the first N items at fold time" pattern precisely because it
  preserves AF in expectation. Don't regress to freebayes' purge
  semantics.

## 4. Findings — backlog to act on

Each finding has a stable id (`Fi` for "freebayes-inspired") so
we can reference them in commits and a one-by-one work plan.

Priority key (same as the samtools review):
- **High** — clear user-visible benefit, low risk.
- **Medium** — usability or hygiene; nice to have.
- **Low** — only matters under conditions we don't currently
  target, or speculative.

### `F1` — Per-read mismatch-fraction filter (default on, single knob)

- **Priority:** Medium
- **Effort:** Small (one config field, one constant, one
  `FilterCounts` field, mismatch counter folded into the
  existing per-read CIGAR walk)
- **Status:** Designed and approved 2026-05-08; implementation
  to follow. Subsumes the original `F2`.

**Observation.** Freebayes drops an alignment **entirely** if
any of four per-read mismatch budgets is exceeded
(`AlleleParser.cpp:2102-2109`):

```cpp
if (ra.alleles.empty()
    || ((float) ra.mismatches / (float) currentAlignment.SEQLEN) > parameters.readMaxMismatchFraction
    || ra.mismatches > parameters.RMU
    || ra.snpCount > parameters.readSnpLimit
    || ra.indelCount > parameters.readIndelLimit) {
    rq.pop_front(); // backtrack
}
```

The `mismatches` counter only increments for mismatches with
`bq >= BQL2` (default 10, `Parameters.cpp`), so low-quality
bases don't poison the budget. All four knobs are off by
default in freebayes — but the `--read-max-mismatch-fraction`
formula is the most useful of the four because it captures the
shape of the failure modes that pass MAPQ + BAQ:

- **Adapter readthrough at the 3' end** — a mismatched tail
  that BAQ can flatten BQ-wise but can't move into a soft clip.
- **Wrong-reference / contamination** — bacterial / organelle
  / PhiX reads aligning at low identity to "best available"
  placement under a uniqueness-permissive aligner.
- **Chimeric reads** where one end has a long mis-located tail.

We have no equivalent in the
[CramMergedReaderConfig](../../src/per_sample_caller/cram_input.rs)
filter cascade today. A read with 30 mismatches in 100 bp
passes our filters as long as MAPQ and length clear; every
mismatch becomes a low-quality SNP allele in the per-position
records and inflates `prodQout` against true alts in
low-coverage cohorts.

**Decision (2026-05-08).** Adopt a **default-on** filter with
**two user-facing knobs** rather than freebayes' four-knob
opt-in set:

- **Config:** two new fields on `CramMergedReaderConfig`:
  - `max_read_mismatch_fraction: Option<f32>` — the threshold.
    `None` = filter disabled; `Some(x)` = drop a read whose
    mismatch fraction exceeds `x`. Default `Some(0.10)`.
  - `mismatch_bq_floor: u8` — only mismatches whose base
    quality clears this floor count toward the fraction. The
    floor exists so the filter doesn't fire on genuinely
    low-quality bases (which the per-base BQ already
    de-emphasises in the likelihood). Default `10`, matching
    freebayes' `BQL2`. `0` disables the floor (every mismatch
    counts).
- **Why two knobs.** The floor is exposed as a calibration
  knob rather than buried as a private constant: it's
  decision-relevant for users with unusually clean or
  unusually noisy lanes, and a future-engineer reading the
  config can see it sitting alongside the threshold. The two
  knobs are semantically distinct (threshold = "how much is
  too much"; floor = "what counts as a real mismatch") and
  combining them into a single value would obscure both.
  Subsumes the originally-proposed `F2`.
- **Default rationale.** Real Illumina at MAPQ 20+ sits at
  ~0.5-1% mismatch rate; 10% is comfortably outside that
  distribution and inside the regime where adapter-runthrough,
  contamination, and chimeric tails live. The BQ floor of 10
  is freebayes' empirically-calibrated `BQL2`; below that
  threshold the base call itself is too noisy to trust as
  evidence of misalignment.
- **Numerator:** count of `M`-op mismatches whose raw base
  quality clears `mismatch_bq_floor`. Note: at the
  `cram_input` stage we have raw BQ, not BAQ-adjusted BQ —
  BAQ runs in a later stage. Filtering on raw BQ is the right
  level for this filter: we're rejecting whole reads, not
  per-base evidence, and we want to do so before paying the
  BAQ cost.
- **Denominator:** the count of `M`-op bases that are
  ATGC-on-both-sides (skipping any position where either the
  read or the reference has an `N`). This differs from
  freebayes' `SEQLEN` denominator, which dilutes the fraction
  with soft-clipped bases. Ours answers "what fraction of the
  *confidently aligned* bases are mismatching?" and is more
  aggressive on heavily-clipped reads — the right direction
  for contamination defence.
- **Where:** drop happens during per-read CIGAR analysis in
  `cram_input` (alongside the existing flag / MAPQ / length
  filters). The read never reaches the walker.
- **Visibility:** new `FilterCounts.high_mismatch_fraction`
  field, surfaced in the run summary alongside the existing
  drop categories.

**Rationale for default-on.** Three reasons (per the design
discussion behind this revision):

1. Filters that are off-by-default don't get used. The whole
   reason `drop_qc_fail` and `drop_duplicate` are on-by-default
   is that nobody would discover them otherwise.
2. Reproducibility — defaults define "what the tool does." A
   safety net that sits in the codebase but never fires is
   worse than no safety net (it implies coverage it doesn't
   deliver).
3. Inverts the burden — instead of "opt-in to safety," it's
   "opt-out for permissive mode." Right direction for a
   defensive filter.

**Rationale for one knob, not four.** Freebayes' separate
`readSnpLimit` / `readIndelLimit` / `RMU` thresholds duplicate
the same signal at different angles; in practice
`max_read_mismatch_fraction` alone catches the failure modes
above. Fewer knobs is easier for users to understand and
harder to misconfigure. If real data later shows a need for an
indel-specific budget (e.g. PCR-stutter regions in homopolymer
contexts), a follow-up adds it then with measured justification.

**Note on roles.** Strictly speaking, dropping high-mismatch
reads is "the aligner's job" — a permissive aligner shouldn't
emit them. But in practice BWA/minimap2/etc. emit alignments
at a wide quality range and rely on downstream filtering;
MAPQ alone doesn't catch high-mismatch but unique placements.
This filter sits at the same `cram_input` BAM-input layer as
the existing MAPQ / flag filters, on the same defence-in-depth
grounds: the layer's job is to enforce a baseline of input
quality the walker can rely on.

**Risk.** Default-on filters can quietly drop more than the
user expects. Mitigations:

1. The threshold (10%) is conservative enough that
   normal Illumina at MAPQ 20+ is essentially never affected.
2. The drop counter is surfaced in the run summary so a
   user looking at run stats can immediately see how many
   reads it fired on. Consistently high counts indicate the
   threshold needs revisiting; consistently zero counts mean
   the filter is paying for itself only in worst-case
   insurance.
3. A user with a known-noisy lane can opt out via
   `max_read_mismatch_fraction = None`.

### `F3` — Stage 1 indel left-alignment for placement canonicalisation

- **Priority:** Medium → High once we have real samples to test
- **Effort:** Medium (port `freebayes/src/LeftAlign.cpp`'s
  `stablyLeftAlign`, run it during `cram_input` BAM-ingest
  before decomposition; ~400 lines)

**Observation.** Freebayes left-aligns indels at registration
time, before allele construction
(`AlleleParser.cpp:2055-2058`):

```cpp
if (parameters.leftAlignIndels) {
    int length = currentAlignment_end_position - currentAlignment.POSITION + 1;
    stablyLeftAlign(currentAlignment,
                    currentSequence.substr(currentSequencePosition(currentAlignment), length));
}
```

Default is **on** (`Parameters.cpp:426`; the off switch is
`-O / --dont-left-align-indels` at line 655). The algorithm
(`LeftAlign.cpp:25-200`) shifts each indel left as far as the
preceding reference bases match the inserted/deleted sequence,
canonicalising placements within homopolymers and short tandem
repeats. The CIGAR is mutated in place, so the downstream
`registerAlignment` decomposition sees the canonical placement.

We don't do this. Our walker takes the BAM's CIGAR verbatim and
emits indel events at whatever anchor the aligner picked. The
consequence is a real correctness issue under our allele-equality
rule:

> Allele equality (so two reads contribute to the same entry) is
> byte-for-byte identity of `seq` under a fixed REF-span.
> ([pileup_walker.md:140](../specs/pileup_walker.md#L140))

If two reads support the same homopolymer-context insertion but
the aligner placed them at different anchors (this is the common
case in a 5-bp poly-A tract — bwa-mem will scatter the anchors
across the run depending on which side of the read the insertion
sits), they end up in **different `PileupRecord`s**. The
per-allele scalars get split across two adjacent records, the
cohort merger sees a 1-vs-1 frequency spread instead of 2-vs-0,
and the variant looks weaker than it is. BAQ does not fix this —
BAQ adjusts BQ, it does not move CIGAR ops.

**Decision (2026-05-08).** Adopt as **always-on, no
configuration knob**:

- Port the `stablyLeftAlign` algorithm to a private function in
  `cram_input` that mutates the CIGAR in place, leaving the read
  sequence and BQ array untouched (the bytes don't change — only
  their CIGAR annotation does). Single function, no module
  boundary.
- Run it during the per-read decode in `cram_input`, immediately
  after `record_buf_to_mapped_read`, **before** the F1
  mismatch-fraction filter (so F1 sees post-canonicalisation
  CIGARs).
- **No CLI / API knob.** Unlike F1's defensive filter (which has
  a real opt-out case for known-noisy lanes), left-alignment is
  pure normalisation: it never throws away data, and the only
  cost of "turning it off" is fragmented allele scalars at
  homopolymer-context indels. The performance gain from skipping
  it is minimal. Adding a knob would mean an API/CLI surface and
  a documentation entry for a switch that no caller has a
  legitimate reason to flip.

**Orientation invariance.** BAM/CRAM stores `seq` and `qual` in
forward-reference orientation (a reverse-mapped read has its
sequence reverse-complemented before storage), and CIGAR is in
the same forward-reference order. The left-alignment algorithm
operates on these forward-reference inputs and never branches on
`is_reverse_strand`. Two reads of the same biological indel on
opposite strands therefore *converge* on the same canonical
anchor, which is the point of the filter — accidentally
mirroring the shift direction by strand would defeat it. Tests
explicitly cover the forward-vs-reverse-strand convergence case.

**Scope of the algorithm we implement.** Only single-indel
left-alignment within the immediately preceding `M`-class op.
Specifically:

- For each indel in the CIGAR, walk leftward one base at a time
  while the shift condition holds (insertion: rotated last
  inserted base equals the base just before the current anchor;
  deletion: base just before the deletion equals the last
  deleted base) AND the preceding `M`-class op has length to
  give. Each base of shift takes one base from the preceding
  `M` op and gives it to the following `M`-class op.
- Bounded by the read's reference span — an indel cannot shift
  past the read's `alignment_start`. If the homopolymer extends
  to the read's left edge, the shift terminates with the indel
  at first-CIGAR-op position; the walker's existing
  first/last-op-indel rejection rule
  ([decompose.rs:106](../../src/per_sample_caller/pileup/decompose.rs#L106),
  [decompose.rs:120](../../src/per_sample_caller/pileup/decompose.rs#L120))
  drops it on the next stage. That's the right outcome —
  flanking evidence on only one side is the same condition the
  rule already protects against.
- **Not implemented:** freebayes' "merge neighbouring indels"
  pass after left-alignment. That handles the rare case of two
  close indels canonicalising into adjacency. Real Illumina
  data rarely produces those patterns at variant-call quality;
  if they show up in real samples, we add the merging pass
  then. Adjacent-indel CIGARs are left untouched.

**Rationale.** Restores the convergence property our
allele-equality rule needs to be lossless. Without it, the
"per-allele scalars are sufficient" property breaks down at
exactly the regions where short-read alignment is most ambiguous
— which are also exactly the regions Stage 3's DUST filter
already drops. The two filters compose: DUST drops
non-canonicalisable repetitive regions; left-alignment fixes the
canonicalisable ones.

**Alternative considered and rejected.** *"Live with it; trust
the aligner."* This is freebayes' Option A
([calling_pipeline_architecture.md §"Per-read likelihood
quality: BAQ vs. PairHMM vs. trust-the-aligner"](../specs/calling_pipeline_architecture.md))
that we already rejected for SNP confidence. The same
argument — that aligner artefacts dominate at low coverage —
applies to indel anchor choice. Picking BAQ for SNPs but
then trusting the aligner's indel anchors is internally
inconsistent.

**Risk.** Low. The algorithm is small (~80 lines), the inputs
and outputs are well-defined (CIGAR-in, CIGAR-out, read bytes
unchanged), and the tests pin both the homopolymer-anchor-scatter
case and the forward-vs-reverse-strand convergence
explicitly. The "merge neighbouring indels" gap is acknowledged
and a follow-up if real data demands it.

### `F4` — Document the `--useMinIndelQuality` choice and contrast with freebayes' default

- **Priority:** Low (documentation only)
- **Effort:** Tiny (a paragraph in
  [calling_pipeline_architecture.md §"Indel BQ proxy"](../specs/calling_pipeline_architecture.md))

**Observation.** Our spec says:

> for a supporting read we take the window of `l + 2` BAQ-adjusted
> base qualities centred on the indel ... and use the **minimum**
> BQ in that window as the per-read contribution to the indel
> allele's `Σ max(ln_BQ, ln_MQ)` scalar. This matches freebayes'
> `--useMinIndelQuality` mode.
> ([calling_pipeline_architecture.md:170-179](../specs/calling_pipeline_architecture.md))

Freebayes' **default** is **not** `useMinIndelQuality`. The
default branch (`AlleleParser.cpp:1659-1673` for deletions,
`1735-1749` for insertions) does:

```cpp
qual  = sumQuality(qualstr);
qual += ln2phred(log((long double) L / (long double) l));
qual /= harmonicSum(l);
```

i.e. summed in-window BQ, with a `log(L/l)` scaling correction
and a `harmonicSum(l)` denominator that penalises long indels
(`harmonicSum(1) = 1`, `harmonicSum(5) ≈ 2.28`). The min-mode
branch fires only with `--use-min-indel-quality`
(`Parameters.cpp:644`).

We have a deliberate reason to pick min over freebayes'
default — the min combines cleanly into our `Σ max(ln_BQ, ln_MQ)`
scalar without a length-dependent reweight, and it matches
bcftools' default indel quality treatment (which is the
calibration target for the BAQ-based pipeline). But the spec
currently leaves this as an undefended preference; a future
engineer reading "matches freebayes' `--useMinIndelQuality`
mode" will reasonably ask "why didn't we pick the default
mode then?" without finding an answer.

**Proposal.** Add a one-paragraph note to
[calling_pipeline_architecture.md §"Indel BQ proxy"](../specs/calling_pipeline_architecture.md)
spelling out:

1. Freebayes' default is summed-with-harmonic, not min;
2. We deliberately pick min because our per-allele scalar is a
   sum of `max(ln_BQ, ln_MQ)` and a min-aggregation slots
   into that sum without a per-allele length-dependent
   correction (the harmonic denominator would have to be
   carried as a per-allele field to undo at merge time);
3. bcftools' indel quality is min-based, so the calibration is
   well-studied;
4. If a future calibration study shows the harmonic
   correction is doing real work, the choice can be revisited
   — it is not committed to in code, only in the scalar sum.

**Rationale.** Cheap insurance against future re-litigation.

**Risk.** None — documentation only.

### `F5` — Define behaviour when reference base is `N`

- **Priority:** Low (rare in practice but the rule should be explicit)
- **Effort:** Small (one branch in `CigarCursor::Match` event
  emission; one test)

**Observation.** Freebayes treats `N` in the reference as a
forced mismatch (`AlleleParser.cpp:1465`):

```cpp
if (b != sb || sb == "N") {  // when the reference is N, we should always call a mismatch
```

In our walker the ref base is fetched at fold time
([walker.rs](../../src/per_sample_caller/pileup/walker.rs))
and goes into `alleles[0].seq` verbatim, including `N`.
Allele equality is byte identity, so a read whose `M` op aligns
an `A` to a reference `N` will produce a Match event with
`base = b'A'` that votes against an `alleles[0]` whose REF byte
is `b'N'` — i.e. it would land in a non-REF allele bucket with
`seq = b"A"` while the REF bucket carries `seq = b"N"`. This is
arguably the right thing (we record what the read showed; the
merger can decide), but it's not specified anywhere.

The complementary case — `N` in the *read* — is also unspecified.
Freebayes records `N` read bases as `ALLELE_NULL` and filters
them downstream (`AlleleParser.cpp:1532-1547`). Our walker today
puts `N` into the allele's `seq` like any other base, so an
`N` read base creates an allele with `seq = b"N"`, which goes
through to Stage 5 as a regular allele. Stage 5 will treat it
as a distinct allele from the SNPs at the same position — which
is correct for the scalar accumulation but wastes a bucket on
sequencing-error `N`s.

**Proposal.**

1. Add a paragraph to
   [pileup_walker.md §"Read decomposition"](../specs/pileup_walker.md)
   that documents the current behaviour explicitly: `N` ref
   bases pass through verbatim into REF; `N` read bases pass
   through into their own allele bucket (or into a Match event
   if the ref also says `N`).
2. Decide whether to filter `N` read bases at decomposition
   time (skip the Match event when `seq == N`) — strictly
   preserves the scalar sums for the real alleles at the
   position, at the cost of losing the depth-counter
   contribution. The freebayes precedent is "drop them at
   genotyping time, not registration time", so the conservative
   choice is to keep them and let Stage 5 decide.
3. Add a test pinning either choice so a future refactor cannot
   silently change it.

**Rationale.** Cheap to specify, currently ambiguous; the
question will come up the first time someone runs the tool on
a draft assembly with `N` patches in the reference.

**Risk.** None for option 1 (doc only). Option 2 changes
observable depth at `N` positions — but those positions are
already low-confidence by definition, and Stage 3's DUST
filter drops most of them anyway.

## 5. Findings explicitly *rejected*

For the record, so we don't relitigate them:

- **`R1` — Per-position-saturating coverage purge (`skipCoverage`).**
  Rejected; see §3 `A4` and
  [pileup_samtools_comparison_2026-05-07.md `S5`](pileup_samtools_comparison_2026-05-07.md#s5--per-column-depth-cap-originally-framed-as-per-allele).
  We adopted samtools' per-column truncation, which preserves
  AF in expectation. Freebayes' purge biases AF in
  intermittently-saturating regions and additionally permanently
  poisons the position via `coverageSkippedPositions` — every
  future read that touches the position also gets dropped, even
  if the burst was a single instant of pile-up.

- **`R2` — Whole-chromosome reference loaded eagerly with no eviction.**
  Rejected. `AlleleParser.cpp:673-691` loads the entire
  chromosome into `currentSequence` on chromosome change and
  never evicts the previous one until the next chromosome is
  loaded — which happens to behave like our
  `ChromBoundaryRefFetcher` ([ref_fetcher.rs](../../src/per_sample_caller/ref_fetcher.rs))
  for *single-chromosome* runs but accumulates O(genome) memory
  if a process drives multiple chromosomes through the same
  `Repository` without explicit clears. We picked
  [pileup_samtools_comparison_2026-05-07.md `S6`](pileup_samtools_comparison_2026-05-07.md#s6--chromosome-boundary-eviction-in-the-production-refbasefetcher)'s
  per-chromosome eviction precisely to avoid the latent
  memory blow-up. Don't import freebayes' "load and forget"
  pattern.

- **`R3` — Per-read post-decomposition `clumpAlleles` re-flow.**
  Rejected. `AlleleParser.cpp:1013-1095` runs *after* the
  per-CIGAR-op decomposition and merges adjacent
  non-reference alleles (separated by ≤ `maxComplexGap`
  reference bases of match) into a single
  `ALLELE_COMPLEX`, then redistributes flanking match bases
  between neighbouring alleles to maintain VCF anchoring
  (the `subtractFromEnd` / `addToStart` dance at lines
  1075-1089). The intent — turn `M2 X1 M5 X1 M2` into one
  complex allele covering both X1s with five M bases between
  them — is a deliberate signal-aggregation step ("a
  haplotype is more informative than its parts").

  Our open-record model achieves the same effect at
  *fold* time: when two reads share a non-REF observation
  inside a single open record's footprint, the record's REF
  span is widened and both reads' allele literals are
  rewritten against the wider span
  ([calling_pipeline_architecture.md §"Overlapping events extend
  the anchor REF"](../specs/calling_pipeline_architecture.md)).
  The trigger is different: ours fires when events *overlap*,
  freebayes' fires when events are *near* (within
  `maxComplexGap` reference bases). Two SNPs separated by 5
  matching bases on the same read would be one
  `ALLELE_COMPLEX` in freebayes and two adjacent
  `PileupRecord`s in our walker.

  The samtools study's §6 covered the analogous decision
  for samtools' per-position callback model. The same
  argument applies here: our scalar set + phase chain id
  lets Stage 5 reconstruct the same compound-haplotype
  evidence at merge time, without committing to a single
  `maxComplexGap` window during `.psf` writing. Stage 5
  has the cohort context to decide which
  near-but-not-overlapping events deserve to be co-called;
  Stage 1 does not. Pulling the freebayes pattern forward
  would lock the choice in at a stage that can't make it
  well.

- **`R4` — Soft-clipped bases as `ALLELE_NULL` observations.**
  Rejected. Freebayes records soft clips as
  `ALLELE_NULL` (`AlleleParser.cpp:1776-1798`), filtered out
  by every downstream consumer. They contribute nothing to the
  five scalars (REF or ALT) — the only purpose is to keep
  per-read bookkeeping symmetric with the original read length.
  Our scalar accumulation does not need that symmetry, so we
  drop soft-clipped bases at decomposition time
  ([decompose.rs:137-140](../../src/per_sample_caller/pileup/decompose.rs#L137-L140))
  without recording any event. Strictly less bookkeeping for
  identical scalar output.

- **`R5` — Reference allele quality = mapping quality.**
  Rejected. Freebayes' `ALLELE_REFERENCE` entries take their
  `qual` field from `alignment.MAPPINGQUALITY` rather than
  from the per-base BQs of the matched stretch
  (`AlleleParser.cpp:1484`, `1621`). The base qualities ride
  along as `qualstr` for use by downstream filters but the
  headline `qual` is MQ. This makes sense in freebayes'
  no-BAQ design — without BAQ, base qualities are
  uninformative about *placement* confidence, so MQ is the
  only signal available for the "this stretch is correctly
  reference" claim. With BAQ
  ([calling_pipeline_architecture.md §"Why BAQ in-process"](../specs/calling_pipeline_architecture.md)),
  base qualities *are* informative about local placement
  confidence (BAQ is precisely the per-base
  `P(correctly aligned)`), so feeding them as
  `max(ln_BQ_BAQ, ln_MQ)` per base is strictly more
  information than MQ alone. Don't import.

## 6. Why we keep per-record closure (and not freebayes' open-allele queue)

Freebayes' equivalent of our open-record set is
`registeredAlleles` (a `vector<Allele*>`) plus
`registeredAlignments` (a `map<long unsigned int,
deque<RegisteredAlignment>>` keyed by alignment end position).
Alignments stay alive until the position passes their end
coordinate; alleles stay alive until
`position + referenceLength < currentPosition`. There is no
equivalent of our "the record at position Q closes once the
walker is past `Q + ref_span`" rule. The output unit is a
per-position genotype call, not a per-record observation
summary, so freebayes does not need an explicit close boundary
on its open state — it just lets state age out under per-position
sweeps.

This is the same design tension the samtools review's §6
addressed in a different shape: samtools' callback model also
has no record concept. The argument we made there — that our
records carry compound-haplotype state (the per-allele
`chain_slots`, the haplotype-extended REF span) that a
position-shaped output cannot represent — applies
verbatim against freebayes too. Even though freebayes builds
haplotype alleles internally via `buildHaplotypeAlleles`, that
machinery runs at *evaluation* time on a window of registered
alleles, not at *write* time on a stream. Our `.psf` is the
write-time artefact, and the explicit close boundary is what
lets us stream it without back-patching while preserving
phase-chain reach.

If the closure rule is ever in doubt, the proof of safety is
spelled out in
[calling_pipeline_architecture.md §"Open-record closure and the
read-span filter"](../specs/calling_pipeline_architecture.md):
because reads arrive coordinate-sorted and event anchors sit
at `M`-positions of their reads, a future event with anchor
`S ≥ walker + 1` can merge into a record at `Q` only if
`S < Q + s`, which forces `Q + s > walker`. Once
`Q + s ≤ walker`, no future event can reach back. Freebayes'
"sweep stale alleles each position" model is correct given its
output shape; ours is correct given ours. They are not
substitutable.

## 7. Suggested execution order

Each finding is independent and can be its own small commit.
None of them block each other; pick whichever is most useful
on the current sample backlog.

1. `F4` — document the `--useMinIndelQuality` choice
   (paragraph-only, ~15 minutes; gets the spec consistent
   with itself).
2. `F5` — document `N`-base behaviour and add a pinning
   test (~30 minutes; clears an open question).
3. `F1` — per-read mismatch-fraction filter in
   `CramMergedReaderConfig` (single knob, default-on at 10%,
   drop counter surfaced in run summary; subsumes the old
   `F2` BQ-floor refinement).
4. `F3` — Stage 1 indel left-alignment. The biggest of the
   four and the only one that materially changes behaviour at
   defaults; justifies its own design pass before
   implementation, ideally tied to a homopolymer-context test
   sample so the win is measurable.
