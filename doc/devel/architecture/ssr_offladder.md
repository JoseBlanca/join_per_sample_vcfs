# SSR off-ladder wiring — Stage 1 (`ssr-pileup`)

**Status:** first draft, 2026-06-17, branch `pileup-offladder`. A discussion
starter, in the style of the other SSR architecture docs (proposals **flagged to
argue with**, not conclusions). It designs how off-ladder allele evidence is
generated, scored, aggregated, and written **end-to-end** in Stage 1 — the gap
the Stage-1 remaining-work checklist
([ssr_stage1_remaining.md](../implementation_plans/ssr_stage1_remaining.md) §3)
parked and that Stage 2's candidate-set assembly (spec §5.1) depends on.

Settle this with the user **before** any implementation plan or code.

Cross-refs: spec [ssr_genotyping.md](../specs/ssr_genotyping.md) §4.2/§4.3/§5.1;
architecture [ssr_genotyping_architecture.md](ssr_genotyping_architecture.md)
§5/§6/§8.2/§10; [ssr_pileup.md](ssr_pileup.md) §2 (revision)/§5/§6/§11;
[ssr_shared_types.md](ssr_shared_types.md) §2/§4; the four source files
[`candidate_generation.rs`](../../src/ssr/pileup/candidate_generation.rs),
[`read_analysis.rs`](../../src/ssr/pileup/read_analysis.rs),
[`locus_record.rs`](../../src/ssr/pileup/locus_record.rs),
[`registry_ssr.rs`](../../src/psp/registry_ssr.rs).

---

## 0. The headline: two grounding facts that reshape the spec's design

Before the layer-by-layer design, two facts established from a read of the tree
(2026-06-17) — both *correct the documents this task points at*, so they are
called out first.

### 0.1 The all-CSR revision deleted the storage shape the spec's off-ladder
mirrored

Spec §4.3 designs off-ladder as **six columns** —
`offl_seqs`(dict)/`offl_counts`/`offl_weight` (a confident *tally*, the
off-ladder analog of `hist_*`) **plus** `offl_amb_offsets`/`offl_amb_idx`/
`offl_amb_logliks` (the ambiguous per-read CSR, the analog of `amb_*`). That
mirror only makes sense because the on-ladder side *also* had two regimes:
confident reads → a locus histogram (`hist_*`), ambiguous reads → per-read CSR
(`amb_*`).

The **realign-everything + all-CSR revision (2026-06-15)** removed the histogram
entirely ([locus_record.rs](../../src/ssr/pileup/locus_record.rs) module doc;
[ssr_pileup.md](ssr_pileup.md) §2 revision note). Every spanning read is now
realigned and stored as **one per-read profile**; there is no confident/ambiguous
split and no `hist_*`. The as-built `.ssr.psp` schema
([registry_ssr.rs](../../src/psp/registry_ssr.rs)) has only `amb-lengths` /
`amb-logliks` (the per-read CSR) — no histogram.

**So off-ladder must follow the same collapse.** There is no `offl_counts`
tally and no separate `offl_amb_*`; off-ladder is just **additional
candidate(s) inside the same per-read profile**, stored all-CSR. This doc
designs *that* shape, not the spec's six-column one. (Spec §4.3 is owed an
amendment to the all-CSR reality — already tracked as a deferred item in
[ssr_stage1_remaining.md](../implementation_plans/ssr_stage1_remaining.md) §3;
this extends that amendment to the off-ladder columns.)

### 0.2 There are **no** off-ladder columns today — `ssr_cohort.md` does not exist

The task brief flags a (future) `ssr_cohort.md` §3 grounding line claiming
"off-ladder columns present in the schema but empty in v1." Confirmed
**inaccurate, twice over**: (a) there is no `ssr_cohort.md` in
`doc/devel/architecture/` at all — the Stage-2 consumer is specified only in
**spec §5.1** (the cohort candidate-set union `A_ℓ`); (b) the schema has **zero**
off-ladder columns, not "present but empty." This task *creates* them. Flag for
correction when `ssr_cohort.md` is eventually written.

---

## 1. What "wiring off-ladder" means, in one picture

The four layers and the single new datum threaded through them — a read's
**off-ladder leg**: at most one read-derived non-rung candidate, sequence-keyed,
jointly normalized with that read's on-ladder rungs.

```
 triage_read ──► analyze_read ───────────► aggregate ──────► registry_ssr (.ssr.psp)
 (unchanged)     (§3 isolate + §4 score)    (§5)               (§6)
   │               │                           │                  │
 region +        Viterbi+traceback vs the    prune+renormalize    new off-ladder
 observed_count  best rung → tract bytes;    over on∪off, split   columns (all-CSR,
                 build_offladder; forward-   into the two legs    Bytes-keyed seq)
                 score it (Allele-agnostic)  of one read's profile
```

Everything rests on one representational decision (already made in
[ssr_shared_types.md](ssr_shared_types.md) §2/§4 and built in
[types.rs](../../src/ssr/types.rs)): an off-ladder allele's identity **is its
verbatim tract sequence** wrapped in a `NormalizedSeq`, `Allele::OffLadder`. This
doc does not reopen that; it wires it through the four Stage-1 layers.

---

## 2. The prior-art anchor (why sequence-keyed, no assembly)

From a read of the vendored sources:

- **HipSTR** keys alleles by **explicit sequence** (`HapBlock.alt_seqs_`,
  `seq_set_`); two reads/samples with byte-identical repeat-region sequence are
  the same allele regardless of motif parity. It *assembles* candidate
  haplotypes across reads (`HaplotypeGenerator::gen_candidate_seqs`).
- **GangSTR** keys alleles by **integer repeat count** (`Locus::allele1/allele2`
  are `int`); an interrupted tract and a clean tract of the same length collapse
  to one allele. Off-ladder detail is lost.

Our design is **HipSTR's sequence-identity without HipSTR's cross-read
assembly** (spec §4.2 "we do not assemble haplotypes across reads"): an
off-ladder allele is captured only when **a single read supports its tract
cleanly enough to emit it**. That is the precision-first sensitivity limit
already documented in spec §4.2 — and it is what makes Stage-1 off-ladder a
*per-read* operation, not a locus-level assembly step.

---

## 3. Layer 1 — isolating the observed tract (one alignment + traceback)

`build_offladder(locus, observed_tract)` needs the read's **observed tract
bytes** — the repeat region *with the flanks stripped*. The pair-HMM **forward**
([pair_hmm.rs](../../src/ssr/pileup/pair_hmm.rs)) gives a *likelihood* (a sum over
all alignments), **not an alignment**, so it does not hand us a boundary. Three
ways to get one were weighed (discussion, 2026-06-17):

1. **CIGAR coordinate mapping** (`ref_to_read` over the read's CIGAR) — *rejected*:
   it re-trusts the indel placement realign-everything distrusts, and fails
   exactly on the impure / interior-indel reads that *are* the off-ladder sources.
2. **Isolated flank anchoring** (match the catalog flanks against the read on
   their own) — *rejected*: a short flank can match by chance, and matching it in
   isolation throws away the global context that pins the boundary.
3. **A single Viterbi (max-path) + traceback against the read's best rung** —
   **chosen.** The full read-to-haplotype alignment anchors the flank/tract
   boundary *globally* (the whole read constrains it, so chance flank-similarity
   can't hijack it), reuses our existing 3-state model, and the tract falls out
   as a byproduct.

triage is **unchanged** — it still classifies coverage from the footprint and
hands `analyze_read` the `region` (flanks + tract) and `observed_count`. Tract
isolation moves **into the scoring layer** (§4): it needs the alignment, which is
produced there.

**The mechanism.** Every rung haplotype is `left_flank + motif×L + right_flank`,
so its flank-boundary columns are known a-priori — the L/tract boundary is
haplotype column `left_flank.len()`, the tract/R boundary is column
`hap.len() − right_flank.len()`. Run Viterbi + traceback of the read region
against the **argmax rung** (the highest-forward-scoring candidate, already
computed in §4); the traceback maps each haplotype column to a read position. The
read bases between the two boundary columns' read positions are the **observed
tract** — interior indels and all. We extract everything *between the anchored
flanks*, so where an internal indel is placed never changes the extracted bytes;
the clean-flank guarantee (Stage 0) makes the flank columns align as matches, so
the boundary read positions are unambiguous in the common case.

**Why the argmax rung, and why one alignment suffices.** The flanks dominate the
boundary placement, so *which* rung we trace against barely moves it — but the
choice must be **fixed** for determinism, and the argmax is the natural,
already-computed one. Because the tract is read off the flank anchors (not the
rung's middle), aligning against a rung whose tract is "wrong" still extracts the
right tract: a 12-unit-plus-1bp read traced against rung 12 places the extra base
as an interior insertion, and the bytes between the anchored flanks are still the
full 13-base tract.

**Determinism — the tie-break rule (pinned here so it isn't lost).** The
traceback must be reproducible, because the boundary read positions determine the
stored tract bytes, and two reads of one molecule must agree (cross-sample
identity, §7). **The rule: on equal Viterbi cell scores, prefer the predecessor
state in the fixed priority `Match > Deletion > Insertion`.** This is purely a
determinism device — the extracted tract bytes are invariant to *internal* indel
placement (we slice between the flanks regardless), so the rule's only job is to
pin the flank-boundary-column → read-position mapping reproducibly. In the
clean-flank common case the boundary column is reached by a `Match` and no tie
arises; the rule covers the pathological remainder.

**This collapses the old O1 and O2 into one decision.** The verbatim off-ladder
key is canonical *because* the tract is "deterministically delimited"
(candidate_generation module doc) — and "deterministically delimited" **is** this
traceback with its fixed tie-break. So the same decision that isolates the tract
also licenses the verbatim `NormalizedSeq` key, and **no separate
`normalize_alleles` step is needed** (the types-doc §4 "reuse the kernel" line is
superseded for the full-tract form — owed an amendment). The remaining
sensitivity limit: a read whose flank runs off its end (allele ≥ read length, spec
§1.4) has no extractable tract and earns no off-ladder leg — the documented
precision-first recall trade, same family as the no-assembly rule (spec §4.2).

**The gate, restated** (re-derived for the realign-everything world, where the
old fast/slow `SlowReason` no longer exists): **a read earns an off-ladder
candidate iff (a) the traceback extracts a tract** (both flank boundaries land
inside the read — no flank ran off the end) **AND (b) that tract is not a pure
motif tiling** (`build_offladder` returns `None` for a tiling — that is a rung).

> **Cost note (deferred, not an open design question).** v1 runs the
> Viterbi+traceback for every spanning read (one extra short DP per read — the
> realign-everything path is already the slow path; per the user, contorting the
> design to avoid it would be premature optimization). A later optimization can
> skip it when the argmax rung's forward score already indicates a clean fit;
> gated on *measuring* it first, like the `count_repeats` fast path.

---

## 4. Layer 2 — `analyze_read`: the mixed candidate set + scoring

**Off-ladder *forward scoring* needs no new scorer** — `score_candidates`
([pair_hmm.rs](../../src/ssr/pileup/pair_hmm.rs)) is already **`Allele`-agnostic**:
it returns `Vec<(Allele, f64)>`, scoring each `CandidateAllele.candidate_seq` by
the forward and reading back its allele key verbatim. An off-ladder
`CandidateAllele` (built by `build_offladder` as `left_flank + tract +
right_flank`, structurally identical to a rung) scored by `forward` comes back as
`(Allele::OffLadder(seq), loglik)`, exactly like a rung. **The one new DP is the
`viterbi_align` traceback (§3)** the tract isolation needs.

**Proposal — `analyze_read` change:**

```rust
TriageResult::Spanning(spanning) => {
    candidates.clear();
    build_rungs(locus, spanning.observed_count, window, candidates);
    let region_seq  = &read.seq[spanning.region.clone()];
    let region_qual = &read.qual[spanning.region.clone()];

    // Forward over the rungs → Qᵣ + (by argmax) the best rung.
    let mut scores =
        score_candidates(region_seq, region_qual, candidates, hmm_scratch, model);

    // Off-ladder: isolate the tract by traceback against the best rung (§3),
    // build the candidate, score it by the same forward.
    if let Some(tract) =
        isolate_tract(region_seq, region_qual, candidates, &scores, locus, hmm_scratch, model)
    {
        if let Some(c) = build_offladder(locus, tract) {
            let off = forward(region_seq, region_qual, &c.candidate_seq, hmm_scratch, model);
            scores.push((c.allele, off)); // appended last
        }
    }
    ReadOutcome::Spanning(scores)
}
```

`isolate_tract` picks the argmax of `scores` (the best rung), runs
`viterbi_align` of `region_seq` against that rung's haplotype, and returns the
read-coordinate tract slice between the two flank-boundary columns — or `None` if
a flank ran off the read end. `pair_hmm.rs` gains `viterbi_align` (max +
backpointer matrix + traceback, with the `Match > Deletion > Insertion`
tie-break) beside `forward`; both share the emission/transition model, but
`viterbi_align` keeps a full backpointer matrix for the one traced haplotype
(relaxing the forward's rolling-two-rows discipline — a tiny, bounded matrix,
§3 cost note). **Reference implementations to start from:** HipSTR's
`HapAligner::align_seq_to_hap` (same domain — read-to-STR-haplotype alignment with
traceback), our own `src/baq/probaln.rs` (in-tree banded pair-HMM with the
matrix/decoding machinery), and GangSTR's `expansion_aware_realign`.

Two correctness notes on the forward, neither requiring code:

- **LCP optimization stays correct.** `score_candidates` shares the forward over
  the candidates' longest common prefix. The off-ladder candidate is *not* in
  that set (it is scored by a standalone `forward` after the fact), so the LCP
  logic is untouched. (Were it ever added to the set, it shares the left flank, so
  the LCP stays ≥ `left_flank.len()` — never a wrong score, only less reuse.)
- **Candidate order.** `build_rungs` emits ascending by length; appending the
  off-ladder score **last** keeps the on-ladder run ascending, which the
  aggregator relies on (§5). Its position is irrelevant once the aggregator splits
  by `Allele` variant.

`build_offladder` / `normalize_offladder` are wired from production for the first
time (only their own tests call them today). The `Qᵣ` returned by `analyze_read`
is now a mixed dense vector: `(OnLadder, ll)×(2w+1)` plus, when earned, one
`(OffLadder, ll)`.

---

## 5. Layer 3 — `aggregate`: joint pruning/renormalization + the record field

### 5.1 The in-memory record field shape

Today: `SsrLocusRecord.spanning: Vec<Vec<(u16, f32)>>` — one on-ladder profile
per spanning read. We must add the off-ladder leg.

Each read earns **at most one** off-ladder candidate (§3/§4), and a read's two
legs are **jointly normalized** (§5.2) — they are two parts of one read's
distribution. The shape must keep that coupling explicit. Three options:

| option | shape | verdict |
|---|---|---|
| **A: parallel vec** | keep `spanning`, add `spanning_offladder: Vec<Vec<(NormalizedSeq, f32)>>` (same length, usually-empty rows) | least churn; coupling is positional/implicit |
| **B: per-read struct (recommended)** | `spanning: Vec<SpanningProfile>` where `SpanningProfile { on_ladder: Vec<(u16,f32)>, off_ladder: Vec<(NormalizedSeq,f32)> }` | coupling is structural; names the noun the task asks for; touches the container mapping |
| **C: one mixed vec** | `spanning: Vec<Vec<(Allele, f32)>>` | smallest type, but re-mixes the two storage regimes the container must split anyway, and pushes the on/off split into every consumer |

**Recommendation: B.** A read's profile *is* a single object with two legs; a
struct says so. The `off_ladder` vec is `0..=1` long in v1 but a `Vec` keeps the
door open (and the container CSR handles it uniformly). The user's naming note
(nouns) gives the type name — `SpanningProfile` (or, for the off-ladder leg
alone, `OffLadderProfile`).

> **Open question O3 — A vs B.** B is cleaner but ripples into `registry_ssr`'s
> `SsrLocusRecord` (the container mirror) and the `append`/`next_record` mapping.
> A is a smaller diff. Recommendation B for clarity; flag the extra surface.

### 5.2 Joint pruning and renormalization

`prune_and_renormalize` currently (a) drops candidates > `AMB_LL_DROP` nats below
the per-read max, (b) renormalizes survivors to log-probs summing to 1, and (c)
`debug_assert!(false)` on any `OffLadder` (the "not yet wired" tripwire).

**The change: prune and renormalize over the *union* of on- and off-ladder
candidates, then split the survivors into the two legs.** This is the load-bearing
statistical invariant:

> A read's stored mass sums to 1 **across both allele kinds together** — the
> off-ladder candidate competes with the rungs for the same probability budget.
> The `logsumexp` denominator includes the off-ladder leg; the `AMB_LL_DROP`
> cutoff is applied to the joint max.

Consequences:

- A read whose best explanation is off-ladder keeps the off-ladder leg at high
  mass and prunes far rungs — correct.
- A clean on-ladder read whose off-ladder candidate (if any) scores far below the
  best rung prunes the off-ladder leg to empty — correct, and the common case.
- **Pruning can empty the off-ladder leg but never the on-ladder leg** (the best
  candidate, which always survives, may be either kind; but at least one
  candidate survives, so a profile is never wholly empty unless the read had no
  candidates — which spanning reads always have, ≥1 rung).

The split is by `Allele` variant: `OnLadder { units }` → the `(u16, f32)` leg,
`OffLadder(seq)` → the `(NormalizedSeq, f32)` leg. Order within each leg is
preserved (rungs stay ascending).

> **Determinism / value change.** This **changes the on-ladder `f32` values**
> versus today's on-ladder-only output, because the renormalization denominator
> now includes off-ladder mass wherever a read earns an off-ladder leg. That is
> intended (off-ladder was previously unwired and silently dropped). The
> **across-thread-count byte-identity invariant still holds** (§7): the
> computation is per-read and order-preserving; nothing about it depends on thread
> count.

---

## 6. Layer 4 — container schema (`registry_ssr`): the on-disk format change

This is an **on-disk `.ssr.psp` format change** (new columns). The format is
pre-alpha, no backwards-compat (no production `.ssr.psp` exist — Stage 1 has no
driver yet). The change is additive and rides **existing wire codecs** — no new
block-level encode/decode primitives — keeping faith with §10.4's "a table + a
record mapping, no new wire code."

### 6.1 What the container already gives us

- **CSR ragged lists** (`encode/decode_list_column_csr`) — used by `amb-lengths`
  (u16) / `amb-logliks` (f32), grouped per profile via shared `amb_offsets`.
- **Variable-length byte sequences** — `ColumnPayload::Bytes { length_column }`:
  a flat byte slab chunked by a paired per-element length column, plus the CSR
  sibling `decode_bytes_split_csr`. **This is exactly how the SNP schema stores
  `allele-seq` (chunked by `allele-seq-len`)** — the prior art for variable-length
  sequence keys. Off-ladder sequences ride this, not a new codec.
- `MAX_ALLELE_SEQ_LEN` — the existing per-sequence byte cap the Bytes decoder
  enforces; off-ladder sequences reuse it (a tract+nothing is far under it).

### 6.2 Proposed columns (per-profile, all-CSR, sequence-keyed)

Mirroring `amb-*` but for the off-ladder leg. New `SsrColumnKey` tags in the
per-profile range `[0x12, 0x1F]`; new `ColumnDef`s. The off-ladder leg per
profile has **0 or 1 entries** in v1, so these CSRs are mostly-empty (zstd
collapses them to ~nothing — consistent with the all-CSR "let zstd dedup"
decision, [locus_record.rs](../../src/ssr/pileup/locus_record.rs) module doc).

| column | tag | payload | grouping | holds |
|---|---|---|---|---|
| `offl-logliks` | `0x12` | `List<F32>` (CSR) | per profile (shares an `offl_offsets`, len = n_profiles+1) | each profile's off-ladder log-prob(s) — 0 or 1 |
| `offl-seq-len` | `0x13` | `List<Varint>` (CSR) | per profile (same `offl_offsets`) | byte length of each off-ladder sequence — 0 or 1 |
| `offl-seq` | `0x14` | `Bytes { length_column: "offl-seq-len" }` | per off-ladder *entry* | the verbatim tract bytes, chunked by `offl-seq-len` |

Decode chains cleanly: `offl-logliks` (CSR) gives the per-profile offsets whose
last value is the total off-ladder entry count `E`; `offl-seq-len` (CSR, same
offsets) gives `E` byte-lengths; `offl-seq` (Bytes) is chunked by those `E`
lengths. No per-locus dict (the spec's `offl_seqs` dict is superseded by the
all-CSR "store verbatim, zstd dedups" choice — §0.1).

> **Open question O4 — the one codec wrinkle: the `Bytes` length column is a
> *CSR list*, not a flat scalar.** The SNP `allele-seq-len` is a flat per-allele
> scalar (one length per allele). Here `offl-seq-len` is a per-profile CSR list
> (0/1 lengths per profile). The `Bytes` chunking only needs the **flattened**
> length stream (`E` lengths in order), which the CSR data slab already is — so
> the Bytes decoder consumes `offl-seq-len`'s flat data and the CSR offsets are
> used only to regroup logliks/lengths into profiles. This is the single piece to
> get right in the implementation plan; it is a wiring detail, not a new codec.
> **Alternative considered:** flatten off-ladder into a per-locus entry stream
> with a `n-offl` per-record grouping count (like `n-spanning` groups profiles).
> Rejected — it decouples the off-ladder leg from its read/profile, and Stage 2's
> per-read likelihood needs the leg tied to the same read's on-ladder profile
> (spec §5.1, `P(read|a) = Σ_L Qᵣ(L)·S_θ`). Per-profile grouping is mandatory.

> **Open question O5 — required vs optional columns.** Recommendation:
> **required** (always emitted, like `amb-*`), so the manifest invariant
> (`decode_block`'s B1 check) and the decoder stay unconditional. An all-empty
> off-ladder column at a block with no off-ladder zstd's to a few bytes. The
> decoder already tolerates unknown/optional columns (skips them), but making
> these required avoids a "some files have it, some don't" matrix.

### 6.3 Writer / decoder / structural-checks

- `SsrBlock` gains `offl_offsets: Vec<u32>` (shared by the two CSR lists),
  `offl_logliks_data: Vec<f32>`, `offl_seq_len_data: Vec<u32>`,
  `offl_seq_bytes: Vec<u8>`; `append` walks each profile's off-ladder leg exactly
  as it walks the on-ladder leg today.
- `SsrDecoder` mirrors with reuse buffers; `next_record` reassembles each
  profile's `off_ladder: Vec<(NormalizedSeq, f32)>` from the slabs.
- **Structural checks** (extending the M2/M3 family already in `registry_ssr`):
  `offl_offsets.len() == n_profiles + 1`; the `offl-seq-len` CSR offsets agree
  with `offl-logliks`'; `sum(offl-seq-len data) == offl-seq.len()`; each sequence
  ≤ `MAX_ALLELE_SEQ_LEN`. A NaN off-ladder loglik is rejected on write (the
  existing `NonFiniteLoglik` path generalizes; `-inf` permitted as a legit
  `log(0)`).

### 6.4 Format version

Bump the **ssr `schema_version`** (additive columns). The shared
**`container_version` is unchanged** — no new wire codec. Both `SsrLocusRecord`
types (the in-memory `pileup` one and the container mirror) gain the off-ladder
field; the Stage-1 driver's name→chrom_id adapter
([ssr_stage1_remaining.md](../implementation_plans/ssr_stage1_remaining.md) §2.3)
carries it across, exactly as it carries `spanning`.

---

## 7. Cross-sample identity & invariants

- **Cross-sample identity (the property Stage 2 §5.1 consumes).** An off-ladder
  allele's key is its verbatim tract `NormalizedSeq`. Two reads of the same
  molecule — in the same sample or across samples — isolate the same tract bytes
  (Stage 0's clean-flank guarantee makes the tract deterministically delimited),
  so they produce **byte-identical keys** → `Allele::OffLadder` `==` and `Hash`
  are equal → Stage 2's candidate-set union (`A_ℓ`, spec §5.1) folds them into
  **one** cohort allele. This is already unit-asserted
  (`same_tract_two_reads_gives_an_identical_key`); wiring does not weaken it. The
  spec §5.1 off-ladder union ("unioned across samples by their normalized key")
  is precisely what these columns feed.
- **Determinism (the load-bearing Stage-1 invariant).** `.ssr.psp` must be
  byte-identical across `--threads ∈ {1, N}`
  ([ssr_stage1_remaining.md](../implementation_plans/ssr_stage1_remaining.md) §2.5;
  [ssr_pileup.md](ssr_pileup.md) §8.4). Off-ladder wiring preserves it: candidate
  generation, scoring, and joint renormalization are **pure per-read functions**
  with fixed candidate ordering; nothing depends on thread count or read-arrival
  order beyond what determinism already pins. (As noted in §5.2, the *values*
  change vs the current on-ladder-only output — intended — but split-invariance is
  untouched.)
- **Dense writing unchanged.** One record per covered locus; a locus with no
  off-ladder simply has all-empty off-ladder CSR rows. No new no-call/sparse
  behavior.

---

## 8. Testing (spec §7 anti-tautology; Bucket-1 on the critical path)

- **Unit — `viterbi_align` + tract isolation:** the traceback against a rung
  extracts the correct tract bytes (clean read → the clean tiling; interior-indel
  read → the indel-bearing tract; a read with a flank off the end → `None`);
  determinism of the `Match > Deletion > Insertion` tie-break (two equal-scoring
  alignments resolve identically). Reuse / cross-check against the HipSTR/BAQ
  reference patterns where useful.
- **Unit — candidate generation:** already present for `build_offladder` /
  `normalize_offladder` (verbatim key, pure-tiling → `None`).
- **Unit — analyze_read:** an interior-indel read earns an off-ladder candidate
  and scores it highest; a clean read earns none (off-ladder leg empty); a read
  whose flank runs off the end earns none (the §3 gate).
- **Unit — aggregate:** joint renormalization sums to 1 across both legs; a read
  bimodal between a rung and an off-ladder sequence keeps both, normalized; the
  off-ladder-prunes-to-empty case.
- **Round-trip — `registry_ssr`:** extend `ssr_round_trips_through_the_container`
  with profiles carrying off-ladder legs (including a locus with none, a profile
  with one, mixed); assert per-record equality. Plus the structural-rejection
  tests (mismatched offsets, oversized sequence, NaN loglik).
- **Determinism:** the §2.5 byte-identity gate covers off-ladder once a driver
  exists; until then, an aggregate-level property test (shuffled outcome order →
  identical record, since order is by read not thread) guards the unit.
- **Anti-tautology:** the simulator's off-ladder injection (a 12 vs 12+1bp cohort,
  spec §7) is the end-to-end guard — kept separate from the caller's code; it
  belongs to the simulator pass but the column shape here must support it (two
  distinct alleles, reconciled by key, not collapsed).

---

## 9. Decisions vs open questions (for discussion)

**Proposed decisions (argue with these):**

- **D1** Off-ladder is **per-read profile entries**, all-CSR — *not* the spec's
  six-column histogram+CSR mirror (superseded by the all-CSR revision, §0.1).
- **D2** A read **earns** an off-ladder candidate iff (a) the traceback extracts a
  tract (no flank ran off the read end) **and** (b) that tract is non-tiling — the
  realign-everything re-derivation of the old `SlowReason` gate (§3).
- **D3** The tract is isolated by **one Viterbi + traceback against the argmax
  rung** — not CIGAR coordinate mapping (re-trusts indel placement), not isolated
  flank anchoring (chance flank match) (§3). *(decided in discussion, was O1.)*
- **D4** Traceback tie-break is the fixed priority **`Match > Deletion >
  Insertion`**; this deterministic delimitation is what licenses the **verbatim**
  `NormalizedSeq` key, so **no separate `normalize_alleles` step** is needed (the
  types-doc §4 "reuse the kernel" line is superseded for the full-tract form and
  owed an amendment) (§3). *(decided in discussion, was O2; collapses into D3.)*
- **D5** Off-ladder *forward scoring* reuses `forward`/`score_candidates`
  unchanged (already `Allele`-agnostic); the off-ladder score is appended last.
  The **only new DP** is `viterbi_align` for D3's traceback (§4).
- **D6** Joint prune/renormalize over on∪off, then split by variant (§5.2).
- **D7** Schema rides the existing `Bytes`/CSR codecs; ssr `schema_version` bump,
  `container_version` unchanged (§6).

**The open questions I want your input on (load-bearing):**

- **O3 — record field shape (§5.1).** `SpanningProfile` struct (B, recommended,
  more surface) vs a parallel `spanning_offladder` vec (A, smaller diff)?
- **O4 — the CSR-length-column Bytes wrinkle (§6.2).** Confirm the
  per-profile-CSR + flat-Bytes-chunking mapping (recommended) over a decoupled
  per-locus off-ladder entry stream (rejected for breaking the read↔leg
  coupling)?
- **O5 — required vs optional off-ladder columns (§6.2).** Required/always-emitted
  (recommended) vs optional?

*(O1 the tract-isolation gate and O2 verbatim-vs-kernel are now resolved →
D3/D4.)*

**Deliverable for this step:** this draft. Stop here; iterate to agreement,
then an implementation plan, then code (incremental: triage → analyze_read →
aggregate → schema, each landed green before the next, per the user's working
style).
