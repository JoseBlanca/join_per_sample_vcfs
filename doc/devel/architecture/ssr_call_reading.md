# SSR Stage 2 — `ssr-call` reading & merge (architecture sketch)

**Status:** draft, 2026-06-19, branch `ssr-cohort`. **A discussion starter, not a
settled design** — the first of three sketches that together cover Stage 2
(`ssr-call`), the cohort caller. This one is **Phase 1: read the N per-sample
`.ssr.psp` files, decompress lazily, k-way-merge by locus into the `CohortLocus`
work-item, and feed it through the producer→queue→worker→writer topology.** The
two companions:

- [ssr_call_parameters.md](ssr_call_parameters.md) — Phase 2: the pre-pass that
  freezes `ε` (per sample group) and builds the stutter-shape / stutter-level / `G₀`
  priors (the "fixed parameters + priors").
- [ssr_call_genotyping.md](ssr_call_genotyping.md) — Phase 3: candidate assembly,
  the HipSTR likelihood, the per-locus EM, the `F` loop, and the VCF.

It is the *how-the-code-is-wired* companion to spec
[ssr_cohort_mark2.md §4.1](../specs/ssr_cohort_mark2.md) (reading & orchestration —
**intent settled 2026-06-19**; this doc proposes the module/struct shape that §4.1
deferred). Where the two disagree on intent, the spec wins; on code layout, this
doc wins. It refines the Mark-1-era `cohort/` sketch in
[ssr_genotyping_architecture.md §4/§6](ssr_genotyping_architecture.md) for the
Mark-2 model.

---

## 1. What this phase produces

```
 sample_1.ssr.psp ─┐  per-sample      ┌──────────── producer thread ────────────┐
 sample_2.ssr.psp ─┤  SampleEvidence  │ for each locus in catalog order:         │
       ...         ├─►  Cursor (lazy ─► │   ask EVERY cursor → evidence_at(locus)  │ ─► bounded queue
 sample_N.ssr.psp ─┤  +prefetch pool)  │   gather the Some(..) responses          │    of CohortLocus
                   │                   │   ≥1 present → push ONE CohortLocus       │   (one at a time)
 .ssr.catalog ─────┘  (motif, frame)   │   all absent → SPARSE-OMIT (drop)        │       │
                                       └──────────────────────────────────────────┘       ▼
                                                                          EM worker pool ─► writer
                                                                          (Phases 2-3)      (reorder by
                                                                                             locus seq → VCF)
```

The deliverable of Phase 1 is the **`CohortLocus`** and the machinery that streams
them in catalog order with bounded memory. Phases 2–3 consume it.

---

## 2. The work-item — `CohortLocus`

Spec §4.1's "analysis unit." The catalog **frame** + the present samples' evidence:

```rust
// shape sketch — not final
struct CohortLocus {
    // --- from the catalog (Stage 0), the coordinate frame, NOT an allele claim ---
    locus: LocusId,            // (chrom, start, end) — the cohort-match key
    motif: Motif,              // the stutter unit (Phase 3 S2)
    ref_tract: Box<[u8]>,      // reference tract + flanks: the Δ frame + align frame
    // --- present samples only (sparse); the two vecs are PARALLEL, same length ---
    present: Vec<u32>,             // cohort sample indices (into 0..N), ascending; absent omitted
    samples: Vec<SampleEvidence>,  // evidence[k] belongs to sample present[k]  (§4.1)
}

struct SampleEvidence {
    seq_counts: Vec<(Box<[u8]>, u32)>, // distinct repeat-region seq -> read count, sorted by bytes
    qc: SsrQc,                          // depth, mapped_reads, n_filtered, n_low_quality, n_border_off_end
}
```

This is exactly the Stage-1 `SsrLocusObs`
([locus_tally.rs](../../../src/ssr/pileup/locus_tally.rs)) per present sample, plus
the catalog frame, plus the absent/present vector. **Tiny per locus** (a handful of
distinct sequences each) — so the cost of this phase is **decompression**, not
assembly (spec §4.1), which is what the worker pool exists to overlap.

> **Q-R1 — RESOLVED (sparse, struct-of-arrays).** `CohortLocus` holds only the
> present samples: `samples: Vec<SampleEvidence>` paired with a parallel
> `present: Vec<u32>` of cohort sample indices (into the fixed `0..N` file order),
> ascending. No `Option`/`None` gaps. The index→name (`SampleId`) map is identical
> for every locus, so it lives **once on the driver** (VCF column headers), not
> copied per locus. The per-individual `F_i` reduce (Phase 3) iterates the present
> vector and attributes evidence `k` to individual `present[k]` — more direct than
> scanning N `Option`s.

---

## 3. The per-sample reader — `SampleEvidenceCursor`

Spec §4.1: *a coordinate cursor with a one-block decode cache.* **One per input
file (N total)** — it knows only *its own* sample. **Responsibility:** own the
per-sample cache + position and decide *when* to refill; it **delegates** the actual
zstd+column decode to the reused columnar layer (`BlockColumnReader`, table below).
The cross-sample gathering is the merger's job (§4), not the cursor's — the cursor's
product is one sample's `SampleEvidence`, never the locus-wide `CohortLocus`.

**State (per sample):**

- the **block index** in memory (cheap — `BlockIndexEntry`s, already loaded by
  `PspReader`) + a forward within-block position;
- **two block slots** — `current` (being consumed) and `next` (prefetched ahead, see
  *Decompression* below) — each holding a decoded block *or* a pending decode future;
- a handle to the **shared decode pool** (the same pool given to every cursor);
- **`held`** — the cursor's *next stored locus*, decoded and buffered ahead (`None`
  once the sample is exhausted). The cursor "always holds the next locus";
- **`last_query`** — the last locus the merger asked for: a **monotonic guard**
  (the differentiator below).

**The contract.** `evidence_at(q)` is called by the merger with loci in **strictly
ascending catalog order**. `held` is preloaded in `new` and re-filled after every hit
(via `advance()`), so the cursor always has its next stored locus ready — or `None`:

Decode can fail lazily on a block refill, so the real return is
`Result<Option<SampleEvidence>, _>` (the `Option` is Present/Absent; the `Result` is a
decode/contract error):

```rust
fn evidence_at(&mut self, q: LocusId) -> Result<Option<SampleEvidence>, SsrCohortReadError> {
    // monotonic guard — q must strictly ascend. The merger validates catalog sortedness,
    // so a rewind can't come from on-disk data; this is belt-and-suspenders (debug only).
    debug_assert!(self.last_query.map_or(true, |p| q > p), "out-of-order / rewind");
    self.last_query = Some(q);

    match self.held {                  // held = next stored locus; None once exhausted
        None               => Ok(None),                       // exhausted → Absent forever after
        Some(f) if q == f  => { let e = self.take_held(); self.advance()?; Ok(Some(e)) }
        Some(f) if q  < f  => Ok(None),                        // forward, sample lacks q → Absent
        Some(_) /* q > f */ => Err(LocusNotInCatalog { .. }),  // sample has a locus the catalog lacks
    }
}
```

The three outcomes:

1. **`q == held`** → `Some(evidence)`, then `advance()` to preload the next stored
   locus.
2. **`q > held`** (asked *past* our front) → **hard error** (`LocusNotInCatalog`): given
   the merger walks a *sorted* catalog (it validates this — §4), the only way it asks
   past a locus this sample holds is that the sample contains a locus the catalog does
   not — i.e. it was genotyped against a **different catalog**. A typed error, never a
   panic, never silent — the "same catalog for all samples" violation is a hard error.
3. **`q < held`** (asked for something *before* our front) → **`None` (Absent)**:
   this sample simply has no record at `q`.

**Why the `last_query` guard is what makes case 3 unambiguous.** Without it, `q < held`
has *two* meanings — "this sample genuinely lacks `q`" (Absent) and "the merger went
backwards" (a bug) — and the cursor can't tell them apart from the `q`-vs-`held`
comparison alone. The guard splits them: a rewind is `q <= last_query`, caught at the
**top**, *before* the match. So by the time control reaches the `q < held` arm, `q` is
guaranteed a forward query — which can only mean Absent. Cost: one `Option<LocusId>`
field and one comparison.

**Dense or sparse — the same code.** Today's catalog is non-sparse (every sample
stores every catalog locus), so the `q < held` arm never fires (a forward query always
lands `q == held`) and the guard merely catches rewinds. If the catalog ever becomes
sparse, the *same* `evidence_at` already returns the correct `None`. The contract is a
**strict superset** — build it dense-first and it is already sparse-correct, no rewrite.
Keep the guard even in the pure-dense design: cheap insurance that turns an out-of-order
merger bug into a loud panic instead of silently-wrong evidence.

**`advance()` — refill `held`, crossing blocks as needed.** Called from `new` (preload
the first stored locus) and after each hit. It moves the within-block position to the
next stored record and decodes it into `held`. When that record lies in a later block it
promotes `next → current` and immediately fires a fresh prefetch (below). No more records
anywhere → `held = None`. zstd and blocks are the cursor's private business; the merger
never sees them.

**Decompression — shared pool + prefetched futures.** Decode (zstd + columnar) is the
heavy work, and we do *not* run it inline on the merger thread. Each cursor holds a
handle to **one decode pool shared across all N cursors**, and each block slot is either
a decoded block or a **pending future** — concretely a `oneshot`/crossbeam receiver that
`recv()` returns instantly if the worker already finished, or **parks on if not** (no
async runtime — fits the existing rayon/crossbeam stack):

```rust
enum BlockSlot {
    Decoded(DecodedBlock),
    Pending(Receiver<DecodedBlock>),  // the "future"; recv() = ready-now or park
    Empty,
}
// cursor: { pool: Arc<DecodePool>, current: BlockSlot, next: BlockSlot, ..index, pos, held, last_query }
```

The cursor **self-double-buffers**: the moment it starts consuming `current`, it submits
the *next* block's decode to the pool and stows the receiver in `next` as `Pending`. When
`advance()` crosses into that block it `recv()`s `next` — **usually already `Decoded`**,
because the prefetch overlapped the downstream EM work and the other cursors — promotes it
to `current`, and fires the next prefetch. The blocking `recv()` is the **rare degraded
path** (pool saturated), not the norm. The decode task is a **stateless pure transform**
(compressed bytes + column set → decoded columns); the worker writes only the `oneshot`
sender, never the cursor — so there is no shared-mutable cursor state and no lock. Decode
**reuses its target buffer** rather than reallocating ([[feedback_scratch_memory]]), and
decodes *all* columns (no light/heavy two-phase, §6).

> **No "loci per block" assumption.** This works whether a block holds one locus or
> thousands: the prefetch overlaps decode with the *pipeline* (EM + sibling cursors via
> the bounded queue, §5), not with the current block's interior — so fine-grained blocks
> don't collapse it. Memory cost is `current + next` = **N × 2 blocks** resident (vs.
> N × 1, §7), tunable by prefetch depth.

**Invariant Stage 1 must honour (spec §4.1):** a locus lives **wholly inside one
block** — its writer must never split a locus across a block boundary. That is what
lets the cursor treat a locus as atomic and skip all the SNP straddler machinery
(§6).

| reused plumbing | path | role |
|---|---|---|
| `PspReader` + block index | [src/psp/reader.rs](../../../src/psp/reader.rs), [index.rs](../../../src/psp/index.rs) | open file, load `BlockIndexEntry`s, seek a block; **interval-keyed `last_pos`** (arch §10) is what the covering-block test needs |
| `BlockColumnReader` / columnar decode | [src/psp/reader.rs](../../../src/psp/reader.rs) | decompress one block's columns |
| SSR typed decode | [src/psp/registry_ssr.rs](../../../src/psp/registry_ssr.rs) (`SsrColumnKey`, `BlockDecoder`) | columns → per-locus `(seq,count)` + QC → `SsrLocusObs` |
| same-catalog check | header md5 ([common.rs](../../../src/pop_var_caller/common.rs) `FastaVerify` convention) | spec §4.1 hard error on catalog-md5 mismatch, checked once at open |

**What SSR does *not* reuse from the SNP cohort reader**
([sample_reader.rs](../../../src/var_calling/sample_reader.rs)): `next_two_phase`,
the light/heavy `TwoPhaseSegment` fold, the straddler/safe-gap/watermark logic, and
position streaming. SSR's unit is *given by the catalog* (not discovered), atomic,
and wholly in-block — so the cursor is much thinner than `SamplePspReader`. We reuse
the **columnar layer beneath** it (the generic core), not the SNP cohort reader.

---

## 4. The k-way merger (the producer)

Spec §4.1. Walks the **catalog** (the authoritative ordered master locus list; each
sample file is a coordinate-ordered subset), and for each locus:

- calls `evidence_at(locus)` on **every** cursor (monotonic, so each cursor advances at
  most forward);
- gathers the `Present` responses; **≥1 present → build a `CohortLocus`**; **every
  cursor absent → drop it, never emit** (sparse-omit);
- pulls the catalog frame (motif, ref tract) for that locus from the catalog reader;
- pushes the `CohortLocus` onto the **bounded queue** — **one locus at a time**,
  tagged with a monotonic locus sequence number for the writer to reorder by.

> **Keep it simple — no batching yet.** We do *not* accumulate loci into batches of
> `K`. One `CohortLocus` per queue item is the starting point; the bounded-queue
> back-pressure and the per-locus EM cost already balance the pipeline. If profiling
> later shows queue/channel traffic dominates or workers are starved by per-item
> latency, batching is an easy, localized add (group `K` loci per item; the SNP path's
> `target-variants-per-chunk` is the precedent) — deferred, not designed in now.

It is a **plain coordinate k-way merge**, not position streaming — loci are
disjoint atomic intervals in a shared order, so "merge" is just "ask all cursors for
the same key and advance." The catalog is the merge driver; the per-sample files are
followers.

> **Q-R2 — RESOLVED: catalog-driven.** The merge iterates *all* catalog loci in their
> authoritative order and sparse-omits the empty ones — simplest, deterministic, and
> the catalog is already ordered so no k-way heap is needed. The alternative (heap-merge
> of the N files' present loci, visiting only covered loci) is rejected: the empty-locus
> visits are near-free anyway (each is N cheap `evidence_at`s answering `None`, mostly
> off blocks already in hand), so they don't justify the heap. On today's non-sparse
> catalog (every sample stores every locus) there are no empty visits at all.

---

## 5. Execution topology

Mirrors the SNP cohort pipeline: **main producer + W workers + 1 writer** (spec
§4.1).

```
 Cursors (N) ─prefetch decode──►┐
   ▲ evidence_at                 ├─► shared worker pool ── decode tasks (priority) + EM tasks
   │                             │         │
 merger (1 producer) ────────────┘         │ EM(locus)
   k-way by catalog locus ─► bounded queue ─┘ ─► writer (1): reorder by locus seq → catalog VCF
   builds one CohortLocus      one locus at a time, back-pressure
```

- **Decompression is offloaded to the pool, not run on the producer thread.** Cursors
  submit their next-block decodes as **prefetched futures** to a **shared decode pool**
  (§3) — so the heavy zstd+columnar work overlaps the EM and the sibling cursors
  through the bounded queue, and the merger thread only ever `recv()`s an
  already-(usually)-decoded block. **No "loci per block" assumption** — this holds for
  one-locus blocks as for thousand-locus blocks (§3).
- **One *shared* pool, with decode-priority — *not* a rigid decode/EM split.** §5 used
  to call for a dedicated EM pool walled off from decode (to dodge the SNP path's
  decode-starvation, [compact_samples / pipeline balance], [worker bucket profile]).
  Prefetch dissolves that: because the future gives decode **lead time**, decode is
  latency-insensitive and need not race the EM — so it can ride the *same* pool as a
  higher-priority task, which avoids the steep mis-split penalty of a fixed partition
  ([thread-budget split sweep]). Giving decode tasks priority is the cheap insurance
  that a burst can't sit behind queued EM work longer than its lead time.
- **Output ordering is trivial.** Loci are independent and catalog-ordered, so the
  writer reorders finished loci by **locus sequence number** (a `BTreeMap` keyed on
  seq, exactly as the SNP writer does on chunks) → VCF in catalog order.
- **Bounded queue depth** caps peak resident loci → the memory knob.

> **Open Q-R3 — queue depth.** Queue depth trades throughput for peak RSS — a tuning
> knob to set against a real run, not pre-defaulted. (Batching `K` loci per item is
> *not* in scope yet — §4; if it's ever added, its size becomes the second knob here.)

---

## 6. Why SSR drops three SNP reader layers (spec §4.1)

The SNP unit (a variant group) is **discovered from the data**; the SSR unit (a
locus) is **given by the catalog**. So SSR omits:

- **no light/heavy two-phase decode** — there is no "is this row variable?" fold;
  every covered locus is genotyped, the SSR schema is already flat, so decode *all*
  columns;
- **no straddler / safe-gap / watermark** — a locus is atomic and wholly in one
  block (the Stage-1 writer invariant, §3);
- **the merge is a plain coordinate k-way merge**, not position streaming.

This is why the cursor (§3) and merger (§4) are markedly simpler than their SNP
analogs — and why the sketch is short.

---

## 7. Bounded memory & determinism

- **Lockstep ⇒ bounded working set.** Every cursor is asked for the same locus in
  order, so each pins only the block(s) covering the front → resident set ≈
  **N × 2 decompressed blocks** (`current` + one prefetched `next`, §3; N × 1 if
  prefetch is off), the cohort-scaling property (spec §4.1). `--block-window-bp` is
  the RSS lever exactly as on the SNP path ([block window memory lever]); prefetch
  depth is the other.
- **Lockstep is approximate today.** Byte/count-gridded Stage-1 blocks make refills
  land at slightly different boundaries (~1–2 blocks/sample). Moving the Stage-1
  writer to **genomically-aligned blocks** (a future task, mirroring the SNP win)
  makes it exact. *Not a Phase-1 blocker* — flag it (memory thesis: keep the knob
  visible).
- **Determinism.** Reading + merge are a pure function of the inputs + catalog
  order — no per-thread state on this path, so the `CohortLocus` stream is identical
  across thread counts. The shared decode pool does **not** change this: decode is a
  pure transform and the merger consumes blocks strictly in catalog order, so *when* a
  block decodes (and on which thread) affects only timing, never content. (The producer
  is single-threaded by design; the workers' determinism is Phases 2–3's concern,
  carried by per-sample-group-frozen `ε` + per-locus stutter shape `θ_locus`.)

---

## 8. The two-pass question (the one real cross-phase wrinkle)

Phase 2 (the pre-pass) needs to stream loci to estimate `ε` (per sample group) / the
stutter-shape (cohort, per period) + stutter-level (per sample group) priors *before*
Phase 3 can genotype any locus. So the merge stream (§4) is consumed **twice**: once
for the pre-pass, once for genotyping.

**RESOLVED: re-read** — reset all cursors and run the merge again for Phase 3. The
rejected alternative was to **cache** the `CohortLocus` stream in memory between
passes (decompress once, pay RAM); caching the whole stream is **too memory-intensive**
and breaks the cohort-scaling thesis, so we keep it simple and pay decompression
twice. It's cheap to do so: spec §4.4's pre-pass touches only a **representative
stratified subset** of loci, far fewer than genotyping, so the first pass is the
small one. The merge layer must therefore be **cheap to restart** — cursors re-seek
from their block index, no global state to rebuild (the cursor's `new`/`advance` and
the shared decode pool re-spin from scratch).

> This is the **module-boundary handoff to [ssr_call_parameters.md](ssr_call_parameters.md)**:
> the pre-pass is a *consumer* of the same merger, run first. Whether it shares the
> producer/queue topology (§5) or runs as a simpler serial sweep over a subset is an
> open question that doc owns.

---

## 9. Proposed module layout (sketch)

```
src/ssr/cohort/            # Stage 2 (ssr-call) — all three phases live here
  mod.rs
  types.rs       # CohortLocus, SampleEvidence, LocusId  (Phase 1 spine; shared with 2-3)
  reader.rs      # SampleEvidenceCursor: block-index + current/next slots + prefetch + evidence_at(locus)  (§3)
  decode_pool.rs # shared decode pool + DecodedBlock futures (§3, §5)
  merge.rs       # k-way catalog-driven merger → one CohortLocus at a time  (§4)
  driver.rs      # producer + bounded queue + worker pool + writer wiring  (§5);
                 #   hosts the two-pass orchestration (§8)
  # ... parameters.rs / candidate_set.rs / em.rs / vcf_out.rs land via docs 2-3
```

(`cohort/` keeps the name from [ssr_genotyping_architecture.md §4](ssr_genotyping_architecture.md);
the CLI subcommand is `ssr-call`.)

---

## 10. Open items (Phase-1 agenda)

- **Q-R1** — ~~dense `Vec<Option>` vs sparse~~ **RESOLVED**: sparse SoA — `samples: Vec<SampleEvidence>` + parallel `present: Vec<u32>` of cohort indices (§2).
- **Q-R2** — ~~catalog vs union-of-files driving~~ **RESOLVED**: catalog-driven (§4) —
  iterate all catalog loci, sparse-omit empties; no k-way heap.
- **Q-R3** — queue depth default (§5), measured. (Batching of `K` loci is **deferred**,
  not designed in — §4; one `CohortLocus` per queue item to start.)
- **Q-R4** — ~~producer-thread vs fanned~~ **RESOLVED (design)**: decode is offloaded
  to a **shared decode-priority pool** as **prefetched futures**, self-double-buffered
  in the cursor (§3); decode rides the same pool as the EM because prefetch makes it
  latency-insensitive (§5). *Build* gated on profiling: ship a simple version + deep
  queue first, turn prefetch on when the decode burst bites (interacts with Q-R6).
- **Q-R5** — ~~re-read vs cache~~ **RESOLVED**: re-read (caching the stream is too
  memory-intensive, §8). *Still open (parameters doc):* whether the pre-pass reuses
  this producer/queue topology or runs a simpler serial sweep over the subset.
- **Q-R6** — genomically-aligned Stage-1 blocks to make lockstep exact (§7) — a
  Stage-1-writer follow-up, tracked, not a Phase-1 blocker.

**Settled upstream (spec §4.1), not re-opened here:** absent = no-data-for-this-
sample; sparse-omit all-absent loci; lockstep bounded memory; writer reorders by
locus seq; SSR drops the three SNP reader layers (§6). (The earlier "dedicated EM
pool, no shared decode" stance is **superseded** by the prefetch design, Q-R4/§5.)
