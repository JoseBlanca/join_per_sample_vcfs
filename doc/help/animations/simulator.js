// Simulator for the pileup walker.
//
// This is a JavaScript reimplementation of the algorithm in
// src/per_sample_caller/pileup/. It is a simplified port — the
// production Rust code has more error handling, more edge cases
// (BAQ-adjusted BQ from upstream, column depth caps, slot exhaustion,
// chromosome boundaries, indel-overlap mate collapse, etc.), and
// owns its data more carefully. The simulator's job is to produce a
// faithful step-by-step trace for the three teaching scenarios; the
// architecture.md doc lists the simplifications.
//
// The simulator runs the algorithm phase by phase (admit → process →
// expire → close → advance) and snapshots state after every visible
// change, so the viewer can step through it.
//
// Each Rust submodule maps onto one block of methods below:
//
//   walker.rs           Simulator.run, phaseAdmit, phaseProcess,
//                       phaseExpire, phaseClose, phaseAdvance,
//                       resolveMateOverlap
//   active_set.rs       admit (inside phaseAdmit), expire (inside
//                       phaseExpire)
//   cigar_cursor.rs     buildCursor, cursorEventsAt,
//                       activeNextEventPos
//   open_record.rs      foldContributor, foldAtRecord (open + widen
//                       + fold), close (inside phaseClose)
//   slot_allocator.rs   allocateSlot, releaseSlot
//
// Where a step has a clear analogue in the Rust source, the snapshot's
// code_ref points to that file:line.

class Simulator {
  constructor(scenario) {
    this.scenario = scenario;
    this.refStart = scenario.reference.start;
    this.refSeq = scenario.reference.seq;
    this.refEnd = scenario.reference.start + scenario.reference.seq.length - 1;

    // Pending reads: shallow copies preserve scenario data.
    this.pending = scenario.reads.map(r => ({
      ...r,
      cigar: r.cigar.map(o => ({ ...o })),
      bq: [...r.bq],
    }));

    this.active = []; // [{ read_id, read, cursor, slot_id }]
    this.openRecords = new Map(); // pos -> open record
    this.closed = [];
    this.nextReadId = 0;

    // Slot allocator state — see slot_allocator.rs:21+
    this.slots = {
      nextFresh: 0,
      free: [],
      pendingFree: [],
      pendingMates: {}, // qname -> { slot, read_id }
      refcounts: {}, // slot_id -> count
    };

    // Walker starts at the first read's alignment_start so the first
    // step admits something rather than firing an empty advance.
    this.walkerPos = this.pending.length > 0
      ? this.pending[0].alignment_start
      : 1;

    this.trace = [];
    this.snapshot(
      "init",
      `Walker initialised. ${this.pending.length} read(s) pending. ` +
      `walker_pos starts at ${this.walkerPos} (the first read's alignment_start).`,
      "walker.rs:34",
      "walker"
    );
  }

  // ============================================================
  // Top-level loop                                  walker.rs:24
  // ============================================================

  run() {
    let safety = 0;
    while (this.pending.length > 0 || this.active.length > 0) {
      this.phaseAdmit();
      this.phaseProcess();
      this.phaseExpire();
      this.phaseClose();
      this.phaseAdvance();
      if (++safety > 500) {
        console.error("Simulator safety limit hit");
        break;
      }
    }
    this.snapshot(
      "done",
      `Walker finished. All reads consumed; all open records closed and emitted.`,
      "walker.rs:96",
      "walker"
    );
    return this.trace;
  }

  // ============================================================
  // Phase 1 — admit                              active_set.rs:89
  // ============================================================

  phaseAdmit() {
    while (
      this.pending.length > 0 &&
      this.pending[0].alignment_start <= this.walkerPos
    ) {
      const read = this.pending.shift();
      const readId = this.nextReadId++;
      const cursor = this.buildCursor(read);
      const allocated = this.allocateSlot(read, readId);
      const slotId = allocated.slot;
      this.active.push({ read_id: readId, read, cursor, slot_id: slotId });

      let slotMsg;
      if (allocated.reused) {
        slotMsg =
          `Reuse slot ${slotId} (first mate "${read.qname}" already in flight; ` +
          `slot refcount stays at 2).`;
      } else if (read.has_mate) {
        slotMsg =
          `Allocate slot ${slotId}; refcount pre-bumped to 2 anticipating second mate. ` +
          `Register "${read.qname}" in pending_mates.`;
      } else {
        slotMsg = `Allocate slot ${slotId} (refcount = 1, solo read).`;
      }

      this.snapshot(
        "admit",
        `Admit "${read.qname}" (alignment_start=${read.alignment_start} ≤ walker_pos=${this.walkerPos}). ` +
        `Build CIGAR cursor (precomputed offsets for ${read.cigar.length} op(s)). ${slotMsg}`,
        "active_set.rs:89",
        "active_set",
        { highlightReadId: readId }
      );
    }
  }

  buildCursor(read) {
    // Precompute (ref_pos, read_pos) at each op start, plus a sentinel
    // after the last op. See cigar_cursor.rs's `OpOffset` table.
    let refPos = read.alignment_start;
    let readPos = 0;
    const offsets = [];
    for (const op of read.cigar) {
      offsets.push({ ref_pos: refPos, read_pos: readPos });
      if ("MD=XN".includes(op.op)) refPos += op.len;
      if ("MIS=X".includes(op.op)) readPos += op.len;
    }
    offsets.push({ ref_pos: refPos, read_pos: readPos }); // sentinel
    return { ops: read.cigar, offsets };
  }

  allocateSlot(read, readId) {
    // Second mate? Reuse first mate's slot.
    if (read.has_mate && this.slots.pendingMates[read.qname] !== undefined) {
      const existing = this.slots.pendingMates[read.qname];
      delete this.slots.pendingMates[read.qname];
      return { slot: existing.slot, reused: true };
    }
    // Otherwise mint a fresh (or recycled) slot.
    const slot = this.slots.free.length > 0
      ? this.slots.free.shift()
      : this.slots.nextFresh++;
    if (read.has_mate) {
      this.slots.refcounts[slot] = 2;
      this.slots.pendingMates[read.qname] = { slot, read_id: readId };
    } else {
      this.slots.refcounts[slot] = 1;
    }
    return { slot, reused: false };
  }

  releaseSlot(slotId) {
    this.slots.refcounts[slotId] -= 1;
    if (this.slots.refcounts[slotId] === 0) {
      this.slots.pendingFree.push(slotId);
    }
  }

  // ============================================================
  // Phase 2 — process_position                walker.rs:200ish
  // ============================================================

  phaseProcess() {
    // 2a. Query each active cursor for events anchored at walker_pos.
    const contributors = [];
    for (const a of this.active) {
      const events = this.cursorEventsAt(a, this.walkerPos);
      if (events.length > 0) contributors.push({ active: a, events });
    }
    if (contributors.length === 0) return;

    const evtSummary = contributors
      .map(c =>
        `${c.active.read.qname} → ${c.events.map(e => this.eventStr(e)).join(", ")}`
      )
      .join("; ");

    this.snapshot(
      "cursor_query",
      `cigar_cursor.events_at(walker_pos=${this.walkerPos}) on each active read. ${evtSummary}.`,
      "cigar_cursor.rs:190",
      "cigar_cursor",
      { contributors: contributors.map(c => c.active.read_id) }
    );

    // 2b. Mate-overlap resolution.
    this.resolveMateOverlap(contributors);

    // 2c. Fold each contributor into its target record.
    for (const c of contributors) {
      this.foldContributor(c);
    }
  }

  cursorEventsAt(active, walkerPos) {
    const events = [];
    const { ops, offsets } = active.cursor;

    for (let i = 0; i < ops.length; i++) {
      const op = ops[i];
      const opRef = offsets[i].ref_pos;
      const opRead = offsets[i].read_pos;

      if (op.op === "M" || op.op === "=" || op.op === "X") {
        const local = walkerPos - opRef;
        if (local >= 0 && local < op.len) {
          const base = active.read.seq[opRead + local];
          // Read-N is silently skipped (F5 commit).
          if (base === "N") continue;
          // Adaptor boundary check (G1 commit).
          if (this.inAdaptor(active.read, walkerPos)) continue;
          events.push({
            type: "M",
            ref_pos: walkerPos,
            base,
            bq: active.read.bq[opRead + local],
            read_pos: opRead + local,
          });
        }
      } else if (op.op === "D") {
        // First-CIGAR indel is rejected upstream — skip if i === 0.
        if (i === 0) continue;
        const anchor = opRef - 1;
        if (anchor === walkerPos) {
          // BQ proxy: min over the BAQ window (l + 2). Simplified:
          // take the two flanking BAQ values.
          const flank = offsets[i].read_pos;
          const lo = Math.max(0, flank - 1);
          const hi = Math.min(active.read.bq.length - 1, flank);
          let bq_min = Infinity;
          for (let k = lo; k <= hi; k++) {
            bq_min = Math.min(bq_min, active.read.bq[k]);
          }
          events.push({
            type: "D",
            anchor,
            deleted_len: op.len,
            bq_proxy: bq_min,
          });
        }
      } else if (op.op === "I") {
        if (i === 0 || i === ops.length - 1) continue;
        const anchor = opRef - 1;
        if (anchor === walkerPos) {
          const insSeq = active.read.seq.slice(
            offsets[i].read_pos,
            offsets[i].read_pos + op.len
          );
          const lo = Math.max(0, offsets[i].read_pos - 1);
          const hi = Math.min(
            active.read.bq.length,
            offsets[i].read_pos + op.len + 1
          );
          let bq_min = Infinity;
          for (let k = lo; k < hi; k++) {
            bq_min = Math.min(bq_min, active.read.bq[k]);
          }
          events.push({
            type: "I",
            anchor,
            inserted_seq: insSeq,
            bq_proxy: bq_min,
          });
        }
      }
      // S, H, N produce no events at any walker_pos; skipped.
    }

    return events;
  }

  inAdaptor(read, refPos) {
    if (read.adaptor_boundary == null) return false;
    return read.is_reverse_strand
      ? refPos <= read.adaptor_boundary
      : refPos >= read.adaptor_boundary;
  }

  // ----- Mate-overlap resolution -------- walker.rs:394
  resolveMateOverlap(contributors) {
    // Group by slot_id; any slot with ≥2 contributors is a mate pair
    // overlapping at this position.
    const bySlot = new Map();
    for (const c of contributors) {
      if (!bySlot.has(c.active.slot_id)) bySlot.set(c.active.slot_id, []);
      bySlot.get(c.active.slot_id).push(c);
    }

    let resolved = false;
    let detail = "";

    for (const [slotId, members] of bySlot) {
      if (members.length < 2) continue;
      // Simplified scenarios: only the M+M case (no indel mate collapse).
      const m1 = members[0];
      const m2 = members[1];
      const e1 = m1.events.find(e => e.type === "M");
      const e2 = m2.events.find(e => e.type === "M");
      if (!e1 || !e2) continue;

      // Deterministic keeper: is_first_mate, then smaller alignment_start.
      let keeper, loser;
      if (m1.active.read.is_first_mate && !m2.active.read.is_first_mate) {
        keeper = m1; loser = m2;
      } else if (!m1.active.read.is_first_mate && m2.active.read.is_first_mate) {
        keeper = m2; loser = m1;
      } else if (m1.active.read.alignment_start <= m2.active.read.alignment_start) {
        keeper = m1; loser = m2;
      } else {
        keeper = m2; loser = m1;
      }

      const kE = keeper.events.find(e => e.type === "M");
      const lE = loser.events.find(e => e.type === "M");

      if (kE.base === lE.base) {
        const summed = Math.min(kE.bq + lE.bq, 200);
        const oldKBq = kE.bq;
        const oldLBq = lE.bq;
        kE.bq = summed;
        lE.bq = 0;
        detail =
          `Slot ${slotId}: "${keeper.active.read.qname}" and "${loser.active.read.qname}" agree on base ${kE.base}. ` +
          `Keep ${keeper.active.read.qname} with BQ = min(${oldKBq} + ${oldLBq}, 200) = ${summed}; ` +
          `zero ${loser.active.read.qname}'s BQ at this position.`;
      } else {
        // Disagreement: keep higher BQ (after the deterministic
        // tiebreak above), scale by 0.8.
        if (lE.bq > kE.bq) {
          // Keeper picks the higher-BQ side regardless of first-mate.
          [keeper, loser] = [loser, keeper];
        }
        const k = keeper.events.find(e => e.type === "M");
        const l = loser.events.find(e => e.type === "M");
        const oldKBq = k.bq;
        k.bq = Math.round(0.8 * k.bq);
        l.bq = 0;
        detail =
          `Slot ${slotId}: bases disagree (${k.base} vs ${l.base}). ` +
          `Keep higher-BQ side scaled to 0.8 × ${oldKBq} = ${k.bq}; zero the other side.`;
      }
      resolved = true;
    }

    if (resolved) {
      this.snapshot(
        "mate_overlap",
        `Mate-pair overlap detected at walker_pos=${this.walkerPos}. ${detail}`,
        "walker.rs:394",
        "walker",
        {}
      );
    }
  }

  // ----- Fold ---------------- open_record.rs:494
  foldContributor(c) {
    // events_at(walker_pos) returns events all anchored at walker_pos,
    // so they all target the same record.
    this.foldAtRecord(c.active, this.walkerPos, c.events);
  }

  foldAtRecord(active, anchorPos, events) {
    // 1. Open the record if not already open.
    let rec = this.openRecords.get(anchorPos);
    if (!rec) {
      const refBase = this.refBase(anchorPos);
      rec = {
        pos: anchorPos,
        ref_span: 1,
        alleles: [{
          seq: refBase,
          num_obs: 0, q_sum: 0, fwd: 0, placed_left: 0, placed_start: 0,
          slots: new Set(),
        }],
        folded_reads: new Map(),
      };
      this.openRecords.set(anchorPos, rec);
      this.snapshot(
        "fold_open",
        `Open new record at pos=${anchorPos}. Fetch REF base "${refBase}" from FASTA. ` +
        `alleles[0] (REF) initialised with 0 observations.`,
        "open_record.rs:288",
        "open_record",
        { highlightRecordPos: anchorPos }
      );
    }

    // 2. Determine required REF span from event footprints.
    let neededSpan = rec.ref_span;
    for (const e of events) {
      if (e.type === "D") neededSpan = Math.max(neededSpan, 1 + e.deleted_len);
      // M and I don't extend ref_span (I has zero ref-side footprint
      // beyond the anchor; M has footprint 1).
    }

    if (rec.ref_span < neededSpan) {
      const oldSpan = rec.ref_span;
      const oldRef = rec.alleles[0].seq;
      const newBases = this.refSlice(rec.pos + rec.ref_span, neededSpan - rec.ref_span);
      for (const a of rec.alleles) {
        a.seq = a.seq + newBases;
      }
      rec.ref_span = neededSpan;
      this.snapshot(
        "fold_widen",
        `Widen record at pos=${rec.pos}: ref_span ${oldSpan} → ${neededSpan} ` +
        `(needed by deletion footprint). Fetch "${newBases}" from FASTA, append to every allele's seq. ` +
        `REF: "${oldRef}" → "${rec.alleles[0].seq}".`,
        "open_record.rs:281",
        "open_record",
        { highlightRecordPos: rec.pos, widened: true }
      );
    }

    // 3. Build this contributor's haplotype string under [pos, pos+ref_span).
    let allZeroBQ = events.every(e => {
      if (e.type === "M") return e.bq === 0;
      return (e.bq_proxy ?? 0) === 0;
    });
    if (allZeroBQ) {
      this.snapshot(
        "fold_skip",
        `Skip ${active.read.qname}'s contribution at pos=${rec.pos}: all event BQs zeroed by mate overlap.`,
        "open_record.rs:494",
        "open_record",
        { highlightRecordPos: rec.pos }
      );
      return;
    }

    const dEvent = events.find(e => e.type === "D");
    const iEvent = events.find(e => e.type === "I");
    const mByPos = new Map();
    for (const e of events) if (e.type === "M") mByPos.set(e.ref_pos, e);

    let hap = "";
    for (let pos = rec.pos; pos < rec.pos + rec.ref_span; pos++) {
      // Skip positions swallowed by a deletion (anchor < pos ≤ anchor+len).
      if (dEvent && pos > dEvent.anchor && pos <= dEvent.anchor + dEvent.deleted_len) {
        continue;
      }
      const m = mByPos.get(pos);
      if (m) {
        hap += m.base;
      } else {
        // No event from this read at this ref position — read implies REF.
        // (Real algorithm uses events_overlapping; a position with no Match
        // event from this read is one the read doesn't observe directly.
        // For our scenarios, fall back to REF.)
        hap += this.refBase(pos);
      }
    }
    if (iEvent) {
      // Insertion: anchor base at position 0; insert after.
      hap = hap.charAt(0) + iEvent.inserted_seq + hap.slice(1);
    }

    // 4. Compute BQ contribution.
    let bqContribution = 0;
    for (const e of events) {
      if (e.type === "M") bqContribution += e.bq;
      else bqContribution += e.bq_proxy;
    }

    // 5. Find or create allele bucket.
    let alleleIdx = rec.alleles.findIndex(a => a.seq === hap);
    let alleleAdded = false;
    if (alleleIdx === -1) {
      rec.alleles.push({
        seq: hap,
        num_obs: 0, q_sum: 0, fwd: 0, placed_left: 0, placed_start: 0,
        slots: new Set(),
      });
      alleleIdx = rec.alleles.length - 1;
      alleleAdded = true;
    }

    // 6. Update scalars.
    const a = rec.alleles[alleleIdx];
    a.num_obs += 1;
    a.q_sum += bqContribution;
    if (!active.read.is_reverse_strand) a.fwd += 1;
    if (active.read.alignment_start < rec.pos) a.placed_left += 1;
    if (active.read.alignment_start === rec.pos) a.placed_start += 1;
    a.slots.add(active.slot_id);
    rec.folded_reads.set(active.read_id, alleleIdx);

    // 7. Describe what happened.
    const events_str = events.map(e => this.eventStr(e)).join(" + ");
    let desc;
    if (alleleAdded) {
      if (alleleIdx === 0) {
        desc =
          `Fold "${active.read.qname}" — ${events_str}. ` +
          `Haplotype "${hap}" matches REF (allele 0). ` +
          `num_obs ${a.num_obs - 1} → ${a.num_obs}.`;
      } else if (hap.length < rec.ref_span) {
        desc =
          `Fold "${active.read.qname}" — ${events_str}. ` +
          `Haplotype "${hap}" (length ${hap.length} < ref_span ${rec.ref_span}) → DELETION allele. ` +
          `Create new allele bucket idx=${alleleIdx}. num_obs=${a.num_obs}.`;
      } else if (hap.length > rec.ref_span) {
        desc =
          `Fold "${active.read.qname}" — ${events_str}. ` +
          `Haplotype "${hap}" (length ${hap.length} > ref_span ${rec.ref_span}) → INSERTION allele. ` +
          `Create new allele bucket idx=${alleleIdx}. num_obs=${a.num_obs}.`;
      } else {
        desc =
          `Fold "${active.read.qname}" — ${events_str}. ` +
          `Haplotype "${hap}" differs from REF "${rec.alleles[0].seq}" → SNP / MNP allele. ` +
          `Create new allele bucket idx=${alleleIdx}. num_obs=${a.num_obs}.`;
      }
    } else {
      desc =
        `Fold "${active.read.qname}" — ${events_str}. ` +
        `Haplotype "${hap}" matches existing allele idx=${alleleIdx}. ` +
        `num_obs ${a.num_obs - 1} → ${a.num_obs}.`;
    }

    this.snapshot(
      "fold",
      desc,
      "open_record.rs:494",
      "open_record",
      { highlightRecordPos: rec.pos, highlightAlleleIdx: alleleIdx }
    );
  }

  // ============================================================
  // Phase 3 — expire                          active_set.rs:134
  // ============================================================

  phaseExpire() {
    const expired = [];
    this.active = this.active.filter(a => {
      if (a.read.alignment_end < this.walkerPos) {
        expired.push(a);
        this.releaseSlot(a.slot_id);
        return false;
      }
      return true;
    });
    if (expired.length === 0) return;
    const summary = expired
      .map(a =>
        `"${a.read.qname}" (alignment_end=${a.read.alignment_end} < walker_pos=${this.walkerPos}; ` +
        `slot ${a.slot_id} refcount → ${this.slots.refcounts[a.slot_id]})`
      )
      .join(", ");
    this.snapshot(
      "expire",
      `Expire ${expired.length} read(s): ${summary}.`,
      "active_set.rs:134",
      "active_set",
      {}
    );
  }

  // ============================================================
  // Phase 4 — close                           open_record.rs:142
  // ============================================================

  phaseClose() {
    const sorted = [...this.openRecords.entries()].sort((x, y) => x[0] - y[0]);
    const closed = [];
    for (const [pos, rec] of sorted) {
      if (pos + rec.ref_span <= this.walkerPos) {
        const finalRec = {
          pos: rec.pos,
          ref_span: rec.ref_span,
          alleles: rec.alleles.map(a => ({
            seq: a.seq,
            num_obs: a.num_obs,
            q_sum: a.q_sum,
            fwd: a.fwd,
            placed_left: a.placed_left,
            placed_start: a.placed_start,
            slots: [...a.slots].sort((u, v) => u - v),
          })),
        };
        this.closed.push(finalRec);
        this.openRecords.delete(pos);
        closed.push(finalRec);
      } else {
        // BTreeMap iteration is in-order; once we hit a non-closeable
        // record, no later record can satisfy the rule earlier.
        break;
      }
    }
    if (closed.length === 0) return;
    const summary = closed
      .map(r =>
        `pos=${r.pos} ref_span=${r.ref_span} (${r.pos}+${r.ref_span} ≤ ${this.walkerPos}); ` +
        `${r.alleles.length} allele(s)`
      )
      .join("; ");
    this.snapshot(
      "close",
      `Close ${closed.length} record(s): ${summary}. Send each through SyncSender to the encoder thread.`,
      "open_record.rs:142",
      "open_record",
      {}
    );
  }

  // ============================================================
  // Phase 5 — advance                            walker.rs:88
  // ============================================================

  phaseAdvance() {
    const nextActive = this.activeNextEventPos();
    const nextPending = this.pending.length > 0
      ? this.pending[0].alignment_start
      : Infinity;
    let nextPos;
    if (nextActive !== null && nextActive <= nextPending) {
      nextPos = nextActive;
    } else if (nextPending !== Infinity) {
      nextPos = nextPending;
    } else {
      // No future events from active reads, no pending reads. We may
      // still have records to close or reads to expire — advance by 1.
      nextPos = this.walkerPos + 1;
    }
    if (nextPos <= this.walkerPos) nextPos = this.walkerPos + 1;
    const oldPos = this.walkerPos;
    this.walkerPos = nextPos;
    const isJump = nextPos !== oldPos + 1;
    this.snapshot(
      "advance",
      isJump
        ? `Advance walker_pos: ${oldPos} → ${this.walkerPos} (jump; no active read has events between).`
        : `Advance walker_pos: ${oldPos} → ${this.walkerPos}.`,
      "walker.rs:88",
      "walker",
      {}
    );
  }

  activeNextEventPos() {
    let min = null;
    for (const a of this.active) {
      const { ops, offsets } = a.cursor;
      for (let i = 0; i < ops.length; i++) {
        const op = ops[i];
        const opRef = offsets[i].ref_pos;
        if (op.op === "M" || op.op === "=" || op.op === "X") {
          for (let k = 0; k < op.len; k++) {
            const p = opRef + k;
            if (p > this.walkerPos && (min === null || p < min)) min = p;
          }
        } else if (op.op === "D" && i > 0) {
          const anchor = opRef - 1;
          if (anchor > this.walkerPos && (min === null || anchor < min)) {
            min = anchor;
          }
        } else if (op.op === "I" && i > 0 && i < ops.length - 1) {
          const anchor = opRef - 1;
          if (anchor > this.walkerPos && (min === null || anchor < min)) {
            min = anchor;
          }
        }
      }
    }
    return min;
  }

  // ============================================================
  // Helpers
  // ============================================================

  refBase(pos) {
    const i = pos - this.refStart;
    return i >= 0 && i < this.refSeq.length ? this.refSeq[i] : "?";
  }

  refSlice(start, len) {
    const i = start - this.refStart;
    return this.refSeq.slice(i, i + len);
  }

  eventStr(e) {
    if (e.type === "M") return `Match{base=${e.base}, bq=${e.bq}}`;
    if (e.type === "D") return `Del{anchor=${e.anchor}, deleted_len=${e.deleted_len}, bq=${e.bq_proxy}}`;
    if (e.type === "I") return `Ins{anchor=${e.anchor}, seq="${e.inserted_seq}", bq=${e.bq_proxy}}`;
    return e.type;
  }

  // ============================================================
  // Snapshot recorder
  // ============================================================

  snapshot(phase, description, codeRef, module, extra = {}) {
    this.trace.push({
      phase,
      walker_pos: this.walkerPos,
      pending_count: this.pending.length,
      active: this.active.map(a => ({
        read_id: a.read_id,
        qname: a.read.qname,
        alignment_start: a.read.alignment_start,
        alignment_end: a.read.alignment_end,
        cigar_str: a.read.cigar.map(o => `${o.len}${o.op}`).join(""),
        seq: a.read.seq,
        bq: [...a.read.bq],
        slot_id: a.slot_id,
        is_reverse: a.read.is_reverse_strand,
        is_first_mate: a.read.is_first_mate,
        has_mate: a.read.has_mate,
      })),
      open_records: [...this.openRecords.entries()]
        .sort((x, y) => x[0] - y[0])
        .map(([_, r]) => ({
          pos: r.pos,
          ref_span: r.ref_span,
          alleles: r.alleles.map(al => ({
            seq: al.seq,
            num_obs: al.num_obs,
            q_sum: al.q_sum,
            fwd: al.fwd,
            placed_left: al.placed_left,
            placed_start: al.placed_start,
            slots: [...al.slots].sort((u, v) => u - v),
          })),
        })),
      closed_records: this.closed.map(r => ({
        pos: r.pos,
        ref_span: r.ref_span,
        alleles: r.alleles.map(a => ({ ...a, slots: [...a.slots] })),
      })),
      slot_state: {
        next_fresh: this.slots.nextFresh,
        free: [...this.slots.free],
        pending_free: [...this.slots.pendingFree],
        pending_mates: Object.keys(this.slots.pendingMates),
        refcounts: { ...this.slots.refcounts },
      },
      description,
      code_ref: codeRef,
      module,
      extra,
    });
  }
}

window.Simulator = Simulator;
