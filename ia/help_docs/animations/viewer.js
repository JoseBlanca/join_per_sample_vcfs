// Viewer for the pileup walker animation.
//
// Loads a scenario, runs the simulator to produce a trace, and
// renders one snapshot at a time. The user steps forward/backward
// through the trace; each step shows the state after one atomic
// algorithmic action.

const $ = (sel) => document.querySelector(sel);
const $$ = (sel) => [...document.querySelectorAll(sel)];

class Viewer {
  constructor() {
    this.scenarios = window.SCENARIOS;
    this.scenarioKey = Object.keys(this.scenarios)[0];
    this.scenario = null;
    this.trace = [];
    this.stepIdx = 0;
    this.playing = false;
    this.playTimer = null;
    this.speed = 700;

    this.populateScenarios();
    this.bindUI();
    this.loadScenario(this.scenarioKey);
  }

  populateScenarios() {
    const picker = $("#scenario-picker");
    for (const [key, scn] of Object.entries(this.scenarios)) {
      const opt = document.createElement("option");
      opt.value = key;
      opt.textContent = scn.name;
      picker.appendChild(opt);
    }
  }

  bindUI() {
    $("#scenario-picker").addEventListener("change", (e) =>
      this.loadScenario(e.target.value)
    );
    $("#prev").addEventListener("click", () => this.step(-1));
    $("#next").addEventListener("click", () => this.step(1));
    $("#reset").addEventListener("click", () => this.goto(0));
    $("#end").addEventListener("click", () => this.goto(this.trace.length - 1));
    $("#play").addEventListener("click", () => this.togglePlay());
    $("#speed").addEventListener("input", (e) => {
      this.speed = parseInt(e.target.value, 10);
      $("#speed-display").textContent = `${this.speed}ms`;
    });
    $("#speed-display").textContent = `${this.speed}ms`;
    document.addEventListener("keydown", (e) => {
      const tag = (e.target.tagName || "").toLowerCase();
      if (tag === "input" || tag === "select" || tag === "textarea") return;
      if (e.key === "ArrowRight") this.step(1);
      else if (e.key === "ArrowLeft") this.step(-1);
      else if (e.key === " ") {
        e.preventDefault();
        this.togglePlay();
      }
    });
  }

  loadScenario(key) {
    this.stopPlay();
    this.scenarioKey = key;
    $("#scenario-picker").value = key;

    const scn = this.scenarios[key];
    this.scenario = scn;
    this.refStart = scn.reference.start;
    this.refSeq = scn.reference.seq;

    const summary = scn.summary;
    const teaches = (scn.teaches || [])
      .map((t) => `<li>${this.escapeHtml(t)}</li>`)
      .join("");
    $("#scenario-description").innerHTML =
      `<p>${this.escapeHtml(summary)}</p>` +
      (teaches ? `<p><strong>What this scenario shows:</strong></p><ul>${teaches}</ul>` : "");

    const sim = new Simulator(scn);
    this.trace = sim.run();
    this.stepIdx = 0;
    this.render();
  }

  step(delta) {
    this.goto(this.stepIdx + delta);
  }

  goto(idx) {
    idx = Math.max(0, Math.min(idx, this.trace.length - 1));
    this.stepIdx = idx;
    this.render();
  }

  togglePlay() {
    if (this.playing) this.stopPlay();
    else this.startPlay();
  }

  startPlay() {
    if (this.stepIdx >= this.trace.length - 1) this.stepIdx = 0;
    this.playing = true;
    $("#play").textContent = "⏸";
    this.tick();
  }

  stopPlay() {
    this.playing = false;
    $("#play").textContent = "▶";
    if (this.playTimer) {
      clearTimeout(this.playTimer);
      this.playTimer = null;
    }
  }

  tick() {
    if (!this.playing) return;
    if (this.stepIdx >= this.trace.length - 1) {
      this.stopPlay();
      return;
    }
    this.step(1);
    this.playTimer = setTimeout(() => this.tick(), this.speed);
  }

  // -------------------------------------------------- rendering

  render() {
    const snap = this.trace[this.stepIdx];
    if (!snap) return;

    $("#step-counter").textContent = `${this.stepIdx + 1} / ${this.trace.length}`;
    $("#phase-name").textContent = snap.phase;
    $("#walker-display").textContent = snap.walker_pos;
    $("#code-ref").innerHTML =
      `<code>${this.escapeHtml(snap.code_ref)}</code>`;
    $("#description").innerHTML =
      `<p>${this.escapeHtml(snap.description)}</p>`;

    $$(".modules li").forEach((li) => {
      li.classList.toggle("active", li.dataset.module === snap.module);
    });

    this.renderRefAxis(snap);
    this.renderReads(snap);
    this.renderOpenRecords(snap);
    this.renderClosedRecords(snap);
  }

  renderRefAxis(snap) {
    const container = $("#ref-axis");
    container.innerHTML = "";

    const grid = document.createElement("div");
    grid.className = "pos-grid";
    grid.style.gridTemplateColumns = `repeat(${this.refSeq.length}, var(--cell-w))`;

    for (let i = 0; i < this.refSeq.length; i++) {
      const pos = this.refStart + i;
      const cell = document.createElement("div");
      cell.className = "ref-cell";
      if (pos === snap.walker_pos) cell.classList.add("walker-pos");
      cell.innerHTML =
        `<div class="pos-num">${pos}</div>` +
        `<div class="ref-base">${this.refSeq[i]}</div>`;
      grid.appendChild(cell);
    }
    container.appendChild(grid);

    // Walker pointer underneath the ref axis (separate row so it
    // doesn't overlap the bases).
    const ptrRow = document.createElement("div");
    ptrRow.className = "pos-grid walker-row";
    ptrRow.style.gridTemplateColumns = `repeat(${this.refSeq.length}, var(--cell-w))`;
    const col = snap.walker_pos - this.refStart + 1;
    if (col >= 1 && col <= this.refSeq.length) {
      const ptr = document.createElement("div");
      ptr.className = "walker-pointer";
      ptr.style.gridColumn = String(col);
      ptr.innerHTML = `<span class="walker-arrow">▲</span> walker_pos = ${snap.walker_pos}`;
      ptrRow.appendChild(ptr);
    } else {
      const note = document.createElement("div");
      note.className = "walker-pointer-out";
      note.style.gridColumn = `1 / -1`;
      note.textContent = `walker_pos = ${snap.walker_pos} (off the visible reference window)`;
      ptrRow.appendChild(note);
    }
    container.appendChild(ptrRow);
  }

  renderReads(snap) {
    const container = $("#reads");
    container.innerHTML = "";

    // Always render every read in the scenario — never just the active
    // ones. State (pending / active / expired) is shown visually so the
    // viewer follows each read through its full lifecycle.
    const reads = this.scenario.reads;
    if (reads.length === 0) {
      const e = document.createElement("div");
      e.className = "empty";
      e.textContent = "Scenario has no reads.";
      container.appendChild(e);
      return;
    }

    // Build a lookup from (qname, alignment_start) → active snapshot
    // entry; that's enough to disambiguate mate pairs since both mates
    // share a qname but have different alignment_starts.
    const activeByKey = new Map();
    for (const a of snap.active) {
      activeByKey.set(`${a.qname}|${a.alignment_start}`, a);
    }

    // Reads are admitted in coordinate order. The simulator records
    // pending_count per snapshot, so the first (total - pending_count)
    // reads have been admitted by this step. That distinguishes
    // "pending" (not admitted yet) from "expired" (admitted, then
    // dropped) without depending on walker_pos heuristics.
    const numAdmitted = reads.length - snap.pending_count;

    for (let i = 0; i < reads.length; i++) {
      const r = reads[i];
      const key = `${r.qname}|${r.alignment_start}`;
      const active = activeByKey.get(key);

      let state;
      if (active) state = "active";
      else if (i >= numAdmitted) state = "pending";
      else state = "expired";

      const wrapper = document.createElement("div");
      wrapper.className = `read-wrapper read-${state}`;

      // Header: name, strand, mate role, slot (only when active), CIGAR,
      // ref range, and a state badge.
      const strandSym = r.is_reverse_strand ? "← (rev)" : "(fwd) →";
      const mateSym = r.has_mate
        ? r.is_first_mate
          ? '<span class="mate-tag mate-1">first mate</span>'
          : '<span class="mate-tag mate-2">second mate</span>'
        : "";
      const slotLabel = active
        ? `<span class="read-slot">slot=${active.slot_id}</span>`
        : `<span class="read-slot read-slot-na">slot=—</span>`;
      const cigarStr = r.cigar.map((o) => `${o.len}${o.op}`).join("");
      const stateBadge =
        state === "active"
          ? '<span class="state-badge state-active">in active set</span>'
          : state === "pending"
          ? '<span class="state-badge state-pending">pending</span>'
          : '<span class="state-badge state-expired">expired</span>';

      const header = document.createElement("div");
      header.className = "read-header";
      header.innerHTML =
        `<span class="read-name">${this.escapeHtml(r.qname)}</span>` +
        `<span class="read-strand">${strandSym}</span>` +
        mateSym +
        slotLabel +
        `<span class="read-cigar">CIGAR: <code>${cigarStr}</code></span>` +
        `<span class="read-range">ref [${r.alignment_start}, ${r.alignment_end}]</span>` +
        stateBadge;
      wrapper.appendChild(header);

      const grid = document.createElement("div");
      grid.className = "read-row pos-grid";
      grid.style.gridTemplateColumns = `repeat(${this.refSeq.length}, var(--cell-w))`;

      // For pending/expired reads we still want to show what the read
      // looks like — same cell shape, just dimmed via the wrapper's
      // .read-{state} class.
      const readForCells = active
        ? active
        : { ...r, cigar_str: cigarStr, bq: r.bq, seq: r.seq };
      const cells = this.readCells(readForCells);
      for (const [pos, cell] of cells) {
        const c = document.createElement("div");
        c.className = `read-cell ${cell.cls}`;
        if (state === "active" && pos === snap.walker_pos) {
          c.classList.add("at-walker");
        }
        c.style.gridColumn = String(pos - this.refStart + 1);
        c.innerHTML = cell.html;
        grid.appendChild(c);
      }

      wrapper.appendChild(grid);
      container.appendChild(wrapper);
    }
  }

  readCells(r) {
    // Translate read.cigar_str + read.seq + read.bq into a list of
    // (ref_pos, {cls, html}) per ref column. M/=/X cells show the
    // base + bq; D cells show a dash; I shows as a sup marker on the
    // previous cell.
    const cells = [];
    const ops = this.parseCigar(r.cigar_str);
    let refPos = r.alignment_start;
    let readPos = 0;
    for (let i = 0; i < ops.length; i++) {
      const op = ops[i];
      if (op.op === "M" || op.op === "=" || op.op === "X") {
        for (let k = 0; k < op.len; k++) {
          const refBase = this.refSeq[refPos - this.refStart];
          const readBase = r.seq[readPos];
          const isSnp = readBase !== refBase && readBase !== "N";
          cells.push([
            refPos,
            {
              cls: isSnp ? "cell-snp" : readBase === "N" ? "cell-n" : "cell-match",
              html: `<span class="cell-base">${readBase}</span><span class="cell-bq">${r.bq[readPos]}</span>`,
            },
          ]);
          refPos++;
          readPos++;
        }
      } else if (op.op === "D") {
        for (let k = 0; k < op.len; k++) {
          cells.push([
            refPos,
            {
              cls: "cell-del",
              html: `<span class="cell-base">−</span>`,
            },
          ]);
          refPos++;
        }
      } else if (op.op === "I") {
        const prevPos = refPos - 1;
        const last = cells.find((c) => c[0] === prevPos);
        if (last) {
          last[1].html += `<sup class="ins-marker">+${op.len}</sup>`;
        }
        readPos += op.len;
      } else if (op.op === "S") {
        readPos += op.len;
      }
      // H, N don't appear in our scenarios; ignore.
    }
    return cells;
  }

  parseCigar(s) {
    const ops = [];
    let len = "";
    for (const ch of s) {
      if (/\d/.test(ch)) len += ch;
      else {
        ops.push({ op: ch, len: parseInt(len, 10) });
        len = "";
      }
    }
    return ops;
  }

  renderOpenRecords(snap) {
    const container = $("#open-records");
    container.innerHTML = "";

    if (snap.open_records.length === 0) {
      const e = document.createElement("div");
      e.className = "empty";
      e.textContent = "No open records.";
      container.appendChild(e);
      return;
    }

    const grid = document.createElement("div");
    grid.className = "pos-grid records-grid";
    grid.style.gridTemplateColumns = `repeat(${this.refSeq.length}, var(--cell-w))`;

    const highlightedPos = snap.extra ? snap.extra.highlightRecordPos : null;
    const widened = !!(snap.extra && snap.extra.widened);

    for (const rec of snap.open_records) {
      const card = document.createElement("div");
      card.className = "rec-card open";
      const startCol = rec.pos - this.refStart + 1;
      card.style.gridColumn = `${startCol} / span ${rec.ref_span}`;

      const isHigh = rec.pos === highlightedPos;
      if (isHigh) card.classList.add("highlighted");
      if (isHigh && widened) card.classList.add("widened");

      const allelesHtml = rec.alleles
        .map((a, idx) => {
          const tag = idx === 0 ? "REF" : `ALT${idx}`;
          const isAH =
            isHigh &&
            snap.extra &&
            snap.extra.highlightAlleleIdx === idx;
          return `
            <div class="allele${isAH ? " highlighted" : ""}">
              <span class="allele-tag">${tag}</span>
              <code class="allele-seq">${this.escapeHtml(a.seq)}</code>
              <span class="scalars">
                <span title="num_obs">obs=${a.num_obs}</span>
                <span title="q_sum (simplified Σ BQ)">q=${a.q_sum}</span>
                <span title="forward-strand count">fwd=${a.fwd}</span>
                <span title="placed_left">pl=${a.placed_left}</span>
                <span title="placed_start">ps=${a.placed_start}</span>
                <span title="contributing slot ids">slots={${a.slots.join(",")}}</span>
              </span>
            </div>
          `;
        })
        .join("");

      card.innerHTML =
        `<div class="rec-head">` +
        `<span class="rec-pos">pos=${rec.pos}</span>` +
        `<span class="rec-span">ref_span=${rec.ref_span}</span>` +
        `<span class="rec-close">closes when walker_pos ≥ ${rec.pos + rec.ref_span}</span>` +
        `</div>` +
        `<div class="rec-alleles">${allelesHtml}</div>`;

      grid.appendChild(card);
    }

    container.appendChild(grid);
  }

  renderClosedRecords(snap) {
    const container = $("#closed-records");
    container.innerHTML = "";

    if (snap.closed_records.length === 0) {
      const e = document.createElement("div");
      e.className = "empty";
      e.textContent = "No records emitted yet.";
      container.appendChild(e);
      return;
    }

    for (const rec of snap.closed_records) {
      const card = document.createElement("div");
      card.className = "rec-card closed";
      const refSeq = rec.alleles[0].seq;
      const variantBits = rec.alleles
        .map((a, idx) => {
          const tag = idx === 0 ? "REF" : `ALT${idx}`;
          return `<span class="closed-allele"><strong>${tag}</strong>=<code>${this.escapeHtml(a.seq)}</code> obs=${a.num_obs} q=${a.q_sum}</span>`;
        })
        .join(" ");
      card.innerHTML =
        `<span class="rec-head">pos=${rec.pos} (ref_span=${rec.ref_span})</span>` +
        `<span class="closed-alleles">${variantBits}</span>`;
      container.appendChild(card);
    }
  }

  escapeHtml(s) {
    return String(s)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  }
}

window.addEventListener("DOMContentLoaded", () => {
  window.viewer = new Viewer();
});
