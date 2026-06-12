# Stage 0 — `ssr-catalog` implementation sketch

**Status:** sketch for alignment, 2026-06-12. Files + main structs + main
function signatures — the *shape* we'll implement, not full pseudocode. Built on
the architecture ([ssr_catalog.md](../architecture/ssr_catalog.md)), the shared
types ([ssr_shared_types.md](../architecture/ssr_shared_types.md)), and the spec
(§3). We iterate on this until we agree, then code it.

Grounding facts already checked in the tree:
- **trf-mod runs per contig via temp files** (not pipes — §3): it prints BED to
  **stdout even with a file input** (`trf_print_bed(stdout)`), so we give it a
  temp input FASTA and redirect its stdout to a temp `.bed`. No stdin/stdout
  pipes ⇒ no deadlock, clean `wait()` + exit-status check.
- `tempfile = "3"` is already a dep, and
  [cram_files.rs](../../src/pileup/per_sample/cram_files.rs) `build_fasta` is a
  precedent: per-task `tempfile::tempdir()` + write FASTA, auto-cleaned on drop.
- We have `noodles-bgzf`, `noodles-csi`, `noodles-fasta`. FASTA side exposes
  `ContigList` + a `MultiChromRefFetcher` trait ([src/fasta/](../../src/fasta/)).
- `Locus` / `Motif` / `ContigId` come from `ssr::types` (shared types doc §5).

---

## 1. File layout

```
src/ssr/catalog/
├── mod.rs          # CatalogArgs, run() orchestrator, CatalogError
├── trf.rs          # locate + spawn trf-mod, parse its BED → TrfRecord
├── postprocess.rs  # TrfRecord[] + contig bytes → Vec<Locus> (the §4 pipeline)
└── io.rs           # CatalogHeader, CatalogWriter (+ CatalogReader for Stage 1/2)
```
Subcommand wiring: a `CatalogArgs` (clap) dispatched from the existing
subcommand router (beside `pileup`/`var-calling`).

---

## 2. `mod.rs` — orchestrator

```rust
pub struct CatalogArgs {
    reference: PathBuf,
    output: PathBuf,            // catalog .bed.gz (index written beside it)
    num_chroms_in_parallel: usize,  // §8.3 — DOP cap AND a RAM knob
    flank_bp: u32,              // §5 — placeholder default; pinned at Stage 1
    trf_mod_path: Option<PathBuf>,  // §2.4 discovery override
    temp_dir: PathBuf,          // §3 — root for per-contig trf-mod temp files;
                                //   default CWD-relative + DISK-backed (NOT /tmp/tmpfs)
    min_purity: f32,            // §4 purity floor (accuracy knob)
    min_score: i32,             // §4 early TRF-score accept-gate
}

pub fn run(args: CatalogArgs) -> Result<(), CatalogError>;
```

`run()` flow (the fan-out/collect of §8):

```
1. trf = trf::locate_trf_mod(args.trf_mod_path)?          // §2.4
2. contigs: ContigList   = read from reference (noodles-fasta)
   reference_md5: String = compute (reuse pop_var_caller::common md5)
3. header = CatalogHeader { reference, reference_md5,
                            trf_mod_version: trf::version(&trf)?,
                            params, flank_bp, tool_version, date }
4. per-contig build, IN PARALLEL, ORDER-PRESERVING:
     contigs.par_iter()                     // pool sized to num_chroms_in_parallel
        .map(|c| {
            let seq = load full contig bytes (held resident this task only)
            let recs = trf::run_on_contig(&trf, c.name(), &seq, &args.temp_dir)?;
            Ok(postprocess::build_loci(recs, c.id(), &seq, &params))  // Vec<Locus>, start-sorted
        })
        .collect::<Result<Vec<Vec<Locus>>>>()?   // index order == contig order  ⇒ ordered collector for free
5. let mut w = io::CatalogWriter::new(bgzf(args.output), header)?;
   for per_contig in results { for locus in per_contig { w.write_locus(&locus)?; } }
   w.finish()?;                              // flush + finalize bgzf
6. io::write_index(&args.output)?;   // bgzip is sorted ⇒ csi/tabix index
```

- **Ordered collector = `par_iter().collect::<Vec<_>>()`** — rayon preserves
  input order, so writing the Vec in sequence *is* contig order (the bgzip/tabix
  sort requirement, §8.2). No manual reorder buffer.
- **Determinism** (§8.4): the collected Vec is order-invariant to pool size →
  byte-identical catalog regardless of `--num-chroms-in-parallel`.
- **Memory** (§8.3): the collected `Vec<Vec<Locus>>` is the (small) row buffer;
  the large contig `seq` is dropped when each task returns. Concurrent resident
  sequences ≈ `num_chroms_in_parallel`.
- A custom rayon `ThreadPool` of `num_chroms_in_parallel` threads bounds
  concurrent `trf-mod` processes.

```rust
pub enum CatalogError { TrfModNotFound, TrfSpawn(std::io::Error), TrfParse{line, reason},
                        Fasta(..), Io(std::io::Error), Md5Mismatch{..} }
```
*(Note: with a sibling module named `io`, spell std I/O as `std::io::Error`
explicitly — the one small tax of the `io.rs` name.)*

---

## 3. `trf.rs` — locate, run, parse

```rust
/// §2.4 layered discovery: override → sibling of our exe → PATH → error.
pub fn locate_trf_mod(override_path: Option<&Path>) -> Result<PathBuf, CatalogError>;

/// `trf-mod -v` → version string for the catalog header (pinned + recorded).
pub fn version(bin: &Path) -> Result<String, CatalogError>;

/// One parsed BED row. Coords already 0-based half-open (trf_print_bed emits
/// first-1 .. last). We keep only what the catalog needs + sanity fields.
pub struct TrfRecord {
    start: u32, end: u32,     // 0-based half-open, from the BED
    period: u16,              // authoritative period (may differ from pattern len!)
    frac_match: f32,          // already in [0,1] — sanity-check only (purity recomputed §4)
    score: i32,               // early accept-gate (§4)
    pattern: Box<[u8]>,       // TRF consensus — sanity-check; motif comes from the tract (§4)
    // copy_num / frac_gap / entropy parsed-and-dropped (not needed)
}

/// Run trf-mod on one contig via temp files under `temp_root`, parse the BED.
/// Per contig — called from a parallel task, so collision-free by construction.
pub fn run_on_contig(bin: &Path, name: &str, seq: &[u8], temp_root: &Path)
    -> Result<Vec<TrfRecord>, CatalogError>;

/// Split one BED line on '\t'; assert 10 columns; parse; validate col-1 ctg == name.
fn parse_bed_line(line: &[u8], expect_ctg: &str) -> Result<TrfRecord, CatalogError>;
```

`run_on_contig` mechanics (temp files, no pipes):

1. **Unique per-task dir** — `tempfile::Builder::new().prefix("ssr-catalog-")
   .tempdir_in(temp_root)?`. The random dir name is what makes parallel tasks
   collision-free (the Rust equivalent of Python's `NamedTemporaryFile`); fixed
   `input.fa` / `output.bed` names *inside* it are then fine.
2. Write the contig to `input.fa` (`>{name}\n{seq}\n`).
3. `Command::new(bin).arg("input.fa").current_dir(&task_dir)
   .stdout(File::create("output.bed"))` → spawn → **`child.wait()` + assert
   `status.success()`** (this is the clean termination control we wanted; a
   crashed/killed trf-mod surfaces as a per-contig `CatalogError`, not a hang).
   `current_dir(task_dir)` corrals any stray files trf writes into the disposable
   dir, not the CWD.
4. Parse `output.bed` (skip any banner / "No TRs found" noise lines).
5. The `TempDir` drops at end of the call → `input.fa` + `output.bed` removed.

- **No stdin/stdout pipes** ⇒ no deadlock, and termination is an explicit
  exit-status check. **`temp_root` is CWD-relative + disk-backed by default**, so
  the big `input.fa` does not land on RAM-backed `/tmp` (tmpfs).
- BED column order is **source-confirmed** (arch §3); the parser hard-asserts the
  field count/types and fails loudly, so an upstream format change is caught.

---

## 4. `postprocess.rs` — the §4 pipeline (per contig)

**We DROP compounds + bundles (no split — GangSTR-style, arch §4).** That kills
the change-point problem *and* the inner-flank case: with `bundle_threshold ≥
flank_bp`, every surviving locus is isolated, so flanks are clean by construction.
No `SubLocus`, no overlap recording.

```rust
pub struct PostProcessParams {
    min_purity: f32, min_score: i32, flank_bp: u32,
    bundle_threshold: u32,   // ≥ flank_bp ⇒ surviving loci have clean flanks
    // per-period copy-number floors {1:10,2:5,3:4,4:3,5:3,6:3} (GangSTR minimal_trim)
}

/// Order: period≤6 → drop compound-motif → drop bundles → end-trim (+copy floor)
/// → recompute purity & filter → embed ref_seq. Start-sorted Locus per contig.
pub fn build_loci(recs: Vec<TrfRecord>, contig: ContigId, contig_seq: &[u8],
                  p: &PostProcessParams) -> Vec<Locus>;
```

Internal steps (each a small fn, in order):

```rust
// 1. filter raw TRF calls: drop period > 6 (+ an early TRF score floor).
fn filter_trf_calls(recs, p) -> Vec<TrfRecord>;

// 2. drop loci whose motif is internally periodic (ATAT=(AT)²) —
//    GangSTR minimal_trim::is_compound. Fundamental period survives via TRF dedup.
fn drop_compound_motifs(recs: Vec<TrfRecord>) -> Vec<TrfRecord>;

// 3. drop bundles: any locus within `bundle_threshold` bp of another → drop the
//    whole cluster (GangSTR remove_bundles). Streaming over start-sorted recs.
fn drop_bundles(recs: Vec<TrfRecord>, bundle_threshold: u32) -> Vec<TrfRecord>;

// 4. end-trim partial motifs to clean whole-motif boundaries (GangSTR minimal_trim);
//    then apply per-period copy-number floor. Adjusts start/end.
fn end_trim(rec: TrfRecord, contig_seq, p) -> Option<TrfRecord>;

// 5. motif = first `period` bases of the (trimmed) tract (verbatim/phase-faithful,
//    types §5; len(motif)==period despite TRF pattern-vs-period drift, arch §3).
fn motif_from_tract(contig_seq, start, period) -> Motif;

// 6. recompute purity from the trimmed tract vs a perfect motif tiling — spec §3.2
//    definition exactly (NOT TRF frac_match). Then apply the purity floor.
fn recompute_purity(tract: &[u8], motif: &Motif) -> f32;

// 7. embed ref_seq: tract + flank_bp each side, clamped at contig ends,
//    upper-cased. Returns (ref_bytes, ref_bytes_start) for the Locus.
fn embed_ref_seq(contig_seq, start, end, flank_bp) -> (Box<[u8]>, u32);
```

Each surviving record → a `Locus { chrom, start, end, motif, purity_fraction,
ref_bytes, ref_bytes_start }`. **Imperfect single-motif loci are kept** (the floor
is a degeneracy cutoff, not perfect-only — our one divergence from GangSTR's
perfect-only `remove_messy`).

---

## 5. `io.rs` — header, writer (+ reader for Stage 1/2)

```rust
pub struct CatalogHeader {           // the ## metadata block
    reference: PathBuf, reference_md5: String,
    trf_mod_version: String, params: PostProcessParams, flank_bp: u32,
    tool_version: String, date: String,
}

/// bgzip-wrapped TSV writer. Writes ## header + '#'-column header on new(),
/// then one tab row per locus. Rows arrive already (contig, start)-sorted.
pub struct CatalogWriter<W: Write> { inner: bgzf::Writer<W>, ... }
impl<W: Write> CatalogWriter<W> {
    pub fn new(sink: W, header: CatalogHeader) -> Result<Self, CatalogError>;
    pub fn write_locus(&mut self, locus: &Locus) -> Result<(), CatalogError>;
    pub fn finish(self) -> Result<(), CatalogError>;
}

/// Build the coordinate index from the finished bgzip TSV.
pub fn write_index(catalog_path: &Path) -> Result<(), CatalogError>;

// Serialization (one place, shared by writer + reader):
fn locus_to_row(&Locus) -> String;          // chrom start end motif purity_fraction ref_seq_start ref_seq
fn row_to_locus(&str) -> Result<Locus>;
```

**Read side (consumed by Stage 1/2 — sketched here for the round-trip; detailed
design belongs to those passes):**
```rust
pub struct CatalogReader { ... }
impl CatalogReader {
    pub fn open(path: &Path) -> Result<Self>;
    pub fn header(&self) -> &CatalogHeader;           // for md5 / param checks
    pub fn iter(&mut self) -> impl Iterator<Item = Result<Locus>>;     // whole catalog
    pub fn query(&mut self, chrom, start, end) -> impl Iterator<Item = Result<Locus>>; // --regions
}
```

> **To settle in code (index):** tabix `.tbi` needs adding `noodles-tabix`;
> `.csi` is buildable from `noodles-csi` we already depend on. Both are
> region-queryable. Lean **`.csi`** (no new dep) unless `.tbi` interop is wanted.

---

## 6. What this sketch deliberately leaves vague

- **Exact parser field handling** — pinned by the source layout + golden tests
  (arch §3), written when coding.
- **`flank_bp` / `bundle_threshold` values** — placeholder defaults here
  (`bundle_threshold ≥ flank_bp` is the invariant); `flank_bp` pinned at Stage 1.
- **`is_compound` / `end_trim` thresholds** — port the GangSTR `minimal_trim.py`
  constants (0.8 compound threshold; copy-number floors), confirm against its
  behaviour when coding.
- **CLI flag names/defaults** beyond the core set.

*(Compound-splitting — the previously-vague piece — is **gone**: we drop, not
split, per arch §4.)*

---

## 7. Round-trip & tests (the safety net we build alongside)

- **TRF-mod golden:** run pinned `trf-mod` on `TRF-mod/t/*.fasta`, snapshot the
  BED, assert `parse_bed_line` extracts the expected `TrfRecord`s.
- **Locus round-trip:** `row_to_locus(locus_to_row(l)) == l`.
- **Determinism:** same catalog bytes across `--num-chroms-in-parallel ∈ {1, N}`.
- **Pipeline units:** `recompute_purity` on perfect/interrupted tracts;
  `drop_bundles` on a synthetic cluster (and that survivors keep clean flanks);
  `drop_compound_motifs` on `ATAT`; `end_trim` boundary cleanup; `embed_ref_seq`
  clamping at a contig end; period≤6 filtering.
- Catalog-accuracy harness (randomized-seq FP, tomato markers) is spec §3.4 —
  separate, not part of this module's unit tests.
