//! S2 — the ephemeral per-locus spill file.
//!
//! When the paralog filter is on, the main caller pass cannot write the VCF
//! directly (a locus's keep/drop verdict needs global quantities not known
//! until the whole genome has streamed — arch §1). Instead each scored locus is
//! **spilled** to a temporary file, read back once to calibrate (S4) and once
//! more to write the surviving calls (S5). This module is that store: a
//! self-contained, length-framed binary stream — written once, read twice,
//! deleted when the run ends (success *or* failure).
//!
//! Framing (settled 2026-07-01 over the columnar `.psp` container): the spill's
//! payload is dominated by **per-sample** vectors (GT/GQ/AD per sample) which the
//! `.psp` container's per-record / per-allele column model does not fit, so the
//! spill is its own simple format — a `u32` little-endian length prefix followed
//! by one record's bytes, `BufWriter`/`BufReader`-buffered. Each record carries
//! the **whole** [`PosteriorRecord`] verbatim, so the write pass reconstructs a
//! bit-identical record and feeds it back through the unchanged VCF writer — the
//! byte-identity safety net holds by construction, with no dependence on which
//! record fields the writer happens to read.
//!
//! The paralog score's coverage input (per-sample centred-window GC + mean
//! coverage) rides **on** each record spill entry ([`ParalogSpillRecord::window_coverage`]),
//! gathered by the caller from the psp windowed columns. The sibling window-spill
//! stream (Approach A, S6c) that once carried it, joined by tile key, is being
//! retired (its read side is gone; the write side follows in M6).
//!
//! Memory-flat by construction: the writer encodes into one reused scratch
//! buffer and the reader decodes one record at a time into a reused byte buffer
//! — neither retains a genome-wide structure (arch §9).

use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

use tempfile::{NamedTempFile, TempPath};

use crate::pileup_record::AlleleSupportStats;
use crate::var_calling::per_group_merger::{CompoundConstituent, MergedAllele};
use crate::var_calling::posterior_engine::{EmDiagnostics, PosteriorRecord, RecordLocus};
use crate::var_calling::types::LocusWindowCoverage;

/// One locus on the **record spill**: the full caller output, carried verbatim,
/// plus the locus's stored paralog likelihood ratio.
///
/// The [`PosteriorRecord`] is stored in full so the write pass (S5) can
/// reconstruct it exactly and emit byte-identical VCF.
///
/// `paralog_lr` is the locus's H1-vs-H2 likelihood ratio, computed **once** and
/// stored here so the write pass reads it instead of re-scoring (the
/// inline-scoring architecture,
/// `doc/devel/architecture/hidden_paralog_inline_scoring.md`). `f64::NAN` marks a
/// locus that was not scored (not a biallelic SNP, no usable samples, or no
/// window) — which is *kept*, matching the existing "non-finite LR is never
/// flagged" rule. (Transitional: until inline scoring lands, the sink writes
/// `NAN` here and the calibrate/write passes still recompute the LR.)
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ParalogSpillRecord {
    /// The caller's emitted record for this locus, carried in full.
    pub record: PosteriorRecord,
    /// The stored paralog likelihood ratio (`f64::NAN` = unscored → kept).
    pub paralog_lr: f64,
    /// Per-sample centred-window GC + coverage at this locus, gathered from the
    /// psp windowed columns by the caller ([`CalledChunk::window_coverage`]).
    /// The score reads it directly — no window join. `NaN` per sample = no
    /// covered record there (that sample is skipped). Empty (both vectors) for a
    /// locus whose chunk carried no window coverage (e.g. the write pass's
    /// re-emitted survivors, which never re-score).
    pub window_coverage: LocusWindowCoverage,
}

/// One analysis window on the **window spill** (Approach A, S6c): the window's
/// tile key, its shared reference GC, and each sample's mean depth over the
/// window. The calibrate / write passes join this to each locus record by tile
/// key, so the paralog score's coverage inputs are supplied without carrying
/// per-sample window data in the (per-locus, all-samples) record spill.
///
/// Written in coordinate (tile) order by the producer's window builder; read
/// back in lockstep with the record spill (both coordinate-ordered).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct WindowSpillRecord {
    /// Chromosome id (matches the record spill's `locus.chrom_id`).
    pub chrom_id: u32,
    /// Window tile index, `(pos − 1) / window_bp`.
    pub tile: u32,
    /// The window's shared reference GC fraction.
    pub gc: f32,
    /// Per cohort sample: mean depth over the window, or `None` if the sample
    /// had no covered position in it.
    pub depths: Vec<Option<f32>>,
}

/// Failure modes of the spill read/write path.
#[derive(Debug, thiserror::Error)]
pub enum SpillError {
    /// Underlying file / stream I/O failed.
    #[error("spill I/O error")]
    Io(#[from] io::Error),
    /// A record's encoded payload exceeded the `u32` length frame. Per-locus
    /// records are far under 4 GiB, so this only fires on corruption / a bug.
    #[error("spill record length {len} exceeds the u32 frame limit")]
    RecordTooLarge { len: usize },
    /// A record's bytes ended before all fields were decoded — a truncated
    /// spill (partial final write) or a framing bug.
    #[error("spill record truncated: needed {needed} more byte(s) decoding {field}")]
    Truncated { field: &'static str, needed: usize },
    /// A decoded field held a value outside its representable range (e.g. a
    /// length that overflows `usize` on this target).
    #[error("spill record field {field} out of range")]
    OutOfRange { field: &'static str },
}

/// An owned, ephemeral spill file that deletes itself on drop.
///
/// Built on [`tempfile::TempPath`], so the backing file is removed when this
/// value is dropped — on the normal path, on an early `?` return, and on an
/// unwinding panic alike (arch §5: deleted on success *and* failure). Create it
/// in a project-local scratch directory, write it once via [`writer`], then
/// read it back (twice) via [`reader`].
///
/// [`writer`]: ParalogSpill::writer
/// [`reader`]: ParalogSpill::reader
pub(crate) struct ParalogSpill {
    /// The temp file's path; its `Drop` unlinks the file.
    temp_path: TempPath,
}

impl ParalogSpill {
    /// Create a fresh, empty spill file inside `dir` (a project-local scratch
    /// directory — never the system temp per the project's scratch policy).
    pub(crate) fn create_in(dir: &Path) -> Result<Self, SpillError> {
        let temp_path = NamedTempFile::new_in(dir)?.into_temp_path();
        Ok(Self { temp_path })
    }

    /// The backing file's path (for diagnostics / size reporting).
    pub(crate) fn path(&self) -> &Path {
        &self.temp_path
    }

    /// Open a buffered writer positioned at the start, truncating any previous
    /// contents. Call once; the caller streams records through
    /// [`ParalogSpillWriter::append`] and then [`ParalogSpillWriter::finish`].
    pub(crate) fn writer(&self) -> Result<ParalogSpillWriter<BufWriter<File>>, SpillError> {
        let file = OpenOptions::new()
            .write(true)
            .truncate(true)
            .open(&self.temp_path)?;
        Ok(ParalogSpillWriter::new(BufWriter::with_capacity(
            64 * 1024,
            file,
        )))
    }

    /// Open a buffered one-pass reader from the start. May be called more than
    /// once (the calibrate pass, then the write pass) — each returns an
    /// independent cursor over the same file.
    pub(crate) fn reader(&self) -> Result<ParalogSpillReader<BufReader<File>>, SpillError> {
        let file = File::open(&self.temp_path)?;
        Ok(ParalogSpillReader::new(BufReader::with_capacity(
            64 * 1024,
            file,
        )))
    }

    /// Open a buffered [`WindowSpillWriter`] over this file (truncating). Used
    /// when the spill holds [`WindowSpillRecord`]s (the window spill).
    pub(crate) fn window_writer(&self) -> Result<WindowSpillWriter<BufWriter<File>>, SpillError> {
        let file = OpenOptions::new()
            .write(true)
            .truncate(true)
            .open(&self.temp_path)?;
        Ok(WindowSpillWriter::new(BufWriter::with_capacity(
            64 * 1024,
            file,
        )))
    }
}

/// Streaming append-only writer over any sink. Encodes each record into one
/// reused scratch buffer, then emits a `u32` length frame followed by the
/// payload — so writing holds exactly one record's bytes at a time.
pub(crate) struct ParalogSpillWriter<W: Write> {
    sink: W,
    /// Reused per-record encode buffer (grows to the largest record, then
    /// stays — one record's worth, never the whole stream).
    scratch: Vec<u8>,
    records_written: u64,
}

impl<W: Write> ParalogSpillWriter<W> {
    pub(crate) fn new(sink: W) -> Self {
        Self {
            sink,
            scratch: Vec::new(),
            records_written: 0,
        }
    }

    /// Append one locus record: encode into the scratch buffer, then write its
    /// length frame + payload to the sink.
    pub(crate) fn append(&mut self, record: &ParalogSpillRecord) -> Result<(), SpillError> {
        self.scratch.clear();
        encode_record(record, &mut self.scratch);
        write_frame(&mut self.sink, &self.scratch)?;
        self.records_written += 1;
        Ok(())
    }

    /// Records appended so far. (Exercised by the round-trip tests and a
    /// diagnostic accessor for callers; allowed so the non-test lib build
    /// doesn't flag it before a production caller reads it.)
    #[allow(dead_code)]
    pub(crate) fn records_written(&self) -> u64 {
        self.records_written
    }

    /// Flush the sink and return it (the caller then drops or fsyncs it).
    pub(crate) fn finish(mut self) -> Result<W, SpillError> {
        self.sink.flush()?;
        Ok(self.sink)
    }
}

/// One-pass reader that yields records from a spill stream, holding exactly one
/// record's bytes at a time in a reused buffer.
pub(crate) struct ParalogSpillReader<R: Read> {
    source: R,
    /// Reused per-record decode buffer (sized to the current record).
    buf: Vec<u8>,
}

impl<R: Read> ParalogSpillReader<R> {
    pub(crate) fn new(source: R) -> Self {
        Self {
            source,
            buf: Vec::new(),
        }
    }

    /// Read the next record, or `None` at a clean end of stream. A stream that
    /// ends *inside* a frame (partial length prefix or short payload) is a
    /// [`SpillError::Truncated`]/[`SpillError::Io`], not a clean end.
    pub(crate) fn next_record(&mut self) -> Option<Result<ParalogSpillRecord, SpillError>> {
        match read_frame(&mut self.source, &mut self.buf) {
            Ok(false) => return None, // clean EOF at a record boundary
            Ok(true) => {}
            Err(e) => return Some(Err(e)),
        }
        Some(decode_record(&self.buf))
    }
}

impl<R: Read> Iterator for ParalogSpillReader<R> {
    type Item = Result<ParalogSpillRecord, SpillError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.next_record()
    }
}

/// Streaming append-only writer for the window spill (same length-framing as the
/// record spill, one [`WindowSpillRecord`] per frame).
pub(crate) struct WindowSpillWriter<W: Write> {
    sink: W,
    scratch: Vec<u8>,
}

impl<W: Write> WindowSpillWriter<W> {
    pub(crate) fn new(sink: W) -> Self {
        Self {
            sink,
            scratch: Vec::new(),
        }
    }

    /// Append one window record (tile-ordered by the caller).
    pub(crate) fn append(&mut self, record: &WindowSpillRecord) -> Result<(), SpillError> {
        self.scratch.clear();
        encode_window(record, &mut self.scratch);
        write_frame(&mut self.sink, &self.scratch)
    }

    /// Flush and return the sink.
    pub(crate) fn finish(mut self) -> Result<W, SpillError> {
        self.sink.flush()?;
        Ok(self.sink)
    }
}

/// Encode a window record: `chrom_id`, `tile`, `gc`, then the per-sample
/// `Option<f32>` depth vector (tag byte + `f32`).
fn encode_window(record: &WindowSpillRecord, out: &mut Vec<u8>) {
    put_u32(out, record.chrom_id);
    put_u32(out, record.tile);
    put_f32(out, record.gc);
    put_len(out, record.depths.len());
    for d in &record.depths {
        match d {
            Some(depth) => {
                put_bool(out, true);
                put_f32(out, *depth);
            }
            None => put_bool(out, false),
        }
    }
}

/// One-pass reader for the window spill, one record at a time. **Test-only**: the
/// production score reads its coverage from the record spill now, so nothing
/// reads the window spill at runtime — this reader survives only to round-trip
/// the still-live writer in its unit tests, until M6 retires the window spill.
#[cfg(test)]
pub(crate) struct WindowSpillReader<R: Read> {
    source: R,
    buf: Vec<u8>,
}

#[cfg(test)]
impl<R: Read> WindowSpillReader<R> {
    pub(crate) fn new(source: R) -> Self {
        Self {
            source,
            buf: Vec::new(),
        }
    }

    /// The next window record, or `None` at a clean end of stream.
    pub(crate) fn next_record(&mut self) -> Option<Result<WindowSpillRecord, SpillError>> {
        match read_frame(&mut self.source, &mut self.buf) {
            Ok(false) => None,
            Ok(true) => Some(decode_window(&self.buf)),
            Err(e) => Some(Err(e)),
        }
    }
}

/// Decode a window record (mirror of [`encode_window`]). Test-only alongside
/// [`WindowSpillReader`].
#[cfg(test)]
fn decode_window(buf: &[u8]) -> Result<WindowSpillRecord, SpillError> {
    let mut d = Decoder::new(buf);
    let chrom_id = d.u32("window.chrom_id")?;
    let tile = d.u32("window.tile")?;
    let gc = d.f32("window.gc")?;
    let n = d.len("window.depths.len")?;
    let mut depths = Vec::with_capacity(n);
    for _ in 0..n {
        depths.push(if d.bool("window.depths.tag")? {
            Some(d.f32("window.depths.value")?)
        } else {
            None
        });
    }
    Ok(WindowSpillRecord {
        chrom_id,
        tile,
        gc,
        depths,
    })
}

/// Read a 4-byte length frame. `Ok(false)` if the stream is cleanly at EOF
/// before any byte of the frame; `Ok(true)` if the frame was filled;
/// [`SpillError::Truncated`] if EOF arrived mid-frame.
fn read_frame_len<R: Read>(source: &mut R, out: &mut [u8; 4]) -> Result<bool, SpillError> {
    let mut filled = 0;
    while filled < out.len() {
        match source.read(&mut out[filled..]) {
            Ok(0) => {
                if filled == 0 {
                    return Ok(false); // clean boundary EOF
                }
                return Err(SpillError::Truncated {
                    field: "length frame",
                    needed: out.len() - filled,
                });
            }
            Ok(n) => filled += n,
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(SpillError::Io(e)),
        }
    }
    Ok(true)
}

/// Write one length-framed payload: a `u32` LE length prefix + the bytes. Shared
/// by the record spill and the window spill.
fn write_frame<W: Write>(sink: &mut W, payload: &[u8]) -> Result<(), SpillError> {
    let len = u32::try_from(payload.len())
        .map_err(|_| SpillError::RecordTooLarge { len: payload.len() })?;
    sink.write_all(&len.to_le_bytes())?;
    sink.write_all(payload)?;
    Ok(())
}

/// Read one length-framed payload into `buf` (cleared first). `Ok(false)` at a
/// clean boundary EOF; `Ok(true)` with `buf` holding the payload; a mid-frame
/// EOF is [`SpillError::Truncated`]. Grows `buf` only as bytes actually arrive
/// (`take` + `read_to_end`), so a corrupt length can't balloon RAM (arch §9).
fn read_frame<R: Read>(source: &mut R, buf: &mut Vec<u8>) -> Result<bool, SpillError> {
    let mut len_bytes = [0u8; 4];
    if !read_frame_len(source, &mut len_bytes)? {
        return Ok(false);
    }
    let len = u32::from_le_bytes(len_bytes) as usize;
    buf.clear();
    let read = source.by_ref().take(len as u64).read_to_end(buf)?;
    if read != len {
        return Err(SpillError::Truncated {
            field: "record payload",
            needed: len - read,
        });
    }
    Ok(true)
}

// ---------------------------------------------------------------------------
// Encoding — little-endian primitives + the record layout.
//
// The encode/decode pair is the whole contract: it must round-trip every field
// bit-for-bit (f64/f32 via `to_le_bytes`, so NaN / ±inf survive — `qual_phred`
// can be +inf). Symmetry is enforced by the round-trip proptest below.
// ---------------------------------------------------------------------------

fn put_u32(out: &mut Vec<u8>, v: u32) {
    out.extend_from_slice(&v.to_le_bytes());
}
fn put_u64(out: &mut Vec<u8>, v: u64) {
    out.extend_from_slice(&v.to_le_bytes());
}
fn put_len(out: &mut Vec<u8>, v: usize) {
    put_u64(out, v as u64);
}
fn put_f32(out: &mut Vec<u8>, v: f32) {
    out.extend_from_slice(&v.to_le_bytes());
}
fn put_f64(out: &mut Vec<u8>, v: f64) {
    out.extend_from_slice(&v.to_le_bytes());
}
fn put_bool(out: &mut Vec<u8>, v: bool) {
    out.push(v as u8);
}
fn put_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    put_len(out, bytes.len());
    out.extend_from_slice(bytes);
}

fn put_support(out: &mut Vec<u8>, s: &AlleleSupportStats) {
    // Destructured (not field-by-field on `s`) so a new `AlleleSupportStats`
    // field is a compile error here — the same refactor-safety net the
    // struct-literal decode gives every other type in this codec. (In-crate, so
    // the `#[non_exhaustive]` on the struct does not block this.)
    let AlleleSupportStats {
        num_obs,
        q_sum,
        fwd,
        placed_left,
        placed_start,
        mapq_sum,
        mapq_sum_sq,
    } = *s;
    put_u32(out, num_obs);
    put_f64(out, q_sum);
    put_u32(out, fwd);
    put_u32(out, placed_left);
    put_u32(out, placed_start);
    put_u32(out, mapq_sum);
    put_u64(out, mapq_sum_sq);
}

fn put_allele(out: &mut Vec<u8>, a: &MergedAllele) {
    put_bytes(out, &a.seq);
    put_bool(out, a.is_compound);
    put_len(out, a.constituents.len());
    for c in &a.constituents {
        put_len(out, c.record_idx);
        put_len(out, c.local_allele_idx);
    }
}

fn encode_record(spill: &ParalogSpillRecord, out: &mut Vec<u8>) {
    let r = &spill.record;
    // locus
    put_u32(out, r.locus.chrom_id);
    put_u32(out, r.locus.start);
    put_u32(out, r.locus.end);
    // alleles
    put_len(out, r.alleles.len());
    for a in &r.alleles {
        put_allele(out, a);
    }
    out.push(r.ploidy);
    put_len(out, r.n_samples);
    put_len(out, r.n_genotypes);
    // allele_frequencies
    put_len(out, r.allele_frequencies.len());
    for &p in &r.allele_frequencies {
        put_f64(out, p);
    }
    // compound_frequencies (Option<f64>)
    put_len(out, r.compound_frequencies.len());
    for opt in &r.compound_frequencies {
        match opt {
            Some(f) => {
                put_bool(out, true);
                put_f64(out, *f);
            }
            None => put_bool(out, false),
        }
    }
    // posteriors
    put_len(out, r.posteriors.len());
    for &p in &r.posteriors {
        put_f64(out, p);
    }
    // best_genotype (usize)
    put_len(out, r.best_genotype.len());
    for &g in &r.best_genotype {
        put_len(out, g);
    }
    // gq_phred
    put_len(out, r.gq_phred.len());
    for &q in &r.gq_phred {
        put_f64(out, q);
    }
    put_f64(out, r.qual_phred);
    // scalars
    put_len(out, r.scalars.len());
    for s in &r.scalars {
        put_support(out, s);
    }
    // other_scalars
    put_len(out, r.other_scalars.len());
    for s in &r.other_scalars {
        put_support(out, s);
    }
    // chain_anchor_flags
    put_len(out, r.chain_anchor_flags.len());
    for &b in &r.chain_anchor_flags {
        put_bool(out, b);
    }
    // diagnostics
    put_u32(out, r.diagnostics.iterations);
    put_f64(out, r.diagnostics.final_max_delta_p);
    put_bool(out, r.diagnostics.converged);
    // the stored paralog LR (NaN = unscored)
    put_f64(out, spill.paralog_lr);
    // per-sample window coverage (gc then coverage; NaN = sample absent)
    put_len(out, spill.window_coverage.gc.len());
    for &g in &spill.window_coverage.gc {
        put_f32(out, g);
    }
    put_len(out, spill.window_coverage.coverage.len());
    for &c in &spill.window_coverage.coverage {
        put_f32(out, c);
    }
}

/// A forward cursor over a record's bytes with range-checked reads.
struct Decoder<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Decoder<'a> {
    fn new(buf: &'a [u8]) -> Self {
        Self { buf, pos: 0 }
    }

    fn take(&mut self, n: usize, field: &'static str) -> Result<&'a [u8], SpillError> {
        let end = self
            .pos
            .checked_add(n)
            .ok_or(SpillError::OutOfRange { field })?;
        if end > self.buf.len() {
            return Err(SpillError::Truncated {
                field,
                needed: end - self.buf.len(),
            });
        }
        let slice = &self.buf[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    fn u32(&mut self, field: &'static str) -> Result<u32, SpillError> {
        let b = self.take(4, field)?;
        Ok(u32::from_le_bytes(b.try_into().expect("4 bytes")))
    }
    fn u64(&mut self, field: &'static str) -> Result<u64, SpillError> {
        let b = self.take(8, field)?;
        Ok(u64::from_le_bytes(b.try_into().expect("8 bytes")))
    }
    fn len(&mut self, field: &'static str) -> Result<usize, SpillError> {
        usize::try_from(self.u64(field)?).map_err(|_| SpillError::OutOfRange { field })
    }
    fn f32(&mut self, field: &'static str) -> Result<f32, SpillError> {
        let b = self.take(4, field)?;
        Ok(f32::from_le_bytes(b.try_into().expect("4 bytes")))
    }
    fn f64(&mut self, field: &'static str) -> Result<f64, SpillError> {
        let b = self.take(8, field)?;
        Ok(f64::from_le_bytes(b.try_into().expect("8 bytes")))
    }
    fn u8(&mut self, field: &'static str) -> Result<u8, SpillError> {
        Ok(self.take(1, field)?[0])
    }
    fn bool(&mut self, field: &'static str) -> Result<bool, SpillError> {
        Ok(self.u8(field)? != 0)
    }
    fn bytes(&mut self, field: &'static str) -> Result<Vec<u8>, SpillError> {
        let n = self.len(field)?;
        Ok(self.take(n, field)?.to_vec())
    }
}

fn get_support(d: &mut Decoder) -> Result<AlleleSupportStats, SpillError> {
    let num_obs = d.u32("scalar.num_obs")?;
    let q_sum = d.f64("scalar.q_sum")?;
    let fwd = d.u32("scalar.fwd")?;
    let placed_left = d.u32("scalar.placed_left")?;
    let placed_start = d.u32("scalar.placed_start")?;
    let mapq_sum = d.u32("scalar.mapq_sum")?;
    let mapq_sum_sq = d.u64("scalar.mapq_sum_sq")?;
    Ok(AlleleSupportStats::new(
        num_obs,
        q_sum,
        fwd,
        placed_left,
        placed_start,
        mapq_sum,
        mapq_sum_sq,
    ))
}

fn get_allele(d: &mut Decoder) -> Result<MergedAllele, SpillError> {
    let seq = d.bytes("allele.seq")?;
    let is_compound = d.bool("allele.is_compound")?;
    let n = d.len("allele.constituents.len")?;
    let mut constituents = Vec::with_capacity(n);
    for _ in 0..n {
        let record_idx = d.len("constituent.record_idx")?;
        let local_allele_idx = d.len("constituent.local_allele_idx")?;
        constituents.push(CompoundConstituent {
            record_idx,
            local_allele_idx,
        });
    }
    Ok(MergedAllele {
        seq,
        is_compound,
        constituents,
    })
}

fn decode_record(buf: &[u8]) -> Result<ParalogSpillRecord, SpillError> {
    let mut d = Decoder::new(buf);
    let locus = RecordLocus {
        chrom_id: d.u32("locus.chrom_id")?,
        start: d.u32("locus.start")?,
        end: d.u32("locus.end")?,
    };
    let n_alleles = d.len("alleles.len")?;
    let mut alleles = Vec::with_capacity(n_alleles);
    for _ in 0..n_alleles {
        alleles.push(get_allele(&mut d)?);
    }
    let ploidy = d.u8("ploidy")?;
    let n_samples = d.len("n_samples")?;
    let n_genotypes = d.len("n_genotypes")?;

    let allele_frequencies = decode_f64_vec(&mut d, "allele_frequencies")?;

    let n_cf = d.len("compound_frequencies.len")?;
    let mut compound_frequencies = Vec::with_capacity(n_cf);
    for _ in 0..n_cf {
        compound_frequencies.push(if d.bool("compound_frequencies.tag")? {
            Some(d.f64("compound_frequencies.value")?)
        } else {
            None
        });
    }

    let posteriors = decode_f64_vec(&mut d, "posteriors")?;

    let n_bg = d.len("best_genotype.len")?;
    let mut best_genotype = Vec::with_capacity(n_bg);
    for _ in 0..n_bg {
        best_genotype.push(d.len("best_genotype.value")?);
    }

    let gq_phred = decode_f64_vec(&mut d, "gq_phred")?;
    let qual_phred = d.f64("qual_phred")?;

    let scalars = decode_support_vec(&mut d, "scalars")?;
    let other_scalars = decode_support_vec(&mut d, "other_scalars")?;

    let n_ca = d.len("chain_anchor_flags.len")?;
    let mut chain_anchor_flags = Vec::with_capacity(n_ca);
    for _ in 0..n_ca {
        chain_anchor_flags.push(d.bool("chain_anchor_flags.value")?);
    }

    let diagnostics = EmDiagnostics {
        iterations: d.u32("diagnostics.iterations")?,
        final_max_delta_p: d.f64("diagnostics.final_max_delta_p")?,
        converged: d.bool("diagnostics.converged")?,
    };
    let paralog_lr = d.f64("paralog_lr")?;
    let window_coverage = LocusWindowCoverage {
        gc: decode_f32_vec(&mut d, "window_coverage.gc")?,
        coverage: decode_f32_vec(&mut d, "window_coverage.coverage")?,
    };

    Ok(ParalogSpillRecord {
        paralog_lr,
        window_coverage,
        record: PosteriorRecord {
            locus,
            alleles,
            ploidy,
            n_samples,
            n_genotypes,
            allele_frequencies,
            compound_frequencies,
            posteriors,
            best_genotype,
            gq_phred,
            qual_phred,
            scalars,
            other_scalars,
            chain_anchor_flags,
            diagnostics,
        },
    })
}

fn decode_f64_vec(d: &mut Decoder, field: &'static str) -> Result<Vec<f64>, SpillError> {
    let n = d.len(field)?;
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        v.push(d.f64(field)?);
    }
    Ok(v)
}

fn decode_f32_vec(d: &mut Decoder, field: &'static str) -> Result<Vec<f32>, SpillError> {
    let n = d.len(field)?;
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        v.push(d.f32(field)?);
    }
    Ok(v)
}

fn decode_support_vec(
    d: &mut Decoder,
    field: &'static str,
) -> Result<Vec<AlleleSupportStats>, SpillError> {
    let n = d.len(field)?;
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        v.push(get_support(d)?);
    }
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn support(num_obs: u32, q_sum: f64) -> AlleleSupportStats {
        AlleleSupportStats::new(num_obs, q_sum, num_obs / 2, 1, 2, 30, 900)
    }

    fn allele(seq: &[u8], is_compound: bool) -> MergedAllele {
        MergedAllele {
            seq: seq.to_vec(),
            is_compound,
            constituents: if is_compound {
                vec![
                    CompoundConstituent {
                        record_idx: 3,
                        local_allele_idx: 1,
                    },
                    CompoundConstituent {
                        record_idx: 7,
                        local_allele_idx: 0,
                    },
                ]
            } else {
                Vec::new()
            },
        }
    }

    /// A representative biallelic-SNP spill record for a 2-sample cohort, with
    /// non-trivial values in every field (incl. `qual_phred = +inf`, a compound
    /// allele) so the round-trip exercises all branches.
    fn spill_record(pos: u32) -> ParalogSpillRecord {
        ParalogSpillRecord {
            // A finite, pos-varying LR so the round-trip's `assert_eq` exercises
            // the field (NaN would break equality — that path is covered by the
            // dedicated nan_lr round-trip test below).
            paralog_lr: (pos as f64) * 0.25 - 5.0,
            // Distinct finite per-sample window coverage so the round-trip's
            // `assert_eq` exercises the field (NaN handling covered separately).
            window_coverage: LocusWindowCoverage {
                gc: vec![0.31 + pos as f32 * 0.001, 0.44],
                coverage: vec![100.0 + pos as f32, 205.5],
            },
            record: PosteriorRecord {
                locus: RecordLocus {
                    chrom_id: 2,
                    start: pos,
                    end: pos,
                },
                alleles: vec![allele(b"A", false), allele(b"T", true)],
                ploidy: 2,
                n_samples: 2,
                n_genotypes: 3,
                allele_frequencies: vec![0.7, 0.3],
                compound_frequencies: vec![None, Some(0.15)],
                posteriors: vec![0.8, 0.15, 0.05, 0.1, 0.7, 0.2],
                best_genotype: vec![0, 1],
                gq_phred: vec![42.0, 99.0],
                qual_phred: f64::INFINITY,
                scalars: vec![
                    support(20, -3.5),
                    support(0, 0.0),
                    support(11, -1.2),
                    support(9, -0.9),
                ],
                other_scalars: vec![support(1, -0.1)],
                chain_anchor_flags: vec![false, true, false, false],
                diagnostics: EmDiagnostics {
                    iterations: 12,
                    final_max_delta_p: 1e-7,
                    converged: true,
                },
            },
        }
    }

    /// N records round-trip through an in-memory stream bit-for-bit (including
    /// `+inf` QUAL, a `None` window, and a compound allele's constituents).
    #[test]
    fn round_trips_records_through_a_stream() {
        let records: Vec<ParalogSpillRecord> = (0..4).map(|i| spill_record(100 + i)).collect();
        let mut writer = ParalogSpillWriter::new(Cursor::new(Vec::new()));
        for r in &records {
            writer.append(r).expect("append");
        }
        assert_eq!(writer.records_written(), 4);
        let bytes = writer.finish().expect("finish").into_inner();

        let mut reader = ParalogSpillReader::new(Cursor::new(bytes));
        let read_back: Vec<ParalogSpillRecord> = std::iter::from_fn(|| reader.next_record())
            .map(|r| r.expect("decode"))
            .collect();
        assert_eq!(read_back, records);
    }

    /// The `paralog_lr = NaN` unscored sentinel round-trips (an `assert_eq`
    /// can't check it — `NaN != NaN` — so assert `is_nan()` explicitly). This is
    /// the value the sink stores today for every locus.
    #[test]
    fn nan_lr_round_trips_as_unscored() {
        let mut rec = spill_record(42);
        rec.paralog_lr = f64::NAN;
        let mut buf = Vec::new();
        encode_record(&rec, &mut buf);
        let back = decode_record(&buf).expect("decode");
        assert!(back.paralog_lr.is_nan(), "NaN LR sentinel must survive");
        assert_eq!(back.record, rec.record, "the record round-trips unchanged");
    }

    /// A `f64::INFINITY` QUAL survives the round-trip exactly (a Phred passed
    /// through when a sample is certainly variant — arch note on `qual_phred`).
    #[test]
    fn infinite_qual_round_trips() {
        let mut buf = Vec::new();
        let rec = spill_record(1);
        encode_record(&rec, &mut buf);
        let back = decode_record(&buf).expect("decode");
        assert!(back.record.qual_phred.is_infinite());
        assert_eq!(back, rec);
    }

    /// A `Read` wrapper counting bytes pulled from the underlying stream, to
    /// prove the reader streams one record at a time rather than slurping.
    struct CountingReader<R> {
        inner: R,
        read: &'static AtomicUsize,
    }
    impl<R: Read> Read for CountingReader<R> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = self.inner.read(buf)?;
            self.read.fetch_add(n, Ordering::SeqCst);
            Ok(n)
        }
    }

    /// The reader is lazy: after one `next_record`, only the first record's
    /// frame has been consumed from the source (not the whole stream), and the
    /// decode buffer holds a single record's bytes — memory flat in record
    /// count.
    #[test]
    fn reader_consumes_one_record_at_a_time() {
        static READ: AtomicUsize = AtomicUsize::new(0);
        READ.store(0, Ordering::SeqCst);

        let records: Vec<ParalogSpillRecord> = (0..5).map(|i| spill_record(200 + i)).collect();
        let mut writer = ParalogSpillWriter::new(Cursor::new(Vec::new()));
        for r in &records {
            writer.append(r).unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();
        let total_len = bytes.len();
        // One record's on-disk size (frame + payload), all records equal shape.
        let one_record_len = total_len / records.len();

        let mut reader = ParalogSpillReader::new(BufReader::with_capacity(
            8,
            CountingReader {
                inner: Cursor::new(bytes),
                read: &READ,
            },
        ));
        let first = reader.next_record().expect("some").expect("ok");
        assert_eq!(first, records[0]);
        // Far less than the whole stream has been pulled — the reader did not
        // materialise all records to yield the first.
        let consumed = READ.load(Ordering::SeqCst);
        assert!(
            consumed < total_len,
            "reader pulled {consumed} of {total_len} bytes for the first of {} records",
            records.len()
        );
        // The reused decode buffer holds exactly one record's payload.
        assert!(reader.buf.len() < one_record_len + 4);
    }

    /// Window-spill records round-trip through a stream (including a `None`
    /// per-sample depth and a `NaN`-free GC).
    #[test]
    fn window_records_round_trip() {
        let records = vec![
            WindowSpillRecord {
                chrom_id: 0,
                tile: 2,
                gc: 0.41,
                depths: vec![Some(18.3), None, Some(0.0)],
            },
            WindowSpillRecord {
                chrom_id: 1,
                tile: 0,
                gc: 1.0,
                depths: vec![None, None],
            },
        ];
        let mut writer = WindowSpillWriter::new(Cursor::new(Vec::new()));
        for r in &records {
            writer.append(r).unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();
        let mut reader = WindowSpillReader::new(Cursor::new(bytes));
        let back: Vec<WindowSpillRecord> = std::iter::from_fn(|| reader.next_record())
            .map(|r| r.unwrap())
            .collect();
        assert_eq!(back, records);
    }

    /// A stream truncated inside a record's payload surfaces as an error, not a
    /// silent short record or a clean end.
    #[test]
    fn truncated_payload_is_an_error() {
        let mut writer = ParalogSpillWriter::new(Cursor::new(Vec::new()));
        writer.append(&spill_record(1)).unwrap();
        let mut bytes = writer.finish().unwrap().into_inner();
        bytes.truncate(bytes.len() - 3); // chop the tail of the payload
        let mut reader = ParalogSpillReader::new(Cursor::new(bytes));
        let result = reader.next_record().expect("frame present");
        assert!(result.is_err(), "truncated payload must error");
    }

    /// A stream that ends inside the 4-byte length prefix is a `Truncated`
    /// error on the frame, not a clean end and not a panic.
    #[test]
    fn truncated_length_prefix_is_an_error() {
        let mut writer = ParalogSpillWriter::new(Cursor::new(Vec::new()));
        writer.append(&spill_record(1)).unwrap();
        let mut bytes = writer.finish().unwrap().into_inner();
        bytes.truncate(2); // only half of the length prefix survives
        let mut reader = ParalogSpillReader::new(Cursor::new(bytes));
        match reader.next_record() {
            Some(Err(SpillError::Truncated { field, .. })) => assert_eq!(field, "length frame"),
            other => panic!("expected Truncated(length frame), got {other:?}"),
        }
    }

    /// An empty stream yields a clean `None` at the first record boundary.
    #[test]
    fn empty_stream_yields_none() {
        let mut reader = ParalogSpillReader::new(Cursor::new(Vec::<u8>::new()));
        assert!(reader.next_record().is_none());
    }

    /// A corrupt length prefix far larger than the remaining bytes surfaces as
    /// `Truncated` without pre-allocating the claimed size (memory-flatness).
    #[test]
    fn corrupt_huge_length_does_not_balloon() {
        // A 12-byte stream: a length prefix claiming ~4 GiB, then 8 payload bytes.
        let mut bytes = u32::MAX.to_le_bytes().to_vec();
        bytes.extend_from_slice(&[0u8; 8]);
        let mut reader = ParalogSpillReader::new(Cursor::new(bytes));
        match reader.next_record() {
            Some(Err(SpillError::Truncated { field, .. })) => assert_eq!(field, "record payload"),
            other => panic!("expected Truncated(record payload), got {other:?}"),
        }
        // The buffer grew only to the bytes actually present, not to ~4 GiB.
        assert!(reader.buf.len() <= 8);
    }

    proptest::proptest! {
        /// Round-trip over randomly-shaped records: variable allele / sample /
        /// genotype counts, arbitrary finite floats, optional windows and
        /// compound frequencies. Any encode/decode asymmetry (a field written
        /// in the wrong order or width) fails here where a single hand-built
        /// example would miss it.
        #[test]
        fn round_trips_arbitrary_records(
            chrom_id in 0u32..30,
            start in 1u32..1_000_000,
            n_alleles in 1usize..4,
            n_samples in 0usize..6,
            ploidy in 1u8..5,
            seed in proptest::prelude::any::<u64>(),
        ) {
            // Derive deterministic-but-varied field values from the seed so the
            // floats and flags differ across cases without another strategy.
            let f = |k: u64| ((seed ^ k).wrapping_mul(2654435761) % 10_000) as f64 / 9_999.0;
            let n_genotypes = n_alleles * usize::from(ploidy) + 1; // any valid shape; not semantically checked
            let alleles: Vec<MergedAllele> = (0..n_alleles)
                .map(|a| allele(if a == 0 { b"A" } else { b"C" }, a % 2 == 1))
                .collect();
            let rec = ParalogSpillRecord {
                paralog_lr: f(700) * 20.0 - 10.0,
                window_coverage: LocusWindowCoverage {
                    gc: (0..n_samples).map(|s| f(s as u64 + 800) as f32).collect(),
                    coverage: (0..n_samples).map(|s| f(s as u64 + 900) as f32 * 300.0).collect(),
                },
                record: PosteriorRecord {
                    locus: RecordLocus { chrom_id, start, end: start },
                    alleles,
                    ploidy,
                    n_samples,
                    n_genotypes,
                    allele_frequencies: (0..n_alleles).map(|a| f(a as u64)).collect(),
                    compound_frequencies: (0..n_alleles)
                        .map(|a| if a % 2 == 0 { None } else { Some(f(a as u64 + 100)) })
                        .collect(),
                    posteriors: (0..n_samples * n_genotypes).map(|i| f(i as u64 + 200)).collect(),
                    best_genotype: (0..n_samples).map(|s| (s + seed as usize) % n_genotypes).collect(),
                    gq_phred: (0..n_samples).map(|s| f(s as u64 + 300) * 99.0).collect(),
                    qual_phred: f(400) * 1000.0,
                    scalars: (0..n_samples * n_alleles)
                        .map(|i| support(i as u32, -f(i as u64 + 500)))
                        .collect(),
                    other_scalars: Vec::new(),
                    chain_anchor_flags: (0..n_samples * n_alleles).map(|i| i % 3 == 0).collect(),
                    diagnostics: EmDiagnostics {
                        iterations: (seed % 100) as u32,
                        final_max_delta_p: f(600) * 1e-3,
                        converged: seed % 2 == 0,
                    },
                },
            };
            let mut buf = Vec::new();
            encode_record(&rec, &mut buf);
            let back = decode_record(&buf).expect("decode");
            proptest::prop_assert_eq!(back, rec);
        }
    }

    /// The RAII owner creates a file in the given dir, round-trips through it,
    /// and deletes it on drop (both success and, by construction, failure).
    #[test]
    fn spill_file_round_trips_and_is_deleted_on_drop() {
        let dir = tempfile::tempdir().unwrap();
        let path = {
            let spill = ParalogSpill::create_in(dir.path()).expect("create spill");
            let path = spill.path().to_path_buf();
            assert!(path.exists(), "spill file exists after create");

            let records: Vec<ParalogSpillRecord> = (0..3).map(|i| spill_record(300 + i)).collect();
            let mut writer = spill.writer().expect("writer");
            for r in &records {
                writer.append(r).unwrap();
            }
            writer.finish().unwrap();

            let mut reader = spill.reader().expect("reader");
            let read_back: Vec<ParalogSpillRecord> = std::iter::from_fn(|| reader.next_record())
                .map(|r| r.unwrap())
                .collect();
            assert_eq!(read_back, records);
            // A second reader re-reads the same file (the spill is read twice).
            let second = spill.reader().expect("reader 2").count();
            assert_eq!(second, 3);
            path
        };
        assert!(
            !path.exists(),
            "spill file is deleted when ParalogSpill drops"
        );
    }
}
