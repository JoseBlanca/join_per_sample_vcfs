//! Streaming writer that turns a [`PileupRecord`] stream into a v1.0
//! `.psp` file.
//!
//! Three public entry points:
//!
//! - [`PspWriter::new`] — emits the framed TOML header to the sink
//!   and prepares the in-memory state for accepting records.
//! - [`PspWriter::write_record`] — validates one record and buffers
//!   it into the open block; auto-flushes when the projected
//!   uncompressed payload reaches [`DEFAULT_TARGET_BLOCK_BYTES`] or
//!   the next record crosses a chromosome boundary.
//! - [`PspWriter::finish`] — flushes any open block, writes the
//!   block index and the trailer (with its XXH3-64 index checksum),
//!   and returns the sink.
//!
//! Every per-record / per-block invariant from spec
//! §"Header-binary consistency" and §"Block sizing / Block
//! invariants" is enforced at write_record / flush time, so a sink
//! that received a sequence of `Ok` returns is *necessarily* a
//! reader-valid `.psp`. The reader's checks (in step 6) duplicate
//! these — that is by design, because the file may be read on a
//! machine that did not write it.

use std::collections::BTreeSet;
use std::io::Write;

use super::block::{
    BlockHeader, ColumnManifestEntry, encode_block_header, encode_bytes_concat, encode_list_column,
    encode_scalar_column, encode_varint_column, zstd_compress,
};
use super::errors::PspWriteError;
use super::header::{WriterHeader, build_header_bytes};
use super::index::{BlockIndexEntry, checksum_index, encode_index};
use super::registry::{
    Cardinality, ColumnDef, ElementType, MAX_ALLELE_SEQ_LEN, Shape, V1_0_COLUMNS,
};
use super::trailer::{Trailer, encode_trailer};
use crate::per_sample_caller::pileup::{PileupRecord, SlotId};

/// Target uncompressed bytes per block. Spec §"Block sizing". Not
/// CLI-exposed; the writer flushes around this value.
pub const DEFAULT_TARGET_BLOCK_BYTES: usize = 16 * 1024 * 1024;

/// Streaming `.psp` writer over any `Write` sink.
pub struct PspWriter<W: Write> {
    sink: W,
    header: WriterHeader,
    /// Bytes written to `sink` so far. Updated after every flush so
    /// the block index records absolute offsets.
    sink_offset: u64,
    /// One entry per emitted block.
    index_entries: Vec<BlockIndexEntry>,
    /// Open block — `None` when the writer is between blocks (just
    /// after `new` or just after a flush).
    block: Option<BlockAccumulator>,
    /// Running active phase-chain slot set. Updated as records'
    /// `new_chains` / `expired_chains` markers are applied;
    /// snapshotted at block start.
    active_slots: BTreeSet<SlotId>,
    /// Last admitted record's `(chrom_id, pos)`, for monotonicity
    /// enforcement across `write_record` calls.
    last_locus: Option<(u32, u32)>,
    /// Running count of records the caller has handed us. Used in
    /// error messages so the offending input is identifiable.
    records_seen: u64,
}

impl<W: Write> PspWriter<W> {
    /// Validate the header, frame it (magic + length prefix + TOML +
    /// sentinel), and write it to `sink`. After this returns, the
    /// writer is ready to accept records.
    ///
    /// **Buffering.** Each block flush issues one `write_all` for the
    /// block header plus one per column payload, and `finish` adds
    /// one for the index and one for the 32-byte trailer. Against a
    /// real-file sink each `write_all` becomes one `write(2)`
    /// syscall. Wrap a `File` in
    /// `BufWriter::with_capacity(64 * 1024, file)` (or larger) before
    /// passing it here. In-memory sinks (`io::sink()`,
    /// `Cursor<Vec<u8>>`) need no wrapper.
    pub fn new(mut sink: W, header: WriterHeader) -> Result<Self, PspWriteError> {
        let header_bytes = build_header_bytes(&header)?;
        sink.write_all(&header_bytes)
            .map_err(|e| PspWriteError::Io {
                context: "file header",
                source: e,
            })?;
        let sink_offset = header_bytes.len() as u64;
        Ok(Self {
            sink,
            header,
            sink_offset,
            index_entries: Vec::new(),
            block: None,
            active_slots: BTreeSet::new(),
            last_locus: None,
            records_seen: 0,
        })
    }

    /// Append one record. Returns the number of bytes pushed to the
    /// sink as a side effect of any auto-flush this call triggered
    /// (zero if no flush).
    pub fn write_record(&mut self, record: &PileupRecord) -> Result<u64, PspWriteError> {
        let record_index = self.records_seen;
        self.records_seen += 1;
        self.validate_record(record_index, record)?;

        // Should we flush before appending? Flush triggers:
        // - chrom_id change (blocks never cross chromosomes)
        // - projected uncompressed size at or past the target
        let pre_flush_bytes = self.sink_offset;
        let should_flush = match &self.block {
            None => false,
            Some(b) => {
                b.chrom_id != record.chrom_id || b.projected_bytes >= DEFAULT_TARGET_BLOCK_BYTES
            }
        };
        if should_flush {
            self.flush_block()?;
        }
        let flushed_bytes = self.sink_offset - pre_flush_bytes;

        // Start a new block if none is open.
        if self.block.is_none() {
            // The first block on a chromosome always has an empty
            // active-slot snapshot — Stage 1 never produces a chain
            // that spans a contig boundary (reads don't), so by the
            // time we cross a boundary every slot has expired. The
            // walker is expected to honour this and we double-check
            // here: if active_slots is non-empty AND we're starting
            // the first block on a new chromosome, that's a writer
            // upstream bug.
            let snapshot = if self.is_first_block_on_chrom(record.chrom_id) {
                if !self.active_slots.is_empty() {
                    return Err(PspWriteError::PhaseChainMarkerInconsistency {
                        record_index,
                        reason: format!(
                            "active_slots non-empty ({}) at first record of chromosome {}",
                            self.active_slots.len(),
                            record.chrom_id
                        ),
                    });
                }
                Vec::new()
            } else {
                self.active_slots.iter().copied().collect()
            };
            self.block = Some(BlockAccumulator::new(record.chrom_id, record.pos, snapshot));
        }

        self.apply_record_to_block(record_index, record)?;
        self.last_locus = Some((record.chrom_id, record.pos));

        Ok(flushed_bytes)
    }

    /// Flush any open block, write the block index, and write the
    /// trailer. Consumes the writer; the returned sink is positioned
    /// at the end of a complete `.psp` file.
    ///
    /// **End-of-stage discipline for `BufWriter`-wrapped file sinks.**
    /// `BufWriter::drop` may swallow flush errors, which for a
    /// billions-of-records file can silently truncate the trailer.
    /// The production caller should do, in order:
    /// ```ignore
    /// let buf = writer.finish()?;       // PSP-level flush
    /// let file = buf.into_inner()?;     // surface BufWriter errors
    /// file.sync_all()?;                 // durability for downstream stages
    /// ```
    /// `sync_all` is end-of-stage only — never per-block.
    pub fn finish(mut self) -> Result<W, PspWriteError> {
        if self.block.is_some() {
            self.flush_block()?;
        }
        let index_offset = self.sink_offset;
        let index_bytes = encode_index(&self.index_entries);
        let index_byte_length = index_bytes.len() as u64;
        let index_checksum = checksum_index(&index_bytes);
        self.sink
            .write_all(&index_bytes)
            .map_err(|e| PspWriteError::Io {
                context: "block index",
                source: e,
            })?;
        self.sink_offset += index_byte_length;

        let trailer = Trailer {
            index_offset,
            index_byte_length,
            n_blocks: self.index_entries.len() as u64,
            index_checksum,
        };
        let trailer_bytes = encode_trailer(&trailer);
        self.sink
            .write_all(&trailer_bytes)
            .map_err(|e| PspWriteError::Io {
                context: "file trailer",
                source: e,
            })?;
        self.sink_offset += trailer_bytes.len() as u64;
        self.sink.flush().map_err(|e| PspWriteError::Io {
            context: "final flush",
            source: e,
        })?;
        Ok(self.sink)
    }

    /// `true` if no block has been emitted yet OR the most recent
    /// emitted block was on a different chromosome.
    fn is_first_block_on_chrom(&self, chrom_id: u32) -> bool {
        match self.index_entries.last() {
            None => true,
            Some(last) => last.chrom_id != chrom_id,
        }
    }

    fn validate_record(
        &self,
        record_index: u64,
        record: &PileupRecord,
    ) -> Result<(), PspWriteError> {
        // Chromosome id known.
        let n_chroms = self.header.chromosomes.len();
        if record.chrom_id as usize >= n_chroms {
            return Err(PspWriteError::UnknownChromId {
                record_index,
                chrom_id: record.chrom_id,
                n_chroms: n_chroms as u32,
            });
        }
        let chrom = &self.header.chromosomes[record.chrom_id as usize];

        // pos in [1, chrom.length].
        if record.pos == 0 || record.pos > chrom.length {
            return Err(PspWriteError::PosOutOfRange {
                record_index,
                chrom_id: record.chrom_id,
                pos: record.pos,
                chrom_length: chrom.length,
            });
        }

        // Monotonic locus.
        if let Some((prev_chrom, prev_pos)) = self.last_locus {
            let regression = record.chrom_id < prev_chrom
                || (record.chrom_id == prev_chrom && record.pos <= prev_pos);
            if regression {
                return Err(PspWriteError::OutOfOrderRecord {
                    record_index,
                    prev_chrom,
                    prev_pos,
                    this_chrom: record.chrom_id,
                    this_pos: record.pos,
                });
            }
        }

        // At least one allele.
        if record.alleles.is_empty() {
            return Err(PspWriteError::InvalidRecord {
                record_index,
                reason: "record has zero alleles (REF always present)".to_string(),
            });
        }

        // Per-allele rules.
        for (i, allele) in record.alleles.iter().enumerate() {
            let len = allele.seq.len();
            if len == 0 || (len as u64) > MAX_ALLELE_SEQ_LEN {
                return Err(PspWriteError::InvalidRecord {
                    record_index,
                    reason: format!(
                        "allele {i} sequence length {len} outside [1, {MAX_ALLELE_SEQ_LEN}]"
                    ),
                });
            }
            for (j, &b) in allele.seq.iter().enumerate() {
                if !matches!(b, b'A' | b'C' | b'G' | b'T' | b'N') {
                    return Err(PspWriteError::InvalidRecord {
                        record_index,
                        reason: format!("allele {i} byte {j} = {b:#04x} (only A/C/G/T/N allowed)"),
                    });
                }
            }
            if !allele.support.q_sum.is_finite() {
                return Err(PspWriteError::InvalidRecord {
                    record_index,
                    reason: format!("allele {i} q_sum is non-finite ({})", allele.support.q_sum),
                });
            }
            // chain_slots ascending + distinct (defensive).
            for w in allele.chain_slots.windows(2) {
                if w[0] >= w[1] {
                    return Err(PspWriteError::InvalidRecord {
                        record_index,
                        reason: format!("allele {i} chain_slots not strictly ascending"),
                    });
                }
            }
        }

        // new / expired chains ascending + distinct.
        for (name, slots) in [
            ("new_chains", &record.new_chains),
            ("expired_chains", &record.expired_chains),
        ] {
            for w in slots.windows(2) {
                if w[0] >= w[1] {
                    return Err(PspWriteError::InvalidRecord {
                        record_index,
                        reason: format!("{name} not strictly ascending"),
                    });
                }
            }
        }

        Ok(())
    }

    /// Apply the record's content to the open block's per-column
    /// buffers, update the running active-slot set, and validate
    /// phase-chain marker consistency against that set.
    fn apply_record_to_block(
        &mut self,
        record_index: u64,
        record: &PileupRecord,
    ) -> Result<(), PspWriteError> {
        // Phase-chain marker consistency: expired must be currently
        // active; new must not be currently active.
        for &slot in &record.expired_chains {
            if !self.active_slots.contains(&slot) {
                return Err(PspWriteError::PhaseChainMarkerInconsistency {
                    record_index,
                    reason: format!("expired chain slot {slot} not currently active"),
                });
            }
        }
        for &slot in &record.new_chains {
            if self.active_slots.contains(&slot) {
                return Err(PspWriteError::PhaseChainMarkerInconsistency {
                    record_index,
                    reason: format!("new chain slot {slot} already active"),
                });
            }
        }
        // Apply markers to running set.
        for &slot in &record.expired_chains {
            self.active_slots.remove(&slot);
        }
        for &slot in &record.new_chains {
            self.active_slots.insert(slot);
        }
        // Now per-allele chain_slots must all be in the running set.
        for (i, allele) in record.alleles.iter().enumerate() {
            for &slot in &allele.chain_slots {
                if !self.active_slots.contains(&slot) {
                    return Err(PspWriteError::PhaseChainMarkerInconsistency {
                        record_index,
                        reason: format!(
                            "allele {i} references chain slot {slot} not in active set"
                        ),
                    });
                }
            }
        }

        // Push to block buffers.
        let block = self.block.as_mut().expect("block open by construction");
        block.append_record(record);

        Ok(())
    }

    /// Encode the current block to wire bytes, compress its columns,
    /// emit header + payloads to the sink, and append an index
    /// entry.
    fn flush_block(&mut self) -> Result<(), PspWriteError> {
        let block = self
            .block
            .take()
            .expect("flush_block called with no open block");
        let block_offset = self.sink_offset;
        let n_records = block.delta_pos.len() as u32;
        let n_total_alleles = block.allele_seq_len.len() as u32;

        // Encode every column, compress each, build the manifest.
        let mut payloads: Vec<(u16, Vec<u8>, u32)> = Vec::with_capacity(V1_0_COLUMNS.len());
        for column_def in V1_0_COLUMNS {
            let payload = encode_column(column_def, &block)?;
            let uncompressed_len = payload.len() as u32;
            let compressed = zstd_compress(&payload).map_err(|e| PspWriteError::Io {
                context: "zstd compression of column payload",
                source: e,
            })?;
            payloads.push((column_def.tag, compressed, uncompressed_len));
        }
        let manifest: Vec<ColumnManifestEntry> = payloads
            .iter()
            .map(|(tag, compressed, uncompressed_len)| ColumnManifestEntry {
                tag: *tag,
                compressed_len: compressed.len() as u32,
                uncompressed_len: *uncompressed_len,
            })
            .collect();

        let block_index = self.index_entries.len() as u64;
        let header = BlockHeader {
            chrom_id: block.chrom_id,
            first_pos: block.first_pos,
            n_records,
            n_total_alleles,
            active_chain_slots: block.snapshot_active_slots.clone(),
            manifest,
        };
        let mut header_bytes = Vec::new();
        encode_block_header(&header, &mut header_bytes).map_err(|source| {
            PspWriteError::BlockEmission {
                block_index,
                source,
            }
        })?;

        self.sink
            .write_all(&header_bytes)
            .map_err(|e| PspWriteError::Io {
                context: "block header",
                source: e,
            })?;
        let mut written = header_bytes.len() as u64;
        for (_, compressed, _) in &payloads {
            self.sink
                .write_all(compressed)
                .map_err(|e| PspWriteError::Io {
                    context: "block column payload",
                    source: e,
                })?;
            written += compressed.len() as u64;
        }
        self.sink_offset += written;

        self.index_entries.push(BlockIndexEntry {
            chrom_id: block.chrom_id,
            first_pos: block.first_pos,
            last_pos: block.last_pos,
            n_records,
            block_offset,
        });
        Ok(())
    }
}

// ---------------------------------------------------------------------
// BlockAccumulator
// ---------------------------------------------------------------------

struct BlockAccumulator {
    chrom_id: u32,
    first_pos: u32,
    last_pos: u32,
    /// Snapshot of the active-slot set at the moment this block
    /// opened. Sits verbatim in the block header.
    snapshot_active_slots: Vec<SlotId>,
    // Per-record columns.
    delta_pos: Vec<u64>,
    n_alleles: Vec<u64>,
    new_chain_slots: Vec<Vec<SlotId>>,
    expired_chain_slots: Vec<Vec<SlotId>>,
    // Per-allele columns.
    allele_seq_len: Vec<u64>,
    allele_seq_bytes: Vec<u8>,
    allele_obs_count: Vec<u32>,
    allele_q_sum_log: Vec<f64>,
    allele_fwd_count: Vec<u32>,
    allele_placed_left_count: Vec<u32>,
    allele_placed_start_count: Vec<u32>,
    allele_chain_slots: Vec<Vec<SlotId>>,
    /// Rough projection of the uncompressed byte total. Used to
    /// decide when to auto-flush.
    projected_bytes: usize,
}

impl BlockAccumulator {
    fn new(chrom_id: u32, first_pos: u32, snapshot_active_slots: Vec<SlotId>) -> Self {
        Self {
            chrom_id,
            first_pos,
            last_pos: first_pos,
            snapshot_active_slots,
            delta_pos: Vec::new(),
            n_alleles: Vec::new(),
            new_chain_slots: Vec::new(),
            expired_chain_slots: Vec::new(),
            allele_seq_len: Vec::new(),
            allele_seq_bytes: Vec::new(),
            allele_obs_count: Vec::new(),
            allele_q_sum_log: Vec::new(),
            allele_fwd_count: Vec::new(),
            allele_placed_left_count: Vec::new(),
            allele_placed_start_count: Vec::new(),
            allele_chain_slots: Vec::new(),
            projected_bytes: 0,
        }
    }

    fn append_record(&mut self, record: &PileupRecord) {
        let is_first = self.delta_pos.is_empty();
        let delta = if is_first {
            0u64
        } else {
            (record.pos - self.last_pos) as u64
        };
        self.delta_pos.push(delta);
        self.n_alleles.push(record.alleles.len() as u64);
        self.new_chain_slots.push(record.new_chains.clone());
        self.expired_chain_slots.push(record.expired_chains.clone());

        for allele in &record.alleles {
            self.allele_seq_len.push(allele.seq.len() as u64);
            self.allele_seq_bytes.extend_from_slice(&allele.seq);
            self.allele_obs_count.push(allele.support.num_obs);
            self.allele_q_sum_log.push(allele.support.q_sum);
            self.allele_fwd_count.push(allele.support.fwd);
            self.allele_placed_left_count
                .push(allele.support.placed_left);
            self.allele_placed_start_count
                .push(allele.support.placed_start);
            self.allele_chain_slots.push(allele.chain_slots.clone());
        }

        self.last_pos = record.pos;

        // Rough size projection — per record, plus per-allele +
        // per-byte. Doesn't need to be precise; the target is a soft
        // cap.
        let per_record = 1 // delta-pos varint typical
            + 1 // n-alleles varint typical
            + (1 + 2 * record.new_chains.len()) // new_chain marker
            + (1 + 2 * record.expired_chains.len()); // expired_chain marker
        let per_allele: usize = record
            .alleles
            .iter()
            .map(|a| {
                1                  // allele-seq-len varint typical
                + a.seq.len()      // allele-seq bytes
                + 4 + 8 + 4 + 4 + 4 // the five scalars
                + (1 + 2 * a.chain_slots.len())
            })
            .sum();
        self.projected_bytes += per_record + per_allele;
    }
}

// ---------------------------------------------------------------------
// Column encoders — dispatched by (cardinality, shape, element-type)
// ---------------------------------------------------------------------

/// Encode one column from the block accumulator into its
/// uncompressed wire form. Dispatches by tag — each tag corresponds
/// to a specific accumulator field with a specific element type.
fn encode_column(def: &ColumnDef, block: &BlockAccumulator) -> Result<Vec<u8>, PspWriteError> {
    let mut out = Vec::new();
    match def.tag {
        0x01 => encode_varint_column(&block.delta_pos, &mut out),
        0x02 => encode_varint_column(&block.n_alleles, &mut out),
        0x03 => encode_varint_column(&block.allele_seq_len, &mut out),
        0x04 => encode_bytes_concat(
            &[&block.allele_seq_bytes[..]],
            // The bytes-column payload is the concatenation of every
            // allele's sequence bytes; the per-allele lengths live
            // in the 0x03 length column. We've kept the bytes
            // concatenated as we appended records, so the payload
            // is just `allele_seq_bytes` as-is.
            &mut out,
        ),
        0x10 => encode_scalar_column(&block.allele_obs_count, &mut out),
        0x11 => {
            // q_sum was finite-validated at write_record time.
            encode_scalar_column(&block.allele_q_sum_log, &mut out)
        }
        0x12 => encode_scalar_column(&block.allele_fwd_count, &mut out),
        0x13 => encode_scalar_column(&block.allele_placed_left_count, &mut out),
        0x14 => encode_scalar_column(&block.allele_placed_start_count, &mut out),
        0x20 => {
            let slices: Vec<&[u16]> = block.new_chain_slots.iter().map(|v| v.as_slice()).collect();
            encode_list_column(&slices, &mut out);
        }
        0x21 => {
            let slices: Vec<&[u16]> = block
                .expired_chain_slots
                .iter()
                .map(|v| v.as_slice())
                .collect();
            encode_list_column(&slices, &mut out);
        }
        0x22 => {
            let slices: Vec<&[u16]> = block
                .allele_chain_slots
                .iter()
                .map(|v| v.as_slice())
                .collect();
            encode_list_column(&slices, &mut out);
        }
        unknown => {
            // V1_0_COLUMNS is the source of truth; if encode_column
            // is called with a tag not in this match, the registry
            // has been extended without updating the dispatch — a
            // programmer error.
            return Err(PspWriteError::InvalidRecord {
                record_index: 0,
                reason: format!("internal: encode_column has no dispatch for tag {unknown:#x}"),
            });
        }
    }
    // Verify the encoded size matches the schema-predicted size for
    // fixed-width scalars and the bytes column (with the length
    // column also known). This is a writer-side self-check: a
    // mis-match means a bug in the writer.
    debug_assert!(verify_encoded_size(def, block, &out));
    Ok(out)
}

fn verify_encoded_size(def: &ColumnDef, block: &BlockAccumulator, encoded: &[u8]) -> bool {
    let n_records = block.delta_pos.len();
    let n_total_alleles = block.allele_seq_len.len();
    match (def.cardinality, def.shape, def.element_type) {
        (Cardinality::PerRecord, Shape::Scalar, Some(ElementType::Varint)) => true, // variable
        (Cardinality::PerAllele, Shape::Scalar, Some(ElementType::Varint)) => true, // variable
        (Cardinality::PerRecord, Shape::Scalar, Some(e)) if e.fixed_byte_width().is_some() => {
            encoded.len() == n_records * e.fixed_byte_width().unwrap()
        }
        (Cardinality::PerAllele, Shape::Scalar, Some(e)) if e.fixed_byte_width().is_some() => {
            encoded.len() == n_total_alleles * e.fixed_byte_width().unwrap()
        }
        (_, Shape::List, _) => true, // variable
        (_, Shape::Bytes, _) => {
            // For tag 0x04, encoded length equals the sum of
            // allele_seq_len.
            encoded.len() as u64 == block.allele_seq_len.iter().sum::<u64>()
        }
        _ => true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::per_sample_caller::pileup::{AlleleObservation, AlleleSupportStats};
    use crate::per_sample_caller::psp::header::{
        ChromosomeEntry, ParameterValue, ParsedHeader, WriterProvenance, parse_header_bytes,
    };
    use crate::per_sample_caller::psp::index::decode_index;
    use crate::per_sample_caller::psp::trailer::{TRAILER_BYTES, decode_trailer};
    use std::collections::BTreeMap;
    use std::io::Cursor;

    // ---------- Fixture builders ---------------------------------

    fn writer_header(n_chroms: usize) -> WriterHeader {
        let chromosomes = (0..n_chroms)
            .map(|i| ChromosomeEntry {
                name: format!("chr{}", i + 1),
                length: 1_000_000,
                md5: "0".repeat(32),
            })
            .collect();
        let mut params = BTreeMap::new();
        params.insert("min-mapq".to_string(), ParameterValue::Integer(30));
        WriterHeader {
            format_version: (1, 0),
            sample: "sample".to_string(),
            reference: "ref.fa".to_string(),
            created: "2026-05-13T10:00:00Z".parse().unwrap(),
            chromosomes,
            writer: WriterProvenance {
                tool: "test".to_string(),
                version: "0.0.1".to_string(),
                subcommand: "per-sample".to_string(),
                input_crams: vec!["a.cram".to_string()],
                input_fasta: "ref.fa".to_string(),
                parameters: params,
            },
        }
    }

    fn support(num_obs: u32, q_sum: f64) -> AlleleSupportStats {
        AlleleSupportStats {
            num_obs,
            q_sum,
            fwd: num_obs / 2,
            placed_left: 0,
            placed_start: num_obs,
        }
    }

    fn allele(seq: &[u8], num_obs: u32, q_sum: f64, chain_slots: &[u16]) -> AlleleObservation {
        AlleleObservation {
            seq: seq.to_vec(),
            support: support(num_obs, q_sum),
            chain_slots: chain_slots.to_vec(),
        }
    }

    fn record(
        chrom_id: u32,
        pos: u32,
        alleles: Vec<AlleleObservation>,
        new_chains: Vec<u16>,
        expired_chains: Vec<u16>,
    ) -> PileupRecord {
        PileupRecord {
            chrom_id,
            pos,
            new_chains,
            expired_chains,
            alleles,
        }
    }

    // ---------- new() emits a valid framed header -----------------

    #[test]
    fn new_emits_valid_framed_header() {
        let header = writer_header(2);
        let cursor = Cursor::new(Vec::new());
        let writer = PspWriter::new(cursor, header.clone()).expect("new should succeed");
        // Drop the writer to extract the sink (or just call finish to get the trailer too).
        let sink = writer.finish().expect("finish should succeed");
        let bytes = sink.into_inner();
        // Parse the header back.
        let (parsed, consumed): (ParsedHeader, usize) =
            parse_header_bytes(&bytes).expect("header parse round-trip");
        assert_eq!(parsed.sample, "sample");
        assert_eq!(parsed.chromosomes.len(), 2);
        assert!(consumed > 0);
    }

    // ---------- finish() with zero records -----------------------

    #[test]
    fn finish_with_zero_records_writes_header_empty_index_trailer() {
        let header = writer_header(1);
        let cursor = Cursor::new(Vec::new());
        let writer = PspWriter::new(cursor, header).unwrap();
        let bytes = writer.finish().unwrap().into_inner();

        // Last 32 bytes are the trailer.
        let trailer_bytes: &[u8; TRAILER_BYTES] =
            bytes[bytes.len() - TRAILER_BYTES..].try_into().unwrap();
        let trailer = decode_trailer(trailer_bytes).expect("trailer should decode");
        assert_eq!(trailer.n_blocks, 0);
        assert_eq!(trailer.index_byte_length, 0);
        // Index lives between header-end and trailer-start, here
        // empty — index_offset points to trailer-start.
        let expected_index_offset = (bytes.len() - TRAILER_BYTES) as u64;
        assert_eq!(trailer.index_offset, expected_index_offset);

        let index = decode_index(&[], trailer.n_blocks).unwrap();
        assert!(index.is_empty());
    }

    // ---------- write_record validation ---------------------------

    #[test]
    fn rejects_unknown_chrom_id() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                5,
                100,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("chrom_id 5 with 1 chromosome should fail");
        assert!(matches!(err, PspWriteError::UnknownChromId { .. }));
    }

    #[test]
    fn rejects_pos_zero() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                0,
                0,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("pos 0 should fail");
        assert!(matches!(err, PspWriteError::PosOutOfRange { .. }));
    }

    #[test]
    fn rejects_pos_beyond_contig_length() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                0,
                10_000_000,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("pos beyond contig length should fail");
        assert!(matches!(err, PspWriteError::PosOutOfRange { .. }));
    }

    #[test]
    fn rejects_out_of_order_positions() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        let err = writer
            .write_record(&record(
                0,
                50,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("going backwards should fail");
        assert!(matches!(err, PspWriteError::OutOfOrderRecord { .. }));
    }

    #[test]
    fn rejects_chrom_regression() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(2)).unwrap();
        writer
            .write_record(&record(
                1,
                100,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        let err = writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("chrom regression should fail");
        assert!(matches!(err, PspWriteError::OutOfOrderRecord { .. }));
    }

    #[test]
    fn rejects_zero_alleles() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(0, 100, vec![], vec![], vec![]))
            .expect_err("zero alleles should fail");
        assert!(matches!(err, PspWriteError::InvalidRecord { .. }));
    }

    #[test]
    fn rejects_non_acgtn_allele_byte() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"X", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("non-ACGTN should fail");
        assert!(matches!(err, PspWriteError::InvalidRecord { .. }));
    }

    #[test]
    fn rejects_oversized_allele_sequence() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let seq = b"A".repeat((MAX_ALLELE_SEQ_LEN + 1) as usize);
        let err = writer
            .write_record(&record(
                0,
                100,
                vec![allele(&seq, 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("seq > cap should fail");
        assert!(matches!(err, PspWriteError::InvalidRecord { .. }));
    }

    #[test]
    fn rejects_nan_q_sum() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 10, f64::NAN, &[])],
                vec![],
                vec![],
            ))
            .expect_err("NaN q_sum should fail");
        assert!(matches!(err, PspWriteError::InvalidRecord { .. }));
    }

    #[test]
    fn rejects_phase_chain_marker_inconsistency_expired_not_active() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 1, -1.0, &[])],
                vec![],
                vec![5], // 5 is not in the active set yet
            ))
            .expect_err("expired-not-active should fail");
        assert!(matches!(
            err,
            PspWriteError::PhaseChainMarkerInconsistency { .. }
        ));
    }

    #[test]
    fn rejects_phase_chain_marker_inconsistency_double_open() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 1, -1.0, &[3])],
                vec![3],
                vec![],
            ))
            .unwrap();
        let err = writer
            .write_record(&record(
                0,
                101,
                vec![allele(b"A", 1, -1.0, &[3])],
                vec![3], // 3 is already active
                vec![],
            ))
            .expect_err("re-opening an active slot should fail");
        assert!(matches!(
            err,
            PspWriteError::PhaseChainMarkerInconsistency { .. }
        ));
    }

    #[test]
    fn rejects_allele_chain_slot_not_in_active_set() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 1, -1.0, &[99])], // 99 not opened
                vec![],
                vec![],
            ))
            .expect_err("allele references unknown slot");
        assert!(matches!(
            err,
            PspWriteError::PhaseChainMarkerInconsistency { .. }
        ));
    }

    // ---------- One block round-trip via writer alone ----------

    #[test]
    fn one_record_one_block_finishes_cleanly() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 10, -2.5, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        let bytes = writer.finish().unwrap().into_inner();

        // Trailer at tail says one block.
        let trailer: &[u8; TRAILER_BYTES] =
            bytes[bytes.len() - TRAILER_BYTES..].try_into().unwrap();
        let trailer = decode_trailer(trailer).unwrap();
        assert_eq!(trailer.n_blocks, 1);

        // Index decodes to one entry with the expected coordinates.
        let index_bytes = &bytes[trailer.index_offset as usize
            ..(trailer.index_offset + trailer.index_byte_length) as usize];
        let entries = decode_index(index_bytes, trailer.n_blocks).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].chrom_id, 0);
        assert_eq!(entries[0].first_pos, 100);
        assert_eq!(entries[0].last_pos, 100);
        assert_eq!(entries[0].n_records, 1);

        // Index checksum verifies.
        let stored = trailer.index_checksum;
        let computed = checksum_index(index_bytes);
        assert_eq!(stored, computed);
    }

    /// Chromosome boundary forces a flush.
    #[test]
    fn chrom_change_triggers_flush() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(2)).unwrap();
        writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 1, -1.0, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        writer
            .write_record(&record(
                0,
                200,
                vec![allele(b"C", 1, -1.0, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        // Switch chromosome — the writer should flush block 0 and
        // start block 1.
        writer
            .write_record(&record(
                1,
                1,
                vec![allele(b"G", 1, -1.0, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        let trailer: &[u8; TRAILER_BYTES] =
            bytes[bytes.len() - TRAILER_BYTES..].try_into().unwrap();
        let trailer = decode_trailer(trailer).unwrap();
        assert_eq!(
            trailer.n_blocks, 2,
            "chrom change should produce two blocks"
        );
    }

    /// Active-slot snapshot on subsequent (same-chrom) blocks
    /// matches the running set at flush time. Tested indirectly via
    /// a write_record sequence that closes block 1's last open slot
    /// just before block 2 opens.
    #[test]
    fn active_slot_snapshot_carries_across_blocks() {
        // Build a writer with a small target so a single record
        // doesn't flush; we'll force a flush by chromosome change.
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(2)).unwrap();
        // chrom 0: open slot 7, no per-allele reference; expire it.
        writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 1, -1.0, &[7])],
                vec![7],
                vec![],
            ))
            .unwrap();
        writer
            .write_record(&record(
                0,
                101,
                vec![allele(b"C", 1, -1.0, &[])],
                vec![],
                vec![7],
            ))
            .unwrap();
        // Move to chrom 1; active set is empty at this point.
        writer
            .write_record(&record(
                1,
                1,
                vec![allele(b"G", 1, -1.0, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        // Smoke test: file is non-zero and ends with a valid trailer.
        let trailer: &[u8; TRAILER_BYTES] =
            bytes[bytes.len() - TRAILER_BYTES..].try_into().unwrap();
        assert!(decode_trailer(trailer).is_ok());
    }

    /// Phase-chain still active when chromosome switches: writer
    /// catches the upstream bug at first-record time on the new
    /// chromosome.
    #[test]
    fn rejects_chain_active_at_chromosome_boundary() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(2)).unwrap();
        // Open a slot on chrom 0 and never close it.
        writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 1, -1.0, &[5])],
                vec![5],
                vec![],
            ))
            .unwrap();
        // Try to start chrom 1 — writer notices the active set is
        // non-empty and rejects.
        let err = writer
            .write_record(&record(
                1,
                1,
                vec![allele(b"G", 1, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("chain active across chrom boundary should fail");
        assert!(matches!(
            err,
            PspWriteError::PhaseChainMarkerInconsistency { .. }
        ));
    }
}
