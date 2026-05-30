//! Chunk-driven stream of [`PerPositionPileups`] for the
//! contamination side-pass.
//!
//! [`estimate_contamination`](crate::var_calling::contamination_estimation::estimate_contamination)
//! consumes any
//! `Iterator<Item = Result<PerPositionPileups, _>>`. Before this
//! module landed, that iterator came from
//! [`PerPositionMerger`](crate::var_calling::per_position_merger::PerPositionMerger)
//! reading per-sample PSP record iterators directly. The stream
//! here builds the same iterator shape on top of the chunk loader's
//! columnar output:
//!
//! 1. For each chromosome, walk its BP range one chunk at a time
//!    (driven by `chunk_genomic_span` + `target_variants_per_chunk`).
//! 2. After each `load_chunk_from_iters` call, take the union of
//!    `chunk.per_sample[s].positions` — these are the cohort-wide
//!    kept positions for this chunk (the loader's variant filter
//!    already dropped monomorphic ones).
//! 3. For each kept position, build a [`PerPositionPileups`] by
//!    looking up the per-sample record (or `None`) at that position
//!    and materialising it.
//!
//! **Why this is byte-identical to the streaming pipeline.** The
//! streaming `PerPositionMerger` emits every position any sample
//! has a record at; the side-pass's step 1a filter rejects every
//! position that is monomorphic in the cohort (one distinct allele
//! type only). The chunk loader's variant filter drops positions
//! that are monomorphic in the cohort, with the
//! grouping-simulation expansion keeping homref positions inside
//! MNP / DEL / INS reach — the same positions the streaming
//! pipeline emits. Step 1a still runs on the chunk-emitted
//! positions and still rejects every position that has only one
//! cohort allele type, so the set of positions that contribute
//! observations to the side-pass is identical.

use std::io::{Read, Seek};

use thiserror::Error;

use crate::pileup_record::PileupRecord;
use crate::psp::header::ParsedChromosome;
use crate::psp::{PspReadError, PspReader};
use crate::var_calling::cohort_block::columns::{MaterialisedChunk, SampleColumns};
use crate::var_calling::cohort_block::loader::{
    ChunkLoadError, ChunkLoadExtent, ChunkLoadScratch, load_chunk_from_iters,
};
use crate::var_calling::per_position_merger::PerPositionPileups;

/// Errors surfaced by [`ChunkedPositionStream`]. Wraps the chunk
/// loader's error type so the side-pass sees a single surface.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum ContaminationStreamError {
    /// Surfaced when the chunk loader fails on a per-sample PSP
    /// iterator or detects a cohort-shape invariant violation.
    #[error("chunk load: {0}")]
    ChunkLoad(#[from] ChunkLoadError<PspReadError>),
}

/// Chunk-driven stream of cohort-wide per-position pileups.
///
/// The stream owns the PSP readers + the chunk loader scratch +
/// the materialised chunk. Each `next()` either returns the next
/// kept position's pileups or advances to the next chunk (and, if
/// needed, the next chromosome). Memory peak is bounded by one
/// chunk's worth of columnar data.
pub struct ChunkedPositionStream<R: Read + Seek> {
    psp_readers: Vec<PspReader<R>>,
    chromosomes: Vec<ParsedChromosome>,
    /// Cohort size; equals `psp_readers.len()`.
    n_samples: usize,
    /// Index into `chromosomes` of the chromosome the next chunk
    /// will load from. Increments when the current chromosome's
    /// BP range is exhausted.
    chrom_idx: usize,
    /// BP position (1-based, inclusive) where the next chunk on
    /// the current chromosome will start.
    chrom_cursor: u32,
    /// Loader scratch — reused across every chunk in every
    /// chromosome.
    chunk_scratch: ChunkLoadScratch,
    /// Current loaded chunk.
    chunk: MaterialisedChunk,
    /// Empty per-sample carryover. The side-pass does not use
    /// carryover semantics; the loader signature requires it so
    /// we pass empty slots and never touch them.
    carryover: Vec<SampleColumns>,
    /// Cohort-wide sorted+dedup union of `chunk.per_sample[s].positions`
    /// for the current chunk. Each entry corresponds to one
    /// `PerPositionPileups` the stream will emit.
    chunk_positions: Vec<u32>,
    /// Cursor into `chunk_positions`. Advances on every emitted
    /// pileup; reset to `0` on each chunk load.
    chunk_position_cursor: usize,
    /// First-attempt BP span passed to the loader.
    chunk_genomic_span: u32,
    /// Soft lower bound on the post-filter variant count per chunk.
    /// `0` disables the loader's extension loop.
    target_variants_per_chunk: u32,
    /// Latched after every chromosome is consumed; subsequent
    /// `next()` calls short-circuit to `None`.
    is_exhausted: bool,
}

impl<R: Read + Seek + Send> ChunkedPositionStream<R> {
    /// Build a stream over `psp_readers` covering every chromosome
    /// in `chromosomes`. The readers are taken by value because the
    /// stream needs `&mut` access to each one across every chunk.
    pub fn new(
        psp_readers: Vec<PspReader<R>>,
        chromosomes: Vec<ParsedChromosome>,
        chunk_genomic_span: u32,
        target_variants_per_chunk: u32,
    ) -> Self {
        let n_samples = psp_readers.len();
        let chunk_scratch = ChunkLoadScratch::with_n_samples(n_samples);
        let chunk = MaterialisedChunk::with_n_samples(n_samples);
        let carryover: Vec<SampleColumns> =
            (0..n_samples).map(|_| SampleColumns::empty()).collect();
        Self {
            psp_readers,
            chromosomes,
            n_samples,
            chrom_idx: 0,
            chrom_cursor: 1,
            chunk_scratch,
            chunk,
            carryover,
            chunk_positions: Vec::new(),
            chunk_position_cursor: 0,
            chunk_genomic_span,
            target_variants_per_chunk,
            is_exhausted: false,
        }
    }

    /// Load chunks until one with at least one kept position is
    /// found, or the entire cohort is exhausted. Returns `Ok(true)`
    /// when `self.chunk_positions` is non-empty and the cursor was
    /// reset to `0`; `Ok(false)` when the stream is fully drained.
    fn advance_to_next_chunk(&mut self) -> Result<bool, ContaminationStreamError> {
        loop {
            if self.chrom_idx >= self.chromosomes.len() {
                self.is_exhausted = true;
                return Ok(false);
            }
            let chrom = &self.chromosomes[self.chrom_idx];
            let chrom_id = self.chrom_idx as u32;
            let chrom_length = chrom.length;
            let chrom_one_past_end = chrom_length.saturating_add(1);

            if chrom_length == 0 || self.chrom_cursor >= chrom_one_past_end {
                // Current chromosome consumed; advance.
                self.chrom_idx += 1;
                self.chrom_cursor = 1;
                continue;
            }

            let range_start = self.chrom_cursor;
            // First-attempt span; the loader extends up to max_span
            // when `target_variants_per_chunk > 0`.
            let initial_span = self.chunk_genomic_span.max(1);
            let remaining_span = chrom_one_past_end.saturating_sub(range_start);
            let max_span = remaining_span.max(initial_span);
            // Open PSP iterators for the full remaining chromosome
            // range so the loader's internal extension loop can keep
            // pulling beyond `initial_span` without exhausting the
            // iterator early.
            let psp_inclusive_end = chrom_length;

            let iters: Vec<_> = self
                .psp_readers
                .iter_mut()
                .map(|r| r.region_records(chrom_id, range_start, psp_inclusive_end))
                .collect();

            let _stats = load_chunk_from_iters(
                &mut self.chunk_scratch,
                &mut self.chunk,
                ChunkLoadExtent {
                    chrom_id,
                    range_start,
                    initial_span,
                    target_variants: self.target_variants_per_chunk,
                    max_span,
                },
                iters,
                &mut self.carryover,
            )?;

            // The loader sets `chunk.range = range_start..attempt_end`
            // where `attempt_end` is the BP boundary the loader
            // actually stopped at. Advance the cursor past it.
            self.chrom_cursor = self.chunk.range.end;

            // Build the chunk's cohort-wide kept-position list as
            // the sorted+dedup union of every sample's positions.
            // The loader has already filtered out positions whose
            // groups carry no variant signal, so this is exactly the
            // set the side-pass needs.
            self.chunk_positions.clear();
            for sample in &self.chunk.per_sample {
                self.chunk_positions
                    .extend_from_slice(sample.positions.as_slice());
            }
            self.chunk_positions.sort_unstable();
            self.chunk_positions.dedup();
            self.chunk_position_cursor = 0;

            if !self.chunk_positions.is_empty() {
                return Ok(true);
            }
            // Chunk had zero kept positions. Loop and load the next
            // one. The loop terminates because `chrom_cursor`
            // advances on every iteration (the loader always pulls
            // at least one BP of range).
        }
    }
}

impl<R: Read + Seek + Send> Iterator for ChunkedPositionStream<R> {
    type Item = Result<PerPositionPileups, ContaminationStreamError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_exhausted {
            return None;
        }
        if self.chunk_position_cursor >= self.chunk_positions.len() {
            match self.advance_to_next_chunk() {
                Ok(true) => {}
                Ok(false) => return None,
                Err(e) => return Some(Err(e)),
            }
        }
        let pos = self.chunk_positions[self.chunk_position_cursor];
        self.chunk_position_cursor += 1;
        let chrom_id = self.chunk.chrom_id;

        let mut per_sample: Vec<Option<PileupRecord>> = Vec::with_capacity(self.n_samples);
        for sample in &self.chunk.per_sample {
            let entry = match sample.binary_search_position(pos) {
                Ok(row_idx) => Some(sample.materialise_record(chrom_id, row_idx)),
                Err(_) => None,
            };
            per_sample.push(entry);
        }
        Some(Ok(PerPositionPileups {
            chrom_id,
            pos,
            per_sample,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup_record::{AlleleObservation, AlleleSupportStats};
    use crate::psp::header::{ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance};
    use crate::psp::writer::PspWriter;
    use std::collections::BTreeMap;
    use std::io::Cursor;

    /// Per-sample input bundle. Records must be sorted by
    /// `(chrom_id, pos)` ascending; the PSP writer rejects
    /// out-of-order input.
    #[derive(Clone)]
    struct SampleFixture {
        sample_name: &'static str,
        records: Vec<PileupRecord>,
    }

    fn record(
        chrom_id: u32,
        pos: u32,
        ref_seq: &[u8],
        alt_seq: &[u8],
        ref_obs: u32,
        alt_obs: u32,
    ) -> PileupRecord {
        PileupRecord::new(
            chrom_id,
            pos,
            vec![
                AlleleObservation::new(
                    ref_seq.to_vec(),
                    AlleleSupportStats::new(ref_obs, 0.0, ref_obs, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    alt_seq.to_vec(),
                    AlleleSupportStats::new(alt_obs, 0.0, alt_obs, 0, 0, 0, 0),
                    Vec::new(),
                ),
            ],
        )
    }

    fn ref_only(chrom_id: u32, pos: u32, ref_obs: u32) -> PileupRecord {
        PileupRecord::new(
            chrom_id,
            pos,
            vec![AlleleObservation::new(
                b"A".to_vec(),
                AlleleSupportStats::new(ref_obs, 0.0, ref_obs, 0, 0, 0, 0),
                Vec::new(),
            )],
        )
    }

    fn chrom(name: &str, length: u32) -> ParsedChromosome {
        ParsedChromosome {
            name: name.to_string(),
            length,
            md5: "0".repeat(32),
        }
    }

    /// Build an in-memory PSP file holding `records` for `sample_name`
    /// against `chromosomes`. Returns the byte buffer ready to be
    /// wrapped in a `PspReader<Cursor<Vec<u8>>>`.
    fn make_psp_bytes(
        sample_name: &str,
        chromosomes: &[ParsedChromosome],
        records: Vec<PileupRecord>,
    ) -> Vec<u8> {
        let mut params = BTreeMap::new();
        params.insert("min-mapq".to_string(), ParameterValue::Integer(30));
        let header = WriterHeader {
            format_version: (1, 0),
            sample: sample_name.to_string(),
            reference: "ref.fa".to_string(),
            created: "2026-05-13T10:00:00Z".parse().unwrap(),
            chromosomes: chromosomes
                .iter()
                .map(|c| ChromosomeEntry {
                    name: c.name.clone(),
                    length: c.length,
                    md5: c.md5.clone(),
                })
                .collect(),
            writer: WriterProvenance {
                tool: "test".to_string(),
                version: "0.0.1".to_string(),
                subcommand: "per-sample".to_string(),
                input_crams: vec!["a.cram".to_string()],
                input_fasta: "ref.fa".to_string(),
                parameters: params,
            },
        };
        let cursor = Cursor::new(Vec::new());
        let mut writer = PspWriter::new(cursor, header).expect("PspWriter::new");
        for rec in &records {
            writer.write_record(rec).expect("write_record");
        }
        let sink = writer.finish().expect("finish");
        sink.into_inner()
    }

    fn build_stream(
        samples: Vec<SampleFixture>,
        chromosomes: Vec<ParsedChromosome>,
        chunk_genomic_span: u32,
        target_variants_per_chunk: u32,
    ) -> ChunkedPositionStream<Cursor<Vec<u8>>> {
        let psp_readers: Vec<PspReader<Cursor<Vec<u8>>>> = samples
            .into_iter()
            .map(|s| {
                let bytes = make_psp_bytes(s.sample_name, &chromosomes, s.records);
                PspReader::new(Cursor::new(bytes)).expect("PspReader::new")
            })
            .collect();
        ChunkedPositionStream::new(
            psp_readers,
            chromosomes,
            chunk_genomic_span,
            target_variants_per_chunk,
        )
    }

    #[test]
    fn empty_cohort_yields_no_positions() {
        let chromosomes = vec![chrom("chr1", 1000)];
        let samples = vec![
            SampleFixture {
                sample_name: "s0",
                records: Vec::new(),
            },
            SampleFixture {
                sample_name: "s1",
                records: Vec::new(),
            },
        ];
        let mut stream = build_stream(samples, chromosomes, 100, 0);
        assert!(stream.next().is_none());
    }

    #[test]
    fn single_chunk_emits_variant_positions_only() {
        // Cohort has ALT at pos 10 (sample 0) and pos 30 (sample 1).
        // Pos 20 is pure REF in the cohort — filtered out.
        let chromosomes = vec![chrom("chr1", 100)];
        let samples = vec![
            SampleFixture {
                sample_name: "s0",
                records: vec![record(0, 10, b"A", b"T", 3, 4), ref_only(0, 20, 5)],
            },
            SampleFixture {
                sample_name: "s1",
                records: vec![ref_only(0, 20, 5), record(0, 30, b"A", b"G", 3, 4)],
            },
        ];
        let mut stream = build_stream(samples, chromosomes, 200, 0);
        let positions: Vec<u32> = std::iter::from_fn(|| stream.next())
            .map(|r| r.expect("stream OK").pos)
            .collect();
        assert_eq!(positions, vec![10, 30]);
    }

    #[test]
    fn multi_chunk_spans_chromosome_in_order() {
        // Variants at pos 10, 50, 90 on a 100-bp chromosome with
        // chunk_genomic_span=30. Expected chunks: [1..31), [31..61),
        // [61..91), [91..101). Each contains one variant.
        let chromosomes = vec![chrom("chr1", 100)];
        let samples = vec![SampleFixture {
            sample_name: "s0",
            records: vec![
                record(0, 10, b"A", b"T", 3, 4),
                record(0, 50, b"A", b"T", 3, 4),
                record(0, 90, b"A", b"T", 3, 4),
            ],
        }];
        let mut stream = build_stream(samples, chromosomes, 30, 0);
        let positions: Vec<u32> = std::iter::from_fn(|| stream.next())
            .map(|r| r.expect("stream OK").pos)
            .collect();
        assert_eq!(positions, vec![10, 50, 90]);
    }

    #[test]
    fn multi_chromosome_emits_in_chromosome_order() {
        let chromosomes = vec![chrom("chr1", 100), chrom("chr2", 100)];
        let samples = vec![SampleFixture {
            sample_name: "s0",
            records: vec![
                record(0, 20, b"A", b"T", 3, 4),
                record(1, 30, b"A", b"T", 3, 4),
                record(1, 80, b"A", b"T", 3, 4),
            ],
        }];
        let mut stream = build_stream(samples, chromosomes, 200, 0);
        let pairs: Vec<(u32, u32)> = std::iter::from_fn(|| stream.next())
            .map(|r| {
                let p = r.expect("stream OK");
                (p.chrom_id, p.pos)
            })
            .collect();
        assert_eq!(pairs, vec![(0, 20), (1, 30), (1, 80)]);
    }

    #[test]
    fn per_sample_some_and_none_at_each_position() {
        // Sample 0 has a record at pos 10 only; sample 1 has a
        // record at pos 20 only. Both positions are variant.
        let chromosomes = vec![chrom("chr1", 100)];
        let samples = vec![
            SampleFixture {
                sample_name: "s0",
                records: vec![record(0, 10, b"A", b"T", 3, 4)],
            },
            SampleFixture {
                sample_name: "s1",
                records: vec![record(0, 20, b"A", b"T", 3, 4)],
            },
        ];
        let mut stream = build_stream(samples, chromosomes, 200, 0);

        let p10 = stream.next().expect("first emit").expect("OK");
        assert_eq!(p10.pos, 10);
        assert!(p10.per_sample[0].is_some());
        assert!(p10.per_sample[1].is_none());

        let p20 = stream.next().expect("second emit").expect("OK");
        assert_eq!(p20.pos, 20);
        assert!(p20.per_sample[0].is_none());
        assert!(p20.per_sample[1].is_some());

        assert!(stream.next().is_none());
    }

    #[test]
    fn target_variants_extends_chunk_until_target_reached() {
        // Variants every 100bp from pos 10 to pos 910. With
        // initial_span=50 and target=4, the loader should extend the
        // first chunk until at least 4 kept positions are accumulated.
        let chromosomes = vec![chrom("chr1", 1000)];
        let positions: Vec<u32> = (1..=10).map(|i| i * 100 - 90).collect();
        let s0_records: Vec<PileupRecord> = positions
            .iter()
            .map(|&p| record(0, p, b"A", b"T", 3, 4))
            .collect();
        let samples = vec![SampleFixture {
            sample_name: "s0",
            records: s0_records,
        }];
        let mut stream = build_stream(samples, chromosomes, 50, 4);
        let emitted: Vec<u32> = std::iter::from_fn(|| stream.next())
            .map(|r| r.expect("OK").pos)
            .collect();
        assert_eq!(emitted, positions);
    }
}
