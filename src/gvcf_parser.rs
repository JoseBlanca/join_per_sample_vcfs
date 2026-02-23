//! High-performance parser for GVCF (Genomic VCF) files.
//!
//! This module provides an efficient streaming iterator for parsing GVCF files with
//! optimizations for common genomic data patterns.
//!
//! # Features
//!
//! - **Zero-allocation fast paths**: Common genotypes like `0/0` and `0|0` are handled
//!   without allocations
//! - **Reusable buffers**: Genotype and phase buffers are reused across variants to
//!   minimize allocations
//! - **Large I/O buffers**: 256KB buffer size for better throughput on large files
//! - **LRU caching**: FORMAT field indices are cached to avoid repeated parsing
//! - **Header processing**: Samples are extracted and available immediately after
//!   iterator construction
//!
//! # Example
//!
//! ```no_run
//! use join_per_sample_vcfs::gvcf_parser::GVcfRecordIterator;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Parse a gzipped GVCF file
//! let parser = GVcfRecordIterator::from_gzip_path("sample.g.vcf.gz")?;
//!
//! // Access sample names
//! let samples = parser.samples();
//! println!("Found {} samples", samples.len());
//!
//! // Iterate over variants
//! for result in parser {
//!     let record = result?;
//!
//!     // Access variant information
//!     println!("{}:{} {} -> {:?}",
//!         record.chrom, record.pos,
//!         record.alleles[0], &record.alleles[1..]);
//!
//!     // Access genotypes for each sample
//!     for sample_idx in 0..record.n_samples() {
//!         let genotypes = record.sample_genotypes(sample_idx);
//!         let is_missing = record.is_missing(sample_idx);
//!         let is_phased = record.is_phased(sample_idx);
//!
//!         if !is_missing {
//!             println!("  Sample {}: {:?} (phased: {})",
//!                 sample_idx, genotypes, is_phased);
//!         }
//!     }
//! }
//! # Ok(())
//! # }
//! ```

use crate::errors::VcfParseError;
use crate::utils_magic::file_is_gzipped;
use flate2::read::MultiGzDecoder;
use lru::LruCache;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::num::NonZeroUsize;
use std::path::Path;

/// Result type alias for VCF parsing operations.
pub type VcfResult<T> = std::result::Result<T, VcfParseError>;

/// Special ALT allele indicating non-reference positions in GVCF files.
const NON_REF: &str = "<NON_REF>";

/// Default number of variants to buffer for efficient iteration.
const DEF_N_VARIANTS_IN_BUFFER: usize = 100;

/// LRU cache size for FORMAT field parsing. Typically sufficient for VCF files
/// with a small number of alternating FORMAT patterns.
const GT_FORMAT_LRU_CACHE_SIZE: usize = 5;

/// BufReader capacity for I/O operations (256KB).
/// Larger buffer reduces system calls and improves throughput on large files.
const BUFREADER_CAPACITY: usize = 256 * 1024;

/// Number of standard VCF columns before sample data begins.
/// Standard columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT.
const VCF_STANDARD_COLUMNS: usize = 9;

/// Default ploidy assumption for variants (diploid).
const DEFAULT_PLOIDY: u8 = 2;

/// Marker value for missing alleles in genotype data.
const MISSING_ALLELE: i8 = -1;

/// A single variant record from a GVCF file.
///
/// # Genotype Storage
///
/// Genotypes are stored as a flat `Vec<i8>` where each sample's alleles are consecutive.
/// For a sample at index `i` with ploidy `p`:
/// - Alleles are at indices `[i*p .. (i+1)*p)`
/// - Missing alleles are represented as `-1`
/// - For diploid (p=2): `genotypes[i*2]` and `genotypes[i*2+1]` are the two alleles
///
#[derive(Debug)]
pub struct GVcfRecord {
    /// Chromosome/contig name
    pub chrom: String,
    /// Genomic position (1-based)
    pub pos: u32,
    /// Alleles: first is reference, rest are alternate alleles
    pub alleles: Vec<String>,
    /// Length of the reference allele in base pairs
    pub ref_allele_len: u8,
    /// Variant quality score (QUAL field). `NaN` if missing.
    pub qual: f32,
    /// Index of GT field in FORMAT (cached for internal use)
    gt_index: usize,
    /// Index of PL field in FORMAT if present (cached for internal use)
    pl_index: Option<usize>,
    /// Genotypes for all samples, stored as a flat vector.
    /// For sample `i` with ploidy `p`: `genotypes[i*p..(i+1)*p]`
    /// Missing genotypes are represented as `-1`
    pub genotypes: Vec<i8>,
    /// Phase information for each sample: `true` = phased (`|`), `false` = unphased (`/`)
    pub phase: Vec<bool>,
    /// Ploidy of this variant (typically `2` for diploid organisms)
    pub ploidy: u8,
}

/// Parses the FORMAT field and returns (gt_index, pl_index).
/// Returns an error if GT is not found. PL is optional.
fn parse_format_field(format: &str) -> VcfResult<(usize, Option<usize>)> {
    let mut gt_index = None;
    let mut pl_index = None;

    for (idx, field) in format.split(':').enumerate() {
        if field == "GT" {
            gt_index = Some(idx);
        } else if field == "PL" {
            pl_index = Some(idx);
        }

        // Early exit if we found both
        if gt_index.is_some() && pl_index.is_some() {
            break;
        }
    }

    let gt_idx = gt_index.ok_or(VcfParseError::MissingGtFieldInFormat {
        line: format.to_string(),
    })?;

    Ok((gt_idx, pl_index))
}

/// Fast genotype parser optimized for common cases.
/// Returns (alleles, is_phased) where alleles contains parsed allele indices (-1 for missing).
#[inline]
fn parse_genotype_fast(gt_str: &str, expected_ploidy: u8) -> VcfResult<(Vec<i8>, bool)> {
    let bytes = gt_str.as_bytes();
    let len = bytes.len();

    // Fast path: check for common diploid patterns first
    if expected_ploidy == 2 && len == 3 {
        match (bytes[0], bytes[1], bytes[2]) {
            (b'0', b'/', b'0') => return Ok((vec![0, 0], false)),
            (b'0', b'|', b'0') => return Ok((vec![0, 0], true)),
            (b'.', b'/', b'.') => return Ok((vec![MISSING_ALLELE, MISSING_ALLELE], false)),
            (b'.', b'|', b'.') => return Ok((vec![MISSING_ALLELE, MISSING_ALLELE], true)),
            _ => {}
        }
    }

    // Fast path: haploid reference
    if expected_ploidy == 1 && len == 1 && bytes[0] == b'0' {
        return Ok((vec![0], false));
    }

    // General parser for other cases
    let mut alleles = Vec::with_capacity(expected_ploidy as usize);
    let mut phased = false;
    let mut current: i8 = 0;
    let mut parsing_number = false;
    let mut seen_dot = false;

    for &byte in bytes {
        match byte {
            b'0'..=b'9' => {
                current = (current as u8)
                    .checked_mul(10)
                    .and_then(|v| v.checked_add(byte - b'0'))
                    .ok_or_else(|| VcfParseError::AlleleIndexOverflow {
                        gt: gt_str.to_string(),
                    })? as i8;
                parsing_number = true;
                seen_dot = false;
            }
            b'.' => {
                alleles.push(MISSING_ALLELE);
                seen_dot = true;
                parsing_number = false;
            }
            b'|' => {
                if parsing_number {
                    alleles.push(current);
                    current = 0;
                    parsing_number = false;
                }
                phased = true;
            }
            b'/' => {
                if parsing_number {
                    alleles.push(current);
                    current = 0;
                    parsing_number = false;
                } else if seen_dot {
                    // Handle cases like "./." or ".|."
                    seen_dot = false;
                }
            }
            _ => {}
        }
    }

    // Push the last allele if we were parsing a number
    if parsing_number {
        alleles.push(current);
    }

    Ok((alleles, phased))
}

/// Validates that the number of sample fields matches the expected count.
fn validate_sample_count(sample_fields: &str, expected: usize) -> VcfResult<()> {
    let actual = sample_fields.split('\t').count();
    if actual != expected {
        return Err(VcfParseError::RuntimeError {
            message: format!("Expected {} samples, found {}", expected, actual),
        });
    }
    Ok(())
}

/// Detects the ploidy of a variant from sample genotypes.
///
/// Scans through samples to find the first informative genotype (not 0/0, 0|0, or .)
/// and uses its allele count to determine ploidy. Falls back to `last_ploidy` if
/// no informative genotypes are found.
fn detect_ploidy(sample_fields: &str, gt_index: usize, last_ploidy: u8) -> u8 {
    for sample in sample_fields.split('\t') {
        if let Some(gt_str) = sample.split(':').nth(gt_index) {
            // Skip genotypes that don't help with ploidy detection:
            // - 0/0 and 0|0 (common reference genotypes)
            // - "." (missing with unknown ploidy)
            if matches!(gt_str, "0/0" | "0|0" | ".") {
                continue;
            }

            if let Ok((alleles, _)) = parse_genotype_fast(gt_str, last_ploidy) {
                if !alleles.is_empty() {
                    return alleles.len() as u8;
                }
            }
        }
    }
    last_ploidy
}

/// Parses genotypes for all samples in a VCF line.
fn parse_genotypes(
    sample_fields: &str,
    n_samples: usize,
    gt_index: usize,
    last_ploidy: u8,
    genotypes_buffer: &mut Vec<i8>,
    phase_buffer: &mut Vec<bool>,
) -> VcfResult<(Vec<i8>, Vec<bool>, u8)> {
    validate_sample_count(sample_fields, n_samples)?;

    // Detect ploidy from first informative genotype
    let ploidy = detect_ploidy(sample_fields, gt_index, last_ploidy);
    let total_alleles = n_samples * ploidy as usize;

    // Resize buffers if needed
    if genotypes_buffer.len() < total_alleles {
        genotypes_buffer.resize(total_alleles, 0);
    }
    if phase_buffer.len() < n_samples {
        phase_buffer.resize(n_samples, false);
    }

    // Zero-fill the genotypes buffer, most alleles are reference alleles
    genotypes_buffer[0..total_alleles].fill(0);

    // Parse each sample
    for (sample_idx, sample_field) in sample_fields.split('\t').enumerate() {
        // Get GT field using nth() to avoid collecting into Vec
        let gt_str = match sample_field.split(':').nth(gt_index) {
            Some(s) => s,
            None => {
                // Missing GT field - mark as missing
                let start = sample_idx * ploidy as usize;
                for i in 0..ploidy as usize {
                    genotypes_buffer[start + i] = MISSING_ALLELE;
                }
                phase_buffer[sample_idx] = false;
                continue;
            }
        };

        // Fast path: check for common patterns
        if ploidy == 2 {
            match gt_str {
                "0/0" => {
                    phase_buffer[sample_idx] = false;
                    // Already zero-filled, skip
                    continue;
                }
                "0|0" => {
                    phase_buffer[sample_idx] = true;
                    // Already zero-filled, skip
                    continue;
                }
                "./." => {
                    genotypes_buffer[sample_idx * 2] = MISSING_ALLELE;
                    genotypes_buffer[sample_idx * 2 + 1] = MISSING_ALLELE;
                    phase_buffer[sample_idx] = false;
                    continue;
                }
                ".|." => {
                    genotypes_buffer[sample_idx * 2] = MISSING_ALLELE;
                    genotypes_buffer[sample_idx * 2 + 1] = MISSING_ALLELE;
                    phase_buffer[sample_idx] = true;
                    continue;
                }
                "." => {
                    // Single "." means missing with unknown ploidy - treat as diploid missing
                    genotypes_buffer[sample_idx * 2] = MISSING_ALLELE;
                    genotypes_buffer[sample_idx * 2 + 1] = MISSING_ALLELE;
                    phase_buffer[sample_idx] = false;
                    continue;
                }
                _ => {}
            }
        } else if ploidy == 1 && gt_str == "." {
            // Haploid missing
            genotypes_buffer[sample_idx] = MISSING_ALLELE;
            phase_buffer[sample_idx] = false;
            continue;
        }

        // General path: parse genotype
        let (alleles, phased) = parse_genotype_fast(gt_str, ploidy)?;

        // Validate ploidy consistency
        if alleles.len() != ploidy as usize {
            return Err(VcfParseError::RuntimeError {
                message: format!(
                    "Inconsistent ploidy: expected {}, got {} in sample {}",
                    ploidy,
                    alleles.len(),
                    sample_idx
                ),
            });
        }

        // Copy alleles to buffer
        let start = sample_idx * ploidy as usize;
        for (i, &allele) in alleles.iter().enumerate() {
            genotypes_buffer[start + i] = allele;
        }
        phase_buffer[sample_idx] = phased;
    }

    // Return owned copies of the exact size needed
    let genotypes = genotypes_buffer[0..total_alleles].to_vec();
    let phase = phase_buffer[0..n_samples].to_vec();

    Ok((genotypes, phase, ploidy))
}

impl GVcfRecord {
    fn from_line(
        line: &str,
        format_cache: &mut LruCache<String, (usize, Option<usize>)>,
        n_samples: usize,
        last_ploidy: u8,
        genotypes_buffer: &mut Vec<i8>,
        phase_buffer: &mut Vec<bool>,
    ) -> VcfResult<Self> {
        let mut fields = line.splitn(10, '\t');
        let chrom = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields)?;

        let pos = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields)?;
        fields.next(); // ID
        let ref_allele = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields)?;
        let alt_alleles = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields)?;
        let qual_str = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields)?;
        fields.next(); // FILTER
        fields.next(); // INFO
        let format_str = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields)?;

        // Get sample fields (the 10th field from splitn contains all remaining sample data)
        let sample_fields = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields)?;

        let pos = pos
            .parse::<u32>()
            .map_err(|_| VcfParseError::GVCFLineNotEnoughFields)?;

        let ref_allele_len = ref_allele.len() as u8;

        let qual = if qual_str == "." {
            f32::NAN
        } else {
            qual_str
                .parse::<f32>()
                .map_err(|_| VcfParseError::GVCFLineNotEnoughFields)?
        };

        let alleles: Vec<String> = std::iter::once(ref_allele)
            .chain(alt_alleles.split(','))
            .filter(|allele| allele != &NON_REF)
            .map(str::to_string)
            .collect();

        // Check LRU cache first
        let (gt_index, pl_index) = if let Some((gt_idx, pl_idx)) = format_cache.get(format_str) {
            // Cache hit: return cached indices
            (*gt_idx, *pl_idx)
        } else {
            // Cache miss: parse and insert into cache
            let (gt_idx, pl_idx) = parse_format_field(format_str)?;
            format_cache.put(format_str.to_string(), (gt_idx, pl_idx));
            (gt_idx, pl_idx)
        };

        // Parse genotypes for all samples
        let (genotypes, phase, ploidy) = parse_genotypes(
            sample_fields,
            n_samples,
            gt_index,
            last_ploidy,
            genotypes_buffer,
            phase_buffer,
        )?;

        Ok(GVcfRecord {
            chrom: chrom.to_string(),
            pos,
            alleles,
            ref_allele_len,
            qual,
            gt_index,
            pl_index,
            genotypes,
            phase,
            ploidy,
        })
    }

    // ============================================================================
    // Core Accessors
    // ============================================================================

    /// Returns the chromosome/contig name.
    #[inline]
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    /// Returns the genomic position (1-based).
    #[inline]
    pub fn position(&self) -> u32 {
        self.pos
    }

    /// Returns all alleles (reference first, then alternates).
    #[inline]
    pub fn alleles(&self) -> &[String] {
        &self.alleles
    }

    /// Returns the variant quality score. May be `NaN` if not available.
    #[inline]
    pub fn quality(&self) -> f32 {
        self.qual
    }

    // ============================================================================
    // Genotype Accessors
    // ============================================================================

    /// Returns the genotypes for a specific sample as a slice.
    ///
    /// For a diploid sample, this returns a slice of length 2 containing the two allele indices.
    /// Missing genotypes are represented as `-1`.
    ///
    /// # Panics
    ///
    /// Panics if `sample_idx >= self.n_samples()`.
    #[inline]
    pub fn sample_genotypes(&self, sample_idx: usize) -> &[i8] {
        let start = sample_idx * self.ploidy as usize;
        let end = start + self.ploidy as usize;
        &self.genotypes[start..end]
    }

    /// Returns whether the specified sample is phased.
    ///
    /// Phased genotypes use `'|'` separator, unphased use `'/'`.
    #[inline]
    pub fn is_phased(&self, sample_idx: usize) -> bool {
        self.phase[sample_idx]
    }

    /// Returns whether the specified sample has a missing genotype.
    ///
    /// A genotype is missing if any allele is `-1` (represented as `'.'` in VCF).
    #[inline]
    #[must_use]
    pub fn is_missing(&self, sample_idx: usize) -> bool {
        self.sample_genotypes(sample_idx)
            .iter()
            .any(|&a| a == MISSING_ALLELE)
    }

    /// Returns whether the specified sample is homozygous reference (all alleles are 0).
    #[inline]
    #[must_use]
    pub fn is_hom_ref(&self, sample_idx: usize) -> bool {
        self.sample_genotypes(sample_idx).iter().all(|&a| a == 0)
    }

    /// Returns whether the specified sample is homozygous (all alleles are the same).
    #[inline]
    #[must_use]
    pub fn is_homozygous(&self, sample_idx: usize) -> bool {
        let gts = self.sample_genotypes(sample_idx);
        if gts.is_empty() {
            return false;
        }
        let first = gts[0];
        gts.iter().all(|&a| a == first)
    }

    // ============================================================================
    // Metadata
    // ============================================================================

    /// Returns the number of samples in this record.
    #[inline]
    pub fn n_samples(&self) -> usize {
        self.phase.len()
    }

    /// Returns the ploidy of this variant.
    #[inline]
    pub fn ploidy(&self) -> u8 {
        self.ploidy
    }

    // ============================================================================
    // Utility Methods
    // ============================================================================

    /// Returns the genomic span covered by this variant as (start, end) inclusive.
    ///
    /// For SNPs, start == end. For indels, the span covers the maximum allele length.
    pub fn get_span(&self) -> VcfResult<(u32, u32)> {
        let max_allele_len = self.alleles.iter().map(|allele| allele.len()).max().ok_or(
            VcfParseError::RuntimeError {
                message: "There should be at least one allele".to_string(),
            },
        )?;
        if max_allele_len == 1 {
            Ok((self.pos, self.pos))
        } else {
            Ok((self.pos, self.pos + max_allele_len as u32 - 1))
        }
    }

    // ============================================================================
    // Internal/Testing Methods
    // ============================================================================

    /// Returns the index of GT in the FORMAT field (for internal use and testing).
    #[doc(hidden)]
    pub fn gt_index(&self) -> usize {
        self.gt_index
    }

    /// Returns the index of PL in the FORMAT field if present (for internal use and testing).
    #[doc(hidden)]
    pub fn pl_index(&self) -> Option<usize> {
        self.pl_index
    }
}

/// Streaming iterator for parsing GVCF files.
pub struct GVcfRecordIterator<B: BufRead> {
    /// Underlying buffered reader
    reader: B,
    /// Reusable line buffer
    line: String,
    /// Internal buffer of parsed records
    buffer: VecDeque<GVcfRecord>,
    /// LRU cache for FORMAT field parsing: maps format_string -> (gt_index, pl_index).
    format_cache: LruCache<String, (usize, Option<usize>)>,
    /// Sample names extracted from the #CHROM header line
    samples: Vec<String>,
    /// Reusable buffer for genotypes to minimize allocations
    genotypes_buffer: Vec<i8>,
    /// Reusable buffer for phase information
    phase_buffer: Vec<bool>,
    /// Last seen ploidy to optimize buffer allocation
    last_seen_ploidy: u8,
}

impl<B: BufRead> GVcfRecordIterator<B> {
    fn new(mut reader: B) -> VcfResult<Self> {
        let mut line = String::new();

        // Read the first line to start header processing
        reader.read_line(&mut line).map_err(VcfParseError::from)?;

        let mut iterator = GVcfRecordIterator {
            reader,
            line,
            buffer: VecDeque::new(),
            format_cache: LruCache::new(NonZeroUsize::new(GT_FORMAT_LRU_CACHE_SIZE).unwrap()),
            samples: Vec::new(),
            genotypes_buffer: Vec::new(),
            phase_buffer: Vec::new(),
            last_seen_ploidy: DEFAULT_PLOIDY,
        };

        // Process the header and populate samples
        iterator.process_header()?;

        Ok(iterator)
    }

    /// Parses sample names from the #CHROM header line.
    ///
    /// The #CHROM line has 9 standard columns followed by sample names:
    /// `#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [sample1] [sample2] ...`
    fn parse_samples_from_header(&self) -> VcfResult<Vec<String>> {
        // Verify this is the #CHROM line
        if !self.line.starts_with("#CHROM") {
            return Err(VcfParseError::BrokenHeader);
        }

        let fields: Vec<&str> = self.line.trim().split('\t').collect();

        // VCF spec requires standard columns plus at least 1 sample
        if fields.len() < VCF_STANDARD_COLUMNS + 1 {
            return Err(VcfParseError::BrokenHeader);
        }

        let samples: Vec<String> = fields
            .iter()
            .skip(VCF_STANDARD_COLUMNS)
            .map(|s| s.to_string())
            .collect();

        Ok(samples)
    }

    /// Processes the VCF header, parsing sample names and preparing for variant reading.
    /// This should be called during construction to initialize the iterator.
    fn process_header(&mut self) -> VcfResult<()> {
        // Read through ## header lines until we find the #CHROM line
        loop {
            if self.line.starts_with("##") {
                self.line.clear();
                match self.reader.read_line(&mut self.line) {
                    Ok(0) => return Err(VcfParseError::BrokenHeader),
                    Ok(_) => {
                        if !self.line.starts_with("##") {
                            break;
                        }
                    }
                    Err(error) => return Err(VcfParseError::from(error)),
                }
            } else {
                break;
            }
        }

        // At this point, self.line contains the #CHROM line - parse samples
        self.samples = self.parse_samples_from_header()?;

        Ok(())
    }

    pub fn fill_buffer(&mut self, n_items: usize) -> VcfResult<usize> {
        let mut n_items_added: usize = 0;
        let n_samples = self.samples.len();

        while self.buffer.len() < n_items {
            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => break, // EOF
                Ok(_) => {
                    match GVcfRecord::from_line(
                        &self.line,
                        &mut self.format_cache,
                        n_samples,
                        self.last_seen_ploidy,
                        &mut self.genotypes_buffer,
                        &mut self.phase_buffer,
                    ) {
                        Ok(record) => {
                            // Update last seen ploidy for next record
                            self.last_seen_ploidy = record.ploidy;
                            self.buffer.push_back(record);
                            n_items_added += 1;
                        }
                        Err(err) => return Err(err),
                    }
                }
                Err(err) => {
                    return Err(VcfParseError::from(err));
                }
            }
        }
        Ok(n_items_added)
    }

    /// Returns a reference to the next variant without consuming it.
    ///
    /// If the buffer is empty, it will attempt to fill it first.
    /// Returns `None` if the iterator is exhausted (no more variants).
    pub fn peek_variant(&mut self) -> VcfResult<Option<&GVcfRecord>> {
        if self.buffer.is_empty() {
            self.fill_buffer(DEF_N_VARIANTS_IN_BUFFER)?;
        }
        Ok(self.buffer.front())
    }

    /// Returns the sample names from the VCF header.
    ///
    /// Samples are parsed from the #CHROM line during construction and are
    /// always available.
    pub fn samples(&self) -> &[String] {
        &self.samples
    }
}

impl<R: Read> GVcfRecordIterator<BufReader<R>> {
    pub fn from_reader(reader: R) -> VcfResult<Self> {
        let buf_reader = BufReader::with_capacity(BUFREADER_CAPACITY, reader);
        Ok(GVcfRecordIterator::new(buf_reader)?)
    }
}
impl<R: Read> GVcfRecordIterator<BufReader<MultiGzDecoder<R>>> {
    pub fn from_gzip_reader(reader: R) -> VcfResult<Self> {
        let gz_decoder = MultiGzDecoder::new(reader);
        let buf_reader = BufReader::with_capacity(BUFREADER_CAPACITY, gz_decoder);
        Ok(GVcfRecordIterator::new(buf_reader)?)
    }
}
impl GVcfRecordIterator<BufReader<MultiGzDecoder<File>>> {
    pub fn from_gzip_path<P: AsRef<Path>>(path: P) -> VcfResult<Self> {
        if !file_is_gzipped(&path).map_err(|_| VcfParseError::MagicByteError)? {
            return Err(VcfParseError::VCFFileShouldBeGzipped);
        }
        let file = File::open(&path)?;
        let gz_decoder = MultiGzDecoder::new(file);
        let buf_reader = BufReader::with_capacity(BUFREADER_CAPACITY, gz_decoder);
        Ok(GVcfRecordIterator::new(buf_reader)?)
    }
}

impl<R: BufRead> Iterator for GVcfRecordIterator<R> {
    type Item = VcfResult<GVcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.buffer.len() == 0 {
            match self.fill_buffer(DEF_N_VARIANTS_IN_BUFFER) {
                Err(error) => return Some(Err(error)),
                Ok(n_items_added) => match n_items_added {
                    0 => return None,
                    _ => (),
                },
            }
        }

        if let Some(variant) = self.buffer.pop_front() {
            return Some(Ok(variant));
        } else {
            return Some(Err(VcfParseError::RuntimeError {
                message: "The buffer should contain something".to_string(),
            }));
        }
    }
}
