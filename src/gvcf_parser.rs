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
//! }
//! # Ok(())
//! # }
//! ```

use crate::errors::VcfParseError;
use crate::utils_magic::file_is_gzipped;
use ahash::RandomState; // fast hasher
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

/// Pre-computed string patterns for common genotypes at a given ploidy.
///
/// These patterns are generated once and reused to enable fast-path matching
/// in `parse_genotypes` without hardcoding ploidy-specific strings.
struct CommonGenotypePatterns {
    ref_unphased: String,
    ref_phased: String,
    missing_unphased: String,
    missing_phased: String,
}

impl CommonGenotypePatterns {
    fn new(ploidy: u8) -> Self {
        let n = ploidy as usize;
        let ref_unphased = (0..n).map(|_| "0").collect::<Vec<_>>().join("/");
        let ref_phased = (0..n).map(|_| "0").collect::<Vec<_>>().join("|");
        let missing_unphased = (0..n).map(|_| ".").collect::<Vec<_>>().join("/");
        let missing_phased = (0..n).map(|_| ".").collect::<Vec<_>>().join("|");
        Self {
            ref_unphased,
            ref_phased,
            missing_unphased,
            missing_phased,
        }
    }
}

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
    pub n_samples: usize,
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
fn parse_genotype(gt_str: &str, expected_ploidy: u8) -> VcfResult<(Vec<i8>, bool)> {
    let bytes = gt_str.as_bytes();

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

/// Parses genotypes for all samples in a VCF line.
fn parse_genotypes(
    sample_fields: &str,
    n_samples: usize,
    gt_index: usize,
    ploidy: u8,
    patterns: &CommonGenotypePatterns,
) -> VcfResult<(Vec<i8>, Vec<bool>)> {
    let total_num_alleles = n_samples * ploidy as usize;

    // Zero-fill the genotypes buffer, most alleles are reference alleles
    let mut genotypes: Vec<i8> = vec![0; total_num_alleles];
    let mut phase: Vec<bool> = vec![false; total_num_alleles];

    // Parse each sample
    let mut last_idx: usize = 0;
    for (sample_idx, sample_field) in sample_fields.split('\t').enumerate() {
        last_idx = sample_idx;
        if sample_idx >= n_samples {
            return Err(VcfParseError::RuntimeError {
                message: format!(
                    "More than the expected number of samples {n_samples} found in line"
                ),
            });
        }

        // Get GT field using nth() to avoid collecting into Vec
        let gt_str = match sample_field.split(':').nth(gt_index) {
            Some(s) => s,
            None => {
                // Missing GT field - mark as missing
                let start = sample_idx * ploidy as usize;
                for i in 0..ploidy as usize {
                    genotypes[start + i] = MISSING_ALLELE;
                }
                phase[sample_idx] = false;
                continue;
            }
        };

        // Fast path: hardcoded diploid reference (most common case)
        if ploidy == 2 {
            if gt_str == "0/0" {
                continue;
            } else if gt_str == "0|0" {
                phase[sample_idx] = true;
                continue;
            } else if gt_str == "./." {
                let start = sample_idx * ploidy as usize;
                genotypes[start] = MISSING_ALLELE;
                genotypes[start + 1] = MISSING_ALLELE;
                continue;
            }
        } else {
            // Fast path: check for common patterns (any ploidy)
            if gt_str == patterns.ref_unphased {
                continue;
            } else if gt_str == patterns.ref_phased {
                phase[sample_idx] = true;
                continue;
            } else if gt_str == patterns.missing_unphased {
                let start = sample_idx * ploidy as usize;
                for i in 0..ploidy as usize {
                    genotypes[start + i] = MISSING_ALLELE;
                }
                phase[sample_idx] = false;
                continue;
            } else if gt_str == patterns.missing_phased {
                let start = sample_idx * ploidy as usize;
                for i in 0..ploidy as usize {
                    genotypes[start + i] = MISSING_ALLELE;
                }
                phase[sample_idx] = true;
                continue;
            }
        }

        if gt_str == "." {
            let start = sample_idx * ploidy as usize;
            for i in 0..ploidy as usize {
                genotypes[start + i] = MISSING_ALLELE;
            }
            phase[sample_idx] = false;
            continue;
        }

        // General path: parse genotype
        let (alleles, phased) = parse_genotype(gt_str, ploidy)?;

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
            genotypes[start + i] = allele;
        }
        phase[sample_idx] = phased;
    }

    let samples_seen = last_idx + 1;
    if samples_seen > n_samples {
        return Err(VcfParseError::RuntimeError {
            message: format!(
                "More samples than expected: expected {n_samples}, got {samples_seen}"
            ),
        });
    }

    Ok((genotypes, phase))
}

impl GVcfRecord {
    fn from_line(
        line: &str,
        format_cache: &mut LruCache<String, (usize, Option<usize>), RandomState>,
        n_samples: usize,
        ploidy: u8,
        patterns: &CommonGenotypePatterns,
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
        let (genotypes, phase) =
            parse_genotypes(sample_fields, n_samples, gt_index, ploidy, patterns)?;

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
            n_samples,
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
    // Metadata
    // ============================================================================

    /// Returns the number of samples in this record.
    #[inline]
    pub fn n_samples(&self) -> usize {
        self.n_samples
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
    vars_buffer: VecDeque<GVcfRecord>,
    /// LRU cache for FORMAT field parsing: maps format_string -> (gt_index, pl_index).
    format_cache: LruCache<String, (usize, Option<usize>), RandomState>,
    /// Sample names extracted from the #CHROM header line
    samples: Vec<String>,
    ploidy: u8,
    /// Pre-computed patterns for fast-path genotype matching
    common_gt_patterns: CommonGenotypePatterns,
}

impl<B: BufRead> GVcfRecordIterator<B> {
    fn new(mut reader: B) -> VcfResult<Self> {
        let mut line = String::new();

        // Read the first line to start header processing
        reader.read_line(&mut line).map_err(VcfParseError::from)?;

        let mut iterator = GVcfRecordIterator {
            reader,
            line,
            vars_buffer: VecDeque::new(),
            format_cache: LruCache::with_hasher(
                NonZeroUsize::new(GT_FORMAT_LRU_CACHE_SIZE).unwrap(),
                RandomState::new(),
            ),
            samples: Vec::new(),
            ploidy: DEFAULT_PLOIDY,
            common_gt_patterns: CommonGenotypePatterns::new(DEFAULT_PLOIDY),
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

        while self.vars_buffer.len() < n_items {
            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => break, // EOF
                Ok(_) => {
                    match GVcfRecord::from_line(
                        &self.line,
                        &mut self.format_cache,
                        n_samples,
                        self.ploidy,
                        &self.common_gt_patterns,
                    ) {
                        Ok(record) => {
                            self.vars_buffer.push_back(record);
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
        if self.vars_buffer.is_empty() {
            self.fill_buffer(DEF_N_VARIANTS_IN_BUFFER)?;
        }
        Ok(self.vars_buffer.front())
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
        if self.vars_buffer.len() == 0 {
            match self.fill_buffer(DEF_N_VARIANTS_IN_BUFFER) {
                Err(error) => return Some(Err(error)),
                Ok(n_items_added) => match n_items_added {
                    0 => return None,
                    _ => (),
                },
            }
        }

        if let Some(variant) = self.vars_buffer.pop_front() {
            return Some(Ok(variant));
        } else {
            return Some(Err(VcfParseError::RuntimeError {
                message: "The buffer should contain something".to_string(),
            }));
        }
    }
}
