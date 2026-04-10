//! This module provides an efficient streaming iterator for parsing GVCF files with
//! optimizations for common genomic data patterns.
//!
//! # Example
//!
//! ```no_run
//! use merge_per_sample_vcfs::gvcf_parser::VarIterator;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Parse a gzipped GVCF file
//! let parser = VarIterator::from_gzip_path("sample.g.vcf.gz")?;
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

use crate::decompression_pool::{DecompressionPool, PooledReader};
use crate::errors::VcfParseError;
use crate::threaded_reader::ThreadedReader;
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
/// Other FORMAT fields (PL, DP, GQ, etc.) are accessible lazily via
/// [`gt_field_index`](Variant::gt_field_index) and
/// [`get_gt_field_by_index`](Variant::get_gt_field_by_index).
#[derive(Debug)]
pub struct Variant {
    /// Chromosome/contig name
    pub chrom: String,
    /// Genomic position (1-based)
    pub pos: u32,
    /// Alleles: first is reference, rest are alternate alleles
    pub alleles: Vec<String>,
    /// Variant quality score (QUAL field). `NaN` if missing.
    pub qual: f32,
    /// Genotypes for all samples, stored as a flat vector.
    /// For sample `i` with ploidy `p`: `genotypes[i*p..(i+1)*p]`
    /// Missing genotypes are represented as `-1`
    pub genotypes: Vec<i8>,
    /// Phase information for each sample: `true` = phased (`|`), `false` = unphased (`/`)
    pub phases: Vec<bool>,
    /// FORMAT field names from the VCF line (e.g., `["GT", "PL", "DP", "GQ"]`).
    /// Empty for programmatically constructed variants.
    pub gt_format_fields: Vec<String>,
    /// Raw sample genotype field strings, one per sample.
    /// Each string contains colon-separated values matching `gt_format_fields`.
    /// Empty for programmatically constructed variants.
    pub sample_gt_fields: Vec<String>,
    pub n_samples: usize,
}

/// Parses the FORMAT field and returns (gt_index, field_names).
/// Returns an error if GT is not found.
fn parse_format_field(format: &str) -> VcfResult<(usize, Vec<String>)> {
    let fields: Vec<String> = format.split(':').map(String::from).collect();
    let gt_index = fields.iter().position(|f| f == "GT").ok_or_else(|| {
        VcfParseError::MissingGtFieldInFormat {
            line: format.to_string(),
        }
    })?;
    Ok((gt_index, fields))
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
    sample_gt_fields: &[String],
    gt_index: usize,
    ploidy: u8,
    patterns: &CommonGenotypePatterns,
) -> VcfResult<(Vec<i8>, Vec<bool>)> {
    let n_samples = sample_gt_fields.len();
    let total_num_alleles = n_samples * ploidy as usize;

    // Zero-fill the genotypes buffer, most alleles are reference alleles
    let mut genotypes: Vec<i8> = vec![0; total_num_alleles];
    let mut phases: Vec<bool> = vec![false; n_samples];

    // Parse each sample
    for (sample_idx, sample_field) in sample_gt_fields.iter().enumerate() {
        // Get GT field using nth() to avoid collecting into Vec
        let gt_str =
            sample_field
                .split(':')
                .nth(gt_index)
                .ok_or_else(|| VcfParseError::MissingGtField {
                    sample: sample_idx.to_string(),
                    line: sample_field.to_string(),
                })?;

        // Fast path: hardcoded diploid reference (most common case)
        if ploidy == 2 {
            if gt_str == "0/0" {
                continue;
            } else if gt_str == "0|0" {
                phases[sample_idx] = true;
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
                phases[sample_idx] = true;
                continue;
            } else if gt_str == patterns.missing_unphased {
                let start = sample_idx * ploidy as usize;
                for i in 0..ploidy as usize {
                    genotypes[start + i] = MISSING_ALLELE;
                }
                phases[sample_idx] = false;
                continue;
            } else if gt_str == patterns.missing_phased {
                let start = sample_idx * ploidy as usize;
                for i in 0..ploidy as usize {
                    genotypes[start + i] = MISSING_ALLELE;
                }
                phases[sample_idx] = true;
                continue;
            }
        }

        if gt_str == "." {
            let start = sample_idx * ploidy as usize;
            for i in 0..ploidy as usize {
                genotypes[start + i] = MISSING_ALLELE;
            }
            phases[sample_idx] = false;
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

        // Copy alleles to genotypes vector
        let start = sample_idx * ploidy as usize;
        for (i, &allele) in alleles.iter().enumerate() {
            genotypes[start + i] = allele;
        }
        phases[sample_idx] = phased;
    }

    Ok((genotypes, phases))
}

impl Variant {
    /// Returns the end position of the variant's span (1-based, inclusive).
    pub fn get_var_end(&self) -> VcfResult<u32> {
        let ref_allele = self.alleles.first().ok_or_else(|| VcfParseError::RuntimeError {
            message: format!(
                "Variant at {}:{} has no alleles",
                self.chrom, self.pos
            ),
        })?;
        let ref_len: u32 = ref_allele.len().try_into().map_err(|_| VcfParseError::RuntimeError {
            message: format!(
                "Reference allele length overflow at {}:{}",
                self.chrom, self.pos
            ),
        })?;
        if ref_len == 0 {
            return Err(VcfParseError::RuntimeError {
                message: format!(
                    "Empty reference allele at {}:{}",
                    self.chrom, self.pos
                ),
            });
        }
        self.pos
            .checked_add(ref_len - 1)
            .ok_or_else(|| VcfParseError::RuntimeError {
                message: format!(
                    "Position overflow computing variant end at {}:{}",
                    self.chrom, self.pos
                ),
            })
    }

    /// Creates a new Variant with the given fields.
    /// Sets gt_format_fields and sample_gt_fields to empty (no lazy field access).
    pub fn new(
        chrom: String,
        pos: u32,
        alleles: Vec<String>,
        genotypes: Vec<i8>,
        phases: Vec<bool>,
        n_samples: usize,
    ) -> Self {
        Variant {
            chrom,
            pos,
            alleles,
            qual: f32::NAN,
            genotypes,
            phases,
            gt_format_fields: Vec::new(),
            sample_gt_fields: Vec::new(),
            n_samples,
        }
    }

    fn from_line(
        line: &str,
        gt_format_cache: &mut LruCache<String, (usize, Vec<String>), RandomState>,
        n_samples: usize,
        ploidy: u8,
        patterns: &CommonGenotypePatterns,
    ) -> VcfResult<Self> {
        let mut fields = line.splitn(10, '\t');
        let chrom = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields { line: line.to_string() })?;

        let pos = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields { line: line.to_string() })?;
        let _id = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields { line: line.to_string() })?;
        let ref_allele = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields { line: line.to_string() })?;
        let alt_alleles = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields { line: line.to_string() })?;
        let qual_str = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields { line: line.to_string() })?;
        let _filter = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields { line: line.to_string() })?;
        let _info = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields { line: line.to_string() })?;
        let gt_format_str = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields { line: line.to_string() })?;

        // Get sample fields (the 10th field from splitn contains all remaining sample data)
        let sample_fields = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields { line: line.to_string() })?;

        let pos: u32 = pos.parse().map_err(|_| VcfParseError::InvalidPosition {
            value: pos.to_string(),
            line: line.to_string(),
        })?;

        let qual: f32 = if qual_str == "." {
            f32::NAN
        } else {
            qual_str
                .parse()
                .map_err(|_| VcfParseError::InvalidQuality {
                    value: qual_str.to_string(),
                    line: line.to_string(),
                })?
        };

        let alleles: Vec<String> = std::iter::once(ref_allele)
            .chain(alt_alleles.split(','))
            .filter(|allele| allele != &NON_REF && *allele != ".")
            .map(str::to_string)
            .collect();

        // Check LRU cache first
        let (gt_index, gt_format_fields) =
            if let Some((gt_idx, fields)) = gt_format_cache.get(gt_format_str) {
                (*gt_idx, fields.clone())
            } else {
                let (gt_idx, fields) = parse_format_field(gt_format_str)?;
                gt_format_cache.put(gt_format_str.to_string(), (gt_idx, fields.clone()));
                (gt_idx, fields)
            };

        // Split sample fields once, reused for both GT parsing and lazy field access
        let sample_gt_fields: Vec<String> = sample_fields.split('\t').map(String::from).collect();
        if sample_gt_fields.len() != n_samples {
            return Err(VcfParseError::MalformedLine {
                reason: format!(
                    "expected {} samples, found {}",
                    n_samples,
                    sample_gt_fields.len()
                ),
                line: line.to_string(),
            });
        }

        // Parse genotypes for all samples
        let (genotypes, phases) = parse_genotypes(&sample_gt_fields, gt_index, ploidy, patterns)?;

        Ok(Variant {
            chrom: chrom.to_string(),
            pos,
            alleles,
            qual,
            genotypes,
            phases,
            gt_format_fields,
            sample_gt_fields,
            n_samples,
        })
    }

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
    // FORMAT Field Access
    // ============================================================================

    /// Returns the index of a FORMAT field by name, or `None` if not present.
    ///
    /// The returned index can be passed to [`get_gt_field_by_index`](Variant::get_gt_field_by_index)
    /// for efficient repeated access across variants with the same FORMAT layout.
    pub fn gt_field_index(&self, field: &str) -> Option<usize> {
        self.gt_format_fields.iter().position(|f| f == field)
    }

    /// Returns the raw string value for a FORMAT field (by index) for each sample.
    ///
    /// Returns `"."` for samples where the field is missing.
    pub fn get_gt_field_by_index(&self, idx: usize) -> Vec<&str> {
        self.sample_gt_fields
            .iter()
            .map(|s| s.split(':').nth(idx).unwrap_or("."))
            .collect()
    }
}

/// Streaming iterator for parsing GVCF files.
pub struct VarIterator<B: BufRead> {
    /// Underlying buffered reader
    reader: B,
    /// Reusable line buffer
    line: String,
    /// Internal buffer of parsed records
    vars_buffer: VecDeque<Variant>,
    /// LRU cache for FORMAT field parsing: maps format_string -> (gt_index, field_names).
    gt_format_cache: LruCache<String, (usize, Vec<String>), RandomState>,
    /// Sample names extracted from the #CHROM header line
    samples: Vec<String>,
    ploidy: u8,
    /// Pre-computed patterns for fast-path genotype matching
    common_gt_patterns: CommonGenotypePatterns,
}

impl<B: BufRead> VarIterator<B> {
    fn new(mut reader: B) -> VcfResult<Self> {
        let mut line = String::new();

        // Read the first line to start header processing
        reader.read_line(&mut line).map_err(VcfParseError::from)?;

        let mut iterator = VarIterator {
            reader,
            line,
            vars_buffer: VecDeque::new(),
            gt_format_cache: LruCache::with_hasher(
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
            return Err(VcfParseError::BrokenHeader {
                line: self.line.trim().to_string(),
            });
        }

        let fields: Vec<&str> = self.line.trim().split('\t').collect();

        // VCF spec requires standard columns plus at least 1 sample
        if fields.len() < VCF_STANDARD_COLUMNS + 1 {
            return Err(VcfParseError::BrokenHeader {
                line: self.line.trim().to_string(),
            });
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
                    Ok(0) => {
                        return Err(VcfParseError::BrokenHeader {
                            line: "EOF reached before #CHROM header line was found"
                                .to_string(),
                        });
                    }
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

    fn fill_buffer(&mut self, n_items: usize) -> VcfResult<usize> {
        let mut n_items_added: usize = 0;
        let n_samples = self.samples.len();

        while self.vars_buffer.len() < n_items {
            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => break, // EOF
                Ok(_) => {
                    let record = Variant::from_line(
                        &self.line,
                        &mut self.gt_format_cache,
                        n_samples,
                        self.ploidy,
                        &self.common_gt_patterns,
                    )?;
                    self.vars_buffer.push_back(record);
                    n_items_added += 1;
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
    pub fn peek_variant(&mut self) -> VcfResult<Option<&Variant>> {
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

impl<R: Read> VarIterator<BufReader<R>> {
    pub fn from_reader(reader: R) -> VcfResult<Self> {
        let buf_reader = BufReader::with_capacity(BUFREADER_CAPACITY, reader);
        VarIterator::new(buf_reader)
    }
}
impl<R: Read> VarIterator<BufReader<MultiGzDecoder<R>>> {
    pub fn from_gzip_reader(reader: R) -> VcfResult<Self> {
        let gz_decoder = MultiGzDecoder::new(reader);
        let buf_reader = BufReader::with_capacity(BUFREADER_CAPACITY, gz_decoder);
        VarIterator::new(buf_reader)
    }
}
/// Validates that a file is gzipped and returns a `MultiGzDecoder` for it.
fn open_gzip_file<P: AsRef<Path>>(path: P) -> VcfResult<MultiGzDecoder<File>> {
    if !file_is_gzipped(&path).map_err(|e| VcfParseError::MagicByteError {
        path: path.as_ref().to_string_lossy().to_string(),
        reason: e.to_string(),
    })? {
        return Err(VcfParseError::VCFFileShouldBeGzipped {
            path: path.as_ref().to_string_lossy().to_string(),
        });
    }
    let file = File::open(&path)?;
    Ok(MultiGzDecoder::new(file))
}

impl VarIterator<BufReader<MultiGzDecoder<File>>> {
    pub fn from_gzip_path<P: AsRef<Path>>(path: P) -> VcfResult<Self> {
        let decoder = open_gzip_file(path)?;
        let buf_reader = BufReader::with_capacity(BUFREADER_CAPACITY, decoder);
        VarIterator::new(buf_reader)
    }
}
impl VarIterator<BufReader<ThreadedReader>> {
    pub fn from_gzip_path_threaded<P: AsRef<Path>>(path: P) -> VcfResult<Self> {
        let decoder = open_gzip_file(path)?;
        let threaded = ThreadedReader::new(decoder);
        let buf_reader = BufReader::with_capacity(BUFREADER_CAPACITY, threaded);
        VarIterator::new(buf_reader)
    }
}
impl VarIterator<BufReader<PooledReader>> {
    /// Create a `VarIterator` that uses a shared `DecompressionPool` for
    /// gzip decompression, instead of spawning a dedicated thread per file.
    pub fn from_gzip_path_pooled<P: AsRef<Path>>(
        path: P,
        pool: &DecompressionPool,
    ) -> VcfResult<Self> {
        let decoder = open_gzip_file(path)?;
        let pooled = pool.register(decoder);
        let buf_reader = BufReader::with_capacity(BUFREADER_CAPACITY, pooled);
        VarIterator::new(buf_reader)
    }
}

impl<R: BufRead> Iterator for VarIterator<R> {
    type Item = VcfResult<Variant>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.vars_buffer.is_empty() {
            match self.fill_buffer(DEF_N_VARIANTS_IN_BUFFER) {
                Err(error) => return Some(Err(error)),
                Ok(n_items_added) => match n_items_added {
                    0 => return None,
                    _ => (),
                },
            }
        }

        if let Some(variant) = self.vars_buffer.pop_front() {
            Some(Ok(variant))
        } else {
            unreachable!("Buffer was just filled or the no more variants branch was dealt with")
        }
    }
}
