use crate::errors::VcfParseError;
use crate::utils_magic::file_is_gzipped;
use flate2::read::MultiGzDecoder;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

pub type VcfResult<T> = std::result::Result<T, VcfParseError>;

const NON_REF: &str = "<NON_REF>";
const DEF_N_VARIANTS_IN_BUFFER: usize = 100;

#[derive(Debug)]
pub struct GVcfRecord {
    pub chrom: String,
    pub pos: u32,
    pub alleles: Vec<String>,
    pub ref_allele_len: u8,
    pub qual: f32,
    gt_index: usize,
    pl_index: Option<usize>,
}

/// Parses the FORMAT field and returns (gt_index, pl_index).
/// Returns an error if GT is not found. PL is optional.
fn parse_format_field(format: &str) -> VcfResult<(usize, Option<usize>)> {
    let fields: Vec<&str> = format.split(':').collect();

    let gt_index = fields.iter().position(|&field| field == "GT").ok_or(
        VcfParseError::MissingGtFieldInFormat {
            line: format.to_string(),
        },
    )?;

    let pl_index = fields.iter().position(|&field| field == "PL");

    Ok((gt_index, pl_index))
}

impl GVcfRecord {
    fn from_line(
        line: &str,
        format_cache: &mut Box<Option<(String, usize, Option<usize>)>>,
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

        if alt_alleles == NON_REF {
            return Err(VcfParseError::InvariantgVCFLine);
        }

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

        // Check cache first
        let (gt_index, pl_index) = if let Some((cached_format, cached_gt_index, cached_pl_index)) =
            format_cache.as_ref()
        {
            if format_str == cached_format.as_str() {
                (*cached_gt_index, *cached_pl_index)
            } else {
                // Format changed, parse and update cache
                let (gt_idx, pl_idx) = parse_format_field(format_str)?;
                **format_cache = Some((format_str.to_string(), gt_idx, pl_idx));
                (gt_idx, pl_idx)
            }
        } else {
            // First time, parse and cache
            let (gt_idx, pl_idx) = parse_format_field(format_str)?;
            **format_cache = Some((format_str.to_string(), gt_idx, pl_idx));
            (gt_idx, pl_idx)
        };

        Ok(GVcfRecord {
            chrom: chrom.to_string(),
            pos: pos,
            alleles: alleles,
            ref_allele_len: ref_allele_len,
            qual: qual,
            gt_index: gt_index,
            pl_index: pl_index,
        })
    }

    /// Returns the index of GT in the FORMAT field (for internal use and testing)
    pub fn gt_index(&self) -> usize {
        self.gt_index
    }

    /// Returns the index of PL in the FORMAT field if present (for internal use and testing)
    pub fn pl_index(&self) -> Option<usize> {
        self.pl_index
    }

    pub fn get_span(self: &GVcfRecord) -> VcfResult<(u32, u32)> {
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
}

#[derive(Debug, PartialEq, Eq)]
enum VcfSection {
    Header,
    Body,
}

pub struct GVcfRecordIterator<B: BufRead> {
    reader: B,
    line: String,
    section: VcfSection,
    buffer: VecDeque<GVcfRecord>,
    /// Cache for FORMAT field: (format_string, gt_index, pl_index)
    /// Boxed to allow mutation from from_line
    format_cache: Box<Option<(String, usize, Option<usize>)>>,
}

impl<B: BufRead> GVcfRecordIterator<B> {
    fn new(reader: B) -> Self {
        GVcfRecordIterator {
            reader: reader,
            line: String::new(),
            section: VcfSection::Header,
            buffer: VecDeque::new(),
            format_cache: Box::new(None),
        }
    }
    fn process_header_and_first_variant(&mut self) -> Option<VcfResult<GVcfRecord>> {
        loop {
            if self.line.starts_with("##") {
                self.line.clear();
                match self.reader.read_line(&mut self.line) {
                    Ok(0) => return Some(Err(VcfParseError::BrokenHeader)),
                    Ok(_) => {
                        if !self.line.starts_with("##") {
                            break;
                        }
                    }
                    Err(error) => return Some(Err(VcfParseError::from(error))),
                }
            }
        }
        self.line.clear();
        match self.reader.read_line(&mut self.line) {
            Ok(0) => None, // EOF
            Ok(_) => {
                self.section = VcfSection::Body;
                Some(GVcfRecord::from_line(&self.line, &mut self.format_cache))
            }
            Err(error) => return Some(Err(VcfParseError::from(error))),
        }
    }
    pub fn fill_buffer(&mut self, n_items: usize) -> VcfResult<usize> {
        let mut n_items_added: usize = 0;
        while self.buffer.len() < n_items {
            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => break, // EOF
                Ok(_) => {
                    if self.section == VcfSection::Header {
                        let result = self.process_header_and_first_variant();
                        if let Some(Ok(record)) = result {
                            self.buffer.push_back(record);
                            n_items_added += 1;
                        }
                    } else {
                        match GVcfRecord::from_line(&self.line, &mut self.format_cache) {
                            Ok(record) => {
                                self.buffer.push_back(record);
                                n_items_added += 1;
                            }
                            Err(VcfParseError::InvariantgVCFLine) => continue, // skip
                            Err(err) => return Err(err),
                        }
                    }
                }
                Err(err) => {
                    return Err(VcfParseError::from(err));
                }
            }
        }
        Ok(n_items_added)
    }

    pub fn peek_items_in_buffer(&self) -> impl Iterator<Item = &GVcfRecord> {
        self.buffer.iter()
    }
}

impl<R: Read> GVcfRecordIterator<BufReader<R>> {
    pub fn from_reader(reader: R) -> Self {
        let buf_reader = BufReader::new(reader);
        GVcfRecordIterator::new(buf_reader)
    }
}
impl<R: Read> GVcfRecordIterator<BufReader<MultiGzDecoder<R>>> {
    pub fn from_gzip_reader(reader: R) -> Self {
        let gz_decoder = MultiGzDecoder::new(reader);
        let buf_reader = BufReader::new(gz_decoder);
        GVcfRecordIterator::new(buf_reader)
    }
}
impl GVcfRecordIterator<BufReader<MultiGzDecoder<File>>> {
    pub fn from_gzip_path<P: AsRef<Path>>(path: P) -> VcfResult<Self> {
        if !file_is_gzipped(&path).map_err(|_| VcfParseError::MagicByteError)? {
            return Err(VcfParseError::VCFFileShouldBeGzipped);
        }
        let file = File::open(&path)?;
        let gz_decoder = MultiGzDecoder::new(file);
        let buf_reader = BufReader::new(gz_decoder);
        Ok(GVcfRecordIterator::new(buf_reader))
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
