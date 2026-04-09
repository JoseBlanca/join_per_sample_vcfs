use std::collections::HashSet;
use std::io::BufRead;

use rayon::prelude::*;

use crate::errors::VcfParseError;
use crate::gvcf_parser::{Variant, VarIterator, VcfResult};

/// Result of peeking a single iterator during the parallel scan.
enum PeekResult {
    /// Iterator is exhausted.
    Exhausted,
    /// Variant is on a different chromosome (not the current one).
    DifferentChrom,
    /// Variant is on the current chromosome.
    OnCurrentChrom {
        pos: u32,
        span_end: u32,
        is_variable: bool,
    },
    /// Ordering violation detected.
    Error(VcfParseError),
}

/// Metadata about a `VariantIterator`
#[derive(Debug, Clone)]
pub struct VariantIteratorInfo {
    pub samples: Vec<String>,
}

/// A bin of overlapping variants collected from multiple per-sample VCFs.
///
/// The bin covers a contiguous genomic region where all contained variants
/// have reference-allele spans that overlap with the bin's span.
#[derive(Debug)]
pub struct OverlappingVariantGroup {
    pub chrom: String,
    /// Start position of the bin (1-based, inclusive).
    pub start: u32,
    /// End position of the bin (1-based, inclusive).
    pub end: u32,
    /// All variant records in this bin, in consumption order.
    pub variants: Vec<Variant>,
    /// Index into `VariantGroupIterator::iter_info` for each variant,
    /// parallel to `variants`.
    pub source_var_iter_idxs: Vec<usize>,
}

impl OverlappingVariantGroup {
    pub fn span(&self) -> (&str, u32, u32) {
        (&self.chrom, self.start, self.end)
    }
}

/// Computes the reference-allele span end for a variant (1-based, inclusive).
#[inline]
fn ref_span_end(record: &Variant) -> u32 {
    record.pos + record.ref_allele_len as u32 - 1
}

/// Iterator that merges multiple per-sample VCF streams, grouping
/// overlapping variants into bins based on reference-allele span overlap.
///
/// Variants from all input VCFs are consumed in genomic order. Two variants
/// overlap when they are on the same chromosome and one's start position
/// falls within the other's reference-allele span.
pub struct VariantGroupIterator<B: BufRead> {
    vcf_iters: Vec<VarIterator<B>>,
    iter_info: Vec<VariantIteratorInfo>,
    sorted_chromosomes: Vec<String>,
    current_chrom_idx: usize,
    chroms_seen: HashSet<String>,
    last_group_start: Option<u32>,
    done: bool,
}

impl<B: BufRead + Send> VariantGroupIterator<B> {
    /// Creates a new `VariantGroupIterator`.
    ///
    /// # Arguments
    /// * `vcf_iters` - Per-sample VCF iterators to merge.
    /// * `sorted_chromosomes` - The expected chromosome processing order.
    pub fn new(
        vcf_iters: Vec<VarIterator<B>>,
        sorted_chromosomes: Vec<String>,
    ) -> VcfResult<Self> {
        let mut seen_samples = HashSet::new();
        let mut iter_info = Vec::with_capacity(vcf_iters.len());
        for iter in &vcf_iters {
            for sample in iter.samples() {
                if !seen_samples.insert(sample) {
                    return Err(VcfParseError::RuntimeError {
                        message: format!("Sample '{}' appears in multiple VCF files", sample),
                    });
                }
            }
            iter_info.push(VariantIteratorInfo {
                samples: iter.samples().to_vec(),
            });
        }

        Ok(Self {
            vcf_iters,
            iter_info,
            sorted_chromosomes,
            current_chrom_idx: 0,
            chroms_seen: HashSet::new(),
            last_group_start: None,
            done: false,
        })
    }

    pub fn iter_info(&self) -> &[VariantIteratorInfo] {
        &self.iter_info
    }

    fn fail(&mut self, error: VcfParseError) -> Option<VcfResult<OverlappingVariantGroup>> {
        self.done = true;
        Some(Err(error))
    }

    /// Phase A: find the next chromosome span seed:
    /// (chromosome, earliest start, current end).
    ///
    /// The parallel scan peeks all iterators concurrently using rayon,
    /// then reduces the results on the main thread.
    fn compute_next_span_seed(&mut self) -> VcfResult<Option<(String, u32, u32)>> {
        loop {
            if self.current_chrom_idx >= self.sorted_chromosomes.len() {
                self.done = true;
                return Ok(None);
            }

            let chrom = &self.sorted_chromosomes[self.current_chrom_idx];
            let chroms_seen = &self.chroms_seen;
            let last_group_start = self.last_group_start;

            // Parallel peek: each iterator is peeked on a rayon thread.
            let peek_results: Vec<PeekResult> = self
                .vcf_iters
                .par_iter_mut()
                .map(|iter| {
                    let r = match iter.peek_variant() {
                        Ok(Some(r)) => r,
                        Ok(None) => return PeekResult::Exhausted,
                        Err(e) => return PeekResult::Error(e),
                    };

                    if r.chrom.as_str() != chrom.as_str() {
                        if chroms_seen.contains(&r.chrom) {
                            return PeekResult::Error(VcfParseError::RuntimeError {
                                message: format!(
                                    "Out-of-order chromosome '{}' at position {}: chromosome already processed",
                                    r.chrom, r.pos
                                ),
                            });
                        }
                        return PeekResult::DifferentChrom;
                    }

                    if let Some(lbs) = last_group_start {
                        if r.pos < lbs {
                            return PeekResult::Error(VcfParseError::RuntimeError {
                                message: format!(
                                    "Out-of-order position on '{}': {} < last bin start {}",
                                    chrom, r.pos, lbs
                                ),
                            });
                        }
                    }

                    PeekResult::OnCurrentChrom {
                        pos: r.pos,
                        span_end: ref_span_end(r),
                        is_variable: r.alleles.len() > 1,
                    }
                })
                .collect();

            // Reduce: find earliest position and check all_non_variable.
            let mut earliest_pos: Option<u32> = None;
            let mut earliest_end: u32 = 0;
            let mut all_vars_in_same_pos_and_not_variable = true;

            for result in &peek_results {
                match result {
                    PeekResult::Error(_) => {
                        // Extract the error by consuming the vec.
                        // Find the first error and return it.
                        for result in peek_results {
                            if let PeekResult::Error(e) = result {
                                return Err(e);
                            }
                        }
                        unreachable!();
                    }
                    PeekResult::Exhausted | PeekResult::DifferentChrom => continue,
                    &PeekResult::OnCurrentChrom {
                        pos,
                        span_end,
                        is_variable,
                    } => {
                        let peek_pos = pos;
                        let p_end = span_end;

                        match earliest_pos {
                            None => {
                                earliest_pos = Some(peek_pos);
                                earliest_end = p_end;
                                if is_variable {
                                    all_vars_in_same_pos_and_not_variable = false;
                                }
                            }
                            Some(ep) if peek_pos < ep => {
                                earliest_pos = Some(peek_pos);
                                earliest_end = p_end;
                                all_vars_in_same_pos_and_not_variable = !is_variable;
                            }
                            Some(ep) if peek_pos == ep => {
                                if p_end > earliest_end {
                                    earliest_end = p_end;
                                }
                                if is_variable {
                                    all_vars_in_same_pos_and_not_variable = false;
                                }
                            }
                            _ => {
                                all_vars_in_same_pos_and_not_variable = false;
                            }
                        }
                    }
                }
            }

            match earliest_pos {
                Some(pos) => {
                    if all_vars_in_same_pos_and_not_variable {
                        // Skip: consume non-variable variants in parallel.
                        self.vcf_iters.par_iter_mut().for_each(|iter| {
                            if let Ok(Some(r)) = iter.peek_variant() {
                                if r.chrom.as_str() == chrom.as_str() {
                                    iter.next();
                                }
                            }
                        });
                        self.last_group_start = Some(pos);
                        continue;
                    }
                    return Ok(Some((chrom.clone(), pos, earliest_end)));
                }
                None => {
                    self.chroms_seen.insert(chrom.clone());
                    self.current_chrom_idx += 1;
                    self.last_group_start = None;
                }
            }
        }
    }

    /// Phase B: consume all records overlapping the current bin span.
    fn group_overlapping_vars(
        &mut self,
        current_chrom: &str,
        mut group_end: u32,
    ) -> VcfResult<(Vec<Variant>, Vec<usize>, u32)> {
        let mut records: Vec<Variant> = Vec::new();
        let mut source_indices: Vec<usize> = Vec::new();

        loop {
            let mut added_any = false;

            for i in 0..self.vcf_iters.len() {
                loop {
                    let peek_pos = {
                        let Some(r) = self.vcf_iters[i].peek_variant()? else {
                            break;
                        };

                        if r.chrom.as_str() != current_chrom {
                            if self.chroms_seen.contains(&r.chrom) {
                                return Err(VcfParseError::RuntimeError {
                                    message: format!(
                                        "Out-of-order chromosome '{}' at position {}: chromosome already processed",
                                        r.chrom, r.pos
                                    ),
                                });
                            }
                            break;
                        }

                        if let Some(lbs) = self.last_group_start {
                            if r.pos < lbs {
                                return Err(VcfParseError::RuntimeError {
                                    message: format!(
                                        "Out-of-order position on '{}': {} < last bin start {}",
                                        current_chrom, r.pos, lbs
                                    ),
                                });
                            }
                        }

                        r.pos
                    };

                    if peek_pos > group_end {
                        break;
                    }

                    let record =
                        self.vcf_iters[i]
                            .next()
                            .ok_or_else(|| VcfParseError::RuntimeError {
                                message:
                                    "Iterator ended unexpectedly while consuming overlapping record"
                                        .to_string(),
                            })??;

                    let new_end = ref_span_end(&record);
                    if new_end > group_end {
                        group_end = new_end;
                    }

                    records.push(record);
                    source_indices.push(i);
                    added_any = true;
                }
            }

            if !added_any {
                break;
            }
        }

        Ok((records, source_indices, group_end))
    }

    fn chromosome_has_remaining_variants(&mut self, current_chrom: &str) -> VcfResult<bool> {
        for i in 0..self.vcf_iters.len() {
            match self.vcf_iters[i].peek_variant()? {
                Some(r) if r.chrom == current_chrom => return Ok(true),
                _ => {}
            }
        }
        Ok(false)
    }

    fn get_next_variant_group(&mut self) -> Option<VcfResult<OverlappingVariantGroup>> {
        if self.done {
            return None;
        }

        let (current_chrom, group_start, group_end) = match self.compute_next_span_seed() {
            Ok(Some(v)) => v,
            Ok(None) => return None,
            Err(e) => return self.fail(e),
        };

        let (overlapping_vars, source_indices, group_end) =
            match self.group_overlapping_vars(&current_chrom, group_end) {
                Ok(v) => v,
                Err(e) => return self.fail(e),
            };

        self.last_group_start = Some(group_start);

        let any_remaining = match self.chromosome_has_remaining_variants(&current_chrom) {
            Ok(v) => v,
            Err(e) => return self.fail(e),
        };

        if !any_remaining {
            self.chroms_seen.insert(current_chrom.clone());
            self.current_chrom_idx += 1;
            self.last_group_start = None;
        }

        Some(Ok(OverlappingVariantGroup {
            chrom: current_chrom,
            start: group_start,
            end: group_end,
            variants: overlapping_vars,
            source_var_iter_idxs: source_indices,
        }))
    }
}

impl<B: BufRead + Send> Iterator for VariantGroupIterator<B> {
    type Item = VcfResult<OverlappingVariantGroup>;

    fn next(&mut self) -> Option<Self::Item> {
        self.get_next_variant_group()
    }
}
