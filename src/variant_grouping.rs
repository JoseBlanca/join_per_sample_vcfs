use std::collections::HashSet;
use std::io::BufRead;

use rayon::prelude::*;

use crate::errors::VcfParseError;
use crate::gvcf_parser::{VarIterator, Variant, VcfResult};

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

/// Seed for a new variant group: chromosome, start position, and initial end.
struct ChromSpan {
    chrom: String,
    start: u32,
    end: u32,
}

/// Collected overlapping variants for a group being built.
struct CollectedOverlaps {
    variants: Vec<Variant>,
    source_indices: Vec<usize>,
    group_end: u32,
}

/// Metadata about a `VariantIterator`
#[derive(Debug, Clone)]
pub struct VarIteratorInfo {
    pub samples: Vec<String>,
}

/// A bin of overlapping variants collected from multiple per-sample VCFs.
///
/// The bin covers a contiguous genomic region where all contained variants
/// have reference-allele spans that overlap with the bin's span.
#[derive(Debug)]
pub struct OverlappingVarGroup {
    pub chrom: String,
    /// Start position of the bin (1-based, inclusive).
    pub start: u32,
    /// End position of the bin (1-based, inclusive).
    pub end: u32,
    /// All variant records in this bin, in consumption order.
    pub variants: Vec<Variant>,
    /// Index into `VarGroupIterator::iter_info` for each variant,
    /// parallel to `variants`.
    pub source_var_iter_idxs: Vec<usize>,
}

impl OverlappingVarGroup {
    pub fn span(&self) -> (&str, u32, u32) {
        (&self.chrom, self.start, self.end)
    }
}

/// Validates that a record respects chromosome and position ordering constraints.
///
/// A chromosome is considered "already processed" if it appears in
/// `sorted_chromosomes` at an index before `current_chrom_idx`.
fn validate_var_order(
    var: &Variant,
    current_chrom: &str,
    sorted_chromosomes: &[String],
    current_chrom_idx: usize,
    last_group_start: Option<u32>,
) -> Result<(), VcfParseError> {
    if var.chrom != current_chrom {
        let already_processed = sorted_chromosomes[..current_chrom_idx]
            .iter()
            .any(|c| c == &var.chrom);
        if already_processed {
            return Err(VcfParseError::RuntimeError {
                message: format!(
                    "Out-of-order chromosome '{}' at position {}: chromosome already processed",
                    var.chrom, var.pos
                ),
            });
        }
        return Ok(());
    }

    if let Some(last_group_start) = last_group_start {
        if var.pos < last_group_start {
            return Err(VcfParseError::RuntimeError {
                message: format!(
                    "Out-of-order position on '{}': {} < last var group start {}",
                    current_chrom, var.pos, last_group_start
                ),
            });
        }
    }

    Ok(())
}

/// Iterator that merges multiple per-sample VCF streams, grouping
/// overlapping variants into bins based on reference-allele span overlap.
///
/// Variants from all input VCFs are consumed in genomic order. Two variants
/// overlap when they are on the same chromosome and one's start position
/// falls within the other's reference-allele span.
pub struct VarGroupIterator<B: BufRead> {
    vcf_iters: Vec<VarIterator<B>>,
    iter_info: Vec<VarIteratorInfo>,
    sorted_chromosomes: Vec<String>,
    current_chrom_idx: usize,
    last_group_start: Option<u32>,
    done: bool,
}

impl<B: BufRead + Send> VarGroupIterator<B> {
    /// Creates a new `VarGroupIterator`.
    ///
    /// # Arguments
    /// * `vcf_iters` - Per-sample VCF iterators to merge.
    /// * `sorted_chromosomes` - The expected chromosome processing order.
    pub fn new(vcf_iters: Vec<VarIterator<B>>, sorted_chromosomes: Vec<String>) -> VcfResult<Self> {
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
            iter_info.push(VarIteratorInfo {
                samples: iter.samples().to_vec(),
            });
        }

        Ok(Self {
            vcf_iters,
            iter_info,
            sorted_chromosomes,
            current_chrom_idx: 0,
            last_group_start: None,
            done: false,
        })
    }

    pub fn iter_info(&self) -> &[VarIteratorInfo] {
        &self.iter_info
    }

    fn fail(&mut self, error: VcfParseError) -> Option<VcfResult<OverlappingVarGroup>> {
        self.done = true;
        Some(Err(error))
    }

    /// Advances to the next chromosome, resetting position tracking.
    fn advance_chromosome(&mut self) {
        self.current_chrom_idx += 1;
        self.last_group_start = None;
    }

    /// Phase A: find the next chromosome span seed.
    ///
    /// The parallel scan peeks all iterators concurrently using rayon,
    /// then reduces the results on the main thread.
    ///
    /// This method only peeks iterators — it does not consume the seed
    /// records. The caller must pass the seed to `group_overlapping_vars`,
    /// which consumes them as part of its overlap loop.
    fn compute_next_span_seed_and_skip_non_variable(&mut self) -> VcfResult<Option<ChromSpan>> {
        loop {
            if self.current_chrom_idx >= self.sorted_chromosomes.len() {
                self.done = true;
                return Ok(None);
            }

            let chrom = &self.sorted_chromosomes[self.current_chrom_idx];
            let sorted_chromosomes = &self.sorted_chromosomes;
            let current_chrom_idx = self.current_chrom_idx;
            let last_group_start = self.last_group_start;

            // Parallel peek: each iterator is peeked on a rayon thread.
            let peek_results: Vec<PeekResult> = self
                .vcf_iters
                .par_iter_mut()
                .map(|var_iter| {
                    let var = match var_iter.peek_variant() {
                        Ok(Some(var)) => var,
                        Ok(None) => return PeekResult::Exhausted,
                        Err(e) => return PeekResult::Error(e),
                    };

                    if let Err(e) = validate_var_order(
                        var,
                        chrom,
                        sorted_chromosomes,
                        current_chrom_idx,
                        last_group_start,
                    ) {
                        return PeekResult::Error(e);
                    }

                    if var.chrom != *chrom {
                        return PeekResult::DifferentChrom;
                    }

                    match var.get_var_end() {
                        Ok(span_end) => PeekResult::OnCurrentChrom {
                            pos: var.pos,
                            span_end,
                            is_variable: var.alleles.len() > 1,
                        },
                        Err(e) => PeekResult::Error(e),
                    }
                })
                .collect();

            // Check for errors first, then reduce to find earliest position.
            if let Some(error) = peek_results
                .iter()
                .position(|r| matches!(r, PeekResult::Error(_)))
            {
                // Extract the error by consuming the vec at the known index.
                let PeekResult::Error(e) = peek_results.into_iter().nth(error).unwrap() else {
                    unreachable!()
                };
                return Err(e);
            }

            let mut earliest_pos: Option<u32> = None;
            let mut earliest_end: u32 = 0;
            // we are skipping groups that are completely not variable in order to avoid further processing them
            // this is an optimization for performance. In most cases most genomic position won't be variable
            let mut skip_group = true;

            for result in &peek_results {
                let &PeekResult::OnCurrentChrom {
                    pos,
                    span_end,
                    is_variable,
                } = result
                else {
                    continue;
                };

                match earliest_pos {
                    None => {
                        earliest_pos = Some(pos);
                        earliest_end = span_end;
                        if is_variable {
                            skip_group = false;
                        }
                    }
                    Some(ep) if pos < ep => {
                        earliest_pos = Some(pos);
                        earliest_end = span_end;
                        skip_group = !is_variable;
                    }
                    Some(ep) if pos == ep => {
                        if span_end > earliest_end {
                            earliest_end = span_end;
                        }
                        if is_variable {
                            skip_group = false;
                        }
                    }
                    _ => {
                        skip_group = false;
                    }
                }
            }

            match earliest_pos {
                Some(pos) => {
                    if skip_group {
                        // All peeked variants are at the same position and non-variable;
                        // consume them in parallel.
                        self.vcf_iters.par_iter_mut().for_each(|iter| {
                            if let Ok(Some(r)) = iter.peek_variant() {
                                if r.chrom == *chrom {
                                    debug_assert_eq!(r.pos, pos);
                                    iter.next();
                                }
                            }
                        });
                        self.last_group_start = Some(pos);
                        continue;
                    }
                    return Ok(Some(ChromSpan {
                        chrom: chrom.clone(),
                        start: pos,
                        end: earliest_end,
                    }));
                }
                None => {
                    self.advance_chromosome();
                }
            }
        }
    }

    /// Phase B: consume all records overlapping the current bin span.
    fn group_overlapping_vars(
        &mut self,
        current_chrom: &str,
        mut group_end: u32,
    ) -> VcfResult<CollectedOverlaps> {
        let mut variants: Vec<Variant> = Vec::new();
        let mut source_indices: Vec<usize> = Vec::new();

        loop {
            let mut added_any = false;

            for i in 0..self.vcf_iters.len() {
                loop {
                    let peek_pos = {
                        let Some(r) = self.vcf_iters[i].peek_variant()? else {
                            break;
                        };

                        validate_var_order(
                            r,
                            current_chrom,
                            &self.sorted_chromosomes,
                            self.current_chrom_idx,
                            self.last_group_start,
                        )?;

                        if r.chrom != current_chrom {
                            break;
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

                    let new_end = record.get_var_end()?;
                    if new_end > group_end {
                        group_end = new_end;
                    }

                    variants.push(record);
                    source_indices.push(i);
                    added_any = true;
                }
            }

            if !added_any {
                break;
            }
        }

        Ok(CollectedOverlaps {
            variants,
            source_indices,
            group_end,
        })
    }

    fn chromosome_has_remaining_variants(&mut self, current_chrom: &str) -> VcfResult<bool> {
        for iter in &mut self.vcf_iters {
            if let Some(r) = iter.peek_variant()? {
                if r.chrom == current_chrom {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }

    fn get_next_var_group(&mut self) -> Option<VcfResult<OverlappingVarGroup>> {
        if self.done {
            return None;
        }

        let group_span_seed = match self.compute_next_span_seed_and_skip_non_variable() {
            Ok(Some(g_span_seed)) => g_span_seed,
            Ok(None) => return None,
            Err(e) => return self.fail(e),
        };

        let overlaps =
            match self.group_overlapping_vars(&group_span_seed.chrom, group_span_seed.end) {
                Ok(v) => v,
                Err(e) => return self.fail(e),
            };

        self.last_group_start = Some(group_span_seed.start);

        let any_remaining = match self.chromosome_has_remaining_variants(&group_span_seed.chrom) {
            Ok(v) => v,
            Err(e) => return self.fail(e),
        };

        if !any_remaining {
            self.advance_chromosome();
        }

        Some(Ok(OverlappingVarGroup {
            chrom: group_span_seed.chrom,
            start: group_span_seed.start,
            end: overlaps.group_end,
            variants: overlaps.variants,
            source_var_iter_idxs: overlaps.source_indices,
        }))
    }
}

impl<B: BufRead + Send> Iterator for VarGroupIterator<B> {
    type Item = VcfResult<OverlappingVarGroup>;

    fn next(&mut self) -> Option<Self::Item> {
        self.get_next_var_group()
    }
}

impl<B: BufRead + Send> std::iter::FusedIterator for VarGroupIterator<B> {}
