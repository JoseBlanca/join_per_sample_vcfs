use std::env;
use std::process;

use join_per_sample_vcfs::genotype_merging::analyze_groups;
use join_per_sample_vcfs::gvcf_parser::VariantIterator;
use join_per_sample_vcfs::variant_grouping::VariantGroupIterator;

fn print_usage(program: &str) {
    eprintln!(
        "Usage: {} --chroms <chrom1,chrom2,...> <vcf1.gz> <vcf2.gz> ...",
        program
    );
    eprintln!();
    eprintln!("Merges per-sample gVCF files into a single multi-sample VCF.");
    eprintln!();
    eprintln!("Options:");
    eprintln!("  --chroms <list>     Comma-separated chromosome names in processing order");
    eprintln!("  --threads <N>       Number of threads for parallel merging (default: all cores)");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  <vcf.gz>         One or more gzipped gVCF files (one per sample)");
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    let program = &args[0];

    if args.len() < 2 {
        print_usage(program);
        process::exit(1);
    }

    // Parse arguments
    let mut chroms: Option<Vec<String>> = None;
    let mut vcf_paths: Vec<String> = Vec::new();
    let mut n_threads: Option<usize> = None;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--chroms" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: --chroms requires a value");
                    process::exit(1);
                }
                chroms = Some(args[i].split(',').map(|s| s.to_string()).collect());
            }
            "--threads" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: --threads requires a value");
                    process::exit(1);
                }
                let t: usize = args[i].parse().unwrap_or_else(|_| {
                    eprintln!("Error: --threads must be a positive integer");
                    process::exit(1);
                });
                if t == 0 {
                    eprintln!("Error: --threads must be at least 1");
                    process::exit(1);
                }
                n_threads = Some(t);
            }
            "--help" | "-h" => {
                print_usage(program);
                process::exit(0);
            }
            arg => {
                vcf_paths.push(arg.to_string());
            }
        }
        i += 1;
    }

    if vcf_paths.is_empty() {
        eprintln!("Error: at least one VCF file is required");
        process::exit(1);
    }

    let sorted_chromosomes = chroms.ok_or("Error: --chroms is required")?;

    // Configure rayon thread pool
    if let Some(t) = n_threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .map_err(|e| format!("Failed to set thread count: {}", e))?;
    }

    // Parse each gVCF file
    let mut vcf_iters = Vec::with_capacity(vcf_paths.len());
    for path in &vcf_paths {
        let iter = VariantIterator::from_gzip_path(path)?;
        vcf_iters.push(iter);
    }

    // Collect all sample names for the header
    let all_samples: Vec<String> = vcf_iters
        .iter()
        .flat_map(|iter| iter.samples().iter().cloned())
        .collect();

    // Create the variant group iterator
    let grouper = VariantGroupIterator::new(vcf_iters, sorted_chromosomes)?;
    let iter_info = grouper.iter_info().to_vec();

    // Output VCF header
    println!("##fileformat=VCFv4.2");
    print!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for sample in &all_samples {
        print!("\t{}", sample);
    }
    println!();

    // Stream: group -> merge -> write (one variant at a time, no bulk collection)
    let ploidy = 2; // diploid assumption
    for result in analyze_groups(grouper, &iter_info) {
        let variant = result?;

        let ref_allele = &variant.alleles[0];
        let alt_alleles = if variant.alleles.len() > 1 {
            variant.alleles[1..].join(",")
        } else {
            ".".to_string()
        };

        print!(
            "{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT",
            variant.chrom, variant.pos, ref_allele, alt_alleles
        );

        for sample_idx in 0..variant.n_samples {
            let gt_start = sample_idx * ploidy;
            let gt = &variant.genotypes[gt_start..gt_start + ploidy];
            let gt_str: Vec<String> = gt
                .iter()
                .map(|&a| {
                    if a < 0 {
                        ".".to_string()
                    } else {
                        a.to_string()
                    }
                })
                .collect();
            let phase_sep = if variant.phase[sample_idx] { "|" } else { "/" };
            print!("\t{}", gt_str.join(phase_sep));
        }
        println!();
    }

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}
