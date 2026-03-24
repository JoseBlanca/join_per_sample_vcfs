use std::env;
use std::process;
use std::thread;

use join_per_sample_vcfs::decompression_pool::DecompressionPool;
use join_per_sample_vcfs::gvcf_parser::VariantIterator;
use join_per_sample_vcfs::pipeline::merge_alleles_and_genotypes;

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

    // Create decompression pool (defaults to number of CPUs)
    let decompression_threads = n_threads.unwrap_or_else(|| {
        thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4)
    });
    let pool = DecompressionPool::new(decompression_threads);

    // Parse each gVCF file using the shared pool
    let mut vcf_iters = Vec::with_capacity(vcf_paths.len());
    for path in &vcf_paths {
        let iter = VariantIterator::from_gzip_path_pooled(path, &pool)?;
        vcf_iters.push(iter);
    }

    let stdout = Box::new(std::io::stdout().lock());
    merge_alleles_and_genotypes(vcf_iters, sorted_chromosomes, stdout)?;

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}
