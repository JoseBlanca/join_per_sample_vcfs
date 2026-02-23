use join_per_sample_vcfs::errors::VcfParseError;
use join_per_sample_vcfs::gvcf_parser::GVcfRecordIterator;
use std::fs::File;
use std::io::BufReader;

const SAMPLE_GVCF: &str = "##
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
20\t14370\t.\tG\t<NON_REF>\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t1|2:48:1:51,51\t3|4:48:8:51,51\t5/60:43:5:.,.
20\t17330\t.\tT\tA,<NON_REF>\t.\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t.|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t17331\t.\tA\tG,T,<NON_REF>\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t17332\t.\tT\t<NON_REF>\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t17333\t.\tGTC\tG,GTCT,<NON_REF>\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
20\t17334\t.\tGTC\tG,GTCT,<NON_REF>\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t.:35:4\t0/2:17:2\t./1:40:3";

#[test]
fn test_gvcf_parsing() {
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let parser = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");
    let mut n_variants: u32 = 0;
    for record in parser {
        match record {
            Ok(_variant) => {
                n_variants += 1;
            }
            Err(error) => {
                //Fail test
                panic!("Unexpected error: {}", error);
            }
        }
    }
    assert_eq!(n_variants, 6);
}

#[test]
fn test_peek_variant() {
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let mut var_iterator =
        GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");

    // Peek should return the first variant without consuming it
    let peeked = var_iterator.peek_variant().unwrap().unwrap();
    assert_eq!(peeked.pos, 14370);

    // Peeking again should return the same variant
    let peeked_again = var_iterator.peek_variant().unwrap().unwrap();
    assert_eq!(peeked_again.pos, 14370);

    // Consuming via next should return the same variant
    let consumed = var_iterator.next().unwrap().unwrap();
    assert_eq!(consumed.pos, 14370);

    // After consuming, peek should show the next variant
    let next_peeked = var_iterator.peek_variant().unwrap().unwrap();
    assert_eq!(next_peeked.pos, 17330);
}

#[test]
fn test_peek_variant_exhausted() {
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let mut var_iterator =
        GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");

    // Consume all variants
    while var_iterator.next().is_some() {}

    // Peek on exhausted iterator should return None
    let peeked = var_iterator.peek_variant().unwrap();
    assert!(peeked.is_none());
}

#[test]
fn test_gzip_reader() {
    let file = File::open("tests/data/sample.g.vcf.gz").expect("Problem opening test file");
    let records = GVcfRecordIterator::from_gzip_reader(file).expect("Failed to create parser");

    let mut n_variants: u32 = 0;
    for record in records {
        match record {
            Ok(_variant) => {
                n_variants += 1;
            }
            Err(error) => {
                //Fail test
                panic!("Unexpected error: {}", error);
            }
        }
    }
    assert_eq!(n_variants, 63);
}
#[test]
fn test_gzip_path() {
    let path = "tests/data/sample.g.vcf.gz";
    let records = GVcfRecordIterator::from_gzip_path(path).expect("Problem opening test file");

    let mut n_variants: u32 = 0;
    for record in records {
        match record {
            Ok(_variant) => {
                n_variants += 1;
            }
            Err(error) => {
                //Fail test
                panic!("Unexpected error: {}", error);
            }
        }
    }
    assert_eq!(n_variants, 63);
}



#[test]
fn test_g_vcf_record() {
    // Test get_span with parsed records from actual VCF data
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t10\t.\tA\tC\t30\tPASS\t.\tGT\t0/1
chr1\t10\t.\tAT\tA\t30\tPASS\t.\tGT\t0/1
chr1\t10\t.\tA\tATT\t30\tPASS\t.\tGT\t0/1";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 3);

    // SNP: should span a single position
    assert!(matches!(variants[0].get_span(), Ok((10, 10))));

    // Deletion: AT -> A, should span 2 positions
    assert!(matches!(variants[1].get_span(), Ok((10, 11))));

    // Insertion: A -> ATT, should span 3 positions
    assert!(matches!(variants[2].get_span(), Ok((10, 12))));
}

#[test]
fn test_ref_allele_len() {
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let parser = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");

    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    // First variant: 20\t17330\t.\tT\tA,<NON_REF>\t.\tq10
    assert_eq!(variants[1].ref_allele_len, 1);
    assert_eq!(variants[1].alleles[0], "T");
    assert!(variants[1].qual.is_nan());

    // Second variant: 20\t17331\t.\tA\tG,T,<NON_REF>\t67\tPASS
    assert_eq!(variants[2].ref_allele_len, 1);
    assert_eq!(variants[2].alleles[0], "A");
    assert_eq!(variants[2].qual, 67.0);

    // Third variant: 20\t17333\t.\tGTC\tG,GTCT,<NON_REF>\t50
    assert_eq!(variants[4].ref_allele_len, 3);
    assert_eq!(variants[4].alleles[0], "GTC");
    assert_eq!(variants[4].qual, 50.0);

    // Fourth variant: 20\t17334\t.\tGTC\tG,GTCT,<NON_REF>\t50
    assert_eq!(variants[5].ref_allele_len, 3);
    assert_eq!(variants[5].alleles[0], "GTC");
    assert_eq!(variants[5].qual, 50.0);
}

#[test]
fn test_sample_parsing() {
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let parser = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");

    // Samples should be available immediately after construction
    let samples = parser.samples();
    assert_eq!(samples.len(), 3);
    assert_eq!(samples[0], "NA00001");
    assert_eq!(samples[1], "NA00002");
    assert_eq!(samples[2], "NA00003");
}

#[test]
fn test_sample_parsing_with_buffer() {
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let mut parser = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");

    // Samples should be available immediately after construction
    let samples = parser.samples();
    assert_eq!(samples.len(), 3);
    assert_eq!(samples[0], "NA00001");
    assert_eq!(samples[1], "NA00002");
    assert_eq!(samples[2], "NA00003");

    // Fill buffer should still work normally
    parser.fill_buffer(2).unwrap();
}

#[test]
fn test_sample_parsing_single_sample() {
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t10\t.\tA\tC\t30\tPASS\t.\tGT\t0/1";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");

    // Samples should be available immediately after construction
    let samples = parser.samples();
    assert_eq!(samples.len(), 1);
    assert_eq!(samples[0], "Sample1");
}

#[test]
fn test_genotype_parsing_basic() {
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3
chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT\t0/0\t0/1\t1/1
chr1\t200\t.\tG\tT\t40\tPASS\t.\tGT\t0|0\t0|1\t1|1";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 2);

    // First variant: unphased genotypes
    let v1 = &variants[0];
    assert_eq!(v1.ploidy, 2);
    assert_eq!(v1.n_samples(), 3);
    assert_eq!(v1.sample_genotypes(0), &[0, 0]);
    assert_eq!(v1.sample_genotypes(1), &[0, 1]);
    assert_eq!(v1.sample_genotypes(2), &[1, 1]);
    assert!(!v1.is_phased(0));
    assert!(!v1.is_phased(1));
    assert!(!v1.is_phased(2));

    // Second variant: phased genotypes
    let v2 = &variants[1];
    assert_eq!(v2.sample_genotypes(0), &[0, 0]);
    assert_eq!(v2.sample_genotypes(1), &[0, 1]);
    assert_eq!(v2.sample_genotypes(2), &[1, 1]);
    assert!(v2.is_phased(0));
    assert!(v2.is_phased(1));
    assert!(v2.is_phased(2));
}

#[test]
fn test_genotype_parsing_missing() {
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3
chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT\t./.\t0/1\t.|.";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 1);
    let v = &variants[0];

    // Sample 0: missing unphased
    assert_eq!(v.sample_genotypes(0), &[-1, -1]);
    assert!(v.is_missing(0));
    assert!(!v.is_phased(0));

    // Sample 1: normal
    assert_eq!(v.sample_genotypes(1), &[0, 1]);
    assert!(!v.is_missing(1));

    // Sample 2: missing phased
    assert_eq!(v.sample_genotypes(2), &[-1, -1]);
    assert!(v.is_missing(2));
    assert!(v.is_phased(2));
}

#[test]
fn test_genotype_parsing_multiallelic() {
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2
chr1\t100\t.\tA\tC,G,T\t30\tPASS\t.\tGT\t0/1\t2/3";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 1);
    let v = &variants[0];

    assert_eq!(v.sample_genotypes(0), &[0, 1]);
    assert_eq!(v.sample_genotypes(1), &[2, 3]);
}

#[test]
fn test_genotype_api_methods() {
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./.";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    let v = &variants[0];

    // is_hom_ref
    assert!(v.is_hom_ref(0));
    assert!(!v.is_hom_ref(1));
    assert!(!v.is_hom_ref(2));
    assert!(!v.is_hom_ref(3));

    // is_homozygous
    assert!(v.is_homozygous(0)); // 0/0
    assert!(!v.is_homozygous(1)); // 0/1
    assert!(v.is_homozygous(2)); // 1/1
    assert!(v.is_homozygous(3)); // ./.

    // is_missing
    assert!(!v.is_missing(0));
    assert!(!v.is_missing(1));
    assert!(!v.is_missing(2));
    assert!(v.is_missing(3));
}

#[test]
fn test_genotype_with_format_fields() {
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2
chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT:GQ:DP\t0/1:30:10\t1|1:40:15";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 1);
    let v = &variants[0];

    assert_eq!(v.sample_genotypes(0), &[0, 1]);
    assert!(!v.is_phased(0));

    assert_eq!(v.sample_genotypes(1), &[1, 1]);
    assert!(v.is_phased(1));
}

#[test]
fn test_genotype_performance_with_sample_gvcf() {
    // Use the existing SAMPLE_GVCF data which has genotypes
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let parser = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 6);

    // Check first variant: 20\t17330\t.\tT\tA,<NON_REF>\t.\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t.|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
    let v1 = &variants[1];
    assert_eq!(v1.n_samples(), 3);
    assert_eq!(v1.ploidy, 2);

    // Sample 0: .|0 (phased, missing first allele)
    assert_eq!(v1.sample_genotypes(0), &[-1, 0]);
    assert!(v1.is_phased(0));
    assert!(v1.is_missing(0));

    // Sample 1: 0|1 (phased heterozygous)
    assert_eq!(v1.sample_genotypes(1), &[0, 1]);
    assert!(v1.is_phased(1));
    assert!(!v1.is_missing(1));

    // Sample 2: 0/0 (unphased homozygous ref)
    assert_eq!(v1.sample_genotypes(2), &[0, 0]);
    assert!(!v1.is_phased(2));
    assert!(v1.is_hom_ref(2));

    // Check second variant has consistent ploidy
    let v2 = &variants[2];
    assert_eq!(v2.ploidy, 2);
    assert_eq!(v2.n_samples(), 3);
}

#[ignore = "too slow"]
#[test]
fn test_performance() {
    let path = "/home/jose/analyses/g2psol/source_data/TS.vcf.gz";
    let records = GVcfRecordIterator::from_gzip_path(path).expect("Problem opening test file");
    println!("{}", path);

    let mut n_variants: u32 = 0;
    for record in records {
        match record {
            Ok(_variant) => {
                n_variants += 1;
            }
            Err(error) => {
                //Fail test
                panic!("Unexpected error: {}", error);
            }
        }
    }
    println!("Num.loci: {}", n_variants);
}