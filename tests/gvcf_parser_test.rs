use merge_per_sample_vcfs::gvcf_parser::{VarIterator, Variant};
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
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
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
    let mut var_iterator = VarIterator::from_reader(reader).expect("Failed to create parser");

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
    let mut var_iterator = VarIterator::from_reader(reader).expect("Failed to create parser");

    // Consume all variants
    while var_iterator.next().is_some() {}

    // Peek on exhausted iterator should return None
    let peeked = var_iterator.peek_variant().unwrap();
    assert!(peeked.is_none());
}

#[test]
fn test_gzip_reader() {
    let file = File::open("tests/data/sample.g.vcf.gz").expect("Problem opening test file");
    let records = VarIterator::from_gzip_reader(file).expect("Failed to create parser");

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
    let records = VarIterator::from_gzip_path(path).expect("Problem opening test file");

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
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
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
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");

    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    // First variant: 20\t17330\t.\tT\tA,<NON_REF>\t.\tq10
    assert_eq!(variants[1].alleles[0].len(), 1);
    assert_eq!(variants[1].alleles[0], "T");
    assert!(variants[1].qual.is_nan());

    // Second variant: 20\t17331\t.\tA\tG,T,<NON_REF>\t67\tPASS
    assert_eq!(variants[2].alleles[0].len(), 1);
    assert_eq!(variants[2].alleles[0], "A");
    assert_eq!(variants[2].qual, 67.0);

    // Third variant: 20\t17333\t.\tGTC\tG,GTCT,<NON_REF>\t50
    assert_eq!(variants[4].alleles[0].len(), 3);
    assert_eq!(variants[4].alleles[0], "GTC");
    assert_eq!(variants[4].qual, 50.0);

    // Fourth variant: 20\t17334\t.\tGTC\tG,GTCT,<NON_REF>\t50
    assert_eq!(variants[5].alleles[0].len(), 3);
    assert_eq!(variants[5].alleles[0], "GTC");
    assert_eq!(variants[5].qual, 50.0);
}

#[test]
fn test_non_ref_not_included_in_alleles() {
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t100\t.\tG\t<NON_REF>\t30\tPASS\t.\tGT\t0/0
chr1\t200\t.\tT\tA,<NON_REF>\t40\tPASS\t.\tGT\t0/1";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 2);

    // ALT=<NON_REF> only: alleles should contain just the ref allele
    assert_eq!(variants[0].alleles, vec!["G"]);

    // ALT=A,<NON_REF>: alleles should contain ref + A, but not <NON_REF>
    assert_eq!(variants[1].alleles, vec!["T", "A"]);
}

#[test]
fn test_dot_alt_not_included_in_alleles() {
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t100\t.\tG\t.\t30\tPASS\t.\tGT\t0/0
chr1\t200\t.\tT\tA,.\t40\tPASS\t.\tGT\t0/1";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 2);

    // ALT=. only: alleles should contain just the ref allele
    assert_eq!(variants[0].alleles, vec!["G"]);

    // ALT=A,.: alleles should contain ref + A, but not "."
    assert_eq!(variants[1].alleles, vec!["T", "A"]);
}

#[test]
fn test_sample_parsing() {
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");

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
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");

    // Samples should be available immediately after construction
    let samples = parser.samples();
    assert_eq!(samples.len(), 3);
    assert_eq!(samples[0], "NA00001");
    assert_eq!(samples[1], "NA00002");
    assert_eq!(samples[2], "NA00003");

    // Iteration should still work normally after header parsing
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();
    assert!(!variants.is_empty());
}

#[test]
fn test_sample_parsing_single_sample() {
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t10\t.\tA\tC\t30\tPASS\t.\tGT\t0/1";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");

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
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 2);

    // Variant 1: chr1:100, unphased
    let v1 = &variants[0];
    assert_eq!(v1.pos, 100);
    assert_eq!(v1.genotypes.len(), 6); // 3 samples * ploidy 2
    // Sample 1: 0/0
    assert_eq!(v1.genotypes[0], 0);
    assert_eq!(v1.genotypes[1], 0);
    assert!(!v1.phases[0]);
    // Sample 2: 0/1
    assert_eq!(v1.genotypes[2], 0);
    assert_eq!(v1.genotypes[3], 1);
    assert!(!v1.phases[1]);
    // Sample 3: 1/1
    assert_eq!(v1.genotypes[4], 1);
    assert_eq!(v1.genotypes[5], 1);
    assert!(!v1.phases[2]);

    // Variant 2: chr1:200, phased
    let v2 = &variants[1];
    assert_eq!(v2.pos, 200);
    // Sample 1: 0|0
    assert_eq!(v2.genotypes[0], 0);
    assert_eq!(v2.genotypes[1], 0);
    assert!(v2.phases[0]);
    // Sample 2: 0|1
    assert_eq!(v2.genotypes[2], 0);
    assert_eq!(v2.genotypes[3], 1);
    assert!(v2.phases[1]);
    // Sample 3: 1|1
    assert_eq!(v2.genotypes[4], 1);
    assert_eq!(v2.genotypes[5], 1);
    assert!(v2.phases[2]);
}

#[test]
fn test_genotype_parsing_missing() {
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3
chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT\t./.\t0/1\t.|.";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 1);
    let v = &variants[0];
    // Sample 1: ./.  → missing, unphased
    assert_eq!(v.genotypes[0], -1);
    assert_eq!(v.genotypes[1], -1);
    assert!(!v.phases[0]);
    // Sample 2: 0/1 → het, unphased
    assert_eq!(v.genotypes[2], 0);
    assert_eq!(v.genotypes[3], 1);
    assert!(!v.phases[1]);
    // Sample 3: .|. → missing, phased
    assert_eq!(v.genotypes[4], -1);
    assert_eq!(v.genotypes[5], -1);
    assert!(v.phases[2]);
}

#[test]
fn test_genotype_parsing_multiallelic() {
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2
chr1\t100\t.\tA\tC,G,T\t30\tPASS\t.\tGT\t0/1\t2/3";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 1);
    let v = &variants[0];
    // Sample 1: 0/1
    assert_eq!(v.genotypes[0], 0);
    assert_eq!(v.genotypes[1], 1);
    assert!(!v.phases[0]);
    // Sample 2: 2/3
    assert_eq!(v.genotypes[2], 2);
    assert_eq!(v.genotypes[3], 3);
    assert!(!v.phases[1]);
}

#[test]
fn test_genotype_with_format_fields() {
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2
chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT:GQ:DP\t0/1:30:10\t1|1:40:15";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 1);
    let v = &variants[0];
    // Sample 1: 0/1, unphased
    assert_eq!(v.genotypes[0], 0);
    assert_eq!(v.genotypes[1], 1);
    assert!(!v.phases[0]);
    // Sample 2: 1|1, phased
    assert_eq!(v.genotypes[2], 1);
    assert_eq!(v.genotypes[3], 1);
    assert!(v.phases[1]);

    // Verify FORMAT field access
    let gq_idx = v.gt_field_index("GQ").expect("GQ field should exist");
    let gq_values = v.get_gt_field_by_index(gq_idx);
    assert_eq!(gq_values[0], "30");
    assert_eq!(gq_values[1], "40");

    let dp_idx = v.gt_field_index("DP").expect("DP field should exist");
    let dp_values = v.get_gt_field_by_index(dp_idx);
    assert_eq!(dp_values[0], "10");
    // Last sample's DP must be exactly "15", not "15\n" (M1 regression)
    assert_eq!(dp_values[1], "15");
}

#[test]
fn test_genotype_performance_with_sample_gvcf() {
    // Use the existing SAMPLE_GVCF data which has genotypes
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
    let variants: Vec<_> = parser.filter_map(Result::ok).collect();

    assert_eq!(variants.len(), 6);

    // Check first variant: 20\t17330\t.\tT\tA,<NON_REF>\t.\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t.|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
    let v1 = &variants[1];
    assert_eq!(v1.n_samples, 3);

    // Check second variant has consistent ploidy
    let v2 = &variants[2];
    assert_eq!(v2.n_samples, 3);
}

#[test]
fn test_parse_rejects_allele_index_out_of_range() {
    // Allele index 128 exceeds i8::MAX (127) and must be rejected, not silently wrapped
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT\t128/0";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
    let results: Vec<_> = parser.collect();

    assert_eq!(results.len(), 1);
    assert!(
        results[0].is_err(),
        "Allele index 128 should produce an error, not wrap to -128"
    );
}

#[test]
fn test_parse_rejects_invalid_characters() {
    // "0/a" contains a non-digit, non-separator, non-dot character
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT\t0/a";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
    let results: Vec<_> = parser.collect();
    assert_eq!(results.len(), 1);
    assert!(results[0].is_err(), "GT '0/a' should be rejected");
}

#[test]
fn test_parse_rejects_consecutive_separators() {
    // "0//1" has consecutive separators with no allele between them
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT\t0//1";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
    let results: Vec<_> = parser.collect();
    assert_eq!(results.len(), 1);
    assert!(results[0].is_err(), "GT '0//1' should be rejected");
}

#[test]
fn test_parse_rejects_trailing_junk() {
    // "0/1x" has a trailing non-digit character
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT\t0/1x";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
    let results: Vec<_> = parser.collect();
    assert_eq!(results.len(), 1);
    assert!(results[0].is_err(), "GT '0/1x' should be rejected");
}

#[test]
fn test_parse_rejects_leading_separator() {
    // "/0/1" starts with a separator
    let vcf_data = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT\t/0/1";

    let reader = BufReader::new(vcf_data.as_bytes());
    let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
    let results: Vec<_> = parser.collect();
    assert_eq!(results.len(), 1);
    assert!(results[0].is_err(), "GT '/0/1' should be rejected");
}

#[test]
fn test_get_span_returns_error_on_position_overflow() {
    // Position near u32::MAX with an allele long enough to overflow
    let variant = Variant::new(
        "chr1".to_string(),
        u32::MAX - 5,
        vec!["AAAAAAAAAAAA".to_string()], // 12 bases, pos + 11 overflows
        vec![0, 0],
        vec![false],
        1,
    );
    let result = variant.get_span();
    assert!(
        result.is_err(),
        "get_span should return error on position overflow, not wrap"
    );
}
