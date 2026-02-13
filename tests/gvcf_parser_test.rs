use join_per_sample_vcfs::errors::VcfParseError;
use join_per_sample_vcfs::gvcf_parser::{GVcfRecord, GVcfRecordIterator};
use std::fs::File;
use std::io::BufReader;

const SAMPLE_GVCF: &str = "##
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
20\t14370\t.\tG\t<NON_REF>\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t1|2:48:1:51,51\t3|4:48:8:51,51\t5/6000:43:5:.,.
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
    assert_eq!(n_variants, 4);
}

#[test]
fn test_buffer() {
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let mut var_iterator = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");
    assert!(matches!(var_iterator.fill_buffer(3), Ok(3)));
    assert!(matches!(var_iterator.fill_buffer(1), Ok(0)));
    assert!(matches!(var_iterator.fill_buffer(4), Ok(1)));
    assert!(matches!(var_iterator.fill_buffer(5), Ok(0)));
    let variant = var_iterator.next().unwrap().unwrap();
    assert_eq!(variant.pos, 17330);
    let buffered_items = var_iterator.peek_items_in_buffer();
    let poss = [17331, 17333, 17334];
    for (expected_pos, variant) in poss.iter().zip(buffered_items) {
        assert_eq!(&variant.pos, expected_pos);
    }
}
#[test]
fn test_buffer2() {
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let mut var_iterator = GVcfRecordIterator::from_reader(reader).expect("Failed to create parser");
    assert!(matches!(var_iterator.fill_buffer(2), Ok(2)));
    let variant = var_iterator.next().unwrap().unwrap();
    assert_eq!(variant.pos, 17330);
    let variant = var_iterator.next().unwrap().unwrap();
    assert_eq!(variant.pos, 17331);

    let buffered_items: Vec<&GVcfRecord> = var_iterator.peek_items_in_buffer().collect();
    assert_eq!(buffered_items.len(), 0);

    let variant = var_iterator.next().unwrap().unwrap();
    assert_eq!(variant.pos, 17333);

    assert!(matches!(var_iterator.fill_buffer(2), Ok(0)));

    let buffered_items = var_iterator.peek_items_in_buffer();
    let poss = [17334];
    for (expected_pos, variant) in poss.iter().zip(buffered_items) {
        assert_eq!(&variant.pos, expected_pos);
    }
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
    assert_eq!(n_variants, 0);
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
    assert_eq!(n_variants, 0);
}

#[ignore = "slow to run"]
#[test]
fn test_performance() {
    let path = "/home/jose/analyses/g2psol/source_data/TS.vcf.gz";
    let records = GVcfRecordIterator::from_gzip_path(path).expect("Problem opening test file");
    println!("{}", path);

    let mut n_variants: u32 = 0;
    let mut n_invariants: u32 = 0;
    for record in records {
        match record {
            Ok(_variant) => {
                n_variants += 1;
            }
            Err(VcfParseError::InvariantgVCFLine) => {
                n_invariants += 1;
            }
            Err(error) => {
                //Fail test
                panic!("Unexpected error: {}", error);
            }
        }
    }
    println!("Num. variant loci: {n_variants}");
    println!("Num. invariant loci: {n_invariants}");
    println!("Num.loci: {}", n_invariants + n_variants);
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
    assert_eq!(variants[0].ref_allele_len, 1);
    assert_eq!(variants[0].alleles[0], "T");
    assert!(variants[0].qual.is_nan());

    // Second variant: 20\t17331\t.\tA\tG,T,<NON_REF>\t67\tPASS
    assert_eq!(variants[1].ref_allele_len, 1);
    assert_eq!(variants[1].alleles[0], "A");
    assert_eq!(variants[1].qual, 67.0);

    // Third variant: 20\t17333\t.\tGTC\tG,GTCT,<NON_REF>\t50
    assert_eq!(variants[2].ref_allele_len, 3);
    assert_eq!(variants[2].alleles[0], "GTC");
    assert_eq!(variants[2].qual, 50.0);

    // Fourth variant: 20\t17334\t.\tGTC\tG,GTCT,<NON_REF>\t50
    assert_eq!(variants[3].ref_allele_len, 3);
    assert_eq!(variants[3].alleles[0], "GTC");
    assert_eq!(variants[3].qual, 50.0);
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
