use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::Path;
use std::time::Duration;

use join_per_sample_vcfs::gvcf_parser::GVcfRecordIterator;

fn bench_parse_gz_end_to_end(c: &mut Criterion) {
    let path = Path::new("/home/jose/analyses/g2psol/source_data/TS.vcf.gz");
    assert!(path.exists(), "Benchmark file not found: {path:?}");

    c.bench_function("gvcf parse TS.vcf.gz (end-to-end)", |b| {
        b.iter(|| {
            let records =
                GVcfRecordIterator::from_gzip_path(path).expect("Problem opening test file");

            let mut n_variants: u32 = 0;
            for record in records {
                match record {
                    Ok(_variant) => n_variants += 1,
                    Err(error) => panic!("Unexpected error: {error}"),
                }
            }
            black_box(n_variants)
        })
    });
}

fn config() -> Criterion {
    Criterion::default()
        .sample_size(10)
        .measurement_time(Duration::from_secs(1))
}

criterion_group! {
    name = benches;
    config = config();
    targets = bench_parse_gz_end_to_end
}

criterion_main!(benches);
