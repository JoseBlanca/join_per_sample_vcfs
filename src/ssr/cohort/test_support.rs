//! Shared in-memory fixture builders for the cohort tests — catalogs and `.ssr.psp`
//! bytes/readers, so the merge and driver tests don't each reinvent them.

use std::io::Cursor;

use crate::psp::PspReader;
use crate::psp::header::{ChromosomeEntry, ParameterValue, WriterHeader};
use crate::psp::registry_ssr::{SsrKind, SsrLocusRecord};
use crate::psp::test_fixtures::writer_header;
use crate::psp::writer::PspWriter;
use crate::ssr::catalog::CatalogParams;
use crate::ssr::catalog::io::{CatalogHeader, CatalogReader, CatalogWriter};
use crate::ssr::types::{Locus, Motif};

/// The reference md5 the fixtures bind catalog ↔ `.ssr.psp` on.
pub(crate) const REF_MD5: &str = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";

/// The catalog-binding header parameter (mirrors Stage-1's writer).
pub(crate) const CATALOG_MD5_PARAM: &str = "catalog_reference_md5";

pub(crate) fn obs(pairs: &[(&[u8], u32)]) -> Vec<(Box<[u8]>, u32)> {
    pairs
        .iter()
        .map(|(s, c)| (s.to_vec().into_boxed_slice(), *c))
        .collect()
}

/// A catalog locus: CA(3) tract at 0-based `[start, start+6)`, 18-byte ref frame
/// (6 G + CACACA + 6 T) anchored at `start - 6`.
pub(crate) fn loc(chrom: &str, start: u32) -> Locus {
    Locus::new(
        chrom.into(),
        start,
        start + 6,
        Motif::new(b"CA").unwrap(),
        1.0,
        (*b"GGGGGGCACACATTTTTT").into(),
        start - 6,
    )
    .unwrap()
}

pub(crate) fn catalog_bytes(ref_md5: &str, loci: &[Locus]) -> Vec<u8> {
    let header = CatalogHeader {
        tool_version: "0".into(),
        reference: "ref.fa".into(),
        reference_md5: ref_md5.into(),
        trf_mod_version: "1".into(),
        params: CatalogParams::default(),
        date: "2026-06-21".into(),
    };
    let mut w = CatalogWriter::new(Cursor::new(Vec::<u8>::new()), &header).unwrap();
    for l in loci {
        w.write_locus(l).unwrap();
    }
    w.finish().unwrap().into_inner()
}

pub(crate) fn catalog_reader(ref_md5: &str, loci: &[Locus]) -> CatalogReader<Cursor<Vec<u8>>> {
    CatalogReader::new(Cursor::new(catalog_bytes(ref_md5, loci))).unwrap()
}

/// SSR `.psp` header naming `chroms` and binding `cat_md5`.
pub(crate) fn ssr_header(chroms: &[&str], cat_md5: &str) -> WriterHeader {
    let mut header = writer_header(chroms.len());
    header.chromosomes = chroms
        .iter()
        .map(|name| ChromosomeEntry {
            name: (*name).to_string(),
            length: 100_000,
            md5: "0".repeat(32),
        })
        .collect();
    header.writer.parameters.insert(
        CATALOG_MD5_PARAM.to_string(),
        ParameterValue::String(cat_md5.to_string()),
    );
    header
}

/// Like [`ssr_header`] but omits the catalog-binding parameter — for the
/// missing-md5 rejection test.
pub(crate) fn ssr_header_without_md5(chroms: &[&str]) -> WriterHeader {
    let mut header = ssr_header(chroms, REF_MD5);
    header.writer.parameters.remove(CATALOG_MD5_PARAM);
    header
}

/// One container record (1-based coords) — a sample's evidence at 0-based
/// `[start, start+6)` becomes container `[start+1, start+7)`.
pub(crate) fn rec(
    chrom_id: u32,
    start_0based: u32,
    observed: Vec<(Box<[u8]>, u32)>,
) -> SsrLocusRecord {
    SsrLocusRecord {
        chrom_id,
        start: start_0based + 1,
        end: start_0based + 7,
        depth: 30,
        n_filtered: 0,
        mapped_reads: 30,
        n_low_quality: 0,
        n_border_off_end: 0,
        observed,
    }
}

pub(crate) fn ssr_psp(header: WriterHeader, records: &[SsrLocusRecord]) -> Vec<u8> {
    let mut w = PspWriter::<_, SsrKind>::new_ssr(Cursor::new(Vec::<u8>::new()), header).unwrap();
    for r in records {
        w.write_locus(r).unwrap();
    }
    w.finish().unwrap().into_inner()
}

pub(crate) fn reader(bytes: Vec<u8>) -> PspReader<Cursor<Vec<u8>>> {
    PspReader::new(Cursor::new(bytes)).unwrap()
}
