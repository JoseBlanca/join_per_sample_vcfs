//! Per-sample reading (appendix §A) — module home for `SamplePspReader`
//! and `SamplePspChunk`.
//!
//! *(today: `var_calling::column_span_reader`, backed by `psp::PspReader`)*
//!
//! All per-sample responsibility lives here:
//!
//! - `SamplePspReader<R: Read + Seek>` — one per sample, created from
//!   `(psp::PspReader, region/bed segments)`. **No dust** (the producer is
//!   the sole dust consumer, §2.3). `!Send`. Hands out one `SamplePspChunk`
//!   per psp segment (the natural read unit).
//! - `SamplePspChunk` — one sample's columns for **one psp segment**.
//!   `new()` decodes the *light* fold columns eagerly (`positions`,
//!   `nonref_obs`, `ref_spans`); the typed getters (`take_seq` /
//!   `take_chain_ids` / `take_fixed`) decode their *heavy* column(s) lazily,
//!   for the `keep` rows only, moving the result out.
//!
//! Phase 1 builds this with a *simple* decode (decode the needed columns
//! correctly); the column-selective skip-the-rest optimisation is Phase 5.

// Phase 1.
