//! Errors raised by the pileup walker. Each variant carries enough
//! context (qname, chromosome, position) to point at the offending
//! input from the error alone, per `ia/specs/design_principles.md`
//! principle 6 ("typed errors at module boundaries").

use thiserror::Error;

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum WalkerError {
    #[error(
        "out-of-order read: qname='{qname}' at (chrom_id={chrom_id}, pos={pos}) regresses from \
         (chrom_id={prev_chrom_id}, pos={prev_pos})"
    )]
    OutOfOrder {
        qname: String,
        prev_chrom_id: u32,
        prev_pos: u32,
        chrom_id: u32,
        pos: u32,
    },

    #[error(
        "read decoded zero reference bases: qname='{qname}' at (chrom_id={chrom_id}, pos={pos})"
    )]
    ZeroRefSpan {
        qname: String,
        chrom_id: u32,
        pos: u32,
    },

    #[error(
        "active-read cap exceeded (cap={cap}) at chrom_id={chrom_id} pos={pos}; \
         consider raising --max-active-reads or pre-filtering this region"
    )]
    ActiveReadsExhausted { cap: u32, chrom_id: u32, pos: u32 },

    #[error(
        "phase-chain id space exhausted at chrom_id={chrom_id} pos={pos}: \
         this .psp file has reached 2^64 unique read identifiers, the per-file \
         limit imposed by the u64 chain id encoding"
    )]
    ChainIdSpaceExhausted { chrom_id: u32, pos: u32 },

    #[error(
        "pending-mates map exceeded its defensive cap (cap={cap}) at \
         chrom_id={chrom_id} pos={pos}; this usually means the input has \
         a pathologically high orphan-first-mate rate (every paired read \
         flagged FirstOfPair with no SecondOfPair ever arriving). Verify \
         the BAM's pairing flags or pre-filter the affected region."
    )]
    PendingMatesExhausted { cap: usize, chrom_id: u32, pos: u32 },

    #[error(
        "open record reference span exceeded MAX_RECORD_SPAN: anchor (chrom_id={chrom_id}, \
         pos={pos}) reached span {span} (cap={cap}); upstream filter should have rejected the \
         underlying read"
    )]
    RecordTooWide {
        chrom_id: u32,
        pos: u32,
        span: u32,
        cap: u32,
    },

    #[error("FASTA fetch failed at chrom_id={chrom_id} for [{start}, {start_plus_len}): {source}")]
    Fasta {
        chrom_id: u32,
        start: u32,
        start_plus_len: u32,
        #[source]
        source: std::io::Error,
    },

    #[error(
        "internal invariant violated: {detail} (qname='{qname}' chrom_id={chrom_id} pos={pos})"
    )]
    Internal {
        detail: String,
        qname: String,
        chrom_id: u32,
        pos: u32,
    },

    #[error("malformed PreparedRead at qname='{qname}' (chrom_id={chrom_id}, pos={pos}): {reason}")]
    MalformedRead {
        reason: String,
        qname: String,
        chrom_id: u32,
        pos: u32,
    },
}
