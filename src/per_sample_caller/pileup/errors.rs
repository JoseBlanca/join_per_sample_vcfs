//! Errors raised by the pileup walker. Each variant carries enough
//! context (qname, chromosome, position) to point at the offending
//! input from the error alone, per `ia/specs/design_principles.md`
//! principle 6 ("typed errors at module boundaries").

use thiserror::Error;

#[derive(Error, Debug)]
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
        "phase-chain slot pool exhausted (cap={cap}) at chrom_id={chrom_id} pos={pos}; \
         consider raising --max-active-chain-slots"
    )]
    SlotExhausted { cap: u32, chrom_id: u32, pos: u32 },

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

    #[error("walker output channel closed before all records were drained: {context}")]
    ChannelClosed { context: String },
}
