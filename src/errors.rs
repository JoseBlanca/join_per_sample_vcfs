use thiserror::Error;

#[derive(Error, Debug)]
pub enum VcfParseError {
    #[error("Invalid allele '{allele}'")]
    InvalidAllele { allele: String },

    #[error("Insufficient columns in VCF line: '{line}'")]
    NotEnoughColumns { line: String },

    #[error("Insufficient columns in CHROM header line")]
    NotEnoughColumnsInChromLine,

    #[error("Invalid position value '{value}' in line: '{line}'")]
    InvalidPosition { value: String, line: String },

    #[error("Invalid quality value '{value}': {line}")]
    InvalidQuality { value: String, line: String },

    #[error("Missing GT field in sample '{sample}' in line '{line}'")]
    MissingGtField { sample: String, line: String },

    #[error("FORMAT column (#8) not found in line '{line}'")]
    FormatColumnNotFound { line: String },

    #[error("GT field not found in FORMAT column in line '{line}'")]
    MissingGtFieldInFormat { line: String },

    #[error("Not possible to extract ploidy from line '{line}'")]
    ErrorFindingPloidy { line: String },

    #[error("Inconsistent ploidies found in line '{line}'")]
    InconsistentPloidies { line: String },

    #[error("Observed ({observed}) and given ({given}) ploidies are different line '{line}'")]
    DifferentObservedPloidy {
        line: String,
        observed: usize,
        given: usize,
    },

    #[error("I/O error: {source}")]
    Io {
        #[from]
        source: std::io::Error,
    },

    #[error("I/O error creating the ThreadPool to decompress the VCF file")]
    ThreadPoolError,

    #[error("I/O error opening path: '{path}'")]
    PathError { path: String },

    #[error("Magic byte error for file '{path}': {reason}")]
    MagicByteError { path: String, reason: String },

    #[error("Gzip in stdin is not supported")]
    GzipInStdinNotSupported,

    #[error("VCF file should be gzipped: '{path}'")]
    VCFFileShouldBeGzipped { path: String },

    #[error("VCF file should be bgzipped")]
    VCFFileShouldBeBGzipped,

    #[error("gVCF line has not enough fields to be a variant")]
    GVCFLineNotEnoughFields,

    #[error("Malformed VCF line: {reason}. Line: '{line}'")]
    MalformedLine { reason: String, line: String },

    #[error("Broken VCF header at line: '{line}'")]
    BrokenHeader { line: String },

    #[error("Malformed header fields and sample definition line")]
    MalformedHeader,

    #[error("Allele index overflow in genotype '{gt}': allele number exceeds maximum allowed (127)")]
    AlleleIndexOverflow { gt: String },

    #[error("RuntimeError: {message}")]
    RuntimeError { message: String },

    #[error("Variant group has no alternative alleles after merging")]
    NotVariable,
}
