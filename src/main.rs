//! `pop_var_caller` binary entry point. Parses the top-level CLI and
//! dispatches to the subcommand orchestrator. Subcommand logic lives
//! in `pop_var_caller/cli.rs`; this file is intentionally thin.

use std::error::Error;
use std::process;

// Opt-in `mimalloc` global allocator (`--features alloc-mimalloc`, off by
// default). The cohort `var-calling` path is allocation-churn-heavy (per-allele
// record building across the worker pool); on a T=8 / N=50 tomato cohort
// mimalloc cut peak RSS ~9% and shaved wall time vs the system allocator,
// byte-identical. Off by default so production builds carry no vendored-C dep.
#[cfg(feature = "alloc-mimalloc")]
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

use clap::Parser;
use pop_var_caller::pop_var_caller::{
    Cli, PopVarCallerCommand, run_estimate_contamination, run_pileup, run_psp_to_pileup,
    run_ssr_catalog, run_ssr_pileup, run_var_calling,
};

// Walk the `std::error::Error::source()` chain joining messages with
// `: ` so the user sees the leaf cause, not just the outermost wrapper.
// Several of our errors (and noodles') use the `thiserror` `: {source}`
// pattern that already embeds the child's message in the parent's
// Display; skip a level when its message is already a substring of the
// previous one to avoid `invalid record: invalid record: invalid value`
// noise.
fn format_error_chain(err: &(dyn Error + 'static)) -> String {
    let mut out = err.to_string();
    let mut prev = out.clone();
    let mut cur = err.source();
    while let Some(e) = cur {
        let msg = e.to_string();
        if !prev.contains(&msg) {
            out.push_str(": ");
            out.push_str(&msg);
        }
        prev = msg;
        cur = e.source();
    }
    out
}

fn main() {
    let cli = Cli::parse();
    let result = match cli.cmd {
        PopVarCallerCommand::Pileup(args) => run_pileup(&args).map_err(|e| format_error_chain(&e)),
        PopVarCallerCommand::PspToPileup(args) => {
            run_psp_to_pileup(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::EstimateContamination(args) => {
            run_estimate_contamination(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::VarCalling(args) => {
            run_var_calling(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::SsrCatalog(args) => {
            run_ssr_catalog(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::SsrPileup(args) => {
            run_ssr_pileup(&args).map_err(|e| format_error_chain(&e))
        }
    };
    if let Err(msg) = result {
        eprintln!("error: {msg}");
        process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::format_error_chain;
    use std::error::Error;
    use std::fmt;

    #[derive(Debug)]
    struct TestErr {
        msg: &'static str,
        src: Option<Box<TestErr>>,
    }

    impl fmt::Display for TestErr {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            f.write_str(self.msg)
        }
    }

    impl Error for TestErr {
        fn source(&self) -> Option<&(dyn Error + 'static)> {
            self.src.as_deref().map(|e| e as &(dyn Error + 'static))
        }
    }

    fn err(msg: &'static str, src: Option<TestErr>) -> TestErr {
        TestErr {
            msg,
            src: src.map(Box::new),
        }
    }

    #[test]
    fn joins_distinct_messages_with_colon() {
        let leaf = err("duplicate tag: PL", None);
        let mid = err("invalid read group", Some(leaf));
        let top = err("CRAM input", Some(mid));
        assert_eq!(
            format_error_chain(&top),
            "CRAM input: invalid read group: duplicate tag: PL"
        );
    }

    #[test]
    fn skips_level_already_embedded_in_parent() {
        // thiserror's `: {source}` pattern bakes the child's Display
        // into the parent's. Walking the chain naively would then
        // re-print the child — the dedup must drop it.
        let leaf = err("inner cause", None);
        let outer = err("outer prefix: inner cause", Some(leaf));
        assert_eq!(format_error_chain(&outer), "outer prefix: inner cause");
    }

    #[test]
    fn surfaces_deeper_message_when_intermediate_was_a_duplicate() {
        // Mirrors the real noodles chain: an outer error whose
        // Display embeds a generic middle ("invalid record") while
        // the actual useful detail lives one level deeper.
        let leaf = err("duplicate tag: PL", None);
        let mid = err("invalid record", Some(leaf));
        let top = err("failed to open CRAM 'x.cram': invalid record", Some(mid));
        assert_eq!(
            format_error_chain(&top),
            "failed to open CRAM 'x.cram': invalid record: duplicate tag: PL"
        );
    }
}
