//! `pop_var_caller` binary entry point. Parses the top-level CLI and
//! dispatches to the subcommand orchestrator. Subcommand logic lives
//! in `pop_var_caller/cli.rs`; this file is intentionally thin.

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
use pop_var_caller::error_render::format_error_chain;
use pop_var_caller::pop_var_caller::{
    Cli, PopVarCallerCommand, run_estimate_contamination, run_pileup, run_psp_to_pileup,
    run_ssr_call, run_ssr_catalog, run_ssr_pileup, run_var_calling,
};

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
        PopVarCallerCommand::SsrCall(args) => {
            run_ssr_call(&args).map_err(|e| format_error_chain(&e))
        }
    };
    if let Err(msg) = result {
        eprintln!("error: {msg}");
        process::exit(1);
    }
}
