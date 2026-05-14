//! `pop_var_caller` binary entry point. Parses the top-level CLI and
//! dispatches to the subcommand orchestrator. Subcommand logic lives
//! in `pop_var_caller/cli.rs`; this file is intentionally thin.

use std::process;

use clap::Parser;
use merge_per_sample_vcfs::pop_var_caller::{Cli, PopVarCallerCommand, run_pileup};

fn main() {
    let cli = Cli::parse();
    let result = match cli.cmd {
        PopVarCallerCommand::Pileup(args) => run_pileup(&args).map_err(|e| format!("{e}")),
    };
    if let Err(msg) = result {
        eprintln!("error: {msg}");
        process::exit(1);
    }
}
