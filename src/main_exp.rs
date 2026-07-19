//! `pop_var_caller_exp` binary entry point — ng's experiment command
//! surface. Parses the top-level CLI and dispatches to the subcommand
//! driver; the logic lives in `pop_var_caller_exp/`, so this file is
//! intentionally thin, the way `src/main.rs` is.
//!
//! Error rendering is a plain Display here; Milestone A2 hoists
//! `format_error_chain` into the library and both binaries share it.

use std::process;

use clap::Parser;
use pop_var_caller::pop_var_caller_exp::{Cli, PopVarCallerExpCommand, run_typed_regions};

fn main() {
    let cli = Cli::parse();
    let result = match cli.cmd {
        PopVarCallerExpCommand::TypeRegions(args) => {
            run_typed_regions(&args).map_err(|e| e.to_string())
        }
    };
    if let Err(msg) = result {
        eprintln!("error: {msg}");
        process::exit(1);
    }
}
