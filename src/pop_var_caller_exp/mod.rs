//! `pop_var_caller_exp` binary namespace — ng's experiment command
//! surface. A second binary in the same crate, kept apart from the
//! production `pop_var_caller` CLI so ng's experiment knobs never grow it
//! (see `doc/devel/ng/spec/typed_regions_cli.md` §2). The library both
//! binaries link is the same one; only what a user can *invoke* is split.
//!
//! Layout mirrors [`crate::pop_var_caller`]: [`cli`] owns the top-level
//! `Parser` plus the subcommand enum, and one module per subcommand owns
//! its `Args`, its `run_*`, and its `#[non_exhaustive]` error enum.

pub mod cli;
pub mod typed_regions;

pub use cli::{Cli, PopVarCallerExpCommand};
pub use typed_regions::{TypedRegionsArgs, TypedRegionsCliError, run_typed_regions};
