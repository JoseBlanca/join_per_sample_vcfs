#!/usr/bin/env bash
# Pre-commit checks: fmt, clippy, tests, bench-compile. Runs inside the
# dev container via scripts/dev.sh.
#
# Manual:       ./scripts/precommit-check.sh
# As git hook:  ln -s ../../scripts/precommit-check.sh .git/hooks/pre-commit
# Bypass once:  git commit --no-verify
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

"$SCRIPT_DIR/dev.sh" bash -c '
set -euo pipefail
step() { printf "\n\033[1;34m==> %s\033[0m\n" "$1"; }

step "1/4  cargo fmt --check"
cargo fmt --check

step "2/4  cargo clippy --all-targets --all-features -- -D warnings"
cargo clippy --all-targets --all-features -- -D warnings

step "3/4  cargo test"
cargo test

step "4/4  cargo bench --no-run  (compile-only, catches bench bitrot)"
cargo bench --no-run

printf "\n\033[1;32mAll checks passed.\033[0m\n"
'
