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

step "1/5  dependency guard: production must not depend on ng"
# ng is a from-scratch experiment caller; production (everything outside
# src/ng/) must never `use crate::ng`. The exp binary
# (src/pop_var_caller_exp/) is a driver for ng, so it is allowed to depend
# on ng and is excluded too — the dependency points the safe way (see
# doc/devel/ng/spec/typed_regions_cli.md T7). The excluded paths are
# load-bearing: widen this list and the guard becomes a comment.
guard_hits=$(grep -rn "use crate::ng" src --include="*.rs" \
    | grep -Ev "^src/(ng|pop_var_caller_exp)/" || true)
if [ -n "$guard_hits" ]; then
    printf "\033[1;31mguard failed: production code depends on ng:\033[0m\n%s\n" "$guard_hits"
    exit 1
fi

step "2/5  cargo fmt --check"
cargo fmt --check

step "3/5  cargo clippy --all-targets --all-features -- -D warnings"
cargo clippy --all-targets --all-features -- -D warnings

step "4/5  cargo test"
cargo test

step "5/5  cargo bench --no-run  (compile-only, catches bench bitrot)"
cargo bench --no-run

printf "\n\033[1;32mAll checks passed.\033[0m\n"
'
