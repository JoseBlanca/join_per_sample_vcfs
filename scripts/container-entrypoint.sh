#!/usr/bin/env bash
# Container entrypoint: refresh @anthropic-ai/claude-code when the user is
# starting a Claude session (interactive bash with no extra args, or running
# `claude` directly). One-off invocations like `./scripts/dev.sh cargo test`
# skip the update so they stay fast.
set -euo pipefail

should_update=false
if [ "$#" -eq 1 ] && [ "$1" = "bash" ]; then
    should_update=true
elif [ "${1:-}" = "claude" ]; then
    should_update=true
fi

if $should_update; then
    npm install -g --silent @anthropic-ai/claude-code@latest \
        || echo "warning: could not update claude-code (offline?)" >&2
fi

exec "$@"
