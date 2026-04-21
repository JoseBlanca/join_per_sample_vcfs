#!/usr/bin/env bash
# Launch Claude Code inside the project's development container.
#
# The container is the sandbox: Claude can freely read/write/execute within
# the mounted project directory, but cannot touch anything else on the host
# (no home dir, no other projects, no ssh keys). This makes
# --dangerously-skip-permissions safe to use here.
#
# Claude credentials and session state are mounted from the host so that the
# in-container agent shares your login and per-project auto-memory.
#
# Usage:
#   ./scripts/claude.sh                # start a new Claude session
#   ./scripts/claude.sh --continue     # continue the most recent session
set -euo pipefail

IMAGE="${IMAGE:-merge_per_sample_vcfs-dev}"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

if ! podman image exists "$IMAGE"; then
    podman build -t "$IMAGE" -f "$PROJECT_DIR/Containerfile" "$PROJECT_DIR"
fi

# Mount Claude host state if present. The container runs as root (UID 0
# mapped to the host user by rootless podman), so the mount target is /root.
CLAUDE_MOUNTS=()
[[ -f "$HOME/.claude.json" ]] && CLAUDE_MOUNTS+=(-v "$HOME/.claude.json:/root/.claude.json:z")
[[ -d "$HOME/.claude"      ]] && CLAUDE_MOUNTS+=(-v "$HOME/.claude:/root/.claude:z")

exec podman run --rm -it \
    -v "$PROJECT_DIR:$PROJECT_DIR:z" \
    -w "$PROJECT_DIR" \
    "${CLAUDE_MOUNTS[@]}" \
    -e ANTHROPIC_API_KEY \
    -e CARGO_TARGET_DIR="$PROJECT_DIR/target-container" \
    "$IMAGE" \
    claude --dangerously-skip-permissions "$@"
