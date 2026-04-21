#!/usr/bin/env bash
# Open a shell (or run a command) inside the project's development container.
#
# Usage:
#   ./scripts/dev.sh                 # interactive bash
#   ./scripts/dev.sh cargo test      # run a one-off command
#
# The project directory is bind-mounted at the same absolute path inside the
# container as on the host. Keeping paths identical avoids surprises with
# tools that embed cwd into state files (Claude auto-memory, cargo target
# fingerprints, editor configs).
set -euo pipefail

IMAGE="${IMAGE:-merge_per_sample_vcfs-dev}"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

if ! podman image exists "$IMAGE"; then
    podman build -t "$IMAGE" -f "$PROJECT_DIR/Containerfile" "$PROJECT_DIR"
fi

exec podman run --rm -it \
    -v "$PROJECT_DIR:$PROJECT_DIR:z" \
    -w "$PROJECT_DIR" \
    -e CARGO_TARGET_DIR="$PROJECT_DIR/target-container" \
    "$IMAGE" \
    "$@"
