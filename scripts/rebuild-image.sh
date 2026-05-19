#!/usr/bin/env bash
# Rebuild the project's development container image.
#
# `dev.sh` and `claude.sh` only build the image when it doesn't already exist,
# so changes to the Containerfile (or to Cargo.toml/Cargo.lock, which the
# image pre-fetches) won't be picked up automatically. Run this script
# whenever you've edited any of those and want the next `dev.sh` /
# `claude.sh` invocation to use the updated image.
#
# Usage:
#   ./scripts/rebuild-image.sh              # rebuild, reusing layer cache
#   ./scripts/rebuild-image.sh --no-cache   # full rebuild from scratch
#   ./scripts/rebuild-image.sh --pull       # also refresh the FROM base image
#
# Any extra arguments are forwarded to `podman build`.
set -euo pipefail

IMAGE="${IMAGE:-pop-var-caller-dev}"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

exec podman build \
    -t "$IMAGE" \
    -f "$PROJECT_DIR/Containerfile" \
    "$@" \
    "$PROJECT_DIR"
