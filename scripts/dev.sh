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
#
# Supported runtimes (picked in this order):
#   1. podman          — Linux dev box (rootless OK).
#   2. container       — Apple's container CLI on macOS
#                        (https://github.com/apple/container).
set -euo pipefail

IMAGE="${IMAGE:-pop-var-caller-dev}"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

if command -v podman >/dev/null 2>&1; then
    RUNTIME=podman
elif command -v container >/dev/null 2>&1; then
    RUNTIME=container
    # Apple container needs its apiserver running; fail loudly if it isn't.
    if ! container system status >/dev/null 2>&1; then
        echo "Error: Apple 'container' apiserver is not running." >&2
        echo "Start it once per boot with:  container system start" >&2
        exit 1
    fi
else
    echo "Error: neither 'podman' nor 'container' (Apple) found on PATH." >&2
    exit 1
fi

image_exists() {
    case "$RUNTIME" in
        podman)    podman image exists "$IMAGE" ;;
        container) container image inspect "$IMAGE" >/dev/null 2>&1 ;;
    esac
}

if ! image_exists; then
    "$RUNTIME" build -t "$IMAGE" -f "$PROJECT_DIR/Containerfile" "$PROJECT_DIR"
fi

# Podman needs the :z SELinux relabel suffix on RHEL/Fedora hosts; Apple
# container uses virtio-fs and has no equivalent (the suffix is rejected).
MOUNT_SPEC="$PROJECT_DIR:$PROJECT_DIR"
[[ "$RUNTIME" == podman ]] && MOUNT_SPEC="$MOUNT_SPEC:z"

# Opt-in extra read-only bind mount for data that lives outside the
# project tree (e.g. reference genomes at $HOME/genomes). Same path
# inside the container as on the host, matching the convention used
# for the project mount. Single path only; rerun with a different
# value if you need a different mount.
EXTRA_MOUNTS=()
if [[ -n "${DEV_EXTRA_MOUNT:-}" ]]; then
    if [[ ! -d "$DEV_EXTRA_MOUNT" ]]; then
        echo "Error: DEV_EXTRA_MOUNT='$DEV_EXTRA_MOUNT' is not a directory." >&2
        exit 1
    fi
    EXTRA_OPTS="ro"
    [[ "$RUNTIME" == podman ]] && EXTRA_OPTS="${EXTRA_OPTS},z"
    EXTRA_MOUNTS+=(-v "${DEV_EXTRA_MOUNT}:${DEV_EXTRA_MOUNT}:${EXTRA_OPTS}")
fi

# Only request a TTY when stdin actually is one; -t fails when invoked
# from non-interactive contexts (CI, editor-spawned tooling).
TTY_FLAGS=(-i)
[[ -t 0 ]] && TTY_FLAGS+=(-t)

# Apple container runs the workload in a Linux VM with a small default
# memory cap (~2 GiB) that OOM-kills rustc/ld on this project. Bump
# memory + CPUs; podman inherits the host's resources directly and
# needs neither flag.
RESOURCE_FLAGS=()
if [[ "$RUNTIME" == container ]]; then
    RESOURCE_FLAGS+=(--memory "${DEV_MEM:-16g}" --cpus "${DEV_CPUS:-8}")
fi

# Empty-array expansions trip `set -u` on bash 3.2 (macOS default).
# The `${arr[@]+"${arr[@]}"}` form expands to zero args when the
# array is empty rather than erroring. TTY_FLAGS is always populated
# with at least `-i` so doesn't need the guard.
#
# HOME is propagated from the host so scripts that compute paths via
# `Path.home()` resolve to the same location as on the host (used by
# DEV_EXTRA_MOUNT consumers like the perf scripts, which default the
# reference path to $HOME/genomes/...).
exec "$RUNTIME" run --rm "${TTY_FLAGS[@]}" ${RESOURCE_FLAGS[@]+"${RESOURCE_FLAGS[@]}"} \
    -v "$MOUNT_SPEC" ${EXTRA_MOUNTS[@]+"${EXTRA_MOUNTS[@]}"} \
    -w "$PROJECT_DIR" \
    -e CARGO_TARGET_DIR="$PROJECT_DIR/target-container" \
    -e HOME="$HOME" \
    "$IMAGE" \
    "$@"
