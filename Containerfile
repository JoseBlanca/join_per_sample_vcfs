# Development container for merge_per_sample_vcfs.
#
# Serves two purposes:
#   1. Reproducible build/dev environment for the team.
#   2. Sandbox for running Claude Code with broad permissions, isolated
#      from the host filesystem outside the project directory.
#
# Rust version is pinned to match the host toolchain used during development.
FROM docker.io/library/rust:1.95-bookworm

# System dependencies:
#   - build-essential / pkg-config: required by a few cargo crates that
#     may compile native code.
#   - git / ca-certificates: needed by cargo for git-based deps and HTTPS.
#   - bcftools / tabix: inspect and compare VCF output against reference tools.
#   - curl: bootstraps the NodeSource apt repo for Node.js.
#   - linux-perf: sampling profiler; backs cargo-flamegraph. In-container
#     sampling typically needs the host's perf_event_paranoid relaxed or
#     the container run with --cap-add=SYS_ADMIN.
#   - hyperfine: statistical CLI benchmarking for reference comparisons.
#   - valgrind: callgrind/cachegrind for instruction-level profiling.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        pkg-config \
        git \
        ca-certificates \
        bcftools \
        tabix \
        curl \
        linux-perf \
        hyperfine \
        valgrind \
    && curl -fsSL https://deb.nodesource.com/setup_22.x | bash - \
    && apt-get install -y --no-install-recommends nodejs \
    && rm -rf /var/lib/apt/lists/*

# Rust components not included in the base image's minimal profile.
RUN rustup component add rustfmt clippy

# cargo-flamegraph: wraps perf to render flamegraphs of release builds.
# Installed as its own layer so it stays cached across Cargo.toml changes.
RUN cargo install flamegraph --locked

# Claude Code CLI, used when the container hosts an agent session.
RUN npm install -g @anthropic-ai/claude-code

# Pre-warm the cargo registry cache with this project's dependencies so that
# cold-start builds don't re-download every crate. Copy only the manifests
# (not the source) so this layer is cached until a dep changes. cargo fetch
# requires a parseable manifest plus src/main.rs, but it would also demand
# a file for every [[bench]] target — so we strip the [[bench]] sections
# from the manifest copy before fetching. Only the populated registry under
# $CARGO_HOME survives into the final image; the stripped manifest, stub
# src/, and /build itself are byproducts of this layer.
WORKDIR /build
COPY Cargo.toml Cargo.lock ./
RUN mkdir -p src \
    && echo 'fn main() {}' > src/main.rs \
    && awk 'BEGIN{skip=0} /^\[\[bench\]\]/{skip=1; next} /^\[/{skip=0} !skip' Cargo.toml > Cargo.toml.fetch \
    && mv Cargo.toml.fetch Cargo.toml \
    && cargo fetch \
    && rm -rf src

WORKDIR /work
COPY scripts/container-entrypoint.sh /usr/local/bin/container-entrypoint.sh
ENTRYPOINT ["/usr/local/bin/container-entrypoint.sh"]
CMD ["bash"]
