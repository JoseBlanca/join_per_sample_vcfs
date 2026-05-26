# Development container for pop_var_caller.
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
#   - samtools / bcftools / tabix: inspect and compare BAM/VCF output
#     against reference tools.
#   - freebayes: reference variant caller; one of the three we compare
#     against in benchmarks/tomato1/.
#   - openjdk-17-jre-headless: Java runtime for GATK (installed below).
#   - unzip / wget: needed to fetch and unpack the GATK release.
#   - python3-psutil: lets the perf experiment scripts measure
#     subprocess RSS without needing uv inside the container.
#   - curl: bootstraps the NodeSource apt repo for Node.js.
#   - linux-perf: sampling profiler; backs cargo-flamegraph. In-container
#     sampling needs the host's kernel.perf_event_paranoid relaxed
#     (typically to 1). In rootless podman, --cap-add=SYS_ADMIN/PERFMON
#     does NOT help — the syscall still reaches the host kernel as the
#     unprivileged invoking user, so the sysctl is the only knob.
#   - hyperfine: statistical CLI benchmarking for reference comparisons.
#   - valgrind: callgrind/cachegrind for instruction-level profiling.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        pkg-config \
        git \
        ca-certificates \
        samtools \
        bcftools \
        tabix \
        freebayes \
        openjdk-17-jre-headless \
        unzip \
        wget \
        python3-psutil \
        curl \
        linux-perf \
        hyperfine \
        valgrind \
    && curl -fsSL https://deb.nodesource.com/setup_22.x | bash - \
    && apt-get install -y --no-install-recommends nodejs \
    && rm -rf /var/lib/apt/lists/*

# GATK 4: HaplotypeCaller + GenotypeGVCFs reference caller. Distributed
# as a Java fat-jar with a shell wrapper. Installs to /opt/gatk-<ver>/
# with a stable /opt/gatk -> /opt/gatk-<ver>/ symlink so the runner
# scripts (default `GATK_BIN=/opt/gatk/gatk`) stay valid across version
# bumps. Bump GATK_VERSION when you want to upgrade.
ARG GATK_VERSION=4.6.0.0
RUN wget -q "https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip" -O /tmp/gatk.zip \
    && unzip -q /tmp/gatk.zip -d /opt/ \
    && ln -s "/opt/gatk-${GATK_VERSION}" /opt/gatk \
    && rm /tmp/gatk.zip

# Rust components not included in the base image's minimal profile.
RUN rustup component add rustfmt clippy

# cargo-flamegraph: wraps perf to render flamegraphs of release builds.
# Installed as its own layer so it stays cached across Cargo.toml changes.
RUN cargo install flamegraph --locked

# samply: sampling profiler that produces a profile viewable in the Firefox
# profiler UI; complementary to flamegraph. On Linux it still uses
# perf_event_open under the hood, so the same host-side
# kernel.perf_event_paranoid constraint as cargo-flamegraph applies.
RUN cargo install samply --locked

# cargo-show-asm: emit annotated assembly for a single function. Used to
# verify codegen-level claims (bounds-check elision, autovectorization)
# without disassembling the whole binary.
RUN cargo install cargo-show-asm --locked

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

# GATK's `gatk` wrapper script uses `#!/usr/bin/env python`, but
# Debian ships only `python3` (no plain `python` symlink). Installed
# as its own layer so adding it doesn't invalidate the (slow) cargo
# install + GATK download layers above. If the Containerfile is ever
# rewritten, fold this back into the main apt-get block.
RUN apt-get update \
    && apt-get install -y --no-install-recommends python-is-python3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work
COPY scripts/container-entrypoint.sh /usr/local/bin/container-entrypoint.sh
ENTRYPOINT ["/usr/local/bin/container-entrypoint.sh"]
CMD ["bash"]
