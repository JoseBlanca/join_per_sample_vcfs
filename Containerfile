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

# coz: causal profiler — built for the cohort var-calling pipeline's
# producer->workers->writer shape, where "speed up which line actually
# improves end-to-end throughput?" is the question wall-clock and self-time
# profiles can't answer. Crucially it samples on PERF_TYPE_SOFTWARE /
# PERF_COUNT_SW_TASK_CLOCK, so it works inside the Apple-`container` Linux VM
# where the PMU is NOT virtualized (hardware counters, `perf c2c`, and
# `perf sched` off-CPU all return <not supported> there). Built from source
# pinned to a verified commit; CMake fetches LIEF + libelfin. Installs
# libcoz.so + the `coz` CLI to /usr/local. cmake is the only added build dep
# (build-essential / pkg-config are already present above).
RUN apt-get update \
    && apt-get install -y --no-install-recommends cmake \
    && git clone --recursive https://github.com/plasma-umass/coz /tmp/coz \
    && git -C /tmp/coz checkout 10630c542bc6d8a24595fd8283bea33bef892016 \
    && cmake -S /tmp/coz -B /tmp/coz/build -DCMAKE_BUILD_TYPE=Release \
    && cmake --build /tmp/coz/build -j \
    && cmake --install /tmp/coz/build \
    && ldconfig \
    && rm -rf /tmp/coz /var/lib/apt/lists/*

# TRF-mod (lh3): Tandem Repeats Finder with a BED-like output format — the
# SSR caller's Stage 0 (`ssr-catalog`) detection engine, shelled out per
# contig (see doc/devel/architecture/ssr_catalog.md §2). Single-file C build
# (`make -f compile.mak`, only build-essential needed), pinned to the commit
# the project vendors under TRF-mod/ so the catalog header's recorded
# `trf_mod_version` is reproducible. Installs to /usr/local/bin/trf-mod (on
# PATH for the catalog's layered binary discovery). Own layer so it does not
# invalidate the slow cargo-install / GATK layers above.
ARG TRF_MOD_COMMIT=3e891db310124f7e5f7a630a1c006650be9d1f3a
RUN git clone https://github.com/lh3/TRF-mod /tmp/trf-mod \
    && git -C /tmp/trf-mod checkout "${TRF_MOD_COMMIT}" \
    && make -C /tmp/trf-mod -f compile.mak \
    && install -m 0755 /tmp/trf-mod/trf-mod /usr/local/bin/trf-mod \
    && rm -rf /tmp/trf-mod

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
