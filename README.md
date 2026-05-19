# merge_per_sample_vcfs

Multi-sample SNP caller built around a six-stage pipeline:
per-sample pileup → `.psp` artefact → DUST filter → variant grouping →
per-group merger → posterior engine.

> The project was originally a gVCF merger and is mid-pivot to the
> caller described above. The gVCF-merger code has been removed; the
> tree now contains only the new pipeline. The CLI surface (`pileup`,
> `psp-to-pileup`, `var-calling`, `var-calling-from-bam`,
> `estimate-contamination`) is shipped via the `pop_var_caller`
> binary in [src/main.rs](src/main.rs).

## Authoritative documentation

- **Pipeline architecture (spec):**
  [doc/devel/specs/calling_pipeline_architecture.md](doc/devel/specs/calling_pipeline_architecture.md)
- **Design principles:**
  [doc/devel/specs/design_principles.md](doc/devel/specs/design_principles.md)
- **Current implementation status, per stage:**
  [PROJECT_STATUS.md](PROJECT_STATUS.md)
- **Per-feature implementation plans:**
  [doc/devel/implementation_plans/](doc/devel/implementation_plans/)
- **Implementation and review reports:**
  [doc/devel/reports/](doc/devel/reports/)

## Development container

A [Containerfile](Containerfile) is provided for reproducible
development with [podman](https://podman.io/). The image pins Rust
and includes `bcftools` / `tabix` for inspecting VCF output, plus the
[Claude Code](https://claude.com/claude-code) CLI for in-container
agent sessions.

```
./scripts/dev.sh                 # interactive shell inside the container
./scripts/dev.sh cargo test      # run a one-off command
./scripts/dev.sh cargo build --release
```

The first invocation builds the image (~4 min); subsequent runs reuse
the cached image.

### Running Claude Code in the container

```
./scripts/claude.sh              # new Claude Code session
./scripts/claude.sh --continue   # resume last session
```

The container acts as the sandbox: Claude can freely
read/write/execute inside the mounted project directory but cannot
touch anything else on the host, so `--dangerously-skip-permissions`
is passed by default. Host `~/.claude.json` and `~/.claude/` are
mounted so login and per-project auto-memory persist.

### Design notes

- The project directory is bind-mounted at the **same absolute path**
  inside and outside the container. This keeps paths consistent for
  tools that embed `cwd` into state (notably Claude auto-memory).
- Container builds write to `target-container/` instead of `target/`
  so host and container don't invalidate each other's incremental
  compilation state.
- Rootless podman maps container UID 0 to the host user, so files
  created inside the container appear with the correct ownership on
  the host.
