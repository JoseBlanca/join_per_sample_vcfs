# Pop Var Caller

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

## Memory tuning

The cohort `var-calling` driver opens one `PspReader` per sample per
chromosome worker, each holding one decoded block live, so peak heap
scales as

```
peak_RSS ≈ ~5 GB + n_threads × N × per_block × 7
```

where `N` is the number of samples, `n_threads` is the rayon worker
count (capped at `n_chromosomes`), and `per_block` is roughly
`block-target-bytes × 7` (the writer's projected-byte threshold
under-counts the live in-RAM footprint of decoded columns by about
that factor).

The PSP writer's `--block-target-bytes` knob on `pop_var_caller
pileup` is the lever for this. Lower values cut peak heap; higher
values shrink the .psp on disk. Wall time is essentially flat across
the entire range. Sweep on tomato1 (N=18, T=4, real per-sample PSPs):

| `--block-target-bytes` | peak cohort RSS | cohort .psp size | wall  |
|------------------------|----------------:|-----------------:|------:|
| 16 MiB                 |         2501 MB |           183 MB | 10.4s |
|  4 MiB                 |          757 MB |           189 MB | 10.3s |
|  1 MiB (default)       |          261 MB |           206 MB | 10.6s |
|  512 KiB               |          161 MB |           233 MB | 10.6s |
|  256 KiB               |          108 MB |           246 MB | 10.9s |
|   64 KiB               |           59 MB |           321 MB | 10.8s |

Picking a value:

- **Single-sample / archival** workflows pay only the on-disk axis —
  dial up (4 MiB, 16 MiB) for the smallest .psp.
- **Large-cohort joint genotyping** (N in the hundreds to thousands)
  pays the cohort-step memory axis — start from the default and dial
  down if the formula above projects past the memory budget.

The knob is exposed only on `pileup` (it's a writer-side setting);
the cohort subcommands read whatever block size each input was
written with, no extra flag needed. The accepted range is 16 KiB
to 16 MiB — values outside that are rejected at parse time. Below
16 KiB zstd loses meaningful compression context; above 16 MiB sits
the legacy hardcoded default that motivated this knob in the first
place. The full per-stage rationale is in
[doc/devel/specs/per_sample_pileup_format.md §"Block sizing"](doc/devel/specs/per_sample_pileup_format.md).

## License

Licensed under the MIT License — see [LICENSE](LICENSE) for the full
text.
