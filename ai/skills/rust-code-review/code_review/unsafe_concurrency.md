# Unsafe & concurrency checklist

**Purpose.** Every `unsafe` block is justified and audited; every shared-state primitive is correctly used; every async lock or blocking call is examined for the cancellation, deadlock, and starvation hazards that humans and LLMs both miss.

**Triggers.** Search the in-scope code for: `unsafe`, `Arc`, `Mutex`, `RwLock`, atomic types, channel types (`mpsc`, `broadcast`, `oneshot`, `crossbeam-channel`), `async fn`, `.await`, `tokio::spawn`, `std::thread::spawn`, `tokio::select!`, `static mut`, `*mut`/`*const`, `Cell`, `RefCell`, FFI handles, `unsafe impl Send/Sync`.

**Skip when.** None of the above appear anywhere in the in-scope code. The orchestrator should not dispatch this category for pure-data, single-threaded, panic-free code; if dispatched anyway, write `No findings.`

## Rules

- **`unsafe` is a last resort.** Each `unsafe` block requires a justification explaining why the safe alternative is insufficient (benchmarked perf gap, FFI, or an `std` API that is itself `unsafe`). "Faster" without a measurement is a finding.
- **Every `unsafe` block has a `// SAFETY:` comment** that names (a) every precondition the called function requires, (b) the specific reason each precondition holds here, referencing the code that establishes it, and (c) for `unsafe fn` and `unsafe impl`, the caller obligations — documented in `///` docs as `# Safety`. "Pointer is valid" is not a safety comment.
- **`Send`/`Sync` are correct and minimal.** Manual `unsafe impl Send`/`Sync` requires a safety comment proving the type contains no trait-violating shared mutability. Auto-derived `Send`/`Sync` on types containing `*mut T`, `Cell`, `RefCell`, or FFI handles is a finding — explicit, justified impls (positive or negative) are required.
- **No `static mut`.** Mutable globals use `Mutex`, `RwLock`, `OnceLock`, or atomics. Never roll your own initialization with `unsafe`.
- **Inside `async` code:** no non-`Send` guard (`std::sync::MutexGuard`, `RwLockReadGuard`/`WriteGuard`, `RefMut`, `Ref`) held across `.await`. Use `tokio::sync` primitives or scope the `std` lock to a sync block that drops the guard before awaiting. No blocking I/O (`std::fs`, `std::net`, `std::thread::sleep`, blocking drivers) — use the runtime's primitives or `spawn_blocking`. Holding a non-`Send` guard across `.await` is Blocker.
- **Lock ordering.** When multiple locks are acquired anywhere in the crate, the crate-level docs declare a total ordering. Every acquisition site is reviewed against it — inconsistent ordering across call sites is a deadlock and a Blocker. Prefer designs that hold one lock at a time.
- **Atomics: every memory ordering is justified** in a comment naming the synchronization relationship (which load pairs with which store, what data is protected). Default `SeqCst`; weaker orderings require explicit argument. `Relaxed` used for synchronization (rather than independent counters) is a Blocker until proven correct.
- **Channel choice is justified.** Bounded vs. unbounded (unbounded is a memory hazard under load and requires justification), `mpsc` vs. `broadcast` vs. `oneshot`, sync vs. async. Sender and receiver drop semantics are handled on both ends.
- **Thread and async task panics are not silently ignored.** `JoinHandle::join` (threads) and `JoinHandle` / spawned-task results (async) are awaited and surfaced as typed errors. Detached spawns require justification.
- **Cancellation safety.** Functions used in `tokio::select!` arms or that may be cancelled are documented as cancellation-safe or not. State that must survive cancellation is held to allow resumption or rollback.
