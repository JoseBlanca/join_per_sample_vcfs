# Data layout & cache locality checklist

**Purpose.** Reduce CPU stalls from L2/L3 misses by packing the data the hot loop touches into the smallest possible footprint per access, and by avoiding cross-thread cache-line invalidation on shared atomics.

**Triggers.** Hot loops over `Vec<Struct>`; structs with many fields where loops read only one or two; atomics or `Mutex` accessed by multiple threads; `pub struct` whose layout matters for callers; very large structs (> ~64 bytes) on the hot path.

**Skip when.** Code does not iterate over collections at scale and has no shared-state concurrency. If dispatched anyway, write `No findings.`

## Rules

- **The L1 (~1 ns) ↔ main-memory (~100 ns) gap is the budget.** Refactors that get a hot loop into L1-resident data can dwarf any micro-tweak. When a profile shows high cache-miss rates, structural changes to data are the right tool.
- **Cache lines are typically 64 bytes; assume 128 on x86_64 due to the spatial prefetcher.** Hot-loop structs that span more than one cache line cost double on every iteration when the loop touches the trailing fields. Filed for any `Vec<S>` where `mem::size_of::<S>() > 64` and the loop reads more than one field of `S`.
- **Project hot fields into a cache struct (data-oriented design).** When sorters / filters / loops repeatedly read the same handful of fields off a larger struct, build a `RowCache { f1, f2, f3 }` of just those fields, refresh it once per update, and run the loop over the cache. The mechanism: collapse the loop's working set so it fits in L1 and the per-iteration body stops chasing pointers. Cross-reference: `allocations.md` (the cache projection also eliminates the per-comparison clone).
- **Array-of-Structs vs Struct-of-Arrays is a function of access pattern.**
  - AoS (`Vec<S>`) is correct when the hot loop reads several fields of each element together.
  - SoA (`struct Grid { a: Vec<A>, b: Vec<B> }`) is correct when each loop reads exactly one field; SoA halves cache footprint and unlocks autovectorization.
  - A `Vec<S>` whose hot loops read only `s.a` is filed Likely with the SoA experiment.
- **Order struct fields large-to-small only when `repr(C)` requires it.** `repr(Rust)` already reorders for minimum padding. Findings here are limited to `repr(C)` structs that gratuitously waste padding without an FFI / serialization reason.
- **`#[repr(transparent)]` newtypes are free.** Use them for domain types over primitives in hot loops; they cost nothing and let the compiler optimize as if it were the underlying primitive.
- **Niche optimization is real.** `Option<&T>`, `Option<NonZeroU32>`, `Option<Box<T>>` are the same size as their inner type. Forcing `Option<u32>` (8 bytes due to the discriminant) where `Option<NonZeroU32>` (4 bytes) would do is a finding for hot, large arrays.
- **False sharing kills atomics-heavy code.** Two thread-private atomics that happen to land on the same cache line cause every write on one to invalidate the line on the other CPU's L1. Wrap shared, frequently-written atomics in `crossbeam_utils::CachePadded<T>` to put each on its own cache line. **Pad to 128 bytes on x86_64** because of the prefetcher. Filed for any per-thread counter or shard array of atomics not cache-padded. Cross-reference: `concurrency.md` on per-thread counters.
- **Inline small keys in maps.** A `HashMap<u64, V>` is dramatically smaller and faster than `HashMap<String, V>` if the key has a natural integer encoding. Hashes of small fixed-size types vectorize. The choice of *hasher* lives in `hot_loops.md`.
- **Pointer chases fragment the cache.** `Vec<Box<Node>>` requires a pointer chase per element. If the elements are uniform, store inline (`Vec<Node>`); if heterogeneous and the type set is closed, use an enum.
