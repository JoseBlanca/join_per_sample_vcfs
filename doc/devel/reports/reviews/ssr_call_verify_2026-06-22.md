# Adversarial re-verification — the seven "couldn't break" claims of `ssr_call_2026-06-21.md`

**Status:** independent skeptic's re-pass, 2026-06-22, branch `ssr-cohort`. Target:
[ssr_call_2026-06-21.md §3](ssr_call_2026-06-21.md) ("Things I tried to break but
couldn't"), re-verified against the **current** (post-amendment) docs and the actual
source — not the review's snapshot. Sources read in full: spec
[ssr_cohort_mark2.md](../../specs/ssr_cohort_mark2.md); architecture
[ssr_call_reading.md](../../architecture/ssr_call_reading.md),
[ssr_call_parameters.md](../../architecture/ssr_call_parameters.md),
[ssr_call_genotyping.md](../../architecture/ssr_call_genotyping.md);
[src/var_calling/posterior_engine.rs](../../../src/var_calling/posterior_engine.rs);
HipSTR `StutterAlignerClass.cpp`. Citations inline. I tried to refute each claim; what
survives my own refutation is reported as **holds**, what survives only with a carve-out
as **caveat**, what does not as **breaks**.

## Verdict table

| # | Claim | Verdict |
|---|---|---|
| 3 | Slip-placement marginalization is a faithful HipSTR borrow | **BREAKS for impure alleles** (holds for pure) — **RESOLVED 2026-06-22 (Fix 1: explicit `Σ_v` placement sum)** |
| 1 | Pre-pass reduce determinism is bit-identical across threads | **CAVEAT** — integer stats + `F_i`/level fixed-point hold; the M4-new floats were *asserted* order-independent but not established — **RESOLVED 2026-06-22 (verify-fix #1: fixed-point `ℓ_pen`/BIC + clustering tie-breaks)** |
| 7b | The engine's exact-AF convolution is reusable per-allele for m4 QUAL | **CAVEAT** — kernel reusable, but Step-1 hard-codes allele-0, and Σ-over-collapses is an approximation, not "exact/properly normalized" — **RESOLVED 2026-06-22 (verify-fix #7b: name the generalization + scope the approximation)** |
| 4 | Absorbing-trap guards are sound (`π=0` floor, `F=0.99` ceiling) | **CAVEAT** — `F` ceiling holds solidly; the geometric `G₀` pseudocount had *no hard floor* against f64 underflow at large \|Δ\| — **RESOLVED 2026-06-22 (verify-fix #4: explicit `G₀` floor / log-domain)** |
| 5 | Cohort recurrence defeats within-sample PCR jackpotting | **HOLDS** (headline) — **CLARIFIED 2026-06-22 (verify-fix #5): dedup is upstream (markdup) + Stage 1 skips flag `0x400`; depth duplicate-free by construction, no in-caller dedup; residual unmarked jackpotting out of scope** |
| 2 | Per-locus EM is a pure deterministic function of its inputs | **HOLDS** — conditional on #1's inputs being deterministic |
| 6 | Engine already takes a per-sample `F` vector; `F_i` loop needs no prior-side change | **HOLDS** |
| 7a | Engine generalizations (pseudocount swap, class-path bypass) are low-invasiveness | **HOLDS** — pinned to three exact paths |

---

## 3 — BREAKS for impure alleles: substitution-only in-tract `align` does *not* reproduce HipSTR's placement marginalization

> **RESOLVED 2026-06-22 (Fix 1 — restore the explicit placement sum).** `align(obs |
> cand ⊕ Δ)` is now defined as `Σ_{v ∈ placements(cand,Δ)} Pr(v) · align_subst(obs | v)`:
> an **explicit sum over the placement-distinct variants** of `cand ⊕ Δ` (the repeat runs
> separated by interruptions), with a **uniform position prior** and equal-LL placements
> within a homogeneous run **collapsed** (HipSTR's `upstream_match_lengths_` trick), so the
> term count is **≈ (#interruptions + 1), not (#copies)** — a pure allele is one term (the
> closed-form / exact-match fast path, no cost), a singly-interrupted impure allele ~2–3.
> The in-tract substitution-only rule (C1) applies *within each fixed variant `v`*. This is
> the faithful HipSTR borrow; a single substitution-only forward against one committed
> `cand⊕Δ` (which is what the genotyping doc had collapsed to) could not marginalize
> placement for impure alleles. Applied to spec §6 (`Σ_v` formula + in-tract rule + cost +
> §7 reachability), genotyping §3/§4/§8/§9 (new Q-G3b) + settled list, parameters §1.
> Determinism preserved: `Σ_v` is a pure, catalog-run-ordered reduction. Open (calibration,
> simulator): multi-interruption / compound placement enumeration; whether the position
> prior should ever be informed.

**Verdict: breaks** (for impure/interrupted alleles; holds for pure). This is the angle
the prompt flagged as most likely to break, and under C1 it does.

**What HipSTR actually does.** `StutterAlignerClass::align_pcr_deletion_reverse`
([StutterAlignerClass.cpp:106-150](../../../../pop_var_caller/HipSTR/src/SeqAlignment/StutterAlignerClass.cpp))
marginalizes *which copy* a slip lands on by an **explicit position sum**: the loop over
`i` (lines 127-142) walks every placement of the whole-period block-size change, accrues a
running log-prob, and `fast_log_sum_exp(log_probs_)` (line 149) sums them under a **uniform
position prior** `log_prior = -int_log(block_len_+D+1)` (line 112). Crucially, **within a
fixed `D` the scoring is substitution-only** — every term is `base_log_correct` /
`base_log_wrong` (match vs mismatch, lines 119-120, 129-130); there is no per-base gap
competing inside the block. So HipSTR's in-block model is *already* "substitutions at fixed
Δ" **plus a sum over placements**. C1 correctly borrowed the first half; the marginalization
**is** the second half.

**Why a single `align(obs | cand⊕Δ)` term loses it.** The genotyping doc writes the
likelihood as `Σ_Δ S_θ(Δ) × align(obs | candidate ⊕ Δ)` ([genotyping §4](../../architecture/ssr_call_genotyping.md)),
one `align` term per Δ, defined as a **banded forward against a single sequence** with
"in-tract substitutions only, per-base gaps confined to the flanks." A forward sums over
alignment *paths against one reference sequence*; it does **not** sum over *different
reference sequences*. For a **pure** allele every placement of a slip yields the *same*
string `cand⊕Δ`, so the position sum is degenerate and the single term is exactly right.
For an **impure** allele, deleting copy *j* vs copy *k* moves the interruption relative to
the deletion → **different strings**, which is precisely the case HipSTR's position sum
exists to average. The design keeps impure alleles first-class
([genotyping §2/§3](../../architecture/ssr_call_genotyping.md): "Impure (off-lattice) major
peaks are kept, first-class").

**The internal contradiction.** [genotyping §3](../../architecture/ssr_call_genotyping.md)
says placement is "marginalized inside the §6 alignment (HipSTR), so the data picks it,"
**but** the same section defines `⊕` as "B is A after k whole motif units added/removed
within A's own repeat context, **interruption structure preserved**" — i.e. `⊕` produces
**one canonical placement** (interruption-preserving), not a set. A substitution-only
forward against that one committed sequence cannot also "let the data pick" the placement.
The two sentences describe incompatible operations. The original review even credited this
as verified in its §1 ("slip-placement marginalized inside `align` — verified") and again in
the §3 "couldn't break" list — but that verification predates reconciling the borrow with
C1's *single-sequence, gapless-in-tract* `align` for **impure** alleles.

**Why it matters.** For an interrupted allele the committed-placement approximation scores
a mis-placed interruption as substitution **mismatches** under `ε`, inflating the apparent
per-base error and biasing the slip likelihood exactly for the alleles (interrupted ones)
the design went out of its way to keep first-class. It also silently re-routes signal that
should be stutter into `ε`/outlier mass.

**Concrete fix (pick one, state it):**
1. **Restore the sum for impure alleles.** For an impure candidate, define `align(obs |
   cand⊕Δ)` as `Σ_placement (1/n_placement) · align(obs | cand⊕Δ@placement)` over the
   placement-distinct sequences, with HipSTR's uniform position prior. This *is* the
   faithful borrow; it costs a placement loop (bounded by the tract's unit count, cheap, and
   HipSTR's `upstream_match_lengths_` trick collapses equal-LL runs). The "single align
   term" simplification then holds only as the pure-allele special case.
2. **Drop the marginalization claim for impure alleles.** Keep the committed
   interruption-preserving placement, and rewrite §3 to say placement is marginalized
   **only for pure alleles** (where it is degenerate); for impure alleles a *fixed*
   interruption-preserving placement is used, an explicit modeling approximation with the
   bias noted above. Then "the data picks it" must be struck for impure alleles.

The docs currently *claim* (1) while `⊕`'s definition implements toward (2). That is the
break.

---

## 1 — CAVEAT: the integer/fixed-point reduces hold; the M4-introduced float reduces are asserted-deterministic, not established

> **RESOLVED 2026-06-22 (verify-fix #1 — extend the fixed-point/tie-break discipline to the
> decision floats).** The three M4-introduced quantities are now *specified* order-independent,
> not assumed: (i) **`ℓ_pen`** (the burn-in plateau stop) is summed by the **same fixed-point
> integer accumulation** (`scale → round → i128`) as `F_i`/level, so `Δℓ_pen/|ℓ_pen| < tol`
> compares byte-identical values and the **stop iteration** cannot drift across thread counts
> (no `i128` overflow — per-locus normalizers are `reads × O(few nats)`; precision `2⁻⁴⁰` nats
> ≪ any `tol`); (ii) the **BIC comparisons** (ε-freeze, group split/merge, group count) are
> differences of two such fixed-point `ℓ_pen` values against a *deterministic* penalty → a
> deterministic discrete outcome; (iii) the **clustering** keeps its no-random-init plus now
> **deterministic tie-breaks** — strict `<` thresholds with a defined equal-handling rule and a
> total tie-break on the **sample catalog index** — so a borderline `distance ≈ threshold`
> resolves identically. Result: the pre-pass is byte-identical across thread counts at its
> *decision layer* (stop iteration, freeze/merge, group count, group membership), not only its
> sufficient statistics. Applied to spec §4.4 (new verify-fix #1 amendment + burn-in stop +
> cluster + schedule), parameters §3/§4/§7.

**Verdict: caveat.** Two layers, opposite outcomes.

**What survives refutation.**
- The sufficient statistics are genuinely order-free: `SlipProfile { down/up: [u64; …] }`
  and `SampleStutterStats { base_match/mismatch: u64, by_length: …u64 }`
  ([parameters §3](../../architecture/ssr_call_parameters.md)) are `u64`; integer addition
  is associative + commutative → bit-identical regardless of thread count. **Holds.**
- The M1 fix for the `F_i` and per-group **level** reduces is sound. Scale `×2⁴⁰` → round →
  sum into per-individual / per-group `i128`
  ([genotyping §5](../../architecture/ssr_call_genotyping.md) lines 252-262). **No overflow
  risk:** each addend ≤ `2⁴⁰ ≈ 1.1e12`; even `10¹²` loci × `2⁴⁰ ≈ 1.1e24` ≪ `i128::MAX ≈
  1.7e38`. The rounding is a *fixed deterministic* operation (same on every thread), so the
  sum is byte-identical across thread counts — which is the whole point; it is not a
  reproducibility hazard. The final divide is a deterministic scalar. **Holds.**

**What does not survive — the new surfaces M3/M4 opened.** The M1 resolution covers the two
named float reduces (`F_i`, level). It does **not** cover the *other* float quantities the
M4 and M3 amendments introduced, and the docs assert their determinism without establishing
it:

1. **The `ℓ_pen` plateau stop.** [parameters §4](../../architecture/ssr_call_parameters.md)
   defines the burn-in stop as `Δℓ_pen/|ℓ_pen| < tol`, where `ℓ_pen` is "the E-step's own
   normalizer (`log Σ`)" summed over the subset. That is a **float sum of per-locus
   log-normalizers over loci that complete out of order** — the identical non-associativity
   the `F_i` reduce has. The doc *asserts* "the reduce is order-independent → the trajectory
   and the `ℓ_pen` plateau are reproducible across thread counts," but `ℓ_pen` is **not**
   reconstructable from the integer sufficient statistics (those drive the M-step, not the
   marginal likelihood), and nothing states it uses the fixed-point trick. A borderline
   `Δℓ_pen/|ℓ_pen| ≈ tol` can therefore **stop one iteration earlier/later** depending on
   thread count, changing the frozen `ε`/seed and hence downstream calls.
2. **The BIC freeze/merge/group-count comparisons.** ε-freeze, group merge, and the number
   of sample groups are all `ℓ_pen(rich) − ℓ_pen(simple)` vs a BIC penalty
   ([parameters §4 / Q-P2 / Q-P6](../../architecture/ssr_call_parameters.md)). Same float
   `ℓ_pen`; a comparison sitting within last-bits of the penalty can **flip a discrete
   structural decision** (freeze vs covariate; one group vs two) across thread counts — a far
   larger swing than a borderline genotype.
3. **The distance-based clustering threshold.** Sample grouping is "distance-based grouping
   of close neighbours in (ε, level) space" with uncertainty-scaled distances
   ([parameters §4](../../architecture/ssr_call_parameters.md)). A sample whose scaled
   distance ≈ the merge threshold can land in different groups if the inputs (`ε`, level)
   carry thread-dependent last bits — changing per-group chemistry for a whole cluster.

**Why it matters.** Each of these feeds the *frozen* parameters and group assignments that
all of Phase 3 consumes; a flip is not a one-bit genotype wobble but a structural change
(different group count, different freeze decision, different stop point) propagating
cohort-wide. The byte-identity-across-threads invariant the design repeatedly claims is
**not yet earned** for the pre-pass's decision layer.

**Concrete fix.** Specify `ℓ_pen` (and the per-model log-likelihoods entering the BIC
comparisons) as a **fixed-order** reduce keyed by locus catalog sequence, or accumulate them
with the same `i128` fixed-point scheme as `F_i`/level (one scalar per locus — cheap). Make
the clustering deterministic under ties: strict thresholds with a defined equal-handling
rule and a total tie-break on sample index. Until then, downgrade "byte-identical across
thread counts" for the pre-pass to "byte-identical *given fixed-order `ℓ_pen` / BIC / cluster
reductions*, TBD" — exactly the downgrade M1 forced on `F_i`, now owed again for the M4
floats.

---

## 7b — CAVEAT: the exact-AF convolution is reusable, but not as-is, and the m4 sum is an approximation

> **RESOLVED 2026-06-22 (verify-fix #7b — name the generalization, scope the approximation).**
> Both specifics are now stated rather than glossed: (a) the m4 QUAL is a **generalization of
> `compute_qual_via_exact_af`, not a call to it unchanged** — the docs now say SSR must
> parametrize the allele-0 collapse (`count[i]`, not `count[0]`), set the per-`i` Beta from
> `G₀` (`α_ref = G₀[i]`, `α_alt = Σ_{j≠i} G₀[j]`), and bypass the REF-anchored `P(K=0)`
> wrapper; the kernel (`convolve_ac_linear`) + Beta-Binomial/`K=0` math are reusable, the
> allele-0 indexing is not (calling it as-is yields only `P(fixed-for-the-modal/REF allele)`).
> (b) "properly normalized / exact" is **relabelled to its honest scope** (option (i)): the
> sum `Σ_i P(fixed-for-i)` is **exact for biallelic** and a **per-collapse-normalizer
> approximation for ≥3 alleles** (each term carries its own binary partition `Z_i`; the
> ref/non-ref Beta lumping is the engine's own multiallelic approximation) — accurate at the
> extremes, loose only mid-range, fine for a QUAL score. The **exact form is named as the
> upgrade**: extract each term's unnormalized `K=0` numerator (`log_p_ac_next[0]` pre-divide)
> and divide the summed numerators by one joint multiallelic `Z`. Applied to spec §4.5 + §9
> settled list, genotyping §5 (Q-G2) / §6 / §9 (Q-G4).

**Verdict: caveat.** The m4 QUAL estimator is feasible and largely as the doc describes, but
two specifics undercut "thin reuse" and "exact, properly normalized."

**Where allele-0 is hard-wired.** `compute_qual_via_exact_af`
([posterior_engine.rs:3010-3095](../../../src/var_calling/posterior_engine.rs)) anchors on
REF = allele-0 in exactly **one** structural place: Step 1, line 3036,
`let c = (p as u32 - counts[0]) as usize;` — every sample's collapsed likelihood is bucketed
by *non-allele-0* copy count. Everything else is already either parametrized or
allele-agnostic:
- The convolution kernel `convolve_ac_linear`
  ([:3132-3207](../../../src/var_calling/posterior_engine.rs)) takes pre-collapsed
  per-sample taps — **fully reusable**, no allele identity inside.
- The Beta-Binomial prior (Step 3) is parametrized by `alpha_ref`, `alpha_alt` arguments
  (lines 3018-3019) — the caller already supplies them (lines 2327-2328:
  `alpha_ref = pseudocounts[0]`, `alpha_alt = Σ pseudocounts[1..]`).
- The `P(K=0)` readout (Step 4, line 3089) means "all samples have zero non-collapsed
  copies" = fixed-for-the-collapsed-allele — **correct** for fixed-for-i once Step 1 collapses
  on i.

So the m4 reuse needs: (a) generalize Step 1's `counts[0]` → `counts[i]` (a signature/index
parameter, not a rewrite); (b) per term set `alpha_ref = G₀[i]`, `alpha_alt = Σ_{j≠i} G₀[j]`;
(c) loop `i ∈ A_ℓ`, sum. That matches what [genotyping §6 / Q-G2](../../architecture/ssr_call_genotyping.md)
says ("reuse the convolution kernel `|A_ℓ|`×, bypass the REF-anchored P(K=0) wrapper"), so
the doc is *not* contradicted — but "thin reuse" undersells the Step-1 generalization, and an
implementer who calls `compute_qual_via_exact_af` unchanged gets P(fixed-for-the-modal-allele)
only (the existing single REF collapse), **not** the sum. State the Step-1 index parameter
explicitly so it is not missed.

**The "exact, properly normalized" claim overstates.** Each `P(fixed-for-i)` is the posterior
`P(K_i = 0 | data)` computed under its **own** binary-collapse model `M_i` with its **own**
normalizer `Z_{M_i}` (Step 4's `log_z`, line 3081). The terms are summed across `i`, but each
is divided by a *different* `Z_{M_i}`, so `Σ_i P_{M_i}(K_i=0 | data)` is not a coherent
posterior under any single model. Two distinct approximations are folded in:
- For `|A_ℓ| ≥ 3`, the binary collapse lumps all non-i alleles into one Beta category. The
  *prior* aggregation is exact (Dirichlet → Beta), but Step 1 sums genotype likelihoods within
  a bucket **unweighted**, losing the relative π of the lumped alleles. (This approximation
  already exists in the engine's own multi-allelic QUAL — biallelic is exact — so it is
  inherited, not new.)
- The per-i normalizers `Z_{M_i}` differ, so the cross-term sum is only approximately the
  joint `P(monomorphic | data)`. It is accurate at the extremes (clearly fixed → dominant
  term ≈1, rest ≈0; clearly variable → all ≈0) and loosest in the borderline middle.

**Why it (mildly) matters.** For a QUAL *score* the approximation is defensible, but the
docs' "exact (no inclusion–exclusion)... properly normalized" reads as a stronger guarantee
than the math delivers. "No inclusion–exclusion" is correct (the events are mutually
exclusive); "exact" is not, for `≥3` alleles or borderline loci.

**Concrete fix.** Either (i) keep the sum but relabel it "approximate `P(monomorphic)`
(exact for biallelic; per-collapse normalizers for `≥3` alleles)"; or (ii) if exactness is
wanted, extract the **unnormalized** `K=0` numerator per i (`scratch.log_p_ac_next[0]` before
the `log_z` divide, line 3089) and divide the summed numerators by a **single** joint
normalizer — at the cost of computing that joint `Z` (a genuinely multi-allelic AF
marginal, a bigger change than the kernel reuse). For a QUAL field (i) is the honest, cheap
choice.

---

## 4 — CAVEAT: the `F` ceiling holds; the geometric `G₀` pseudocount has no hard floor against f64 underflow

> **RESOLVED 2026-06-22 (verify-fix #4 — give `G₀` the hard floor the `F` ceiling already
> has).** The `F` ceiling needed no change (it is a hard clamp, verified airtight). The `G₀`
> half is fixed: the geometric pseudocount is now **floored explicitly** —
> `G₀[i] = max(p^|Δ|, FLOOR)` for a tiny positive `FLOOR` (any representable `> 0` re-floors
> it above hard zero; the Dirichlet M-step grows `π_i` from real evidence regardless of the
> pseudocount's size, so the floor only has to keep a *zero-evidence* far candidate alive),
> or equivalently `G₀` carried in the **log domain** into the M-step (the cleaner principled
> form). The docs now also state the `|Δ|·(−ln p) ≳ 744` underflow bound and that
> short reads are safe while long-read / steep-decay cohorts cross it. Framed as the `π`
> analogue of `F_CEILING = 0.99` — both keep a multiplicatively-absorbing parameter off its
> trap value. Applied to spec §4.3, parameters §5, genotyping §5.

**Verdict: caveat.** The two guards are not equally robust.

**The `F` ceiling holds — solidly.** The clamp order is raw `F_i` → shrink to cohort mean →
clamp ≤ user cap → clamp ≤ `F_CEILING = 0.99`
([genotyping §5](../../architecture/ssr_call_genotyping.md) line 245; [spec §4.4](../../specs/ssr_cohort_mark2.md)
lines 806-810). No path reaches exactly `1.0`: the raw `F_i` is a mean of autozygous
responsibilities `F·π_i/(F·π_i+(1−F)·π_i^ploidy)`, each `<1` because the *current* `F ≤ 0.99`
and `π_i < 1`; the cohort-mean shrink is a convex combination of values `<1`; and the final
`0.99` clamp is unconditional and last. The engine consumes it as `log_one_minus_f =
safe_ln(1 − f_s)` ([posterior_engine.rs:2170](../../../src/var_calling/posterior_engine.rs));
`1 − 0.99 = 0.01` → finite log. The absorbing trap is genuinely closed.

**The `G₀` floor is not hard.** The pseudocount decays geometrically as `p^|Δ|` in the unit
offset ([parameters §0/§5](../../architecture/ssr_call_parameters.md); [spec §4.3](../../specs/ssr_cohort_mark2.md)
lines 403-412). There is **no explicit `max(p^|Δ|, floor)`** anywhere in the spec — the
"floor" language refers to the pseudocount keeping `π > 0` *in principle*, not to a numeric
clamp. In f64, `p^|Δ|` underflows to **exactly 0.0** once `|Δ|·(−ln p) ≳ 744`:
- Short reads bound it safely. Candidates are "hard-capped at the read length"
  ([spec line 970](../../specs/ssr_cohort_mark2.md)), so for 150 bp reads on a di-repeat
  `|Δ| ≤ ~75`; even a steep `p = 10⁻³` gives `75 × 6.9 = 518` → `e⁻⁵¹⁸ ≈ 10⁻²²⁵`, well above
  the `~10⁻³⁰⁸` f64 floor. Safe.
- **Long reads break it.** A HiFi read spanning a 1 kb di-repeat admits `|Δ| ~ 500`; with
  `p = 10⁻³`, `500 × 6.9 = 3450` → `e⁻³⁴⁵⁰ = 0.0` exactly. The candidate's pseudocount is
  then **0**, re-creating the `π = 0` absorbing trap *numerically* for exactly the
  far-from-mode allele the floor exists to keep recoverable — and because candidates are
  empirical, a far candidate exists only because a read supports it, i.e. precisely when we
  must not trap it. The tool explicitly targets non-model genomes where long-read input is
  plausible.

The review's §1 credited this guard as "sound and correctly mirror[ing]" the `F` ceiling.
The mirror is imperfect: the `F` ceiling is a clamp, the `G₀` floor is an unclamped
exponential.

**Concrete fix.** Floor the pseudocount explicitly — `G₀[i] = max(p^|Δ|, FLOOR)` with a tiny
positive `FLOOR` (e.g. `1e-300`), or, cleaner, carry `G₀` in the **log domain** and feed the
engine `log`-pseudocounts (the Dirichlet M-step at
[posterior_engine.rs:2827-2829](../../../src/var_calling/posterior_engine.rs) sums
`expected_counts[a] + pseudocounts[a]`, so a log-domain swap touches the M-step too — but it
is the principled fix). At minimum, state the `|Δ|·(−ln p)` underflow bound and that
long-read / steep-decay cohorts must clamp.

---

## 5 — HOLDS (headline); residual within-sample exposure is real, acknowledged, and unspecified

> **CLARIFIED 2026-06-22 (verify-fix #5 — dedup is an upstream input contract, not a
> caller concern).** The "unspecified dedup" half is resolved by stating the contract, not by
> adding machinery: PCR/optical duplicates are marked **upstream** (post-mapping `markdup`,
> BAM flag `0x400`), and **Stage 1 skips flagged reads** when it builds the `(seq, count)`
> tallies — so the depth `ssr-call` sees is **duplicate-free by construction**, and there is
> no "assumed deduplicated" unchecked precondition and no in-caller dedup to specify. Any
> residual PCR jackpotting is whatever the upstream marker missed (imperfect without UMIs) —
> an **input-quality matter, explicitly out of scope**. The cohort-recurrence headline is
> unchanged. Applied to spec §2 (input contract) + §4.4 + §9, parameters §3/§4/§7. (The
> overdispersion-discounted-effective-depth idea is therefore **not** adopted — it solved a
> problem the input contract already removes.)

**Verdict: holds, with the residual flagged.** Rung recurrence counts **distinct samples**,
not reads: a length is admitted by recurrence when "seen in ≥ *k* samples"
([genotyping §2](../../architecture/ssr_call_genotyping.md)), so within-sample PCR
duplication cannot manufacture a cohort-recurrent false allele. The headline claim survives
refutation.

**The residual the review punted to m2(b) does not yet have an answer.** Within a single
sample, a jackpot is not gated by sample-level recurrence and still touches three quantities:
- the **clear-maximum prominence** (`> 3 reads above each neighbor`,
  [genotyping §2](../../architecture/ssr_call_genotyping.md)) — a jackpotted minor length can
  clear the floor;
- the **allele-balance / overdispersion** term — designed to catch *systematic concentrated*
  minor signal, so it is the right defense, but it runs on deconvolved responsibilities and
  its overdispersion form is "the documented upgrade," **not v1**
  ([genotyping §4](../../architecture/ssr_call_genotyping.md));
- the **clustering precision** — defined as effective *deduplicated* depth
  ([parameters §3](../../architecture/ssr_call_parameters.md), `effective_depth`), but
  **nothing enforces dedup** and the "overdispersion-discounted effective depth" is listed as
  *still open* (m2(b)). A non-deduplicated jackpotted sample overstates its precision and can
  dominate the distance-based clustering it should be *uncertain* in.

So the answer to "is the overdispersion-discounted effective depth sufficient, and is it
specified anywhere yet?" is: **not specified** — it is a named upgrade in two places and
dedup is assumed, not checked. The claim holds where it was made (cohort recurrence); the
residual is correctly scoped but remains an open obligation, not a closed one.

---

## 2 — HOLDS: per-locus EM is a pure deterministic function of its (now-expanded) inputs

**Verdict: holds**, conditional on #1. The per-locus EM is one locus = one worker
([genotyping §7](../../architecture/ssr_call_genotyping.md)); its float reduces are over that
locus's own samples/genotypes in a fixed array layout, so the summation order is fixed by the
data, not by thread scheduling — the iteration count (non-decreasing penalized log-lik stop)
is therefore deterministic and identical on whichever worker runs it. The reused engine does
no cross-record float accumulation (verified: every record's EM is independent;
[posterior_engine.rs](../../../src/var_calling/posterior_engine.rs) `run_em_for_record` is
self-contained). The inputs now include the per-`(group, period)` shape and the
outer-loop-refined per-group level — but within a sweep these are **fixed** (updated only at
the barrier), so they enter the locus EM as constants and do not break purity.

**The one linkage to state.** "Pure function of its inputs" is only as deterministic as those
inputs. The `F_i` vector, per-group level, and group assignments are produced by the
cross-locus reductions of #1 — so the per-locus EM is byte-identical across threads **iff**
those reductions are (which, per #1, is established for `F_i`/level but not yet for the
clustering that fixes group membership). The claim itself holds; its determinism inherits
#1's caveat.

---

## 6 — HOLDS: the engine's per-sample `F` input matches the SSR `F_i` update

**Verdict: holds.** The engine consumes a per-sample fixation index vector
(`fixation_index_overrides: Option<Vec<f64>>`,
[posterior_engine.rs:237](../../../src/var_calling/posterior_engine.rs)) →
`scratch.log_f_per_sample[s] = safe_ln(f_s)` and `log_one_minus_f_per_sample[s] =
safe_ln(1 − f_s)` ([:2167-2171](../../../src/var_calling/posterior_engine.rs)), feeding the
IBD-mixture prior `F·π_i + (1−F)·π_i^ploidy` (used at
[:2569, :2770](../../../src/var_calling/posterior_engine.rs)). Domain requirement: `f_s ∈
[0,1)` — `f_s = 1` would make `log(1−f) = −∞`. The SSR `F_i` update is the mean
autozygous-branch responsibility ∈ `[0,1]`, then shrunk and clamped ≤ `0.99`
([genotyping §5](../../architecture/ssr_call_genotyping.md)), so the fed-back value is always
`< 1` → prior stays valid. The only wiring note: `resolve_fixation_indices` requires the
override vector length to equal the record's sample count
([:616-619](../../../src/var_calling/posterior_engine.rs)), and `CohortLocus` is sparse
(present samples only), so the per-locus engine call must pass the **present-subset** `F_i`
slice indexed by `present[k]` — which is exactly what [reading Q-R1](../../architecture/ssr_call_reading.md)
specifies. No prior-side engine change is needed; the narrow reuse claim is correct (it would
over-claim only if it implied the engine *estimates* `F`, which it does not, and the docs do
not).

---

## 7a — HOLDS: the pseudocount swap and class-path bypass are concrete and minimal

**Verdict: holds**, now pinned to exact lines. The per-allele pseudocounts are materialized
once per record as an opaque `Vec<f64>` —
`scratch.pseudocounts[a] = pseudocount_for(class, config)`
([posterior_engine.rs:2195](../../../src/var_calling/posterior_engine.rs)) — and consumed by
the M-step as a plain Dirichlet vector,
`p_hat_next[a] = (expected_counts[a] + pseudocounts[a]) / denom`
([:2827-2829](../../../src/var_calling/posterior_engine.rs)), with **no further class
dependency** (the biallelic fast path at
[:2844-2861](../../../src/var_calling/posterior_engine.rs) is allele-count-based, not
class-based). So the SSR swap is the single fill line 2195 → the per-candidate `G₀` vector.

The three paths that key off `AlleleClass`, and how each is handled:
1. **`pseudocount_for(class)`** ([:1737-1744](../../../src/var_calling/posterior_engine.rs))
   → `scratch.pseudocounts` → the Dirichlet M-step **and** the QUAL alphas (`alpha_ref =
   pseudocounts[0]`, `alpha_alt = Σ pseudocounts[1..]`,
   [:2327-2328](../../../src/var_calling/posterior_engine.rs)). **This is the swap point** —
   and note it ties to 7b: swapping in `G₀` makes the *existing* QUAL call compute
   P(fixed-for-the-modal-allele), which is exactly why m4 needs the generalized per-i loop.
2. **`map_to_contam_class(class)`** ([:1358-1365](../../../src/var_calling/posterior_engine.rs))
   → `mixture_contam_class_per_allele`, consumed only inside
   `compute_mixture_log_likelihoods` when the contamination side-pass is active. SSR runs no
   contamination → `c_s = 0` → the mixture is skipped
   ([:1433-1444](../../../src/var_calling/posterior_engine.rs)) → **bypassed**.
3. **`is_compound` / `compound_mask` / `f_hat_compound` / `compound_pseudocount`**
   ([:2196, :2294, :2532, :2890-2898](../../../src/var_calling/posterior_engine.rs)) → the
   compound-allele M-step branch. SSR has no compound alleles → `compound_mask` all false →
   `p[a]` reads `p_hat[a]` not `f_hat_compound[a]` → branch **inert**.

So "bypass contamination/compound" is concrete: contamination is gated off by `c_s = 0`,
compound by an all-false `compound_mask`. The pseudocount-vector generalization is genuinely
low-invasiveness on the math path, as claimed.

---

## Bottom line

Of the seven, **one breaks** (#3, for impure alleles — the placement marginalization is not
reproduced by a single-sequence substitution-only `align`, and `⊕`'s
interruption-preserving definition contradicts the "data picks it" claim), and **three newly
caveat** because the post-review amendments moved the ground:

- **#1** — the M1 fixed-point fix secured `F_i`/level, but M3/M4 introduced new float reduces
  (`ℓ_pen` plateau, BIC freeze/merge/group-count, clustering distances) whose
  order-independence was asserted, not established; a borderline tie could flip a *structural*
  pre-pass decision across thread counts. **Resolved 2026-06-22 (verify-fix #1):** `ℓ_pen` and
  the BIC log-likelihoods now use the same fixed-point `i128` accumulation as `F_i`/level, and
  the clustering takes deterministic tie-breaks on the sample catalog index — so the decision
  layer is byte-identical across thread counts.
- **#7b** — feasible and consistent with the doc's stated bypass, but `compute_qual_via_exact_af`
  hard-codes allele-0 in Step 1 (needs a generalized collapse index), and the sum-over-collapses
  is an approximation, so "exact, properly normalized" overstated. **Resolved 2026-06-22
  (verify-fix #7b):** the docs now name the allele-0 generalization explicitly and relabel the
  estimator "exact for biallelic / per-collapse-normalizer approximation for ≥3 alleles," with
  the single-joint-`Z` exact form named as the upgrade.
- **#4** — the `F` ceiling is airtight, but the geometric `G₀` pseudocount had no hard floor
  and underflowed to exactly 0 at large `|Δ|` (long-read / steep-decay), numerically
  re-creating the very trap it guards against. **Resolved 2026-06-22 (verify-fix #4):** `G₀`
  now floored explicitly (`max(p^|Δ|, FLOOR)`) / carried in the log domain — the `π` analogue
  of the `F` ceiling.

The remaining three (**#2, #6, #7a**) survive my refutation as stated, with #2's determinism
explicitly inheriting #1's caveat. The genuinely-HipSTR borrow is faithful **only for pure
alleles**; everything the design newly *froze* or *bolted on* (the chemistry-group axis, the
likelihood-based pre-pass criteria, the m4 QUAL reuse) is where the load-bearing risk now
concentrates — which is the same conclusion the review reached about Milestone D, now extended
into the determinism and impure-allele machinery the "couldn't break" list had cleared.
