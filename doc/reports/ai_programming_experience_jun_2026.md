# The never-ending toy

Alternative title: Storming the heavens of Scientific software

*An evolving experience report: building scientific software with an AI coding agent.*

---

## The problem we were stuck with

We needed to genotype **2,500 tomato samples** from whole-genome sequencing data.

The hardware: **32 cores, 128 GB RAM**. Mapping was the easy part —
`fastp` + `minimap2`, a couple of months and done. Calling variants
across a cohort that large was where we got stuck:

- **freebayes** doesn't scale to this many samples, in time and in memory.
- **GATK**: per-sample calling was doable — a couple of months. But the
  population-wise gVCF merge:
  
  - would take **about a year**,
  - it **doesn't really merge indels**, and
  - the genome has to be cut into segments because the merge step has a **memory leak** and can't run for long without running out of memory.

So the standard recipe — call each sample, then merge — was blocked.

## The plan we tried

We tried writing a population-wise gVCF merger in **Python** (the language we knew).
The plan had an optional later phase: trying to translate it to **Rust** for speed.

We built the Python merger and we translated it to Rust. And then we hit a wall that had nothing to do
with speed:

**You cannot compute the joint genotype posterior for indels from the data
stored in per-sample gVCFs.** By the time you're merging outputs, the
information you'd need is already gone.

The only way across would be to stop post-processing and compute the statistics from the alignments yourself — a real caller, from BAM.

## The doubt

Which left us staring at a much bigger task than we'd signed up for: write
a multi-sample variant caller from scratch, in Rust, with the statistics
done properly.

The problem:

- Although I can *read* Rust, my writing is very slow.
- I don't have the statistics. I have a feel for Bayesian inference, not
  the ability to derive or code the models.

A biologist who could neither write the language nor code the math, about
to attempt a tool that — historically — took a team. **Should we even try?**

We had managed to write the gVCF merger in about six weeks (serious work
started mid-February 2026; a working Rust merger with EM-derived genotype
posteriors was running by 25 March), so we decided to try this new
direction for a couple of weeks, then evaluate where we were and decide
what to do.

After those two weeks we were impressed by what was already done. In the
first twelve days after the pivot (17–29 April) we had written the
authoritative **six-stage architecture spec** and a set of design-principle
specs, and then implemented *and* code-reviewed the first stage of the new
pipeline — the CRAM input reader (header validation, multi-file
peek-and-scan merge, filter cascade). Within three weeks (6 May) the
pileup walker was working and being validated against samtools. That was
enough to convince us to keep going.

---

## First contact: it's not the chat window

My first instinct about "AI coding" was the chat interface: describe a
problem, get a snippet, paste it, fix it by hand. If that's all you've
seen, you don't yet know what these tools can do.

The thing that changed everything was the **agent loop**: the AI writes
the code, compiles it, runs the tests, *reads its own compiler and test
errors*, and fixes them — without me in the inner loop. I stopped typing
code and started steering the coding: I became a project manager, and my
colleagues were the coding agents.

With the new approach we managed to build a six-stage cohort caller (BAM → pileup →
`.psp` → variant grouping → per-group merge → posterior engine → VCF).
And once that worked, a *second* caller from scratch, for microsatellites
(SSR/STR). Full history: `doc/reports/project_history.md`.

The milestones, dated from the git log (elapsed measured from the pivot,
17 April 2026):

| Date | From pivot | Milestone | Size |
|---|---|---|---|
| 2026-02-28 | — | Fast gVCF parser working | ~1.7k LOC, 10 files |
| 2026-03-25 | — | gVCF merger working (EM posteriors) | ~6k LOC |
| **2026-04-17** | **day 0** | **Pivot: decide to build a caller from BAM** | — |
| 2026-04-24 | +7 days | Six-stage architecture spec lands | — |
| 2026-04-29 | +12 days | First pipeline code (CRAM input reader) | ~9.4k LOC |
| 2026-05-06 | +19 days | Pileup walker working | — |
| 2026-05-13 | +26 days | `.psp` binary format (writer + reader) | — |
| **2026-05-19** | **+32 days** | **Six-stage caller runs end-to-end** (synthetic) | ~45k LOC |
| 2026-05-24 | +37 days | Caller running on real tomato BAMs | ~59.8k LOC, 96 files |
| 2026-06-04 | +48 days | Pipeline re-architected, byte-identical swap | ~70k LOC |
| 2026-06-23 | +67 days | Second caller (SSR/STR) emits a VCF end-to-end | ~80k LOC |

As of today: **~1,073 commits, 136 Rust source files, ~83k lines of
code** — one person, with the assistant model itself moving Opus 4.6 → 4.7
→ 4.8 over the six months.

In the process I have changed a lot about how I use the AI tool to create the scientific tooling that I need.

## How do I use AI for programming

The naive version is: "tell the AI what you want and let it cook".
This is what people call vibe coding, and although the results could be useful,
they won't be the best. What I backed into, feature by feature, was a much
more disciplined loop.

Before coding:

1. **Specification**: what we want to build.
2. **Architecture**: how are we going to build it, how the data is going to flow inside the application and how are we going to structure the code.
3. **Implementation plan**: the steps and milestones to code it.

### The spec is key

Before committing to the full caller, I spent **nearly two weeks writing
specifications, with zero pipeline code** — from the first design notes on
17 April to the first line of the CRAM input reader on 29 April. That
architecture document is still cited by name in commit messages *months*
later.

How a spec gets made, in practice:

- First draft by the AI.
- Point-by-point argument until I understand the questions at hand and adapt them to my needs.
- **Adversarial review** — a literature search *and* reading the source of
  related tools (we keep freebayes, GATK, bcftools, htslib checked out so
  the AI reads them before porting an algorithm; for the SSR caller it
  digested the HipSTR / GangSTR / ConSTRain papers and source into a
  research report that *fed* the design).
- Point-by-point discussion of everything that review turns up.

AI is very useful at that adversarial pass. One spec for the
microsatellite cohort caller was hardened against a **13-finding
adversarial review** — thirteen real problems in *my own design*, found
and fixed before a line of code existed.

### Coding

1. **Coding**, following the plan.
2. Mechanical **code review**, then fixing what it finds.
3. **Performance analysis**, then tuning.

All these actions are governed by skills.

Across ~1,073 commits, raw implementation
is only about **30–34%** of them. Review (correctness fixes), performance,
and docs/specs together **equal or exceed** new feature work; and in the
later months, docs alone nearly doubled to ~18%. We spend as much effort
reviewing and specifying as writing new code.

The coding follows some recommendations:

- **TDD, red/green.** Tests first, then the code that turns them green.
- **Strict with the linter**, lint suggestions are near-mandatory.
- Rust's static, strong typing turned out to be a force multiplier *for
  the AI*: the compiler hands it a tight, unambiguous signal on every
  iteration. In Python we'd have had to add type annotations and a type
  checker just to recover that guardrail.
- In Rust the compiler kills whole classes of error before
  anything runs, by taking advantage of the strong types.
- **Byte-identity oracles for rewrites.** When we re-architected the entire
  pipeline — **+21,700 / −31,100 lines in a single day** — it was allowed
  into production only because an oracle proved the new pipeline emits the
  *exact same VCF* as the old one. Trust came from a proof, not from me
  squinting at a diff.


### Skills: codify the best practices

We wrote reusable **skills** for code review, for applying the
fixes, for performance review, for feature implementation.

Skills:

- are not optional.
- Don't hand-write them: have the AI draft them via deep research, then
  revise.
- They're living documents: update them whenever the process breaks.
- There are skills to write and revise skills.

## QA is paramount and the limiting factor

If the AI is writing most
of the code and I'm not reading every line: *how do I trust any of it?*

"I don't read every line" is worrying, but what I eventually realized is that even when I did write the code I had to do QA, and that was the critical bit: you should never trust what you write. So what's new is *who* writes the code — not what you're allowed to trust, which is still nothing until QA says otherwise.

Our **trust stack**:

- **Validation against reference tools**, samtools, GATK, freebayes, on
  real data.
- **Simulators with known ground truth**, so we can score against an answer
  we control.

At the end of the day, we should never trust the code we write; we need other QA measures.

## What I could suddenly afford: experiments

In software we usually code what we think will be the best algorithm for the task at hand, but on many occasions there are other algorithms.
Ideally we should implement all those algorithms, test them with real data, and choose what to do only after having empirical evidence; but if coding one algorithm is hard enough, coding several is often too time-consuming, so we just pick the option that sounds best, or the one a respected paper used, because
building all the candidates to actually compare them is too expensive. AI
coding changes that arithmetic. When fully implementing a candidate drops
from weeks to an afternoon, building them all and deciding only after seeing the data becomes a possibility.

Just an example case: `doc/reports/workshop_ai_fast_coding_idea_testing.md`

For the microsatellite genotyper we needed an HMM alignment algorithm. We had three credible ways to compute it: the literature standard (HipSTR), an incremental fix to what we already had, and our own in-house first implementation. We built all three
behind one interface and raced them on simulated data with known answers.

- Our own elegant idea collapsed: it just created bad alignments in many cases, and no whiteboard argument would have settled it as cleanly as one run did.
- The literature choice won; it ran ~100× faster than our first implementation.

Cost: ~1,400 lines committed, of which ~310 were written, tested, and
deliberately deleted. That "wasted" 310 lines is the
entire point: in a slower world, that cost is exactly what stops you from
testing the idea at all, so you follow your intuition or argue instead.

The same freedom shows up as a willingness to throw work
away: the project pivoted off its own name twice and deleted a whole
first-generation module once a better one existed. *Cheap to build makes it
cheap to discard*.

## AI is not good at everything?

This is just my current and evolving opinion about a series of very rapidly evolving tools.

AI is **very good at:**

- Writing focused code to a clear spec, architecture doc, and implementation plan — especially with the review loop on top.
- Writing code I *can't* write myself. The Bayesian and HMM models are the
  clearest case: they are not just faster to write, but things that were out of my
  reach. For a domain scientist, this is the genuinely new power.
- Understanding large codebases, ours and other people's.
- Finding bugs, much better than I do.

AI is quite **bad at:**

- **High-level architecture and data flow.** It produces *an*
  architecture, rarely a *good* one.
- **Naming and readability.** Its vocabulary tends to be both jargony and generic. Left alone, it writes code that works but reads badly.

## Where it actually got us

- **The pileup, validated.** About **three weeks** after deciding to build
  a caller (pivot 17 April → samtools comparison 7 May), Stage-1 pileup
  matched samtools to **~95%** on SNPs and alleles on real data — our
  walker slightly stricter, most likely because of BAQ.
- **The caller, compared.** The end-to-end six-stage caller ran on
  synthetic fixtures on **19 May** — about **two weeks after the pileup
  worked** — and first runs against GATK on real tomato data followed the
  next day (20 May). Head-to-head tuning against GATK / freebayes /
  samtools has continued through June (e.g. F1 vs GATK 0.285 → 0.317 after
  a MAPQ filter). This is real-data *tuning*, not yet a published or
  production bar.
- **A whole second caller.** The SSR/STR genotyper went from the first
  design docs (9 June) to an **end-to-end VCF in about two weeks** (23
  June). It is **synthetic/simulator-validated only** — real-data
  validation is the open item, exactly where the SNP caller stood a month
  earlier.

---

## Lessons I think I've learned

### Coding is not the bottleneck

Coding is not the slow step. The slow
steps are now **specification and the discussions around it**, and
**validation**.
The limiting factor has moved from "can I build the tool?" to "can I prove the tool is right?".

### We can create more

Building the original GATK took statisticians, Java engineers, and geneticists, ~9–12 core people
over many months, which means real money, a managed project, and funding
to match. One person plus an AI did a comparable *scope* in a couple of
months. Put to numbers
(`doc/reports/cost_estimate_first_working_version.md`): the current
two-caller codebase would conventionally be a **~2–3 person-year,
€0.9M–€3.5M** project; actual cost so far is roughly **€11k–€14k** including
AI usage — a **~75–400× compression**, large enough to survive every
conservative adjustment.

**The speed-up isn't uniform.** Parsers, plumbing, scaffolding, plan
  docs: huge. Genuine design, numerical-accuracy debugging, surprises on
  real data: much less.

### More quality or more slop?

My prediction is that the scientific community will produce both **much better specialized tools** *and* much **more slop and cargo-cult science**.

The *same* tools that let you build a caller in a
month let you ship confident-looking garbage even faster.
Under publish-or-perish that's a real hazard: if you don't understand what
you're doing, you can now produce much more cargo-cult science, which will now be renamed slop science.
I prefer the former name because this is not a new phenomenon; we just have better tools to produce more bad science.

When anyone can spin up a specialized
one-person tool, **how do we find the trustworthy ones**, and *how will the scientific community react*?

### Could AI create high-quality software?

TigerBeetle (https://tigerbeetle.com/blog/) is a piece of software created with one objective: build extremely high-quality software for a critical mission.
Every line is hand-crafted and revised, both manually and with automated test suites.

Could an AI spec-based development create something of such extreme quality? The TigerBeetle project leader doubts it, and he is likely right; but even he says that AI will be a boon for improving quality.

Of course, creating something like TigerBeetle is a labor of love that requires an extreme dedication to only one goal: high quality.

### I like to build, not to code

I've found out that what I find exhilarating is not writing code, but assembling ideas to build tools. And this new way of working is a joy: you can now build almost anything you dare to imagine. This is an infinite toy.
