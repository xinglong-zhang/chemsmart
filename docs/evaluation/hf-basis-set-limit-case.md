# W1: a paper case the harness can decide end to end

## Why this case

Codex's campaign report closed with a rule: no further provider call until a
case has either a complete numerical answer key or a complete expert-authored
canonical YAML–CLI DAG answer. Every case tried before it failed the first
condition for the same reason — the systems were too large to run, so no
number could be produced locally and no numerical claim could be graded.

This case removes that obstacle rather than working around it. PySCF is a
library backend whose executable is a Python interpreter, so a small-molecule
Hartree–Fock protocol runs to completion on a laptop in seconds. That makes an
answer key obtainable **through ChemSmart itself**, not only from the paper.

**Source.** A. Halkier, W. Klopper, T. Helgaker, P. Jørgensen, P. R. Taylor,
*Basis-set convergence of the energy in molecular Hartree–Fock calculations*,
[Chem. Phys. Lett. **302** (1999) 437](https://doi.org/10.1016/S0009-2614(99)00179-7).

The paper fixes water at its experimental equilibrium geometry
(r(O–H) = 0.9572 Å, ∠H–O–H = 104.52°), computes Hartree–Fock total energies in
cc-pVXZ, and establishes that the Hartree–Fock basis-set error decays
exponentially in the cardinal number rather than as an inverse power — so the
limit must be reached by an exponential extrapolation.

## The answer key exists, and ChemSmart produces it

Four `chemsmart run pyscf sp` jobs, one per basis, each driven by a project
YAML whose entire scientific content is `ab_initio: hf`, `basis: cc-pvXZ`,
`freq: false`, `density_fit: false`:

| basis | ChemSmart / PySCF | paper (Dalton) | difference |
|---|---|---|---|
| cc-pVDZ | −76.026798697 | −76.026799 | +0.30 µE_h |
| cc-pVTZ | −76.057168515 | −76.057168 | −0.51 µE_h |
| cc-pVQZ | −76.064835339 | −76.064835 | −0.34 µE_h |
| cc-pV5Z | −76.067090832 | −76.067091 | +0.17 µE_h |

Every deviation is inside the paper's own printed precision. The
three-parameter exponential limit through (D,T,Q) is −76.067424 E_h, against
the near-Hartree–Fock limit of −76.067488 E_h the paper quotes — 0.064 mE_h.

The case is therefore decidable, and the four required quantities and their
tolerances were sealed before the first provider call.

## What the first live arm showed

`deepseek-v4-flash`, `agent plan`, 14 provider turns, 40 successful tool calls,
2 rejections, terminal state `planned`.

**The fundamentals held.** The model authored four PySCF project YAML
documents that are substantively identical to the sealed reference, bound the
approved geometry by digest, compiled four canonical `chemsmart run … pyscf …
sp` commands that passed safe preview, and refused to report any number
without execution.

**Three harness defects surfaced, each general.**

### 1. A rejected argument named the field but not the value or the rule

    tool argument analysis_nodes[].outputs[].output_id does not match its
    pattern

Eleven analysis nodes, no index, no value, no pattern. An earlier commit had
added exactly those three facts to one dataclass validator; that repair was
local, and the generic argument validator that guards every tool still had the
defect. It now indexes array paths and states the value it received against
the rule it broke, with rendered values bounded so a refusal cannot become an
echo channel.

This also repairs the audit trail: a `tool_failed` event records only
`arguments_sha256`, so nothing durable could say what the model sent. The
message is recorded, so a self-describing message makes the rejection
reconstructable.

### 2. The operation that owns a convention was unreachable, so the model rebuilt it

The model produced the basis-set limit from fifteen primitive nodes:

    E_CBS = (E2*E4 - E3^2) / (E2 - 2*E3 + E4)

ChemSmart already owns that convention as `exponential_cbs_limit`. Two things
kept the model away. The operation enum carried no descriptions, so twenty-four
bare names had to be read from spelling alone; and the operation accepted only
a single three-element input while the extraction plane yields one scalar per
calculation, so a model holding three separate energies had no way in.

This is the failure mode the project exists to prevent, stated in its own
terms: the model re-derived a computational-chemistry convention instead of
selecting the canonical one. It is also not cosmetic — the closed form is
valid only for equally spaced cardinals, nothing checked that, and its scale
factors of 2 and −1 enter the result as numbers the model chose.

`exponential_cbs_limit` now accepts three scalar energies as well as a series,
and every operation declares what it computes and what it expects, pinned to
the operation set so a new operation cannot ship undescribed.

### 3. A derived value could not say which of its digits the model supplied

The receipt closure proved which measurements a value depended on and said
nothing about constants introduced on top of them. From the same executed
cc-pVTZ and cc-pVQZ jobs, the three-point fit gives −76.067424 E_h and the
two-point form with a model-chosen α = 5.34 gives −76.067245 E_h; both were
equally well evidenced. Each expression output now carries the sorted set of
literal values, scale factors, power exponents and extrapolation exponents
that actually reach it, attributed to the node that introduced them, following
the DAG and bound into the receipt digest.

## The paired arms

Identical task, identical workspace, identical provider
(`deepseek-fallback` / `deepseek-v4-flash`). Arm B differs only by the three
repairs above. Hypotheses were written down before arm B was dispatched.

| | arm A | arm B |
|---|---|---|
| terminal state | `planned` | `complete` |
| tool calls (ok / rejected) | 40 / 2 | 29 / **0** |
| provider turns | 14 | 10 |
| input tokens | 829,712 | 435,932 |
| output tokens | 46,443 | 28,864 |
| reasoning tokens | 28,005 | 18,784 |
| expression nodes for the limit | **15** | **1** |
| operation used | hand-built closed form | `exponential_cbs_limit` |
| model-authored constants in `e_hf_cbs` | scale factors 2 and −1 | none |
| analysis nodes | 9 | 6 |

The guard held: arm B's four project YAML documents are identical in
scientific content to arm A's and to the sealed reference — `ab_initio: hf`,
`basis: cc-pvXZ`, `freq: false`, `density_fit: false` — and both arms bound
the same approved geometry digest `202250f8d1e9a7eb…`. The manipulation moved
the post-processing route and nothing else.

Honest limits on this pair. One run per arm, so this is an observation, not an
effect size. Three commits differ rather than one, though only the operation
vocabulary can plausibly change a plan. The token reduction is confounded with
the shorter DAG and with two fewer rejection-and-retry cycles, so it should
not be read as an independent efficiency claim. And arm B's calculation plan
declares `required_observables: ['energy']` where arm A declared the four
named energies — a small regression in the calculation half that the
manipulation did not target and did not fix.

## What this case does not yet show

- **No number has been graded.** `agent plan` does not execute. The numerical
  answer key is sealed and unused until an execution arm runs.
- **`agent run` is blocked on approval construction.** The approval file binds
  `workflow_sha256`, the digest of a workflow a fresh session would have to
  reproduce byte-for-byte, and nothing in the harness or the documentation
  emits an approval from a completed plan. This is the structural reason the
  campaign has never closed a numerical case, and it is not fixed here.
- **One run per arm.** The paired arms below are observations, not effect
  sizes.
