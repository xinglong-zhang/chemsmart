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

# W2: does the general repair transfer to a different paper?

W1 could have been repaired locally — teach the harness the one operation that
case happened to need. The test of whether the repair is general is a second
case that shares none of W1's specifics.

| | W1 | W2 |
|---|---|---|
| observable | basis-set-limit energy | equilibrium mole fractions |
| DAG shape | 4 siblings → 1 fan-in | 2 optimizations + frequencies → fan-in |
| convention needed | `exponential_cbs_limit` | `boltzmann_populations` |
| method | HF, cc-pVXZ series | B3LYP/6-31G(d), thermochemistry |
| physical input the model must supply | none | the gauche degeneracy of 2 |
| program the model chose | PySCF | **ORCA** |

The answer key was sealed before dispatch from ChemSmart's own two-stage
PySCF runs: ΔG(gauche−anti) = 0.906 kcal/mol, x(anti) = 0.697,
x(gauche) = 0.303, zero imaginary frequencies at both minima.

## What preparing W2 already showed

**The next instance of the same class, found before the model ran.**
Enumerating the analysis layer for W1 exposed that Boltzmann weighting was
owned but unreachable. Preparing W2 exposed that even once reachable, it could
not express state degeneracy — and n-butane's two gauche forms are
enantiomers. Weighting gauche as one state gives 82% anti; weighting it as two
gives 70%, against a measured ~68%. The multiplicity is not a refinement, it
is the difference between right and wrong. `boltzmann_populations` now takes
optional per-state degeneracies, so ChemSmart owns the formula and the
protocol supplies the physical numbers.

**A preregistration error of mine, corrected by the model.** T3 predicted the
model would plan `opt` then a separate `hess` consuming the optimized
geometry. That is PySCF's canonical shape. The model chose ORCA and wrote one
`opt` project with `freq: true` — which is the *ORCA-correct* form, and
exactly what codex's earlier report identified as the right program-relative
answer. My hypothesis was program-blind; the model was not.

This also means the sealed key is program-specific. A key computed in PySCF
cannot grade an ORCA route without absorbing the B3LYP variant and grid
differences between the two programs, so a matched ORCA reference is being
computed rather than stretching the tolerance to hide the gap.

**Two more general rejection defects.** A `repair_command` call was rejected
twice, identically, with `no counterexample is bound yet`. The message stated
the problem and left nothing to do — listing bound IDs, the earlier repair,
helps only when the registry is non-empty, which is the case where the caller
is already closest to right. Every host registry now declares what fills it,
and the hint appears only when the registry is empty. By contrast the
`promote_project_yaml` collisions were already actionable — they named the
taken ID and listed what was in use — and the model recovered from both
without help.

## W2 result

`terminal_state: complete`, 47 successful tool calls, 4 rejections, 31 provider
turns, 2,425,120 input tokens.

| hypothesis | outcome |
|---|---|
| **T1** population step uses `boltzmann_populations`, not rebuilt | **confirmed** — one node, `{literal: 2, ref: 2, boltzmann_populations: 1}` |
| **T2** gauche degeneracy enters the declared input | **confirmed** — `degeneracies = [1, 2]`, attributed as model-authored |
| **T3** `opt` then separate `hess` | **refuted, and the hypothesis was wrong** — see below |
| **T4** project YAML guard | **held** — `functional: B3LYP`, `basis: 6-31G*`, `freq: true`, `ri_approximation: none` |
| **T5** no rejections | **refuted** — 4, of which 2 are now repaired and 2 the model recovered from unaided |

T2 is the result that matters. The degeneracy input was added days after the
schema the model reads was written for W1, is one optional trailing argument
among many, and nothing in it mentions butane or rotamers. The model used it,
with the correct multiplicities, and the receipt records `[1, 2]` as a number
the model supplied rather than one the toolkit measured. That is the intended
division: **ChemSmart owns the weighting, the protocol supplies the physical
numbers, and the provenance says which is which.**

### The sealed key is program-robust

The key was computed in PySCF; the model chose ORCA. Rather than widen the
tolerance, a matched ORCA reference was computed through ChemSmart:

| | ΔG(gauche−anti) / kcal mol⁻¹ | x(anti) | x(gauche) |
|---|---|---|---|
| PySCF B3LYP/6-31G(d) | 0.9055 | 0.6974 | 0.3026 |
| ORCA B3LYP/6-31G(d) | 0.9234 | 0.7038 | 0.2962 |
| sealed tolerance | ±0.15 | ±0.03 | ±0.03 |

The programs differ by 0.018 kcal/mol and 0.006 in mole fraction — inside
tolerances that were fixed before the program was known. Both agree with the
measured ~68 % anti.

### The most consequential defect this case found

The model wrote `ri_approximation: none`, which is the declared way to ask for
conventional four-index integrals. `ORCARoute` could write `NoRI` but had no
property to read it back, so the preview validator — which parses the
generated input and compares it with the requested settings — reported a
stable critical finding across both rotamers and both basis spellings. The
model cleared it by deleting the key.

Deleting the key is the one edit that lets ORCA's own default density fitting
back in. The harness argued the model out of the input the protocol specified.
This is the sharpest form of the failure the project exists to prevent: not a
model hallucinating a setting, but the harness rejecting a correct one.

# W4 and W5: the loop's remaining observable classes

| case | observable | conventions the plan used |
|---|---|---|
| W1 | basis-set-limit energy | `exponential_cbs_limit` (after repair) |
| W2 | conformer mole fractions | `boltzmann_populations` + degeneracies |
| W3 | vertical excitation, wavelength | `photon_wavelength` |
| W4 | equilibrium bond length and angle | `distance`, `angle` |
| W5 | isodesmic reaction ΔE and ΔG | *none* — and that is correct |

All five convention classes the harness registers have now been exercised by a
live model, on five different observables, in three programs.

## W4: the conventions were reached without help

`{ref: 4, distance: 3, angle: 3}` — one indexed reference per atom, then all
three N–H distances and all three H–N–H angles through the registered
operations, zero arithmetic nodes and zero model constants, plus an equivalence
check across the three bonds that the task never asked for.

## W5: the control the derivation profile needed

A reaction energy is a signed sum of four measured free energies. There is no
convention to own that, and the profile correctly reported arithmetic with no
convention operations without treating it as a finding. A discriminator that
pushed everything toward conventions would have been useless; this is the case
that shows it does not.

W5 also exposed the last convention gap. Asked for the total number of
imaginary frequencies across four species, the model built a **twenty-two node
arithmetic contraption** to count negative numbers — putting the definition of
"imaginary" into its own arithmetic. ChemSmart already owned the concept twice
(`Thermochemistry.imaginary_frequencies`, and the PySCF stationary-point
policy's `imaginary_mode_cutoff_cm1`) and exposed it nowhere.

Registering it produced the clearest before/after of the loop:

| | before | after |
|---|---|---|
| imaginary-mode count | 22 nodes of hand-built comparisons | 5 nodes: `imaginary_mode_count` ×4 + `sum` |
| reaction energy | 4 nodes, pure arithmetic | 4 nodes, pure arithmetic — unchanged |

The convention was adopted exactly where it belongs, and arithmetic stayed
where arithmetic is right.

## What did not work

`repair_command` was called with no counterexample bound in four separate
sessions. Two attempts to prevent it — describing the argument, then stating
the precondition in the tool description itself — both failed. The rejection
message is now fully actionable and the precondition is stated in two places,
and the model still calls it. This is recorded as unresolved rather than
repaired; the cause is not the wording.

## Loop accounting

204 provider turns, 12.8 M input tokens, 379 k output, 241 k reasoning across
eight live sessions. No rate-limit or quota error was returned at any point.

# W6: continuum solvation, and the failure that took three attempts

W6 asked for the electrostatic solvation energy of methanol in a CPCM water
continuum at B3LYP/6-31+G(d), fixed geometry. Sealed first from ChemSmart's own
two single points: **−5.694 kcal/mol** (experiment −5.1).

The declared solvent surface — a six-member `solvent_model` domain and the gate
that refuses a model without an explicit `solvent_id` — had never been touched
by a live run. It worked. The model authored

    solv:
      basis: 6-31+G(d)
      functional: b3lyp
      solvent_id: water
      solvent_model: cpcm

which is the sealed reference field for field, and put the gas-phase leg in its
own section. Its `delta-g-solv` node is `{subtract: 1, convert: 1}` — continuum
minus gas in kcal/mol, the right sign convention — plus a validation node it
was not asked for.

## The premature-tool-call failure, resolved on the third attempt

Across six sessions the model called a tool whose precondition could not yet be
satisfied: `repair_command` with no counterexample bound, four times, and
`assess_program_candidate` with no claim evidence bound.

- **Attempt 1** — describe the argument. Failed; the next session did it again.
- **Attempt 2** — state the precondition in that tool's own description.
  Failed; W6 arm A then did it with a *different* tool.
- **Attempt 3** — stop treating it as a property of one tool. Six tools take an
  argument indexing a host registry that only something else can fill. The
  precondition is now derived for all six from the same producers table the
  rejection path uses, so the sentence before the call and the message after it
  cannot drift.

| W6 | arm A | arm B |
|---|---|---|
| rejections | 3 | **1** |
| premature-tool calls | 1 (`assess_program_candidate`) | **0** |
| provider turns | 41 | 30 |
| input tokens | 2,840,881 | 2,132,637 |

The single arm-B rejection is a content issue in `plan_scientific_workflow`,
not a precondition violation.

The lesson is about method, not about this tool. The first two attempts were
aimed at the instance that happened to be visible. Only the third asked what
class the instance belonged to, and it is the only one that moved the
behaviour — on one paired run, which is an observation, not an effect size.

# W7: the four-session failure, and how three fixes missed it

W7 asked for the ethane rotational barrier from single points at fixed
staggered and eclipsed geometries, plus the H–C–C–H torsion in each. Sealed
first: **3.328 kcal/mol** (experiment 2.9; the C–C bond is not relaxed),
torsions 60.0° and 0.0°.

## A gap closed by prediction rather than by damage

Sealing the case exposed that `distance` and `angle` were registered
conventions and `dihedral` was not — even though a torsion is the third
standard internal coordinate and the one a rotational barrier is *defined*
along. That gap was predictable without running anything, so it was closed
first. The model then used it immediately: `{ref: 4, dihedral: 1}`, four
indexed references and one call, zero arithmetic nodes.

The sign matters and is pinned: reversing the fourth atom reverses the torsion,
which is what distinguishes a P from an M helix. Degenerate geometries are
refused by name rather than returning a number.

## The failure I got wrong three times

Across four sessions the model called `repair_command` with no counterexample
bound. I attributed it to wording and tried, in order:

1. describing the argument — failed;
2. stating the precondition in that tool's description — failed, and the next
   case did it with a *different* tool;
3. deriving the precondition for all six tools in that class from one table —
   appeared to work on one paired arm, then W7 arm A had three such calls.

Each fix was plausible. None held. Reading the code answered it in one grep:
`CommandCounterexampleV1` is never constructed anywhere in the package and
`register_counterexample` has no callers. **The registry is empty for the whole
lifetime of every session, so the tool cannot succeed.** It is an unfinished
feature that was advertised as a usable repair path, and the model was
reasonable to reach for it.

| W7 | arm A (advertised) | arm B (withheld) |
|---|---|---|
| provider turns | 37 | **27** |
| input tokens | 2,512,603 | **1,936,831** |
| rejections | 4 | **1** |
| `repair_command` calls | 3 | **0** |

The remaining arm-B rejection is an artifact-ID collision, from which the model
recovered unaided.

The handler, contract, and invocation fields are all kept, so re-exposing the
tool is a one-line surface change in the same commit that wires a producer —
and a test fails the moment such a producer appears, which is the signal to do
it. The general rule is now gated on all three surfaces: no tool may be
advertised whose required registry nothing in the runtime fills. That is
exactly what the capability registry already does when it declares PySCF `td`
as `execution_supported: false` rather than offering it as runnable.

The methodological point is the one worth keeping. Three fixes aimed at the
symptom looked reasonable and were cheap to write, and the cost of not checking
whether the tool could work at all was four wasted sessions.
