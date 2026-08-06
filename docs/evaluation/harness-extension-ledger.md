# Harness extension ledger

What was added to the ChemSmart agent harness across the W1–W10 adversarial
campaign, which experiment produced each addition, and the general effect each
one was defined to obtain.

The companion narrative is
[hf-basis-set-limit-case.md](hf-basis-set-limit-case.md), which reports the
cases in the order they ran. This file is the cross-cut: it groups the
extensions by the effect they were meant to have, so that a reader can ask
"what class of future failure does this cover?" rather than "what happened in
W5?".

## How to read an entry

Every row binds four things, and no row may be filled from prose:

| field | bound to |
|---|---|
| experiment | a run directory under `scratchpad/runs/`, with its `events.jsonl` |
| observation | an event payload, a rejection message, or a node count from that log |
| definition | a commit and the symbol it introduced |
| gate | the test that fails if the definition is removed or drifts |

Where an entry says *"closed by prediction"*, no live failure preceded it: the
gap was found while sealing a case, and the entry says so rather than implying
a run that never happened.

## The campaign, measured

Sixteen live `agent plan` sessions against `deepseek-fallback` /
`deepseek-v4-flash`, ten distinct cases, six of them paired-arm.

| | |
|---|---|
| live sessions | 16 |
| provider turns | 512 |
| input tokens | 40,492,101 |
| output / reasoning tokens | 851,755 / 547,711 |
| successful tool calls | 792 |
| rejections | 40 |
| rate-limit or quota errors | **0** |

Three of the ten cases are built from genuinely fetched article text (W1, W8,
W9); the other seven are protocols constructed for the benchmark and anchored
to established experimental values, and their task files say so. Every case
was sealed from ChemSmart's own execution of its protocol **before** dispatch.

## The three mechanisms

Every extension below traces to one of three failure classes. Naming them is
the point of the ledger: an extension aimed at a class transfers, an extension
aimed at an instance does not, and the campaign has one clean demonstration of
each.

**M1 — The harness holds a precise typed fact and spends it on refusal instead
of on telling the caller.** Found at twelve seams. The refusal was always
correct; it was never *usable*. This class is why a rejection message and the
sentence that prevents the rejection must come from the same table.

**M2 — A convention ChemSmart owns, but does not expose, gets rebuilt by the
model out of arithmetic.** This is the failure the project exists to prevent,
stated in its own terms. Its signature is a large arithmetic DAG where a single
named operation belongs, and it is measurable: node count before versus after.

**M3 — The harness is wrong about itself.** It advertises a tool nothing can
satisfy, writes a setting it cannot read back, or refuses a project shape that
is correct. This class is the most damaging because the model is behaving
correctly when it hits it, and its recovery is to abandon a capability that
actually works.

---

# A. Refusal that repairs

**General effect.** A refusal must carry everything needed to fix it: the
indexed path, the value received, the rule broken, and — when a registry is
empty — what fills it. The target property is that *no rejection requires a
second rejection to interpret*.

| experiment | observation | defined | commit |
|---|---|---|---|
| W1 arm A | `analysis_nodes[].outputs[].output_id does not match its pattern` — eleven nodes, no index, no value, no pattern | generic argument validator rewritten: indexed array paths, quoted value, quoted rule, union types, output bounded at 120 chars | `647cc47a` |
| W2 | `repair_command` rejected twice, identically, with `no counterexample is bound yet` — nothing to act on | `REGISTRY_PRODUCERS`; an empty registry now names what would fill it, and the hint appears only when it is empty | `fd0ffa90` |
| W2 | a `tool_failed` event records only `arguments_sha256`, so nothing durable said what the validator objected to | `_public_validator_findings()` — up to 8 findings recorded, values ≤120 chars, explicit truncation marker | `baed41f3` |
| W3 | a field that does not apply had no way to say so and was rejected for being absent | `_nullable_positive_number()` → `{"type": ["number","null"], …}` | `ee9cab0a` |
| W4 | a Gaussian excited-state job hit `AttributeError: 'NoneType' object has no attribute 'merge'` three times; the session concluded the program could not preview that job at all | `require_jobtype_settings(...)`, used by all 20 per-jobtype accessors in `settings/gaussian.py` and `settings/orca.py` | `40262735` |
| W9 | a thermochemistry selector was refused by the extraction tool that does not own it | `QUANTITIES_FROM_ANOTHER_TOOL` — a refused selector names the tool that does own the quantity | `da9aa0a5` |
| W9 | `analysis output unit must not be empty` appeared in three separate cases, always on a count, naming neither the output nor the value a dimensionless quantity should carry | one shared `_UNIT_HINT` naming the output and all three acceptable spellings, used by both planes that raise it | `c41b3a68` |
| pre-W1 | a rejected DAG node named neither the node nor the offending value | node id and value quoted in `workflows.py` | `bc039526`, `bfdabbd1`, `053bc81a` |

**Why this class kept recurring.** Each fix was local to the validator that
happened to be visible. `647cc47a` is the one that generalised: an earlier
commit had already added exactly those three facts to *one* dataclass
validator, and the generic validator guarding every tool still had the defect.

---

# B. Declaration before the call

**General effect.** Any rule the harness will enforce *after* a call must be
stated *before* it, derived from the same table that enforces it, so the two
cannot drift. The target property is the elimination of guess–reject–retry
cycles, which are the campaign's dominant token cost.

| experiment | observation | defined | commit |
|---|---|---|---|
| W1 arm A | the operation enum carried no descriptions; 24 bare names had to be read from spelling alone | `OPERATION_DESCRIPTIONS`, pinned to `_OPERATIONS` by a runtime check so a new operation cannot ship undescribed | `8172d368` |
| W3, W4 | required arguments were named but not explained | `ARGUMENT_DESCRIPTIONS` keyed by **argument name**, not by tool — 93/93 required arguments described, applied to all three exposure surfaces | `b0c397f4` |
| W4 | identifier patterns were discoverable only by violating them | `IDENTIFIER_ARGUMENTS` / `OPTIONAL_IDENTIFIER_ARGUMENTS`, with the class **probed against the live validators** in a test rather than asserted | `6e2bea65` |
| W5 → W6 | `repair_command` and `assess_program_candidate` called before their registries could be filled, across six sessions | preconditions derived for all six registry-reading tools from the same producers table the rejection path uses | `1b3ab67f`, `0135584b` |
| goal 4 | every waiting node reported the identical string `"await producer outputs"` | per-node `waiting_on`, `unsatisfied_inputs` (`producer_node_id.producer_output_id`), `feeds`, and a `next_action` naming them | `3d804661` |
| W8 | `project_section_names` was declared on `ProgramCapabilityV1` and projected for **no** program | the loader projection filled in | `7f184da3` |

**A negative result inside this class, kept.** W8's paired arm isolated
`7f184da3`, and the arm carrying it cost **2.7× the input tokens** of the arm
without it (11.6 M vs 4.3 M) with four rejections of kinds not seen before.
One run per arm and a single additive field, so no causal reading is available
— but the projection did not demonstrably help, and that is recorded as it
happened rather than smoothed.

---

# C. Convention ownership

**General effect.** Where ChemSmart owns a computational-chemistry convention,
the model must be able to *select* it rather than rebuild it from primitives.
This is the FUNDAMENTALS claim made measurable: "prevents hallucinations"
becomes the falsifiable statement *no convention the hub owns is re-derived in
the analysis DAG*, and its metric is expression-node count.

The selectable convention vocabulary **doubled** over the campaign, 6 → 12.

| experiment | observation | defined | commit |
|---|---|---|---|
| W1 arm A | basis-set limit built from **15 primitive nodes**, `E_CBS = (E2·E4 − E3²)/(E2 − 2E3 + E4)`, with scale factors 2 and −1 entering as model-chosen numbers | `exponential_cbs_limit` accepts three scalars **or** a series; the closed form's equal-spacing precondition stops being silent | `9a1b4b46`, `8172d368` |
| W1 → W2 | Boltzmann weighting was owned in `analysis/aggregation.py` and selectable nowhere | `boltzmann_populations`, `boltzmann_average`, `correlation_inverse_power_cbs_limit`, `linear_fit_intercept` registered; plus a **drift gate** that fails if any public convention in `aggregation.py` is neither exposed nor excluded with a reason pointing at a real operation | `a4bdec95` |
| W2 prep | n-butane's two gauche forms are enantiomers; weighting gauche as one state gives 82 % anti, as two gives 70 %, against a measured ~68 % — the multiplicity is the difference between right and wrong | optional per-state `degeneracies` on both Boltzmann operations, so ChemSmart owns the formula and the protocol supplies the physical numbers | `ba08cbf4` |
| W5 | asked for the total imaginary-mode count across four species, the model built a **22-node arithmetic contraption** to count negative numbers, putting the *definition* of "imaginary" into its own arithmetic | `imaginary_mode_count(frequencies, cutoff_cm1)`, the concept ChemSmart already owned twice and exposed nowhere | `2420858d` |
| W7 prep — **closed by prediction** | sealing the case exposed that `distance` and `angle` were conventions and `dihedral` was not, though a torsion is the third internal coordinate and the one a rotational barrier is *defined* along | `dihedral`, signed, with degenerate geometries refused by name rather than returning a number | `719ac2aa` |

**Before/after, the clearest measurement in the campaign:**

| | before | after |
|---|---|---|
| basis-set limit (W1) | 15 nodes, 2 model-authored scale factors | 1 node, 0 model constants |
| imaginary-mode count (W5) | 22 nodes of hand-built comparisons | 5 nodes: `imaginary_mode_count` ×4 + `sum` |
| isodesmic reaction energy (W5) | 4 nodes, pure arithmetic | 4 nodes, pure arithmetic — **unchanged** |

The third row is the control. A discriminator that pushed everything toward
conventions would be useless; a reaction energy *is* a signed sum, there is no
convention that owns it, and the profile correctly reported arithmetic with no
convention operations without treating it as a finding.

**Transfer, which is the actual claim.** `imaginary_mode_count` (registered
after W5) and `dihedral` (registered before W7) were both used unprompted in
**W9**, on a different molecule at a different level, in single nodes. Neither
case mentioned the operation. That is the evidence that these are class-level
repairs and not case-level ones.

---

# D. Derivation attribution

**General effect.** A reported number must be able to say which of its digits
came from a measurement and which the model typed in. Without this, C is not
measurable and transparency is asserted rather than auditable.

| experiment | observation | defined | commit |
|---|---|---|---|
| W1 | from the *same* executed cc-pVTZ and cc-pVQZ jobs, the three-point fit gives −76.067424 E_h and the two-point form with a model-chosen α = 5.34 gives −76.067245 E_h — both **equally well evidenced** by the receipt closure | `ModelAuthoredConstantV1` with roles `extrapolation_exponent`, `literal_value`, `power_exponent`, `scale_factor`; each expression output carries the sorted set of constants that actually reach it, attributed to the introducing node, following the DAG and bound into the receipt digest | `b8e86a9a` |
| W1 | nothing distinguished "arithmetic because arithmetic is right" from "arithmetic because the convention was unreachable" | the derivation profile: `CONVENTION_OPERATIONS`, `ARITHMETIC_OPERATIONS`, `_PLUMBING_OPERATIONS`, and a per-output count | `4588b8f4` |

**The check that this does not over-report.** W9's scaling factor 0.9135 comes
from the paper, not from a measurement, and *should* be recorded as
model-supplied. That expectation was written into the answer key before
dispatch precisely as a control, and the profile behaved: a reader can see
exactly which digit came from the article and which from the calculation.

---

# E. The harness must not be wrong about itself

**General effect.** The harness must not advertise what it cannot complete,
must not write a setting it cannot read back, and must not refuse a project
shape that is correct. Every failure in this class ends with the model
abandoning a capability that works — the most expensive outcome available,
because the model was right.

| experiment | observation | defined | commit |
|---|---|---|---|
| **W10** | a gas-phase B3LYP/aug-cc-pVTZ single point written as `gas: {functional, basis}` produced a **0-byte `.inp`**, `ValueError`, and a durable event carrying `findings: []`, three times, including on a minimal project. The session concluded ORCA could not materialise a single point in this environment and **switched program** | `sp` inherits `gas:` when no `solv:` exists, without inheriting its `freq`; plus a gate that derives section→jobtype routing from the loader and sweeps the whole capability registry | `afc20200` |
| W2 | the model wrote `ri_approximation: none` — the declared way to ask for conventional four-index integrals. `ORCARoute` could write `NoRI` and had no property to read it back, so the preview validator reported a stable critical finding and **the model cleared it by deleting the key** — the one edit that lets ORCA's default density fitting back in | `ri_approximation` read-back property on `ORCARoute` | `2b2c6720` |
| W7 | `repair_command` called with no counterexample across **four sessions**; three fixes aimed at wording all failed. One grep answered it: `CommandCounterexampleV1` is never constructed and `register_counterexample` has no callers — **the registry is empty for the lifetime of every session** | the tool withheld from the exposed surface (handler and contract kept); a gate on all three surfaces: no tool may be advertised whose required registry nothing in the runtime fills | `0bfb8c87` |
| W1 | `agent run` needs an approval binding `workflow_sha256`, and nothing emitted an approval from a completed plan | `WorkflowApprovalRequestV1`, `build_workflow_approval_request()`, `approve_workflow_request()` | `57c04ea9` |
| W2 | a compiled command could not be attributed to its DAG node | `command_compiled` now carries `node_id`, `program`, `jobtype`, `execution_target`, `charge`, `multiplicity`, `display_command`, and two binding digests | `72a7e933` |
| pre-W1 | the section-shape gate refused per-jobtype sections the loader really reads | the gate unions in `capability.jobtypes` rather than maintaining a second list | `457fbabb` |

**The methodological result of this family**, and the most useful thing the
campaign produced: in W7, three plausible, cheap fixes aimed at the visible
symptom all failed, and cost four sessions. Only the fourth attempt asked what
*class* the instance belonged to. In W10 the same discipline was applied first
— the model's narrative ("ORCA fails here") was treated as an interpretation to
audit, the receipt was read instead (`size_bytes: 0`, `exception_class:
ValueError`), and the defect was reproduced and root-caused in the loader
before anything was written.

---

# F. Scale honesty

**General effect.** The harness must refuse a request it can no longer serve,
rather than serving it more and more expensively.

| experiment | observation | defined | commit |
|---|---|---|---|
| W8 arm B | context grows on **every** turn without exception (51/51, 78/78 measured), 16 k → 249 k; total cost is near-quadratic in turn count (`total/N²` within a 3× band across 10–79-turn runs); the worst arm spent **11.6 M input tokens** | `estimate_request_input_tokens()` and a pre-request guard emitting `TURN_BLOCKED` with rule `budget.context_would_be_exceeded` | `c4d56de2` |

**This corrects an earlier claim of mine.** Token deltas reported between arms
in W1–W6 were substantially measuring *turn count*, not efficiency. They should
not be read as independent efficiency claims, and the tables above no longer do.

---

# Gates

Each family is held by a test that fails on removal or drift. These are the
instruments; the extensions are what they protect.

| gate | fails when |
|---|---|
| `test_actionable_argument_rejections.py` | a rejection omits path, value, or rule |
| `test_convention_exposure.py` | a public convention in `aggregation.py` is neither exposed nor excluded with a reason naming a real operation |
| `test_operation_vocabulary.py` | an operation ships without a description |
| `test_argument_descriptions.py` | a required argument ships undescribed |
| `test_identifier_constraints.py` | the declared identifier class disagrees with the live validators |
| `test_reachable_tool_surface.py` | a tool is advertised whose registry nothing fills — **and when a producer appears**, which is the signal to re-expose `repair_command` |
| `test_declared_sections_feed_their_jobtypes.py` | a declared project section does not reach a job type that reads it |
| `test_model_authored_constants.py` | a model-typed constant escapes attribution |
| `test_context_growth_guard.py` | the budget guard stops refusing |

Suite state: **2332 passed**, 10 failed, all ten pre-existing and none in the
files above — verified by running them against committed `HEAD` in a clean
checkout.

---

# W11: the first graded number, and why the gap was not what I said it was

W11 re-ran the W9 protocol as an **execution** arm. Information was controlled
to the methodology section plus `water.xyz`; the study's results were withheld
and only its calibrated factor 0.9135 was given, as the protocol states it.

The plan was right. 26 turns, 1,618,784 input tokens, terminal state
`complete`, two rejections — both of them repairs from this campaign firing
live, and the model recovered from each without help:

| rejection | family |
|---|---|
| `artifact_id must be a lower-case public identifier … it is empty` | A |
| `analysis output 'n_imaginary' (count) declares no unit; a dimensionless quantity … uses '1'` | A |

It chose HF/6-31G(d) with `density_fit: false`, planned `opt` then `hess`,
bound the hessian to the optimizer's `optimized_geometry` as a producer edge,
and used `imaginary_mode_count` — registered after W5 — unprompted for the
third case running.

**Graded against the sealed key**, executing its own planned nodes with its own
project YAMLs:

| quantity | agent | sealed | tolerance | |
|---|---|---|---|---|
| `zpve_harmonic_kj` | 60.1445 | 60.145 | ±0.5 | **PASS** |
| `zpve_scaled_kj` | 54.9420 | 54.942 | ±0.5 | **PASS** |
| `n_imaginary` | 0 | 0 | 0 | **PASS** |

Frequencies 1826.0168, 4055.5119, 4173.8509 cm⁻¹, matching the sealed
independent reproduction to every printed digit.

**Two honest limits.** This was deterministic host-side execution of the
agent's own planned nodes, not `agent run --approval-file`; and the sealed key
came from ChemSmart's own execution of the same protocol, so close agreement is
expected *given the same level of theory*. What the case tests is that the
agent chose the level, the two-stage workflow, the geometry handoff and the
post-processing correctly from method text and an xyz alone. It did.

## The gap was worse than recorded, and is now closed

I had recorded the reason no number was ever graded as "no approval file is
emitted from a plan session". Reading the execution path end to end shows that
was wrong. `execute_program_node` refuses when its frozen approval is `None`
— *"legacy V1 approval is preview-only; Runtime V2 frozen approval is required
for execution"* — and `frozen_workflow_approval=` was passed at **no call site
in the package**. No approval file could have driven execution whatever it
contained.

That is the W7 `repair_command` finding again — a path advertised that the
runtime cannot let succeed — on the most important capability in the harness.
Found the same way: by asking whether the thing could work at all before
theorising about why it did not. The frozen body now travels in the approval
file beside the V1 approval and the composition threads it through; deriving it
in-session was rejected because it would make its own `plan_sha256` check
compare the plan against itself.

## The convention the case earned

The model's ZPVE was `sum → scale 0.5 → convert(cm⁻¹→kJ/mol) → scale 0.9135`.
The `convert` node **cannot run**: a wavenumber is `FREQUENCY` and a molar
energy is `ENERGY`, so the unit engine refuses, correctly — `h·c·N_A` is a
spectroscopic convention, not a unit identity. The session predicted the
refusal in its own report and proposed hand-multiplying by `h·c·N_A`, which
would have moved a physical constant into model arithmetic.

`harmonic_zero_point_energy` now owns both the factor of one half and the
conversion, reproducing the sealed value to eight significant figures.
Selectable conventions: **13**.

**My own drift gate missed this.** It sweeps `analysis/aggregation.py`; this
convention was owned in `analysis/thermochemistry.py`. Widening the sweep is
the next repair, and it is listed below rather than claimed as done.

---

# What is not closed

Stated plainly, because the value of the ledger depends on it.

1. **`agent run` has not itself executed a case.** The wiring defect is fixed
   and one case is graded, but the end-to-end CLI path additionally requires a
   fresh session to reproduce byte-identical project artifacts, which has not
   been demonstrated. Nine sealed keys remain ungraded.
2. **The convention drift gate has a known blind spot**: it sweeps
   `analysis/aggregation.py` only, and `analysis/thermochemistry.py` owns
   conventions too. That blind spot is exactly what let the zero-point energy
   stay unexposed through nine cases.
3. **Goal 4 is named, not wired.** `select_task_dependency_context()` implements
   the dependency-context policy with a budget, and is called only from the
   experiments module — never from the live loop. Family B addresses the
   *frontier*; the preceding-task context itself is still not delivered.
4. **One run per arm.** Every paired comparison in this ledger is an
   observation, not an effect size.
5. **Gaussian declares zero `project_parameter_domains`** where ORCA declares
   five. A declaration gap, not a verified absence.
6. **An adjacent asymmetry, reported rather than repaired.** A `solv:`-only
   project gives `sp` `freq: True` from the defaults. Changing it would move
   existing projects, unlike the `gas:` case which was an unconditional
   failure, so it is recorded here and the gate skips it explicitly.
7. **`project_section_names` remains a half-fact.** It lists which sections are
   accepted, not which job type reads which. `afc20200` removed the need to
   know for the case that failed; the general question is open.

---

# W10, in full

The case that produced family E's newest entry, since it postdates the
narrative report.

Vertical ionization energy of water, B3LYP/aug-cc-pVTZ, both states at the
neutral geometry. Sealed before dispatch: E(neutral) = −76.46614403 E_h,
E(cation) = −75.99546614 E_h, **IP = 12.808 eV** (experiment 12.62), tolerance
±0.15 eV. It closed a coverage gap in the benchmark design: every prior case
was a neutral closed-shell singlet.

| | |
|---|---|
| turns / input tokens | 24 / 1,420,086 |
| tool calls | 41 successful, **0 rejected** — the first zero-rejection session |
| identities bound | `{charge: 0, multiplicity: 1}` **and** `{charge: 1, multiplicity: 2}` |
| stage order | `pyscf-sp-neutral → pyscf-sp-cation → extract ×2 → expr-ip → render-claims` |

The session recorded four alternatives it rejected with reasons, six
assumptions, four diagnostics, and four uncertainties — including that PySCF
resolves the requested `b3lyp` literal to the Gaussian-convention alias
`b3lypg`/`vwn3_gaussian`, cited to the host's own functional-resolution
receipt and explicitly *not* claimed as component-level materialisation.

That record is also what surfaced the defect: it is the only reason a 0-byte
`.inp` three levels down in a preview receipt became a reproducible bug rather
than a session that quietly used a different program.
