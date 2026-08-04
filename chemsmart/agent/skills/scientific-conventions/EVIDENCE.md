# Evidence for `scientific-conventions`

Every knowledge item in `SKILL.md` is admitted only after a paired dry run on a
constructed benchmark case showed the targeted defect change. No chemistry
engine was executed for any of these observations.

Paired runs use the same frozen task, the same approved geometry artifact, and
the same provider profile (`alibaba-token-plan` / `qwen3.8-max`,
`reasoning_effort=xhigh`). The arms differ only in `CHEMSMART_AGENT_SKILLS`.

## KB-01 — CH₂ adiabatic singlet–triplet splitting

Task: determine `T0` for methylene with the design hidden (no method, basis,
charge, multiplicity, stage list, or DAG shape supplied).

| Arm | Session | `result_sha256` | Event-stream head |
|---|---|---|---|
| skills off | `live-20260804T081201341496Z-8d520ac033-67097e60` | `42fd43f5a371f8afb5a648d7742b76a2582c82b5faa38ca8165d6cbaa6f76569` | `ed5616153165ce3a5a7d94255e549e3c2fa9c37290f04f0a54ee904ab766804e` |
| skills on | `live-20260804T092617008368Z-8d520ac033-6473953e` | `5bbb172da625cd11e59ec1de4f3756429ddf5fb6f0c6b31003bf97351d5aaa37` | `e8b73f3659fa1bbd24239570d96b9093eb3cee48149ee4df58959dbd36af5094` |

Both arms reached `terminal_state: complete` and produced the same four-node
DAG (two sibling OPTs, two HESS nodes left unresolved behind their own producer
edge), so the comparison is not confounded by one arm failing to plan.

### Items this pair admits

| Item | Skills off | Skills on |
|---|---|---|
| Term-value orientation | `T0 = E_triplet − E_singlet` — inverted, so the literature `+9.0` would print negative | `T0 = E(singlet) − E(triplet)`, "positive if the triplet is the ground state" |
| Established state assignment | "the identity of the ground state cannot be reported in this session" | "the triplet `X ³B₁` is the **established ground state**", attributed to the skill and labelled advisory prior knowledge, with the calculation left to confirm or contradict it |
| Convergence for difference quantities | `scf_tol: null`, `scf_maxiter: null` | `scf_tol: 1.0e-09`, `scf_maxiter: 200`, `defgrid3`, with the reason stated |

The skills-on arm called `consult_domain_skill` and named the skill in its
public answer, so the attribution is observable rather than inferred.

Both arms independently identified that the loader exposes no post-HF or
multireference method, so no item claims to close that gap.

## KB-02 — vertical and adiabatic first ionization energy

Task: report both the vertical and the adiabatic first ionization energy of a
small closed-shell molecule, design hidden. This case probes §2 (adiabatic
versus vertical is a geometry convention), which is a structural claim: the
correct workflow shape differs by definition alone, so the DAG either has it or
does not and no numeric judgement is involved.

| Arm | Session | `result_sha256` | Event-stream head |
|---|---|---|---|
| skills off | `live-20260804T093753430200Z-1503fa1ddd-994955bd` | `6d2ec93b4ebe9e32c3f0ae8fdd43a52d17a827c8341a62bfd6393a982954a021` | `a745b1f3442fc56d1d39e0fefa244051746f7c30aeaf775d8382451a153ea2ed` |
| skills on | `live-20260804T094700392450Z-1503fa1ddd-140cffba` | `443232cf78330e7419ff9e7aaeccb0bd1c30cbe9b3ab7c646cce3fbdf495d4bf` | `3af663ce48ecacb6d34c20df65377fe2add4e32b44e65b86ba1a9933edfafe05` |

Both arms reached `terminal_state: complete`.

### What the pair admits

| | skills off | skills on |
|---|---|---|
| Geometry the vertical quantity is evaluated at | a final-state single point at the **supplied initial** geometry, with **no data edge** from the relaxed initial-state optimization | a final-state single point carrying an explicit data edge **from the relaxed initial-state optimization** |
| Adiabatic branch | a final-state optimization not tied to the relaxed initial state | relaxed initial state → final-state optimization → final-state energy, as an explicit producer chain |
| Data edges | 2 | 4 |

Without the skill the vertical quantity was not bound to the relaxed
initial-state geometry its definition requires. With the skill the binding is
explicit in the graph. This admits §2 as a **general** item: the claim is about
the adiabatic/vertical distinction for any state-to-state quantity, and the case
is one probe of it, not its scope.

## Interference and validity notes

An earlier skills-on attempt (`live-20260804T091515229869Z-...`) terminated on a
provider `rate_limited` response before previewing. It is retained as an
inconclusive observation and is **not** used as evidence for any item; the
re-run above is the admitted arm.

The KB-02 skills-on arm was in flight when `SKILL.md` was revised from 0.1.0 to
0.2.0 for generality. The adiabatic/vertical section is textually equivalent
across both revisions, so the item this pair admits is unaffected; no other
section is claimed from this pair.

## Generality

The skill states principles, not facts about the probe systems. The probe
systems are deliberately absent from `SKILL.md`, and
`tests/agent/test_domain_skills.py` fails the build if any of them is named
there — overfitting the knowledge to its own benchmark is a regression, not a
refinement. Each admitted item is a rule that generates the probe's answer
rather than the answer itself:

| Item | Probe | The general rule admitted |
|---|---|---|
| §1 direction of a difference | KB-01 | every difference quantity states which term was subtracted; a term value is measured upward from the ground state |
| §2 adiabatic vs vertical | KB-02 | the distinction is a geometry convention that changes the workflow shape, for any state-to-state quantity |
| §4 established assignments | KB-01 | settled literature facts may be stated as prior knowledge and labelled as such; the generative principles are given, not a list of systems |
| §6 convergence for differences | KB-01 | a difference is more convergence-sensitive than either absolute, so the settings are part of what was computed |

## Items carried without a paired run

- **§3 (which energy terms are included), §5 (reference-method limitation), §7
  (thermochemistry), §8 (units)** have no paired dry run behind them yet. §7 is
  additionally motivated by an observed defect — an executed run silently used
  1 atm for a request that named 1 bar — but an observed defect is not the
  paired evidence this file requires.
- The `KB-03` standard-state case is constructed but unrun.

These sections are definitional rather than behavioural claims, and none of them
asserts an accuracy target. They are recorded here as unadmitted so the
distinction between confirmed and carried knowledge stays visible.
