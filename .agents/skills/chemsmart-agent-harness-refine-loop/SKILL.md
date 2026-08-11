---
name: chemsmart-agent-harness-refine-loop
description: Run bounded, behavior-first, self-improving experiments on the ChemSmart computational-chemistry agent using difficult scientific tasks, real model calls, direct expert reading, and the smallest general CLI, YAML, parser, DAG, or analysis improvement. Use for iterative agent research without answer-specific hard-coding or deterministic scientific grading.
---

# ChemSmart Agent Harness Refine Loop

## Mission

Act as a challenging senior computational scientist and automation-agent
researcher. Study how the model actually solves computational chemistry through
ChemSmart. Prefer an informative failed experiment to an unrun plan, and a
sharper experiment to another speculative gate.

ChemSmart project YAML, CLI subcommands, program adapters, parsers, scientific
DAGs and typed analysis form the control layer. Keep the model free to discover
a chemically valid route that differs from yours.

## Iteration bound

Require one positive integer iteration count:

```bash
python scripts/validate_loop_request.py --iterations 1
```

Stop after that many completed observe-change-observe cycles. This is the only
deterministic parser in the skill.

## Self-improving experiment cycle

1. **Challenge:** choose a difficult, source-complete problem that fixes only
   chemically necessary identity, state, conditions, artifacts and requested
   observable.
2. **Compress:** write a short model-visible task. Do not reveal an answer,
   preferred DAG, prior failure or repair hint.
3. **Predict:** privately state competing hypotheses about how the current
   agent may succeed or fail.
4. **Observe:** run the real authorised model through the ordinary ChemSmart
   surface before adding a rule or test.
5. **Reconstruct:** read the answer and tool trajectory directly. Understand
   the route the model discovered before comparing it with your own.
6. **Diagnose:** distinguish chemistry, mathematics, model reasoning, program
   capability, environment, parser and missing affordance.
7. **Improve:** change the smallest general ChemSmart layer or project skill.
   Never encode a paper, molecule, DOI, expected number or private route.
8. **Re-observe:** run a fresh real task and judge the scientific behavior.
9. **Internalize:** retain the general lesson in source, documentation or the
   appropriate skill; remove guidance disproved by the observation.
10. **Challenge transfer:** when overfitting is plausible, use a materially
    different chemistry case. Transfer is evidence, not paperwork.

The loop improves itself by revising its hypotheses and guidance from real
observations. It must not autonomously edit outside the user's scope, increase
its privileges, run paid resources without authority or manufacture a success
from tests.

## Human scientific evaluation

Do not use code to grade the agent answer. Read it as a computational scientist
and assess:

- molecular identity, atom order, geometry, charge, multiplicity and state;
- method and program suitability;
- quantities, signs, dimensions, units and physical conventions;
- necessary scientific causality in the DAG;
- distinction between archived evidence, inference, plan, preview, analysis
  and execution; and
- whether limitations matter to the requested claim.

Accept algebraically equivalent post-processing, different valid DAG
decompositions and program-native alternatives. Tool count, tokens, exact graph
shape, hashes and runtime labels are not scientific oracles.

Use **Improvement**, **Partial improvement**, or **No improvement** as defined
in `AGENTS.md`. The main scientist owns this judgment.

## Reporting

Report the discovered route, strong scientific decisions, consequential
errors, general capability learned, and improvement disposition. Keep
mechanical checks secondary and state exactly what was previewed, executed,
analysed or inferred.
