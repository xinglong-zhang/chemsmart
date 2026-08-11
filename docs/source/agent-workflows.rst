Computational Agent Workflows
#############################

CHEMSMART's agent is a provider-neutral computational-chemistry researcher that
operates through the same project YAML, Click commands, program adapters,
runners and result readers as a human user.

The model proposes scientific intent. CHEMSMART owns native input generation,
command construction, program execution and structured analysis. A model must
not repair a workflow by writing native Gaussian, ORCA, xTB or PySCF execution
text.

Agent commands
==============

Plan and safely preview a workflow:

.. code-block:: bash

   chemsmart agent plan \
     --provider PROFILE \
     --task-file task.md \
     --secret-file /secure/path/api.env \
     --workspace ./agent-workspace

Plan and optionally execute reviewed nodes:

.. code-block:: bash

   chemsmart agent run \
     --provider PROFILE \
     --task-file task.md \
     --secret-file /secure/path/api.env \
     --workspace ./agent-workspace \
     --approval-file workflow-approval.json

Without ``--approval-file``, ``agent run`` remains preview-only. Deterministic
execution of an already approved workflow uses ``chemsmart agent execute`` and
requires no provider credential.

Use a disposable workspace containing only approved molecular, project and
result artifacts. Store credentials outside the workspace and repository.

Challenge-driven scientist persona
===================================

The agent should behave as a senior computational scientist:

- ask for scientifically consequential missing facts;
- choose defensible programs and methods rather than mimic a hidden answer;
- preserve molecular identity, geometry, state, conditions and units;
- design causal calculation and post-processing DAGs;
- react to real engine and parser failures by revising the scientific plan;
- distinguish planned, previewed, executed, analysed and inferred work; and
- explain limitations in proportion to their effect on the requested claim.

Creative routes are welcome. A different algebraic decomposition, a
program-native alternative, or a composite method may be superior to the path
anticipated by the task designer.

Self-improving research cycle
=============================

Self-improvement is driven by observation:

1. give the agent a concise, difficult and scientifically complete task;
2. run the real model through the normal CHEMSMART surface;
3. read its answer and tool trajectory as a computational scientist;
4. classify the consequential obstacle as chemistry, model reasoning, program
   capability, environment, parser or tool affordance;
5. improve the smallest general CHEMSMART layer or project skill;
6. rerun a relevant task and compare the scientific behavior; and
7. retain the lesson only when it generalises beyond one paper or molecule.

Self-improvement does not mean autonomous privilege escalation, hidden prompt
mutation or learning benchmark answers. Repository edits, paid provider calls,
engine execution and scheduler use remain within user authority.

Human scientific evaluation
===========================

Read the visible answer directly. Judge molecular and electronic-state
preservation, method suitability, numerical transformations, signs, units,
physical conventions and necessary causal dependencies. Accept any route that
is chemically and mathematically sound.

Use these outcomes:

- **Improvement:** correct requested result or equivalent coherent plan.
- **Partial improvement:** central chemistry improved, but an important
  quantity, condition, dependency or limitation remains unresolved.
- **No improvement:** invented data, altered scientific problem, invalid
  chemistry or mathematics, or an unperformed action presented as executed.

Tool-call count, token use, exact DAG shape and runtime labels are operating
observations, not scientific oracles.

Project skills
==============

The maintained project-local skills are under ``.agents/skills``:

- ``chemsmart-agent-harness``;
- ``chemsmart-agent-harness-refine-loop``;
- ``chemsmart-cli-program-hub``;
- ``chemsmart-scientific-workflow``;
- ``chemsmart-evidence-audit``;
- ``pyscf-v2-14-0``; and
- ``gpu4pyscf-v1-8-0``.

They are operational guidance, not a substitute for the checked-out source or
the live CLI. Update them when a real general lesson makes their guidance
wrong, incomplete or suppressive.
