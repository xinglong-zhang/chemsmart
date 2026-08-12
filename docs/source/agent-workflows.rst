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

For a user-authorized local research episode, ``--execution-envelope`` keeps
planning, ChemSmart preview, engine execution, result registration and typed
analysis in one model conversation. The envelope contains operating authority,
not a method or answer DAG:

.. code-block:: yaml

   schema_version: chemsmart.bounded-execution-envelope.v1
   mode: continuous-local
   allowed_program_engines:
     gaussian: [cpu]
     orca: [cpu]
     pyscf: [cpu]
     xtb: [cpu]
   resources:
     execution_target: run
     cores: 8
     memory_gb: 28
     gpu_count: 0
     scratch_policy: server
     node_timeout_seconds: 3600
   episode_wall_time_seconds: 5400
   postprocess_reserve_seconds: 600
   max_engine_calls: 8
   scratch_root: /workspace/chemsmart-bench/scratch

Use it instead of ``--approval-file``:

.. code-block:: bash

   chemsmart agent run \
     --provider PROFILE \
     --task-file task.md \
     --secret-file /secure/path/api.env \
     --workspace ./agent-workspace \
     --identity-manifest identities.yaml \
     --execution-envelope envelope.yaml

Every calculation still passes through the project loader, live Click schema,
native-input writer, preview verifier, program runner and result validator. A
future-geometry edge may consume a validated optimized structure. When the DAG
explicitly declares a different charge and multiplicity for its consumer,
CHEMSMART preserves the coordinates and atom order while rebinding the target
electronic state; this supports four-point reorganization-energy workflows
without pretending that the producer optimized the second state.

Use a disposable workspace containing only approved molecular, project and
result artifacts. Store credentials outside the workspace and repository.

When a task uses one or more named geometries, pass an approved molecular
input manifest with ``--identity-manifest``. Each entry names a
workspace-relative XYZ file and its exact SHA-256, approved molecular names,
geometry role, coordinate units, charge, multiplicity and source locator.
CHEMSMART validates the exact bytes and atom order before publishing the
path-free identity and state records to the model. Filenames and XYZ comments
remain non-authoritative, and the manifest must not contain an expected value,
answer DAG or native program input.

.. code-block:: yaml

   schema_version: chemsmart.approved-molecular-input-manifest.v1
   inputs:
     - input_id: water-initial
       identity_id: water-initial
       approved_names: [water, H2O]
       geometry_file: inputs/water.xyz
       geometry_sha256: <sha256-of-exact-xyz-bytes>
       coordinate_units: angstrom
       geometry_role: neutral optimization start
       charge: 0
       multiplicity: 1
       source_locator: https://example.org/coordinate-record
       source_record_sha256: <sha256-of-source-record>
       state_source_locator: benchmark-case:water#initial-state

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

Current analysis boundary
=========================

Cross-program quasi-harmonic thermochemistry remains available through the
typed analysis operation and the existing ``chemsmart run thermochemistry``
CLI. The CLI command does not yet load project YAML or materialize as a
result-artifact analysis node in a scientific workflow. Therefore the agent
surface has not been removed in this release: doing so would delete the only
model-callable route used by thermochemistry-dependent workflows. The general
follow-up is a canonical YAML-backed analysis command node whose structured
receipt can re-enter dimensional expressions; it is not a benchmark-specific
repair.

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
