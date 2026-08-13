Computational Agent Workflows
#############################

ChemSmart Agent turns a scientific request into project YAML, compiled
ChemSmart CLI operations, a safe preview, and typed result analysis.  The
model proposes the chemistry; ChemSmart owns native-input generation,
execution, validation, and result parsing.

Supported scope
===============

Version 3.1.4 exposes the following product-capable engine/job pairs through
project YAML, command compilation, safe preview, exact human approval, and the
normal ChemSmart runners:

* Gaussian CPU: ``irc``, ``link``, ``opt``, ``sp``, ``td``, and ``ts``.
* ORCA CPU: ``irc``, ``neb``, ``opt``, ``sp``, ``td``, and ``ts``.
* PySCF CPU: ``hess``, ``opt``, and ``sp``; CPU ``td`` is preview-only.
* GPU4PySCF: PySCF ``hess``, ``opt``, and ``sp`` with the ``gpu`` engine.
* xTB CPU: ``hess``, ``opt``, and ``sp``.

These are product capabilities, not claims that an executable is available on
the current machine.  Each approved run still requires the environment probe
for its exact program and engine to report available.  Gaussian requires a
separately licensed installation and GPU4PySCF requires a compatible CUDA
stack. Typed analysis is available only for quantities actually parsed and
registered from a completed result. NCIPLOT and other CLI leaves without an
Agent engine/job declaration remain outside this surface.

Configure a provider
====================

The Runtime is provider-neutral.  Version 3.1.4 registers Alibaba Token Plan
and DeepSeek adapters.  It has no default model: every selected profile must
name its model explicitly.  Credentials belong in a separate assignment file,
not in ``agent.yaml``.

The following example defines both registered adapters.  Select either profile
with ``active`` or ``--provider`` and replace the model placeholders with model
IDs available to your accounts.

.. code-block:: yaml

   active: alibaba
   fallback: []
   providers:
     alibaba:
       type: openai
       api_key_env: ALIBABA_TOKEN_PLAN_KEY
       model: REPLACE_WITH_YOUR_ALIBABA_MODEL_ID
       context_tokens: REPLACE_WITH_MODEL_CONTEXT_LIMIT
       max_output_tokens: REPLACE_WITH_MODEL_OUTPUT_LIMIT
       base_url: https://token-plan.ap-southeast-1.maas.aliyuncs.com/compatible-mode/v1
       reasoning_effort: max
       preserve_thinking: true
       transport_deadlines:
         connect_seconds: 15
         first_event_seconds: 90
         inter_event_seconds: 90
         absolute_turn_seconds: 300
     deepseek:
       type: openai
       api_key_env: DEEPSEEK_API_KEY
       model: REPLACE_WITH_YOUR_DEEPSEEK_MODEL_ID
       context_tokens: REPLACE_WITH_MODEL_CONTEXT_LIMIT
       max_output_tokens: REPLACE_WITH_MODEL_OUTPUT_LIMIT
       base_url: https://api.deepseek.com
       reasoning_effort: max
       preserve_thinking: true

Replace every model and limit placeholder with values supported by the
selected account and model. Both token-limit fields are required positive
integers, and
``max_output_tokens`` cannot exceed ``context_tokens``.  They are explicit
profile facts: ChemSmart binds them into the profile digest, Runtime
capabilities, and the provider request.  For ``qwen3.8-max``, use
``reasoning_effort: xhigh``.

``fallback`` must be empty.  Runtime does not switch providers inside a
session, so a non-empty fallback list is rejected.  After a classified
provider failure, start a new, explicitly attributed attempt with
``--provider OTHER_PROFILE``.  Private reasoning, response text, and
credentials are not written to the public session record.

The assignment file contains only the label selected by the active profile:

.. code-block:: text

   ALIBABA_TOKEN_PLAN_KEY=your-secret-value

Plan and safe preview
=====================

Use a workspace containing only user-approved molecular, project, and result
artifacts.  The command below contacts the selected provider but never starts
a chemistry engine:

.. code-block:: bash

   chemsmart agent plan \
     --provider PROFILE \
     --provider-config /secure/path/agent.yaml \
     --task-file task.md \
     --secret-file /secure/path/provider.env \
     --workspace /absolute/path/task-workspace

``chemsmart agent run`` has the same preview-only authority.  It is retained
for workflows that also request an execution review; passing
``--approval-file`` to a provider session fails closed.

Exact execution approval
========================

Real execution is a separate, three-command chain.  First provide explicit
program and resource bounds.  These bounds are review input, not permission.

.. code-block:: yaml

   schema_version: chemsmart.bounded-execution-envelope.v1
   mode: bounded-local
   allowed_program_engines:
     orca: [cpu]
     pyscf: [cpu]
     xtb: [cpu]
   resources:
     execution_target: run
     cores: 4
     memory_gb: 8
     gpu_count: 0
     scratch_policy: server
     node_timeout_seconds: 1800
   episode_wall_time_seconds: 3600
   postprocess_reserve_seconds: 300
   max_engine_calls: 4
   scratch_root: /absolute/path/chosen-by-the-user

Create an inert review packet after project compilation and safe preview:

.. code-block:: bash

   chemsmart agent plan \
     --provider PROFILE \
     --provider-config /secure/path/agent.yaml \
     --task-file task.md \
     --secret-file /secure/path/provider.env \
     --workspace /absolute/path/task-workspace \
     --execution-envelope /absolute/path/resources.yaml \
     --review-file /absolute/path/execution-review.json

Read the complete review.  It binds the molecular and electronic state,
project and input bytes, scientific DAG, compiled operations, data handoffs,
resources, and observed program environments.  Then record exactly one human
decision using the displayed ``review_sha256``:

.. code-block:: bash

   chemsmart agent approve \
     --review-file /absolute/path/execution-review.json \
     --reviewed-sha256 FULL_REVIEW_SHA256 \
     --decision approve \
     --actor HUMAN_ID \
     --output /absolute/path/execution-approval.json

``--decision`` also accepts ``deny``, ``revise``, or ``quit``; those decisions
create no execution grant.  There is no permanent, session-wide, or
command-prefix approval.  A change to any reviewed input requires a new review
and approval.

Finally execute the one-shot bundle without a model or provider credential:

.. code-block:: bash

   chemsmart agent execute \
     --approval-file /absolute/path/execution-approval.json \
     --workspace /absolute/path/task-workspace \
     --run-directory /absolute/path/run

The workspace keeps a durable consumption ledger, so the same bundle cannot
be executed twice.  If pre-launch validation or execution stops, review and
approve the corrected workflow again.

Terminal interface
==================

Install the optional interface and open the same approval chain in a terminal:

.. code-block:: bash

   python -m pip install 'chemsmart[agent-tui]'
   chemsmart agent tui \
     --provider PROFILE \
     --provider-config /secure/path/agent.yaml \
     --secret-file /secure/path/provider.env \
     --workspace /absolute/path/task-workspace \
     --execution-envelope /absolute/path/resources.yaml \
     --review-file /absolute/path/execution-review.json

The TUI displays the exact review and offers explicit approve, deny, revise,
and quit-review actions.  Execution still requires re-entering the full
approval-file SHA-256.  See :doc:`agent-tui` for the command reference.

Molecular identity
==================

When a task uses named geometries, pass ``--identity-manifest``.  Every entry
binds a workspace-relative XYZ file to its exact SHA-256, approved name,
geometry role, coordinate units, charge, multiplicity, and source.  Filenames
and XYZ comments are not identity evidence.

.. code-block:: yaml

   schema_version: chemsmart.approved-molecular-input-manifest.v1
   inputs:
     - input_id: water-initial
       identity_id: water-initial
       approved_names: [water, H2O]
       geometry_file: inputs/water.xyz
       geometry_sha256: SHA256_OF_EXACT_XYZ_BYTES
       coordinate_units: angstrom
       geometry_role: neutral optimisation start
       charge: 0
       multiplicity: 1
       source_locator: SOURCE_RECORD
       source_record_sha256: SHA256_OF_SOURCE_RECORD
       state_source_locator: SOURCE_RECORD_STATE

Provider waits and deadlines
============================

Provider profiles may bind connect, first-event, inter-event, and absolute-turn
deadlines as shown above.  The absolute deadline cannot be extended by stream
heartbeats or partial bytes.  Synchronous platform DNS resolution is the
documented exception and must be checked on the target host.

During an approved program run, the Runtime makes no provider call while the
normal ChemSmart runner waits synchronously.  Public wait and wake events
record that interval; they are not a scheduler, sleep service, or parallel
agent.  Process timeout and termination remain the runner's responsibility.

Interpreting evidence
=====================

Keep the following states distinct in reports: proposed, planned,
materialised, previewed, approved, executing, engine-complete, parsed,
scientifically validated, and interpreted.  A provider answer, fake preview,
fixture parser, or engine return code alone is not a scientific result.

Check molecular identity and atom order, charge and multiplicity, method and
program semantics, convergence, stationary-point evidence, geometry handoff,
physical conditions, signs, dimensions, and units before using a value in a
scientific claim.

Analyse completed results
=========================

Analysis is a normal Agent mode and does not require execution approval. Put
one or more supported completed results in an otherwise clean workspace and
ask for the quantities or comparison you need. ChemSmart discovers normally
terminated Gaussian and ORCA outputs, receipt-bound xTB outputs, structured
PySCF HDF5 results, and XYZ geometry/trajectory artifacts. For example::

   chemsmart agent plan \
     --provider PROFILE \
     --provider-config /secure/path/agent.yaml \
     --secret-file /secure/path/provider.env \
     --workspace /absolute/path/completed-results \
     --task "Extract the final energy and frequencies, diagnose the stationary point, and report the result with explicit units."

The Agent may inspect registered artifacts, extract semantic quantities,
derive RRHO or explicitly parameterised quasi-harmonic thermochemistry,
evaluate unit-aware expression DAGs (including energy differences, CBS
extrapolation, Boltzmann populations, distances, angles and dihedrals), and
bind numerical claims to the resulting typed evidence. Supply
``--analysis-completion-file`` only when an application needs a stricter,
host-authored list of mandatory quantities or claims; ordinary existing-result
analysis does not require one.

After ``agent execute`` completes, start this analysis mode over the completed
result workspace. The deterministic executor itself stops at validated engine
receipts and never contacts a provider; scientific interpretation is a new,
explicit Agent request rather than an implicit post-run model call.
