Computational Agent Workflows
#############################

ChemSmart Agent turns a scientific request into project YAML, compiled
ChemSmart CLI operations, a safe preview, and typed result analysis.  The
model proposes the chemistry; ChemSmart owns native-input generation,
execution, validation, and result parsing.

Supported scope
===============

The current release provides broad planning, project-YAML generation, live CLI
compilation, and safe preview for these Agent job families:

* Gaussian CPU: ``irc``, ``link``, ``modred``, ``opt``, ``scan``, ``sp``,
  ``td``, and ``ts``.
* ORCA CPU: ``irc``, ``modred``, ``neb``, ``opt``, ``scan``, ``sp``, ``td``,
  and ``ts``.
* PySCF CPU: ``hess``, ``opt``, ``sp``, and preview-only ``td``.
* GPU4PySCF: PySCF ``hess``, ``opt``, and ``sp`` with the ``gpu`` engine.
* xTB CPU: ``hess``, ``opt``, and ``sp``.

Release-qualified real Agent execution is narrower:

* PySCF CPU ``sp``, ``opt``, and ``hess``;
* xTB CPU ``sp``, ``opt``, and ``hess``; and
* ORCA CPU single-points, optimization/frequency, transition-state,
  excited-state, relaxed coordinate scans, and serial producer-to-consumer
  DAGs.

A relaxed scan states its driven coordinate on the workflow node rather than in
project YAML: which atoms, which coordinate type, and the range and number of
points.  The same specification is rendered into each program's own idiom, so
one physical description reaches ORCA's absolute endpoints and Gaussian's
increment without the caller writing either.  A completed scan's surface is
read into typed quantities through the ordinary analysis layer.

The first constrained optimisation of a scan imposes the driven coordinate on
the geometry supplied, so a range beginning far from that geometry's current
value may be refused by the program before any optimisation runs.

ORCA ``irc``, ``neb``, and ``modred`` remain planning and preview paths until
the selected target is qualified.  Gaussian Agent execution is not claimed in
this release; Gaussian support covers project YAML, generated native input,
safe preview, and typed analysis of user-supplied completed outputs.  GPU4PySCF
remains a configuration and preview path until a compatible GPU target is
qualified.

These boundaries do not alter the wider human ``chemsmart run`` and
``chemsmart sub`` CLI.  They also do not imply that an executable is available
on the current machine.  Every approved CPU run still needs an observed program
environment and sufficient target resources.  Gaussian requires a separately
licensed installation and GPU4PySCF requires a compatible CUDA stack.

Configure a provider
====================

The runtime is provider-neutral.  The current release registers Alibaba Token Plan
and DeepSeek adapters.  It has no default model: every selected profile names
its model, endpoint, context and output limits, reasoning setting, and
credential label explicitly.  Credentials belong in a separate assignment
file rather than ``agent.yaml``.

The following example defines one profile for each registered adapter.  Replace
the placeholders with values supported by the selected account and model.

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
       reasoning_effort: REPLACE_WITH_SUPPORTED_REASONING_VALUE
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
       reasoning_effort: REPLACE_WITH_SUPPORTED_REASONING_VALUE
       preserve_thinking: true

Both token limits are required positive integers and
``max_output_tokens`` cannot exceed ``context_tokens``.  ``fallback`` must be
empty because ChemSmart does not switch providers inside one session.  After a
provider failure, start a new, explicitly attributed attempt with another
profile.

The separate assignment file contains the label selected by the active
profile::

   ALIBABA_TOKEN_PLAN_KEY=your-secret-value

ChemSmart parses this file as data; it does not source it into the shell.
Provider reasoning, response text, and credentials are not scientific
execution evidence.

Plan and safe preview
=====================

Use a workspace containing the molecular, project, and completed-result
artifacts relevant to the task.  This command contacts the selected provider
but cannot start a chemistry engine:

.. code-block:: bash

   chemsmart agent plan \
     --provider PROFILE \
     --provider-config /secure/path/agent.yaml \
     --task-file task.md \
     --secret-file /secure/path/provider.env \
     --workspace /absolute/path/task-workspace

``chemsmart agent plan`` may create and validate project YAML, compile live
Click commands, inspect generated artifacts, and execute fake previews.  It
never accepts a model-created execution grant; real execution belongs to the
approval chain below.

Approved execution in the terminal interface
============================================

Real execution starts from explicit program and resource limits.  The envelope
below permits up to four serial CPU engine calls on the listed programs.  It is
a bound on a possible run, not permission to run one.

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
     memory_gb: 4
     gpu_count: 0
     scratch_policy: server
     node_timeout_seconds: 1800
   episode_wall_time_seconds: 3600
   postprocess_reserve_seconds: 300
   max_engine_calls: 4
   scratch_root: /absolute/path/chosen-by-the-user

The envelope supplies the explicit per-run cores, memory, and GPU allocation,
overriding ordinary server defaults.  The human review displays that allocation;
the selected program, operating system, and any external scheduler must still be
able to provide it.  Open the terminal interface with the envelope:

.. code-block:: bash

   chemsmart agent tui \
     --provider PROFILE \
     --provider-config /secure/path/agent.yaml \
     --secret-file /secure/path/provider.env \
     --workspace /absolute/path/task-workspace \
     --execution-envelope /absolute/path/resources.yaml

After planning and safe preview, the interface displays every planned stage,
including a release-unsupported stage that remains scientifically necessary.
It marks that stage deferred, gives its reason, and displays the molecule and
state, effective project setting, ChemSmart CLI operation, data handoff,
program environment, and resource bound for every executable stage.  Enter
``/approve`` once to execute only those reviewed executable stages; deferred
stages remain unapproved and are not launched.  Use ``/deny`` or ``/revise``
without launching an engine.  The provider is disconnected before execution.

The human does not retype a hash or create an approval-file token.  Internal
receipts and content digests remain provenance.  The pending workflow is
consumed before launch, so a failed or completed run requires a fresh plan and
human review before another execution attempt.

Use ``--review-file /absolute/path/review.json`` only when an audit application
needs an exported copy of the displayed workflow.  The file is optional and is
not the execution authority for the terminal interface.

Molecular identity
==================

When a task uses named geometries, ``--identity-manifest`` can bind each
workspace-relative XYZ file to its approved name, geometry role, coordinate
units, charge, multiplicity, and source.  Filenames and XYZ comments alone are
not molecular-state evidence.

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

The content fields in this optional manifest preserve input provenance.  They
are separate from the human ``/approve`` action and are never retyped in the
terminal interface.

Causal data handoffs
====================

Downstream calculations consume typed artifacts rather than guessed filenames.
An optimization may provide its validated final geometry to a Hessian,
single-point, or excited-state node while preserving atom order, charge, and
multiplicity.  A project that already requests frequencies produces them on
that same scientific node; a duplicate Hessian is not required unless the
scientific protocol asks for an independent calculation.

For an ORCA transition-state-to-IRC workflow, the IRC node has two scientific
inputs: the final TS geometry and the final ORCA Hessian.  ChemSmart selects the
unique Hessian whose atom order, geometry, frequencies, and single
consequential imaginary mode agree with the validated TS output.  This
materialisation capability does not by itself claim a completed IRC run on an
unqualified target.

Provider waits and deadlines
============================

Provider profiles may set connect, first-event, inter-event, and absolute-turn
deadlines.  The absolute deadline is not extended by stream heartbeats or
partial bytes.  During an approved chemistry run the executor makes no provider
call; process timeout and termination remain the ChemSmart runner's
responsibility.

Interpreting evidence
=====================

Keep these states distinct: proposed, planned, materialised, previewed,
approved, executing, engine-complete, parsed, scientifically validated, and
interpreted.  Provider prose, a fake preview, a parser example, or a successful
process exit alone is not a scientific result.

Before using a value, check molecular identity and atom order, charge and
multiplicity, method and program semantics, convergence, stationary-point
evidence, geometry handoff, physical conditions, signs, dimensions, and units.

When a node fails
=================

A node whose engine run did not succeed reports its terminal state, the wrapper
and child exit statuses, the validator findings, and a bounded quotation of what
the program itself printed about the failure.

That quotation is the program's own text, not a ChemSmart claim.  It is evidence
that a diagnosis exists and what it said; it never establishes readiness,
validity, or what to do next.  URLs, absolute and home-relative paths, e-mail
addresses, and credential-like assignments are removed from it.

A run that did not terminate normally still yields no scientific quantities.
Typed extraction continues to require a normally terminated result, so a failed
run can be read for its reason and not for its numbers.

Re-running an approved workflow
===============================

``chemsmart agent review`` re-presents a stored execution review so the same
workflow can be decided on again::

   chemsmart agent review \
     --review-file review.json \
     --workspace /path/to/workspace

With no ``--decision`` it reports only: the review digest, the workflow, whether
every approved artifact is still present under the workspace, and which approval
identities for this review have already been consumed.  Adding
``--decision approve --actor NAME`` records a new human decision and writes a new
one-shot bundle.

This does not reuse a spent approval.  Approvals remain one-shot and bound to the
exact request digest; replay obtains the *current* decision that the approval
chain requires for a launch, over an unchanged displayed plan.

Differences between the approval and the present workspace are displayed rather
than silently accepted, and they change no enforcement: the environment and
command comparison that runs immediately before the first dispatch still refuses
a launch whose facts have drifted from the reviewed ones.  Approved input bytes
that are no longer present under the workspace are refused before a decision is
offered, because resolving them would fail before anything could run.

Domain-knowledge skills
=======================

The Agent can consult short, advisory skill documents that describe general
computational-chemistry practice.  They are surfaced by name in the planning
prompt and fetched on request; they are advisory only and never establish
readiness, approval, terminal state, or an accuracy claim, and they never
replace a typed host receipt.

Released skills:

* ``scientific-conventions`` — how computed quantities are conventionally
  defined and reported: the direction of every difference quantity, adiabatic
  versus vertical geometry, which energy terms are included, and thermochemistry
  standard states.
* ``method-adequacy`` — whether a chosen method, basis set, solvation model or
  conformer sample can answer the question being asked: which errors cancel in a
  comparison and which do not, and how to state the resulting uncertainty.
* ``typed-analysis-contract`` — how the typed analysis layer expects intent to be
  expressed: identifiers, units, declared quantity kinds, and evidence
  references.

Set ``CHEMSMART_SKILL_ROOT`` to add or override skills from a directory of your
own; set ``CHEMSMART_AGENT_SKILLS=0`` to remove the skill index and the
consultation tool entirely.

Analyse completed results
=========================

Analysis is a normal Agent mode and does not require an execution envelope or
approval.  Put supported completed results in a task workspace and request the
quantities or comparison needed.  ChemSmart discovers:

* normally terminated Gaussian and ORCA native outputs;
* validated xTB result folders, including relocated archives whose missing
  original paths are reported as provenance limitations;
* structured PySCF HDF5 results; and
* XYZ geometries and trajectories.

For example::

   chemsmart agent plan \
     --provider PROFILE \
     --provider-config /secure/path/agent.yaml \
     --secret-file /secure/path/provider.env \
     --workspace /absolute/path/completed-results \
     --task "Extract the final energy and frequencies, diagnose the stationary point, and report the result with explicit units."

The Agent may extract energies, structures, frequencies, dipoles, excited
states, spin evidence, solvent treatment, and auxiliary-basis roles.  Text
metadata remains typed text and is not assigned a fictitious physical unit.
It may derive RRHO or explicitly parameterised quasi-harmonic thermochemistry
and evaluate unit-aware expression DAGs for, among other operations:

* energy differences and CBS extrapolation;
* Boltzmann populations or averages from one energy vector or separate scalar
  state energies, with optional degeneracies;
* harmonic ZPE and imaginary-mode counts;
* distances, angles, dihedrals, centres of mass, inertia, rotational constants,
  and connectivity changes; and
* photon wavelengths and dimensional unit conversions.

For a composite high-level electronic energy plus lower-level quasi-harmonic
thermochemistry, use
``quasi_harmonic_thermal_gibbs_correction = G_qh(T) - E_electronic``.  Do not
substitute the harmonic ``thermal_gibbs_correction`` while describing the
result as quasi-harmonic.

Supply ``--analysis-completion-file`` only when another application needs a
host-authored list of mandatory quantities or claims.  Ordinary
existing-result analysis does not require it.

A typed analysis chain planned with a workflow travels verbatim in the review
packet and the approval bundle, and the single approval covers it.  After
every approved calculation node validates, ``agent run`` runs the chain
provider-free -- extraction, thermochemistry, expressions, validation
verdicts, and claim rendering -- and writes ``analysis/
completed-analysis-report.md`` into the run directory.  A failed validation
verdict is a completed determination, not an execution failure.  Scientific
interpretation and the recorded decision remain a session act: start a new
analysis request over the completed result workspace to record them; this
does not rerun the engine or extend the earlier execution decision.  A
workflow approved without an analysis chain executes exactly as before.

Two further affordances round out multi-stage work.
``compose_molecular_arrangement`` places two identity-bound geometries into
one arrangement at an explicit atomic contact; the host owns the placement
mathematics and the composed bytes with full parent lineage, the electronic
state is bound explicitly afterwards, and the consuming stage is a new
workflow.  A validated frequency-bearing ORCA producer may feed its Hessian
to an ORCA transition-state search through the live ``--inhess-filename``
option (selection rule ``validated_producer_orca_hessian``); the starting
Hessian may carry any imaginary-mode count, and the observed count is
recorded in the handoff receipt.
