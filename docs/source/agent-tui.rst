ChemSmart Agent terminal interface
==================================

The terminal interface is the human-facing route from a scientific request to
project YAML, compiled ChemSmart CLI operations, a safe preview, and optional
approved execution.  It does not create a second chemistry engine or ask the
model to write native program input.

Install and open it
-------------------

Install the optional interface dependency alongside ChemSmart::

   python -m pip install 'chemsmart[agent-tui]'

Planning and existing-result analysis need a provider profile, its separate
secret assignment file, and a task workspace::

   chemsmart agent tui \
     --provider PROFILE \
     --provider-config /path/to/agent.yaml \
     --secret-file /path/to/provider.env \
     --workspace /path/to/task

Use ``--plain`` for inline rendering in a simple SSH terminal.  Alibaba Token
Plan and DeepSeek are the registered version-3.1.4 adapters.  The selected
provider, endpoint, model, reasoning setting, and credential label remain user
configuration.

Add an execution envelope when the session may run a chemistry engine::

   chemsmart agent tui \
     --provider PROFILE \
     --provider-config /path/to/agent.yaml \
     --secret-file /path/to/provider.env \
     --workspace /path/to/task \
     --execution-envelope /path/to/resources.yaml

``--review-file`` may additionally export the prepared workflow for an audit.
It is optional and is not an approval token.  A review export requires an
execution envelope because analysis-only and preview-only sessions have no
executable workflow to export.

Plan and safe preview
---------------------

Enter a scientific request at the prompt.  ChemSmart then:

#. preserves the supplied molecular identity, geometry role, charge,
   multiplicity, and coordinate units;
#. creates loader-validated project YAML;
#. compiles the project through the live ``chemsmart run`` surface;
#. generates and inspects program-native artifacts in safe-preview mode; and
#. presents the scientific DAG, effective settings, findings, and typed
   evidence in the transcript.

Planning does not launch a chemistry engine.  Native Gaussian, ORCA, xTB, or
PySCF artifacts remain outputs of ChemSmart rather than a model-authored API.

Approve a displayed workflow
----------------------------

When an execution envelope is present and every executable calculation node
reaches review readiness, the interface shows:

* every planned stage, with any release-unsupported stage marked deferred and
  its reason shown;
* molecule, charge, and multiplicity for every executable node;
* effective project settings;
* the human-readable ChemSmart CLI operation;
* producer-to-consumer data handoffs;
* program and engine observations; and
* core, memory, GPU, timeout, and engine-call limits.

Review that display and enter ``/approve`` once.  The command takes no
arguments.  It removes the workflow from the pending state, disconnects the
provider, and runs the displayed executable partition through the normal
ChemSmart executor.  A deferred stage stays unapproved and is named again in
the execution result.  No hash or approval-file token is required from the
human.  Internal receipts remain available as provenance and mutation
evidence.

The run is written below
``WORKSPACE/.chemsmart-agent/executions/tui-...``.  A failed or completed run
cannot be started again by repeating ``/approve``; create and review a new
workflow instead.

Commands
--------

``/capabilities``
   Show the declared Agent program, engine, job, preview, and execution matrix.
   A declared route is not proof that the selected host has its engine or that
   the route has been release-qualified there.

``/status``
   Show the current phase, workspace, task, workflow awaiting review, and
   whether human approval is pending.

``/approve``
   Approve and run the displayed executable stages once.  It accepts no
   arguments and is available only after their complete safe preview and the
   full-plan human review display.

``/deny``
   Decline the displayed workflow without launching an engine.

``/revise``
   Decline the displayed workflow and return to the prompt for a revised
   scientific request.

``/help``
   Show the in-terminal command guide.

``/quit``
   Exit when no host operation is running.  ``/exit`` is an alias.

Result analysis
---------------

Analysis does not require an execution envelope or approval.  Put supported
completed results in the workspace and ask for the quantities, comparisons,
or interpretation needed.  ChemSmart can register Gaussian and ORCA native
outputs, validated xTB result folders, structured PySCF HDF5 results, and XYZ
geometries or trajectories.  It then exposes semantic quantities and
unit-aware expressions rather than asking the model to copy numbers from a
terminal transcript.

After an approved run completes, enter a new scientific request over the
completed result workspace.  This analysis request may contact the selected
provider, but it does not extend the consumed execution decision or rerun an
engine.

Version-3.1.4 execution boundary
--------------------------------

Release-qualified Agent execution includes PySCF CPU ``sp``, ``opt``, and
``hess``; xTB CPU ``sp``, ``opt``, and ``hess``; and ORCA CPU
single-point, optimization/frequency, transition-state, excited-state, and
serial DAG workflows.  ORCA ``irc`` and ``neb`` remain planning and preview
paths until the selected target is qualified.

Gaussian ``sp``, ``opt``, ``ts``, ``irc``, ``td``, and ``link`` are supported
for project YAML, native-input generation, safe preview, and analysis of
user-supplied completed outputs; this release does not claim Gaussian Agent
execution.  GPU4PySCF ``sp``, ``opt``, and ``hess`` are configuration and
preview paths until a compatible GPU target is qualified.  PySCF CPU ``td`` is
preview-only.

What the interface does not do
------------------------------

The interface does not invent coordinates, native input, server credentials,
or hidden resource defaults.  A provider turn cannot invoke ``/approve``.  The
interface does not submit an unreviewed job, offer an ``always allow`` mode, or
claim that a configured engine is installed.  Program availability and
requested resources must be observed on the selected target before a real
workflow is shown for approval.
