ChemSmart Agent terminal interface
==================================

The terminal interface is the human-facing route from a scientific request to
project YAML, compiled ChemSmart CLI operations, a safe preview, one explicit
approval, provider-free execution, and the completed-analysis report.  It is a
view and controller over the ordinary approval chain -- never a second
permission engine -- and it does not create a second chemistry engine or ask
the model to write native program input.

Install and open it
-------------------

Install the optional interface dependency alongside ChemSmart::

   python -m pip install 'chemsmart[agent-tui]'

The interface is validated on Textual 8.2.  The minimal invocation is two
commands::

   cd /path/to/task
   chemsmart agent

The workspace defaults to the current directory; credentials resolve from
an exported environment variable or the managed key store that
``chemsmart config agent`` maintains; a
workspace's ``execution-envelope.yaml`` and ``identity-manifest.yaml`` are
discovered by convention.  Every default and discovery is announced in the
startup banner -- nothing is silent.  Each may also be given explicitly::

   chemsmart agent tui \
     --provider PROFILE \
     --provider-config /path/to/agent.yaml \
     --workspace /path/to/task

Use ``--plain`` for inline rendering in a simple SSH terminal; it disables the
mouse, replaces the ChemSmart theme with the terminal's own ANSI palette, and
renders the same panels inline.  The selected provider, endpoint, model,
reasoning setting, and credential label remain user configuration; the
interface imposes no default model.

Add an execution envelope when the session may run a chemistry engine::

   chemsmart agent tui \
     ... \
     --execution-envelope /path/to/resources.yaml

``--review-file`` may additionally export the prepared workflow for an audit.
It is optional and is not an approval token.

Watch a plan happen
-------------------

Enter a scientific request at the composer.  From the moment the session's
append-only event stream exists, the interface tails it and every tool call
becomes one humanized line: a pending phrase while the call runs, a
one-sentence summary once it settles (finished lines go muted; a refusal shows
its public message in red).  The project settings are repeated as readable
YAML at the end of the session; the complete canonical payload of every call
remains in the run's ``events.jsonl``.

The interface displays no digests.  Content hashes are provenance, preserved
byte-for-byte in the durable records (event streams, review files, bundles,
and the report file) and behind the ``--json`` flags of the headless
commands; a human panel shows molecules, units, conditions, and words.

Planning can be cancelled: the first ``esc`` arms the interrupt (the footer
says so and disarms after five seconds), the second stops the session at the
next provider-turn boundary.  A cancelled session is still an observed,
well-formed run -- its transcript and event stream are complete -- and it
never leaves a workflow awaiting approval.  A running chemistry engine is
never killable from the interface.

Review what the one approval covers
-----------------------------------

When an execution envelope is present and every executable calculation node
reaches review readiness, the interface displays the complete plan:

* every planned stage, with any release-unsupported stage marked deferred and
  its reason shown;
* molecule, charge, and multiplicity for every executable node;
* the full lineage of a composed molecular arrangement -- both parent
  structures by name, the chosen contact atoms, and the requested versus
  achieved distance;
* the typed analysis chain planned with the workflow (per-node kind, inputs,
  outputs, conditions, and blocked reasons) -- or an explicit statement that
  none is planned;
* per-node project settings as readable YAML (the canonical fallback when a
  node's YAML is unresolvable) and the human-readable ChemSmart CLI
  operation;
* producer-to-consumer data handoffs, program and engine observations, and
  core, memory, GPU, timeout, and engine-call limits; and
* the decision panel.  The full review record is kept in the run evidence;
  no digest is displayed and nothing is ever retyped by the human.

Review that display and enter ``/approve`` once.  It removes the workflow from
the pending state, disconnects the provider, writes a mermaid projection of
the approved DAG (``workflow.mmd``) into the run directory, and runs the
displayed executable partition through the normal ChemSmart executor.
``/deny`` declines; ``/revise`` declines and hands the request back into the
composer for editing.  A deferred stage stays unapproved and is named again in
the execution result.  A failed or completed run cannot be started again by
repeating ``/approve``; create and review a new workflow instead.

Watch the calculation run
-------------------------

While the approved workflow executes, a jobs panel sits between the transcript
and the composer showing what is really running: the node, its program and
molecule, the exact ``chemsmart run`` command (receipt bookkeeping scrubbed),
elapsed wall-clock time, and the queue.  ``/dag`` toggles the workflow panel,
where every node -- calculation and analysis alike -- wears its live status::

   [ ] queued   [•] running   [»] validating   [✓] validated
   [✗] failed   ⊘ deferred (with its reason)

Both panels are fed by tailing the run's own append-only ``events.jsonl`` --
the same evidence stream an audit reads.

When the run finishes, the interface renders the per-node execution table,
the settled analysis chain, and the analysis results -- conditions, claims
with values and units, and any decision prose -- projected from the report
file, whose full provenance stays on disk untouched.  A typed analysis chain approved with the workflow executes
provider-free in the same run; scientific interpretation and the recorded
decision remain a subsequent explicit session act.

Commands
--------

``/approve``
   Approve and run the displayed executable stages once.  Available only
   after their complete safe preview and the full-plan review display.

``/deny`` / ``/revise``
   Decline the displayed workflow.  ``/revise`` also returns the request to
   the composer for editing.

``/status``
   Phase, workspace, task, workflow awaiting review, and pending authority.

``/capabilities``
   The declared Agent program/engine/job matrix.  A declared route is not
   proof that this host has its engine.

``/skills`` and ``/skill <id>``
   List the consultable domain skills; tag the next request with one.  The
   tag is a visible line in the request -- the session still consults the
   skill through its own tool, so the receipt chain is preserved.

``/dag``
   Toggle the workflow panel (works before approval too).

``/resume``
   Restore this workspace's previous session and any pending review (see
   below).

``/report [n]`` and ``/runs``
   Open the completed-analysis report -- the latest, or run ``n`` from the
   ``/runs`` listing of this workspace's executions and replays.

``/export``
   Save the transcript to a file in the workspace.

``/help`` and ``/quit``
   The generated command guide; exit when no host operation is running
   (``/exit`` is an alias).

``ctrl+p`` opens a command palette over the same registry; ``pageup`` /
``pagedown`` scroll the transcript; ``ctrl+c`` quits safely.  A mistyped
slash command answers with the nearest real one.

Resume a workspace
------------------

``chemsmart agent resume`` (or ``/resume`` in a running terminal) restores a
workspace from its durable records: the last session's task and how it
ended, every run with its workflow and report, and the previous session's
closing prose.  A workflow that was reviewed but never decided is
re-presented for one fresh decision -- refused with the reason when the
files it was approved against have moved.  Every terminal approval is
durable (a decision record, a one-shot bundle, and a consumption ledger in
the workspace), so re-approving the same reviewed workflow is refused even
across restarts; ``chemsmart agent review`` records a deliberate second
decision.  A review that was approved but whose run never started is not
re-decided: the greeting prints the exact provider-free launch its one-shot
approval still awaits.

Result analysis
---------------

Analysis of existing completed results does not require an execution envelope
or approval.  Put supported completed results in the workspace and ask for
the quantities, comparisons, or interpretation needed.  ChemSmart registers
Gaussian and ORCA native outputs, validated xTB result folders, structured
PySCF HDF5 results, and XYZ geometries or trajectories, exposing semantic
quantities and unit-aware expressions rather than asking the model to copy
numbers from a terminal transcript.

Execution boundary
------------------

Release-qualified Agent execution includes PySCF CPU ``sp``, ``opt``, and
``hess``; xTB CPU ``sp``, ``opt``, and ``hess``; and ORCA CPU single-point,
optimization/frequency, transition-state, excited-state, relaxed coordinate
scans, and serial DAG workflows.  ORCA ``irc``, ``neb``, and constrained
optimisation (``modred``) remain planning and preview paths until the
selected target is qualified.

Gaussian ``sp``, ``opt``, ``ts``, ``irc``, ``td``, ``link``, ``scan``, and
``modred`` are supported for project YAML, native-input generation, safe
preview, and analysis of user-supplied completed outputs; this release does
not claim Gaussian Agent execution.  GPU4PySCF ``sp``, ``opt``, and ``hess``
are configuration and preview paths until a compatible GPU target is
qualified.  PySCF CPU ``td`` is preview-only.

What the interface does not do
------------------------------

The interface does not invent coordinates, native input, server credentials,
or hidden resource defaults.  A provider turn cannot invoke ``/approve``.
There is no ``always allow`` mode and no standing grant; every workflow needs
its own displayed review and its own single human decision.  The interface
does not submit an unreviewed job and does not claim that a configured engine
is installed: program availability and requested resources must be observed
on the selected target before a real workflow is shown for approval.
