ChemSmart Agent terminal interface
==================================

The terminal interface presents the same provider-neutral Runtime V2 used by
``chemsmart agent plan`` and the same deterministic executor used by
``chemsmart agent execute``.  It is a human interface, not another chemistry
engine or approval system.

Install and open it
-------------------

Install the optional interface dependency alongside ChemSmart::

   python -m pip install 'chemsmart[agent-tui]'

Then select an existing provider profile, secret assignment file, and task
workspace. Planning and existing-result analysis need no execution bounds::

   chemsmart agent tui \
     --provider alibaba-production \
     --secret-file /path/to/provider.env \
     --workspace /path/to/task

Add an envelope and review destination only when the session may prepare real
execution for human approval::

   chemsmart agent tui \
     --provider alibaba-production \
     --secret-file /path/to/provider.env \
     --workspace /path/to/task \
     --execution-envelope /path/to/bounds.json \
     --review-file /path/to/task-review.json

Use ``--plain`` for inline rendering without a mouse or animations, such as
inside a simple SSH terminal.  The supported provider adapters are Alibaba
Token Plan and DeepSeek; the selected model and credentials remain user
configuration.

Plan and safe-preview chain
---------------------------

Enter a scientific request at the prompt.  The interface asks Runtime V2 to:

#. preserve the user-approved molecular identity and coordinate bytes;
#. create canonical ChemSmart project YAML;
#. compile the YAML into the live ChemSmart CLI;
#. inspect and safely preview the compiled operations; and
#. show the resulting scientific DAG, findings, and receipt digests.

The planning step never launches a chemistry engine.  Canonical YAML and
action-grade tool results are shown in the transcript; native input remains a
downstream product of ChemSmart.  The execution envelope supplies explicit
resource and program bounds but grants no authority.  When the preview is
ready, the TUI displays the exact inert review artifact and its SHA-256.

Explicit approval chain
-----------------------

Execution uses four distinct states:

#. The host writes and the TUI displays a ``WorkflowExecutionReviewV1``.  The
   review is inert and carries no grant.
#. The user reviews its complete scientific plan, materialized CLI operations,
   resources, environments, and SHA-256 outside the model turn.
#. The human chooses ``/approve ACTOR OUTPUT``, ``/deny ACTOR``,
   ``/revise ACTOR``, or ``/quit-review ACTOR``.  The same human-only resolver
   used by ``chemsmart agent approve`` records the exact decision.  The model
   and provider runtime cannot invoke this resolver.
#. Approval creates one exact, one-shot bundle and displays its file SHA-256.
   ``/execute APPROVAL_FILE_SHA256 RUN_DIRECTORY`` requires the complete file
   digest to be retyped before calling the provider-free deterministic
   executor.

``/deny`` forgets the current review and approval while retaining the safe
preview.  An approval cannot be changed into an “always allow” or session-wide
permission.  It cannot authorize a different plan, resource allocation,
workspace, environment, or artifact body.

Commands
--------

``/capabilities``
   Show the live program, engine and job matrix. This reports implemented
   product paths; the selected target still needs its own environment probe.

``/help``
   Show the in-terminal guide.

``/status``
   Show the task specification, review digest, approval file digest, and
   current state.

``/request PATH``
   Read and display another exact inert execution-review packet for the same
   task and workspace.  The review configured at launch is displayed
   automatically after planning.

``/approve ACTOR OUTPUT_PATH``
   Record explicit human approval of the displayed review and create one
   exact, one-shot bundle.  It does not execute anything.

``/approval PATH``
   Read a one-shot user-created approval bundle and prove that it binds the
   displayed review.

``/execute SHA256 RUN_DIRECTORY``
   Reconfirm the full approval file SHA-256 and invoke deterministic execution.
   This stage makes no model call.

After a completed execution, put its result artifacts in the task workspace
and enter a new scientific request to inspect and analyse them. This begins a
new provider-backed analysis turn; it does not extend the consumed execution
grant or rerun the engine.

``/deny ACTOR``
   Record denial of an unapproved review without creating a grant. A bundle
   already written to disk is not revocable; simply do not execute it and
   create a new review for revised work.

``/revise ACTOR``
   Record that an unapproved workflow requires revision, without creating a
   grant.

``/quit-review ACTOR``
   Record that review ended without a grant while leaving the TUI open.

``/quit``
   Exit when no host operation is running.

What the interface does not do
------------------------------

The interface does not generate coordinates, native input, server credentials,
or hidden resource defaults.  A provider/model turn cannot create an approval;
only an explicit human TUI or CLI action can.  The interface does not submit
unreviewed jobs or offer autonomous execution.  The product-capable execution
matrix is declared by ChemSmart's program registry and includes the documented
Gaussian, ORCA, xTB, PySCF CPU, and GPU4PySCF pairs.  That matrix does not claim
that this machine has an engine installed: the exact environment must pass its
host probe before the review can become executable. Gaussian remains dependent
on an operator-provided licensed executable and GPU4PySCF on a compatible CUDA
environment.
