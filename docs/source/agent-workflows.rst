###############################
 Computational Agent Workflows
###############################

ChemSmart Agent turns a scientific request into project YAML, compiled ChemSmart CLI operations, a safe preview, and
typed result analysis. The model proposes the chemistry; ChemSmart owns native-input generation, execution, validation,
and result parsing.

Approved execution stops an engine and its whole process tree through POSIX process groups and signals, so the Agent
layer runs on Linux and macOS. The human CLI remains supported on Windows.

*****************
 Supported scope
*****************

The current release provides broad planning, project-YAML generation, live CLI compilation, and safe preview for these
Agent job families:

-  Gaussian CPU: ``irc``, ``link``, ``modred``, ``opt``, ``scan``, ``sp``, ``td``, and ``ts``.
-  ORCA CPU: ``irc``, ``modred``, ``neb``, ``opt``, ``scan``, ``sp``, ``td``, and ``ts``.
-  PySCF CPU: ``hess``, ``opt``, ``sp``, and preview-only ``td``.
-  GPU4PySCF: PySCF ``hess``, ``opt``, and ``sp`` with the ``gpu`` engine.
-  xTB CPU: ``hess``, ``opt``, and ``sp``.

Release-qualified real Agent execution is narrower:

-  PySCF CPU ``sp``, ``opt``, and ``hess``;
-  xTB CPU ``sp``, ``opt``, and ``hess``; and
-  ORCA CPU single-points, optimization/frequency, transition-state, excited-state, relaxed coordinate scans, and serial
   producer-to-consumer DAGs.

A relaxed scan states its driven coordinate on the workflow node rather than in project YAML: which atoms, which
coordinate type, and the range and number of points. The same specification is rendered into each program's own idiom,
so one physical description reaches ORCA's absolute endpoints and Gaussian's increment without the caller writing
either. A completed scan's surface is read into typed quantities through the ordinary analysis layer.

The first constrained optimisation of a scan imposes the driven coordinate on the geometry supplied, so a range
beginning far from that geometry's current value may be refused by the program before any optimisation runs.

A downstream node may consume a completed ORCA scan's minimum-energy sampled point inside the same approval: declare the
consumer's geometry input as a producer edge from the scan node, and the host carries the point with the lowest recorded
energy (ties resolve to the lowest point index), verifying atom identity and the approved electronic state at handoff.
The displayed review names the rule, so approving the workflow approves exactly that settlement. Carrying any other
point of the surface remains an explicit post-scan choice whose consuming stage is a new workflow. This is the ordinary
escape from a torsional saddle: scan the dihedral, then optimise from the carried well.

Typed analysis chains may select host-owned literature constants by registered name through the ``constant`` expression
operation — an aqueous proton free energy, a standard-state correction, a reference acid's measured pKa — each carrying
its value, unit, and standard-state convention. A ``literal`` you supply is recorded as model-authored; a ``constant``
is host-owned, and an unregistered name is refused when the chain is planned. Named convention operations own their own
mathematics (for example ``gibbs_to_pka`` owns pKa = ΔG/(RT ln 10)). The review and the completed-analysis report render
a Literature-constants table whenever a chain selects one. Workflows composed this way — for example an aqueous pKa from
solvated optimisations, a derived conjugate base, thermochemistry at an explicit 1 mol/L standard state, and a registry
proton constant — need no task-specific code, and completed registered results may feed a later workflow's analysis as
typed inputs. A wavenumber becomes an energy only through ``wavenumber_to_energy`` and ``energy_to_wavenumber``, which
own the h·c·N_A factor, so an exchange coupling declared in cm⁻¹ is answered in cm⁻¹; when a claim carries its declared
id in another dimension, the completion miss names both dimensions and this route, and the expectation row compares the
sign and the host-restated band rather than printing ``not_comparable``. The registry also holds a ferrocene reference
for acetonitrile: a computed absolute Fc⁺/Fc potential with its source's stated accuracy, whose purpose phrase prefers a
ferrocene pair computed at the same level, and the experimental construction it was benchmarked against as one family.

A batch is N enumerated records under the one displayed approval. A workspace chemsmart ``.db`` database is an
inspectable artifact: ``inspect_database_records`` enumerates records with their stored fields as observations (a record
may store no electronic state at all), and ``extract_database_record_geometry`` copies one record's exact coordinates
into a lineage-carrying workspace geometry artifact — database digest, record id, explicit structure selection, with
multi-structure ambiguity refused — after which identity and electronic state are bound explicitly, exactly as for a
derived species; execution never reads the database again. N records are planned as N disconnected sub-DAGs in one
workflow, the review shows one row per record with its bound state, origin, and any stored-versus-bound mismatch
flagged, execution is sequential and record-major under the displayed envelope (episode window and engine-call budget
enforced by the provider-free executor, replays counted), and one record's failure settles that record as typed findings
while the others deliver — the report carries per-record delivery verdicts and deliberately no aggregate quantity. An
interrupted or partial run continues by re-entering its own run directory: each resume is recorded on the approval's
consumption ledger naming the remainder, terminal nodes replay from their receipts without re-execution, a mid-engine
interruption reports as ambiguous pending human reconciliation, and a completed approval refuses re-invocation.

A geometry may be built by typed spatial operations whose arithmetic the host owns. ``edit_molecular_geometry`` sets one
internal coordinate of an identity-bound geometry — a bond length, an angle, or a torsion, the same three coordinates a
scan drives — as a rigid motion: you name the coordinate's atoms, the target value in the coordinate's own unit, and
which side moves, named by one of the coordinate's own atoms because the choice is scientific and the common libraries
disagree on a default. The host measures the coordinate before and after with the same arithmetic the analysis layer's
``distance``/``angle``/``dihedral`` operations use, verifies the requested value was reached, enumerates every atom that
moved, and records close contacts and connectivity changes as observations rather than verdicts. Refusals are structural
only — an axis that is not a perceived bond, a coordinate a rigid motion would tear out of a ring (use a constrained
optimisation or a relaxed scan there), collinear or out-of-range atoms; a requested value is never judged, because the
optimisation that consumes the edited structure is what grades it, and the gap between what was requested and what
relaxation returns is the measurement. ``append_molecular_atom`` is derivation's mirror: it adds one atom placed by a
bond length, an angle, and a dihedral against three anchor atoms — protonation, hydrogenation, capping — leaving parent
atom indices unchanged. Both operations write starting structures with electronic state deliberately unbound (adding a
hydrogen gives a cation or a radical depending on whether it brought an electron), so charge and multiplicity are bound
explicitly afterwards and the consuming stage is a new workflow; the displayed review renders every hop of a built chain
in the order it was performed, so the edit that decides what the molecule is stays on the decision surface.

A built or idealised start carries its builder's symmetry and relaxes to the nearest stationary point of that symmetry,
which is a saddle whenever the true minimum lies lower. Every identity binding and every compiled node therefore states
a point-group estimate found within 0.01 Å from the molecule's own atoms (and within 0.1 Å when the two differ), with
the count of appended atoms placed on the exact 60° torsion lattice or at idealised angles, as a host observation on the
review. ``break_symmetry`` perturbs every atom of an identity-bound geometry by a seed and an amplitude you name, so the
same request gives the same bytes; the receipt carries the largest step actually taken and the point-group estimate
before and after, and the perturbed geometry is a starting structure with no electronic state bound.

A vibrational frequency says how fast a mode moves, never which atoms move in it.
``vibrational_mode_atom_participation`` answers the second question: each atom's share of a mode's squared displacement,
one row per mode summing to one, so a row reads as "how much of mode k does atom i carry". It is derived by the host
from the displacement vectors the program printed and renormalised, which is what lets one selector mean the same thing
for programs whose stored vectors differ — ORCA, Gaussian and xTB print Cartesian displacements at unit norm, while
PySCF returns the same physical displacement scaled by one over the square root of the reduced mass. Dividing by the
row's own total removes that per-mode scalar, the arbitrary sign of an eigenvector, and the program's choice of frame;
it does not remove the program's atomic mass table, which perturbs the vectors themselves by about a tenth of a percent
for C/H/O and about one percent for heavy halogens. The number is an observation — that three hydrogens carry 96% of an
imaginary mode is evidence, and calling it a methyl rotor is your reading of that evidence. Inside a degenerate set the
individual eigenvectors are an arbitrary basis, so ``vibrational_mode_degeneracy_group`` reports which modes share a
frequency within a stated tolerance and lets you see that a mode has company before assigning motion to it. Both are
declared for ORCA ``opt`` and ``ts``, xTB ``hess``, and PySCF. Gaussian is deliberately undeclared: ChemSmart never runs
Gaussian, and its displacement block differs between the default output, ``freq=HPModes`` and ``freq=raman``, so which
block a supplied log contains is not something this reader can yet establish.

ORCA ``irc`` is qualified for approved Agent execution on a qualified target: a transition-state search fed two IRC
runs, each consuming the converged transition state's own geometry and analytic Hessian as role-distinct producer
bindings inside one approved workflow, and the chain executed, validated, and delivered host-rendered claims. The IRC
log prints only the starting structure, so only job-level facts (charge, multiplicity, direction, solvation route, atom
identity) are declared for the jobtype and selector declarations gate extraction; the reaction path lives in the
trajectory artifact, which enters the typed layer as a registered geometry artifact — endpoint connectivity is read from
that artifact by the scientist, not rendered as a host claim. ORCA ``neb`` and ``modred`` remain planning and preview
paths until the selected target is qualified. Gaussian Agent execution is not claimed in this release; Gaussian support
covers project YAML, generated native input, safe preview, and typed analysis of user-supplied completed outputs.
GPU4PySCF remains a configuration and preview path until a compatible GPU target is qualified.

These boundaries do not alter the wider human ``chemsmart run`` and ``chemsmart sub`` CLI. They also do not imply that
an executable is available on the current machine. Every approved CPU run still needs an observed program environment
and sufficient target resources. Gaussian requires a separately licensed installation and GPU4PySCF requires a
compatible CUDA stack.

**********************
 Configure a provider
**********************

The runtime is provider-neutral. The current release registers Alibaba Token Plan, DeepSeek, and OpenAI adapters; an
``anthropic`` profile is accepted as configuration and refuses execution until its adapter is registered. There is no
default model: every selected profile names its model, endpoint, context and output limits, reasoning setting, and
credential label explicitly. Credentials never live in ``agent.yaml``: they resolve from an exported environment
variable matching the profile's key label, or from the managed key store at ``~/.chemsmart/agent/keys.env`` that
``chemsmart config agent`` maintains. The guided setup asks for the provider, the exact model id, the effort, the token
limits, and the credential (hidden input):

.. code::

   chemsmart config agent

The following example defines one profile for each registered adapter. Replace the placeholders with values supported by
the selected account and model.

.. code:: yaml

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

Both token limits are required positive integers and ``max_output_tokens`` cannot exceed ``context_tokens``. A profile
may add ``record_reasoning: true`` to keep each turn's provider reasoning in the private run directory for that
campaign; it never enters the public transcript, and the turn receipt says whether it was kept. ``fallback`` must be
empty because ChemSmart does not switch providers inside one session. After a provider failure, start a new, explicitly
attributed attempt with another profile.

The managed key store contains the label selected by the active profile:

.. code::

   ALIBABA_TOKEN_PLAN_KEY=your-secret-value

ChemSmart parses this file as data; it does not source it into the shell, and an exported environment variable with the
same label takes precedence. Provider reasoning, response text, and credentials are not scientific execution evidence.

***********************
 Plan and safe preview
***********************

Use a workspace containing the molecular, project, and completed-result artifacts relevant to the task. This command
contacts the selected provider but cannot start a chemistry engine:

.. code:: bash

   chemsmart agent plan \
     --provider PROFILE \
     --provider-config /secure/path/agent.yaml \
     --task-file task.md \
     --workspace /absolute/path/task-workspace

``chemsmart agent plan`` may create and validate project YAML, compile live Click commands, inspect generated artifacts,
and execute fake previews. It never accepts a model-created execution grant; real execution belongs to the approval
chain below.

**********************************************
 Approved execution in the terminal interface
**********************************************

Real execution starts from explicit program and resource limits. The envelope below permits up to four serial CPU engine
calls on the listed programs. It is a bound on a possible run, not permission to run one.

.. code:: yaml

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

The envelope supplies the explicit per-run cores, memory, and GPU allocation, overriding ordinary server defaults. The
human review displays that allocation; the selected program, operating system, and any external scheduler must still be
able to provide it. An optional ``max_excursion_calls`` (default 0) grants a separate, displayed line for nodes that
investigate an anomaly the host recorded: such a node is tagged with the anomaly receipt's digest, is charged to that
line and never to ``max_engine_calls``, and may feed no untagged node. Open the terminal interface with the envelope:

.. code:: bash

   chemsmart agent tui \
     --provider PROFILE \
     --provider-config /secure/path/agent.yaml \
     --workspace /absolute/path/task-workspace \
     --execution-envelope /absolute/path/resources.yaml

After planning and safe preview, the interface displays every planned stage, including a release-unsupported stage that
remains scientifically necessary. It marks that stage deferred, gives its reason, and displays the molecule and state,
effective project setting, ChemSmart CLI operation, data handoff, program environment, and resource bound for every
executable stage. Enter ``/approve`` once to execute only those reviewed executable stages; deferred stages remain
unapproved and are not launched. Use ``/deny`` or ``/revise`` without launching an engine. The provider is disconnected
before execution.

The human does not retype a hash or create an approval-file token. Internal receipts and content digests remain
provenance in the durable records; the terminal interface displays none of them, and the headless commands print human
summaries by default with ``--json`` for the exact machine payloads. The pending workflow is consumed before launch, so
a failed or completed run requires a fresh plan and human review before another execution attempt.

Use ``--review-file /absolute/path/review.json`` only when an audit application needs an exported copy of the displayed
workflow. The file is optional and is not the execution authority for the terminal interface.

***********************************
 Goals, scheduler dispatch, waking
***********************************

``chemsmart agent goal`` drives one goal to settlement under one human decision: the plan is reviewed and displayed,
execution is provider-free, and after each run the session reads the typed terminal outcome and may revise its route
within the envelope's budgets. Every entry point -- the ``goal`` command, ``chemsmart agent plan``, and the terminal
interface -- is a view of the same driver.

.. code:: bash

   chemsmart agent goal \
     --task-file /absolute/path/TASK.md \
     --workspace /absolute/path/goal-workspace \
     --execution-envelope /absolute/path/resources.yaml \
     --goal-id GOAL-ID \
     --granted-by HUMAN \
     --dispatch scheduler

With ``--dispatch scheduler`` the approved run is submitted through the current server profile's scheduler (or the one
named by ``--server``) using the same submitters as ``chemsmart sub``. The command records the job it created in
``dispatch.receipt.json`` inside the run directory and parks the goal. The job script runs the provider-free executor in
that run directory, writes its result as ``execution-result.json``, and then runs:

.. code:: bash

   chemsmart agent wake --workspace /absolute/path/goal-workspace --goal GOAL-ID

which rebuilds the driver from the goal's own ledger and settles the goal or wakes its next cycle under the budgets the
human approved. A human may run the same command; ``--wait`` polls the scheduler until the job is over first. A settled
goal and a goal with no parked run are refused. Engine wall time is charged from the run's receipts; queue wait is
recorded beside it and never against the engine budget.

A run that ended in a state a revision can answer -- a stationary point of the wrong order, a convergence failure, a
timeout, a memory limit, a program error -- opens a typed recovery when budget remains, and the next session's wake
context carries a repair menu naming the ordinary route for each such ending. A run that ended in a state no revision
can stand on returns the goal to the human.

The recorded decision may say what was done with each route the menu offered (``menu_route_dispositions``: taken,
rejected or deferred, with the mechanism and the receipts); the host checks that the route was offered and the receipts
are its own, and the next wake shows the dispositions beside the menu it re-offers. A declaration may carry ``role:
diagnostic``, the session's own prediction about the route with a ``failure_update_rule`` and an optional
``method_resolution``; it is scored like any expectation and never owed, and a value inside the resolution prints
``indeterminate``.

A surprise the host detects on a completed node -- a stationary point of the wrong order, a walk to another basin, a
spin expectation value far from the bound state -- is recorded as an anomaly observation with its numbers beneath the
verdict, whether or not it was asked for. A certified delivery that carries one settles ``achieved_with_observations``
and its report names each observation; an excursion that re-runs the sensor marks the observation ``replicated`` or
``refuted``, and the latest receipt speaks for it. A scan whose extremum lies on the grid's edge, and a result that
walked from the goal's original bound geometry even when its immediate input was a reached or displaced structure, are
recorded the same way.

A declared observable is delivered when any cycle of the goal claimed it by id, and a goal whose declared observables
remain undelivered while budget remains is woken once more. A session that cannot reach a declared observable refuses it
in ``record_scientific_decision``'s ``unreachable_observable_ids``, naming the producer it would need and the receipts
that show the gap; the host verifies a selector no envelope program declares or a blocked node in the session's own
plan, and only a verified refusal settles ``unreachable_from_evidence``. The declaration reply says at once when a
meaning names a quantity kind no envelope program declares, with the routes, and refuses nothing. When a goal session is
about to end with declared observables undelivered while budget remains, the host says so once, as an informational
notice that demands nothing, and allows one further turn. Every routed refusal reaches the durable stream as a failure
report naming its gate, invariant, diagnosis, route and cost.

Every result the workspace record or a goal's run streams name is registered at bootstrap under
``<program>-result-<sha16>`` when its bytes still hash to the recorded digest, in this workspace or an earlier one, so
``inspect_run`` and ``extract_result_quantities`` open it by that id; a refused id is answered with the nearest
registered ones. Where the active server profile names an ORCA executable, preflight also runs ORCA's own input check on
the previewed input, bounded to 20 s and stopped at the ``INPUT FILE`` banner; its word (passed, aborted with ORCA's
lines, or not run) appears beside the node on the review and is never an engine call. The wake's budgets show
``host_seconds_spent`` beside the engine lines and the goal's ``input_check_probes`` counts.

*************************
 Guides and capabilities
*************************

The tools the model reads are a stem and leaves. Every session reads the stem: the core tools, the operations that
belong to no family, and the universal rules. A guide is a family unit -- ``structure``, ``scan``, ``constants``,
``cbs``, ``ensemble``, ``spectroscopy``, ``database``, ``crossprogram``, ``recovery``, ``saddle`` -- of extra tools,
extra operations, and a few hundred words of guidance. The host opens a guide when the task text, the workspace, the
planned workflow (its job types, operations, or a second program), or a previous run's ending calls for it, and the
model may open any guide itself with ``open_guide``. Each activation is recorded with the tool-schema digest it
produced.

.. code:: bash

   chemsmart agent capabilities [--kind tool|program_jobtype|operation|...] [--json]

lists every capability of the agent on one ladder -- declared, wired, advertised, tested, qualified -- computed from the
registries that own each kind. An executable program and job type is qualified when a live run stands behind it in the
release record; a claim without a recorded run is displayed as a claim, and a job type the agent can run but cannot read
or judge is displayed as unsupported.

********************
 Molecular identity
********************

When a task uses named geometries, ``--identity-manifest`` can bind each workspace-relative XYZ file to its approved
name, geometry role, coordinate units, charge, multiplicity, and source. Filenames and XYZ comments alone are not
molecular-state evidence.

.. code:: yaml

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

The content fields in this optional manifest preserve input provenance. They are separate from the human ``/approve``
action and are never retyped in the terminal interface.

**********************
 Causal data handoffs
**********************

Downstream calculations consume typed artifacts rather than guessed filenames. An optimization may provide its validated
final geometry to a Hessian, single-point, or excited-state node while preserving atom order, charge, and multiplicity.
A project that already requests frequencies produces them on that same scientific node; a duplicate Hessian is not
required unless the scientific protocol asks for an independent calculation.

For an ORCA transition-state-to-IRC workflow, the IRC node has two scientific inputs: the final TS geometry and the
final ORCA Hessian. ChemSmart selects the unique Hessian whose atom order, geometry, frequencies, and single
consequential imaginary mode agree with the validated TS output. This materialisation capability does not by itself
claim a completed IRC run on an unqualified target.

******************************
 Provider waits and deadlines
******************************

Provider profiles may set connect, first-event, inter-event, and absolute-turn deadlines. The absolute deadline is not
extended by stream heartbeats or partial bytes. During an approved chemistry run the executor makes no provider call;
process timeout and termination remain the ChemSmart runner's responsibility.

***********************
 Interpreting evidence
***********************

Keep these states distinct: proposed, planned, materialised, previewed, approved, executing, engine-complete, parsed,
scientifically validated, and interpreted. Provider prose, a fake preview, a parser example, or a successful process
exit alone is not a scientific result.

Before using a value, check molecular identity and atom order, charge and multiplicity, method and program semantics,
convergence, stationary-point evidence, geometry handoff, physical conditions, signs, dimensions, and units.

*******************
 When a node fails
*******************

A node whose engine run did not succeed reports its terminal state, the wrapper and child exit statuses, the validator
findings, and a bounded quotation of what the program itself printed about the failure.

That quotation is the program's own text, not a ChemSmart claim. It is evidence that a diagnosis exists and what it
said; it never establishes readiness, validity, or what to do next. URLs, absolute and home-relative paths, e-mail
addresses, and credential-like assignments are removed from it.

A run that did not terminate normally still yields no scientific quantities. Typed extraction continues to require a
normally terminated result, so a failed run can be read for its reason and not for its numbers.

*********************************
 Re-running an approved workflow
*********************************

``chemsmart agent review`` re-presents a stored execution review so the same workflow can be decided on again:

.. code::

   chemsmart agent review \
     --review-file review.json \
     --workspace /path/to/workspace

With no ``--decision`` it reports only: the review digest, the workflow, whether every approved artifact is still
present under the workspace, and which approval identities for this review have already been consumed. Adding
``--decision approve --actor NAME`` records a new human decision and writes a new one-shot bundle.

This does not reuse a spent approval. Approvals remain one-shot and bound to the exact request digest; replay obtains
the *current* decision that the approval chain requires for a launch, over an unchanged displayed plan.

Differences between the approval and the present workspace are displayed rather than silently accepted, and they change
no enforcement: the environment and command comparison that runs immediately before the first dispatch still refuses a
launch whose facts have drifted from the reviewed ones. Approved input bytes that are no longer present under the
workspace are refused before a decision is offered, because resolving them would fail before anything could run.

*************************
 Domain-knowledge skills
*************************

The Agent can consult short, advisory skill documents that describe general computational-chemistry practice. They are
surfaced by name in the planning prompt and fetched on request; they are advisory only and never establish readiness,
approval, terminal state, or an accuracy claim, and they never replace a typed host receipt.

Released skills:

-  ``scientific-conventions`` — how computed quantities are conventionally defined and reported: the direction of every
   difference quantity, adiabatic versus vertical geometry, which energy terms are included, and thermochemistry
   standard states.

-  ``method-adequacy`` — whether a chosen method, basis set, solvation model or conformer sample can answer the question
   being asked: which errors cancel in a comparison and which do not, and how to state the resulting uncertainty.

-  ``typed-analysis-contract`` — how the typed analysis layer expects intent to be expressed: identifiers, units,
   declared quantity kinds, and evidence references.

Set ``CHEMSMART_SKILL_ROOT`` to add or override skills from a directory of your own; set ``CHEMSMART_AGENT_SKILLS=0`` to
remove the skill index and the consultation tool entirely.

***************************
 Analyse completed results
***************************

Analysis is a normal Agent mode and does not require an execution envelope or approval. Put supported completed results
in a task workspace and request the quantities or comparison needed. ChemSmart discovers:

-  normally terminated Gaussian and ORCA native outputs;
-  validated xTB result folders, including relocated archives whose missing original paths are reported as provenance
   limitations;
-  structured PySCF HDF5 results; and
-  XYZ geometries and trajectories.

For example:

.. code::

   chemsmart agent plan \
     --provider PROFILE \
     --provider-config /secure/path/agent.yaml \
     --workspace /absolute/path/completed-results \
     --task "Extract the final energy and frequencies, diagnose the stationary point, and report the result with explicit units."

The Agent may extract energies, structures, frequencies, dipoles, excited states, spin evidence, solvent treatment, and
auxiliary-basis roles. Text metadata remains typed text and is not assigned a fictitious physical unit. It may derive
RRHO or explicitly parameterised quasi-harmonic thermochemistry and evaluate unit-aware expression DAGs for, among other
operations:

-  energy differences and CBS extrapolation;
-  Boltzmann populations or averages from one energy vector or separate scalar state energies, with optional
   degeneracies;
-  harmonic ZPE and imaginary-mode counts;
-  distances, angles, dihedrals, centres of mass, inertia, rotational constants, and connectivity changes; and
-  photon wavelengths and dimensional unit conversions.

For a composite high-level electronic energy plus lower-level quasi-harmonic thermochemistry, use
``quasi_harmonic_thermal_gibbs_correction = G_qh(T) - E_electronic``. Do not substitute the harmonic
``thermal_gibbs_correction`` while describing the result as quasi-harmonic.

Supply ``--analysis-completion-file`` only when another application needs a host-authored list of mandatory quantities
or claims. Ordinary existing-result analysis does not require it.

A typed analysis chain planned with a workflow travels verbatim in the review packet and the approval bundle, and the
single approval covers it. After every approved calculation node validates, ``agent run`` runs the chain provider-free
-- extraction, thermochemistry, expressions, validation verdicts, and claim rendering -- and writes ``analysis/
completed-analysis-report.md`` into the run directory. A failed validation verdict is a completed determination, not an
execution failure. Scientific interpretation and the recorded decision remain a session act: start a new analysis
request over the completed result workspace to record them; this does not rerun the engine or extend the earlier
execution decision. A workflow approved without an analysis chain executes exactly as before.

Three further affordances round out multi-stage work. ``compose_molecular_arrangement`` places two identity-bound
geometries into one arrangement at an explicit atomic contact; the host owns the placement mathematics and the composed
bytes with full parent lineage, the electronic state is bound explicitly afterwards, and the consuming stage is a new
workflow. ``derive_molecular_species`` is its mirror, taking an ordered subset of one identity-bound parent's atoms --
the operation underneath homolysis, deprotonation, and fragment extraction. Name either the atoms to remove or the atoms
to keep; the host records both, copies the parent's coordinates unchanged so the result is a starting structure rather
than a relaxed one, and reports whether the result is one species or several separated pieces. Derivation never infers
an electronic state, because removing a hydrogen gives a radical or an anion depending on where its electron went. The
selection rule ``validated_producer_orca_hessian`` declares a validated frequency-bearing ORCA producer as a legal
source for a transition-state search's ``--inhess-filename`` starting Hessian; the starting Hessian may carry any
imaginary-mode count, and the observed count is recorded in the handoff receipt. A workflow consuming this rule freezes
into an approved bundle when the geometry also arrives from the producer; none has yet executed here, so
producer-Hessian seeding is not described as completed Agent execution.
