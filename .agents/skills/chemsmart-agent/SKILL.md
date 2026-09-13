---
name: chemsmart-agent
description: Operate, audit, document, or improve the production ChemSmart computational-chemistry Agent through canonical project YAML, the live CLI, scientific DAGs, safe preview, visible human approval, deterministic execution, and typed result analysis. Use for ChemSmart Agent planning, TUI, approval, provider adapters, program support, or scientific workflow work.
---

# ChemSmart Agent

Treat ChemSmart as the execution authority and the model as the scientific
reasoner. Read `AGENTS.md` before changing the product boundary.

## Start from the live product

1. Inspect the project YAML loader and the exact Click command involved.
2. Establish molecular identity, geometry role, units, charge, multiplicity,
   electronic state, constraints, requested observable, and conditions.
3. Ask for a consequential missing fact instead of guessing it.
4. Use ChemSmart to materialise native input; never author a substitute native
   input or shell execution path in the model layer.

## Keep states honest

Name work as proposed, planned, materialised, previewed, approved, executing,
engine-complete, parsed, scientifically validated, or interpreted. Do not
collapse these states or treat provider text, a fake preview, or a fixture as
engine evidence.

## Use the complete scientific layer

The production surface is broader than command generation. Use live
capability/environment inspection, identity and project-YAML binding, causal
DAG/frontier operations, validated geometry handoff, registered-result
inspection, semantic quantity extraction, RRHO or parameterised qRRHO,
unit-aware expression DAGs, evidence-bound claims, and scientific decisions
when the task calls for them. The supported expression vocabulary includes
CBS extrapolation, Boltzmann populations and averages, harmonic ZPE,
imaginary-mode counts, geometry measurements, centres of mass, inertia,
rotational constants, and connectivity changes.

A frequency says how fast a mode moves, never which atoms move in it.
For that, ORCA opt/ts, xTB hess, and PySCF serve
``vibrational_mode_atom_participation``: each atom's share of a mode's
squared displacement, one row per mode summing to one. It is host-derived
from the displacement vectors the program printed, renormalised so it
means the same thing across programs whose vectors do not (PySCF stores
the same displacement scaled by 1/sqrt(reduced mass)). Read it as an
observation -- "these three atoms carry 96% of this imaginary mode" is
evidence; "this is a methyl rotor" is your claim, not the host's. Inside
a degenerate set the individual eigenvectors are an arbitrary basis, so
consult ``vibrational_mode_degeneracy_group`` before assigning motion to
one mode. Gaussian is deliberately undeclared: we never run it, and its
displacement block comes in variants (``freq=HPModes``, ``freq=raman``)
this reader cannot yet tell apart.

A completed solvated ORCA result carries its own solvation
decomposition: ``solvation_electrostatic_energy``,
``solvation_nonelectrostatic_energy`` and
``solvation_cavity_surface_area`` are declared for ``opt``, ``sp`` and
``ts`` beside ``solvation_model`` and ``solvent``. Read them together --
the terms say what the program applied, which can differ from what the
route asked for. A gas-phase result reports them absent, and a CPCM run
carries no cavity-dispersion term; that absence is how the models differ,
not a failure. Use them to decompose a delivered solvated number, never
to conclude that solvation is the source of a discrepancy: a basis-set or
functional error of the same size looks identical in that decomposition.

Per-atom charges arrive as a positional vector in molecular order, paired
with ``symbols``, and named by scheme. ORCA declares
``mulliken_atomic_charges``, ``loewdin_atomic_charges`` and
``hirshfeld_atomic_charges``; PySCF declares ``mulliken_atomic_charges``.
They are different quantities, not different spellings, and can disagree
by more than a third of an electron on the same atom, so quote the scheme
with the number. The two basis partitions divide a sum over basis
functions and inherit its sensitivity; Hirshfeld divides real space
against a promolecular reference and does not, which is why the
condensed-reactivity literature asks for it. ORCA prints Mulliken and
Loewdin unasked but prints Hirshfeld only when the route says
``Hirshfeld``, which the project escape hatch may carry -- a print
directive changes no method and every token you put there is displayed to
the reviewer on its own. CM5 stays undeclared.

The physical check worth doing on any per-atom vector is that it sums to
the formal charge; a session doing exactly that found a reader returning
spin populations for every open-shell result. Expect the two basis
partitions to close to their printed decimals and a basin partition to
close two orders of magnitude looser on its numerical grid -- both far
from what a dropped atom would cost. When you difference two per-atom
vectors across electronic states, check the sum in the direction where a
wrong reading would change its sign, not the direction where charges and
spins happen to agree.

Every program answers the shared selector vocabulary the same way, and
that now includes PySCF: its structured HDF5 result is a registered
reader with job-type declarations for ``sp``, ``opt`` and ``hess``, so
the capability query reports what it carries and the declaration gate
refuses a selector whose meaning was never audited for that job type.
Excitation energies come back in hartree there and in electronvolts from
the log-parsing programs; the reader states its own native unit and the
arithmetic is canonical either way, so never convert one yourself.
PySCF ``td`` is executable and declares the SCF set beside the excitation
set; a second declaration axis says whose density each value is.

A geometry may cross programs -- an xTB optimisation feeding an ORCA or
PySCF single point is the ordinary multi-program protocol, and the
handoff refuses any change of atom identity or order. Numbers may not
cross so freely. A typed value carries its unit and its dimension and not
the method that produced it, so the arithmetic will subtract a
tight-binding energy from a hybrid-DFT one without complaint. The
displayed analysis chain names the level of theory behind every input
for exactly this reason: mixing levels can be a composite method or a
mistake, and the host shows it rather than guessing which. Read that
column as necessary and not sufficient. A functional keyword is a name,
not a definition, and programs disagree about some names -- ORCA's
``B3LYP`` and PySCF's ``b3lyp`` differ in their local correlation and
gave total energies 0.24 hartree apart under two identical displayed
strings. Compare differences across programs, never total energies, and
say which variant each program means.

An electrode potential is expressible: charge is a dimension, so
potential is energy per charge and ``gibbs_to_redox_potential`` owns
E = -dG/(nF) and the IUPAC sign. Subtract a reference electrode as
ordinary arithmetic against a registered constant so the electrode you
chose stays visible. Read a constant's convention family and its purpose
before composing several by hand: a family says which scale an entry
sits on and nothing about its standard state, so two entries on one
scale can still need the term that bridges them, and where a finished
composed value is registered, prefer it. A chain drawing on two families
is displayed to the reviewer, never refused -- mixing can be deliberate.

Two refusals now happen at planning rather than after an engine runs. A
charge and multiplicity no molecule can have is refused when the state
is bound; every state the arithmetic permits is still admitted, because
choosing among possible states is what the calculation is for. And an
expression node that reads a value no earlier node or analysis input
provides is refused when planned -- expression nodes evaluate in the
order given, so define a value before the node that uses it.

Do not force every task through every layer. A valid route may be analysis-only
or may choose a different causal decomposition. Use only the operations needed
to answer the scientific request, and keep every source quantity and convention
visible.

## Apply the production support boundary

- Gaussian CPU ``sp/opt/ts/irc/td/link`` has project-backed planning, native
  input preview, and typed analysis of supplied completed outputs; do not claim
  Agent execution in this release.
- ORCA CPU planning covers ``sp/opt/ts/irc/td/neb/scan/modred``.
  Release-qualified execution covers single-points, optimization/frequency,
  transition-state, excited-state, relaxed coordinate scans, intrinsic
  reaction coordinates, and serial DAG workflows. Treat ``neb`` and
  constrained optimisation (``modred``) as preview paths until the selected
  target is qualified.
- PySCF CPU ``sp/opt/hess`` and xTB CPU ``sp/opt/hess`` have approved real
  execution paths. PySCF CPU ``td``, excited-root ``opt``
  and the ``mp2``, ``ccsd`` and ``ccsd(t)`` methods are executable and
  fixture-qualified, claimed as completed execution only from sealed cases.
- GPU4PySCF ``sp/opt/hess`` is a PySCF-engine configuration and preview
  surface until a compatible GPU target is qualified. NCIPLOT and other human
  CLI families without an Agent declaration are not Agent execution paths.
- Keep product capability, observed scientific evidence, and current-host
  readiness distinct. A supported Gaussian path does not imply a licensed
  executable is present; a supported GPU path does not imply a compatible GPU
  stack. The live environment probe and human review decide whether the exact
  operation can run here.
- Runtime semantics are provider-neutral. Registered adapters in this release
  are Alibaba Token Plan and DeepSeek, configured entirely by a user profile.
  There is no default model; the profile must explicitly state the selected
  model and its context/output limits.

## Use visible one-shot approval

Planning, YAML work, CLI compilation, safe preview, and result analysis are
non-executing. Before an engine launch:

1. produce the complete project-backed DAG;
2. compile and safely preview every executable node;
3. present every planned stage, marking any release-unsupported stage deferred
   with its reason, and present molecular/state identity, effective YAML,
   ChemSmart CLI operations, data handoffs, environment, and resources for the
   executable stages in the terminal interface;
4. let the human enter ``/approve`` once, or choose ``/deny`` or ``/revise``;
   and
5. hand only the displayed executable partition to the provider-free
   deterministic executor; deferred stages remain visible but unapproved.

The unit of the human decision is the goal. One ``/approve`` covers the
requested observables, every identity and electronic-state binding, the
physical conditions, the envelope with its engine-call, wall-clock, and
revision budgets, and the complete initial plan. Under that grant the
host may admit a revised workflow without a returning human action, but
only when the revision preserves identity, state, and conditions, stays
inside the envelope, and its session read the previous run's typed
terminal outcome (``inspect_run``) -- an evidence-free revision
returns to the human. Every admission is recorded with the composite
actor ``goal-approval:<goal-id>`` beside the human ``granted_by``; the
model never approves. A reviewer may freeze the exact initial plan and
keep the frozen grain for that goal, and a stop file cancels at the
next boundary. One driver runs every goal and the ``goal`` command,
the ``plan`` command, and the terminal interface are views of it; with
``--dispatch scheduler`` the approved run is submitted through the
server profile's scheduler, the goal parks, and the job's own tail
runs ``chemsmart agent wake`` to resume it at the outcome phase. A run
that ended in a repairable state (wrong stationary point, convergence
failure, timeout, memory limit, native error) opens a typed recovery
whose wake carries a repair menu; one that ended unanswerably returns
to the human. Say what you did with each offered route in
``record_scientific_decision``'s ``menu_route_dispositions`` (taken,
rejected, deferred, with the mechanism and the receipts); the next wake
shows them beside the menu. Declare your own prediction about the route
with ``role: diagnostic`` and a ``failure_update_rule``: it is scored
like any expectation and never owed. Every result the workspace record
or a run stream names is registered by digest under
``<program>-result-<sha16>`` and opens by that id in ``inspect_run`` and
``extract_result_quantities``; a 64-hex value is a digest, not an id.
Where the server profile names ORCA, preflight also runs ORCA's own
input check on the previewed input, bounded and never charged, and its
word appears beside the node on the review. Every executed result is judged on the program-neutral
stationary-point rule beside its program's validator, and the
capability receipt's coverage cell says which typed axes and validity
rules the host reaches for that program and jobtype. How each executed node ended is one typed vocabulary
derived from the sealed run stream; a failed output is nameable as a
``failed_result`` for inspection and never admissible for extraction or
geometry handoff.

A typed analysis chain displayed in the review is covered by the same single
approval: after every approved calculation node validates, the executor runs
the chain provider-free and renders a completed-analysis report, while
scientific interpretation and the recorded decision remain a session act. A
composed molecular arrangement (compose_molecular_arrangement) is host-built
from two identity-bound parents, and a derived species
(derive_molecular_species) is host-built from an ordered subset of one
identity-bound parent's atoms -- homolysis, deprotonation, and fragment
extraction are that one operation. Both display their lineage at review,
neither infers an electronic state, and the consuming stage of either is a
new workflow. A downstream node may consume a completed ORCA scan's
minimum-energy sampled point inside one approval via a declared producer
edge (the displayed review names the rule); carrying any other point is an
explicit post-scan binding and a new workflow. Typed chains select
host-owned literature constants by registered name through the
``constant`` operation -- a ``literal`` is recorded as model-authored, a
``constant`` as host-owned with its standard-state convention displayed --
and named convention operations own their mathematics (``gibbs_to_pka``).
Composed workflows such as an aqueous pKa need no task-specific code, and
completed registered results may feed a later workflow's analysis.
A batch is N enumerated records under the same single approval: a
workspace chemsmart .db is inspectable (stored fields are observations,
never bindings -- a record may store no electronic state at all), a
record's geometry is extracted with full lineage and its state bound
explicitly like a derived species, N records plan as N disconnected
sub-DAGs in one workflow, the displayed review carries every record
row, execution is sequential and record-major under the displayed
envelope (episode window and engine-call budget enforced by the
provider-free executor, replays counted), one record's failure settles
typed while the others deliver, per-record verdicts render with
deliberately no aggregate quantity, and an interrupted or partial run
continues by re-entering its own run directory -- each resume recorded
on the consumption ledger, terminal nodes replayed without
re-execution, a mid-engine interruption reported as ambiguous, and a
completed approval refusing re-invocation.
A geometry may be built by typed host-owned spatial operations:
edit_molecular_geometry sets one internal coordinate (bond length,
angle, or torsion) as a rigid motion of a named side -- the moving
side is the model's choice, named by one of the coordinate's own
atoms and never defaulted, and the receipt lists every atom that
moved with the coordinate measured before and after; refusals are
structural only (unbonded axis, ring, collinear atoms) and a
requested value is never judged, because the consuming optimisation
grades it. append_molecular_atom is derivation's mirror, placing one
atom by three internal coordinates against three anchors. Both leave
electronic state unbound for explicit binding, both make starting
structures whose consuming stage is a new workflow, and the review
renders every hop of a built chain root-first so the edit that
decides what the molecule is stays on the decision surface.

Never create an approval on the model's behalf. Never offer permanent,
session-wide, prefix-based, or "always allow" chemistry execution. A revised
scientific input, project, environment, resource request, or DAG requires a new
review. Internal receipts remain provenance; never make a human retype their
digests as a second scientific authority. A reviewed multi-node DAG needs one
approval, not one prompt per node.

## Improve at the owning layer

Classify a failure as scientific reasoning, program limitation, environment,
parser, or missing ChemSmart affordance. Fix the smallest general project
setting, CLI compiler, adapter, runner, parser, typed analysis operation, or
error message. Do not learn a molecule, paper answer, DOI, private DAG, tool
order, or provider-specific workaround.

Use one focused mechanical check after a change. Then prefer one decisive
scientific observation through the real public surface. Report exactly which
program ran, which artifact was parsed, and which scientific validation was
or was not established.

## Review as a computational scientist

Check identity and atom order, electronic state, method semantics, geometry
handoff, convergence, stationary-point evidence, signs, units, physical
conditions, and causal dependencies. Accept any chemically and
mathematically sound route. The final interpretation and publication decision
belong to the human scientist.
