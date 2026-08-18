# ChemSmart Product Charter

## Mission

ChemSmart is the canonical, CLI-first hub through which humans and AI agents
operate computational-chemistry programs. Scientific intent belongs in
readable project YAML and typed scientific DAGs. ChemSmart validates that
intent, materialises program-native inputs, compiles the public CLI, controls
execution, and returns typed scientific evidence.

The model is a computational scientist, not an input-file generator. It may
choose a defensible method, program, decomposition, and interpretation when a
task leaves them open. It must not bypass ChemSmart by inventing native input,
shell commands, execution status, or result values.

## Product boundary for version 3.1.4

The production Agent supports:

- project-YAML creation and validation;
- ChemSmart CLI compilation and safe preview;
- causal scientific workflow planning;
- inspection and typed analysis of supported results; and
- explicitly approved execution on release-qualified CPU paths: ORCA
  single-points, optimization/frequency, transition-state, excited-state,
  relaxed coordinate scans, and serial DAG workflows; PySCF ``sp/opt/hess``;
  and xTB ``sp/opt/hess``.

ORCA ``scan`` is qualified for approved execution: a relaxed torsional profile
ran through the ordinary plan, preview, single human approval, and provider-free
execution path, and its surface is read into typed quantities by the same
analysis layer as any other result. A scan's driven coordinate is carried on the
workflow node, not in project YAML, because it is a fact about this molecule in
this calculation rather than reusable method rationale.

ORCA ``modred`` is declared for planning, preview, and native-input generation
only. Constrained optimisation is expressible and previewable, and no
constrained optimisation has yet run here, so this release does not describe it
as completed Agent execution.

Gaussian ``sp/opt/ts/irc/td/link/scan/modred`` is supported for project YAML,
native-input generation, safe preview, and parsing of user-supplied completed
results; this release does not claim Gaussian Agent execution. GPU4PySCF
``sp/opt/hess`` is a PySCF-engine configuration and preview surface, not a
release-qualified Agent execution path. PySCF CPU ``td`` is likewise
preview-only. ORCA ``irc`` and ``neb`` may be planned and previewed, but require
target-specific qualification before they are described as completed execution.
NCIPLOT and additional human CLI families without an Agent declaration remain
outside the version-3.1.4 Agent execution surface.

Product support never asserts that an engine is installed on the current host.
Every real operation must pass its normal environment probe and appear in the
human review before it can run.

Runtime orchestration is provider-neutral. Version 3.1.4 contains registered
adapters for Alibaba Token Plan and DeepSeek. A user-selected profile supplies
the provider, endpoint, model, reasoning setting, and credential label; source
code and documentation must not impose a default model.

## Authority and approval chain

Planning, YAML validation, CLI compilation, safe preview, and result analysis
do not grant engine authority. Real calculation follows this chain:

1. the Agent produces a project-backed DAG;
2. ChemSmart compiles and safely previews every executable node while retaining
   any scientifically necessary release-unsupported stage as explicit
   non-executable intent;
3. the terminal interface displays the complete plan, marks non-executable
   stages and their reasons, and displays molecular identity, electronic state,
   effective project settings, CLI operations, dependencies, environment, and
   resources for the executable partition;
4. a human enters ``/approve`` once, or chooses ``/deny`` or ``/revise``;
5. the displayed workflow is removed from the pending state before launch;
6. a provider-free executor runs only that reviewed executable partition; and
7. ChemSmart records engine and validation evidence; a subsequent explicit
   analysis request reads the completed result into typed quantities,
   expressions, thermochemistry, claims, and interpretation.

There is no permanent calculation grant, session-wide "always allow", command
prefix allow-list, or model-created approval. A revised molecule, state,
project, environment, resource allocation, or DAG is a new workflow and must be
reviewed again. A multi-node causal workflow needs one human action because the
complete graph is displayed together; a displayed non-executable stage remains
unapproved and unlaunched.

The terminal UI is a view and controller for this chain. It is not a second
permission engine. Internal receipts and content digests preserve provenance
and mutation evidence; they are not hashes or approval-file tokens that a human
must retype in the production TUI.

## Scientific invariants

Before materialisation, establish the facts that determine meaning:

- molecular identity and the role of each geometry;
- coordinate units and atom order;
- charge, multiplicity, electronic state, and constraints;
- requested observable and physical conditions;
- method or program requirements fixed by the question; and
- whether the task requests planning, preview, analysis, or execution.

Ask rather than invent a consequential missing fact. Never infer identity or
state from a filename. Preserve artifact lineage across geometry handoff and
state changes. Keep signs, dimensions, units, standard states, temperature,
pressure or concentration, and thermochemical conventions explicit.

Normal process exit is not scientific validation. Distinguish, in order:

- proposed;
- planned;
- materialised;
- previewed;
- approved;
- executing;
- engine-complete;
- parsed;
- scientifically validated; and
- interpreted.

Only the deterministic host owns these states. Provider text is not execution
evidence, and hidden model reasoning is never scientific evidence.

## Product differentiation

ChemSmart does not compete by maximising autonomy or agent count. Its value is
the separation of flexible scientific reasoning from a reproducible,
multi-program execution authority:

- one public YAML-and-CLI layer instead of model-authored native inputs;
- molecular, electronic-state, artifact, and geometry-lineage preservation;
- preview and one explicit human decision over the displayed scientific and
  resource state;
- provider-independent execution semantics;
- native outputs plus typed, unit-aware analysis rather than transcript-only
  provenance; and
- explicit maturity claims for each program and operation.

Do not force one paper answer, molecule-specific branch, preferred DAG, tool
order, or reporting style. Algebraically equivalent transformations and
scientifically stronger program-native routes are acceptable when their
evidence chain is complete.

## Implementation discipline

- Treat live project loaders and Click commands as the public authority.
- Keep provider protocol code inside registered adapters.
- Use the smallest existing architectural layer that owns a defect.
- Do not create a parallel orchestration, scheduler, or grading system.
- Preserve unrelated working-tree changes; never reset, clean, or overwrite
  user work without explicit authority.
- Do not commit credentials, user configuration, engine binaries, generated
  inputs, outputs, scratch data, private transcripts, or one-off reports.
- Keep controller and program compute environments explicit in user or server
  YAML. Never replace an operator-selected executable implicitly.
- Validate a target host from its actual operating system, architecture,
  scheduler, program builds, and resource limits. No single cloud or server is
  the universal reference.

After a material change, run one focused mechanical check and then prefer a
decisive real scientific observation. Tests verify mechanics; they do not
grade computational-chemistry intelligence. Never claim an engine run from a
fake preview, fixture, parser test, or source inspection.

## Human scientific review

The human scientist owns interpretation and publication. Evaluate whether
identity, state, method, numerical transformations, units, conditions,
dependencies, and limitations are coherent. Accept creative valid routes.
Reject invented data, unperformed actions presented as completed, silent
changes to the scientific problem, and invalid chemistry or mathematics.

Report the route, strong scientific decisions, consequential limitations, the
general ChemSmart capability involved, and exactly what was planned,
previewed, executed, parsed, validated, or inferred.

## Documentation and repository hygiene

User documentation lives under `docs/source` and describes released public
behavior only. It must not contain development diaries, hidden evaluation
rubrics, private infrastructure, future implementation status, or internal
class inventories. `README.md` is a concise human entry point.

This charter and `.agents/skills/chemsmart-agent/SKILL.md` are the two
governance exceptions. Keep them aligned with the live product. The repository
source and CLI win if either instruction becomes stale.
