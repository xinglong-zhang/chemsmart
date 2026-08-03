# Forward agent-harness compatibility

This fork integrates a capability-driven `chemsmart.agent` package with the
v3.1.4 command surface. The core registry remains agent-independent; the agent
adds an evidence overlay rather than redefining program capabilities.

## Consume the canonical capability views

| Canonical name | Exact meaning |
|---|---|
| `PROGRAM_CAPABILITIES` | frozen declaration for Gaussian, ORCA, xTB, NCIPLOT, and PySCF |
| `EXECUTABLE_PROGRAMS` | all five executable Click program commands |
| `PRIMARY_PROGRAMS` / `COMPUTATIONAL_PROGRAMS` | molecular-method programs; excludes the advanced NCIPLOT leaf |
| `PROJECT_PROGRAMS` | programs that support project YAML |
| `PROJECT_REQUIRED_PROGRAMS` | programs for which synthesis must resolve project YAML |
| `PROJECT_PROGRAM_ORDER` | deterministic sequence for project discovery or first-match selection |
| `PROGRAM_CLI_JOBTYPES` | immediate Click children only; NCIPLOT is a leaf and therefore has an empty tuple |
| `PROGRAM_EXECUTION_ENGINES` | scientific execution choice (`cpu`/`gpu`), not a software-program synonym |
| `PROGRAM_PROJECT_OWNED_CLI_PARAMETERS` | CLI option names that an agent compiler may not silently replace instead of using approved project YAML |

Treat `ENGINE_CAPABILITIES`, `EngineCapability`, `engine_capability`,
`PROGRAM_JOBTYPES`, `PROGRAM_ENGINES`, and `PROJECT_OWNED_PARAMETERS` as v2
compatibility aliases only. Keep new code on the canonical `PROGRAM_*`
vocabulary.

## Apply the eventual one-line touch points

| Harness touch point | Registry-backed replacement |
|---|---|
| `cli_schema._PRIMARY_PROGRAMS` | `PRIMARY_PROGRAMS` |
| `command_workflow._COMPUTATIONAL_PROGRAMS` | `COMPUTATIONAL_PROGRAMS` |
| `command_workflow` mandatory-project gates | `requires_project_configuration(program)` |
| `command_workflow` project/program mismatch gate | `supports_project_configuration(program)` |
| `command_workflow._project_owns_parameter` | `project_owns_parameter(program, parameter)` |
| `command_workflow_tools` program extraction and JSON enum | `COMPUTATIONAL_PROGRAMS` and a sorted projection |
| `project_yaml_values.KNOWN_PROGRAMS` | `PROJECT_PROGRAMS`, not all executable programs |
| `workspace_bindings._discover_projects` | iterate `PROJECT_PROGRAM_ORDER` |
| `scientific_task._PROGRAMS` | `COMPUTATIONAL_PROGRAMS` |
| `harness.engine_capabilities` | thin re-export of the compatibility aliases and helpers |
| `harness.command_semantics._SOFTWARE_COMMANDS` | `COMPUTATIONAL_PROGRAMS` |
| `harness.command_workflow_receipt._COMPUTATIONAL_PROGRAMS` | `COMPUTATIONAL_PROGRAMS` |
| `settings.workspace_project.PROJECT_PROGRAMS` | membership view plus `PROJECT_PROGRAM_ORDER` where order affects behavior |
| `ParsedModelCommand` | add typed, serialized `engine: cpu|gpu|None` |
| generated-input invariants | add `generated_pyscf` delegating to `verify_provenance()` |

Build command and option behavior from the live Click schema. In particular,
replace copied Gaussian/ORCA child sets with the live walk; the integrated tree
has top-level `pka`, nested `qmmm`, xTB, and PySCF, which static legacy copies do
not describe.

## Preserve the semantic boundaries

- Do not use `PROGRAM_CLI_JOBTYPES` as the `scientific_task` validation map,
  concrete `JobKind` class registry, submit renderer list, planner vocabulary,
  NLP cue map, or expected-artifact map. Those lists describe implemented agent
  behavior and need explicit PySCF branches before they may broaden.
- Do not use `PROGRAM_PROJECT_OWNED_CLI_PARAMETERS` for settings-object fields
  or model/runtime-owned spec keys. The v2 renderer uses attribute aliases and
  the spec invariant owns additional keys and prefixes; centralize those layers
  separately if they are migrated.
- Runtime registry values cannot synthesize precise `typing.Literal` aliases.
  Centralize a static alias with an agreement test, or use `str` plus a runtime
  validator; never keep untested copies across modules.
- xTB has an explicit, narrow execution declaration. Its richer reader surface
  does not broaden the declared `sp|opt|hess` command surface. Keep method,
  solvation, project, and unsupported-workflow checks deterministic.
- PySCF completion and invariants come from `label.h5` schema and typed
  violations, never Gaussian/ORCA route regexes or the human `.out` log.

## Program-management boundary

Capability declaration, environment observation, substitution assessment,
program binding, preflight, execution, and result verification are distinct
states. Environment evidence must bind the finalized compute interpreter,
configuration identity, and execution phase. Aggregate preflight violations
for every selected molecule before constructing any executable job. A model
may propose a candidate, but it never establishes installation, readiness,
approval, execution, or success.
