# CHEMSMART

[![CI](https://github.com/Hongjiseung-ROK/chemsmart/actions/workflows/main.yml/badge.svg)](https://github.com/Hongjiseung-ROK/chemsmart/actions/workflows/main.yml)

<p align="center">
  <img src="docs/source/_static/chemsmart_logo.png" alt="CHEMSMART logo" width="600">
</p>

CHEMSMART is a CLI-first, project-YAML hub for computational chemistry. Human
researchers and AI agents use the same validated project settings and commands
to prepare, preview, run, submit, and analyse calculations across supported
programs.

Scientific choices remain visible in project YAML. CHEMSMART translates those
choices into program-native inputs, compiles the public CLI, controls execution,
and exposes typed, unit-aware results. Generated native files are downstream
artifacts; they are not a second user or model API.

Current package version: **3.1.4**. Python **3.10** is required.

## Product model

```text
scientific question and molecular artifacts
                  |
                  v
        project YAML + typed DAG
                  |
                  v
       CHEMSMART compile and preview
                  |
       human approval for exact DAG
                  |
                  v
        deterministic CLI execution
                  |
                  v
 native output + typed scientific analysis
```

Planning, YAML generation, CLI compilation, safe preview, and result analysis
do not start a chemistry engine. Real Agent execution requires a one-shot human
approval bound to the exact workflow, molecular state, inputs, project files,
commands, environment, and resources. The approved executor uses no model or
provider credential.

## Human CLI support

The human CLI includes Gaussian, ORCA, xTB, PySCF, GPU4PySCF as a PySCF engine,
NCIPLOT, molecular databases, pKa and thermochemistry analysis, ITERATE, and
PyMOL-based molecular operations. Each external program requires a compatible,
separately licensed installation where applicable.

Use the live command help as the option authority:

```bash
chemsmart --help
chemsmart run --help
chemsmart sub --help
chemsmart run PROGRAM --help
```

## Agent support in version 3.1.4

| Program | Planning and safe preview | Explicitly approved execution | Result analysis |
| --- | --- | --- | --- |
| PySCF CPU | `sp`, `opt`, `hess`; `td` preview | `sp`, `opt`, `hess` | structured HDF5 quantities actually produced |
| GPU4PySCF | `sp`, `opt`, `hess` | `sp`, `opt`, `hess` | the same structured PySCF result path |
| xTB CPU | `sp`, `opt`, `hess` | `sp`, `opt`, `hess` | receipt-bound native quantities, orbitals, frequencies and geometry handoffs; XYZ trajectory analysis uses the shared geometry reader |
| ORCA CPU | `sp`, `opt`, `ts`, `irc`, `td`, `neb` | the same six job families | native energies, structures, frequencies, excited states, spin and trajectory evidence actually produced |
| Gaussian CPU | `sp`, `opt`, `ts`, `irc`, `td`, `link` | the same six job families | normal native outputs supplied by the user, including thermochemistry, excited states, spin and trajectory evidence |

The table describes implemented product paths, not the software installed on a
particular workstation. Before execution, ChemSmart still requires the exact
program/engine environment to be available and included in the human review.
Gaussian is separately licensed, GPU4PySCF requires a compatible CUDA stack,
and neither is implied by installing ChemSmart. NCIPLOT and remaining CLI job
families without an Agent engine/job declaration stay available to humans but
are not Agent execution paths in version 3.1.4.

The Runtime is provider-neutral. Version 3.1.4 ships registered adapters for
Alibaba Token Plan and DeepSeek. The user profile supplies the provider,
endpoint, model, reasoning setting, and credential label; CHEMSMART does not
choose a default model.

Across those programs, the public Agent surface includes:

- live capability and target-environment inspection;
- molecular identity/state binding and loader-validated project YAML;
- causal calculation DAGs, workflow-frontier inspection, independent branches,
  and validated optimized-geometry handoff;
- live CLI compilation, generated-artifact inspection, safe preview, and
  program preflight;
- semantic extraction from Gaussian and ORCA native output, receipt-bound xTB
  output, structured PySCF HDF5, and XYZ geometries or trajectories;
- RRHO and explicitly parameterised quasi-harmonic thermochemistry;
- unit-aware expressions for energy differences, CBS extrapolation,
  Boltzmann populations/averages, harmonic ZPE, imaginary-mode counts,
  distances, angles, dihedrals, centres of mass, inertia, rotor constants, and
  connectivity changes; and
- evidence-bound numerical claims and a human-readable scientific decision.

These layers are general ChemSmart operations. They do not encode a benchmark
answer, molecule-specific route, preferred DAG shape, or forced tool order.

## Installation

Create an isolated environment on the workstation or cluster that will host
the controller:

```bash
git clone https://github.com/Hongjiseung-ROK/chemsmart.git
cd chemsmart
conda env create -f environment.yml
conda activate chemsmart
python -m pip install .
chemsmart config
```

Program executables, compute interpreters, scheduler settings, scratch paths,
memory, and core counts belong in the user's ChemSmart configuration. Inspect
the generated files under `~/.chemsmart` before running an engine. Do not place
credentials in a project or calculation workspace.

Confirm the controller and preview path before enabling external programs:

```bash
chemsmart --version
chemsmart run --fake --no-scratch xtb \
  -p preview -f examples/xtb/water.xyz sp
```

## Project YAML and CLI

Project YAML is the canonical location for methods, basis sets, dispersion,
solvent, convergence, and stage settings. Supported CLI flags may override a
field for one invocation; CHEMSMART still validates and materialises the final
settings.

Example ORCA project:

```yaml
gas:
  functional: PBE0
  dispersion: D3BJ
  basis: def2-TZVP
  aux_basis: def2/J
  defgrid: DEFGRID3
  scf_tol: TightSCF
  freq: true
solv:
  functional: PBE0
  dispersion: D3BJ
  basis: def2-TZVP
  solvent_model: SMD
  solvent_id: water
  freq: false
```

Local execution and scheduler submission use the same program subcommands:

```bash
chemsmart run [RUN_OPTIONS] PROGRAM [PROGRAM_OPTIONS] JOBTYPE
chemsmart sub [SUB_OPTIONS] PROGRAM [PROGRAM_OPTIONS] JOBTYPE
```

Safe examples:

```bash
# Generate and inspect ORCA artifacts without invoking ORCA.
chemsmart run --fake --no-scratch orca \
  -p preview -f examples/xtb/water.xyz sp

# Generate a scheduler script without submitting it.
chemsmart sub --test --fake -s my_server orca \
  -p preview -f examples/xtb/water.xyz opt
```

Correct project YAML or CLI intent and regenerate when an artifact is wrong.
Do not hand-edit a generated native input and then treat it as ChemSmart state.

## Agent workflow

A plan session may create project YAML, compile commands, and safe previews but
cannot launch an engine:

```bash
chemsmart agent plan \
  --provider PROFILE \
  --task-file task.md \
  --secret-file /secure/path/provider.env \
  --workspace /absolute/path/agent-workspace
```

The terminal UI presents the planned DAG, molecular and electronic state,
project YAML, compiled commands, environment, resources, and exact digests.
The user may approve the complete workflow once, deny it, request revision, or
quit. An approved packet is executed separately and deterministically with no
provider call. See [Agent workflows](docs/source/agent-workflows.rst) for the
interactive and non-interactive commands.

## Scientific discipline

- Preserve identity, atom order, coordinate units, charge, multiplicity,
  electronic state, constraints, and physical conditions.
- Distinguish planned, previewed, approved, executed, parsed, validated, and
  interpreted work.
- Treat normal process exit as necessary but not sufficient evidence.
- Use validated optimized geometries for downstream calculations.
- Keep signs, units, standard states, temperature, pressure or concentration,
  and thermochemical conventions explicit.
- Cite CHEMSMART and every external program, method, basis set, and dataset used
  in a publication.

## Documentation

- [Documentation index](docs/source/index.rst)
- [Installation on Linux and macOS](docs/source/installation-linux-macos.rst)
- [HPC installation](docs/source/installation-hpc-cluster.rst)
- [CLI overview](docs/source/cli-overview.rst)
- [Project configuration](docs/source/configuration-project-settings.rst)
- [Server configuration](docs/source/configuration-server-settings.rst)
- [Agent workflows](docs/source/agent-workflows.rst)
- [PySCF](docs/source/pyscf-cli-options.rst)
- [xTB](docs/source/xtb-cli-options.rst)
- [ORCA](docs/source/orca-cli-options.rst)
