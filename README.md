# CHEMSMART

[![CI](https://github.com/Hongjiseung-ROK/chemsmart/actions/workflows/main.yml/badge.svg)](https://github.com/Hongjiseung-ROK/chemsmart/actions/workflows/main.yml)

<p align="center">
  <img src="docs/source/_static/chemsmart_logo.png" alt="CHEMSMART logo" width="600">
</p>

CHEMSMART is a CLI-first, project-YAML hub for computational chemistry. It
gives users and research agents one maintained interface for preparing,
running, submitting, and analysing calculations across supported quantum
chemistry programs.

The same command hub also provides molecular-database assembly, inspection,
query, and export operations, together with ITERATE workflows for systematic
structure generation and calculation campaigns.

Scientific choices live in readable project YAML. CHEMSMART translates those
choices into program-native inputs, controls execution through the same
`chemsmart run` and `chemsmart sub` commands used by human researchers, and
returns structured quantities for reproducible post-processing. A model never
needs to invent Gaussian, ORCA, xTB, or PySCF input syntax.

Current package version: **3.1.4**. Python **3.10** is required.

## Architecture

```text
scientific question
       │
       ▼
project YAML + molecular artifact
       │
       ▼
CHEMSMART CLI and scientific DAG
       │
       ├── Gaussian / ORCA / xTB native input
       ├── PySCF standalone compute script
       └── local runner or scheduler submission
       │
       ▼
program output + structured result reader
       │
       ▼
unit-aware analysis, thermochemistry, and claims
```

The CLI is the operational authority for both users and agents. Generated
native files are downstream execution artifacts, not an alternative user or
model API.

## Program support

| Program | Maintained CHEMSMART surface | Execution requirement |
| --- | --- | --- |
| Gaussian | input generation, local/scheduled jobs, result parsing and analysis | licensed Gaussian installation |
| ORCA | SP, OPT, TS, IRC, scan, NEB, TD/TDA, QRC, pKa and analysis workflows | compatible ORCA installation and MPI runtime when parallel |
| xTB | CPU `sp`, `opt`, and `hess` with optional strict project YAML | xTB executable |
| PySCF | CPU `sp`, `opt`, `hess`, structured HDF5 results; bounded TD preview | dedicated PySCF 2.14 compute interpreter |
| GPU4PySCF | PySCF execution engine, not a separate program | compatible NVIDIA driver, CUDA, CuPy, cuTENSOR and GPU4PySCF stack |
| NCIPLOT | input generation and execution integration | NCIPLOT installation |

The Python package does not distribute external program binaries. The private
CPU research image described below intentionally includes the redistributable,
pinned PySCF and xTB runtimes; it does not include Gaussian, ORCA, GPU4PySCF,
or CUDA.

## Private Ubuntu x86 CPU research image

The reference CPU-server environment is the private ``linux/amd64`` image:

```text
ghcr.io/hongjiseung-rok/chemsmart:3.1.4-cpu
```

It contains ChemSmart 3.1.4, Python 3.10, PySCF 2.14.0, xTB 6.7.1,
geomeTRIC 1.1.1, pyscf-dispersion 1.5.0, NumPy 1.26.4, h5py 3.16.0,
PyMOL open source 3.1.0, Open Babel, and the locked runtime dependencies.
GPU4PySCF, CUDA, Gaussian, and ORCA are intentionally excluded.

Authenticate to private GHCR with a token carrying ``read:packages`` and run
the non-root image with explicit workspace, scratch, temporary storage, and
xTB stack limits:

```bash
printf '%s' "$GHCR_READ_TOKEN" | \
  docker login ghcr.io -u Hongjiseung-ROK --password-stdin
docker pull ghcr.io/hongjiseung-rok/chemsmart:3.1.4-cpu

mkdir -p "$PWD/workspace" "$PWD/scratch"
docker run --rm -it --read-only --ulimit stack=-1:-1 \
  --cpus 8 --memory 32g --tmpfs /tmp:rw,size=4g \
  --mount "type=bind,src=$PWD/workspace,dst=/workspace" \
  --mount "type=bind,src=$PWD/scratch,dst=/scratch" \
  ghcr.io/hongjiseung-rok/chemsmart:3.1.4-cpu \
  chemsmart --help
```

See [the CPU image guide](containers/cpu/README.md) for exact package pins,
CLI examples, read-only operation, configuration and secret mounts, image
qualification, and promotion policy.

## Source installation on an Ubuntu x86 CPU server

The recommended benchmark host is an x86-64 Ubuntu server with Conda or
Mamba. Start from a clean clone of the maintained fork:

```bash
sudo apt-get update
sudo apt-get install -y build-essential ca-certificates curl git \
  libegl1 libgl1 libsm6 libxext6 libxrender1

git clone https://github.com/Hongjiseung-ROK/chemsmart.git
cd chemsmart
make env
conda activate chemsmart
make install-dev
make configure
```

`make env` creates or updates the `chemsmart` environment from
`environment.yml`. `make configure` creates local templates under
`~/.chemsmart`; inspect those files before running any engine. Program paths,
scheduler directives, scratch paths, memory, and core counts are machine-local
configuration and must not be committed.

Verify the controller before installing external engines:

```bash
chemsmart --version
chemsmart --help
chemsmart run --help
chemsmart agent --help
```

Run a program-free xTB preview through the real CLI schema:

```bash
chemsmart run --fake --no-scratch xtb \
  -p test -f examples/xtb/water.xyz sp
```

For private-container deployment, source installation, PySCF/xTB execution,
and staged benchmark work, see
[Ubuntu CPU server installation](docs/source/installation-ubuntu-cpu-server.rst).

## Project YAML

Project YAML is the canonical location for method, basis, dispersion, solvent,
convergence, and stage settings. CLI options may override supported fields for
one invocation, but CHEMSMART remains responsible for validating and
materialising the final program settings.

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

Example stage-specific PySCF CPU project:

```yaml
sp:
  ab_initio: hf
  functional: null
  basis: def2-svp
  density_fit: false
  freq: false
opt:
  ab_initio: hf
  functional: null
  basis: def2-svp
  opt_solver: geometric
  opt_maxsteps: 100
  freq: false
hess:
  ab_initio: hf
  functional: null
  basis: def2-svp
  freq: true
```

Never edit a generated native input to make a project appear valid. Correct
the YAML or CLI intent and regenerate the artifact.

## CLI usage

Local execution and scheduler submission share program subcommands:

```bash
chemsmart run [RUN_OPTIONS] PROGRAM [PROGRAM_OPTIONS] JOBTYPE
chemsmart sub [SUB_OPTIONS] PROGRAM [PROGRAM_OPTIONS] JOBTYPE
```

Examples:

```bash
# Local ORCA preview; does not invoke ORCA
chemsmart run --fake --no-scratch orca \
  -p test -f examples/xtb/water.xyz sp

# Real local xTB calculation after xTB is configured
chemsmart run --no-scratch xtb \
  -p test -f examples/xtb/water.xyz opt

# PySCF fixed-geometry Hessian in a configured CPU compute environment
chemsmart run --no-scratch pyscf \
  -p test -f examples/xtb/water.xyz hess

# Generate and inspect a scheduler script without submitting it
chemsmart sub --test --fake -s my_server orca \
  -p test -f examples/xtb/water.xyz opt
```

Use `--help` at the exact command level to inspect the live option surface.

## Agent workflows

The provider-neutral agent uses the same YAML loaders, CLI compiler, program
adapters, result readers, and scientific analysis operations as a user.

```bash
chemsmart agent plan \
  --provider PROFILE \
  --task-file task.md \
  --secret-file /secure/path/api.env \
  --workspace ./agent-workspace
```

`agent plan` may create project YAML and safe previews but does not run a
chemistry engine. `agent run` stays preview-only unless an approval file is
provided. `agent execute` performs an already approved workflow without a
model or provider credential.

The repository includes a focused set of computational-research skills under
`.agents/skills`. They teach future research agents to work as challenging,
self-improving computational scientists while keeping program control inside
CHEMSMART. They do not encode benchmark answers or replace the live CLI.

See [Agent workflows](docs/source/agent-workflows.rst) for the operational
boundary and the behavior-first research loop.

## Scientific execution discipline

- Preserve molecular identity, atom order, coordinate units, charge,
  multiplicity, electronic state, and constraints.
- Distinguish planned, previewed, executed, parsed, and scientifically
  interpreted work.
- Treat normal process exit as necessary but not sufficient; check
  convergence and the requested properties.
- Use validated optimized geometries for downstream Hessian, SP, IRC, or path
  calculations.
- Keep signs, units, thermochemical conventions, temperature, pressure or
  concentration explicit.
- Prefer a new scientific experiment over a benchmark-specific rule. Improve
  the smallest reusable CLI, YAML, parser, or analysis capability exposed by a
  real failure.

## Development

```bash
conda activate chemsmart
make install-dev
make pre-commit
make test
make docs
```

Focused tests should be run while developing. Full tests and documentation
builds should pass before a release branch is published. See
[CONTRIBUTING.md](CONTRIBUTING.md).

## Documentation

- [Documentation index](docs/source/index.rst)
- [Ubuntu CPU server installation](docs/source/installation-ubuntu-cpu-server.rst)
- [Private CPU research image](containers/cpu/README.md)
- [CLI overview](docs/source/cli-overview.rst)
- [Project configuration](docs/source/configuration-project-settings.rst)
- [Server configuration](docs/source/configuration-server-settings.rst)
- [Agent workflows](docs/source/agent-workflows.rst)
- [PySCF](docs/source/pyscf-cli-options.rst)
- [xTB](docs/source/xtb-cli-options.rst)
- [ORCA](docs/source/orca-cli-options.rst)

## Citation

If CHEMSMART contributes to a publication, cite the toolkit and the external
programs, methods, basis sets, and datasets used by the workflow.

Zhang, X.; Tan, H.; Liu, J.; Li, Z.; Wang, L.; Chen, B. W. J.
*CHEMSMART: Chemistry Simulation and Modeling Automation Toolkit for
High-Efficiency Computational Chemistry Workflows*. arXiv 2025,
arXiv:2508.20042. <https://doi.org/10.48550/arXiv.2508.20042>.

## License

CHEMSMART is distributed under the terms in [LICENSE](LICENSE).
