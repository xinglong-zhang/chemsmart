# ChemSmart 3.1.4 integration status

## Purpose

This is the canonical status record for the PySCF 2.14.0, GPU4PySCF 1.8.0,
xTB 6.7.1, and capability-driven agent integration. The source baseline is
upstream ChemSmart revision `5d486b34501b9c040ad1db30d41a0eaa3850b15f`.
The package version remains `3.1.4`.

The statements below describe bounded integration observations. They are not
claims of general scientific reproduction, GPU execution, release readiness,
or state-of-the-art performance.

## Current milestone ledger

| Milestone | State | Current evidence | Remaining boundary |
|---|---|---|---|
| Program plane | implemented and live-observed | PySCF CPU `sp`, `opt`, and `hess` and xTB CPU `sp` executed through ChemSmart; GPU planning remained blocked on missing CUDA dependencies. | No GPU, Gaussian, ORCA, scheduler, or HPC execution was performed. |
| Agent plane | implemented and live-observed | The public `agent plan` and `agent run` paths reached the real DeepSeek session, typed tool host, command compiler, safe preview, approval-bound execution, result validation, and OPT-to-HESS handoff. | Generality beyond the bounded water cases is not established. |
| Focused validation | bounded checks green | The selected group produced 277 passing tests before the final fixture correction. The corrected xTB provenance test then passed as a targeted check, and read-only Ruff over the four corrected Python files was green. | The complete focused group was not rerun after the fixture correction, and the full repository suite was not run. |
| Local integration | committed release candidate | The program, agent, and test layers are recorded as reviewable local commits. Project-local skill packages were removed in a separate commit because maintained skills belong to the user's global registry. | Remote publication is established only by fetching the public ref after the release push; this document does not infer it from local state. |

## Live observations

### PySCF CPU workflow

One approval-enabled DeepSeek session completed the fixed workflow:

`SP(initial geometry) -> OPT(initial geometry) -> HESS(optimized geometry)`

The exact scientific settings were neutral singlet, gas-phase
B3LYP/def2-SVP, `defgrid2`, no density fitting, no dispersion, no solvent,
SCF tolerance `1e-9`, 100 SCF iterations, geomeTRIC with 100 optimization
steps, four CPU cores, 4 GB memory, and no scratch.

The host, rather than the model, compiled and launched all commands. Observed
results were:

- initial-geometry SP energy: `-76.35814945 Eh`;
- optimized energy: `-76.35831577 Eh`;
- optimized O-H distances: `0.9667046 Angstrom`;
- optimized H-O-H angle: `103.0930 degrees`;
- finite harmonic modes: `1638.6873`, `3791.3872`, and
  `3886.6168 cm^-1`;
- imaginary modes: zero.

The HESS input geometry matched the validated OPT HDF5 geometry. Requested
and applied settings agreed, atom identity and order were preserved, and the
three stages reported convergence. These observations validate this bounded
integration workflow, not an independent scientific reproduction.

### xTB CPU workflow

One approval-enabled DeepSeek session completed a gas-phase neutral-singlet
water GFN2 single point through the same generic execution tool. The run
terminated normally with finite energy `-5.070364532463 Eh` and preserved the
three-atom identity. The result audit now independently verifies receipt,
settings, state, executable, environment, source, project, output artifacts,
and finite energy. A forged status-only receipt is rejected.

### GPU4PySCF boundary

The live model created and validated GPU project and workflow planning
artifacts for `sp`, `opt`, and `hess`. Local environment observation found no
usable CuPy device or GPU4PySCF installation. Command materialization and
execution were therefore blocked, and no CPU fallback occurred. No GPU
calculation was run.

### Method-unspecified planning

The first development episode correctly refused to invent a functional or
basis and retained all nodes as unresolved, but it failed to preserve the
declared SP-to-OPT control order. A general prompt rule was added; the paired
episode then produced the ordered DAG
`SP -> OPT -> HESS` while remaining honestly blocked.

Both episodes exposed a deeper information defect: the model saw the
`ab_initio` field but not its loader-bounded value domain, and therefore
described MP2 and CCSD(T) as runnable alternatives. The capability registry
now exposes deterministic project-parameter domains, including
`ab_initio = hf` for the current PySCF backend. This final metadata change was
verified by direct registry projection but was not subjected to a third paid
model episode.

## Implemented authority boundary

- Models propose typed scientific intent, project sections, workflow nodes,
  explanations, and bounded repairs.
- The registry and live Click schema expose maintained commands and exact
  loader-bounded enumerations where available.
- Project loaders own scientific configuration validity.
- The command compiler owns argv and safe preview.
- Workflow approval owns exact nodes, resources, initial inputs, and typed
  producer edges.
- The execution host owns node order, idempotent in-session replay, fresh node
  directories, process launch, output collection, and terminal state.
- Program validators own convergence, provenance, finite values, electronic
  state, atom order, and optimized-geometry handoff.
- GPU requests never fall back silently to CPU.

## Validation boundary

Confirmed in the current tree or preserved local evidence:

- live CLI help and imports;
- CPU PySCF and xTB conformance and real execution observations;
- GPU-unavailable honest blocking;
- DeepSeek thinking-mode tool continuation;
- exact OPT-to-HESS geometry handoff;
- 277 selected focused tests passing;
- the corrected xTB provenance test passing in a targeted run;
- read-only Ruff passing on the four files corrected after the focused run;
- no raw provider response, calculation artifact, credential, or private
  reasoning file added to the repository;
- no case-insensitive forbidden integration-tool name in public content.

Not yet confirmed:

- a complete post-correction rerun of the focused test group;
- the full repository suite;
- GPU execution or CPU/GPU parity;
- broad chemistry or paper-level generalization;
- broad release readiness beyond the bounded checks and remote-ref verification.

## Canonical implementation map

| Responsibility | Source |
|---|---|
| PySCF CLI, projects, execution, and HDF5 | `chemsmart/cli/pyscf/`, `chemsmart/settings/pyscf.py`, `chemsmart/jobs/pyscf/`, `chemsmart/io/pyscf/` |
| xTB CLI, projects, execution, and result audit | `chemsmart/cli/xtb/`, `chemsmart/settings/xtb.py`, `chemsmart/jobs/xtb/` |
| Program and parameter-domain registry | `chemsmart/settings/capabilities.py` |
| Capability, environment, and conformance contracts | `chemsmart/agent/capabilities.py` |
| Project, command, preview, and preflight tools | `chemsmart/agent/projects.py`, `chemsmart/agent/commands.py`, `chemsmart/agent/preview.py`, `chemsmart/agent/preflight.py` |
| Approval, execution, and geometry handoff | `chemsmart/agent/execution.py`, `chemsmart/agent/tool_runtime.py` |
| Runtime V2 and DeepSeek continuation | `chemsmart/agent/runtime/`, `chemsmart/agent/services/`, `chemsmart/agent/live_session.py` |

## Publication gate

The release candidate is locally committed. Publication must bind a fresh
observation of the public fork's `main` ref, preserve its previous value in a
recoverable local bundle, use an exact lease for the replacement, and fetch the
remote ref again to prove local/remote equality. The full repository suite
remains explicitly unrun and must not be implied by a successful push.
