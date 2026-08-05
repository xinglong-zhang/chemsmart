# Paper-to-Workflow Campaign: Failure Report and Benchmark Reset

## Decision

The paper campaign no longer treats a structurally plausible tool trace as a
benchmark pass. A benchmark case is admissible only after an answer has been
prepared independently of the tested model.

Two primary evaluation modes are permitted:

1. **Small-molecule numerical mode.** ChemSmart executes the approved workflow,
   records host-rendered numerical claims, and compares every required value
   with a private reference value in compatible units and a preregistered
   absolute/relative tolerance.
2. **Canonical YAML–CLI DAG mode.** Intermediate or larger systems are not run.
   A computational-chemistry expert first prepares the correct ChemSmart
   project YAML documents, program-relative CLI job forms, exact input-geometry
   assignments, calculation dependencies, result extractions, mathematical
   analysis, and outputs. The model passes only when its host-normalized answer
   is semantically identical. Node names may differ; scientific content may not.

A case with neither a numerical answer nor a complete canonical YAML–CLI DAG
answer is not issued as a benchmark. It may be retained only as a diagnostic
pilot and cannot contribute a pass.

## Fe-aquo failure communicated to the campaign coordinator

The Fe-aquo P0 episode recovered the broad scientific skeleton of an HS/LS
DLPNO-CBS workflow, but it did not answer the research question.

- Assignment of the two coordinate artifacts to the high-spin and low-spin
  states depended on an assumption rather than separate approved identity
  records.
- The ORCA frequency stage was represented incorrectly. In ChemSmart's ORCA
  interface, harmonic frequencies are requested with `freq: true` in the
  optimization project and are produced by the `opt` job. ORCA does not expose
  a separate `hess` CLI job.
- Publication-level optimization settings were not fully materialized.
- The DLPNO basis-set calculations and CBS analysis were not executed.
- No value of the requested high-spin/low-spin energy difference was produced.

The scientifically accurate outcome is therefore **planned / partially
previewed, final observable blocked**. The historical structural oracle marked
this episode strict-pass because it found the expected tool names and terms.
That label is retired as a primary result and retained only as a diagnostic
record.

The public Fe P0 receipt records 72 provider attempts, 17,707,986 input tokens,
192,318 output tokens, 64,228 reasoning tokens, 127 successful tool calls, five
failed calls, and terminal state `planned`. P1 was interrupted after the
evaluation rule changed; P2 was not dispatched. Neither is a benchmark
observation.

## Earlier campaign observations under the revised interpretation

### QM8 methane development case

The P0/P1/P2 runs recovered the requested excited-state workflow vocabulary and
exposed a real result-tool defect: the model-visible schema omitted registered
Gaussian/PySCF excitation-energy and oscillator-strength selectors. The generic
selector surface was corrected. A paired P2 rerun then retained those
observables, but provider input rose to 3,071,241 tokens. This is evidence of a
tool-surface repair, not an efficiency gain or a numerical benchmark pass.

### QM40 transfer case

All three arms preserved the main B3LYP/6-31G(2df,p) optimization/frequency
intent and kept unsupported LModeA analysis visible. P1 used 1,832,605 input
tokens; P2 used 4,801,660, so dependency projection did not satisfy the planned
efficiency criterion. P0 also produced an incorrect molecular formula in free
prose despite retaining the typed atom order. No numerical observable was
computed. These runs are planning diagnostics only.

The associated sources are the [QM8 article](https://doi.org/10.1063/1.4928757),
the [QM40 article](https://doi.org/10.1038/s41597-024-04206-y), and the
[Fe-aquo DLPNO-CBS article](https://doi.org/10.1021/acs.jctc.9b01109).

## Implemented general changes

### Answer-first benchmark contracts

`paper_workflow_campaign.py` now provides separate private contracts for:

- numerical reference quantities and tolerances;
- small-molecule numerical answer keys;
- expert-prepared project-YAML and workflow-DAG answer keys;
- pre-dispatch eligibility;
- numerical claim grading; and
- semantic YAML–CLI DAG grading.

The public black-box prompt contains none of the reference values or answer
structures. The benchmark-dispatch function cannot render an eligible task
without one of the two answer-key types.

The former `PaperWorkflowGradeV1.strict_pass` remains available for debugging
tool coverage, but its role is explicitly `diagnostic_structure_only`.

### Program-relative DAG comparison

Canonical workflow grading is independent of arbitrary model node names. It
compares:

- project program and exact normalized YAML contents;
- target program and ChemSmart CLI job form;
- exact initial geometry assignment;
- producer output to consumer input edges;
- deterministic analysis operations; and
- declared output semantics.

This means an ORCA `opt` project with `freq: true` can match the expert answer,
while a superficially similar `opt → hess` ORCA DAG fails. A PySCF workflow may
correctly use a separate `hess` node because the reference is program-relative.

### Reusable postprocessing and identity support

The typed dimensional expression engine now exposes two general component-wise
CBS operations:

- two-point exponential SCF extrapolation with explicit basis cardinals and
  exponent; and
- two-point inverse-power correlation extrapolation with explicit basis
  cardinals and exponent.

The operations require separately parsed SCF and correlation energies. They do
not encode the Fe-aquo answer, assume a particular correlated method, or accept
model-authored formulas.

Live paper sessions can now register multiple approved molecular identities,
each bound to its own coordinate bytes and atom order. This permits distinct
state geometries without inferring identity from a filename or XYZ comment.

## Verification and remaining work

A focused suite covering the answer contracts, unit-aware numerical grading,
program-relative ORCA workflow comparison, two-point CBS expressions,
multi-geometry identity, result readers, analysis nodes, and dependency context
passed **79/79** in the ChemSmart environment. An earlier launch using the base
Anaconda interpreter stopped before collection because optional PyMOL was not
installed there; that was an environment-selection error, not a code failure.

The following remain unverified:

- automatic conversion of every live Runtime V2 workflow into
  `PaperWorkflowObservationV1`;
- a small-molecule paper run that reaches all reference values;
- expert preparation and independent review of a complete Fe-aquo YAML–CLI DAG
  answer;
- transfer of the new CBS operations to a second paper; and
- any claim of paper reproduction, generality, production readiness, or SOTA
  performance.

No further provider call should be made until the next case has either a
complete numerical answer key (small molecule) or a complete expert-authored
canonical YAML–CLI DAG answer (intermediate or larger system).
