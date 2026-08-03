---
name: chemsmart-scientific-workflow
description: Plan, command-compile, preflight, execute only with approval, or validate Gaussian, ORCA, xTB, and PySCF workflows in ChemSmart. Use for molecular identity, geometry, charge and multiplicity, method and basis selection, solvent or constraints, ScientificTaskSpec and CommandWorkflowSpec grounding, generated-artifact inspection, execution receipts, convergence, frequencies, units, and physical result checks.
---

# ChemSmart Scientific Workflow

Use this skill to turn a chemistry request into a bounded, reproducible
ChemSmart workflow. Read `AGENTS.md` first. A syntactically valid calculation
is not a chemically adequate result.

## Workflow

1. Make the scientific question and requested observable explicit. Separate a
   planned calculation from a completed result.
2. Identify the molecule from a stable artifact and geometry frame. Record
   coordinate units, charge, multiplicity, stereochemistry, fragments,
   constraints, and relevant conformer or spin assumptions.
3. Record these choices in ScientificTaskSpec. Select a program and job kind
   from the current CLI schema; keep method, basis/ECP, dispersion, solvent,
   temperature/standard state, and resources in an approved project artifact.
   Ask for any material missing choice.
4. Propose CommandWorkflowSpec data, not a native input or executable shell
   string. Let the deterministic compiler resolve live schema options, trusted
   project/artifact references, canonical argv, safe preview, and independent
   semantic round-trip checks.
5. Inspect ChemSmart-generated input only as downstream preview evidence. Do
   not edit .com, .gjf, or .inp files to repair a model proposal; repair the
   typed command intent or approved project YAML instead.
6. Obtain explicit approval before real local execution, scheduler submission,
   retry, cancellation, paid compute, or a material method change.
7. Record native inputs and outputs, executable and environment versions,
   hashes, timestamps, command, working directory, and termination state.
8. Parse values with units and apply job-specific validators before supporting
   a claim. Report missing diagnostics, assumptions, and uncertainty.

## Validation boundary

Check molecular identity, electron count, charge/multiplicity, method and
setting compatibility, SCF/geometry convergence, frequency/stationary-point
requirements, spin diagnostics where relevant, stoichiometry, comparability,
and units. A parser success, an exit code of zero, or a critic opinion alone
does not pass scientific validation.

At this roadmap stage, a fake/test CLI preview can establish only
`previewed`. An archived artifact may support `validated` if its deterministic
receipt is complete. Do not call a new calculation `executed`, `reproduced`,
or SOTA without the separately required engine, environment, and evidence
receipts.

## Use the references

- Read [scientific-task-contract.md](references/scientific-task-contract.md)
  for the required calculation specification and evidence receipt.
- Read [program-validation.md](references/program-validation.md) for
  engine-neutral checks and Gaussian/ORCA/xTB/PySCF boundaries.

## Examples

Use this skill for: “compile an ORCA transition-state intent into a safe
preview,” “check whether an xTB solvation request is fully specified,” or
“audit a Gaussian frequency result before reporting a free energy.”

Do not use it to hand-write Gaussian/ORCA/xTB native input or PySCF generated
scripts, change provider protocols, approve a run, or write a research claim without
`chemsmart-evidence-audit`.
