# Scientific task contract

Define a `ScientificTaskSpec` before the agent can plan execution. At minimum,
record:

- stable molecule and geometry-frame identifier plus a content digest;
- coordinate units, charge, multiplicity, fragments, stereochemistry, and
  constraints;
- requested observable and whether the task is a prediction, optimization,
  frequency, pathway, comparison, or analysis;
- program, job kind, method, basis/ECP, dispersion, solvent, temperature,
  standard-state convention, and resource target;
- expected artifacts, validation rules, uncertainties, and decision boundary.

Derive a CommandWorkflowSpec from the task specification. It identifies the
CLI-schema digest and immutable command nodes with a Click path, canonical
parameter map, trusted project reference, ArtifactBinding IDs/hashes,
charge/multiplicity, execution intent, dependencies, and expected artifact
classes. The model must not select a filesystem path, shell spelling, option
order, or native-engine input text.

Use ChemSmart-generated native input and output files as authoritative engine
records after safe preview or execution. A QCSchema-compatible structured
record should point to those artifacts rather than erase code-specific
settings. A CommandPreflightReceipt links the typed IR, canonical argv, parser
observation, intent comparison, resolved project/input/environment, and
safe-preview artifact hash.

## Receipt minimum

Record command, working directory, input/output hashes, program and executable
version, environment identity, start/end timestamps, exit status, parsed
quantities and units, validator versions/results, approval reference, and
claim-to-evidence links.

State a numerical convention before comparing results: energy reference,
temperature, standard state, unit conversion, and whether the value is an
electronic energy, enthalpy, Gibbs energy, barrier, or derived quantity.
