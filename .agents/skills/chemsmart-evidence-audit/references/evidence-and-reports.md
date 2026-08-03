# Evidence and report reference

## Canonical evidence bundle

Keep a structured manifest with stable IDs and SHA-256 digests for:

- scientific task specification, CommandWorkflowSpec, and geometry artifacts;
- CLI schema, canonical argv/display command, trusted project/artifact
  bindings, safe-preview output, independent parser observation, intent
  comparison, and structured counterexamples;
- native inputs, outputs, logs, parsed values, and units;
- program, executable, environment, and container/version information;
- commands, working directory, timestamps, exit status, and approvals;
- validation receipts, citation records, claim records, and review findings.

Use QCSchema-compatible calculation records and an RO-Crate-compatible bundle
manifest where practical. Keep engine-native artifacts alongside structured
views. Markdown, HTML, tables, and notebooks are derived outputs and must be
regenerable from the manifest.

## Claim taxonomy

| Claim | Minimum support |
| --- | --- |
| Observation | directly locatable artifact or log |
| Computed result | parsed value, unit, native output, and validation receipt |
| Inference | supported result plus stated assumptions and limits |
| Literature statement | verified citation and exact supported proposition |
| Unresolved uncertainty | explicit missing evidence or conflicting result |

Never collapse these labels in a conclusion.

For command preparation, call the state `previewed` only when every compiler
gate and the safe CLI preview are green. An archived artifact can be
`validated` only when the named deterministic receipt is complete. A report
must not imply engine execution, reproduction, or SOTA status from a command
or input file alone.
