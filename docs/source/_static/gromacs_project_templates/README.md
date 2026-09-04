# gmxsmart literature → project.yaml seed dataset v0.1

This is a starter training/benchmark dataset, **not a preset library and not a set of runnable simulations**.

## Core representation

`paper/SI section → normalized MD parameters → project.yaml → unsupported/schema-gap fields`

A partial but traceable pair is more valuable than a complete-looking YAML containing guessed values.

## Contents

- `literature_mapping.xlsx` — source tracker + templates + field mappings + schema gaps + queue.
- `templates/` — one training-target YAML per MD stage.
- `training_pairs.jsonl` — compact machine-readable pairs.
- `example_walkthrough.md` — one case mapped end-to-end.
- `schema_gap_proposal.md` — recurring literature parameters the current schema cannot express.
- `validate_templates.py` — basic YAML/key validation; optionally invokes ChemSmart's parser in-repo.

## Collection rule

Never invent a missing MD parameter. If the source does not state it, omit it and record the missing/unsupported field.  
Do not invent local `structure`, `topology`, or other input paths: those belong to a runnable user project, not to the literature training target.

## Search strings

- `site:rsc.org/suppdata GROMACS "Supplementary Information"`
- `site:acs.figshare.com GROMACS "Supporting Information"`
- `"GROMACS" "TIP3P" "Parrinello-Rahman" "Supporting Information"`
- `"GROMACS" "OPLS-AA" "2 fs" "1 bar" "Supplementary Information"`
