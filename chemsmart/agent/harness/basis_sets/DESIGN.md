# BSE-Backed Basis-Set Harness Design

## Goal

Validate whether a model/user basis-set phrase can be grounded to a concrete,
chemsmart-usable basis-set name before a project YAML or command reaches runtime.

This harness validates name intent and program compatibility. It does not vendor
or reinterpret orbital exponents; Basis Set Exchange remains the source of truth
for basis-set definitions.

## Data

`bse_basis_catalog.json` is generated from the installed
`basis_set_exchange` Python package by:

```bash
python -m chemsmart.agent.harness.basis_sets.build_bse_catalog
```

The catalog stores:

- `basis_sets`: canonical BSE entries with display name, family, role, element
  coverage, function types, auxiliaries, and renderable programs.
- `aliases`: normalized spelling variants to canonical BSE keys.
- `programs.gaussian`: basis names renderable through BSE `gaussian94`.
- `programs.orca`: basis names renderable through BSE `orca`.

As of the generated BSE 0.11 catalog, Gaussian and ORCA both expose 748
renderable basis names. They remain split by program so future BSE/software
drift can be handled without changing the harness interface.

## Harness Verdicts

- `ok`: concrete basis name resolved to a BSE canonical entry and is renderable
  for the target program.
- `warn`: multiple concrete basis names were found; caller should decide whether
  this is a mixed-basis request.
- `ask_user`: phrase is qualitative or family-level, not a concrete basis name.
  Examples: "good Karlsruhe triple-zeta basis", "large diffuse basis".
- `reject`: no concrete BSE-backed basis name can be found, or the name exists
  but is not renderable for the target program.

## Token-Efficient Search Tool

`search_basis_sets(query, program, limit=8, role="any")` is the model-facing
lookup tool. It must be used instead of passing the full BSE catalog into a
prompt.

The tool returns a compact top-k payload:

- normalized query and requested program/role
- `ok | warn | ask_user | reject` verdict
- at most `limit` candidates, clamped to 20
- candidate display name, family, role, match reasons, element-count, and a
  short auxiliary-basis map when available
- `token_policy: top_k_only; full catalog is never returned`

The ranking intentionally handles user language rather than exact strings only:

- family terms: Karlsruhe/Ahlrichs/def2, Dunning/cc, Pople/split-valence/6-31
- quality terms: double/triple/quadruple-zeta, diffuse/augmented, polarized/star
- role terms: JFIT, JKFIT, RIFIT, ADMMFIT
- common spoken forms: "six thirty one star" -> `6-31G*`, "tee zeta" -> TZ

If the result is `warn` or `ask_user`, the agent should preserve ambiguity and
show 2-4 candidates rather than silently choosing a basis.

## Integration Points

### Project YAML

Run basis checks on `gas` and `solv` blocks after YAML parse and before
`write_project_yaml`.

Simple basis:

```yaml
gas:
  functional: m062x
  basis: def2tzvp
```

Mixed basis:

```yaml
gas:
  basis: gen
  heavy_elements: [Br]
  heavy_elements_basis: def2svpd
  light_elements_basis: def2svp
```

For Gaussian, mixed basis must continue to use chemsmart's existing `gen` /
`genecp` semantics. The BSE harness validates `heavy_elements_basis` and
`light_elements_basis`; chemsmart runtime decides `gen` vs `genecp` from the
actual molecule.

### Command / Dry-Run Harness

After `dry_run_input`, inspect route lines:

- Gaussian route: simple basis token or `gen`/`genecp`.
- ORCA route: simple basis token in `!` route.

If route contains `gen`/`genecp`, inspect the settings/YAML evidence rather than
trying to infer basis details from the route line alone.

### Dataset Gate

Use the same resolver to reject training/eval records where the target contains
invented basis names, and mark qualitative basis requests as decline/ask-user
unless the record explicitly maps them to a concrete project policy.

## Non-Goals

- Do not choose "best" basis sets from qualitative user language.
- Do not replace BSE coefficient data.
- Do not infer element-specific ECP suitability purely from basis family names;
  use chemsmart runtime and BSE-backed GenECP generation for that stage.
