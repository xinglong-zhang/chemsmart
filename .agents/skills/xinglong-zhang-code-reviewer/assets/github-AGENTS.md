## Code Review Rules

For CHEMSMART code reviews, read and apply:
`.agents/skills/xinglong-zhang-code-reviewer/SKILL.md`

- Preserve scientific evidence and calculation meaning. Derive
  completion, energies, forces, and frequencies from parsed output.
  Retain units, spin/state identity, and complete linked-job routes.
  Report missing evidence explicitly; never fabricate success or
  silently discard scientific data.

- Preserve object ownership and usable APIs. Molecular operations
  belong on the molecular abstraction, parsed-file operations on the
  output/file abstraction, and input rendering on its writer. Reuse
  compatible owners and public APIs. Keep genuinely independent
  utilities as functions.

- Preserve shared Gaussian/ORCA contracts. Check analogous affected
  backends and identify concrete inconsistencies in behaviour or APIs.

- Establish usable instance state through the constructor or an
  explicit property contract. Reset per-run state before reusing a
  job. Identify the affected caller and consequence when reporting
  stale state or incomplete objects.

Support findings with inspected code and concrete consequences.
Separate correctness defects, OOP/API recommendations, and style
preferences. Respect the review host's priority threshold.

If the skill files are unavailable, report that limitation rather than
claiming to have applied them.

## Python class layout

Place __init__ after the class docstring and class attributes, before
ordinary helper definitions.

Group subsequent methods by responsibility, keeping properties and
public operations easy to find. Keep static helpers near related
private methods.

This is a readability convention. Do not classify method order as a
runtime bug or increase its priority to bypass a review filter.
Determine method binding separately: use self for instance state,
cls for class-level dispatch or construction, and static methods for
appropriate operations requiring neither.
