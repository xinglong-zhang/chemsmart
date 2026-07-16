# Review the current code changes

Perform a rigorous code review of the current changes.

Unless the user specifies another scope, review all staged and unstaged changes relative to the current branch base.

## Review workflow

1. Inspect repository instructions under `.cursor/rules/`.
2. Inspect the current Git status and complete diff.
3. Read the surrounding implementation, not only the changed lines.
4. Trace changed functions to their callers and downstream consumers.
5. Inspect relevant tests, CLI mappings, registries, configuration, exports, and documentation.
6. Run focused tests or static checks when they are available and reasonably scoped.
7. Do not modify files unless the user explicitly asks for fixes.

## Review priorities

Review in this order:

1. Incorrect behavior or logic errors
2. Regressions and backward-compatibility breaks
3. Data loss, corruption, or identity inconsistencies
4. Security and unsafe input handling
5. Error handling and hidden failures
6. Incorrect assumptions about external tools or file formats
7. Missing integration points
8. Edge cases and boundary conditions
9. Test quality and missing coverage
10. Performance and maintainability
11. Documentation inconsistencies

Ignore minor formatting preferences unless they cause ambiguity, violate an established project convention, or indicate a deeper problem.

## CHEMSMART-specific checks

When applicable, verify:

* CLI choices, argument forwarding, and dispatch mappings remain synchronized.
* New options are connected through the complete workflow.
* Gaussian, ORCA, xTB, and related program-specific behavior remain separated correctly.
* Program auto-detection does not silently misclassify outputs.
* Database identity fields and hashing inputs remain stable unless a deliberate migration is included.
* Existing schemas, record identifiers, and exported field names remain compatible.
* Scratch-directory behavior does not lose outputs or obscure debugging information.
* Batch and single-job execution paths behave consistently.
* Exceptions are not broadly suppressed.
* Tests do not depend unnecessarily on execution order or local environment state.

## PyMOL-specific checks

For PyMOL or molecular-rendering changes, verify:

* No atom indices, object names, residue numbers, or specific molecules are hard-coded.
* The implementation does not assume Mn is the only possible metal.
* Metal and coordination-core selections are restricted to the requested parent selection.
* Metal detection reuses the canonical project utility.
* Oxygen remains chemically recognizable as red and nitrogen as blue unless a documented style requirement overrides this.
* Multiple metal centers and structures without metals are handled safely.
* Repeated invocation does not accumulate pseudoatoms, labels, distance objects, or temporary selections.
* The style remains visually distinct from existing styles.
* CLI style choices, rendering dispatch, exports, documentation, and tests are all updated together.
* Optional PyMOL settings fail safely without hiding core programming errors.

## Finding format

Report only actionable findings.

For each finding, use:

### [Severity] Concise title

**Location:** `path/to/file.py:line`

**Problem:** Explain what is wrong.

**Impact:** Explain the concrete failure, regression, or maintenance risk.

**Recommendation:** Describe the smallest appropriate correction.

Use these severities:

* `Critical`: likely data loss, security issue, unusable feature, or broadly incorrect results
* `High`: significant functional defect or backward-compatibility break
* `Medium`: real defect affecting a limited case
* `Low`: minor but concrete correctness or maintainability issue

Do not classify optional enhancements or personal style preferences as defects.

## Final response structure

1. Findings, ordered by severity
2. Questions or assumptions
3. Tests and checks executed
4. Brief overall assessment

If no substantive findings are identified, say:

> No substantive correctness or regression issues were identified in the reviewed scope.

Then list any tests that were not run or behavior that could not be validated.

