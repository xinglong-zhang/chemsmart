# Conservative chemistry reassessment of chemsmart agent (v2)

Date: 2026-05-08
Target code reviewed: `fork/main` at `ee3bbd470331` plus one follow-up local fix on this branch for clearer mixed-workflow failure handling.

## What I re-checked
1. Current `fork/main` code and tests for the four claimed fixes:
   - PR #14: Gaussian IRC jobs promote to `GaussianIRCJobSettings` and write `irc=(...)`.
   - PR #17: Gaussian opt+freq is a **single** `gaussian.opt` job with `freq=True`.
   - PR #15: mixed Gaussian→ORCA workflow includes `extract_optimized_geometry`.
   - PR #16: critic inspects **all** generated dry-run inputs.
2. Fresh `chemsmart agent run --dry-submit` reproductions for the four cases I previously rated NEVER:
   - Q2 opt+freq
   - Q6 IRC
   - Q9 tight SCF
   - Q10 Gaussian→ORCA
3. One remaining operational failure in Q10: if local Gaussian execution fails, the agent previously crashed forward into `extract_optimized_geometry`. I fixed that so it now stops cleanly and reports `run_local failed ...`.

## Verification of the claimed fixes in current code
- **PR #14 confirmed**
  - `tests/agent/test_agent_build_job.py` asserts `"irc=(" in job.settings.route_string` for `gaussian.irc`.
  - Live dry-run now writes `# B3LYP 6-31G* irc=(calcfc,stepsize=10,maxpoints=20)`.
- **PR #17 confirmed**
  - `tests/agent/test_agent_build_job.py` asserts one Gaussian opt job with both `opt` and `freq` in one route line when `freq=True`.
  - Live dry-run now writes `# opt freq B3LYP 6-31G*` in one `.com` file.
- **PR #15 confirmed**
  - `extract_optimized_geometry` exists in `chemsmart/agent/tools.py` and the registry.
  - Live Q10 plans now build the ORCA SP job from the extracted optimized geometry step, not from the original `build_molecule` step.
- **PR #16 confirmed**
  - `chemsmart/agent/core.py` sends `dry_run_inputs` to the critic.
  - `tests/agent/test_agent_planner_critic_mock.py::test_critic_receives_all_dry_run_inputs` passes.

## Fresh dry-submit results for the previously NEVER cases

| Case | Current generated file(s) | Current route line(s) | Chemistry verdict now | Current trust |
|---|---|---|---|---|
| Q2 opt+freq | `h2o_opt_freq.com` | `# opt freq B3LYP 6-31G*` | Correct fix. This now does what was requested in one Gaussian job. | TRUST |
| Q6 IRC | `h2o_irc.com` | `# B3LYP 6-31G* irc=(calcfc,stepsize=10,maxpoints=20)` | The missing IRC keyword bug is fixed. Route is now an actual IRC route. For real chemistry I would still verify that the starting structure is a TS, not a minimum. | VERIFY |
| Q9 tight SCF | `h2o_opt.com` | `# opt B3LYP 6-311+G** scf=tight` | Correct fix. The requested SCF tightening is now preserved. | TRUST |
| Q10 Gaussian→ORCA | `h2o_gaussian_opt.com` under default dry-submit | `# opt B3LYP 6-31G*` | The **wrong-geometry handoff bug is fixed in the plan**, but the final ORCA `.inp` is not generated during default dry-submit in this environment because the workflow correctly waits for a successful Gaussian optimization and extracted geometry. | VERIFY |

## What is still unresolved?
### Chemistry/input-generation issues from my original NEVER list
- **Q2 opt+freq:** resolved.
- **Q6 IRC missing keyword:** resolved.
- **Q9 tight SCF dropped:** resolved.
- **Q10 original-geometry ORCA SP handoff:** resolved at the planning/dataflow level.

### Remaining practical issue
The mixed Gaussian→ORCA case is **not chemistry-wrong anymore**, but it is still harder to validate end-to-end in dry-submit mode because the second input depends on a completed upstream optimization.

In this environment, `--allow-remote-unknown` reached `run_local`, but local Gaussian did not produce an output log. Before my patch, the agent then crashed forward into `extract_optimized_geometry`. After my patch, it now stops cleanly with:

- `Error: run_local failed with returncode 1; see .../h2o_gaussian_opt.stderr`

That is a workflow robustness fix, not a chemistry-route fix.

## Code change made during this reassessment
I made one focused fix for the remaining Q10 execution-path failure:
- if `run_local` returns `ok: false`, the agent now logs a tool error and stops **before** calling `extract_optimized_geometry`
- the CLI now wraps that runtime failure as a clean `ClickException` instead of a Python traceback

### Tests added/updated
- `tests/agent/test_agent_session_metadata.py::test_run_local_failure_stops_before_geometry_extraction`
- `tests/agent/test_agent_cli_run.py::test_agent_cli_run_wraps_runtime_errors_as_click_exceptions`

### Validation run
- `pytest -v tests/agent/test_agent_session_metadata.py tests/agent/test_agent_cli_run.py tests/agent/test_agent_build_job.py tests/agent/test_agent_planner_critic_mock.py`
- **22 passed**

## Updated conservative lab recommendation
**Practical recommendation: CONDITIONAL ADOPT**

Compared with my first report, the most dangerous chemistry mistakes I flagged have been substantially repaired:
- opt+freq is now one correct Gaussian job
- IRC now actually contains an IRC route keyword
- explicit `scf=tight` is preserved
- mixed Gaussian→ORCA workflows now attempt geometry handoff from an optimized structure rather than reusing the original geometry

## What I would trust now
### TRUST
- Plain Gaussian optimizations
- Gaussian opt+freq when the route visibly contains both `opt` and `freq`
- Plain solvated Gaussian optimizations
- Plain ORCA single-point jobs
- Simple route-option preservation cases like `scf=tight`

### VERIFY
- IRC jobs: keyword generation is fixed, but the chemist must still confirm the starting point is a TS structure suitable for IRC.
- TS jobs in general
- Charged/open-shell jobs
- Mixed Gaussian→ORCA workflows, because the downstream input cannot be reviewed until the upstream optimization actually succeeds.

### STILL NOT READY FOR BLIND TRUST
- Any chained workflow where downstream files depend on real upstream outputs that have not yet been produced.
- Unsupported workflows such as RESP through the current agent path, unless separately demonstrated.

## My updated single biggest concern as a chemist
The major “looks normal but computes the wrong thing” failures are much improved. My biggest concern now is **workflow transparency for chained jobs**: for mixed-program sequences, I cannot always inspect the final downstream input until the previous step has really run and produced a usable output structure.

## Bottom line
My original **REJECT** verdict is no longer appropriate for the four re-tested NEVER cases. Three are fixed outright, and the fourth (Gaussian→ORCA) is chemically corrected in its handoff logic but remains operationally dependent on a successful upstream local optimization.

So my updated conservative stance is:
- **simple one-step jobs:** usable
- **opt+freq:** now usable
- **IRC/TS/mixed workflows:** still human-reviewed, but no longer automatically disqualified for the same reasons as before
