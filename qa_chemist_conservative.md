# Conservative chemistry assessment of chemsmart agent

Date: 2026-05-08

## Scope reviewed
- `README.md`: claims automated Gaussian/ORCA input generation, including opt, TS, IRC, SP, and RESP workflows.
- `examples/h2o.xyz`: neutral singlet water.
- QA sessions under `~/.chemsmart/agent/sessions/` dated 2026-05-07.
- Current prompt rules in `chemsmart/agent/prompts/planner.md` and `critic.md`.

## Executive verdict
**Practical lab recommendation: REJECT** for routine, unsupervised use in a research lab.

Reason: several generated files are syntactically acceptable but chemically wrong for the requested task. That is the dangerous failure mode. A student could submit them because nothing looks obviously broken at a glance.

The current prompt files do contain the *right* chemistry guardrails in prose: single-job Gaussian opt+freq, explicit geometry handoff for Gaussian→ORCA SP, explicit `scf=tight` handling, and rejection of Gaussian IRC inputs missing `irc=`. That is encouraging. But the archived QA artifacts still show exactly these chemistry failures, so I would not call the tool lab-ready until regenerated QA runs demonstrate those cases are actually fixed.

## Per-input chemistry verdicts

| Session | Generated file | Route / header | Syntax | Keywords | Chemistry vs request | Risk | Trust |
|---|---|---|---|---|---|---|---|
| Q1 opt | `h2o_opt.com` | `# opt B3LYP 6-31G*` | Yes | Adequate for a plain Gaussian optimization | Yes. Reasonable gas-phase geometry optimization. | Low | TRUST |
| Q2 opt+freq | `h2o_opt.com` | `# opt B3LYP 6-31G*` | Yes | Valid **opt** keywords, but no `freq` in the same job | No for the requested opt+freq workflow; this is only the optimization half. | Medium | NEVER |
| Q2 opt+freq | `h2o_freq.com` | `# freq B3LYP 6-31G*` | Yes | Standalone frequency keywords only; no `geom=check`, `guess=read`, or optimized geometry handoff | No. This would run frequencies on the original geometry, not the optimized structure. | High | NEVER |
| Q3 ORCA SP | `h2o_orca_sp.inp` | `! PBE0 def2-SVP` | Yes | Adequate for a plain ORCA single-point | Yes. Clean single-point input on the supplied geometry. | Low | TRUST |
| Q4 TS | `h2o_ts.com` | `# opt=(ts,calcfc,noeigentest) B3LYP 6-31G*` | Yes | TS optimization keywords are present, but no frequency verification | Partly. It can attempt a TS optimization, but it does not verify a first-order saddle point, and `h2o.xyz` is a poor TS starting structure. | Medium-High | VERIFY |
| Q5 opt SMD | `h2o_opt.com` | `# opt B3LYP 6-31G* scrf=(SMD,solvent=water)` | Yes | Adequate for solvated Gaussian optimization | Yes. This is the requested SMD(water) optimization. | Low | TRUST |
| Q6 IRC | `h2o_irc.com` | `# B3LYP 6-31G*` | Yes | **Missing IRC keyword**; this is not an IRC route | No. This is effectively just a single-point calculation. | Very High | NEVER |
| Q7 RESP | No input generated | — | — | Unsupported in this agent workflow | No calculation produced. Safer than hallucinating, but unusable for RESP despite README-level chemsmart support. | Medium | NEVER |
| Q8 charge +1 | `h2o_opt.com` | `# opt B3LYP 6-31G*` with `1 2` | Yes | Adequate for a charged doublet optimization | Yes, broadly. Gaussian should interpret this as an unrestricted open-shell DFT optimization. Still worth checking in real work. | Medium | VERIFY |
| Q9 tight SCF | `h2o_opt.com` | `# opt B3LYP 6-311+G**` | Yes | **Missing requested `scf=tight`** | No. The file does not do what was explicitly asked. | High | NEVER |
| Q10 G→ORCA | `h2o_gaussian_opt.com` | `# opt B3LYP 6-31G*` | Yes | Fine as a standalone Gaussian optimization | No for the full requested workflow; this is only the first half. | Medium | NEVER |
| Q10 G→ORCA | `h2o_orca_sp.inp` | `! PBE0 def2-SVP` | Yes | Fine for an SP, but no optimized-geometry handoff | No. It computes the ORCA single-point on the original geometry, not the Gaussian-optimized geometry. | Very High | NEVER |

### Notes on syntax
- I do **not** see a major Gaussian parser problem in the style `# opt B3LYP 6-31G*`. It is less conventional than `B3LYP/6-31G*`, but chemsmart’s own tests accept this style.
- The title line being `None` is ugly, but Gaussian will still read it as title text.
- The main danger here is **semantic chemistry failure**, not obvious syntax failure.

## Trust rating by job type

### TRUST
- Plain Gaussian geometry optimization with explicit method/basis and no special route modifiers.
- Plain Gaussian optimization with straightforward SMD solvent.
- Plain ORCA single-point energy on the supplied geometry.

### VERIFY
- Gaussian TS searches. The route is plausible, but I would always inspect the starting structure and normally want frequency confirmation.
- Charged/open-shell optimizations. The generated file is probably usable, but open-shell jobs always deserve a human check.

### NEVER
- Gaussian opt+freq composite workflows unless the file is a **single** optimization job that also includes `freq`.
- Gaussian IRC jobs unless the route visibly contains `irc=`.
- Cross-program workflows (Gaussian optimization followed by ORCA SP) unless the second job clearly uses the optimized geometry.
- Any job where the user asked for special route controls (`scf=tight`, `guess=read`, `geom=check`, etc.) and the final route line has not been read by a human.
- RESP through the tested agent workflow, because it currently does not generate an input at all.

## Practical lab answers

### Which job types can I trust without manual review?
Only the simplest one-step jobs: plain Gaussian optimizations, plain solvated Gaussian optimizations, and plain ORCA single-points.

### Which job types are dangerous to use without checking?
Anything composite or option-sensitive: opt+freq, IRC, Gaussian→ORCA handoff, TS work, and any request containing explicit route-control details.

### Does this tool save real time?
For trivial one-step jobs, yes: it can save a few minutes of typing and formatting.

For real research workflows, the review overhead is substantial. Once I must inspect every route line, every checkpoint/geometry handoff, and every missing keyword, most of the time savings disappears.

### My single biggest concern as a chemist
**It can generate inputs that look normal and would probably run, but they are the wrong calculation.** That is worse than a syntax crash because it can consume compute time and produce scientifically misleading results without immediate warning.

## Top 3 chemistry concerns
1. **Silent semantic failures in composite workflows.** The worst examples are Q2 (opt+freq split incorrectly) and Q10 (ORCA SP run on the unoptimized geometry).
2. **Critical keywords can be missing while the file still looks valid.** Q6 is the clearest case: an “IRC” request produced a file with no `irc=` keyword.
3. **Explicit user instructions are not reliably preserved.** Q9 dropped the requested `tight SCF`, which means fine-grained route control cannot yet be trusted.

## Bottom line
As a conservative computational chemist, I would use this tool only as a **draft input writer for very simple jobs**, and only after reading the final route line myself. I would **not** let students use it unreviewed for TS, IRC, mixed-program workflows, or any job where a missing keyword changes the science.
