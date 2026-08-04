# ChemSmart PySCF/GPU4PySCF 24-Hour Goal

## Fundamental

ChemSmart is the canonical CLI and project-YAML control hub for computational
chemistry programs. Models express scientific intent and propose typed
settings. ChemSmart owns program grammar, generated scripts, command
compilation, permissions, execution, deterministic validation, evidence, and
terminal state.

## Objective

Develop ChemSmart 3.1.4 into a scientifically reliable PySCF/GPU4PySCF
control plane by completing and hardening PySCF SP/OPT/HESS, adding a bounded
closed-shell DFT TDA/TDDFT preview, unifying planning through execution in a
replayable scientific DAG, and improving a Qwen 3.8 Max specialist harness
through hypothesis-driven experiments and cross-case validation.

## Active window

- Activated: `2026-08-04T02:39:49+09:00`
- Deadline: `2026-08-05T02:39:49+09:00`
- Goal thread: `019fb8ec-85bf-7012-9aa0-d6c86d7cd2fd`
- Starting commit: `475f238e9b3bd7209cb78fbd678088e45b322d06`
- Working branch: `integration/chemsmart-v314-pyscf-agent`

The original development checkout is outside the implementation boundary. Its
dirty files must not be modified.

## Operating contract

- Use the named Alibaba Token Plan profile and exact observed
  `qwen3.8-max` model with `reasoning_effort=xhigh`.
- Do not use DeepSeek in this campaign.
- Give every model episode a unique hypothesis, changed factor, expected
  observation, deterministic oracle, and configuration record.
- Keep provider reasoning transient. Persist only public responses, typed tool
  traffic, decisions, artifacts, usage, and grader outcomes.
- Continue implementation and non-engine experiments during provider quota
  waits or while execution approval is absent.
- Never create meaningless requests solely to consume quota.
- Do not execute a chemistry engine without an existing exact approval.
  H2O and NH3 are the only molecules eligible for later approved execution.
- Never permit model-authored scripts, shell commands, readiness, approvals,
  or terminal outcomes.
- Retain a harness change only after a one-factor comparison and transfer
  check show an improvement without a deterministic safety regression.
- Create reviewable local commits. Do not push, open a pull request, publish,
  or claim reproduction, generality, GPU validation, or state of the art.

## Completion evidence

The final report must separate implemented behavior, deterministic checks,
live Qwen observations, rejected changes, unresolved defects, waiting
approvals, and unverified scientific claims. Goal completion requires a clean
accounting of every retained change and does not require an unapproved engine
run.
