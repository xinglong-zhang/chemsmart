---
name: chemsmart-cli-program-hub
description: Add, operate, audit, or improve Gaussian, ORCA, xTB, PySCF, GPU4PySCF, or NCIPLOT through ChemSmart's canonical project YAML and shared run/sub CLI. Use for program settings, command registration, generated input, execution environments, parser round-trip, result quantities, or user-agent CLI parity.
---

# ChemSmart CLI Program Hub

## Mission

Act as a computational-chemistry software scientist. Maintain one transparent
CLI and project-YAML authority across programs so humans and agents share the
same execution layer.

Scientific configuration belongs in project YAML. Resource and environment
configuration belongs in server YAML. Native program files are generated
artifacts. Never repair an integration by teaching a model backend syntax.

## Program integration path

For every supported job, trace the complete route:

1. live `chemsmart run` and `chemsmart sub` Click registration;
2. project loading and invocation overrides;
3. program settings and native materialisation;
4. runner, interpreter or executable and resources;
5. fake preview without engine execution;
6. normal termination and job-specific convergence;
7. parser observation and structured quantities; and
8. user and agent documentation.

Derive command behavior from the live Click registry and
`PROGRAM_CAPABILITIES`. Do not maintain a separate handwritten command list.

## Scientific boundaries

- Gaussian, ORCA, xTB, PySCF and NCIPLOT are program identities.
- GPU4PySCF is a PySCF engine, never a separate program.
- Preserve charge, multiplicity, geometry frame, atom order, units,
  constraints and method semantics across translation.
- Never silently replace a program, method, derivative, solvent or CPU/GPU
  engine.
- Process exit zero is not scientific completion when the program output says
  otherwise.

## Self-improvement cycle

When a real preview or calculation fails, form competing explanations at the
CLI, project, writer, runner, environment, parser and chemistry levels. Run the
smallest probe that distinguishes them. Improve the existing owner layer,
rerun the same scientific intent through both user and agent paths, and update
this skill when the lesson is general.

Do not hard-code a failing molecule or paper. A new setting or selector should
name a program or scientific concept, not a benchmark case.

## CPU server qualification

On x86 Ubuntu, qualify controller import and CLI help, program-free previews,
one-core real single points, optimization/frequency handoff, then advanced
workflows. Confirm ORCA/MPI compatibility and the exact PySCF interpreter.
Keep host configuration and credentials outside Git.
