---
name: chemsmart-evidence-audit
description: Perform a bounded human-readable audit of ChemSmart calculations, result analyses, scientific claims, citations, or release documentation. Use when the user explicitly asks whether a reported quantity is supported by program output, units, conditions, execution state, or source literature; do not use as a default benchmark gate.
---

# ChemSmart Evidence Audit

## Mission

Act as a computational scientist checking whether a claim says what the
calculation actually established. Use the ordinary ChemSmart artifacts and
program semantics; do not build an additional control system.

Apply this skill only when an audit is requested or a scientific claim depends
on ambiguous evidence. It is not a prerequisite for every refine-loop trial.

## Audit method

1. Classify the statement as observation, computed quantity, inference,
   literature comparison or unresolved uncertainty.
2. Identify the molecular identity, geometry frame, state, method, program,
   physical conditions and source result.
3. Verify whether the work was planned, previewed, executed, parsed or
   reproduced.
4. Recompute signs, dimensions, units and simple algebra independently.
5. Check normal termination, convergence and requested properties when a real
   result is claimed.
6. Read the scientific answer as a whole and recognise equivalent valid routes.
7. Report consequential support gaps and proportional limitations in plain
   language.

Do not use exact DAG identity, tool count, tokens, hash scores or runtime labels
as the scientific verdict. Stable IDs may locate ordinary artifacts but do not
replace chemical reasoning.

## Self-improvement cycle

When an audit exposes a recurring general ambiguity, improve the owning result
reader, quantity definition, unit operation, documentation page or skill. Test
the smallest mechanical change and challenge it with another scientific case.
Do not convert a private answer or audit rubric into a production gate.

Revise this skill when a real audit shows that it overstates evidence,
suppresses a correct alternative or misses a material scientific convention.

## Report

Lead with the supported scientific conclusion. Separate strengths,
consequential gaps, minor caveats and unresolved questions. State what would
resolve uncertainty without implying that more paperwork is better science.
