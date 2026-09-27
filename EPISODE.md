# R11 episode behav -- the Agent as scientist, the LLM side (lens: behaviour)

Base SHA: 9185770e30e486e0d1e96e83bf35edc914b8590a (verified with
`git rev-parse HEAD` as the first action, 2026-09-28). Researcher model:
claude-opus-5-5[1m]. Episode id `behav`. Brief sha256 prefixes:
brief-behav.md acefe5b7a612a264, common-behav.md ed89ca7e461f8e3d.

## The question (as currently understood)

Under what conditions does the CHEMSMART Agent generate, notice, pursue or
dissent from a hypothesis without being scripted toward the observation,
and which of those conditions live in the model, which in its context and
tools, and which in the loop's moments? Serves claim C4 of the round map
(the Agent's science moves with the acts and moments it is given, not the
messages it is shown -- beyond one model). Phase I's part of it: is R10's
record (one model, message levers null at ceiling) a fact about the
surface or about deepseek-v4-flash-0731? That needs the same archived
decision points answered by two models, with outcomes the host types.

## What the base tree does (read before any sample, 2026-09-28)

- The brief's reading-turn prior holds: `_system_prompt` renders
  `wake.goal_authority` whenever a goal context exists
  (live_session.py `_system_prompt`, execution_sentence), and
  `driver._reading_context` hands the reading turn a goal context whose
  own authority says budgets are zero and nothing it plans is decided;
  `wake.reading_turn` names inspect_run (core), extract_result_quantities,
  record_analysis_claims and record_scientific_decision (all three
  deferred under host_search; `promotions_from_states` pins nothing for a
  certified delivery). Not changed here: Q29 measured this version.
- Mid-session host messages under a goal: the wave notice, the
  termination notice and the goal terms re-injected every 5 tool turns
  (`_HOST_REINJECTION_TURN_INTERVAL`). None recorded what was callable.
  Cycle 1's authority (`_goal_terms_context`) names
  record_scientific_decision, deferred under host_search.
- A goal-cycle-1 context names each geometry by content id with atom
  count and symbols, never by file name; the task text carries the names.
- `alibaba-token-plan` exposes the catalogue by `host_search` (8 core
  tools: search_capabilities, inspect_program, project_yaml,
  bind_scientific_identity, plan_scientific_workflow,
  inspect_workflow_frontier, compile_command, inspect_run).
- The model sees its own provider record (profile name, model id,
  reasoning effort, enable_thinking, deadlines) in its context.

## Commits so far

- d4e63923 loop: a message the host places in a goal session records the
  calls the model could make when it read it (deliverable 2; witness red
  on a pristine export of the base, green here).
- 1aec4e37, 5e7212fb research: the matched-turn instrument
  `.agents/research/loop/matched_turns.py` (deliverable 1) and its two
  corrections (the credential is never named; a failed real attempt keeps
  its turn).
- 49eeaf8d research: `.agents/research/loop/matched_outcomes.py`, the
  pre-registered outcome rules below.

## Corrected premises so far

- "A second model is a provider-profile edit": true only with an output
  cap the model accepts. qwen3.8-max on the Token Plan answered HTTP 400
  (request_invalid) to a max_output_tokens of 200,000 (the cluster
  deepseek profile's value) with reasoning effort xhigh and with max; it is
  served at 65,536 (leased probe, observed model qwen3.8-max).
- A served response must echo the requested model (DeepSeekProtocolError,
  failing_field "model"); replayed prefix turns cannot carry a label, so
  the instrument tells replayed, real and closing turns apart by position.
- The round guard refuses any command naming the credential file; a real
  sample therefore keeps the real HOME and lets the session's own lease
  resolve the store (the configuration directory is still fenced).

## Baseline -- PRE-REGISTRATION (written before any baseline sample)

Instrument: `.agents/research/loop/matched_turns.py` (sha256 prefix
9eeeeb3cd3c20eb8, commit 5e7212fb), outcomes
`.agents/research/loop/matched_outcomes.py` (336519907d762797, commit
49eeaf8d), runner `scratchpad/behav/run_baseline.py` (sha256 prefix
7ae49bd5b9ee266d). All samples local on this Mac under one
`provider_lease` (one session at a time), interleaved in a fixed order
(for i in 0..11: D1 control-deepseek, D1 repaired-deepseek, D1
control-qwen, D1 repaired-qwen, D2 deepseek, D2 qwen), written to
`baseline/manifest.json` before the first sample. No outcome is read
until all 72 samples end or the stop rule fires; only infrastructure
fields (real turns, INFRA reasons, tokens) are watched while it runs.

Models (each in the configuration its own archive ran in, output cap
excepted): deepseek-v4-flash-0731, reasoning effort xhigh, max output
200,000, default transport deadlines -- the CUHK profile every R10 goal
and matched turn used (read in place, non-secret fields only); file
profile-deepseek.yaml sha256 be2192d59757ca97. qwen3.8-max, reasoning
effort max, enable_thinking true, deadlines 1800/1800/600/900 s -- the
provider record of the public po3-r19 / ino3-r12 sessions -- with max
output 65,536 (the largest tried that the provider accepts); file
profile-qwen-65k.yaml sha256 10b3883ca8835a07. Both files name the
profile `alibaba-token-plan`, context 1,000,000. Consequence stated: the
model arms also differ in reasoning-effort word and thinking switch, as
their archives did; the comparison is of each model as run, not of
weights at equal settings.

**D1 -- a message lever at a ceiling decision point (R10 Q26 CF).** The
p-benzyne task and Q26's scripted prefix (TASK.md 6eedde207fbc243a,
pbenzyne.xyz e7626c652c53d525, envelope 7bbbec82a10a055f with scratch
moved locally, prefix-transcript.json 729ee62dbf2d2f16; copied from CUHK
r10/q26/cf, digests equal), cut 7 (9 archived calls), 3 real turns,
plain planning session with the envelope, as Q26 ran. Trees: control
9297d6ba (package .py digest 406f5d5277d3dbe1), repaired 6b5adf48
(e62ab26648d894ea), exported from git. Stub runs on both: 9 of 9 prefix
calls reproduce the archived statuses, the triplet compile is
preview_failed in both, carrying "ValueError" and "inspect the
generated-input validation findings" on control and the program's
broken-symmetry sentence on repaired. Outcome: the first change after
the refusal, Q26's classes (matched_outcomes d1); primary binary ROUTE.
Prior evidence, not cells: deepseek on CUHK, ROUTE 6/6 and 6/6 (CF),
6/6 and 6/6 (CF-G).

**D2 -- unscripted noticing and hypothesis at a goal's first moment
(po3-r19, cycle 1, from its start).** The chemist's task text lists each
product file beside a SMILES whose ester locant contradicts the file's
name (triazole-ester-at-c4.xyz is the ester-at-C5 isomer), and asks for
"the isomer bearing the ester at C4"; nothing asks for the files to be
checked. The archived qwen3.8-max sessions noticed it (po3-r19 message
15, recorded in its decision's assumptions; po3-r17 by the case's own
record) on the 2026-09-11 tree, where the tools it used were in view and
the host differed; this is a new observation on this tree, not a
replication. Setup: the goal driver's first cycle on the R11 head
(package digest c5bc9af7708f415d), goal id goal-po3-r19, granted_by
claude-researcher-behav-owner-delegated (delegated; nothing is approved
or run), max_revisions 12, the archived envelope (scratch moved), the
four archived geometries (ids reproduce 9237de03..., 023cc0b5...,
e57f5c9d..., 560bf8cb...), task.md 3c257d20000d9d80; cut 0, 4 real
turns. Outcomes (matched_outcomes d2): LOOK, HYP, DECLARE (typed) and
NOTICE (text, secondary; mechanical rule primary, every sample then
hand-read unblinded and disagreements listed).

**N and power.** N = 12 per cell, 72 samples. Two-sided Fisher exact,
alpha 0.05, exact power by enumeration: 0.81 for 1.0 vs 0.5; 0.83 for
0.9 vs 0.3; 0.61 for 0.8 vs 0.3 and 0.7 vs 0.2; 0.42 for 0.9 vs 0.5;
0.28 for 0.3 vs 0.0. At N = 12 a 12/12 arm is separated from 7/12 or
fewer (p 0.037), and a 0/12 arm from 5/12 or more.

**Comparisons.** D1-M (per model): repaired vs control ROUTE. D1-X (per
tree): deepseek vs qwen ROUTE. D2-X: deepseek vs qwen on LOOK, HYP,
NOTICE (DECLARE reported, not tested). Every count is reported with its
denominator; no pooling across points.

**Predictions (priors, written before any sample).** D1 deepseek ROUTE
>= 10/12 in both arms; D1 qwen >= 8/12 control and >= 10/12 repaired.
D2 DECLARE >= 9/12 in both models; LOOK qwen about 6/12, deepseek about
3/12; NOTICE qwen about 4/12, deepseek about 1/12; HYP qwen about 3/12,
deepseek about 1/12. These are guesses; the rules below decide.

**What each result does to C4 (fixed now).**
- D1 null in both models (|repaired - control| <= 2 ROUTE each):
  supports "the message does not move the model" beyond one model,
  narrowed to a refusal whose cause the model can see (claims.yaml C4's
  own counterexample).
- D1 moves qwen (repaired - control >= 5, p <= 0.05) and not deepseek:
  contradicts C4's "not the messages ... beyond one model" for this
  message class; the message null is deepseek's.
- D1 control arms differ between models (p <= 0.05): the model is the
  locus at this point, whatever the message does.
- D2 models differ (p <= 0.05) on LOOK, NOTICE or HYP: unscripted
  noticing or hypothesising at this moment is model-dependent; C4 is
  narrowed to name the model as a factor.
- D2 both models at or below 3/12 on LOOK and NOTICE: the goal's first
  moment does not elicit premise checking in either model; the one arm
  then tests a general act or moment (never a task hint) on this point.
- D2 both at or above 9/12 on LOOK and NOTICE: ceiling; no arm is built
  on this point.

**INFRA and stop rule.** A sample with no real turn, a transport failure
that ends it, a turn deadline, a session error, or a prefix that does not
reproduce the archived statuses is INFRA: reported, never counted, never
re-rolled. The runner stops after three consecutive samples lost to
transport (unrecovered); recovered retries count nothing.

**Cost.** Estimated 9-12 M input tokens (D1 about 3 x 35-45 k per sample;
D2 about 4 turns rising from 9 k); cap for the baseline 20 M.

**D3 (conditional, not issued).** R10 Q32's T point (the partial-scan
route; deepseek 2/6 on the repaired tree, 0/6 on the base) is the one
archived point with a known mid-range rate, but it replays a goal cycle
whose records name `/lustre` paths and needs CUHK and 180-310 k tokens a
turn. Issued only if the baseline leaves budget and the master confirms
whether the doorway counts a key-using non-goal slot job as an Agent.

## Jobs issued

- 2026-09-28: baseline runner (run_baseline.py 7ae49bd5b9ee266d), local,
  one provider lease, 72 samples in the manifest order, pre-registration
  commit 1ace84a3. Pre-launch checks (provider-free): the qwen profile
  loads on the R10 control tree and its D1 prefix reproduces 9 of 9; both
  R10 trees resolve the fenced configuration directory with the real
  HOME.

## Status

- 2026-09-28: kernel, CONDUCT, RSL and lessons, charter topics
  (architecture, goal-grain-recovery-wake, dispatch-excursion), claims,
  REPORT-R10 sections 2 and 4, Q26 and Q32 records, the R10 harnesses read.
- Deliverables 1 and 2 committed (above); stub validation on three trees;
  leased probes: both models served.
- Next: issue the baseline runner under one lease; read nothing but
  infrastructure fields until it ends.
