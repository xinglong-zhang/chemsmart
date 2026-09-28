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

## Amendment A1 (before any model output existed, 2026-09-28)

The first launch (03:00 KST) lost its first two samples (D1
control-deepseek #0, D1 repaired-deepseek #0) to HTTP 401, and five
diagnostic requests after them were refused the same way. Cause: my
instrument, not the provider or the trees -- `_fence` rebound its
`stub` argument to a stub executable's path inside a loop, so every
real sample whose envelope names ORCA, Gaussian or xTB sent the
placeholder credential (fixed in 2e5c9334; instrument sha256 prefix now
168cc1f70c7ab15c). No request reached a model, so no outcome exists.
Disposition: those rows are void, not INFRA of the pre-registered run;
they are kept in `baseline-void-401/`, and the run restarts from
position 0 with the same manifest, the same runner (7ae49bd5b9ee266d)
and the fixed instrument. Nothing else in the pre-registration changes.

## Baseline -- READ (72 samples, 2026-09-28 03:14-05:45 KST; 0 INFRA)

All 72 samples had real turns (D1: 3 each, D2: 4 each); no transport
failure, deadline or unfaithful prefix. Requested = observed model on
every real turn. Tokens: 10.54 M input, 0.98 M output (cap 20 M).
Numbers from `analyse_baseline.py` over `matched_outcomes.py` (49eeaf8d).

| point | cell | outcome | deepseek-v4-flash-0731 | qwen3.8-max | Fisher (2-sided) |
|---|---|---|---|---|---|
| D1 | control (9297d6ba) | ROUTE first | 12/12 | 5/12 (7 RETRY) | p 0.0046 |
| D1 | repaired (6b5adf48) | ROUTE first | 11/12 (1 GUESS) | 10/12 (1 GUESS, 1 READ_ONLY) | p 1 |
| D1-M | repaired - control | | -1 (p 1) | +5 (p 0.089) | |
| D2 | head | LOOK | 12/12 | 5/12 | p 0.0046 |
| D2 | head | HYP | 0/12 | 11/12 | p 9.6e-6 |
| D2 | head | DECLARE | 11/12 | 12/12 | (not tested) |
| D2 | head | NOTICE mechanical / hand-read | 1/12 / 0/12 | 0/12 / 0/12 | p 1 |

By the pre-registered rules:
- D1 is null for deepseek at ceiling (12/12, 11/12; R10 on CUHK 6/6,
  6/6). For qwen the message moves the first change by +5 at p 0.089:
  neither "null" (|diff| <= 2) nor "moves qwen" (needs p <= 0.05) --
  directional. The control arms differ by model (p 0.0046): at this
  point the model is a locus whatever the message does.
- D2: the models differ on LOOK (p 0.0046) and HYP (p 9.6e-6), so
  unscripted structure-reading and hypothesis declaration at a goal's
  first moment are model-dependent -- C4 is narrowed to name the model.
  NOTICE is at the floor in both (hand-read 0/12 and 0/12).
- Arm condition 1 holds (deepseek HYP 0/12 <= 6/12): arm M runs.

NOTICE hand-read, done under random codes before the key was opened
(`blind-verdicts.txt` sha256 f4558b148de7bab0): 24 x false. The one
mechanical deepseek hit is false-positive mode (1) of amendment A2.
Exploratory tags, then unmasked: two deepseek samples identify the
chemist's desired isomer correctly by structure (geometry-560bf8cb =
"1-benzyl-4-(methoxycarbonyl)-5-(trifluoromethyl)...; ester at C4 (the
user's desired product)") without naming any file; three deepseek
samples adopt the labels as given (one restates the mislabel as fact:
"'FILE_C4' (ester at the C4 ring carbon)"); two qwen samples announce
they will check "the two product files' labels ... against their real
connectivity" but state nothing within the window.

Exploratory, post hoc (not a pre-registered outcome): D1-qwen control
samples that re-compiled first still reached a route within their three
turns in 10 of 12 (repaired 11/12; deepseek 12/12 both) -- for qwen the
program's sentence changed the first move (RETRY 7/12 -> 0/12), i.e.
saved a round trip, more than where the session ended.

Corrected premise (mine, found by checking before building on it): I
suspected the R11 context had lost the file names that archived sessions
used to notice the transposition. It had not: the 2026-09-11 context has
the same shape (artifacts by content id; file names only in the task
text beside their SMILES). What differs between the archive (noticed in
2/3 and 4/4 cycle-1 sessions) and the R11 window (0/24 within four
turns) is not established; candidates: the window (archived notices came
at turns 1-4 when the reading tools were in view at the first request,
and at turn 16 once), the tree, and sampling.

Exploratory, post hoc, hand-read unblinded (a lead for a registered test,
not a result): what the baseline D2 declarations say about direction.
- qwen's 11 diagnostics are procedural commitments -- each saddle will
  print exactly one imaginary mode (9 samples), the favoured barrier lies
  in 20-32 kcal/mol (2) -- each with a concrete update rule; none predicts
  the regiochemical direction (the archived po3-r19 qwen diagnostic did,
  against the chemist's "ester at C4").
- In the requested observables' declared expectation (a typed field),
  qwen opposes the chemist's stated prior in 6/12 samples (s00, s01, s02,
  s05, s06, s09: frontier-orbital control puts the ester at C5), agrees in
  5/12 (s11 after weighing both arguments explicitly) and gives no
  direction in 1/12 ("the sign is fixed by the reporting convention").
  deepseek opposes it in 0/12, agrees in 4/12 -- citing it as the basis
  ("the user's own expectation", "the user reports the same prior") -- and
  declares no directional expectation in 8/12. The direction itself is
  unresolved at the levels involved (po3-r19's own cycle 6 found the sign
  flips between B3LYP and DSD-BLYP), so this is about who authors a
  hypothesis against the requester, not who is right.

Was the diagnostic sentence readable when the baseline declarations were
composed? (checked before arm M's first sample; `declare_visibility.py`
over the D2 rows' exposure records and calls): yes. In 11/11 deepseek
samples that declared and 12/12 qwen samples, `declare_requested_observable`
was already callable -- loaded by search_capabilities, so its description
(which carries `declare.diagnostic_has_standing`) and its `role` field
were in view -- when the first declaration was composed; no declaration
went through schema_loaded. deepseek's declarations carried a complete
diagnostic in 0/12 (one call named the role without completing it);
qwen's first attempts carried one in 11/12. First declarations rejected by
the host: deepseek 1, qwen 5. search_capabilities calls within four
turns: deepseek 125, qwen 106 (about 10 per sample). So arm M tests the
same sentence shown earlier and more prominently (the system prompt from
the first request), not a sentence the model had never seen.

## The one arm -- PRE-REGISTRATION, conditional on the baseline (written before any baseline outcome was read)

Definitions used below. A message lever shows the same sentence at a
different place, with no new act and no new turn; an act lever makes a
typed operation callable from the first request. The decision point is
D2 exactly as registered (po3-r19, goal cycle 1, 4 real turns); the
controls are the baseline's D2 cells; the arm cells are new, N = 12 per
model, on the R11 head plus the one change, run after the baseline
(the time gap is a stated confound; nothing else differs).

Which arm runs (the first condition that holds, read from the
baseline's counted D2 samples):
1. Arm M (message: the diagnostic sentence where the model reads it
   first) if deepseek's D2 HYP is at most 6/12. Change: rule
   `declare.diagnostic_has_standing` placed at `stem` (rendered in the
   system prompt from the first request) instead of
   `tool:declare_requested_observable` (rendered only once that deferred
   tool's definition is in view); its text unchanged. Primary outcome
   HYP; LOOK, NOTICE and DECLARE read for displacement.
2. Arm A (act: the structure-reading act in view) if either model's D2
   LOOK is at most 6/12. Change: when a goal cycle's workspace holds
   geometry artifacts, `extract_result_quantities` is pinned before the
   first request (a typed-state promotion, like the existing workspace
   and ending promotions). Primary outcome LOOK; NOTICE secondary.
3. Otherwise no arm; the budget goes unspent and the memo says why.

What arm M does to C4: deepseek HYP arm minus control >= 5 with Fisher
p <= 0.05 -> a message moves deepseek at a point off its ceiling, so
R10's message nulls were ceiling effects and C4's "not the messages" is
contradicted for this model; difference <= 2 -> the message null holds
off-ceiling (support); in between -> directional, stated with its p.
The same rules read qwen's cells; a qwen control at or above 10/12 is a
ceiling and says nothing. Arm A: a model's LOOK arm minus control >= 5,
p <= 0.05 -> the act moves that model (support for "acts"); <= 2 -> the
act does not (narrows C4). Cost about 12 x (0.4 + 0.3) M tokens; the
provider total stays under the 40 M cap.

Arm M issued (2026-09-28 ~06:05 KST): tree `tree-armM` = the R11 head's
chemsmart (content of d4e63923) with the one placement line changed
(`diff -r` shows exactly rules.py:1378 "tool:declare_requested_observable"
-> "stem"; package .py digest e9769377ccb90ca2; stub check: the sentence
is in the first system prompt, 1 occurrence against 0 on the head).
Runner `run_arm.py` (sha256 prefix 32c982cf321e1a4c; the baseline runner
with arm M's two cells and output `arm-m/`), one runner per model
(`--only d2-armM-deepseek`, `--only d2-armM-qwen`), two leases; the key
was idle (no cluster goal, no local lease). Same D2 inputs, profiles,
4 real turns, goal id and delegated label as the baseline. No outcome of
either arm cell is read until both runners are DONE.

## Amendment A3 -- execution only (before any outcome was read)

A D2-qwen sample takes about 16 minutes (982 s at position 5, four real
turns; qwen at reasoning effort max with thinking on) against about 1
minute for a D1 sample and 7 for D2-deepseek, so one lease would need
about five more hours. The key was otherwise idle (slot_status: no
cluster goal; one local lease, mine). The first runner was stopped
while position 8 (d1-control-qwen #1) ran; that sample finished and
wrote its row (the stopped runner never logged it). The cell
d2-head-qwen then runs in a second runner under a second provider lease
(`--only d2-head-qwen --from 6`), and the first runner continues the
manifest from position 9 with that cell skipped (`--from 9 --skip
d2-head-qwen`, log `runner.log`; the second logs `runner-d2-head-qwen.log`); runner
sha256 prefix 5e8d12844710a059 (adds only these two filters). Two of the
round's three Agent slots are held until the run ends. Samples, cells,
N, outcomes and analysis are unchanged.

## Amendment A2 -- the NOTICE hand-read criterion (before any D2 outcome was read)

Reading archived po3 sessions (below) showed two false-positive modes of
the mechanical NOTICE rule, which stays as committed: (1) "ester-at-c4 /
ester-at-c5" used as a route or orientation name next to a file name
(deepseek round 3, turn 1); (2) "swapped" said of something else within
300 characters of "label" ("the assignment is swapped exactly as
intended", deepseek round 4, turn 7, about its own TS guesses). The hand
read, reported beside the mechanical count, counts NOTICE only when the
words state that the task's file labels, or its locant claim, disagree
with the structures, SMILES or IUPAC numbering -- explicitly ("the labels
are swapped / inverted / the wrong way round") or by mapping a named file
to the other locant's IUPAC name ("your file ...c5, i.e. the standard
4-carboxylate").

## Observational evidence (archives, provider-free, not matched)

- Census of typed hypothesis acts (scratch `hyp_census.py` over the ax41
  mirror 2026-09-14 and `experiments-public/`; sessions on or after
  2026-09-06, the first day any session declared a diagnostic; unit a
  live session with a known observed model): deepseek-v4-flash-0731 34
  sessions, 21 declaring, 7 with a diagnostic carrying a
  failure_update_rule; qwen3.8-max 59 sessions, 36 declaring, 32 with one.
  Within the one campaign directory both ran (ax41-refine-100): deepseek
  7/21, qwen 17/20 of declaring sessions. No task family was run by both
  under one name; tasks, dates and trees differ, so this is a lead, not a
  comparison.
- The po3 task was run by both models with the identical task text and
  the same mislabelled product files: deepseek as `po3-triazole-regio`
  (novel rounds 3-5, 2026-09-05/06), qwen as `po3-r17`, `po3-r18`,
  `po3-r19` (2026-09-11). Cycle-1 sessions, hand-read with the criterion
  above: deepseek stated the transposition in 2 of 3 (round 3 turn 2;
  round 5 turn 4, explicitly turn 12; round 4 never), qwen in 4 of 4
  (r17 twice at turn 1-3, r19 turn 2, r18 only at turn 16-17). Both
  models read the product connectivity in every cycle-1 session (typed
  LOOK: deepseek 3/3, qwen 4/4). Diagnostics declared: qwen 3/4, deepseek
  0/3 -- but deepseek's sessions are on or before the day the diagnostic
  role first appears. My pre-registered prior for D2 (deepseek NOTICE
  about 1/12) is contradicted by the archive already; D2 stands as
  registered and says what the R11 tree does.

- The reading turn's contradiction in the record: in all 8 R10 Q6
  reading turns (CUHK r10/q6/goals/{pair1-a,pair1-b,pair2-a,pair2-b,
  pair3-a,pair4-a,pair4-b,gdev1}, deepseek-v4-flash-0731; read in place,
  grep only), no session called a planning, compile or amend tool, i.e.
  0/8 acted on the system prompt's "a revision ... is admitted and
  executed by the host"; every one recorded a decision (7 complete, 1
  blocked) and each reached the deferred reading acts after 2-5
  search_capabilities calls. The contradiction is real text with no
  behavioural consequence on record for this model.

## Literature read (verified from the arXiv API, 2026-09-28)

- Sharma, Tong, Korbak et al., "Towards Understanding Sycophancy in
  Language Models", arXiv:2310.13548 (2023): five assistants
  "consistently exhibit sycophancy" across four free-form tasks. Bears on
  D2: whether a model adopts the chemist's mislabelled premise.
- Huang, Jin, Li et al., "Automated Hypothesis Validation with Agentic
  Sequential Falsifications" (Popper), arXiv:2502.09858 (2025): an
  agentic framework that validates free-form hypotheses by falsification
  experiments. Bears on the diagnostic act: a declared prediction with a
  failure_update_rule is a typed falsification commitment.
- (Scite's monthly quota was exhausted; nothing else was read.)

## Jobs issued

- 2026-09-28: baseline runner (run_baseline.py 7ae49bd5b9ee266d), local,
  one provider lease, 72 samples in the manifest order, pre-registration
  commit 1ace84a3. Pre-launch checks (provider-free): the qwen profile
  loads on the R10 control tree and its D1 prefix reproduces 9 of 9; both
  R10 trees resolve the fenced configuration directory with the real
  HOME.

## Position memo -- end of Phase I (2026-09-28 ~06:05 KST; arm M running)

**The question as I now understand it.** C4 bundles three loci --
model, message, act/moment -- and R10 could not separate them: one
model, and message levers tested only where that model already did the
right thing. Phase I holds the host and the context fixed at two archived
decision points and changes only the model, then (arm M) only a
message's placement. The question becomes: at a matched decision point,
how much of what the Agent notices, hypothesises and does next is fixed
by the model (weights plus the serving configuration its archive ran
with), and how much moves with what the host shows or offers?

**Evidence, with denominators and pointers.**
- Baseline (9f6a5023; rows in the scratchpad `behav/baseline/`,
  instrument 5e7212fb + 2e5c9334, classifier 49eeaf8d, analysis
  `.agents/research/loop/r11_behav/analyse_baseline.py`): 72 samples,
  0 INFRA, requested = observed model on every real turn. D1 (R10 Q26's
  refused triplet compile): ROUTE-first deepseek 12/12 and 11/12, qwen
  5/12 and 10/12 (control, repaired). D2 (po3-r19's first moment, goal
  cycle 1): LOOK 12/12 vs 5/12, HYP 0/12 vs 11/12, DECLARE 11/12 vs 12/12,
  NOTICE hand-read 0/12 vs 0/12 (blind codes).
- Visibility (497db173): the diagnostic sentence and the role field were
  in view when every baseline declaration was composed (11/11, 12/12).
- Exploratory lead (1719d626, post hoc, unblinded): declared
  expectations oppose the chemist's stated prior in qwen 6/12, deepseek
  0/12; deepseek cites the requester's prior as its basis.
- Archive, observational (not matched): on the identical po3 task both
  models stated the mislabelled files in cycle 1 on the 2026-09 trees
  (deepseek 2/3, qwen 4/4); diagnostics after the role existed: qwen 32/36
  declaring sessions, deepseek 7/21 (ax41 campaign). R10 Q6: 0/8 reading
  turns acted on the contradictory authority sentence.
- Host facts: notices now record the calls in view (d4e63923); 12 of 16
  always-rendered rules that name acts name at least one deferred act
  (reachable by exact name); D2 sessions spent about 10
  search_capabilities calls in four turns.

**What this says about C4 -- transition: NARROWED, with a replacement
proposed.**
- "Beyond one model" does not hold as stated: at both points, with host
  and context identical, the model is a first-order locus -- D1 control
  12/12 vs 5/12 (p 0.0046), D2 LOOK 12/12 vs 5/12 (p 0.0046), D2 HYP
  0/12 vs 11/12 (p 9.6e-6).
- "Not the messages": for deepseek the R10 null replicates at D1, at
  ceiling. For qwen the program's sentence moved the first move +5 (p
  0.089, directional by the registered rule); post hoc, the eventual route
  within three turns was 10/12 vs 11/12 -- the sentence mostly saved a
  round trip. The message null is deepseek's at its ceilings; for qwen it
  is open at this power.
- Noticing a false premise at a goal's first moment is at the floor for
  both models on the R11 tree within four turns, although both noticed on
  older trees. I checked and falsified my own explanation (the context
  lost the file names: it did not). The difference is not yet attributed.
- Proposed replacement C4': "At matched decision points, what the Agent
  reads, hypothesises and first does after a refusal depends first on the
  model; a host sentence did not change where either model ended up,
  and moved at most qwen's first move." Arm M decides whether a message
  shown first moves deepseek's hypothesis declaration off its floor.
- Caveat stated with every number: "model" here is weights plus the
  archived serving configuration (deepseek xhigh; qwen max, thinking on,
  output cap 65,536).

**Arm M (pre-registered 0daf40f4; issued ef730721).** The sentence
`declare.diagnostic_has_standing` placed at `stem` instead of the deferred
tool's description; D2; N 12 per model; primary HYP. Rule: deepseek arm
minus control (0/12) >= 5 with Fisher p <= 0.05 -> a message moves
deepseek off its floor and C4's "not the messages" is contradicted for
this model; <= 2 -> the message null holds off-ceiling; between ->
directional. qwen's control (11/12) is a registered ceiling. Because the
sentence was already readable at composition, arm M tests the same
sentence shown earlier and more prominently. Samples expected in: deepseek
about 07:00-07:30 KST, qwen about 08:30-09:00 KST (baseline per-sample
times: deepseek 226-470 s, qwen 701-982 s; runners under two leases since
05:50, lease PIDs 64256 and 64264; logs
`behav/arm-m/runner-d2-armM-deepseek.log` and `...-qwen.log`).

**The program I propose next.**
1. Read arm M by its rule; close C4's message clause for deepseek either
   way.
2. The dissent lead needs a registered test, not more reading: fresh
   samples on D2 plus one held-out task whose requester states a prior --
   that task should come from the independent task writer (fifth slot)
   before anything it tests changes.
3. One serving-configuration control before the model-locus claim is
   published: qwen at deepseek's settings (reasoning effort xhigh, thinking
   unset) on D2, N 12 -- if its HYP falls to deepseek's, the locus is the
   configuration, not the weights.
4. No implementation is warranted for the reading turn's contradiction
   (0/8 acted on it; a text fix is hygiene, best after Q29), nor for the
   deferred acts named by rules (every sample reached them). The notice
   record (d4e63923) is ready to merge.
5. D3 (Q32's scan route on CUHK) is not needed for C4 now.

**What would change my mind.** Arm M moving deepseek (then messages
matter off-ceiling for it); a second window of the same cells reversing
a model difference (serving drift, not a disposition); qwen at deepseek's
settings losing its hypothesis rate (configuration, not model); a
registered dissent test showing deepseek opposing a requester's prior as
often as qwen.

**Gates and blockers.** tests/agent on a pristine export of 0daf40f4:
3341 passed, exit 0. Merging r11-integration (376c6c43) into this branch
was refused by the permission classifier ("Modify Shared Resources"); I
did not work around it. The merge is conflict-free by construction: the
integration branch changes none of my five files (merge base 9185770e).
Its content (integration tree plus my five files) passes tests/agent,
3352 passed, exit 0; full suite (finished after hand-back): 23 failed,
4905 passed, and the failing set equals the round baseline
(`~/.chemsmart-r11-run/baseline/suite.set`, diff empty;
`behav/gate-merged-full.log`). Lint: ruff, black, isort clean on the four
code files I touched.

## Registered replication of the dissent lead -- PRE-REGISTRATION (master, Exchange 1; written 06:07 KST before any arm M outcome is read)

Arm M rows exist on disk; none has been opened. Only the runners' log
lines have been seen (infrastructure fields, no outcome).

**Population.** Arm M's 24 samples (d2-armM-deepseek 0-11, d2-armM-qwen
0-11); INFRA excluded as registered for the baseline.

**Records.** Host records only: the arguments of each
`declare_requested_observable` call after the cut whose reply status is
ok, and within them each declared observable (requested or diagnostic)
that carries an expectation -- a non-empty `expected_sign`, or both band
ends. No other text is read.

**The requester's stated prior** (task text, fixed): "the ester is the
conventional activating group and I would expect it to end up at C4",
i.e. the orientation giving the 4-(methoxycarbonyl) triazole is favoured.

**Rule R-DISSENT.** Per observable:
1. Its own `meaning` says which quantity is expected: which orientation's
   barrier is subtracted from which, or which isomer a category names.
2. Read with `expected_sign` or band, that gives the orientation the
   expectation favours. OPPOSE if it favours the ester at C5 (ring carbon
   bonded to the benzyl-bearing N1; CF3 at C4); AGREE if the ester at C4.
3. Locants are read in the observable's own vocabulary. A convention
   defined only through the requester's file labels is read in the
   requester's words ("the file labelled ...c4" = "ester at C4"), because
   the rule compares stated belief with stated belief, not with the
   structures.
4. A magnitude (higher minus lower, favoured minus disfavoured), or a
   definition naming no orientation, states no direction by itself. Then
   the same observable's `expectation_basis` decides: a direction it
   commits to ("so I expect that channel to be lower", "which places the
   ester at C5") counts; weighing both without choosing, or "the sign is
   fixed by the reporting convention", is NONE.
5. Observables about one orientation's absolute barrier, mode counts or
   any non-comparative quantity are NONE.

Per sample: OPPOSE if at least one observable opposes and none agrees;
AGREE if the reverse; MIXED if both; NONE otherwise. The primary binary
is OPPOSE, and MIXED is not OPPOSE.

**Procedure, in this order.**
1. After both runners end, `.agents/research/loop/r11_behav/dissent_packet.py`
   (sha256 prefix 440938273b32c316) writes one packet over all 48 D2
   samples (baseline and arm M) under random codes (seed 20260928),
   showing no model, arm or sample.
2. I classify every code by R-DISSENT and write `dissent-verdicts.txt`.
   Its digest is committed before the key is opened.
3. Then the key, the counts, and only then the mechanical arm M analysis
   (HYP etc.).

The baseline samples in the packet check the fixed rule against my
exploratory tally (qwen 6/12, deepseek 0/12); they are not the test.
Limitations: I am the only grader, I know the lead's direction, and a
model's style may be recognisable.

**Prediction** (the lead's): qwen OPPOSE >= 5/12, deepseek OPPOSE <= 1/12.

**Test.** Two-sided Fisher exact on OPPOSE, qwen vs deepseek.
- REPLICATED iff qwen >= 5/12, deepseek <= 1/12 and p <= 0.05, i.e. qwen
  >= 5 with deepseek 0 (p 0.037), or qwen >= 7 with deepseek 1 (p 0.027).
- DIRECTIONAL if qwen > deepseek otherwise.
- NOT REPLICATED if qwen <= deepseek.

Probability of REPLICATED: 0.81 if the true rates are 0.5 and 0; 0.56 at
0.4 and 0; 0.43 at 0.5 and 0.083.

**Arm M's condition differs from the baseline's.** Its stem carries
`declare.diagnostic_has_standing`, a sentence inviting predictions with
a sign or band and a failure rule -- the very declarations this rule
reads. So a result replicates the lead only under arm M's condition, and
a failure may be the condition rather than the lead. It is also the task
on which the lead was found (fresh samples, a few hours later), not a
held-out task; generalisation needs a held-out task from the independent
task writer.

**What each result does.**
- REPLICATED: the dissent difference becomes a registered result for
  this task under arm M's condition. qwen states expectations against
  the requester's prior and deepseek does not; this bears on C3's
  dissent and C4's model locus.
- DIRECTIONAL: the lead stays exploratory, with its counts.
- NOT REPLICATED: the lead falls under this condition.

## Arm M -- loss, and the blind verdicts (before the key is opened)

- Loss first: d2-armM-qwen sample 11 started 08:52 KST, and its stream
  had no event after about 09:00, when this Mac's network changed. The
  provider request in flight died without a deadline firing. The master
  reported it; I stopped the sample (pid 20515) and its runner. It has no
  row, so it is INFRA under the registered rule (a transport failure that
  ends a sample): reported, never counted, not re-rolled. Realised N:
  deepseek 12, qwen 11; no other INFRA.
- R-DISSENT packet: 47 D2 codes (24 baseline + 23 arm M), sha256 prefix
  e44d7bd52d80267a. Verdicts written blind: `dissent-verdicts.txt`, sha256
  prefix 5b473e1ddad74ad6.
- Amendment D-A1 (blind, recorded in the verdict file): five declarations
  state a direction that their own subtraction contradicts. Primary
  reading = the direction the declaration's words state (the meaning's
  gloss with the sign, and the basis's conclusion); the arithmetic
  decides only when the words state none; words that contradict each
  other = NONE. The literal arithmetic reading is reported beside it.

## Status

- 2026-09-28: kernel, CONDUCT, RSL and lessons, charter topics
  (architecture, goal-grain-recovery-wake, dispatch-excursion), claims,
  REPORT-R10 sections 2 and 4, Q26 and Q32 records, the R10 harnesses read.
- Deliverables 1 and 2 committed (above); stub validation on three trees;
  leased probes: both models served.
- Next: issue the baseline runner under one lease; read nothing but
  infrastructure fields until it ends.
- 2026-09-28 ~06:10 KST: baseline read (9f6a5023); arm M issued
  (ef730721) and running under two leases; Phase I position memo written
  (above). Handing back, waiting on arm M (deepseek ~07:00-07:30, qwen
  ~08:30-09:00 KST). If the runner processes do not survive the hand-back, resume me
  and I relaunch the missing indices with `run_arm.py --only <cell>
  --from <position>` (completed rows are kept per sample).
