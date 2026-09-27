# R11 episode truth -- does the host sign only words its records support, and does any refusal bury a legitimate scientific choice?

Base SHA: 9185770e30e486e0d1e96e83bf35edc914b8590a (verified with
`git rev-parse HEAD` as the first action, 2026-09-28). Model:
claude-opus-5-5[1m]. Episode id `truth`, lens `truth`. Brief sha256
prefixes: brief-truth.md 21724e90de649a55, common-truth.md
ff6d4647e499b66d (both equal `~/.chemsmart-r11-run/briefs/SHA256-prefixes.txt`).

## The question (as currently understood)

Across every archived goal, does the CHEMSMART host at the round pin sign
only words that its records support -- settlements, verified refusals,
certifications, stationarity verdicts, finding standing, category answers
-- including in the goals where the Agent's conclusion diverged from its
task? And does any refusal or rule suppress a legitimate scientific choice
rather than protect an invariant?

Claims served (`~/.chemsmart-r11-run/CLAIMS.md`): C1 (each enumerated class
of signed word is true of the records it cites, recomputed at signing) and
C3 (evidence-bound dissent and unasked findings have standing; host control
does not suppress legitimate science).

## Corrected premises (before any census)

- P1 (brief): "only about 6 of the 210 commits since `002f91cf` touch the
  signing files". True of the five settlement files (3 commits touch
  `driver.py`/`goal.py`/`analysis_*`/`terminal_states.py`: 9964968e,
  451f5435, 1b388d75). False of the signing *surface*: the functions that
  decide stationarity, the free-energy surface, a refusal's re-read and a
  characterisation moved in R10 Q21/Q27/Q31/Q33 (ebf36450, 3041b7fb,
  6ca94a2f, cfaffd58, b25ea98b, 2bcdc763, 07cc7978 and more under
  `chemsmart/analysis/`). Words in classes W7, W10, W11 below are expected
  to move at the pin for rule-change reasons.
- P2 (mine, from reading `_achieved_word`, driver.py:491): the first word
  becomes `achieved_with_observations` only for host anomalies, falsified
  expectations and answered failed criteria in the delivery it reads. A
  session's *unrequested finding* rides in the reasons under plain
  `achieved`. That matches the charter's definition of the word; whether
  it gives an unasked finding "standing" is a C3 question, not a C1 defect.

## The classes of word the host signs (from the code, one signer each)

| id | word | signer (the one function) | record | how it is checked at the pin |
|---|---|---|---|---|
| W1 | settlement: achieved, achieved_with_observations, unreachable_from_evidence, exhausted, returned_to_human, with reasons | `driver.GoalDriver._settle` (run path); `driver._delivery_settlement` (planning / analysis-only path); `driver.GoalDriver._settle_after_reading`; `driver._typed_error_settlement` | ledger `goal_settled` | RE-SIGNED by the walk-mode settle replay (final step) + reader |
| W2 | held words: recovery_opened, reading_opened (what is owed, what is uncertified) | same signers | ledger rows | RE-SIGNED where it is the final settle-type row; else not |
| W3 | `qualified` capability rows | `driver._record_goal_qualification` | ledger `qualified` | RE-SIGNED by the replay + reader |
| W4 | executor analysis word: completed / partial / "" | `executor.ApprovedWorkflowExecutor._run_analysis_phase` | `run_recorded.analysis_status`, `execution-result.json`, `analysis_only_plan_executed` | RE-SIGNED by the walk for chains that ran no node; reader otherwise |
| W5 | certification: completion status passed/partial with limitation, anomaly, falsified-expectation and failed-criterion ids | `tool_runtime.CommandCompiledToolHostV1._record_toolchain_completion` | `analysis_completion_evaluated` | READER only (signed in-session with host state no replay rebuilds) |
| W6 | expectation verdict: agreed / diverged / indeterminate / not_comparable / agreed_as_approximation | `tool_runtime.CommandCompiledToolHostV1._declared_observable_predictions` | completion `declared_observable_predictions` | READER (arithmetic recomputed from the row) |
| W7 | verified refusal: verified + basis | `tool_runtime.CommandCompiledToolHostV1._verify_unreachable_observables` (+ `refusal_read_against_results`; re-read at settlement by `GoalDriver._refusals_the_results_now_answer`) | decision `unreachable_observables` | READER (basis present in the stream's plans; latest typed word; no later extraction of the named selector); the final cycle's re-read is inside the W1 replay |
| W8 | finding standing (answers / on_the_request / unrequested) and relation truth | `tool_runtime.CommandCompiledToolHostV1._verify_findings` (+ `analysis_claims.evaluate_finding_relation`) | decision `findings` | READER (standing from declared ids; `holds` from the recorded values) |
| W9 | category answer: the word the host read | `tool_runtime.CommandCompiledToolHostV1._categorical_answer_row` | finding `answer[]`, completion `declared_categorical_answers` | READER (word against the extraction receipt it names) |
| W10 | stationary-point order and stationarity | `execution.build_stationary_point_characterisation` (+ `result_quantities.structure_stationarity`) | `stationary_point_characterised` | RE-SIGNED as a pure function over the archived result file, where the file exists |
| W11 | a free energy stands on a stationary structure of a named surface | `result_quantities.structure_stationarity`, `result_quantities.free_energy_surface` | `thermochemistry_derived` | RE-SIGNED (stationarity of the archived source result) where the file exists; reader (source node's terminal state) |
| W12 | node terminal state and result validity | `terminal_states.derive_run_outcome`, program verifiers | `program_result_verified`, node state rows | READER-lite (stationary-point rule from the verification record); not re-signed in Phase I |
| W13 | anomaly standing: unreplicated / replicated / refuted | `execution.anomaly_standing` | `anomaly_observed`, ledger `anomalies_observed` | counted as observations the first word must carry (W1 reader) |
| W14 | sufficiency: met / attested / short / unstated | `delivery.judge_sufficiency` | claim `requirement_assessment` rows | READER if the rows carry the operands; else stated unchecked |
| W15 | a failed criterion answered | `goal.failed_criteria` | derived | folded into W1 and W5 |
| W16 | revision admitted / returned | `goal.admit_revision` | ledger | not in scope (authority, not truth of evidence) |
| W17 | refusal: gate, invariant, route | each gate (`RoutedContractError`, `ContractError`) | `tool_failed` | REFUSAL CENSUS (C3) |

C1 is narrowed up front: "recomputed at signing" is testable by replay only
for W1-W4, W10, W11. For W5-W9 and W12-W14 the words were signed inside a
session with host state no replay rebuilds, so the claim tested there is
"true of the records it cites", by an independent reader.

## PRE-REGISTRATION (written before any archived goal is read or re-signed)

Population. Every goal ledger in `~/.chemsmart-r11-run/ARCHIVES.md` except
sealed and private studies: ax41 mirror 2026-09-14 (campaign + research,
195), the public cases in `experiments-public/` (3), CUHK R8 (17), R9 (14)
and the R10 episode directories named in `instruments/r10_jobs.tsv` (never
`r10/m*`, `*/sealed*`, q29). Deduplicated by (goal_id, goal_sha256); every
directory read is named in the result. CUHK records are read and replayed
in place inside slot jobs; nothing is bulk-copied to the Mac.

Oracle. (i) The independent reader (`word_reader.py`, records only,
imports nothing from `chemsmart`); (ii) for re-signable classes, the same
replay on the commit that produced the record (CUHK: the campaign's own
`code/` directory; ax41: the commit identified from campaign records, or
"unidentified").

Divergence buckets for a word re-signed at the pin, decided mechanically:
- (a) the state differs;
- (b) the state is the same and the reasons' content differs: the sets of
  declared observable ids, 8-hex receipt prefixes, node ids and numbers
  named in the reasons are not equal;
- (c) state and content equal, wording differs;
- identical.
Attribution of (a) and (b): the producing-commit replay must reproduce the
archived word byte for byte (else "record insufficient", never counted
either way); then the change is a rule change, named by the commit that
moved it, and classified by the reader: the archived word was false and
the pin repaired it / the pin's word is false / both defensible. (c) is
reported with a count only.

Control for harness faithfulness: the 10 settled goals with the lowest
sha256 of their `goal_sha256` hex string (5 CUHK, 5 ax41) are replayed on
their producing commits whatever their pin bucket.

Falsifiers.
- F1 (C1 contradicted): the pin signs a word the reader flags and the
  flag survives reading the records it cites, in a class the host computes.
- F2 (harness): a producing-commit replay does not reproduce the archived
  word; that goal's comparison is not evidence.
- F3 (reader): a flagged word is defensible from its records; the reader
  is corrected and the correction reported as loudly as a finding.
- F4 (C3 contradicted): a refusal the Agent met refused a chemically
  legitimate choice with no walkable legal route and it changed a
  delivered answer; or a settlement word buries an evidence-bound dissent
  (verified refusal, answered criterion, falsified expectation or
  diagnostic, unasked finding, self-correction) its records hold.

What counts. C1 supported: every reader flag on a pin word lies in a named,
already-repaired class. Narrowed: residual classes named with counts, or
classes only reader-checkable. Contradicted: F1. Widened: the classes cover
every word a delivery rests on and a delivery's conclusion replays from its
receipts alone. C3 supported: dissent-bearing words carry each marker and
no refusal suppressed a legitimate choice that changed a delivered answer.
Narrowed: suppressions or burials named with counts. Contradicted: F4.

Refusal census scheme (per group of (tool, error class, gate or message
template), never per instance): P protects a CONDUCT section 2 invariant;
R resolves a reference the model got wrong (unknown id, schema); L refuses
a chemically legitimate choice (counted as a suppression only if it changed
a delivered answer); U names a route the host would refuse in that state.

## Census 1 on CUHK -- what it runs (written before submission)

Job `/project/xlzhang/jiseung/r11/truth/census1/job.sh`, one node, 4 cores,
16 GB, 3 h, provider-free. Code: this worktree at b89b6108, whose
`chemsmart/` tree equals the pin's (`git diff 9185770e HEAD -- chemsmart
pyproject.toml` is empty); tree digest 1371776d5dd0704c. Population: the
goal ledgers under r8/{gaussian,integration,orca,pyscf} (12),
r9/{gaussian,master,orca,pyscf,xtb} (14) and r10/{q1..q16,q18..q24,q26,
q27,q28,q30,q32} (73 ledgers, q32's ten replay copies of one goal among
them), deduplicated by goal digest and written to `census1/specs.txt`
before anything is replayed. r10/q6 is included (R10's report discusses
its goals openly; pair4-a is the unasked-finding case); r10/q17 is not
(sealed-study guard). The job: (1) `signed_words.py`; (2) `resign.py` at
the pin; (3) `word_reader.py` over the archived words and over the pin's;
(4) `resign_producing.py`: every goal whose job output names a code tree
that still hashes to its printed digest is replayed on it.

Expected, from the pre-registration: the producing-code replays reproduce
the archived words (else F2 for that goal); bucket (a) and (b) words at
the pin are rule changes named by commit; the reader's flags on the pin's
words lie only in classes Q24 or later rounds repaired (else F1).

## Local census (ax41 mirror + public; provider-free, on this Mac)

190 goals after deduplication (195 ax41 ledgers + 3 public; the in-repo
po3-r19 and ino3-r12 are the mirror's research copies).
- Re-signed at the pin (`resign.py`, walk mode): 9 unsettled, 4 typed
  errors (not replayed), 1 error (a copy of ino1 with no `goal.json`:
  record insufficient), 1 planning goal with no stream; of 175 compared,
  identical 9, wording only 30, content 86, state 50 (corrected: the
  first entry of this record said 10 and 49; the transitions below sum
  to 50).
- State changes (a): achieved -> recovery_opened 35, exhausted ->
  recovery_opened 5, returned_to_human -> recovery_opened 3,
  achieved_with_observations -> recovery_opened 2, achieved ->
  returned_to_human 2, achieved -> achieved_with_observations 2,
  unreachable_from_evidence -> returned_to_human 1. Q24 reported 34
  achieved -> recovery_opened at ec41a57c (historical classes repaired
  before R10); attribution needs producing-commit replays (pending).
- Reader over the archived words: 13 achieved over undelivered declared
  ids (equals Q24's count), 1 over an uncertified delivery, 2 plain
  achieved over what the run found, 13 qualified rows under a word that
  is not true, 4 expectation verdicts that disagree with their arithmetic
  (ino2, not_comparable where the numbers compare), 11 free energies
  from a result its verification did not pass (po3-r19 ts-esterc4 x8,
  po3-triazole-regio x3).
- Reader over the pin's words: 0 achieved over undelivered, uncertified
  or later-refused ids, 0 qualified rows under a false word; 1 plain
  achieved flagged (goal-h4). F3, my reader: goal-h4's five "diverged"
  rows are a declaration with `expected_sign: positive` and band [0, 0]
  for an imaginary-mode count whose value is 0 -- a zero has no sign, and
  the pin now refuses such a declaration where it is written. The
  expectation (a minimum) was met; the word is defensible.

## Census 1 on CUHK -- first read (job 2157057, pin replay and readers done)

89 goals after deduplication (81 settled). Reader over the archived
words: achieved over an uncertified delivery 7 (Q24's six empty-chain
words plus r10/q6 gdev1, which Q24 did not read), over undelivered ids 5,
over a later verified refusal 1 (Q24's g2r), plain achieved hiding what
the run found 2, qualified rows under a false word 10. Reader over the
pin's words: 0 uncertified, 0 later refusal, 0 qualified rows under a
false word; two flags survive and were read against their records:

- r10/q7 g2-scan-modred, plain `achieved` at the pin (bucket identical:
  the archived word reproduces). Cycle 1's run completion passed and
  listed `falsified_expectation:cis-barrier` (8.30 kcal/mol against the
  pre-registered band 2-8, role requested, not declared after evidence)
  and `falsified_expectation:oo160-torsion` (141.5 deg against 100-135).
  The goal settled in cycle 2 on a session completion that re-rendered
  nothing (four rows not_comparable, anomaly list empty) and names those
  numbers only as "delivered in an earlier cycle". The first word hides
  two pre-registered expectations the physics left, for numbers it
  certifies as delivered. Candidate F1 (see the class below).
- r10/q1 gdev1, `achieved_with_observations` at the pin (wording only):
  the declared category `rks-stability-verdict` is answered by a finding
  whose own reason line, signed by the pin, says it "rests on no word the
  host read, so it answers nothing", under "the host completion gate
  certified the delivery". The certificate is the archived completion,
  minted before a category had to be answered by a word the host read;
  at the pin such a finding is refused where it is written
  (`finding.answers_through_a_word_the_host_read`), so production cannot
  make this record. Classed: word changed by a rule the old certificate
  predates (the settlement trusts the completion's delivered list while
  the finding reason recomputes; two organs, one question).

The same mechanism as q7, found locally by the refined reader: ax41
goal-ino3-r17 settles `unreachable_from_evidence` at the pin; cycle 3's
passed session completion recorded three pre-registered expectations as
falsified (E vs Fc -0.494 V against [-0.3, 0.9]; quartet-doublet gap
103 kJ/mol against [5, 90]; spin population 0.002 against [0.05,
0.45]); the final word, signed from the last session's completion, names
none of them. CLASS (candidate): the settlement reads falsified
expectations from the completion its delivery stands on, while it
certifies ids an earlier cycle delivered; the observation an earlier
cycle's certificate recorded for those ids never reaches the goal's word
(the goal-grain rule already covers claims, anomalies and failed
criteria). Owning function: `driver._achieved_word` /
`_AnalysisDelivery` (inputs to the observed set).

Stationarity re-signed at the pin (local, `resign_stationarity.py`): of
370 archived free energies, 195 re-read (190 still on a stationary point,
43 of them delivered through a claim; 5 now refused -- po3-r19
ts-esterc4 x4, po3-triazole-regio ts-c4 -- all on ORCA saddle searches
that printed their own non-convergence, none delivered through a claim
the census can link; the charter records po3-r19's 23.194 kcal/mol as
claimed through an expression that predates recorded input bindings),
170 unread (105 digests name no verified result, 67 files absent). Of 16
characterisations, 11 certified again, 3 refused (po3-r19 ts-esterc4,
po3-triazole-regio ts-c5b: not converged; g5-phosphine planar Hessian
at 0.0138 Eh/Bohr), 2 unread. These are rule changes (R10 Q21/Q33).

Refusals met by the Agent (local, `refusal_census.py`): 826 in 190
groups. Read so far: `gibbs_free_energy` on ORCA opt/ts (150 of the
selector refusals) names its route (a thermochemistry stage with bound
temperature, pressure and standard state): P. The decision refused for a
functional-convention word (VWN3/VWN5/PZ81/PW92) without a resolution
receipt (18 instances, 16 goals): every one of the 18 sessions recorded a
decision later in the same session (0 buried); the message at the pin
still names neither invariant nor route (a wording defect, Behaviour's).
Named-site spin flips: ax41 ino2 returned to the human with no J under
its declared ids (no delivered answer changed); R10 Q15's flips were on
p-benzyne, where R10 Q18's oracles showed GuessMix and FlipSpin reach one
solution to the printed digit. Neither changed a delivered answer.

## Census 2 on CUHK -- what it runs (written before submission)

Job `/project/xlzhang/jiseung/r11/truth/census2/job.sh`, same code tree
and population as census 1, instruments at 2e7f0e23 copied to `tools2/`
(census 1's `tools/` is left untouched while it runs): the refined
reader over the archived words and over census 1's pin words; the
refusal census; the stationarity re-sign at the pin. Expected: the new
W1 check finds q7 g2-scan-modred on the pin's words and possibly further
instances of the same class; refusal groups beyond the local ones are
read and classed P/R/L/U by the scheme above.

## Census 1 -- producing code and eras (read)

Producing-code replays (CUHK, `resign_producing.py`): of 89 goals, the
code a goal's job printed still hashes to its printed digest for 66
(6 code moved, 17 no job output). Of those 66: 47 reproduce the archived
word byte for byte, 2 do not (r10/q6 pair3-a, pair4-a), 12 raise in the
harness on the old tree (older driver lacks a function the harness
calls: record insufficient), 5 unsettled or typed errors. F2 for the two:
both goals ran a reading turn; the harness compares the word held at
`reading_opened` but replays over streams that already hold the
reading's finding, so reading-turn goals are not faithfully replayed by
this harness (a harness limitation, not a host fact). pair4-a's final
archived word is plain `achieved` whose reasons carry "the session's
finding n-ch2cl-anomeric-shortening (not asked for)" with the relations
the host checked (N1-C 1.386 A vs N1-CH3 1.445 A).

CUHK pin changes, attributed: 12 state changes, of which 6 reproduce on
their producing code (g1-hono, g2-hooh, merged-smoke: achieved ->
recovery_opened, Q24's empty-chain repair; g2r: awo -> unreachable,
Q24 41f607d4; r10/q3 g2 and q9 g1: returned_to_human -> an achieved
word, Q16 H2's known repairs), 1 errors on its old tree (r9/orca g5,
whose producing replay raises the same contract error the archive
settled on), 1 code moved (r10/q6 gdev1: achieved over a partial newest
completion, false by the reader, repaired at the pin), 4 have no job
output (r8 goal-irc2, r9 gaussian g1 and g3, r9 infra-smoke: Q24
reproduced three of them on their commits).

Eras (local, 180 goals replayed on ec41a57c = R10 base, 002f91cf = Q24
merged, and the pin): 49 state changes happened between the archive and
the R10 base (34 achieved -> recovery_opened, Q24's historical classes;
5 exhausted -> recovery_opened; 3 returned -> recovery_opened; 2 awo ->
recovery_opened; 2 achieved -> returned; 2 achieved -> awo; 1 unreachable
-> returned), 1 at Q24's merge (e2-acetone, the empty chain), and none
after it: the pin and 002f91cf sign byte-identical words (state and
reasons) for all 175 local goals both replayed. The brief's prior "the
pin may sign what Q24's tree signed" holds for this population.

## Census 3 on CUHK -- what it runs (written before submission)

Job `/project/xlzhang/jiseung/r11/truth/census3/job.sh`, instruments at
213843e6 in `tools3/`: `receipt_refusals.py` over the CUHK population;
the refined reader over archived and pin words; `resign.py` on Q24's
merged tree 002f91cf (uploaded as `code-002f91cf/`) over census 1's
specs, to place each CUHK pin change before or after that tree.
Expected: decision-gate refusals of receipts a recorded run minted
recur in R10 (the gate's run-stream reader knows five receipt kinds);
CUHK words at 002f91cf equal the pin's except for goals whose records
postdate what that tree can read.

## Status

- 2026-09-28: pin verified; briefs, kernel, CONDUCT, RSL, charter topics
  (settlement, goal grain, validity, analysis chain, delivery), CLAIMS,
  ARCHIVES, Q24's record and tools read. Trial: Q24's replay at the pin
  on the public po3-r19 runs in 2.5 s and keeps the state.
- Next: commit the instruments, run the local census (ax41 + public), pack
  and submit the CUHK census job.

## Jobs issued

None.
