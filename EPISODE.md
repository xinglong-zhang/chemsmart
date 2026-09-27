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

## Status

- 2026-09-28: pin verified; briefs, kernel, CONDUCT, RSL, charter topics
  (settlement, goal grain, validity, analysis chain, delivery), CLAIMS,
  ARCHIVES, Q24's record and tools read. Trial: Q24's replay at the pin
  on the public po3-r19 runs in 2.5 s and keeps the state.
- Next: commit the instruments, run the local census (ax41 + public), pack
  and submit the CUHK census job.

## Jobs issued

None.
