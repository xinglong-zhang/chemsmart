---
id: verify-when-signing
paths:
  - "chemsmart/agent/driver.py"
  - "chemsmart/agent/goal.py"
  - "chemsmart/agent/tool_runtime.py"
conditions: "the goal driver's settlement, the planning-time refusal verification and the decision gate, on the trees of 2026-09-24 (c79c39a1) and 2026-09-28 (R11: 254 archived goals re-signed, CUHK and ax41 ledgers); CUHK Gaussian 16; R10 episode Q10 and R11 truth"
evidence:
  - "commit:51ff5466"
  - "commit:a5cb66b7"
  - "commit:4ec8957d"
  - "commit:10f09617"
  - "commit:945186b0"
  - "test:tests/agent/test_a_host_word_is_true_of_what_it_read.py::test_a_run_without_an_analysis_chain_certifies_nothing"
  - "test:tests/agent/test_a_host_word_is_true_of_what_it_read.py::test_a_refusal_made_before_the_run_is_read_again_when_the_goal_settles"
  - "test:tests/agent/test_a_goal_word_names_every_expectation_the_physics_left.py::test_a_goal_word_names_the_expectations_an_earlier_cycle_scored"
  - "test:tests/agent/test_a_goal_word_reads_every_criterion_the_goal_holds.py::test_a_run_word_reads_the_criteria_and_decisions_of_the_whole_goal"
  - "test:tests/agent/test_a_decision_gate_names_what_the_host_minted.py::test_a_woken_decision_may_cite_what_an_earlier_session_minted"
  - "note:CUHK Slurm 2151662 (LG1): the base tree signs unreachable_from_evidence over a log that prints the refused value"
repeat_cost: "six false `achieved` settlements, R9's own pushed smoke goals among them, from a settlement that read a chainless run's empty stream as 'nothing undelivered'; a refusal verified at planning, before any result existed, signed after the run had printed the value; and in R11 three words signed from fewer records than their question covers -- a goal's expectations and its run-path criteria read from one stream (r10/q7 signed plain achieved over two falsified expectations), and a decision gate that read goal runs only and told the Agent 29 receipts the host had minted were not its own"
falsifier: "retire when every host-signed word is recomputed at signing over every record its question covers, by construction; wrong if a word computed from an earlier check, or over a narrower scope, is shown to stay true after the evidence, or the wider records, change"
home: prose
supersedes: []
earned: 2026-09-24
last_verified: "2026-09-28 @ 629b5113"
---
A word the host signs -- a settlement, a verified refusal, a certification, a gate's diagnosis -- is computed from the evidence as it stands when the word is signed, over every record its question covers, and says what it read.
A check made earlier (at planning, in another cycle) describes that moment only, and a signer that reads only what it holds (one stream of a goal, one host's decisions, runs but not sessions) describes only that; carrying either into the word is how a verified word goes false.
Evidence: 51ff5466 (a chainless run certified what it never claimed: 6 false achieved), a5cb66b7 (a planning-time refusal signed over a log that printed the value; LG1, CUHK 2151662), 4ec8957d and 10f09617 (a goal's word read one stream's expectations and criteria), 945186b0 (a gate that read goal runs only called 29 host-minted receipts unminted).
