---
id: verify-when-signing
paths:
  - "chemsmart/agent/driver.py"
  - "chemsmart/agent/goal.py"
  - "chemsmart/agent/tool_runtime.py"
conditions: "the goal driver's run-path settlement and the planning-time refusal verification on the tree of 2026-09-24 (c79c39a1); CUHK Gaussian 16 and archived ax41/CUHK ledgers; R10 episode Q10"
evidence:
  - "commit:51ff5466"
  - "commit:a5cb66b7"
  - "test:tests/agent/test_a_host_word_is_true_of_what_it_read.py::test_a_run_without_an_analysis_chain_certifies_nothing"
  - "test:tests/agent/test_a_host_word_is_true_of_what_it_read.py::test_a_refusal_made_before_the_run_is_read_again_when_the_goal_settles"
  - "note:CUHK Slurm 2151662 (LG1): the base tree signs unreachable_from_evidence over a log that prints the refused value"
repeat_cost: "six false `achieved` settlements, R9's own pushed smoke goals among them, from a settlement that read a chainless run's empty stream as 'nothing undelivered'; and a refusal verified at planning, before any result existed, signed after the run had printed the value"
falsifier: "retire when every host-signed word is recomputed from the evidence at signing time by construction; wrong if a word computed from an earlier check is shown to stay true after the evidence changes"
home: prose
supersedes: []
earned: 2026-09-24
last_verified: "2026-09-28 @ 04c2e5aa"
---
A word the host signs -- a settlement, a verified refusal, a certification -- is computed from the evidence as it stands when the word is signed, and says what it read.
A check made earlier (at planning, in another cycle, over another run's stream) describes that moment only; carrying it into a later word is how a verified word goes false.
Evidence: 51ff5466 (a chainless run certified what it never claimed: 6 false achieved) and a5cb66b7 (a planning-time refusal signed over a log that printed the value; LG1, CUHK 2151662).
