# R11 episode evid -- the evidence surface (lens `evidence`)

Base SHA: 9185770e30e486e0d1e96e83bf35edc914b8590a (verified with
`git rev-parse HEAD` as the first action, 2026-09-28). Model:
claude-opus-5-5[1m]. Episode id `evid`, lens `evidence`. Brief sha256
prefix 73d92f085eb7dd1c (`brief-evid.md`, matches
`SHA256-prefixes.txt`); common brief `common-evid.md` prefix
98d71df1938e5547.

## The question (as currently understood)

Where, in real archived CHEMSMART Agent sessions, does the Agent's
scientific work leave the typed evidence layer -- native program words,
numbers computed in prose, raw text, requests refused for lack of
vocabulary, lost lineage, receipts it could not reuse -- and which of
those exits changed a delivered conclusion? What is the smallest change
that closes the load-bearing exits, or is none warranted?

Claim served: C2 (the Agent does its chemistry in typed evidence, with
provenance, without program-native language). May be narrowed,
contradicted, widened or replaced.

## Falsifiers (pre-registered before any census is run)

- C2 is **supported** if, across the archived authoring and analysis
  calls with stated denominators, exits from the typed layer are rare
  AND reading the transcripts behind each exit shows none changed a
  delivered conclusion (no delivered number unbound by a receipt that
  differs from what a receipt would have given; no refused request whose
  absence changed the answer).
- C2 is **narrowed** if exits exist and some are load-bearing, but they
  concentrate in a small, nameable set of missing typed forms.
- C2 is **contradicted** if delivered conclusions routinely rest on
  numbers the typed layer did not produce (prose arithmetic, raw text,
  native words), i.e. the typed layer is bypassed for the load-bearing
  step.
- C2 is **widened** if the load-bearing exits are in exploration
  (deciding what to compute next) rather than delivery.

## Oracle

Host records (ledgers, run streams, receipts, settlement records) and
the archived transcripts themselves, read through the host's own
functions where a replay is needed, on the commit that produced them.

## What orientation established (read, no detector run yet)

- Q28's census (`~/.chemsmart-r11-run/instruments/q28/`): all 17 copies
  match `SHA256SUMS`; `census.py` is 423e919b..., the corrected
  instrument of CUHK 2153720. Q28's "1,788 authoring calls" = 1,237
  distinct ax41 + 551 CUHK R8-R10 (EPISODE.md at cbdde295); 104 hatch
  calls in 16 of 842 sessions.
- Local corpus (walked with private stores pruned: `sealed*`, `seals`,
  `grading`, `grade-*`, `claude`): ax41 mirror 906 public transcripts =
  660 distinct sessions (sha256), 1,538 event streams, 195 goal ledgers;
  `experiments-public` 9 transcripts, 3 ledgers. Tool replies 15,083, of
  which 1,426 `rejected`.
- Code on 9185770e (priors verified by reading, not yet by records):
  `render_workspace_record` (workspace_record.py 979-1064) drops
  `claim_receipt_sha256`/`source_receipt_sha256` that claim rows store
  (543-548) and omits `finding` and `answer` rows entirely; expression
  inputs and claims resolve only in-session receipts
  (`_typed_quantity_from_receipt` 20000-20057, `record_analysis_claims`
  ~21027, both a bare `ContractError` naming no route), while decisions,
  refusals and dispositions accept a recorded run's receipt
  (`_recorded_run_receipt` 7262-7302).
- What a human reads: the goal CLI prints only settlement and reasons;
  claims and categorical answers are host-rendered (executor
  `completed-/partial-analysis-report.md`); finding statements are
  model prose "recorded as yours and never checked" (only their
  `rests_on` relations are); a session's `final_text` is the model's
  last message unless the task carries an `analysis_completion_policy`,
  in which case loop.py (~676-690) replaces it with the host report.
- Excluded by design: Mac `~/.chemsmart/agent/sessions/` (300 entries,
  May-July 2026, the v8 plan-step agent -- `build_molecule`,
  `recommend_method` -- not the typed layer); CUHK `r10/q29` (sealed
  study being closed), `r10/q3/sealed`, `r10/q6/goals`,
  `r10/q17/sealed*`, never `r10/m*`.

## Census E -- PRE-REGISTRATION (written before any detector has run)

Instrument: `.agents/research/loop/evidence_census.py` (new, stdlib
only, imports nothing from chemsmart; Q28's `census.py` logic adopted
verbatim inside it). Rows stay in scratch (they carry private paths);
numbers come here. Distinct sessions by `transcript_sha256`; each
session labelled with its observed model (`provider_turn_observed`);
sessions with zero provider turns or a `turn_deadline_exceeded` ending
are infrastructure, counted apart and never in a denominator.
Corpus: ax41 mirror (both slices), `experiments-public`, CUHK `r8`,
`r9`, `r10/<q>` named one by one (q29 excluded), Q28's 32 pre-R8
directories (a slot job, read-only).

Detectors and classes:
- X (native, path, shell): Q28's hatch keys and unknown-key refusals
  unchanged; extended to every string leaf of every tool call --
  native input (ORCA `!`/`%block ... end`, Gaussian `#` routes and
  `%chk/%mem/%nproc`, `$`-blocks, PySCF code), paths (absolute,
  `~/`, `../`, native output names), shell (pipes, `&&`, `$(`,
  `bash/cat/grep/sed/python`, `chemsmart run|sub`, `sbatch`); the reply
  (accepted / refused) recorded.
- V (vocabulary refusals): rejected replies classified by message into
  selector (undeclared for jobtype / not in enum / reader does not
  provide / value absent), operation or field not accepted, unit,
  thermochemistry kind, shape. Each requested selector is looked up
  against 9185770e's reader declarations: declared for this program
  and another jobtype, declared for another program, or declared
  nowhere (a gap still open today). Plus the Agent's own declarations
  of a gap: analysis nodes planned `blocked_unsupported` with their
  `blocked_reason`, and `plan_unsupported_external`.
- R (reuse): each cited receipt digest or artifact id is placed as
  minted in-session, shown only by the task/wake message (context), or
  neither; with the reply. Re-extraction of a context-shown artifact
  is the legitimate route and is counted as reuse that succeeded.
- N (numbers): targets are unit-bearing numbers (energy, frequency,
  distance, angle, wavelength, dipole, shift units; not K, %, bare
  integers) with >= 2 significant digits, in delivered prose -- finding
  statements, unreachable-observable statements, and each session's
  last assistant message -- and, reported apart as exploration,
  decision text fields and mid-session assistant text. Classes, first
  match wins: bound (equals a value in a typed tool reply of the same
  session at the displayed precision), context (equals a value in the
  task/wake message only), own-argument (equals a number the model
  itself put in an earlier tool argument), prose-computed (one unit
  factor applied to, or a difference/sum/ratio of, seen values),
  unmatched. Control: the same matcher against a different session's
  seen set (fixed derangement); a class counts as evidence only where
  its real rate is at least twice its control rate.
- S (route shapes, for C5, no prediction): per session the accepted
  plans' (program, jobtype) nodes, analysis kinds and edge kinds; per
  goal the executed (program, jobtype) and handoff/data-edge counts.

Predictions (never tuned after a result):
- P1: sessions carrying any Q28 hatch or free-word key <= 3% of
  distinct sessions; hatch calls whose intent has no typed form on
  9185770e <= 15.
- P2: path, shell or native-input strings outside project sections in
  <= 0.5% of tool calls; none reached a native input or a shell.
- P3: selector refusals >= 200 calls; >= 60% of the distinct requested
  selector names are declared somewhere on 9185770e; names declared
  nowhere <= 10.
- P4: Agent-declared gaps (blocked_unsupported, plan_unsupported_external)
  in <= 5% of sessions that plan analysis; <= 10 distinct gaps, among
  them mode vectors / per-mode participation.
- P5: re-extraction of a context-shown artifact in >= 20 sessions;
  refused expression or claim inputs citing a receipt not minted
  in-session <= 20 calls.
- P6: of delivered-prose target numbers, bound >= 70%, prose-computed
  <= 15%, context-only <= 10%, own-argument + unmatched <= 15%.
- P7: after reading every flagged (prose-computed, context, unmatched)
  number and every V/R refusal in a delivered goal (achieved,
  achieved_with_observations, unreachable_from_evidence), at most 3
  goals where an exit changed a delivered conclusion -- a wrong number
  delivered, or a conclusion resting on a computation or a read the
  typed layer could not do.

What each outcome does to C2:
- supported: P6 bound >= 70%, P7 <= 1 goal, P3 gaps none load-bearing.
- narrowed: P7 finds 2-5 load-bearing goals concentrated in <= 3
  missing typed forms (each named with its smallest typed home and
  what it replaces).
- contradicted: bound < 50% of delivered numbers, or >= 6 goals whose
  delivered conclusion rests on untyped computation.
- widened: exits concentrate in exploration (refused reads or prose
  arithmetic that changed the route taken) rather than in delivery.

## Status

- 2026-09-28: base verified; governance, archives, Q28 records and the
  evidence surfaces read; census E pre-registered (above) before any
  detector run. Next: write the instrument, run it on the local
  corpus, then the CUHK slot job.
