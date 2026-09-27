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

## Jobs issued

- CUHK census E (slot job, 1 core, 8 GB, <= 1 h, reads only):
  `/project/xlzhang/jiseung/r11/evid/census/census.sh` sha256 a00ecf84...;
  `evidence_census.py` = `git show 223c66ed:.agents/research/loop/evidence_census.py`
  sha256 9ab91894...; `declarations-9185770e.json` (its
  `--dump-declarations` on this worktree at 223c66ed) sha256 3986d46b....
  Roots: r8, r9, r10/q1-q24, q26-q28, q30-q36 named one by one (q25 does
  not exist; q29 excluded), Q28's 32 pre-R8 directories.
  Submitted as CUHK Slurm 2157056 (r11-evid-a), pre-registration
  a4d88ca2304f; COMPLETED in 2 min 55 s, 45 MB.
- CUHK census E v2: the same roots and prunes, `census-v2.sh` sha256
  1dfef59d..., instrument `git show 4c38357b:...` sha256 6f005dc9...
  (the grouped-number tokenizer and witness fields; see its commit).
  Submitted as CUHK Slurm 2157060 (r11-evid-a), pre-registration
  d54aee39b738; COMPLETED in 34 s. Its derived outputs (rows,
  sessions, summary; 8 MB) were fetched to scratch for analysis; the
  records themselves were read in place with read-only greps.
- CUHK probe M (slot job, 1 core, 4 GB, <= 30 min, reads only):
  `/project/xlzhang/jiseung/r11/evid/probe-m/job.sh` sha256 e542f405...;
  probe = `git show 23c08541:.agents/research/loop/probe_mode_selection.py`
  sha256 a114fe1b...; code packed from 23c08541 (426 files, tree digest
  1371776d..., `chemsmart/` identical to 9185770e), unpacked at
  `/project/xlzhang/jiseung/r11/evid/code`, verified on the node.

## Census E -- READ (ax41 mirror + experiments-public local; CUHK 2157056/2157060)

Corpus (behavioural sessions; infrastructure -- zero provider turns or
`turn_deadline_exceeded` -- counted apart): ax41 + public 642 of 660
distinct (594 deepseek-v4-flash-0731, 48 qwen3.8-max); CUHK R8-R10 197
of 199 (all deepseek-v4-flash-0731); CUHK pre-R8 96 of 109. Pooled 935
sessions, 24,451 tool calls, 2,389 refusals. The instrument reproduces
Q28 exactly where they overlap (ax41: 1,237 authoring calls, hatches
in 9 sessions; CUHK: 82 hatch calls in 7 sessions -- none after q27).
Instrument corrected once after running (4c38357b: grouped numbers and
Unicode exponents were split; effect on ax41 delivered 4-digit numbers:
bound 1,092/1,252 -> 1,101/1,255); both versions ran on CUHK.

Predictions, as pre-registered (pooled; per-corpus where it differs):
- P1 HOLDS pooled: hatch sessions 18/935 (1.9%); per corpus ax41 1.4%,
  CUHK R8-R10 3.6% (fails alone), pre-R8 2.1% (PySCF numeric scf_tol).
  Untyped-intent hatch calls on 9185770e: Q28's 12 (Hirshfeld x3,
  named-site FlipSpin x8, NoUseSym x1); none added after Q28.
- P2 HOLDS: operative (non-prose) path/shell/native strings outside
  project sections in 13 of 24,451 calls (0.05%); 12 are
  `search_capabilities` query text my patterns mis-tag, 1 is
  `bind_scientific_identity(input_artifact_id="binol.xyz")`, refused
  (`artifact.id_is_registered`). No tool call reached a shell or native
  input outside the project sections Q28 already closed.
- P3 HOLDS: 311 selector refusals; 34 distinct requested names; 32
  declared somewhere on 9185770e; declared nowhere: `stationary_point_gradient`
  and xTB `solvation_shift_energy` (a naming slip for
  `xtb_solvation_shift_energy`). Refusals are dominated by deliberate
  semantic routing: ORCA-printed `gibbs_free_energy` on opt/ts x155
  (-> `derive_thermochemistry`), `spin_square` on closed shells x41 (->
  `multiplicity`); closed since: `solvent` x19, ORCA IRC trajectories
  x10, ORCA sp frequencies x16 (now an ORCA `freq` job type).
  Vocabulary refusals fell from 435 of 1,426 refusals (ax41) to 22 of
  563 (CUHK R8-R10).
- P4 FAILS: Agent-declared gaps (`blocked_unsupported`,
  `plan_unsupported_external`) in 108 of 861 analysis-planning sessions
  (12.5%). Read by reason, they are mostly not vocabulary: missing
  producers 32 workspaces, thermochemistry refused at a non-minimum 27,
  engine/envelope/execution 5, "other" 20 (read: evidence and
  execution). Vocabulary gaps: 23 workspaces over 11 types -- ORCA sp
  frequencies 9, frontier orbitals 3, IRC path geometry 3, atomic spin
  populations 2, and one each for energy->cm-1, scan coverage, NEB
  energies, Gaussian <R^2>, adaptive mode selection, categorical
  claims, per-root excited-state <S^2>; settings gaps: counterpoise
  ghost atoms 3 (open), functional MN15 3 (ORCA 6.1.1 itself refuses
  it); literature constants 7 (closed: `literature_constants.py` now
  carries the aqueous proton, SHE and ferrocene references with
  sources). On 9185770e every vocabulary type is closed except NEB
  energies (NEB is not executable), per-root excited-state <S^2> (a
  program limit), and adaptive mode selection (below).
- P5 HOLDS on its stated definition: re-extraction of a
  context-shown artifact in 301 sessions; expression or claim inputs
  citing a receipt no earlier reply of the session carried, refused:
  10 (all ax41; 0 in CUHK R8-R10). A further 16 same-tool refusals cite
  in-session digests of the wrong kind. Decisions citing digests that
  are not receipts (plan, binding, wake-shown digests) are refused with
  a route naming what the digest is: 54 CUHK, 25 ax41 -- citation
  hygiene, not chemistry.
- P6 HOLDS on its pre-registered definition ("a value in a typed tool
  reply"): of 5,039 unit-bearing numbers in delivered prose (final
  messages, finding and unreachable statements), 4,770 (94.7%) equal a
  host-returned value at the displayed precision (control 12.8%);
  prose-computed 2.0% (control 9.9%: the arithmetic class is not
  evidence), context-only 1.7%, own-argument + unmatched 1.6%. The
  stricter typed-record class is 63.4% (control 3.8%); at >= 4
  significant digits it is 88% (ax41 final messages), 86% (CUHK final
  messages) and 94% (CUHK finding statements: 216/229, control 0).
- P7 (reading): every unmatched delivered number (ax41 final messages
  51; CUHK findings 53 and final messages 38) and the ax41 residue with
  witnesses (149 rows, 101 sessions) read in context. The residue is
  correct prose restatements and unit conversions of host values,
  earlier-cycle claim values shown by the wake without receipts,
  labelled literature priors and expectation bands, and structural
  descriptions of input geometries. Exits that changed a delivered
  conclusion, as pre-registered: 3 goals --
  (a) r10/q17 hc1 `dE-h2occ-U1` (achieved): the final text says
      "-0.214269064 hartree (about -13.46 eV, -5.62 kcal/mol)"; the
      correct conversions are -5.83 eV and -134.46 kcal/mol; the claim
      of record is the correct host value (the sibling S1 wrote
      "-134.5 kcal/mol"). A wrong number delivered beside a right claim.
  (b) ax41 novel-round-3 `ino2-dinickel-exchange` (returned_to_human):
      J in the declared unit cm-1 "stated in prose using 1 kcal/mol =
      349.755 cm^-1" because energy->wavenumber was refused then; all
      six conversions check (e.g. 2.9850 x 349.755 = 1044.0). Closed
      since (`energy_to_wavenumber`).
  (c) CUHK r8/orca `goal-ts` (unreachable_from_evidence): IRC direction
      "delivered in prose from the typed extraction" because claims
      were numeric-only. Closed since (categorical answers).
  Plus one truncation (r10/q23 rt1: -56.50003 for -56.500035850 Eh,
  0.004 kcal/mol). Refusals for missing vocabulary prevented a
  requested observable in the 23 vocabulary-gap workspaces above; all
  but three types are closed on 9185770e.
- benchmark-v2 (pre-R10 baseline, 63 runs): its final messages carry
  almost no computed numbers (0-15 per run, none prose-computed); when
  the host lost a delivery (F7, the <S^2> crash, the arity crash) the
  Agent did not substitute prose numbers -- the referee's "reader told
  nothing".

Two behaviours read beside the counts:
- The Agent refuses to leave the typed layer when the value is in
  sight: r10/q10 lg1 saw Gaussian's printed <R**2> = 18.9424 in host
  text and rejected "reading the value by hand and binding it through a
  literal-bearing quantity expression ... that is fabricated evidence",
  declaring the observable unreachable (the selector gap is closed on
  9185770e).
- It composes selections inside the typed layer: the pre-R8
  `xtb-ir-acetamide-pyscf-stability` series found the PySCF mode nearest
  the strongest xTB band with `subtract` -> `abs` -> `min` and
  `coordinate_at_minimum(absdiff, frequencies)` (r3-r10:
  `coordinate_at_maximum` in 29 of 31 transcripts, `coordinate_at_minimum`
  in 25; corrected from a miscount "38 of 38" first written here).
  Argmax-then-select over aligned vectors is already expressible.

Numbers the model puts into the typed layer (host-recorded): 658 of
4,160 expression outputs (16%, 125 streams, ax41) depend on a
model-authored constant -- exponents, cutoffs, 298.15, R, kT/hc, a
Sackur-Tetrode constant, atomic masses, and literature values (an
experimental band 1746 cm-1; oxidant potentials 0.27 and 0.7 V). The
host names every one in the receipt and treats it as asserted for
sufficiency; the rendered report does not appear to mark a claim whose
value rests on one (open question, not yet checked on a CUHK record).

The one open vocabulary gap with a delivered loss: adaptive mode
selection. `vibrational_mode_atom_participation` is modes x atoms; `ref`
selects a mode's row by static index but never an atom's column, and
`sum/mean/max/min` take no axis. So "the mode with the most C+O
participation" cannot be planned before the Hessian exists. CORRECTED
(an inference first written here as fact, from operation counts): 46
ax41 sessions request participation, but what most did with it is not
established; the one opened (general-round h1) used its many `ref`s for
atom positions. Read instead: 8 transcripts in 7 workspaces index
participation with static `[mode, atom]` refs, mostly the imaginary mode
of a saddle, whose index is known before the Hessian exists. The need
fails exactly where the mode's index is not known in advance: r9/xtb g3
(exhausted), whose last wave had no further cycle (the graph's
`negative.an_analysis_chain_cannot_plan_a_mode_it_has_not_read`).

Retained numbers the Agent could not cite (`evidence_census_report.py
carried-literals`, ax41): of 67 literal values with >= 5 significant
digits that are not conditions or physical constants, 26 re-type a
number an earlier typed reply of the same session carried and 18 carry
one only the task or wake message showed -- earlier cycles' uncertainty
terms in ino3-r13b and ino3-r14a (nodes named `*-carried-v`), a cycle-one
Delta G in sm2-hcn-hnc, and two task-given experimental references. The
host then records them as model-authored. Claims standing on any
model-authored constant: 332 of 2,580 (12.9%, 86 streams; mostly counts,
exponents, cutoffs, 298.15, R). Rendered reports name model-authored
constants only in uncertainty rows (6 of 288 local reports), never beside
a claim value.

Where the three P7 cases live: the goal driver never reads a session's
`final_text` (only `live_session.py`, `loop.py` and the TUI do; the goal
CLI prints settlement and reasons), so in goals the final message is a
transcript surface, not a host-delivered one; all three cases are final
messages. Finding statements are delivered at goal grain and their
numbers are bound through `rests_on`: r10/q18 g1's restated 4.26/2.22
kcal/mol rest on `dst-proj > dst-unproj`, host claims -- so the
`context` class in findings restates bound claims, and the wake's dropped
receipts are a citation gap for expression inputs, not a delivery gap.

Gap typing after one fix (7982c0b3's reader types a "non-stationary"
reason first): ORCA sp frequencies 8 workspaces (not 9); r10/q21 g1-hooh
and r10/q24 g2r are the thermochemistry-at-a-held-dihedral gate. P2 with
`query` read as prose: 1 of 24,451 calls. q28-q36 (after Q28's native
gate): 15 sessions, 614 calls, 46 authoring calls, 0 hatch rows, 0
native-word refusals. The CUHK count of participation requests was
abandoned (a login-node grep over engine outputs; stopped, not counted).
Every figure in this section is recomputed by
`.agents/research/loop/evidence_census_report.py` (7982c0b3) from the
census outputs.

Route shapes (C5 by-product): 538 accepted plans with calculation
nodes -- 462 single-program (ORCA 364, PySCF 68, xTB 29, Gaussian 1),
76 cross-program (PySCF+xTB 23, ORCA+PySCF 19, ORCA+xTB 17,
Gaussian+ORCA 9, three programs 8); 1 node 162, 2-4 nodes 267, 5+ 109.
248 goals executed nodes: 208 single-program, 40 cross-program
(ORCA+PySCF 26); 264 geometry handoffs and 269 data edges recorded.

## POSITION MEMO (Phase I)

**The question, as I now understand it.** Not "is the typed layer
complete" but: when the Agent's chemistry reaches a human, does it pass
through host values, and where it does not, is that because the
vocabulary lacked something, because a receipt could not travel, or
because the Agent chose prose? The delivery surfaces differ by era: in
goals the host delivers claims, categorical answers and finding
statements (bound through `rests_on`); a session's final message is
delivered only when a person runs a session directly (CLI, TUI).

**Evidence** (all model-labelled; deepseek-v4-flash-0731 in 886 of 935
sessions, deepseek-v4.1-flash in 1, qwen3.8-max in 48): census E over
935 behavioural sessions,
24,451 tool calls, 2,389 refusals, ax41 mirror + experiments-public +
CUHK R8-R10 and pre-R8 (Slurm 2157056, 2157060), with a
shuffled-session control; the reading of every unmatched delivered
number and of the Agent's 108 declared-gap sessions by reason.
DEVIATION from the pre-registration: P7 promised to read "every V/R
refusal in a delivered goal"; the 311 selector refusals and the reuse
refusals were read by class and host diagnosis, with a few opened, not
one by one;
benchmark-v2's referee report as the pre-R10 baseline. Instruments:
`evidence_census.py` (223c66ed, 4c38357b, 7982c0b3) and
`evidence_census_report.py` (7982c0b3).

**What it says about C2: NARROWED** (by the pre-registered map: P7 = 3,
P4 failed).
- Supported where the claim is strongest today: delivered numbers are
  host values -- at >= 4 significant digits the typed-record class is
  86% of final-message numbers and 94% of finding-statement numbers
  (216/229; control 0/229); the looser pre-registered class is 94.7% of
  all 5,039 (control 12.8%). Outside project sections, one operative
  exit in 24,451 calls (a file name given as an artifact id), refused;
  inside them, Q28's 109 hatch calls in 18 sessions all predate its gate
  and none follow it (0 in 46 authoring calls, 15 sessions -- a small
  n). When the value was in sight the Agent refused to type it
  by hand (r10/q10 lg1), and it composed selections inside the DAG
  (xtb-ir `coordinate_at_minimum`).
- Narrowed to: the archive's exits that changed a delivered conclusion
  are three final messages -- two where a missing typed form forced
  prose (energy->wavenumber, ino2; categorical answer, r8 goal-ts; both
  forms built since, both restate typed values -- ino2's six
  conversions checked; goal-ts's direction is a word restated from its
  typed extraction, not checked against the trajectory) and one wrong prose
  conversion beside a correct claim (r10/q17 dE-h2occ-U1). Vocabulary
  refusals prevented a requested observable in 23 workspaces over 11
  types; on 9185770e all are closed except NEB energies (not
  executable), per-root excited-state <S^2> (a program limit), and
  adaptive mode selection (one delivered loss, r9/xtb g3).
- Failed prediction P4, and what it hid: the Agent declares a gap in
  12.5% of analysis-planning sessions, but by reason (post hoc) most are
  missing producers and host physics gates, not vocabulary. The
  declaration mechanism is how the Agent keeps scientific intent honest;
  it is not a measure of vocabulary lack.
- The retained-number gap is real and small: receipts cannot be
  expression inputs across runs (10 refusals, all pre-R10 ax41), and the
  Agent's workaround re-types earlier results as literals (18 carried, 5
  workspaces), turning host numbers into model-authored ones. Findings
  are not affected (bound at goal grain).
- Widening stays open: I measured no route change after a refused read
  in exploration; two anecdotes point both ways (r9/xtb g3 blocked;
  xtb-ir composed).

**Program proposed next.**
1. No implementation for the delivery path. Three cases in ~5,000
   delivered numbers, two already closed by typed forms, earn no gate;
   a host check on prose numbers is not proposed.
2. Witness-first candidate, not a commitment: an axis on
   `sum/mean/max/min` (or a column select in `ref`) for matrix
   quantities (participation, populations, bond orders). It replaces a
   wake cycle, not a tool. First probe, provider-free: plan the C=O
   selection against r9/xtb g3's archived xTB Hessian and check the
   typed chain returns the mode the Agent chose by eye (mode 99, C10+O1
   0.969). Whether the model would plan it is Q17's use-not-sight
   question and belongs to a later live check.
3. Receipt reuse across runs: evidence for it is the carried literals
   (ino3-r13b/r14a uncertainty budgets); the smallest home is letting
   `_typed_quantity_from_receipt` resolve a recorded run's receipt as
   decisions already do (`_recorded_run_receipt`, Truth's hotspot).
   Proposed only if Truth finds a sufficiency word changed by a carried
   literal; otherwise leave it.
4. Routed to the master, not ruled on: Gibbs energy refused at a held,
   non-stationary dihedral (r10/q21 g1-hooh, r10/q24 g2r) -- a possible
   C3 suppression, since projected-Hessian thermochemistry at a
   constrained point is a legitimate route; ORCA MDCI launched with more
   MPI ranks than electron pairs (r10/q12 g1-hi; execution contract);
   rendered reports never mark a claim value that rests on a
   model-authored literal (Truth, host words; load 332/2,580 claims, the
   literature-valued ones few); counterpoise ghost atoms have no typed
   form (open settings gap; a designed boundary in benchmark-v2 b7).

**What would change my mind.** A host-delivered claim, answer or finding
whose number disagrees with its receipt on a tree after R10 (narrows
C2 further toward contradiction); the provider-free probe failing to
return the Agent's mode (the axis reduction is not the missing form);
a live check where the Agent, given the form, still does not plan the
selection (use, not sight: no implementation); a second model's exit
rates differing materially from deepseek's (C2 becomes model-specific --
48 qwen sessions are too few to say); a measure of route change after
refused exploratory reads (would decide widening).

## Probe M -- PRE-REGISTRATION (master's instruction after the memo; provider-free only)

Question (the master's): can a typed plan, written before the Hessian
exists, select a vibrational mode by its atom participation and return
mode 99 on R9 xtb g3's archived Hessian? If today's vocabulary
expresses it, no code; if not, the smallest change and what it replaces
(the extra wake cycle). No implementation beyond a witness before
Exchange 1.

Inputs, read in place and read-only on CUHK (`/project/xlzhang/jiseung/r9/xtb/goals/g3/workspace/`):
- route 2: `cycle-2/gf2_hess_gf2min_r2/.../geom-gf2min-cyc1_hess.out`
  (GFN2-xTB Hessian at the GFN2 minimum; the one the Agent read);
- route 1: `cycle-3/gf2_hess_gffgeom_r3/.../gff-relaxed-geometry-r3_hess.out`
  (GFN2-xTB Hessian at the GFN-FF geometry; the last wave, whose
  assignment the goal could not plan).
Atom sets, 0-based, exactly the indices the Agent's own `ref` nodes used
(`[96..98, 1|2|4|10|14|16]` in live-20260921T103342...): lactam C=O
{10, 1}, amide C=O {14, 2}, carboxyl C=O {16, 4}.

Instrument: `.agents/research/loop/probe_mode_selection.py` (committed
before the run). It drives `CommandCompiledToolHostV1.dispatch`
(`extract_result_quantities`, `evaluate_quantity_expression`), the
public tool surface. Its mechanics were checked locally on
`tests/data/XTBTests/outputs/acetaldehyde_hess` (not g3): C=O {2, 0} ->
mode 11, 1798.58 cm-1, share 0.980; runner-up share 0.556 (mode 2, 501.81
cm-1). Every blind construction was refused or reduced to a scalar there.
Tree: this worktree's committed `chemsmart/` (identical to 9185770e),
packed by `pack_code.sh`, digest re-verified on the node; HOME fenced to
the job directory.

Oracle: the Agent's reading of route 2 in its transcript --
mode 99 = 1866.33 cm-1, lactam, 0.649 + 0.320 = 0.969; mode 97 = 1752.42,
amide, 0.634 + 0.346 = 0.980; mode 98 = 1778.55, carboxyl; "an all-real
117-mode spectrum (min 10.68 cm-1)".

Predictions (never tuned after a result):
- M0 (replay gate): route 2 extracts 117 modes x 41 atoms; symbols at
  {10, 1}, {14, 2} and {16, 4} are C and O; the static-index route
  returns 0.969 (3 decimals) and 1866.33 cm-1 for row 98, and 0.980 and
  1752.42 for row 96. If M0 fails, the probe is void (reader drift since
  R9), and nothing after it is read.
- M1 (today's vocabulary): every blind construction is refused or
  non-selecting on both routes. B1 `coordinate_at_maximum` over the
  table and B2 `multiply` by an atom mask raise shape errors; B3 (column
  select in `ref`, as a null index and as an `axis` field) is refused by
  the tool schema; B4 `sum`/`max` reduce to a scalar. Falsified if any
  construction returns the target frequency.
- M2 (the criterion): on route 2 the argmax of summed participation
  selects row 98 (mode 99, 1866.33) for {10, 1}, row 96 (1752.42) for
  {14, 2} and row 97 (1778.55) for {16, 4}. The Agent's by-eye choice is
  exactly this argmax. Falsified if any row differs.
- M3 (route 1, what the goal could not deliver): each set selects one
  mode with share >= 0.85 and runner-up share <= 0.15. If falsified,
  argmax alone is not a safe selector at a non-stationary geometry, and
  the typed form must also carry the margin.
- M4 (the wake-cycle route): the static-index route at route 1's argmax
  rows returns typed receipts equal to the targets (share to 1e-12,
  frequency exact).

Outcome -> proposal:
- M1 and M2 hold: today's vocabulary does not select the mode at plan
  time. The only formal route, reasoned from code and not run, enumerates
  every mode: for cephalexin about 590 expression nodes over at least five
  stages (128 nodes per expression, 64 inputs per node) plus one
  validation node per mode as a selector, 116 of whose verdicts would read
  "failed". I do not count that as expressing it. The smallest change is
  a column select in `ref` (null in `indices` meaning every element along
  that axis). With `add` and `coordinate_at_maximum` it makes the
  selection a four-node plan written before the Hessian, and it replaces
  the wake cycle that reads the table and then plans static refs.
- M1 fails: no code; report the construction.
- M2 fails: report what the Agent chose instead of the argmax.

## Probe M -- READ (CUHK Slurm 2157069, pre-registration d3b4787f1946)

Run: COMPLETED in 12 s on chpc-cn071. It imported the packed tree
(`/project/xlzhang/jiseung/r11/evid/code`), remote digest 1371776d =
local, commit 23c08541. Inputs were hashed in place (route 2 `.out`
6c78dc46..., route 1 `.out` 3017d95a...). Output:
`/project/xlzhang/jiseung/r11/evid/probe-m/probe.jsonl` (27 rows).

ERROR FIRST -- my pre-registration, not the reader. M0 and M2 named
rows 98/96/97 because I read the Agent's "mode 99" as 1-based. The
Agent labelled modes 0-based, as it did atoms: its `ref` nodes read rows
96-99, and on the R9 tree `ref [99,10]` returned 0.6490384615384617 and
`ref [99,1]` 0.3200120192307692 (transcript
live-20260921T103342..., nodes r99c10, r99o1). Today's reader returns
those identical values at the identical rows: lactam row 99, carboxyl
row 98 (C 0.6094017956219107), amide row 97 (C 0.6339746132487108).
So M0 FAILED AS WRITTEN. Its purpose, detecting reader drift, is decided
by the Agent's own receipts, and it shows none. I read M1-M4 after that
diagnosis; whether the gate stands is the master's call.

- M0 in substance: route 2 extracts 117 modes x 41 atoms (lowest 10.68
  cm-1, as the Agent recorded); every set's atoms are C and O; the static
  route returns the Agent's values to the last digit.
- M1 HOLDS on both routes: B1 `coordinate_at_maximum` over the table
  raises "needs one coordinate per value; got 4797 values and 117
  coordinates"; B2 `multiply` by an atom mask raises "multiply accepts
  scalar broadcasting or identical shapes"; B3a/B3b are refused by the
  tool schema ("nodes[0].indices[0] must be integer, but got NoneType
  None"; "supplied ['axis'], which this object does not accept"); B4
  `sum` returns 117.0 and `max` 0.943/0.941, selecting nothing.
- M2 HOLDS in substance: the argmax of summed participation is exactly
  the Agent's by-eye choice -- lactam {10,1} row 99, 1866.33 cm-1, 0.969;
  amide {14,2} row 97, 1752.42, 0.980; carboxyl {16,4} row 98, 1778.55,
  0.962.
- M3 FALSIFIED as written. All six selected shares are >= 0.96, but the
  runner-up shares are 0.36-0.63 on BOTH routes (route 1: 0.364 at
  136.58 cm-1, 0.633 at 71.34, 0.516 at 743.52; route 2: 0.358, 0.358,
  0.494). Low-frequency group torsions and bends load the same C and O,
  so an atom pair's share is not small for non-stretch modes. The argmax
  still separates by >= 0.33. As pre-registered, the typed form must
  also carry the margin. The instability is not specific to the
  non-stationary geometry.
- M4 HOLDS: the static-index route at route 1's argmax rows returns typed
  receipts equal to the targets (shares identical as floats, frequencies
  exact).

What the probe also showed:
- The Agent was right to refuse index reuse across the two Hessians.
  Amide and carboxyl swap rows: route 2 has amide at 97 and carboxyl at
  98; route 1 has carboxyl at 97 (1699.03) and amide at 98 (1714.42). Route
  2's amide index reused on route 1 would have delivered the carboxyl
  stretch as the amide. Selection by participation is right on both.
- Route 1's assignment, the one the goal could not deliver (probe output,
  not a delivery). GFN2-xTB Hessian at the GFN-FF geometry, three
  imaginary modes (-56.06, -35.65, -22.09 cm-1): lactam 1754.65 cm-1
  (share 0.981), amide 1714.42 (0.972), carboxyl 1699.03 (0.964). The
  lactam-amide separation is 40.2 cm-1 against 113.9 at the GFN2 minimum:
  curvature at a point that is not stationary.

Answer to the master's question:
- Today's vocabulary does not select the mode at plan time (M1). The
  only formal route is the per-mode enumeration reasoned from code
  (about 590 expression nodes and 117 validation verdicts for
  cephalexin), not run and not counted.
- The smallest change is a column select in `ref`: null in `indices`
  means every element along that axis. The plan is then written before
  the Hessian exists, in five nodes: `ref(part,[null,10])`,
  `ref(part,[null,1])`, `add`, `coordinate_at_maximum(share, freq)` for
  the frequency, and `max(share)` for its share. It has no model-authored
  literal; the atom indices are identity facts of the bound geometry.
  It replaces the wake cycle between a Hessian and its assignment, which
  g3 did not have. An axis on `sum/mean/max/min` alone does not do it: it
  also needs vector-matrix broadcasting and a literal atom mask
  (model-authored) -- two changes and weaker provenance.
- Size of the change, from the code read: `ref` indexes sequentially
  (`selected = selected[index]` per index, quantity_expressions.py
  ~1648), so a null cannot be a slice in that loop. The change touches
  the evaluator (tuple indexing with a slice for null) and the tool
  schema (`indices` items integer or null), not the schema alone. A null
  index is preferred to an `axis` field because it generalises to any
  axis, and to mixed selections, with no second parameter.
- The margin (M3) is a separate, optional candidate: the host names the
  second extremum beside every `coordinate_at_maximum/minimum` receipt
  as an observation, the pattern it already uses for geometry over
  non-bonded atoms. No new operation. A second observation for
  Exchange 1, not a proposal: every runner-up is a low-frequency group
  motion (42-744 cm-1). A chemist would exclude those by a frequency
  window, and today's vocabulary cannot express that either (no
  elementwise comparison).
- The archive was read, not written: no file in either g3 Hessian
  directory is newer than the job's start.
- Witness: `probe_mode_selection.py` step `candidate_plan`. It is red on
  9185770e (the schema refuses the null index; shown locally on the
  acetaldehyde fixture after 2157069, whose B3a is the one-column form of
  the same refusal on g3). It turns green when it returns these targets:
  - route 2: lactam 1866.33 / 0.9690504807692308, amide 1752.42 /
    0.9799682665608885, carboxyl 1778.55 / 0.9624735196207002;
  - route 1: lactam 1754.65 / 0.9813981398139814, amide 1714.42 /
    0.971859296482412, carboxyl 1699.03 / 0.9640676608948053;
  - acetaldehyde C=O {2,0}: 1798.58 / 0.9803607214428858.
  Nothing is implemented; the decision is Exchange 1's.

## Status

- 2026-09-28: base verified; governance, archives, Q28 records and the
  evidence surfaces read; census E pre-registered before any detector
  run; instrument committed (223c66ed, corrected 4c38357b, 7982c0b3) and
  run on the local corpus and on CUHK (2157056, 2157060); census read;
  position memo written. Master's exchange message received: C2
  "narrowed" logged as provisional. Probe M pre-registered, run (CUHK
  2157069) and read above; witness defined, nothing implemented. Waiting
  for Exchange 1.
