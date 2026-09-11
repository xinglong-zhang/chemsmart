# ChemSmart Development Conduct

This file binds every session that changes this repository, human or
model. It says what a test may exist for, what a gate may be, how large
a change may be, and how a change is proven. `AGENTS.md` is the product
charter; this is the working discipline beneath it. When the two
disagree, the charter wins and this file is wrong.

## 0. The values, in the order they win

When two of these conflict, the earlier one wins, and the conflict is
written down.

1. **The host's word is true.** Every verdict, settlement, match and
   refusal is backed by a receipt and means what the physics means. A
   word that can be false for a real structure is a defect: `failed`
   for a quasi-planar radical whose inversion barrier lies below the
   zero-point level, `achieved` over an undelivered headline, `agreed`
   by dimension alone.
2. **One human decision per goal.** Displayed, one-shot, digest-bound;
   a delegated approval is disclosed on every record it touches.
3. **Replication before belief.** An anomaly, a number or a claim
   stands after it reproduces under a stated perturbation; N runs are
   N observations, and a weak run is never re-rolled.
4. **Freedom of route.** Any chemically valid route, decomposition or
   interpretation is admissible; a gate exists only where language
   cannot compute an invariant, and a refusal names the route.
5. **The anomaly has standing.** What the host detects it records with
   the numbers that tripped it, whether or not it was asked for; the
   model interprets, the human judges, and a gate is never loosened to
   buy a discovery.
6. **Errors first, receipts always, the seal is the method.**

## 1. What a test may exist for

A test exists for exactly one of three reasons:

- **Heartbeat.** The live `chemsmart` CLI compiles and safely previews
  a program and jobtype through project YAML, and a written native
  input reads back to what was requested.
- **Reachability.** A tool the Agent can call reaches the ChemSmart
  function it claims: every advertised parameter is settable, written,
  and read back; every declared selector is requestable; every operation
  in the vocabulary is exposed.
- **Production path.** A release-qualified path runs end to end on
  archived real program output: the approval chain, the provider-free
  executor, the result readers, unit and dimension arithmetic, the goal
  driver's settlements, and one pin per code gate that protects a
  scientific or authority invariant.

A test that pins wording, a single observed stream shape, a private
helper, or a historical defect a general invariant already covers is not
written; if it exists, it is deleted. A test says what it pins with a
``capability(...)`` marker (``kind:id``, a trailing ``*`` for a kind),
which is how the capability ladder learns that a capability is tested. A
defect earns one general test at the invariant it broke, never a test of
the case. Tests verify mechanics; a real observation through the public
surface establishes behaviour, and no test is ever cited as engine
execution.

## 2. What a gate may be

A refusal lives in code only when it protects an invariant that
language cannot compute or that an optimising model would erode:

- molecular identity, atom order, and geometry lineage;
- explicit charge and multiplicity, and the arithmetic of impossible
  states;
- the single human decision: the digest-bound one-shot bundle and the
  equality of recompiled argv with reviewed argv;
- no model-authored native input, path, shell, or status;
- the execution envelope's budgets and the dispatch target;
- units and dimensions across the analysis DAG, and the per-jobtype
  meaning of a selector;
- host-rendered claims and completion receipts;
- credentials;
- the terminal-state vocabulary and the stationary-point rule;
- the observable-regression guard, because deleting the node that
  carries a finding is the cheapest way to clear it.

Everything else is a sentence at the point of use: on the tool whose
argument it governs, in the wake context that carries it, or in a guide
the host opens. A sentence is a registered rule with an id, a placement,
and the provenance that earned it (``chemsmart/agent/rules.py``); prose
that lives nowhere else is not a rule. Before a gate is added, its
reachability from the model is measured by a direct probe; a check the
host already normalises away is not a gate and gets no sentence. A gate
earned by a live loss names that loss in a comment.

A refusal is a first-class output, and writing one is a design act.
Measured over one round: the agent met 470 refusals, retried the
refused tool in 60% of them, recovered within the session in 82%, and
opened a guide in 2% -- the refusal message, not the guide, is where
the model is taught. Three routes changed this round because a refusal
named them. So a refusal states the invariant it protects and, where a
legal route exists, names it; it does not script the science that
follows. And it is never the mechanism by which an unexpected result
is discarded: what contradicts an expectation is delivered as itself,
and a refusal that would bury a finding is a defect in the refusal, not
a fact about the finding. The live tension to watch is the quarantine
of a result typed ``failed``: its numbers are exactly where inversion
transition states and other surprises live, and today only the anomaly
receipt carries them.

## 3. How large a change may be

- One general commit per defect or affordance, at the smallest layer
  that owns it. A repair that rescues exactly one case is the wrong
  repair.
- Never repair while a session is live; a defect exists when a stream
  shows it.
- A deletion is its own commit, so it reverts cleanly.
- A premise stated in a plan is verified against the tree before the
  commit that depends on it; a corrected premise is written into the
  commit message as loudly as the change.

## 4. How a change is proven

- After every change: the fast suite with its exit code checked, then
  `ruff check chemsmart tests`, `black`, `isort`, and the docs linters as
  fixed points.
- Before a window is issued, the witness bank runs green: small probes
  of connected paths, program to artifact to selector to operation to
  claim, each derived from a loss this laboratory actually paid for and
  run through the public tool surface over archived evidence with no
  provider and no engine. A composition discovered inside a
  twelve-hour chemistry run is a composition nobody tested; the seal
  records the bank's report digest. The bank names the tree it actually
  imported: it once read the code under test from one clone and
  reported the harness of another, which is the class of defect it
  exists to catch, in the instrument that catches it. A new repair adds
  its witness, and the witness is shown red on the tree before the
  repair and green after -- a witness that was never red witnesses
  nothing, and one that constructs its own intermediate state witnesses
  nothing either: drive the public entry point and let the host build
  the state. A repair is proven on the path it was built for. An
  independent adversarial audit of the round's own implementation earns
  its cost: the first one found a blocking defect, two disconnected
  wires and three readers disagreeing about one state; the second found
  two more blocking defects over a green suite and a green bank; and an
  independent *scientific* review by a different model, asked to judge
  the round against the charter's values rather than its commits, found
  five more -- including an authority bypass in which
  ``--initial-decision deny`` launched an engine, and a defect the
  previous audit's own repair had introduced hours earlier. Three
  reviewers, three kinds of eye, and the round's own author found none
  of the fifteen.
  Both audits reported one pattern, and it is now a thing to look for
  by name: **the mechanism is right where it is computed and
  unconnected where it is consumed.** A projection that removes
  something needs a reader for what it removed -- subtracting a
  verified refusal was correct and nothing read the difference, so the
  refusal bought `achieved`. Two projections cannot each subtract the
  other's facts; one must inherit explicitly. A value the host resolves
  and then discards leaves a boolean nothing can audit. And in four of
  five findings the test or witness that should have caught it built by
  hand the exact state production fails to produce, which is why a
  witness drives the public entry point and a test that asserts on
  `inspect.getsource` is deleted on sight. One test went further and
  *asserted the defect*, so read what a failing test was pinning before
  changing the code to satisfy it. Two more rules earned by the third
  review: **an explanation string never confers authority** -- a
  mapping that describes both admitted and refused things must not be
  the thing that grants them -- and **a record's existence is not a
  grant**; a gate that asks whether state exists is asking the wrong
  question when the invariant is whether a human decided. A fourth
  review, asked to treat two live defects as probes rather than patch
  requests, found a third blocking defect and a laundering channel
  cheaper than the one the round had just closed -- so **use a defect
  as a probe, not a ticket**: ask what class it belongs to and where
  else that class lives, because the two lines it points at are rarely
  the whole of it. Two more rules it earned: a name rebound inside a
  branch of a long loop corrupts everything read after it, and it is
  mechanically detectable; and **surviving an error and preserving what
  the error interrupted are two different properties** -- a handler
  that returns before the projection keeps the process and loses the
  science.
- A behavioural change to the Agent is followed by a sealed live
  observation on chemically different tasks. N runs are N observations;
  a rate lives inside one contiguous window; a weak run is never
  re-rolled.
- The seal is written and closed before issue, with input digests, the
  runtime configuration, the physics bands never scored, and the
  falsifiers armed. If it leaks, the case is void and reported void.
- Physics outranks the scoreboard: every delivery is read against
  geometry, arithmetic, and the constants registry, and
  achieved-per-contract and wrong-per-physics are both written down.
- Errors are reported as loudly as wins, in the first section. An
  inferred mechanism reported as fact is corrected in place.
- An insight is registered in the memory ledger the day it is read.
- A session-ending mechanism is named by replaying the session's public
  transcript on the commit that ran it, never inferred from the stream;
  the same replay on the repaired tree is the repair's first probe.
- Two host organs that answer one question call one function; a
  frontier that admits what a review refuses is a defect in whichever
  organ grew alone.
- A literature constant is registered only from a source text read in
  the session, and its provenance names what was read; a value recalled
  from memory is not registered, however plausible, and a value whose
  primary text could not be read is absent by name.
- When a test disagrees with the host, the physics is checked before
  the code: a planar triatomic is Cs and never C1, a triatomic with
  equal bonds is C2v about its bisector however it is tilted, and a
  uniform shift of every heavy atom is a translation that alignment
  removes. Two failing tests in one day were wrong tests.
- A planted gem is checked to survive the optimiser before issue, a
  planted false gem is checked to exist, and the task text never asks
  for the observation the change under test is meant to elicit (the
  first STANDING window: a planar cyclohexane relaxed to the
  twist-boat unplanted, MMFF benzenes gave no imaginary modes, and
  both arms reported saddle character because the text asked).

## 5. The invariants every session carries

The claim ladder -- proposed, planned, materialised, previewed,
approved, executing, engine-complete, parsed, scientifically validated,
interpreted -- is owned by the deterministic host alone. Provider text
is never execution evidence. The hub invariant holds: one public
YAML-and-CLI layer, no second approval plane, no LLM grading a step.

## 6. Repository hygiene

- One branch per round; push only on fresh explicit instruction.
- Campaign evidence, licensed media, credentials, private transcripts,
  and generated program inputs never enter Git. `experiments/` is
  untracked scratch and is never touched.
- `/opt/chemsmart` stays stable; work happens in the research clone with
  `PYTHONPATH` set, because the controller environment shadows it.
