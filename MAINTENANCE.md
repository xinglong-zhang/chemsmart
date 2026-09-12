# Maintenance notes: what has worked here, and under what conditions

**This file is evidence, not authority.** `AGENTS.md` is the product
charter and `CONDUCT.md` is the working discipline; where either disagrees
with this file, they win and this file is out of date. Nothing here is a
specification. Each section says what was tried, what happened, and what
the conditions were, so that a later session can judge whether its own
situation resembles ours — and depart from us when it does not.

The failure mode this file exists to avoid is the opposite of
under-documentation: a frozen architecture, defended by a document, that
prevents the science the product exists to do. If a rule below blocks a
legitimate scientific action, the rule is the thing to change.

---

## Adding a quantity the analysis plane can serve

**What we did.** Added `reached_positions` — the structure an optimisation
stopped on, as distinct from the geometry its thermochemistry block
describes.

**What it took, in the order the failures arrived.** An accessor on the
reader was not enough. It also needed: an entry in `SELECTOR_UNITS`; an
entry in `SUPPORTED_SELECTORS` (the request gate, without which nothing
can ask for it); an entry in `_SELECTOR_DIMENSIONS`; and a per-jobtype
declaration, because a declaration is a semantic claim about what the
value means *for the job that ran*. Four separate registries, each of
which failed loudly and separately, one test run at a time.

**What worked well.** The failures were loud and immediate, and each named
the registry it wanted. Threading a new quantity took four iterations of
about a minute each. We would not trade that for a single permissive path.

**What we would watch.** We declared the new selector for `opt` and `ts`
only, and deliberately not for `irc` (whose log prints one structure) or
`scan` (whose last printed structure is a scan point, a different thing).
Choosing the jobtypes is the scientific act; if a future program's log
carries a reached structure for a jobtype we excluded, that exclusion is
the thing to revisit, not the mechanism.

## Adding a numeric policy the host decides science through

**What we did.** Added a 13th capability kind, `policy`, after finding
seven numeric thresholds living as bare module constants with no registry
and no rung — one of which put the bond/no-bond line at 1.120 Å for C–H
and reported a converged formaldehyde as having no C–H bonds.

**Why the existing kinds did not fit.** The ladder's only numeric kind was
`constant`, and all thirteen of its entries are literature values read
from a source text with provenance. A threshold is a decision the host
*makes*, not a value it *reads*. That distinction turned out to be the
whole reason the ladder could stay green over a false delivered vector: it
had no kind that could hold the thing.

**What worked well.** Declaring the policy with its value, unit,
functional form, owning module and consumers made the gap measurable
rather than invisible: one policy is `tested`, six are honestly
`advertised`, and two human-CLI conventions are marked `legacy`. The
`legacy` flag mattered more than expected — it lets a differential oracle
treat a *named* difference between conventions as a scientific
observation rather than a defect to eliminate.

**Condition worth stating.** We only declared policies we could name an
owner for. A policy with no single owning module would have needed the
owner first; the declaration is not a substitute for one.

## Repairs that create an owner

**What we did, twice, with different outcomes.**

In September a repair found a bad X–H tolerance, created
`_molecule_graph` as a single wide-tolerance perceiver, routed two
consumers through it, and its commit message said hydrogens now get the
ordinary buffer "in the one graph every host perception shares". Three
other callers kept the narrow tolerance. Nine days later one of them
delivered a wrong bond-pair vector into a scientific claim.

In the round after, a repair created
`WorkflowExecutionApprovalBundleV1.node_observation_lines()` as the single
reducer of a `tuple[dict, ...]` field — **and in the same commit added an
AST lint forbidding any other consumer from reducing that shape**, with a
planted offender that drives the production scanner.

**What worked.** The second pattern. When a repair creates an owner, the
same commit forbids the bypass mechanically. A comment does not hold, and
a commit message asserting a general property ("every", "the one",
"always") without a check is a claim that reads as done to every later
reader.

**Empirical note on the lint.** Our first version matched by field name
alone and over-captured: a field that is a tuple on one class and a
genuine `dict` on another was flagged. Narrowing it to *names the package
never uses as a mapping* made it correct and kept it general. This is the
same over-capture a marker regex without a left boundary made in the same
round, which suggests the general lesson is to scope a lint by what the
code *does* rather than what a name *looks like*.

## Adding a refusal

**What worked well.** Three properties, each earned by a loss:

- The refusal names the invariant it protects **and a legal route onward**.
  Measured over one round: the agent met 470 refusals, retried the refused
  tool in 60% of them, recovered within the session in 82%, and opened a
  guide in 2%. The refusal message, not the guide, is where the model is
  taught.
- The refusal carries the numbers behind it. When a perception convention
  blocks a geometry edit, the message now states the distance, the cutoff,
  the signed margin, the policy id and the routes out — because
  hypervalent, agostic and proton-transfer cases sit near such a line by
  their nature, and a bare "not bonded" told a session its structure was
  impossible when the host meant its threshold was close.
- The refusal is a **typed** error. A check that accepts any exception
  cannot tell a designed refusal from a defect, and a defect inside a
  refusal path is what ends goals. Assertions on refusals name the host's
  error class.

**What did not work.** Refusals that graded chemistry. Three checks on
composed uncertainty magnitudes fired 13 times with 0 legitimate
derivations prevented, and each was walked past by an equivalent
spelling. They became observations instead, and the delivery improved.

## Writing an oracle for a host-derived value

**Three kinds have been used here, and they catch different things.**

1. **Self-consistency** — declaration against implementation, name against
   state, parameter against round trip. Cheap, and it found real defects.
   **Its limit, measured:** `connectivity` satisfied every one of them with
   every declaration true and still reported that formaldehyde has no C–H
   bonds. A self-consistency oracle cannot see a wrong *value*.
2. **Metamorphic** — a relation that must hold over a transformation of
   the input. Worked where the host can synthesise the input (a bond
   graph over coordinates: perturb the distance, sweep the boundary). Did
   **not** transfer to a closed-form transformation over a program's
   printed block (per-atom mode participation), where the invariants are
   identities of the formula and the real risk — each program's mass
   table — needs the programs to test.
3. **Referential / cross-representation** — two host answers to one
   question must not disagree in silence. This is the one that caught the
   perception defect, and it needed **no chemistry at all**: the generator
   enumerates element pairs from the radii table, takes each declared
   policy's decision boundary, and probes below, at, above and *between*
   boundaries. It was red with 87 generated disagreements and green after.

**Condition that decides which applies.** A transformation that carries a
**threshold** has boundaries to generate from. One that does not (a
normalisation, an index remap) does not, and asking for a generated
adversarial domain there produces noise. The `policy` kind happens to
enumerate exactly the thresholded transformations, which is a convenient
accident worth preserving.

## Receipts, digests and schema growth

**What bit us.** We added a field to a thermochemistry receipt dataclass
and put it in the canonical digest body. Every receipt already written to
disk then failed revalidation, because the `record` an event carries **is**
the digest body and the validator requires the two to hash alike. Three
tests caught it; the run streams this laboratory has produced would have
stopped being evidence.

**What worked instead.** Carrying the new information in a field that is
*already* inside the body — `assumptions`, the free-text channel every
other thermochemistry control narrates itself in. No selection adds no
line, so historical digests are untouched; a real selection changes the
digest, which is what the field was for.

**What also worked.** When a receipt field's permitted contents had to
grow (adding a policy id and per-pair margins to a delivered adjacency),
replacing an exact-set check with an **allow-list** kept the invariant it
actually protected — no perceived *label* rides there — while admitting
measurement provenance. Read what a check protects before widening it.

## Removing a capability

**What we did.** Withdrew a distance-derived `bond_order` from everything
the agent can reach, after measuring that it reads ethane's single C–C as
2.0 and benzene's as 3.0 at the buffer a previous repair had chosen.

**What worked.** Withdrawing rather than recalibrating, because the
quantity was on the wrong side of the boundary: where electrons are is not
a distance. The human-CLI consumers that need an rdkit molecule kept
working through rdkit's own perception, and a test now holds that no
agent-reachable module reads a distance-derived order.

**What to be careful of.** The observable-regression guard refuses a
replan that *removes* a stage, on purpose — deleting the node that carries
a finding is the cheapest way to clear it. When we tried to distinguish
two plans by changing a node id, that guard fired correctly and we changed
the discriminator instead of the guard.

## Changing a test that disagrees with a change

**What worked.** Reading what the test pinned before satisfying it. Two
cases in this round:

- A test asserted all-zeros connectivity on a fixture whose O–H is
  1.200 Å. It held only because of the defect, so the suite had encoded
  the defect as a requirement. We re-derived the fixture from chemistry
  (a 1.200 Å O–H is a proton in transit) and added a second test covering
  the case the first had asserted away.
- A test asserted `freq` declares what `opt` declares minus `converged`.
  Adding a reached-geometry selector to `opt` broke it — correctly, because
  a fixed-geometry frequency job has no second structure. The set
  difference gained a second member with the reason written down.

**What we avoid.** Changing an assertion to match new behaviour without
establishing which of the two is right. One test in an earlier round
*asserted the defect*, which is why this is a standing habit.

## Before a long live window

**What worked, measurably.** Three unrelated four-atom cases, one question
each, run before a twelve-hour window. They cost about 23 minutes and 100
seconds of engine time and found two host defects that 3,183 tests, 67
witnesses and five lint gates had passed over — one of which ended a goal
unsettled. The previous long window took twelve hours to reach its first
defect.

**Conditions that made them useful.** Genuinely different chemistry;
one question each; deliverables that are single selectors or a
subtraction; and a deliberately wrong starting geometry so the relaxation
has work to do. A case that cannot fail teaches nothing.

**What we also learned to do.** Freeze the tree for the window's whole
life and arm that as a falsifier. Editing the clone while a goal is parked
made one window's later cycles a mixed-tree observation, because each
per-node subprocess imports the clone fresh.

## Delegating a review

**What worked.** Two reviewers, kept independent until both returned, each
given the same verified evidence and the same fundamentals verbatim, with
different primary questions — one forensic ("why did prevention fail"),
one architectural ("how should the pattern evolve"). They disagreed on a
central point, which was the most useful thing they produced: it located
the real decision instead of confirming a plan.

**What we insist on afterwards.** Every load-bearing claim verified
against source, behaviour and run record before it enters a synthesis. In
this round both reviewers corrected me and one of their corrections was
itself slightly wrong (a commit id that was a pre-rebase duplicate). Two
of their findings were about defects in work committed hours earlier.
