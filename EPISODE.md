# R11 episode c5 -- lens: composition

- **Base SHA (verified by `git rev-parse HEAD` at start):** `6fe89afda861b82cce32e4ff0221f1d836bdb1cd`
- **Branch:** `worktree-agent-acce5db1d31e6bd00`
- **Claim served:** C5 -- *general typed steps compose into multi-stage and
  cross-program routes without task-specific code.*
- **Radius:** census and audit instruments under `.agents/research/loop/`
  (new files only). No product code, rules or settings; no Agent session,
  provider arm or engine job. CUHK reads only as read-only slot jobs.

## The question as I currently understand it

Two halves, tested separately:

1. *Compose* -- do the records show approved, executed routes in which one
   calculation's typed product (a structure, a Hessian, a scan minimum, a
   reached geometry, a registered result) becomes another calculation's
   input, inside one approval or across cycles, within one program and
   across programs; and analysis chains that combine several calculations
   (thermodynamic cycles, side-by-side comparisons)?
2. *Without task code* -- does any of those routes rely on code that exists
   for one task (a pKa, a PCET scheme, a Fukui index, a CBS protocol), or is
   every step a general typed operation whose domain covers many tasks?

## Corrected premises (written before any counting)

- **P1 (brief prior "all small, <= 9 heavy atoms").** `release.json`'s own
  text contradicts it before any bundle is read: record 59/60 (Gaussian
  `ts` -> `irc` across the host geometry handoff on C17H10N2O2, 21 heavy
  atoms), 61/62 (PySCF `opt` -> three `td` on DANS, 20 heavy atoms), 67-70
  (xTB GFN2 minimum Hessian of C16H17N3O4S, 24 heavy atoms), 71 (ORCA
  `modred` -> `ts` on (R)-BINOL, 22 heavy atoms, across a cycle). To be
  confirmed or refuted from the bundles, not from the narration.
- **P2 (brief prior "xTB -> PySCF (2139548, 2139556, 2139560)").** These are
  scheduler-dispatch job ids inside ONE goal
  (`xtb-ir-acetamide-pyscf-stability-r10`, controller Slurm 2139448; they
  appear in its cycle-1/3/4 `dispatch.receipt.json`), so they are one route,
  not three.
- **P3 (brief prior "12 multi-stage routes, 2 side-by-side comparisons").**
  That is a count of `release.json` narration. The census counts a record
  as citing a multi-stage route only when the cited run's approval bundle
  shows the edge.
- **P4 (orientation, code).** `producer_edge_selection_rule`
  (`chemsmart/agent/execution.py`) admits four producer rules inside one
  approval: `validated_optimized_geometry` (any of gaussian/orca/pyscf/xtb
  to any), `validated_scan_minimum_geometry` (ORCA scans only),
  `validated_final_orca_ts_hessian` and `validated_producer_orca_hessian`
  (ORCA only). The reach limits are per program, not per task -- the
  distinction the census has to keep.

## Definitions (pre-registered)

Unit: an **approved cycle** = one approval bundle
(`.chemsmart-agent/replays/<id>/bundle.json`, any agent-state directory,
`resolution.decision == approve`), de-duplicated by real path. Goal ledgers
add cycle order, the run's `workflow_state` and the settlement word.

- *calc node*: a node of `approved_scientific_plan.nodes`.
- *in-approval data edge*: an `approved_scientific_plan.edges` entry with
  `edge_kind == data`; its rule is the `selection_rule` of the matching
  `frozen_workflow_approval.producer_edge_rules` entry (read, not inferred).
- *executed*: a node reached a terminal state in a host stream of that
  approval (`workflow_node_state_changed` records naming its approval id);
  an edge executed when its consumer reached `engine_complete` or later.
- *cross-cycle lift*: a node whose `molecular_identity.geometry_lineage`
  (or coordinate identity) names a result artifact of an earlier
  calculation (mode displacement, reached geometry, scan point, registered
  result) -- producer program read from the lineage record.
- Classes of an approved cycle (a cycle may carry several):
  - **S** one calc node, no data edge; **P** >= 2 calc nodes, no data edge;
  - **M** >= 1 in-approval data edge, all within one program;
  - **X** >= 1 in-approval data edge joining two programs;
  - **L** a cross-cycle lift (L-intra / L-cross by program);
  - **C** side-by-side comparison: calc nodes of >= 2 programs, no
    cross-program data edge, and one `quantity_expression` whose calc-node
    closure spans >= 2 programs;
  - **A** analysis composition: a `quantity_expression` whose calc-node
    closure holds >= 2 calc nodes.
- *Route shape*: the sorted set of `program:stage -> program:stage [rule]`
  edge triples plus lift kinds; S/P cycles have the shape of their node set.
- *Heavy atoms*: non-H atoms of `molecular_identity.atom_order`; a route's
  size is its largest node.
- *Operations used*: every `expression_nodes[*].operation` and
  `constant_name` of the cycle's toolchain.

Strata, counted and reported separately: ax41 mirror (`2026-09-14` slice,
all bundles), `experiments-public/`, CUHK R8, R9, R10 (goal paths from
`r10_jobs.tsv`; no `m*`, no `sealed*`), and a separate *pre-R8 CUHK*
stratum for the campaigns `release.json` cites (`pyscf-irc-20260920`,
`pyscf-ts-irc-20260921`, `xtb-ir-acetamide-pyscf-stability-r10`).

## Static audit criteria (pre-registered)

An operation or CLI whose name is a task is a **general conversion** when
(a) its inputs are typed by dimension, not by species or task (any free
energy, any temperature, any cardinal pair), (b) it owns only a definition
or convention (RT ln 10, the IUPAC sign and n F, h c / 2 pi k_B, a stated
extrapolation law with explicit parameters), (c) it holds no molecule- or
protocol-specific constant, and (d) its domain covers more than one task.
It is **task code** when its inputs or behaviour are tied to one task's
workflow (it takes the N species of one named cycle, encodes one
protocol's constants, or branches on the task). For task code the verdict
also records whether any Agent route reaches it (from the census:
`real_execution_argv` job words, expression operations) and what the
general form would be.

## Replay at the pin (pre-registered)

Sample: among **executed** cycles, the first (by sorted bundle path) of
every distinct route shape in classes M, X, L; the first of every distinct
(program set, task-named operation set) in classes C and A; plus the named
routes: CUHK `pyscf-irc-20260920` g2 (Slurm 2140679), CUHK
`xtb-ir-acetamide-pyscf-stability-r10`, R10 Q11 G1 (2151912) and G2
(2151911), ax41 `standing-round` `e6-pcet-1`, ax41
`qualification/interop-fukui-path`, ax41 `qualification/pka-agent-path`.
Cap 40 cycles; above the cap, shapes with a cross-program edge first, then
by cycle count.

Levels, each through the host's own functions at the pin (PYTHONPATH =
this worktree, HOME fenced, provider-free, engine-free):

- **R1 load** -- `load_workflow_execution_approval_bundle` on the archived
  bundle (the executor's own loader: plan, toolchain, frozen rules,
  digests).
- **R2 edges** -- `producer_edge_selection_rule` per data edge on the
  parsed plan equals the frozen rule; `admitted_producer_edge_rules`
  admits the edge set.
- **R3 analysis** -- the toolchain parses (every operation, arity and
  analysis kind still in the vocabulary).
- **R4 compile+preview** -- each approved node's archived
  `real_execution_argv`, placeholders bound to the archived artifacts
  (project YAML and geometry found by digest in the workspace; a consumer
  gets the archived handed-off geometry) and `--no-fake` -> `--fake`, run
  through the pin's live Click in a scratch directory: exit 0 and a native
  input written; compared to the archived native input where one exists.

Faithfulness gate: R1-R4 on the producing tree must reproduce the archived
record (R1 loads, R2 equals the frozen rules, R4 writes the archived native
input). Where the producing tree is unavailable or does not reproduce, the
pin result is reported **unanchored** and is not counted either way.

A cycle **replays at the pin** when R1-R3 pass and R4 compiles and
previews every approved node. A failure is classified at its first failing
level as: *drift* (schema/digest/option/writer change with the composition
intact), *refused native field* (the pin's loader refuses a project key as
native), *missing operation* (an operation, selector, rule or jobtype the
route used no longer exists), *task code* (the failure is caused by code
that exists for one task), or *infrastructure* (the archive lacks an
artifact the level needs).

## What would count (pre-registered, before any count)

- **Support**: executed M/X/L routes exist in >= 2 strata and >= 2 program
  pairs, every edge under a registered general rule, every analysis step
  in the general vocabulary; the audit justifies each task-named operation
  as a general conversion and finds no Agent route reaching task code;
  every sampled cycle that is anchored replays at the pin or fails only as
  drift / refused native field, with R2 passing on 100% of loadable cycles.
- **Narrowing**: crossings exist only for structures (no cross-program
  Hessian, scan point or lifted TS inside one approval), or only in one
  direction (cheap -> expensive), or only for small molecules, or R2/R3
  fail for a route class because a general operation or rule is missing
  at the pin; the claim is then restated to the reach the records show.
- **Contradiction**: an executed route's edge, analysis or compile step
  depends on task code (a task jobtype such as `pka` in an Agent argv, a
  task branch in the executor, a protocol constant inside an operation),
  or a route shape that can only be expressed through per-task code.
- **Widening**: >= 3 programs in one route, crossings in both directions,
  crossings on >= 20 heavy atoms, or thermodynamic cycles assembled from
  geometry-origin operations (derive / compose / append) plus general
  analysis.

## Replay sample, fixed before any replay ran

Computed by `composition_census.py --sample` over the de-duplicated census
(local pass + CUHK Slurm 2157197): 90 distinct keys among 316 executed
distinct cycles, 13 named-route cycles, cap 40 (named first, then
cross-program keys, then by cycle count). R4 (compile+preview through the
live CLI) is run only if budget remains after R1-R3 on all 40; if it is
not run, the replay verdict is stated for R1-R3 only.

```
ax41:qualification/interop-fukui-path/workspace2  replay-875de3e8de824e8b  named
ax41:qualification/pka-agent-path/workspace  replay-40849b5c766d48a3  named
ax41:qualification/pka-agent-path/workspace  replay-4af34edeb0d14510  named
ax41:qualification/pka-agent-path/workspace  replay-af0b0423658c4978  named
ax41:standing-round/workspaces/e6-pcet-1  goal-goal-e6-pcet-1-cycle-1  named
cuhk:pyscf-irc-20260920/goals/g2-hono  cycle-1, cycle-2  named
cuhk:r10/q11/goals/g1  cycle-1, cycle-2, cycle-3  named
cuhk:r10/q11/goals/g2  cycle-1  named
cuhk:xtb-ir-acetamide-pyscf-stability-r10  cycle-1, cycle-3  named
ax41:general-round h1-a0-r1 c1 (analysis x28), c4-r3 c1 (shape x13), h2 c1 (analysis x4)
cuhk:r10/q12 g2-nh3 c1 (analysis x2), r10/q15 g1 c2 (analysis x2)
ax41:pyscf-round-2 E4-formic-acid c1, pyscf-round g4-methanol c1, standing-round e6-pcet-2 c1 (shape x1)
cuhk:r10/q2 g1-hono c2, r10/q3 g2 c3, r10/q4 g1 c1, r8/integration smoke c1 + c2 (shape x1)
ax41:general-round c1-r1 c1 (analysis x122), pyscf-round-2 E2-acetone c1 (analysis x33)
ax41:general-round c5-r1 c1 (shape x14), s6 c1 (shape x10), c7-r1 c1 (analysis x9)
ax41:novel-round-3 po2-fluoro-sulfone-gauche c2 (shape x9); cuhk:r10/q28 g2 c2 (analysis x8)
ax41:pyscf-round g1-methylamine c1 (shape x6); general-round s4 c1, w1 c1 (analysis x5)
ax41:novel-round-3 ino2-dinickel-exchange c3 (shape x4); cuhk:r10/q28 g2 c1 (shape x4)
ax41:general-round h1b c1 (shape x3); novel-round-7 ino3-r14b c1 (shape x3)
```

## Position memo (Phase I)

### The question as I now understand it

C5 has three separable parts, each with its own instrument: (a) do
approved, executed routes compose one calculation's typed product into
another's input, within and across programs (census); (b) is every step
of those routes general -- no task branch, no task jobtype, no protocol
constant hidden in an operation (audit); (c) does the composition still
stand at the pin (replay). The brief's candidate narrowing "crossing where
a stage is exclusive" is not what the records show: crossings are
level-ladder (xTB -> DFT, ORCA saddle -> PySCF path) and side-by-side, and
occur where both programs have the stage.

### Evidence (denominators after de-duplication)

Census (`composition_census.py`; local pass + CUHK Slurm 2157197): 463
bundle paths, **458 distinct approved cycles** (ax41 304, public 5, CUHK
R8 21, R9 20, R10 87 from the 63 goal dirs `r10_jobs.tsv` names, pre-R8
21), **316 executed** (>= 1 node validated in a host stream), 240 goals.

- **Multi-stage**: 150 executed cycles carry an in-approval data edge or a
  lift from an earlier result (M 77, X 28, L-intra 49, L-cross 7, overlap
  allowed); 18 of them have a node of >= 20 heavy atoms, all within one
  program (Gaussian ts -> irc on C17H10N2O2; PySCF opt -> td on DANS; ORCA
  modred -> ts on BINOL; xTB on C16H17N3O4S; ORCA ts -> sp on 20-heavy
  saddles). P1 (brief prior) is **refuted for multi-stage routes**.
- **Cross-program**: 34 executed cycles (X or L-cross) in 4 strata (ax41,
  pre-R8, R8, R10). Ordered pairs: orca -> pyscf (37 edges, 10 lifts),
  xtb -> orca (19 edges), pyscf -> orca (7 edges, 2 lifts), xtb -> pyscf
  (4 edges, 1 lift), gaussian -> pyscf (2 edges, 2 lifts), orca ->
  gaussian (2 edges). Both directions for orca <-> pyscf. Two executed
  cycles hold three programs (ax41 `c4-r3`, `interop-fukui-path`: xTB ->
  ORCA and PySCF). Size: 33 of 34 have <= 9 heavy atoms (median 3); one has
  20 (R10 Q4 g1 cycle 1: ORCA ts -> PySCF irc -> ORCA opt; the goal settled
  `unreachable_from_evidence`). P1 **holds for cross-program routes but
  one**.
- **What crosses**: every cross-program in-approval edge is
  `validated_optimized_geometry` (a structure). Hessians cross only ORCA ->
  ORCA (`validated_final_orca_ts_hessian`, 14 edges) and scan minima only
  ORCA -> ORCA (21); `validated_producer_orca_hessian` appears in no bundle.
  A saddle reaches another program's IRC as a lifted structure
  (`reached_geometry` 4, `mode_displacement` 2), never with its Hessian.
- **Admission is not execution on the scheduler target**: CUHK's
  in-approval edges were mostly admitted and not run inside the approval
  that admitted them (gaussian -> gaussian 26 admitted, 0 consumer
  validated in-approval); the same compositions ran across cycles through
  lifts -- 141 lifts in all (91 `reached_geometry`, 50 `mode_displacement`).
- **Analysis composition**: 179 executed cycles hold an expression over >= 2
  calculations; 75 read registered results of earlier workflows (A-reg);
  12 are side-by-side cross-program comparisons (C).
- **Operations on routes**: `gibbs_to_redox_potential` 64 nodes and
  `gibbs_to_pka` 22 (all ax41, 0 on CUHK); the four CBS operations and
  `transition_state_crossover_temperature` on 0 routes; 5 registered
  constants selected. Task job words in approved argv: **0** of all nodes.
  Node kind `aggregate`: 0 (every node is `program_call`).

Audit (`task_code_audit.py` at the pin):

- `gibbs_to_pka` -- general conversion: any free energy, any T > 0; owns
  RT ln 10 (gas constant + hartree->kcal only); its domain is any pK =
  dG/(RT ln 10), so its name is narrower than its domain; assumes the
  caller's dG already carries the solution standard state (built in
  thermochemistry + the constants registry, visible in the chain).
- `gibbs_to_redox_potential` -- general: any dG, n > 0; owns E = -dG/(nF)
  and the IUPAC sign; no constant inside (F is the unit system); the
  reference electrode stays a visible subtraction.
- `exponential_cbs_limit`, `scf_exponential_cbs_limit`,
  `scf_inverse_power_cbs_limit`, `correlation_inverse_power_cbs_limit` --
  general extrapolation laws; cardinals and exponents are explicit and
  required at node construction; the evaluator passes them explicitly, so
  the helper defaults (alpha 3.9, Helgaker p = 3) never reach a route.
- `transition_state_crossover_temperature` -- general spectroscopic
  conversion (h c |nu| / (2 pi k_B)) that owns single-imaginary-mode
  selection.
- `chemsmart/cli/pka.py` with `cli/{gaussian,orca}/pka.py` and
  `jobs/{gaussian,orca}/pka.py` (3,810 lines, plus pKa code in `io/file.py`,
  `io/*/output.py`, `utils/datasets.py`, `jobs/*/settings.py`; 29 test
  files) -- **task code**, in the human hub (upstream, Feb-Mar 2026): an
  8-file proton-exchange/direct pKa scheme with protocol defaults (qRRHO
  Grimme 100 cm^-1, 1 M). **Unreachable from the Agent**: the Agent registry
  declares no `pka` (nor dias/nci/resp/wbi/crest/qrc/traj/userjob) jobtype;
  an overlay admitting (cpu, pka) is refused "support overlay broadens
  declared jobtypes" while the control pairs gaussian:link and orca:neb are
  admitted; `run pka` is `UNKNOWN_PROGRAM`; 0 approved argv carry it. Its
  general form already exists and is qualified: derive -> opt+freq per
  species -> thermochemistry at a stated standard state -> constant ->
  `gibbs_to_pka`.
- Lexicon scan of `chemsmart/agent` + `chemsmart/analysis` code: three
  task-flavoured identifiers (`extrapolate_correlation_helgaker`, a named
  published law; the crossover operation and its constant); no task branch.

Replay at the pin (`replay_composition.py`; 40 pre-registered cycles, 23
local + 17 CUHK Slurm 2157214; pin package digest 5a645356... on both
hosts):

- **R2 composition, 40/40**: every approved plan parses and reproduces its
  archived `plan_sha256` (40/40); producer rules recomputed at the pin equal
  the frozen rules on all 80 frozen data edges (85 edges checked);
  `admitted_producer_edge_rules` admits all 24 frozen edge sets (16 cycles
  have none: lifts or analysis-only); 0 operations missing from the
  vocabulary; 0 of 191 nodes lost its (engine, stage) execution pair.
- **R1 whole bundle, 2/40 load** (both late R10: Q2 g1-hono c2, Q28 g2 c1).
  All 38 refusals fall in the toolchain parse (R1 and R3 give the same
  reason 38/38): "scientific toolchain plan digest mismatch" 25, "claim
  rendering node declares output(s) that no input of the node carries" 12,
  "thermochemistry node requests entropy_method 'grimme' but declares only
  harmonic outputs" 1 -- a later analysis-chain admission contract.
  Classified **drift**: no refused native field, no missing operation, no
  task code.
- **Faithfulness**: 18 of 18 bundles run on their producing trees load and
  parse (local e8909401 x2, 006a65bf x1; CUHK Slurm 2157214, 15 bundles
  on 8 producing code directories found from each goal's own job output),
  so 16 of the 38 pin refusals are anchored as the tree's drift; 22 are
  unanchored (20 ax41 bundles not run on a producing tree, 2 acetamide
  cycles whose job output names no code). R2's pin result is anchored by
  the archived digest itself (40/40 `plan_sha256` reproduced).
- **R4 (compile + preview) not run** -- budget. "Still compiles and
  previews at the pin" is therefore not established.

### What it says about C5

**Narrowed** (with one widening). Supported: general typed steps compose
into executed multi-stage routes (150 cycles, 18 at >= 20 heavy atoms) and
cross-program routes (34 cycles, 6 ordered program pairs, 3 programs in one
approval twice); every edge rides one of three general producer rules;
every analysis step is in the general vocabulary; no Agent route reaches
task code; each task-named operation is a general conversion. Narrowed to:
*a crossing between programs carries a structure, not a Hessian or a scan
point, and is realised small (33/34 <= 9 heavy atoms); on the scheduler
target composition is realised across approvals by lifts rather than inside
one; and at the pin the composition layer replays while the archived
approved bundles do not load, because the analysis-chain contract moved.*
Widened: routes compose across workflows through registered results (75
executed) and through geometry-origin operations (derive, append, compose)
into thermodynamic cycles. Product-wide, "without task-specific code" is
false in the human hub (the pKa stack), which is duplicated authority
beside the Agent's composed pKa route.

### Proposed next program

1. R4 on the same 40 cycles (live CLI `--fake`, pin vs producing tree),
   locally for ax41 and in one CUHK slot job -- the one level that decides
   "still compiles and previews".
2. No implementation for C5 itself is warranted. Two reduction candidates
   for the master/owner, not implemented: re-express the human pKa CLI stack
   as a composed plan (one authority for one question); delete the stale
   `AGGREGATE_OPERATIONS` tuple in `chemsmart/analysis/aggregation.py`
   (no reader; its comment contradicts `workflows.AGGREGATE_NODE_STAGE`).

### What would change my mind

An R4 failure at the pin whose cause is a composed input (a handed-off or
lifted geometry a writer refuses) would contradict "still composes"; an
executed route whose argv carries a task job word, or an operation found to
hold a protocol constant, would contradict "without task code"; an admitted
cross-program Hessian or scan-point edge, or a cross-program route above 20
heavy atoms that validated, would widen the narrowed claim.

## Status

- Phase I, step 1 (pre-registration) committed before any census row.
- Census run: local (ax41 + public) and CUHK (Slurm 2157183, 2157188,
  2157197; read-only slot jobs). The master's duplicate warning is applied:
  rows are counted once per bundle content and goals once per goal digest.
