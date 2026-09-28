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

## Status

- Phase I, step 1 (pre-registration) committed before any census row.
