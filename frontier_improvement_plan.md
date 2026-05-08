# chemsmart frontier improvement plan

_Date: 2026-05-08_
_Basis: `qa_frontier_research.md` frontier review + `fork/main` code inspection only. No implementation in this change._

## Scope and planning assumptions

This document proposes a **design-only** migration plan for three improvements:

1. hierarchical multi-agent planning
2. RAG-backed method advisory
3. self-debug / replan loop

The goal is to raise chemsmart toward the 2024-2025 frontier **without breaking its current strength**: typed, inspectable, HPC-aware execution with a conservative Plan-Approve-Execute gate.

Two constraints shape the plan:

- **Do not replace typed execution with free-form text generation.** Frontier evidence still does not support raw Gaussian/ORCA input generation as the default production path.
- **Preserve the current approval semantics.** Hierarchical planning and self-debug should strengthen, not bypass, dry-run + critic + runtime validation.

---

## Current architecture baseline (for change planning)

### Existing control flow in `fork/main`

The current agent runtime is centered in `chemsmart/agent/core.py`:

- `AgentSession.run()` initializes session state and performs the first planner call (**lines 165-206**).
- `AgentSession._continue_run()` executes all non-risky steps, calls the critic, applies deterministic gates, and then executes risky steps (**lines 208-348**).
- `AgentSession._planner_call()` performs a single planner LLM call and returns a flat `Plan` (**lines 355-375**).
- `AgentSession._critic_call()` performs a single critic LLM call on dry-run inputs (**lines 376-398**).
- `AgentSession._execute_step()` is the single tool execution choke point (**lines 400-456**).
- `AgentSession._apply_deterministic_gates()` and `_block_reason()` provide hard safety gating (**lines 666-766**).

The typed tool surface is defined in:

- `chemsmart/agent/registry.py` (`ToolRegistry.default()`, **lines 59-80**)
- `chemsmart/agent/tools.py`
  - `recommend_method()` at **lines 551-640**
  - execution tools such as `build_job`, `dry_run_input`, `validate_runtime`, `run_local`, `submit_hpc`

This means chemsmart already has the right substrate for frontier improvements: a narrow typed tool layer plus auditable session control.

---

# 1. Hierarchical Multi-Agent architecture transition plan

## 1.1 Target design

### Today

```
User request
  -> planner.md
  -> flat Plan
  -> execute safe tools
  -> critic.md
  -> execute risky tools
```

### Target

```
User request
  -> Root planner agent
      -> subtask graph + task decomposition
  -> Specialist agents (method / workflow / runtime)
      -> subplans + constraints + evidence
  -> Aggregator agent
      -> final flat executable Plan
  -> existing dry-run / critic / runtime gates
  -> execution
```

### Architectural principle

The hierarchical system should be **internal-only** at first. That means:

- user-facing CLI still exposes one `chemsmart agent run ...`
- final executable artifact is still the current flat `Plan`
- all risky actions still pass through the existing critic and deterministic gates

So the migration is **planner-internal refactoring**, not a rewrite of execution.

## 1.2 Proposed agent roles

### A. Root Planner Agent

Responsibilities:

- classify request complexity
- decide whether to stay in `linear` mode or use `hierarchical` mode
- produce a task DAG / ordered subtask list
- assign subtasks to specialist agents

Typical subtask splits:

- structure preparation
- method-selection subproblem
- execution workflow assembly
- runtime / submission readiness review

### B. Method Specialist Agent

Responsibilities:

- consume request + molecule metadata + retrieved evidence
- propose settings intent only, not raw input text
- emit normalized settings recommendation for `build_gaussian_settings` / `build_orca_settings`

### C. Workflow Specialist Agent

Responsibilities:

- map chemistry intent into tool sequence
- choose job kinds, mixed-program handoff pattern, and expected intermediate artifacts
- remain bound to registered tool names

### D. Runtime / HPC Specialist Agent

Responsibilities:

- reason about whether server validation, local run, or preview submit is required
- identify points where `validate_runtime`, `run_local`, or `submit_hpc` should appear
- annotate remote-unknown risks for critic visibility

### E. Aggregator Agent

Responsibilities:

- merge specialist outputs into one final `Plan`
- resolve conflicting proposals conservatively
- guarantee the final plan is linear, JSON-serializable, and valid against current tool schemas

## 1.3 How it integrates with the existing Plan-Approve-Execute pattern

This must be an **extension**, not a replacement.

### Keep unchanged in v1 transition

- `Plan` and `Step` models in `core.py` (**lines 53-62**)
- flat executable step list consumed by `_continue_run()`
- `dry_run_input` -> `validate_runtime` -> critic -> risky-tool gate
- final execution through `_execute_step()`

### Add before the current planner output is finalized

New internal stage:

1. root planner creates `HierarchicalPlanDraft`
2. specialist agents refine their assigned subproblems
3. aggregator converts the hierarchical draft into the current flat `Plan`
4. current `_continue_run()` consumes that flat plan unchanged

### Key design rule

**Only the planner side becomes hierarchical.** The approval and execution side remains flat until the planner stack is stable.

## 1.4 Concrete file and function changes

### A. `chemsmart/agent/core.py`

#### Change 1: replace direct `_planner_call()` usage in `run()`

Current integration point:

- `run()` lines **191-198** call `_planner_call(request)` and immediately persist the returned `Plan`.

Planned change:

- replace `_planner_call()` with `_build_execution_plan()`
- `_build_execution_plan()` decides `linear` vs `hierarchical`
- `linear` mode returns current planner output for rollback safety
- `hierarchical` mode orchestrates root planner -> specialists -> aggregator -> `Plan`

#### Change 2: split planner responsibilities out of `core.py`

Current state:

- `_planner_call()` and `_critic_call()` are embedded in `AgentSession`

Planned change:

- keep `_critic_call()` in `AgentSession` for now
- move planner orchestration into a new module (`orchestration.py` or `planner_runtime.py`)
- `AgentSession` becomes session lifecycle + gatekeeper, not planner brain

#### Change 3: persist hierarchical planning artifacts

Current result artifacts:

- only final step JSONs and `decision_log.jsonl`

Planned addition:

- persist `hierarchical_plan.json`
- persist `shared_canvas.json`
- persist `subagent_outputs/<role>.json`

This preserves explainability.

### B. `chemsmart/agent/registry.py`

#### Change 1: add capability metadata to `ToolSpec`

Current `ToolSpec` fields (**lines 23-29**) are minimal.

Planned addition:

- `capabilities: list[str]`
- `risk_level: Literal["safe", "review", "risky"]`
- `domain_tags: list[str]`

Why:

- root planner and specialists should route tasks using structured tool metadata instead of prompt-only heuristics

#### Change 2: add tool groups

Add helper methods:

- `list_tools_by_capability(capability)`
- `list_risky_tools()`
- `list_method_tools()`

This helps specialist prompts stay compact.

### C. `chemsmart/agent/cli.py`

Planned additive changes:

- new hidden or expert flag: `--planner-mode linear|hierarchical|auto`
- optional debug flag: `--dump-canvas`

Reason:

- staged rollout and rollback
- no default behavioral break on day 1

### D. `chemsmart/agent/prompts/planner.md`

Planned change:

- narrow it to `linear` planner mode only
- do **not** overload it with hierarchical role instructions

Instead create separate prompts for separate roles.

## 1.5 New files / modules to add

### Planning runtime

1. `chemsmart/agent/orchestration.py`
   - top-level coordinator for hierarchical planning
   - owns mode selection (`linear`, `hierarchical`, `auto`)
   - exposes `build_execution_plan(request, registry, provider, session_dir, ...)`

2. `chemsmart/agent/schemas.py`
   - shared Pydantic models for:
     - `SubTask`
     - `HierarchicalPlanDraft`
     - `AgentAssignment`
     - `SpecialistRecommendation`
     - `AggregatedExecutionPlan`
     - `RepairContext` (reused by self-debug loop)

3. `chemsmart/agent/memory.py`
   - shared planning canvas
   - stores decomposition, molecule summary, retrieved evidence summary, unresolved risks, and assumptions

### Specialist agent modules

4. `chemsmart/agent/agents/root_planner.py`
   - request decomposition
   - task graph construction

5. `chemsmart/agent/agents/method_specialist.py`
   - settings intent generation from RAG evidence + molecule metadata

6. `chemsmart/agent/agents/workflow_specialist.py`
   - tool sequence planning and intermediate artifact design

7. `chemsmart/agent/agents/runtime_specialist.py`
   - runtime/scheduler/HPC planning logic

8. `chemsmart/agent/agents/aggregator.py`
   - merges specialist outputs into current flat `Plan`

### Prompts

9. `chemsmart/agent/prompts/root_planner.md`
10. `chemsmart/agent/prompts/method_specialist.md`
11. `chemsmart/agent/prompts/workflow_specialist.md`
12. `chemsmart/agent/prompts/runtime_specialist.md`
13. `chemsmart/agent/prompts/aggregator.md`

### Tests

14. `tests/agent/test_hierarchical_planner.py`
15. `tests/agent/test_shared_canvas.py`
16. `tests/agent/test_aggregated_plan_schema.py`
17. `tests/agent/test_linear_mode_compat.py`

## 1.6 Detailed change order (dependency order)

1. Add shared planning schemas (`schemas.py`)
2. Add shared canvas store (`memory.py`)
3. Add tool capability metadata in `registry.py`
4. Add root planner + specialist modules
5. Add aggregator module that always produces current flat `Plan`
6. Add `orchestration.py`
7. Update `AgentSession.run()` to call orchestration entrypoint instead of `_planner_call()`
8. Add CLI feature flag and default to `linear`
9. Add tests for parity: hierarchical planner must produce same/compatible flat plan on simple requests
10. Only after parity is shown, move default from `linear` to `auto`

---

# 2. RAG-based Method Advisor addition plan

## 2.1 Goal

Upgrade `recommend_method()` from a rule-only YAML matcher to a **hybrid advisor**:

- first consult local authoritative project YAMLs
- then retrieve program-documentation and method-guidance evidence
- then return a structured recommendation with citations and confidence
- still emit typed fields for `build_gaussian_settings` / `build_orca_settings`

## 2.2 Strategy: hybrid, not full replacement

### Current function

- `recommend_method()` in `chemsmart/agent/tools.py` **lines 551-640**

### Planned behavior

Mode order:

1. **project-hit mode**: current YAML rule path remains highest priority
2. **RAG-backed fill mode**: when no exact project hit exists, retrieve evidence and produce a structured recommendation
3. **fallback-decline mode**: if evidence is weak or contradictory, return `match=None` plus retrieved evidence summary and force the planner to ask for user confirmation or use conservative defaults

This preserves safety while broadening coverage.

## 2.3 Recommended corpus sources

### Tier 1: local authoritative chemsmart context (must-have)

1. `~/.chemsmart/gaussian/*.yaml`
2. `~/.chemsmart/orca/*.yaml`
3. `~/.chemsmart/usersettings.yaml` only for context labels, not method ranking
4. chemsmart README and agent prompts as operational policy context

### Tier 2: program documentation (must-have, but license-aware)

5. **ORCA manual / ORCA docs** for method keywords, solvent options, frequency/IRC syntax, auxiliary basis, DLPNO settings
6. **Gaussian documentation** only through **user-supplied local docs or internally authored curated notes**
   - do **not** bundle or redistribute proprietary manual content in the repo
   - index user-local extracted notes under `~/.chemsmart/agent/corpus/gaussian/`

### Tier 3: curated chemistry-method guidance (recommended)

7. internally curated markdown notes summarizing widely used method heuristics, e.g.:
   - B3LYP / PBE0 / M06-2X / ωB97X-D use cases
   - def2-SVP / def2-TZVP / ma-def2 / diffuse basis heuristics
   - SMD vs CPCM/IEFPCM solvent selection notes
   - heavy-element / ECP guidance
   - open-shell caution notes
8. Basis Set Exchange metadata or internally curated basis lookup tables
9. optionally curated benchmark summaries for method choice (not raw papers as first-class retrieval units)

### Tier 4: optional future sources

10. local copies of lab/group SOPs
11. validated previous successful input files from the user’s own projects

## 2.4 Vector DB and indexing strategy

## Recommended DB

**Primary choice: ChromaDB**

Why:

- simple local persistence
- rich metadata filters
- easy file-backed deployment under `~/.chemsmart/agent/rag/`
- low operational burden compared with server DBs

## Recommended retrieval shape

Use **hybrid retrieval** rather than pure vector search:

- vector retrieval for semantic matching
- BM25 keyword retrieval for exact chemistry tokens (`SMD`, `CPCM`, `DLPNO-CCSD(T)`, `IRC`, `ma-def2-TZVP`, etc.)
- rerank merged candidates before prompting the advisor

## Chunking strategy

Chunk by **semantic unit**, not fixed token length:

- one YAML section (`gas`, `solv`, `td`) per chunk
- one ORCA/Gaussian manual subsection per chunk
- one curated heuristic note per chunk

Target chunk size:

- 300-800 tokens
- overlap only where section headers and option lists need continuity

## Metadata fields to store

For each chunk:

- `source_type`: `project_yaml`, `orca_doc`, `gaussian_local_doc`, `heuristic_note`, `historical_input`
- `program`: `gaussian`, `orca`, `general`
- `task_types`: `opt`, `ts`, `freq`, `sp`, `irc`, `scan`
- `charge_regime`: `neutral`, `anion`, `cation`, `unknown`
- `multiplicity_regime`: `closed_shell`, `open_shell`, `unknown`
- `heavy_element_support`: list or bool
- `solvent_keywords`
- `confidence_tier`
- `project_name`
- `citation`
- `license_class`: `redistributable`, `user_local_only`

This enables hard filtering before any LLM sees the context.

## 2.5 Embedding model choice

### Default

- local embedding model via `sentence-transformers`, e.g. `BAAI/bge-small-en-v1.5`

Why:

- avoids extra paid API dependency
- reproducible in offline-ish HPC-adjacent environments
- good enough for documentation retrieval

### Optional premium mode

- OpenAI embeddings only as opt-in backend

Do **not** make hosted embeddings the default because the rest of chemsmart is file-centric and user-local.

## 2.6 Codebase changes

### A. `chemsmart/agent/tools.py`

#### Existing function to refactor

- `recommend_method()` at **lines 551-640**

#### Planned split

Keep `recommend_method()` as the public tool name, but internally route to:

- `_recommend_method_rules(...)`
- `_recommend_method_rag(...)`
- `_recommend_method_hybrid(...)`

The public tool contract should remain stable for planner compatibility.

#### Planned output additions

Extend returned dict with:

- `citations: list[str]`
- `evidence_summary: str`
- `confidence: float`
- `advisory_mode: Literal["rules", "rag", "hybrid", "decline"]`

### B. `chemsmart/agent/registry.py`

If output schema text becomes long, add tool description notes that `recommend_method` is evidence-backed and may return citations.

### C. `chemsmart/agent/core.py`

No major control-flow changes required for first RAG phase, but:

- `render_plan()` should eventually display cited rationale compactly
- session metadata should persist method citations in decision logs for auditability

### D. New modules to add

1. `chemsmart/agent/rag/corpus.py`
   - corpus loaders for YAML, local docs, curated notes

2. `chemsmart/agent/rag/index.py`
   - build/update/load vector index

3. `chemsmart/agent/rag/retriever.py`
   - hybrid retrieval and reranking

4. `chemsmart/agent/rag/advisor.py`
   - transforms retrieved evidence + task context into normalized recommendation output

5. `chemsmart/agent/rag/models.py`
   - Pydantic models for retrieved chunks and advisor outputs

6. `chemsmart/agent/rag/README.md`
   - explains corpus licensing and local-doc setup

### E. New tests

- `tests/agent/test_method_advisor_rules_mode.py`
- `tests/agent/test_method_advisor_rag_mode.py`
- `tests/agent/test_method_advisor_decline_mode.py`
- `tests/agent/test_method_advisor_citations.py`
- `tests/agent/test_rag_index_filters.py`

## 2.7 New dependencies needed

Recommended new dependencies:

- `chromadb`
- `sentence-transformers`
- `rank-bm25`
- `pypdf` or `pymupdf` (for local docs ingestion; choose one)
- `numpy` (likely already present transitively, but pin if needed)

Optional:

- `rapidfuzz` for exact-ish option normalization
- `whoosh` if BM25 implementation should be file-backed rather than in-memory

## 2.8 Deployment notes

- default to **rules mode** if no index exists
- expose `chemsmart agent rag build-index` later as a separate CLI command
- keep RAG opt-in until corpus quality is proven

---

# 3. Self-Debug Loop addition plan

## 3.1 Goal

Add a bounded reflexion-style repair loop that reacts to:

- critic `reject`
- critic `warn` that is non-remote and fixable
- tool execution errors during safe steps
- malformed dry-run inputs

The loop should **repair the plan**, not blindly retry the same plan.

## 3.2 Target behavior

### Repair loop states

1. initial plan generated
2. safe steps executed
3. critic + deterministic gates evaluate plan
4. if fixable failure detected, create `RepairContext`
5. repair planner produces revised plan diff
6. rerun safe steps from the earliest invalidated step
7. re-critic
8. stop on success or bounded termination

## 3.3 Integration points in `core.py` (current line anchors)

### Primary integration point A: blocked return path in `_continue_run()`

Current location:

- `_continue_run()` **lines 276-295** determine block reason and return early when blocked.

Planned insertion:

- before returning blocked result, call `_attempt_repair(...)`
- `_attempt_repair(...)` decides whether the failure class is repairable
- if repair succeeds, replace `self.state.plan`, reset `current_step_index` to the rewind point, and re-enter `_continue_run()`

### Primary integration point B: exception path in `_continue_run()`

Current location:

- exception branch at **lines 338-348**

Planned insertion:

- intercept tool/runtime exceptions for non-risky steps before final raise
- convert exception into `RepairContext`
- only escalate immediately for risky-step failures or repeated identical errors

### Primary integration point C: planner/critic helpers

Current locations:

- `_planner_call()` **lines 355-375**
- `_critic_call()` **lines 376-398**

Planned addition:

- add `_repair_planner_call(repair_context)`
- optionally add `_repair_critic_call(...)` only later; first version should reuse existing critic

### Primary integration point D: deterministic gate output

Current location:

- `_apply_deterministic_gates()` **lines 666-733**

Planned addition:

- classify issues as:
  - `repairable`
  - `user_confirmation_required`
  - `hard_stop`

Examples:

- repairable: malformed route, geometry handoff missing, missing IRC keyword, missing safe-step prerequisite
- hard stop: duplicate submission, server local validation fail due missing executable, unsupported capability

## 3.4 New modules to add

1. `chemsmart/agent/repair.py`
   - repair loop controller
   - issue hashing
   - plan diffing
   - rewind-point calculation

2. `chemsmart/agent/prompts/repair_planner.md`
   - prompt specialized for revising an existing plan under constraints

3. `chemsmart/agent/prompts/repair_critic.md` (optional phase 2)
   - only if separate post-repair validation prompt is needed

4. `tests/agent/test_repair_loop_critic_reject.py`
5. `tests/agent/test_repair_loop_tool_error.py`
6. `tests/agent/test_repair_loop_repeat_issue_stop.py`
7. `tests/agent/test_repair_loop_risky_step_no_retry.py`

## 3.5 Maximum iteration count and termination rules

### Recommended defaults

- `max_replans_per_session = 2`
- `max_total_safe_repair_cycles = 3`
- `max_same_issue_repetitions = 1`

### Termination conditions

Stop repairing if any of the following holds:

1. revised plan hash == previous plan hash
2. normalized issue signature repeats
3. rewound plan still fails at the same step with same error class
4. repair budget exceeded
5. repair would require bypassing deterministic hard gate
6. risky step already executed and failed in a side-effectful way

## 3.6 Infinite-loop prevention mechanisms

1. **plan hash tracking**
   - hash `plan.model_dump()` after normalization
2. **issue signature tracking**
   - hash `(verdict, sorted(issues), failing_step_tool)`
3. **rewind boundary**
   - only re-execute from earliest invalidated step, not whole session unless required
4. **time budget**
   - e.g. stop if repair planning exceeds 2x initial planning wall time
5. **no autonomous override of user gates**
   - repair loop cannot flip `allow_remote_unknown` or `allow_critic_override`
6. **no auto-retry after duplicate submit rejection**
   - duplicate detection remains terminal in the same session

## 3.7 State model changes

Extend `SessionState` in `core.py` (**lines 72-90**) with:

- `repair_iteration: int = 0`
- `repair_history: list[dict[str, Any]] = Field(default_factory=list)`
- `planner_mode: str = "linear"`
- `active_plan_hash: str | None = None`

Persist this in `session.json` so resume works correctly.

## 3.8 Decision-log additions

Add new event kinds:

- `repair_attempt`
- `repair_plan`
- `repair_rejected`
- `repair_budget_exhausted`
- `repair_success`

This is essential for auditability because self-debug loops are otherwise opaque.

---

# 4. Implementation roadmap

## 4.1 Dependency graph

```text
Shared schemas / state extensions
        |
        +--> Hierarchical planner modules
        |
        +--> Self-debug repair context

Tool capability metadata ----> Hierarchical planner routing

RAG corpus/index/retriever ---> Method specialist agent
                           \-> Hybrid recommend_method()

Hierarchical planner --------> final flat Plan generation
Self-debug loop -------------> revised flat Plan regeneration

All of the above -----------> existing critic + deterministic gates
```

## 4.2 Recommended implementation order

### Phase 1 — RAG-backed Method Advisor

**Why first:**

- smallest blast radius
- improves method selection immediately
- gives useful infrastructure to the future method specialist agent
- preserves the existing linear planner

### Phase 2 — Self-Debug Loop

**Why second:**

- uses the existing single-agent planner
- directly improves robustness of current workflows
- forces clean issue typing and state persistence before multi-agent complexity is added

### Phase 3 — Hierarchical Multi-Agent Planner

**Why third:**

- highest architectural complexity
- benefits strongly from already having:
  - evidence-backed method advisor
  - repair infrastructure
  - richer tool metadata

## 4.3 Complexity estimates

| Improvement | Complexity | Why |
|---|---|---|
| RAG-based Method Advisor | **Medium** | new retrieval stack, but limited execution-path disruption |
| Self-Debug Loop | **Medium-High** | touches session state, control flow, retries, and resumability |
| Hierarchical Multi-Agent Planner | **High** | new planner runtime, new prompts, shared canvas, mode flags, and compatibility constraints |

## 4.4 Expected impact on current QA score (current baseline: 10/10)

Because the current baseline is already 10/10, the relevant metric is **stability risk vs frontier capability gain**.

| Improvement | Short-term QA impact during rollout | Steady-state QA impact after stabilization | Frontier capability gain |
|---|---|---|---|
| RAG-based Method Advisor | slight risk: **10 -> 9.5** if evidence quality is noisy | back to **10/10** with better rationale traceability | High on method correctness / user trust |
| Self-Debug Loop | moderate risk: **10 -> 9.0** if loop control is weak | back to **10/10** once bounded and well-tested | High on robustness |
| Hierarchical Multi-Agent Planner | highest risk: **10 -> 8.5/9.0** if enabled too early | back to **10/10** only after parity + feature flags + evals | Highest on autonomy / workflow breadth |

Interpretation:

- **RAG** is the safest frontier-oriented first move.
- **Self-debug** is valuable but must be bounded.
- **Hierarchical planning** should be rolled out last and behind flags, because it is the easiest way to lose the very reliability that currently earns chemsmart its strong QA profile.

---

# 5. Risk assessment

## 5.1 Hierarchical Multi-Agent transition risks

### Potential impact on existing functionality

- planner may emit more complex but less stable plans
- step references may become invalid after aggregation
- prompt count and token cost will increase
- explainability could get worse if subagent outputs are not persisted

### Rollback strategy

- keep `linear` planner path intact
- add `planner_mode=linear|hierarchical|auto`
- allow immediate fallback to current `_planner_call()` behavior

### Staged deployment

1. ship hidden `hierarchical` mode
2. run in shadow mode on saved sessions
3. compare aggregated flat plan vs current linear plan on benchmark prompts
4. only then enable `auto` for complex workflows

## 5.2 RAG Method Advisor risks

### Potential impact on existing functionality

- retrieval noise could lower recommendation quality
- stale or conflicting documents could overrule project YAMLs incorrectly
- licensing issues if proprietary docs are bundled improperly

### Rollback strategy

- preserve current rules path as authoritative default
- make RAG advisory opt-in or `auto` only when exact project rule misses
- allow env flag `CHEMSMART_METHOD_ADVISOR=rules`

### Staged deployment

1. ship corpus/index builder with no runtime use
2. enable RAG only for `match=None` cases
3. add citations to logs
4. only later allow hybrid override where confidence exceeds threshold

## 5.3 Self-Debug Loop risks

### Potential impact on existing functionality

- accidental infinite loops
- repeated safe-step execution may create confusing artifact churn
- repair plans may silently drift away from original user intent

### Rollback strategy

- feature flag: `CHEMSMART_SELF_DEBUG=off|critic_only|full`
- hard max repair counts
- preserve current blocked-return path as default

### Staged deployment

1. critic-reject-only repair mode
2. add safe-step tool-error repair mode
3. never auto-repair post-submit failures in first release
4. require explicit logging of every repair attempt

## 5.4 Cross-cutting safety risks

1. **Token/cost explosion** from multi-agent planning + retrieval
   - mitigate with compact evidence summaries and role-specific prompts
2. **Auditability loss**
   - mitigate with persisted canvas, retrieved citations, repair history
3. **Regressing deterministic safety gates**
   - mitigate by treating `_apply_deterministic_gates()` as non-bypassable
4. **Configuration sprawl**
   - mitigate with conservative defaults and three top-level feature flags only

---

# Recommended final rollout sequence

## Stage 0 — prep

- add shared schemas/state fields
- add feature flags
- add benchmark fixtures and regression prompts

## Stage 1 — method advisory

- ship hybrid `recommend_method()` behind rules-first policy
- no planner behavior change yet

## Stage 2 — bounded repair

- add repair loop for critic-reject before risky steps
- verify no-loop guarantees

## Stage 3 — hierarchical planner (opt-in)

- planner-side only
- aggregator must emit current flat `Plan`
- keep current critic and execution unchanged

## Stage 4 — hierarchical planner (auto mode)

- only after parity on simple workflows and better-than-linear behavior on complex workflows

---

# Summary recommendation

If only one principle is followed, it should be this:

> **Do not let frontier autonomy outrun chemsmart’s typed, auditable execution model.**

The right migration path is:

1. **RAG first** for better evidence-backed method advice
2. **bounded self-debug second** for robustness
3. **hierarchical planning last** behind flags and with plan flattening into the existing approval/execution path

That sequence maximizes frontier gain while protecting the current reliability profile.
