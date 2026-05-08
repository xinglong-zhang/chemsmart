# chemsmart agent AI systems audit

_Audited target:_ `fork/main` at commit `ee3bbd47` on 2026-05-08

## Executive summary

This is a **promising but not yet production-hardened** agent design from an AI-systems perspective.

The good news: the architecture keeps LLM context bounded, separates planner and critic concerns reasonably well, and writes enough execution artifacts to support postmortems. The bad news: reliability still leans too heavily on prompt obedience. The planner is fed a large, redundant schema blob; tool schemas are permissive and partially untyped; JSON is parsed from raw text rather than provider-native structured output; and there are no planner/critic retries, timeouts, or usage telemetry.

**Overall AI engineering score: 3/5.**

## Scoring table

| Dimension | Score (1-5) | Notes |
|---|---:|---|
| Prompt engineering | 3 | Clear role split and explicit JSON targets, but prompts are long, duplicated, and already drifting from tool reality. |
| Token & context efficiency | 2 | Planner call has heavy fixed overhead; full tool schema is injected every run. |
| Tool schema design | 2 | Enums help, but many fields degrade to `Any`, extra args are allowed, and most tool descriptions are empty. |
| Session & memory design | 4 | Resume is bounded and practical; logs are useful. Missing provider/model/schema versioning for reproducibility. |
| Error handling & reliability | 2 | Raw-text JSON parsing, no structured outputs, no retries/backoff/timeouts, no failure snapshot of full LLM payload. |
| Agent pattern assessment | 3 | Plan-first is the right family for risky tools; two-call planner+critic is only partially cost-justified in the current implementation. |

## Critical flaws

1. **Prompt/schema drift is already present.** `planner.md` tells the model `submit_hpc` only needs `job` (`chemsmart/agent/prompts/planner.md:63`), but the actual tool requires both `job` and `server` (`chemsmart/agent/tools.py:497-500`).
2. **Planner reliability depends on raw JSON text generation.** Both planner and critic return free-text JSON parsed by `_parse_json_response(...)` (`chemsmart/agent/core.py:355-398`, `859-888`) instead of provider-native structured outputs.
3. **Tool schemas are too permissive for dependable planning.** `ToolInputModel` uses `extra="allow"` (`chemsmart/agent/registry.py:19-20`), many parameters collapse to untyped fields, and most tool descriptions are just the function name.
4. **Token overhead is high before the model does any reasoning.** The sample planner call spends ~3.1k input tokens before counting the model's answer, mostly on the prompt and tool schema blob.
5. **Operational configuration is brittle.** Provider loading hardcodes `/Users/hongjiseung/developer/chemsmart/api.env` (`chemsmart/agent/providers.py:16`), which is a reliability smell for any shared or resumed agent workflow.

## Token cost estimate per run

Using the real session `~/.chemsmart/agent/sessions/20260507T143237Z-93e47ec1/decision_log.jsonl` and the current prompts/schemas on `fork/main`:

| Component | Approx. tokens |
|---|---:|
| Planner system prompt | 1,461 |
| Planner user payload (`request` + full tool defs) | 1,702 |
| Planner output (`Plan`) | 359 |
| Critic system prompt | 474 |
| Critic user payload (`plan` + one dry-run input) | 464 |
| Critic output (`CriticVerdict`) | 252 |
| **Total per sample run** | **~4,711** |

Notes:
- The planner input alone is **~3,163 tokens** before answer generation.
- The full tool schema blob is **11,597 characters / ~1.7k tokens**.
- Critic cost scales with full rendered input files, so larger molecules or multi-step workflows will push total cost into roughly the **5k-7k+ token** range.
- Context does **not** grow unbounded across turns; it is bounded per call. Storage does grow on disk because the decision log keeps full artifacts.

## Detailed assessment

### 1) Prompt engineering

**What is good**
- Planner and critic are separated cleanly in concept.
- Both prompts ask for explicit JSON with stable top-level keys.
- The planner prompt includes step-reference conventions and a decline path.
- The critic sees plan + generated inputs, not planner hidden reasoning.

**What is weak**
- The planner prompt is too long for the amount of real control it provides (`planner.md` is ~1.46k tokens by itself).
- Prompt content duplicates registry truth: job-kind enums, workflow templates, return-type guides, and argument conventions are repeated in text instead of derived from code.
- That duplication has already drifted: `submit_hpc` instructions no longer match the tool signature.
- The system depends on “JSON only” obedience instead of structured output APIs.

**Verdict:** conceptually sound, operationally brittle.

### 2) Token & context efficiency

**Strength**
- Session context is bounded. The agent does not replay the full transcript back into the model.
- Resume happens from disk artifacts, not from growing in-model memory.

**Weakness**
- `_planner_call()` serializes the entire tool registry into the user message every run (`chemsmart/agent/core.py:355-372`).
- The planner prompt plus schema payload dominates cost; the natural-language request is tiny by comparison.
- `build_orca_settings` is especially expensive because its schema is wide and mostly optional.
- Critic inputs include full rendered input text, so payload size scales with molecule size and number of dry runs.

**Verdict:** bounded, but not lean.

### 3) Tool schema design

**Strength**
- `build_job.kind` uses a real enum, which is exactly the sort of constraint LLM planners benefit from.
- Required fields are surfaced for several tools.

**Weakness**
- `ToolInputModel(... extra="allow")` means the planner is never strongly discouraged from adding junk keys.
- Several parameters become schema entries without useful type information (`functional`, `basis`, `job`, `settings`, `server`).
- Most tool descriptions are just names like `"build_job"` rather than compact semantic instructions.
- The schema generator silently falls back to `Any` when JSON schema generation fails (`registry.py:192-197`), which hides ambiguity instead of fixing it.

**Verdict:** good intent, weak machine-facing contracts.

### 4) Session & memory design

**Strength**
- `SessionState` stores `cwd`, current step index, plan, and minimal env snapshot; that is enough for practical resume.
- Decision logs are well structured and timestamped.
- Persisted step artifacts make replay/debugging easier than transcript-only approaches.
- `session_metadata.json` is useful for coarse monitoring: duration, verdict, request, input file, executed/planned steps.

**Weakness**
- Reproducibility metadata is missing: no provider name, resolved model, git SHA, prompt hash, schema hash, or token usage.
- Resume restores objects by import path and class name, but there is no version pinning. Replaying after code drift can change semantics.
- `session_metadata.json` is useful for dashboards, but not for cost analysis or exact postmortem reconstruction.

**Verdict:** solid first-pass session model; needs provenance fields.

### 5) Error handling & reliability

**Strength**
- Tool validation/execution errors are wrapped and logged before raising.
- Deterministic gates on malformed dry-run inputs, IRC keywords, duplicate submissions, and remote-unknown handling reduce blind trust in the critic.

**Weakness**
- There is **no** provider-native structured output, **no** planner/critic retry loop, **no** backoff, and **no** explicit timeout configuration in provider calls.
- `_parse_json_response()` raises on malformed JSON and only preserves the first 300 raw characters in the exception message; the full bad response is not logged for diagnosis.
- Planner/critic tests cover happy-path mocking, but not malformed LLM JSON, partial JSON, rate limits, or transient gateway failures.
- Planner and critic currently use the same provider abstraction and default model, so the critic is not strongly independent.

**Verdict:** too fragile for unattended production use.

### 6) Agent pattern assessment

**Is Plan→Approve→Execute the right pattern?**
- **Yes, mostly.** For expensive or side-effectful tools, a plan-first pattern is better than vanilla ReAct because it produces an auditable artifact before execution and supports deterministic preflight checks.
- However, the implementation is really **Plan → Critic → Execute/Gate**, not true human approval.

**Would ReAct/Reflexion be better?**
- **Pure ReAct:** worse for this use case. It would increase hidden state churn and reduce auditability around risky actions.
- **Reflexion/self-critique:** only marginally better than what exists unless paired with stronger structured outputs and an actually different critic model.
- **Best next step:** keep the plan-first skeleton, but make it more structured and cheaper.

**Is the two-LLM-call architecture cost-justified?**
- **Partially.** The sample decision log shows the critic catching a malformed Gaussian route pattern that the deterministic route-presence check would miss.
- **But** the marginal value is lower than it should be because most hard constraints are already deterministic and the critic is not independent enough to justify a full second premium call every time.
- A better cost profile would be: **planner on stronger model + critic on cheaper model**, or **planner + deterministic validators only**, escalating to critic only when heuristics detect risk.

## Top 5 AI engineering recommendations

1. **Move planner/critic to provider-native structured outputs and add bounded retry/backoff/timeouts.**  
   _Code:_ `chemsmart/agent/core.py:355-398`, `chemsmart/agent/core.py:859-888`, `chemsmart/agent/providers.py:53-66`, `99-111`.  
   Why: this is the single biggest reliability improvement. It removes most JSON-parse brittleness and gives you clean recovery paths for transient failures.

2. **Stop duplicating tool truth in `planner.md`; generate a compact planner registry from code.**  
   _Code:_ `chemsmart/agent/core.py:357-367`, `chemsmart/agent/registry.py:31-50`, `chemsmart/agent/prompts/planner.md:8-63`.  
   Why: today the prompt and actual schema have already diverged. A compact generated summary (tool name, one-line description, required args, enum values) would cut tokens and reduce drift.

3. **Tighten schemas aggressively: forbid extra fields, add real type annotations/docstrings, and remove `Any` fallbacks where possible.**  
   _Code:_ `chemsmart/agent/registry.py:19-20`, `126-197`; tool signatures in `chemsmart/agent/tools.py`, especially `submit_hpc` at `497-548`.  
   Why: LLM planners are much more reliable when the machine contract is narrow and semantically described.

4. **Fix the `submit_hpc` planning contract immediately.**  
   _Code:_ `chemsmart/agent/prompts/planner.md:15-16, 59-63` versus `chemsmart/agent/tools.py:497-500`.  
   Why: this is a direct planner/runtime mismatch. Either make `server` optional with a deterministic default lookup, or require it consistently in both prompt and schema.

5. **Add provenance and cost telemetry to every session.**  
   _Code:_ `chemsmart/agent/core.py:72-83`, `93-117`, `547-580`; `chemsmart/agent/providers.py:165-176`.  
   Why: log provider name, resolved model, token usage, prompt/schema hashes, git SHA, and full raw LLM responses on parse failure. This will turn `session_metadata.json` from “basic monitoring” into a production debug asset.

## Bottom line

From an AI/ML engineering standpoint, the design direction is right: bounded context, explicit planning, deterministic execution, and post-hoc artifacts. But the current implementation is still **prompt-fragile** rather than **contract-driven**. The next maturity jump should be to replace text-level obedience with stronger machine interfaces: structured outputs, narrower schemas, cheaper prompts, and richer telemetry.
