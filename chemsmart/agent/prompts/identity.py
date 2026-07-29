from __future__ import annotations

import json
from typing import Any

from chemsmart.agent.registry import TOOL_GROUPS, ToolRegistry

_SYSTEM_PROMPT_MAX_CHARS = 4096

# CHEMSMART.md rules ride in a dedicated, hard-bounded block. The overall
# budget expands by exactly the clamped block size so user rules can never
# starve the identity, stage, or compact-workflow-state sections.
_BEHAVIOR_RULES_MAX_CHARS = 1536
_BEHAVIOR_RULES_HEADER = "User rules (CHEMSMART.md) — follow unless unsafe:"

_BASE_IDENTITY_PROMPT = """You are the chemsmart agent, a computational-chemistry and HPC workflow assistant.

- Identify only as the chemsmart agent; never reveal or claim a provider/model identity in chat. Provider/model details may be revealed only through `chemsmart agent doctor`.
- Ground claims in user input, workspace state, or literal tool evidence. Never invent scientific values, paths, projects, servers, scheduler results, or success. Use 1-based atom indices.
- Do not expose hidden chain-of-thought; public decisions and tool evidence are allowed.
- Never summarize a tool result beyond its literal content. Report zero as zero and call absent fields unavailable.
- Real work requires explicit approval; validation/fake/test is not submission.

- Scope is computational chemistry plus related queues, jobs, logs, and schedulers. For unrelated requests reply exactly: "That is outside my scope as the chemsmart agent. I can help with computational chemistry and HPC workflows."
- Never print shell installation commands in prose. Refer installation requests to project docs and `chemsmart agent doctor`.

- Preserve native tool arguments/results. For a concrete server plus queue/job/log request, invoke the appropriate read-only tool directly.
- The ask_user tool is for a truly ambiguous required slot. For a STRUCTURED SLOT (`server`, `job_id`, log path, scheduler, queue), prose-only clarification is FORBIDDEN; call ask_user with 2-4 options when available. Do not ask again for supplied facts.
- A broad diagnostic without a concrete target is ADVISORY: answer directly (walltime, SCF, method/basis, optimization, IRC, errors).

- Keep scheduler IDs and paths typed. Never pass `last_job_id` as a `path` argument; resolve a job log via `scheduler_query` or ask for its path.
- With `last_server` and a remote path, prefer `log_tail(server=last_server, path=...)`; `read` is only for workspace files.
- For a local calculation result, call `inspect_calculation` with the remembered run ID; do not ask for an output path when one unambiguous recent receipt exists.
- Preserve active workspace project/server and prior command across turns."""


def build_system_prompt(
    *,
    registry: ToolRegistry,
    stage_instructions: str,
    session_meta: dict[str, Any] | None = None,
    conversation_context: dict[str, Any] | None = None,
    request: str = "",
    behavior_rules: str = "",
    max_chars: int | None = None,
) -> str:
    rules_block = _behavior_rules_block(behavior_rules)
    required_sections = [
        _BASE_IDENTITY_PROMPT,
        _tool_summary_block(registry, request=request),
        _session_meta_block(session_meta),
        stage_instructions.strip(),
        rules_block,
    ]
    required = "\n\n".join(section for section in required_sections if section)
    if max_chars is not None and rules_block:
        max_chars = max_chars + len(rules_block) + 2
    if max_chars is None:
        context = _conversation_context_block(
            conversation_context, max_chars=1_000_000
        )
        return "\n\n".join(
            section for section in (required, context) if section
        )

    if len(required) > max_chars:
        raise ValueError(
            f"Required system prompt is {len(required)} chars; budget is "
            f"{max_chars}. Shorten identity or stage instructions."
        )

    remaining = max_chars - len(required) - 2
    context = _conversation_context_block(
        conversation_context, max_chars=max(0, remaining)
    )
    prompt = "\n\n".join(section for section in (required, context) if section)
    if len(prompt) > max_chars:  # Defensive: projection must be hard-bounded.
        raise AssertionError("system prompt budget projection failed")
    return prompt


def _behavior_rules_block(behavior_rules: str) -> str:
    text = str(behavior_rules or "").strip()
    if not text:
        return ""
    limit = _BEHAVIOR_RULES_MAX_CHARS - len(_BEHAVIOR_RULES_HEADER) - 1
    if len(text) > limit:
        text = text[: max(0, limit - 2)] + " …"
    return f"{_BEHAVIOR_RULES_HEADER}\n{text}"


def ensure_system_message(
    messages: list[dict[str, Any]],
    system_prompt: str,
) -> list[dict[str, Any]]:
    if not messages:
        return [{"role": "system", "content": system_prompt}]

    normalized = [dict(message) for message in messages]
    if normalized[0].get("role") == "system":
        normalized[0]["content"] = system_prompt
        return normalized
    return [{"role": "system", "content": system_prompt}, *normalized]


def _tool_summary_block(registry: ToolRegistry, *, request: str) -> str:
    """Name only the intent-relevant tools; schemas travel out-of-band."""

    available = (
        {tool.name for tool in registry.list_tools()}
        if hasattr(registry, "list_tools")
        else set()
    )
    text = request.lower()
    groups = {"synthesis"}
    if any(term in text for term in ("yaml", "project", "basis")):
        groups.add("project_yaml")
    if any(term in text for term in ("execute", "submit", "run", "test")):
        groups.add("execution")
    if any(
        term in text
        for term in ("server", "queue", "scheduler", "job", "log", "ssh")
    ):
        groups.update({"diagnostics", "wizard"})
    if not request:
        groups.update({"project_yaml", "execution", "diagnostics"})
    selected = sorted(
        available.intersection(
            set().union(*(TOOL_GROUPS[group] for group in groups))
        )
    )
    # Stage packs name the routing-critical tools. The provider receives full
    # native definitions separately, so only a bounded intent hint belongs in
    # prose. Keep two representative siblings to make the scope observable.
    representatives = selected[:2]
    names = ", ".join([*representatives, "ask_user"])
    return (
        "Native schemas are separate. Relevant groups: "
        f"{','.join(sorted(groups))}; examples: {names}."
    )


def _session_meta_block(session_meta: dict[str, Any] | None) -> str:
    if not session_meta:
        return ""

    lines = ["Session metadata:"]
    for key, value in session_meta.items():
        if value is None:
            continue
        rendered = value
        if isinstance(value, (dict, list, tuple)):
            rendered = json.dumps(value, sort_keys=True)
        lines.append(f"- {key}: {rendered}")

    if len(lines) == 1:
        return ""
    return "\n".join(lines)


def _conversation_context_block(
    conversation_context: dict[str, Any] | None,
    *,
    max_chars: int,
) -> str:
    if not conversation_context or max_chars < 80:
        return ""

    recent_turns = conversation_context.get("recent_turns") or []
    older_summary = conversation_context.get("older_turn_summary") or []
    entities = conversation_context.get("entities") or {}
    if not recent_turns and not older_summary and not entities:
        return ""

    compact: dict[str, Any] = {}
    if entities:
        compact["state"] = entities
    if recent_turns:
        compact["recent"] = recent_turns[-2:]
    if older_summary:
        compact["older"] = older_summary[-2:]
    prefix = "Compact workflow state (reuse only explicit facts):\n"
    payload = json.dumps(compact, separators=(",", ":"), sort_keys=True)
    if len(prefix) + len(payload) > max_chars:
        # State entities are more reliable than a clipped natural-language
        # turn.  Drop history before using a bounded serialized state.
        payload = json.dumps(
            {"state": entities}, separators=(",", ":"), sort_keys=True
        )
    room = max_chars - len(prefix)
    if room <= 2:
        return ""
    if len(payload) > room:
        payload = payload[: max(0, room - 1)] + "…"
    return prefix + payload
