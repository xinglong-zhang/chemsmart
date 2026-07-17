"""Plan rendering, projection, and lightweight request classification."""

from __future__ import annotations

import json
import re
from typing import Any, Literal

from chemsmart.agent.models import Plan, Step
from chemsmart.agent.services.result_codec import preview_value

_INTENT_PATTERNS = {
    "opt": (
        r"\bopt(?:imize|imization|imisation)?\b",
        r"\bgeometry optimi[sz]ation\b",
    ),
    "ts": (r"\btransition state\b", r"\bts\b"),
    "irc": (r"\birc\b", r"\breaction path\b"),
    "sp": (
        r"\bsingle[ -]?point\b",
        r"\bsp\b",
        r"\bsingle point energy\b",
    ),
    "freq": (
        r"\bfrequenc(?:y|ies)\b",
        r"\bfreq\b",
        r"\bvibrational\b",
    ),
    "scan": (r"\bscan\b", r"\bpes\b"),
}
_CHITCHAT_EXACT_PATTERNS = (
    r"(?:hi|hello|hey)(?: there)?[!. ]*",
    r"(?:thanks|thank you|thx)[!. ]*",
    r"good (?:morning|afternoon|evening)[!. ]*",
    r"what can you do(?: for me)?\??",
    r"what do you do\??",
    r"who are you\??",
    r"help\??",
)
_CHITCHAT_IDENTITY_PATTERNS = (
    r"(?:hi|hello|hey)[!.,? ]*(?:what(?:'s| is) your name|who are you|what are you|introduce yourself)[!.,? ]*",
    r"what(?:'s| is) your name[!.,? ]*",
    r"who are you[!.,? ]*",
    r"what are you[!.,? ]*",
    r"introduce yourself[!.,? ]*",
    r"tell me about yourself[!.,? ]*",
    r"which (?:model|llm|ai) (?:are you|do you use)[!.,? ]*",
)
_CHITCHAT_TOKENS = {
    "hello",
    "hey",
    "hi",
    "thanks",
    "thank",
    "you",
    "thx",
}
_PROJECT_YAML_TOOLS = {
    "extract_project_protocol",
    "render_project_yaml",
    "validate_project_yaml",
    "critic_project_yaml",
    "write_project_yaml",
    "read_project_yaml",
    "update_project_yaml",
}


def render_plan(plan: Plan) -> str:
    lines = ["Plan:"]
    if plan.rationale:
        lines.append(f"Rationale: {plan.rationale}")
    if plan.estimated_cost:
        lines.append(f"Estimated cost: {plan.estimated_cost}")
    for index, step in enumerate(plan.steps, start=1):
        args = json.dumps(preview_value(step.args), sort_keys=True)
        lines.append(f"{index}. {step.tool} {args}")
        if step.rationale:
            lines.append(f"   - {step.rationale}")
    return "\n".join(lines)


def synthetic_plan_from_tool_requests(tool_requests: list[Any]) -> Plan:
    steps = [
        Step(tool=request.name, args=dict(request.arguments), rationale="")
        for request in tool_requests
    ]
    return Plan(
        steps=steps,
        rationale="Synthetic plan projected from tool_use_request entries.",
        intent="workflow" if steps else "advisory",
    )


def is_project_yaml_workflow(plan: Plan) -> bool:
    return any(step.tool in _PROJECT_YAML_TOOLS for step in plan.steps)


def classify_intent(request: str) -> str:
    normalized = request.lower()
    matches = [
        intent
        for intent, patterns in _INTENT_PATTERNS.items()
        if any(re.search(pattern, normalized) for pattern in patterns)
    ]
    if not matches:
        return "unknown"
    unique_matches = list(dict.fromkeys(matches))
    if len(unique_matches) == 1:
        return unique_matches[0]
    if _single_implicit_opt_intent(normalized, unique_matches):
        return next(intent for intent in unique_matches if intent != "opt")
    return "composite"


def resolve_plan_intent(
    request: str,
    plan: Plan,
) -> Literal["workflow", "advisory", "chitchat"]:
    if plan.steps:
        return "workflow"
    if is_chitchat_request(request):
        return "chitchat"
    if plan.intent in {"workflow", "advisory", "chitchat"}:
        return plan.intent
    return "advisory"


def is_chitchat_request(request: str) -> bool:
    normalized = re.sub(r"\s+", " ", request).strip().lower()
    if not normalized:
        return False
    if any(
        re.fullmatch(pattern, normalized)
        for pattern in _CHITCHAT_IDENTITY_PATTERNS
    ):
        return True
    if any(
        re.fullmatch(pattern, normalized)
        for pattern in _CHITCHAT_EXACT_PATTERNS
    ):
        return True
    tokens = re.findall(r"[a-z]+", normalized)
    return (
        bool(tokens)
        and len(tokens) <= 3
        and set(tokens).issubset(_CHITCHAT_TOKENS)
    )


def _single_implicit_opt_intent(
    normalized: str,
    unique_matches: list[str],
) -> bool:
    has_composite_marker = any(
        marker in normalized for marker in ("+", " then ", " and ", " after ")
    )
    non_opt = [intent for intent in unique_matches if intent != "opt"]
    return len(non_opt) == 1 and not has_composite_marker


__all__ = [
    "classify_intent",
    "is_chitchat_request",
    "is_project_yaml_workflow",
    "render_plan",
    "resolve_plan_intent",
    "synthetic_plan_from_tool_requests",
]
