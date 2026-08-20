"""Small, evidence-preserving projections for the terminal transcript."""

from __future__ import annotations

from dataclasses import dataclass
import json
import re
from typing import TYPE_CHECKING, Any, Iterable, Mapping

from chemsmart.agent._contracts import canonical_json

_RECEIPT_PLACEHOLDER = re.compile(
    r"<(?P<role>[a-z0-9-]+):sha256=[0-9a-f]{64}>"
)


def human_cli_operation(argv: Iterable[str]) -> str:
    """Show ChemSmart input roles without exposing receipt bookkeeping."""

    return " ".join(
        _RECEIPT_PLACEHOLDER.sub(r"<\g<role>>", token) for token in argv
    )

if TYPE_CHECKING:
    from chemsmart.agent.live_session import LiveAgentSessionResultV1


_VISIBLE_TOOLS = frozenset(
    {
        # Capability, environment, and program-choice evidence.
        "inspect_program_capability",
        "inspect_program_environment",
        "assess_program_candidate",
        "consult_domain_skill",
        # Molecular identity and project-YAML authority.
        "render_project_yaml",
        "promote_project_yaml",
        "bind_scientific_identity",
        "compose_molecular_arrangement",
        "read_project_yaml",
        "validate_project_yaml",
        # Calculation and analysis DAG construction and repair.
        "plan_command_workflow",
        "plan_scientific_workflow",
        "amend_scientific_workflow",
        "inspect_workflow_frontier",
        # Live CLI compilation, safe preview, and preflight.
        "prepare_program_node",
        "synthesize_command",
        "compile_command",
        "inspect_compiled_command",
        "preview_command",
        "preflight_program_node",
        # Existing-result inspection and typed scientific analysis.
        "inspect_calculation_artifact",
        "extract_result_quantities",
        "derive_thermochemistry",
        "evaluate_quantity_expression",
        "evaluate_scientific_validation",
        "record_analysis_claims",
        "record_scientific_decision",
    }
)


@dataclass(frozen=True)
class EvidenceBlockV1:
    title: str
    text: str
    language: str = "json"


def _tool_name_by_call_id(
    transcript: Iterable[Mapping[str, Any]],
) -> dict[str, str]:
    names: dict[str, str] = {}
    for message in transcript:
        if message.get("role") != "assistant":
            continue
        for call in message.get("tool_calls") or ():
            if not isinstance(call, Mapping):
                continue
            function = call.get("function")
            if isinstance(function, Mapping):
                names[str(call.get("id") or "")] = str(
                    function.get("name") or ""
                )
    return names


def _canonical_tool_results(
    transcript: Iterable[Mapping[str, Any]],
) -> tuple[tuple[str, dict[str, Any]], ...]:
    messages = tuple(transcript)
    names = _tool_name_by_call_id(messages)
    results = []
    for message in messages:
        if message.get("role") != "tool":
            continue
        tool_name = names.get(str(message.get("tool_call_id") or ""), "")
        content = message.get("content")
        if not isinstance(content, str):
            continue
        try:
            parsed = json.loads(content)
        except json.JSONDecodeError:
            continue
        if isinstance(parsed, dict):
            results.append((tool_name, parsed))
    return tuple(results)


def _walk_key(value: Any, key: str) -> tuple[Any, ...]:
    found = []
    if isinstance(value, Mapping):
        for name, item in value.items():
            if name == key:
                found.append(item)
            found.extend(_walk_key(item, key))
    elif isinstance(value, (list, tuple)):
        for item in value:
            found.extend(_walk_key(item, key))
    return tuple(found)


def session_evidence_blocks(
    result: LiveAgentSessionResultV1,
) -> tuple[EvidenceBlockV1, ...]:
    """Return the complete public scientific tool chain used by the session.

    The TUI intentionally presents model-visible results from every general
    ChemSmart layer: capability discovery, project YAML, identity, causal DAG,
    CLI compilation, safe preview, registered-result parsing, typed
    thermochemistry, dimensional expressions, claims, and decisions.  Private
    provider reasoning and unparsed transcript text remain excluded.
    """

    blocks: list[EvidenceBlockV1] = []
    seen: set[tuple[str, str]] = set()
    for tool_name, record in _canonical_tool_results(result.public_transcript):
        for rendered_yaml in _walk_key(record, "rendered_yaml"):
            if not isinstance(rendered_yaml, str):
                continue
            marker = ("Project YAML", rendered_yaml)
            if marker not in seen:
                blocks.append(
                    EvidenceBlockV1(
                        title="Canonical project YAML",
                        text=rendered_yaml.rstrip() + "\n",
                        language="yaml",
                    )
                )
                seen.add(marker)
        if tool_name not in _VISIBLE_TOOLS:
            continue
        text = canonical_json(record)
        marker = (tool_name, text)
        if marker in seen:
            continue
        blocks.append(
            EvidenceBlockV1(
                title=tool_name.replace("_", " ").title(), text=text
            )
        )
        seen.add(marker)
    return tuple(blocks)


__all__ = [
    "EvidenceBlockV1",
    "human_cli_operation",
    "session_evidence_blocks",
]
