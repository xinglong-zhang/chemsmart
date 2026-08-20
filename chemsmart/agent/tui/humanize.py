"""One English sentence per tool call, derived only from public channels.

Every visible tool gets an icon, a pending phrase, and a summarizer over the
host's canonical result envelope. The full canonical JSON stays one
expansion away in the transcript; these lines exist so a 200-step session
reads as a narrative instead of a wall of payloads. Nothing here touches a
provider channel: inputs are the append-only runtime events, whose tool
payloads the loop already restricts to public data.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

Summarizer = Callable[[Mapping[str, Any]], str]


@dataclass(frozen=True)
class HumanizedRowV1:
    icon: str
    text: str
    state: str  # running | finished | failed | note


def _inner(payload: Mapping[str, Any]) -> Mapping[str, Any]:
    envelope = payload.get("canonical_result")
    if isinstance(envelope, Mapping):
        result = envelope.get("result")
        if isinstance(result, Mapping):
            return result
        return envelope
    return {}


def _get(record: Mapping[str, Any], *path: str) -> Any:
    value: Any = record
    for key in path:
        if not isinstance(value, Mapping):
            return None
        value = value.get(key)
    return value


def _join(values: Any, limit: int = 4) -> str:
    items = [str(item) for item in (values or ())]
    shown = ", ".join(items[:limit])
    if len(items) > limit:
        shown += f", +{len(items) - limit} more"
    return shown


def _program_stage(record: Mapping[str, Any]) -> str:
    program = record.get("program") or _get(record, "query", "program") or ""
    jobtype = (
        record.get("jobtype")
        or record.get("stage")
        or _get(record, "query", "jobtype")
        or ""
    )
    engine = record.get("engine") or _get(record, "query", "engine") or ""
    parts = [str(part) for part in (program, jobtype) if part]
    label = " ".join(parts)
    if engine:
        label += f" ({engine})"
    return label


def _capability(record: Mapping[str, Any]) -> str:
    status = _get(record, "status") or _get(record, "capability", "status")
    return f"capability {_program_stage(record)} -> {status or 'observed'}"


def _environment(record: Mapping[str, Any]) -> str:
    version = record.get("observed_version") or ""
    status = record.get("status") or "observed"
    label = _program_stage(record) or str(record.get("program") or "")
    text = f"environment {label} -> {status}"
    if version:
        text += f" v{version}"
    return text


def _project_yaml(record: Mapping[str, Any]) -> str:
    artifact = record.get("artifact_id") or record.get("project_artifact_id")
    return f"project YAML {artifact or ''}".rstrip()


def _validate_yaml(record: Mapping[str, Any]) -> str:
    return (
        f"validated {record.get('project_artifact_id') or 'project'} "
        f"for {_program_stage(record)} -> {record.get('status') or 'checked'}"
    )


def _identity(record: Mapping[str, Any]) -> str:
    geometry = record.get("geometry")
    if isinstance(geometry, Mapping):
        formula = geometry.get("formula") or ""
        atoms = geometry.get("atom_count") or ""
        return f"identity bound: {formula} ({atoms} atoms)"
    return "scientific identity bound"


def _compose(record: Mapping[str, Any]) -> str:
    receipt = record.get("composition") or record
    achieved = receipt.get("achieved_contact_distance_angstrom")
    if achieved is None:
        achieved = _get(receipt, "placement", "distance_angstrom")
    label = f"composed arrangement {record.get('artifact_id') or ''}".strip()
    if achieved is not None:
        label += f" at {achieved} A contact"
    return label


def _plan_workflow(record: Mapping[str, Any]) -> str:
    plan = record.get("calculation_plan") or record
    workflow = (
        record.get("workflow_id")
        or _get(plan, "workflow_id")
        or _get(record, "scientific_toolchain_plan", "workflow_id")
        or ""
    )
    nodes = _get(plan, "nodes")
    count = len(nodes) if isinstance(nodes, (list, tuple)) else None
    analysis = _get(record, "scientific_toolchain_plan", "analysis_nodes")
    text = f"workflow {workflow}".strip()
    if count is not None:
        text += f": {count} calculation node(s)"
    if isinstance(analysis, (list, tuple)) and analysis:
        text += f" + {len(analysis)} analysis node(s)"
    return text


def _frontier(record: Mapping[str, Any]) -> str:
    frontier = record.get("workflow_frontier") or record
    actionable = _join(frontier.get("actionable_node_ids"))
    waiting = frontier.get("waiting_node_ids") or ()
    blocked = frontier.get("blocked_node_ids") or ()
    text = f"frontier: actionable [{actionable or 'none'}]"
    if waiting:
        text += f", waiting {len(waiting)}"
    if blocked:
        text += f", blocked {len(blocked)}"
    return text


def _prepared_node(record: Mapping[str, Any]) -> str:
    return f"prepared node {record.get('node_id') or ''} {_program_stage(record)}".rstrip()


def _command(record: Mapping[str, Any]) -> str:
    path = record.get("command_path") or ()
    argv = record.get("argv") or ()
    label = " ".join(str(item) for item in path) or _join(argv, limit=6)
    return f"compiled `chemsmart {label}`" if label else "compiled command"


def _preview(record: Mapping[str, Any]) -> str:
    artifacts = record.get("generated_artifacts") or record.get("artifacts")
    count = len(artifacts) if isinstance(artifacts, (list, tuple)) else None
    label = "safe preview"
    node = record.get("node_id")
    if node:
        label += f" of {node}"
    if count:
        label += f" -> {count} native artifact(s)"
    return label


def _preflight(record: Mapping[str, Any]) -> str:
    return (
        f"preflight {record.get('node_id') or _program_stage(record)}"
        f" -> {record.get('status') or 'checked'}"
    )


def _artifact_inspect(record: Mapping[str, Any]) -> str:
    return (
        f"inspected {record.get('artifact_id') or 'artifact'}"
        f" ({record.get('kind') or record.get('artifact_kind') or 'result'})"
    )


def _extraction(record: Mapping[str, Any]) -> str:
    quantities = record.get("quantities")
    ids = (
        [str(item.get("quantity_id")) for item in quantities]
        if isinstance(quantities, (list, tuple))
        and all(isinstance(item, Mapping) for item in quantities)
        else ()
    )
    return f"extracted {_join(ids) or 'quantities'}"


def _thermochemistry(record: Mapping[str, Any]) -> str:
    temperature = record.get("temperature_k")
    text = "derived thermochemistry"
    if temperature is not None:
        text += f" at {temperature} K"
    return text


def _expression(record: Mapping[str, Any]) -> str:
    value = record.get("value")
    unit = record.get("unit") or ""
    quantity = record.get("quantity_id") or ""
    if value is None:
        outputs = record.get("outputs")
        if isinstance(outputs, (list, tuple)) and outputs:
            first = outputs[0]
            if isinstance(first, Mapping):
                quantity = first.get("quantity_id") or quantity
                value = first.get("value")
                unit = first.get("unit") or unit
    if value is not None:
        return f"evaluated {quantity or 'expression'} = {value} {unit}".rstrip()
    return f"evaluated {quantity or 'expression'}"


def _validation(record: Mapping[str, Any]) -> str:
    verdict = record.get("verdict")
    if verdict is None:
        verdict = record.get("validated")
    label = record.get("validation_id") or record.get("rule_id") or "validation"
    return f"validation {label}: {verdict}"


def _claims(record: Mapping[str, Any]) -> str:
    claims = record.get("claims")
    count = len(claims) if isinstance(claims, (list, tuple)) else None
    return f"recorded {count or ''} claim(s)".replace("  ", " ")


def _decision(record: Mapping[str, Any]) -> str:
    return "recorded the scientific decision"


def _skill(record: Mapping[str, Any]) -> str:
    return f"consulted skill '{record.get('skill_id') or ''}'".replace(" ''", "")


_REGISTRY: dict[str, tuple[str, str, Summarizer]] = {
    "inspect_program_capability": ("⚛", "checking program capability", _capability),
    "inspect_program_environment": ("⚛", "probing program environment", _environment),
    "assess_program_candidate": ("⚛", "assessing program candidate", _capability),
    "consult_domain_skill": ("→", "consulting a domain skill", _skill),
    "render_project_yaml": ("←", "rendering project YAML", _project_yaml),
    "promote_project_yaml": ("←", "promoting project YAML", _project_yaml),
    "read_project_yaml": ("→", "reading project YAML", _project_yaml),
    "validate_project_yaml": ("✱", "validating project YAML", _validate_yaml),
    "bind_scientific_identity": ("⚛", "binding molecular identity", _identity),
    "compose_molecular_arrangement": ("⚛", "composing a molecular arrangement", _compose),
    "plan_command_workflow": ("✱", "planning the command workflow", _plan_workflow),
    "plan_scientific_workflow": ("✱", "planning the scientific workflow", _plan_workflow),
    "amend_scientific_workflow": ("✱", "amending the scientific workflow", _plan_workflow),
    "inspect_workflow_frontier": ("✱", "inspecting the workflow frontier", _frontier),
    "prepare_program_node": ("⚛", "preparing a program node", _prepared_node),
    "synthesize_command": ("←", "synthesizing the command", _command),
    "compile_command": ("←", "compiling the command", _command),
    "inspect_compiled_command": ("→", "inspecting the compiled command", _command),
    "preview_command": ("✱", "running the safe preview", _preview),
    "preflight_program_node": ("✱", "preflighting the node", _preflight),
    "inspect_calculation_artifact": ("→", "inspecting a calculation artifact", _artifact_inspect),
    "extract_result_quantities": ("→", "extracting typed quantities", _extraction),
    "derive_thermochemistry": ("Σ", "deriving thermochemistry", _thermochemistry),
    "evaluate_quantity_expression": ("Σ", "evaluating a quantity expression", _expression),
    "evaluate_scientific_validation": ("✱", "evaluating a validation rule", _validation),
    "record_analysis_claims": ("←", "recording analysis claims", _claims),
    "record_scientific_decision": ("←", "recording the scientific decision", _decision),
    "execute_approved_program_node": ("$", "running the approved engine node", _prepared_node),
}

_GENERIC_ICON = "⚙"


def humanize_tool_started(payload: Mapping[str, Any]) -> HumanizedRowV1:
    tool = str(payload.get("tool") or "tool")
    icon, pending, _summarize = _REGISTRY.get(
        tool, (_GENERIC_ICON, tool.replace("_", " "), None)
    )
    return HumanizedRowV1(icon=icon, text=f"~ {pending}...", state="running")


def humanize_tool_settled(
    kind: str, payload: Mapping[str, Any]
) -> HumanizedRowV1:
    tool = str(payload.get("tool") or "tool")
    icon, _pending, summarize = _REGISTRY.get(tool, (_GENERIC_ICON, "", None))
    if kind == "tool_failed":
        error = str(payload.get("error_class") or "refused")
        message = str(_get(payload, "canonical_result", "message") or "")
        text = f"{tool.replace('_', ' ')} refused ({error})"
        if message:
            text += f": {message}"
        return HumanizedRowV1(icon="✗", text=text, state="failed")
    record = _inner(payload)
    if summarize is not None:
        try:
            text = summarize(record)
        except Exception:  # pragma: no cover - a summary must never crash
            text = tool.replace("_", " ")
    else:
        keys = _join(sorted(record)[:3], limit=3)
        text = tool.replace("_", " ")
        if keys:
            text += f" [{keys}]"
    return HumanizedRowV1(icon=icon, text=text, state="finished")


def humanize_provider_turn(payload: Mapping[str, Any]) -> HumanizedRowV1:
    model = payload.get("observed_model") or payload.get("requested_model")
    return HumanizedRowV1(
        icon="▣", text=f"provider turn · {model}", state="note"
    )


__all__ = [
    "HumanizedRowV1",
    "humanize_provider_turn",
    "humanize_tool_settled",
    "humanize_tool_started",
]
