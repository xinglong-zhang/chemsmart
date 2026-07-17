"""Pure parsing and command-normalization helpers for synthesis."""

from __future__ import annotations

import json
import os
import re
import shlex
from pathlib import Path
from typing import Any

import click

from chemsmart.agent.harness.command_semantics import CommandSemanticResult
from chemsmart.agent.harness.workflow_state import current_workflow_state
from chemsmart.agent.model_command_parser import parse_model_command
from chemsmart.agent.provider_adapter import extract_response_text
from chemsmart.settings.workspace_project import (
    resolve_workspace_project,
    workspace_project_path,
)

JsonDict = dict[str, Any]
_ALLOWED_STATUSES = {"ready", "needs_clarification", "infeasible"}
_ALLOWED_CONFIDENCE = {"low", "medium", "high"}
_FRONTIER_PROVIDER_NAMES = {"openai", "anthropic"}

_THINK_BLOCK = re.compile(r"<think>(.*?)</think>", re.DOTALL | re.IGNORECASE)


def _strip_think(text: str) -> str:
    """Remove a leading <think>…</think> reasoning preamble for parsing.

    Reasoning models (Qwen3-Thinking, DeepSeek) emit their chain-of-thought
    before the JSON. It must be removed before ``json.loads`` and is never
    persisted as SFT training evidence.
    """

    return _THINK_BLOCK.sub("", text).strip()


def _extract_reasoning(response: Any) -> str:
    """Pull provider reasoning only for ephemeral response handling.

    Two conventions: a separate ``reasoning_content`` field (DashScope /
    DeepSeek / OpenAI-reasoning) or an inline ``<think>…</think>`` block in
    the content. Callers must not persist it in the SFT ledger.
    """

    if isinstance(response, dict):
        message: Any = None
        choices = response.get("choices") or []
        if choices and isinstance(choices[0], dict):
            message = choices[0].get("message") or {}
        rc = None
        if isinstance(message, dict):
            rc = message.get("reasoning_content")
        if not rc:
            rc = response.get("reasoning_content")
        if isinstance(rc, str) and rc.strip():
            return rc.strip()
    text = (
        extract_response_text(response)
        if not isinstance(response, str)
        else response
    )
    match = _THINK_BLOCK.search(text or "")
    return match.group(1).strip() if match else ""


def _parse_json_response(response: Any) -> JsonDict:
    if isinstance(response, dict):
        for key in ("parsed", "json", "content"):
            value = response.get(key)
            if isinstance(value, dict):
                return value
    text = _strip_think(extract_response_text(response).strip())
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
    parsed = json.loads(text)
    if not isinstance(parsed, dict):
        raise ValueError("synthesis response must be a JSON object")
    return parsed


def _clarification_from_semantic(
    semantic: CommandSemanticResult | None,
) -> list[str]:
    """Turn a runtime-semantic rejection into clarifying questions.

    When the safe validation run fails because a *required user input* is
    missing (charge/multiplicity being the common one), return targeted
    questions so the agent can ask the user instead of dead-ending.
    """

    if semantic is None:
        return []
    questions: list[str] = list(semantic.missing_info)
    text = " ".join(
        [
            *(issue.message for issue in semantic.issues),
            semantic.stderr_tail or "",
            semantic.stdout_tail or "",
        ]
    ).lower()
    has_charge = "charge" in text and "must be set" in text
    has_mult = "multiplicity" in text and "must be set" in text
    if has_charge and has_mult:
        questions.append(
            "What are the molecular charge and spin multiplicity? "
            "(e.g. charge 0, multiplicity 1 for neutral closed-shell)"
        )
    elif has_charge:
        questions.append("What is the molecular charge? (e.g. 0)")
    elif has_mult:
        questions.append(
            "What is the spin multiplicity? (e.g. 1 for a closed-shell singlet)"
        )
    seen: set[str] = set()
    deduped: list[str] = []
    for question in questions:
        question = str(question).strip()
        if question and question not in seen:
            seen.add(question)
            deduped.append(question)
    return deduped


def _normalize_result(result: JsonDict) -> JsonDict:
    normalized = dict(result)
    status = normalized.get("status")
    confidence = normalized.get("confidence")
    if status not in _ALLOWED_STATUSES:
        raise ValueError(f"invalid synthesis status: {status!r}")
    if confidence not in _ALLOWED_CONFIDENCE:
        raise ValueError(f"invalid confidence: {confidence!r}")
    normalized.setdefault("command", "")
    normalized.setdefault("explanation", "")
    for key in ("missing_info", "alternatives"):
        value = normalized.get(key)
        if value is None:
            normalized[key] = []
        if not isinstance(normalized[key], list):
            raise ValueError(f"{key} must be a list")
    if status == "ready" and not str(normalized.get("command", "")).startswith(
        "chemsmart"
    ):
        raise ValueError("ready command must start with 'chemsmart'")
    return normalized


def _is_v8_spec(result: JsonDict) -> bool:
    intent = result.get("intent")
    return (
        isinstance(intent, str)
        and (intent in {"workflow", "advisory", "decline", "chitchat"})
        and ("jobs" in result or "message" in result)
    )


def _normalize_v8_spec(
    result: JsonDict,
    *,
    default_project: str | None = None,
) -> JsonDict:
    from chemsmart.agent.v8_adapter import adapt

    adapted = adapt(result, validate=True, default_project=default_project)
    intent = str(adapted.get("intent") or result.get("intent") or "")
    if intent != "workflow":
        message = str(adapted.get("message") or result.get("message") or "")
        status = "infeasible" if intent == "decline" else "needs_clarification"
        if intent == "chitchat":
            status = "infeasible"
        return {
            "status": status,
            "command": "",
            "explanation": message or "No executable workflow was requested.",
            "confidence": "high",
            "missing_info": [],
            "alternatives": [],
        }

    commands = [
        command
        for command in adapted.get("commands", [])
        if isinstance(command, str) and command.strip()
    ]
    if not commands:
        raise ValueError(
            "v8 compact SPEC rendered no commands: "
            + "; ".join(map(str, adapted.get("errors") or []))
        )
    if adapted.get("valid") is False:
        raise ValueError(
            "v8 compact SPEC failed chemsmart validation: "
            + "; ".join(map(str, adapted.get("errors") or []))
        )
    return {
        "status": "ready",
        "command": commands[0],
        "explanation": "Prepared chemsmart command from compact SPEC.",
        "confidence": "high" if adapted.get("valid") else "medium",
        "missing_info": [],
        "alternatives": commands[1:],
    }


def resolve_default_project() -> str:
    selected = current_workflow_state().project
    if selected is not None and _workspace_project_exists(selected.name):
        return selected.name
    env_project = os.environ.get("CHEMSMART_AGENT_PROJECT", "").strip()
    if env_project and _workspace_project_exists(env_project):
        return env_project
    try:
        from chemsmart.agent.provider_config import load_active_provider_config

        config = load_active_provider_config()
    except Exception:
        config = None
    if (
        config is not None
        and config.project
        and _workspace_project_exists(config.project.strip())
    ):
        return config.project.strip()
    workspace = resolve_workspace_project()
    if workspace.loaded:
        return workspace.project
    return ""


def _workspace_project_exists(project: str) -> bool:
    return any(
        workspace_project_path(project, program).exists()
        for program in ("gaussian", "orca")
    )


def _is_frontier_provider(provider: Any) -> bool:
    return (
        str(getattr(provider, "name", "")).lower() in _FRONTIER_PROVIDER_NAMES
    )


def _extract_command_from_text(text: str) -> str:
    if "chemsmart" not in text:
        return ""
    for line in text.splitlines():
        stripped = line.strip().strip("`")
        if stripped.startswith("chemsmart "):
            return stripped
    match = re.search(r"(chemsmart\s+.+)", text)
    if not match:
        return ""
    command = match.group(1).strip().strip("`")
    return command


def _request_needs_workspace_project(request: str) -> bool:
    lowered = request.lower()
    if "chemsmart " in lowered and (
        "explain" in lowered or "what" in lowered or "mean" in lowered
    ):
        return False
    job_markers = (
        "gaussian",
        "orca",
        "optimization",
        "optimisation",
        "transition",
        "single point",
        "frequency",
        "td-dft",
        "tddft",
        "scan",
        "irc",
        "neb",
        "dias",
        "wbi",
        "run",
        "submit",
        "prepare",
        "set up",
    )
    return any(marker in lowered for marker in job_markers)


def _build_decision_trace(
    *,
    decision: JsonDict,
    action: str,
    request: str,
    target_command: str,
    default_project: str,
    last_command: str,
) -> JsonDict:
    evidence = decision.get("evidence")
    if not isinstance(evidence, list):
        evidence = []
    rejected = decision.get("rejected_actions")
    if not isinstance(rejected, dict):
        rejected = {}
    trace: JsonDict = {
        "router": "api_frontier_intent_router",
        "action": action,
        "confidence": str(decision.get("confidence") or "medium"),
        "decision_summary": str(decision.get("decision_summary") or ""),
        "target_command": target_command,
        "default_project": default_project,
        "last_command_available": bool(last_command),
        "request_excerpt": request[:240],
        "evidence": [str(item) for item in evidence if str(item).strip()],
        "rejected_actions": {
            str(key): str(value) for key, value in rejected.items()
        },
        "note": (
            "This is a public decision trace, not hidden chain-of-thought. "
            "It records observable routing evidence for user audit."
        ),
    }
    if not trace["decision_summary"]:
        trace["decision_summary"] = _fallback_decision_summary(action)
    if not trace["evidence"]:
        trace["evidence"] = _fallback_decision_evidence(
            action,
            target_command=target_command,
            last_command=last_command,
        )
    if not trace["rejected_actions"]:
        trace["rejected_actions"] = _fallback_rejected_actions(action)
    return trace


def _fallback_decision_summary(action: str) -> str:
    summaries = {
        "synthesize_command": "The user requested a computational job command.",
        "explain_command": "The user asked to interpret an existing command.",
        "critique_command": "The user asked to judge command validity.",
        "repair_command": "The user asked to fix a command or validation failure.",
        "clarify": "The request does not contain enough command or job detail.",
    }
    return summaries.get(
        action, "The request was routed by the API intent gate."
    )


def _fallback_decision_evidence(
    action: str,
    *,
    target_command: str,
    last_command: str,
) -> list[str]:
    evidence = [f"Selected action: {action}."]
    if target_command:
        evidence.append("A chemsmart command is available for this turn.")
    elif last_command:
        evidence.append("A previous chemsmart command is available in memory.")
    else:
        evidence.append("No chemsmart command was available in memory.")
    return evidence


def _fallback_rejected_actions(action: str) -> dict[str, str]:
    rejected: dict[str, str] = {}
    if action != "synthesize_command":
        rejected["synthesize_command"] = (
            "The request was not primarily asking for a new CLI command."
        )
    if action != "explain_command":
        rejected["explain_command"] = (
            "The request was not primarily asking for a command explanation."
        )
    if action != "critique_command":
        rejected["critique_command"] = (
            "The request was not primarily asking for validation judgment."
        )
    if action != "repair_command":
        rejected["repair_command"] = (
            "The request was not primarily asking for command repair."
        )
    return rejected


def _ensure_program_project(command: str, project: str) -> str:
    project = project.strip()
    if not project:
        return command
    try:
        tokens = shlex.split(command)
    except ValueError:
        return command
    program_index = _find_program_index(tokens)
    if program_index is None:
        return command
    updated = list(tokens)
    for index in range(program_index + 1, len(updated)):
        token = updated[index]
        if token in {"-p", "--project"}:
            if index + 1 >= len(updated):
                return command
            if _is_project_placeholder(updated[index + 1]):
                updated[index + 1] = project
                return shlex.join(updated)
            return command
        if token.startswith("--project="):
            value = token.split("=", 1)[1]
            if _is_project_placeholder(value):
                updated[index] = f"--project={project}"
                return shlex.join(updated)
            return command
    updated[program_index + 1 : program_index + 1] = ["-p", project]
    return shlex.join(updated)


def _selected_project_for_command(
    command: str,
    default_project: str,
) -> Any | None:
    """Return the selected workspace project only when it matches the command."""

    selected = current_workflow_state().project
    if (
        selected is None
        or selected.name != default_project
        or not Path(selected.path).is_file()
    ):
        return None
    parsed = parse_model_command(command)
    if parsed.parse_error or parsed.program != selected.program:
        return None
    return selected


def _command_project_is_unresolved(command: str) -> bool:
    try:
        tokens = shlex.split(command)
    except ValueError:
        return False
    program_index = _find_program_index(tokens)
    if program_index is None:
        return False
    for index in range(program_index + 1, len(tokens)):
        token = tokens[index]
        if token in {"-p", "--project"}:
            return index + 1 < len(tokens) and _is_project_placeholder(
                tokens[index + 1]
            )
        if token.startswith("--project="):
            return _is_project_placeholder(token.split("=", 1)[1])
    return True


def _is_project_placeholder(value: str) -> bool:
    return value.strip().lower() in {
        "<project>",
        "<project-name>",
        "<project_name>",
        "{project}",
        "{project_name}",
    }


def _is_project_name_missing(value: str) -> bool:
    normalized = re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()
    if "project" not in normalized:
        return False
    return any(
        marker in normalized
        for marker in (
            "project name",
            "project yaml",
            "selected project",
            "select project",
            "which project",
            "workspace project",
        )
    )


def _find_program_index(tokens: list[str]) -> int | None:
    if not tokens or tokens[0] != "chemsmart":
        return None
    top_index = None
    for index, token in enumerate(tokens[1:], start=1):
        if token in {"run", "sub"}:
            top_index = index
            break
    if top_index is None:
        return None

    index = top_index + 1
    while index < len(tokens):
        token = tokens[index]
        if token in {"gaussian", "orca"}:
            return index
        if token.startswith("-"):
            index += 1 if token in _TOP_LEVEL_FLAGS else 2
            continue
        return None
    return None


_TOP_LEVEL_FLAGS = {
    "--fake",
    "--no-fake",
    "--scratch",
    "--no-scratch",
    "--test",
}


def _program_has_project(tokens: list[str], program_index: int) -> bool:
    return any(
        token in {"-p", "--project"} for token in tokens[program_index + 1 :]
    )


def _validate_tokens_against_schema(
    tokens: list[str], schema: JsonDict
) -> None:
    node = schema
    path_nodes = [schema]
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if token == "--":
            return
        subcommands = node.get("subcommands") or {}
        if not token.startswith("-") and token in subcommands:
            node = subcommands[token]
            path_nodes.append(node)
            index += 1
            continue
        if not token.startswith("-") and subcommands:
            raise click.ClickException(f"unknown subcommand: {token}")
        if token.startswith("-"):
            option_token, has_inline_value = _split_option_token(token)
            option = _find_option(option_token, path_nodes)
            if option is None:
                raise click.NoSuchOption(option_token)
            if bool(option.get("is_flag")) or has_inline_value:
                index += 1
                continue
            nargs = int(option.get("nargs") or 1)
            consume = max(1, nargs)
            if index + consume >= len(tokens):
                raise click.BadOptionUsage(
                    option_token, "option requires a value"
                )
            index += 1 + consume
            continue
        # Positional arguments are legal for many ChemSmart commands; click
        # performs command-path parsing above, while this pass blocks invented
        # options regardless of where the model placed inherited options.
        index += 1


def _split_option_token(token: str) -> tuple[str, bool]:
    if token.startswith("--") and "=" in token:
        option, _value = token.split("=", 1)
        return option, True
    return token, False


def _find_option(option_token: str, nodes: list[JsonDict]) -> JsonDict | None:
    for node in reversed(nodes):
        for option in node.get("options") or []:
            if option_token in (option.get("opts") or []):
                return option
    return None
