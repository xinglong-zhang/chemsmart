"""Interactive natural-language to ChemSmart CLI command synthesis."""

from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
from pathlib import Path
from typing import Any

import click
from rich.console import Console
from rich.panel import Panel

from chemsmart.agent.cli_schema import build_chemsmart_cli_schema
from chemsmart.agent.harness.command_semantics import (
    CommandSemanticResult,
    evaluate_command_semantics,
)
from chemsmart.agent.command_answerer import (
    ComposedAnswer,
    compose_command_answer,
    reason_missing_info,
)
from chemsmart.agent.harness.spec_invariants import check_spec
from chemsmart.agent.kind_disambiguator import disambiguate
from chemsmart.agent.model_command_parser import (
    parse_model_command,
)

# Structural spec-invariant rule ids that are reliable to enforce at synthesis
# time. The query-heuristic rules (decline_contract / required_present / label)
# are intentionally excluded here to avoid false rejects on well-formed specs.
_RELIABLE_SPEC_REJECTS = frozenset({
    "spec.kind.canonical",
    "spec.runtime_owned.leak",
    "spec.settings.allowed",
    "spec.ts.runtime_freq",
    "spec.ts.runtime_route",
    "spec.ts.bad_extra",
    "spec.jobs.empty",
})


def _apply_spec_invariants(
    result: "JsonDict", spec: "JsonDict", query: str
) -> None:
    """Run the deterministic SPEC invariants on the model's v8 SPEC and, on a
    reliable structural violation, downgrade a ``ready`` result so a malformed
    spec (runtime-owned leak, disallowed setting, redundant TS route token, …)
    is never surfaced as an executable command."""
    try:
        issues = check_spec(spec, query)
    except Exception:  # pragma: no cover - never fail synthesis on the gate
        return
    reliable = [i for i in issues if i.rule_id in _RELIABLE_SPEC_REJECTS]
    if not reliable:
        return
    ids = sorted({i.rule_id for i in reliable})
    msgs = [f"{i.rule_id}: {i.message}" for i in reliable]
    result["spec_invariant_failed"] = ids
    if result.get("status") == "ready":
        result["status"] = "infeasible"
        result["explanation"] = (
            (result.get("explanation") or "")
            + " Spec-invariant check rejected: "
            + "; ".join(msgs)
        ).strip()
from chemsmart.agent.prompts.synthesis import build_synthesis_system_prompt
from chemsmart.agent.services.conversation_memory import (
    ConversationMemory,
    ConversationTurn,
)

JsonDict = dict[str, Any]
_ALLOWED_STATUSES = {"ready", "needs_clarification", "infeasible"}
_ALLOWED_CONFIDENCE = {"low", "medium", "high"}
_FRONTIER_PROVIDER_NAMES = {"openai", "anthropic"}
_ROUTER_ACTIONS = {
    "synthesize_command",
    "explain_command",
    "critique_command",
    "repair_command",
    "clarify",
}
_INTENT_ROUTER_SYSTEM_PROMPT = """You are the intent router for the ChemSmart agent UI.

Return ONLY one JSON object with exactly these fields:
{
  "action": "synthesize_command" | "explain_command" | "critique_command" | "repair_command" | "clarify",
  "target_command": "chemsmart ..." | "",
  "confidence": "low" | "medium" | "high",
  "missing_info": ["question", ...],
  "decision_summary": "one short sentence",
  "evidence": ["short observable reason", ...],
  "rejected_actions": {"action_name": "short reason", ...}
}

Choose synthesize_command when the user asks to prepare, set up, run, submit, or change a computational job.
Choose explain_command when the user asks what an existing command does, what it means, or how to read it.
Choose critique_command when the user asks whether a command is valid, safe, correctly grounded, or will work.
Choose repair_command when the user asks to fix a command or recover from a failed validation.
Choose clarify only when neither the current request nor the last command is enough.

Use target_command from the request if it contains one; otherwise use last_command when the user refers to "this command", "that", "it", or the previous result. Do not reveal hidden chain-of-thought."""


class SynthesisSession:
    """Synthesize, validate, confirm, and execute one ChemSmart CLI command."""

    def __init__(
        self,
        provider: Any | None = None,
        schema: JsonDict | None = None,
        max_retries: int = 2,
        debug: bool = False,
        semantic_gate: bool = True,
        semantic_timeout_s: float = 30.0,
        default_project: str | None = None,
        enable_intent_router: bool | None = None,
    ) -> None:
        """Initialize a synthesis session.

        Args:
            provider: Optional provider test double. Defaults to the active
                configured provider.
            schema: Optional precomputed CLI schema.
            max_retries: Number of corrective retries for invalid JSON.
            debug: When True, executed chemsmart subprocesses keep their
                default ``--verbose`` behavior; when False (default), a
                ``--no-verbose`` flag is injected so the child stays quiet.
            semantic_gate: When True, ready run/sub commands must pass the
                safe runtime semantic gate before user confirmation.
            semantic_timeout_s: Timeout for the safe semantic validation run.
            default_project: Runtime-owned Gaussian/ORCA project name to inject
                when a synthesized command omits ``-p/--project``.
            enable_intent_router: Enable the API/frontier-only intent router
                that decides whether a turn should synthesize, explain, critique,
                or repair a command before the command-synthesis schema runs.
        """

        if provider is None:
            from chemsmart.agent import providers

            provider = providers.get_provider()
            if default_project is None:
                default_project = resolve_default_project()
        self.provider = provider
        self.default_project = (default_project or "").strip()
        self.schema = schema or build_chemsmart_cli_schema()
        self.max_retries = max_retries
        self.debug = debug
        self.semantic_gate = semantic_gate
        self.semantic_timeout_s = semantic_timeout_s
        if enable_intent_router is None:
            enable_intent_router = _is_frontier_provider(provider)
        self.enable_intent_router = bool(enable_intent_router)
        self.memory = ConversationMemory()
        self._last_raw_response = ""
        self._last_semantic_result: CommandSemanticResult | None = None
        self._last_decision_trace: JsonDict | None = None

    def synthesize(self, request: str) -> JsonDict:
        """Return a structured synthesis result for ``request``.

        Invalid JSON is retried with a corrective message up to
        ``max_retries`` times.
        """

        messages = self._messages_for_request(request)
        last_error: Exception | None = None
        for attempt in range(self.max_retries + 1):
            response = self.provider.chat(messages)
            self._last_raw_response = _extract_text(response)
            # The local provider pre-adapts the model output into a status/command
            # result; replay the model's raw SPEC (when provided) as assistant history
            # so multi-turn follow-ups see the schema the model actually emits.
            if isinstance(response, dict) and isinstance(
                response.get("raw_plan"), str
            ):
                self._last_raw_response = response["raw_plan"]
            try:
                parsed = _parse_json_response(response)
                if _is_v8_spec(parsed):
                    parsed, _changed = disambiguate(request, parsed)
                    result = _normalize_v8_spec(
                        parsed,
                        default_project=self.default_project or None,
                    )
                    _apply_spec_invariants(result, parsed, request)
                    return self._apply_default_project(result)
                return self._apply_default_project(_normalize_result(parsed))
            except (TypeError, ValueError, json.JSONDecodeError) as exc:
                last_error = exc
                if attempt >= self.max_retries:
                    break
                messages.append(
                    {
                        "role": "assistant",
                        "content": self._last_raw_response,
                    }
                )
                messages.append(
                    {
                        "role": "user",
                        "content": (
                            "Your previous response was invalid. Return ONLY "
                            "valid JSON matching the required schema. Error: "
                            f"{exc}"
                        ),
                    }
                )
        raise ValueError(
            f"LLM failed to return valid synthesis JSON: {last_error}"
        )

    def validate_command(self, command: str) -> tuple[bool, str]:
        """Dry-parse ``command`` and return ``(ok, error)``."""

        stripped = command.strip()
        try:
            tokens = shlex.split(stripped)
        except ValueError as exc:
            return False, str(exc)
        if not tokens or tokens[0] != "chemsmart":
            return False, "command must start with 'chemsmart'"
        if len(tokens) == 1:
            return False, "command must include a chemsmart subcommand"

        from chemsmart.cli.main import entry_point

        remainder = tokens[1:]
        try:
            entry_point.make_context(
                "chemsmart",
                remainder,
                resilient_parsing=True,
            )
            _validate_tokens_against_schema(remainder, self.schema)
        except Exception as exc:
            return False, str(exc)
        return True, ""

    def run_interactive(self, request: str) -> None:
        """Run the interactive synthesis, clarification, validation flow."""

        current_request = request
        result = self.synthesize(current_request)
        while result["status"] == "needs_clarification":
            answers: list[str] = []
            missing_info = result.get("missing_info") or ["clarification"]
            for item in missing_info:
                answer = click.prompt(str(item))
                answers.append(f"{item}: {answer}")
            self._remember_turn(current_request, "; ".join(answers))
            current_request = (
                f"{current_request}\nClarifications: " + "; ".join(answers)
            )
            result = self.synthesize(current_request)

        if result["status"] == "infeasible":
            click.echo(result.get("explanation") or "Request is infeasible.")
            return

        result = self._repair_ready_result(current_request, result)
        if result is None:
            return
        command = str(result.get("command") or "")

        self._remember_turn(
            current_request,
            command,
            assistant_message=self._last_raw_response,
        )
        self.confirm_and_execute(
            command,
            str(result.get("explanation") or ""),
            str(result.get("confidence") or "low"),
        )

    def prepare_command(self, request: str) -> JsonDict:
        """Synthesize and validate a command without prompting or executing.

        This is the non-interactive half of ``run_interactive`` used by the TUI:
        it lets the user chat with the active synthesis provider, then inspect
        the rendered command and runtime semantic gate result before deciding
        whether to move into the full harness workflow.
        """

        self._last_semantic_result = None
        self._last_decision_trace = None
        routed = self._route_frontier_request(request)
        if routed is not None:
            return routed

        result = self.synthesize(request)
        if result["status"] != "ready":
            return self._attach_decision_trace(result)

        repaired = self._repair_ready_result(request, result)
        if repaired is None:
            return self._attach_decision_trace({
                "status": "infeasible",
                "command": "",
                "explanation": (
                    "Synthesized command failed validation and could not be "
                    "repaired."
                ),
                "confidence": "low",
                "missing_info": [],
                "alternatives": [],
            })
        command = str(repaired.get("command") or "")
        self._remember_turn(
            request,
            command,
            assistant_message=self._last_raw_response,
        )
        return self._attach_decision_trace(repaired)

    def _route_frontier_request(self, request: str) -> JsonDict | None:
        if not self.enable_intent_router:
            return None
        decision = self._frontier_intent_decision(request)
        if decision is None:
            return None
        action = str(decision.get("action") or "synthesize_command")
        target_command = _extract_command_from_text(request) or str(
            decision.get("target_command") or ""
        ).strip() or self._last_ready_command()
        self._last_decision_trace = _build_decision_trace(
            decision=decision,
            action=action,
            request=request,
            target_command=target_command,
            default_project=self.default_project,
            last_command=self._last_ready_command(),
        )
        if action == "synthesize_command":
            return None
        if action in {"explain_command", "critique_command", "repair_command"}:
            if not target_command:
                return self._attach_decision_trace({
                    "status": "needs_clarification",
                    "command": "",
                    "explanation": "No previous chemsmart command is available.",
                    "confidence": "medium",
                    "missing_info": [
                        "Paste the chemsmart command you want me to explain or critique."
                    ],
                    "alternatives": [],
                    "action": action,
                })
            target_command = _ensure_program_project(
                target_command,
                self.default_project,
            )
            if self._last_decision_trace is not None:
                self._last_decision_trace["target_command"] = target_command
            parsed_command = parse_model_command(target_command)
            semantic_summary: JsonDict | None = None
            if action in {"critique_command", "repair_command"}:
                self._last_semantic_result = evaluate_command_semantics(
                    target_command,
                    timeout_s=self.semantic_timeout_s,
                )
                semantic_summary = self._last_semantic_result.to_dict()
            # The deterministic parser + semantic gate are the grounding tools;
            # the provider composes the user-facing answer in the user's own
            # language, and the grounding harness rejects any contradiction.
            composed = compose_command_answer(
                self.provider,
                request=request,
                action=action,
                parsed=parsed_command,
                semantic_summary=semantic_summary,
                timeout_s=self.semantic_timeout_s,
            )
            self._record_answer_reasoning(composed)
            return self._attach_decision_trace({
                "status": "informational",
                "command": target_command,
                "explanation": composed.answer,
                "confidence": str(decision.get("confidence") or "high"),
                "missing_info": [],
                "alternatives": [],
                "action": action,
                "project": self.default_project,
                "reasoning": list(composed.reasoning),
                "caveats": list(composed.caveats),
                "grounded": composed.grounded,
            })
        if action == "clarify":
            missing = decision.get("missing_info")
            if not isinstance(missing, list) or not missing:
                missing = []
            summary = str(decision.get("decision_summary") or "")
            reasoning: list[str] = []
            cot = reason_missing_info(
                self.provider,
                request=request,
                hints=tuple(str(item) for item in missing),
                timeout_s=self.semantic_timeout_s,
            )
            if cot is not None:
                if cot.questions:
                    missing = list(cot.questions)
                reasoning = list(cot.reasoning)
                if cot.reasoning and self._last_decision_trace is not None:
                    self._last_decision_trace["reasoning"] = reasoning
            if not missing:
                missing = ["What command or calculation should I work with?"]
            return self._attach_decision_trace({
                "status": "needs_clarification",
                "command": "",
                "explanation": summary,
                "confidence": str(decision.get("confidence") or "medium"),
                "missing_info": missing,
                "alternatives": [],
                "action": action,
                "reasoning": reasoning,
            })
        return None

    def _record_answer_reasoning(self, composed: ComposedAnswer) -> None:
        """Surface the composed answer's public reasoning on the trace."""

        if self._last_decision_trace is None:
            return
        if composed.reasoning:
            self._last_decision_trace["reasoning"] = list(composed.reasoning)
        if composed.caveats:
            self._last_decision_trace["caveats"] = list(composed.caveats)
        self._last_decision_trace["answer_grounded"] = composed.grounded
        if composed.fallback_used and composed.fallback_reason:
            self._last_decision_trace["fallback_reason"] = (
                composed.fallback_reason
            )

    def _attach_decision_trace(self, result: JsonDict) -> JsonDict:
        if not self._last_decision_trace:
            return result
        normalized = dict(result)
        normalized.setdefault("decision_trace", self._last_decision_trace)
        return normalized

    def _frontier_intent_decision(self, request: str) -> JsonDict | None:
        last_command = self._last_ready_command()
        messages = [
            {
                "role": "system",
                "content": _INTENT_ROUTER_SYSTEM_PROMPT,
            },
            {
                "role": "user",
                "content": json.dumps(
                    {
                        "request": request,
                        "last_command": last_command,
                        "default_project": self.default_project,
                    },
                    ensure_ascii=False,
                ),
            },
        ]
        try:
            parsed = _parse_json_response(self.provider.chat(messages))
        except Exception:
            return None
        action = str(parsed.get("action") or "")
        if action not in _ROUTER_ACTIONS:
            return None
        return parsed

    def _last_ready_command(self) -> str:
        for turn in reversed(self.memory.turns):
            command = str(turn.plan_rationale or "").strip()
            if command.startswith("chemsmart "):
                return command
        return ""

    def _repair_ready_result(
        self,
        current_request: str,
        result: JsonDict,
    ) -> JsonDict | None:
        """Validate a ready result and ask the model to repair gate failures."""

        semantic_retried = False
        for attempt in range(self.max_retries + 1):
            if result["status"] != "ready":
                click.echo(
                    result.get("explanation")
                    or "Could not synthesize command."
                )
                return None

            command = str(result.get("command") or "")
            valid, error = self.validate_command(command)
            retry_request: str | None = None
            if not valid:
                retry_request = (
                    f"Original request: {current_request}\n"
                    f"Your command failed validation: {error}\n"
                    "Return a corrected legal command."
                )
            elif self.semantic_gate:
                semantic_result = evaluate_command_semantics(
                    command,
                    timeout_s=self.semantic_timeout_s,
                )
                self._last_semantic_result = semantic_result
                if semantic_result.verdict == "reject":
                    semantic_retried = True
                    click.echo(
                        "Notice: synthesized command failed runtime semantic "
                        "validation; asking the model to repair it."
                    )
                    retry_request = semantic_result.correction_prompt(
                        current_request
                    )
                elif semantic_result.verdict == "warn":
                    click.echo(
                        "Notice: synthesized command passed with runtime "
                        "semantic warnings."
                    )
                else:
                    if semantic_retried:
                        click.echo(
                            "Notice: synthesized command was repaired after "
                            "runtime semantic validation."
                        )
                    return result
            else:
                return result

            if retry_request is None:
                return result
            if attempt >= self.max_retries:
                if self.semantic_gate and self._last_semantic_result:
                    click.echo(
                        "Synthesized command failed runtime semantic "
                        "validation: "
                        + "; ".join(
                            issue.message
                            for issue in self._last_semantic_result.issues
                        )
                    )
                else:
                    click.echo(f"Synthesized command is invalid: {error}")
                return None
            result = self.synthesize(retry_request)
        return None

    def confirm_and_execute(
        self,
        cmd: str,
        explanation: str,
        confidence: str,
    ) -> None:
        """Ask the user to confirm, edit, cancel, or test the command."""

        console = Console()
        panel_text = (
            f"[bold]Command[/bold]\n{cmd}\n\n"
            f"[bold]Explanation[/bold]\n{explanation or '(none)'}\n\n"
            f"[bold]Confidence[/bold] {confidence}"
        )
        console.print(Panel(panel_text, title="ChemSmart synthesis"))
        choice = (
            click.prompt(
                "Run? [Y]es/[E]dit/[N]o/[T]est",
                default="Y",
                show_default=False,
            )
            .strip()
            .upper()
        )

        if choice in {"", "Y", "YES"}:
            subprocess.run(self._quiet_argv(shlex.split(cmd)), check=False)
            return
        if choice in {"N", "NO"}:
            click.echo("cancelled")
            return
        if choice in {"E", "EDIT"}:
            edited = click.edit(cmd)
            if edited is None:
                click.echo("cancelled")
                return
            edited = edited.strip()
            valid, error = self.validate_command(edited)
            if not valid:
                click.echo(f"Edited command is invalid: {error}")
                return
            self.confirm_and_execute(edited, explanation, confidence)
            return
        if choice in {"T", "TEST"}:
            tokens = shlex.split(cmd)
            if len(tokens) > 1 and tokens[1] == "sub":
                tokens.insert(2, "--test")
                subprocess.run(self._quiet_argv(tokens), check=False)
            else:
                click.echo("test only for sub")
            return
        click.echo("cancelled")

    def _quiet_argv(self, tokens: list[str]) -> list[str]:
        """Inject ``--no-verbose`` after ``chemsmart`` unless debug is on.

        The top-level ``chemsmart`` CLI defaults to ``--verbose=True`` which
        prints the ASCII banner and INFO logs to stdout. Quiet mode is the
        user-facing default for synthesis-driven invocations; ``--debug``
        passes the command through unchanged so the child stays loud.
        """
        if self.debug:
            return list(tokens)
        if not tokens or tokens[0] != "chemsmart":
            return list(tokens)
        if "--no-verbose" in tokens or "--verbose" in tokens:
            return list(tokens)
        argv = list(tokens)
        argv.insert(1, "--no-verbose")
        return argv

    def _messages_for_request(self, request: str) -> list[dict[str, str]]:
        messages = [
            {
                "role": "system",
                "content": build_synthesis_system_prompt(self.schema),
            }
        ]
        for turn in self.memory.turns:
            if not turn.assistant_message:
                continue
            messages.append({"role": "user", "content": turn.request})
            messages.append(
                {"role": "assistant", "content": turn.assistant_message}
            )
        messages.append({"role": "user", "content": request})
        return messages

    def _apply_default_project(self, result: JsonDict) -> JsonDict:
        """Attach the active runtime project to ready Gaussian/ORCA commands."""

        if result.get("status") != "ready" or not self.default_project:
            return result
        normalized = dict(result)
        command = str(normalized.get("command") or "")
        normalized["command"] = _ensure_program_project(
            command,
            self.default_project,
        )
        alternatives = normalized.get("alternatives") or []
        if isinstance(alternatives, list):
            normalized["alternatives"] = [
                (
                    _ensure_program_project(str(item), self.default_project)
                    if isinstance(item, str)
                    else item
                )
                for item in alternatives
            ]
        normalized["project"] = self.default_project
        return normalized

    def _remember_turn(
        self,
        request: str,
        result: str,
        assistant_message: str = "",
    ) -> None:
        turn = ConversationTurn(
            turn_index=len(self.memory.turns) + 1,
            request=request,
            plan_rationale=result,
            assistant_message=assistant_message,
            status="completed",
        )
        self.memory.turns.append(turn)


def _parse_json_response(response: Any) -> JsonDict:
    if isinstance(response, dict):
        for key in ("parsed", "json", "content"):
            value = response.get(key)
            if isinstance(value, dict):
                return value
    text = _extract_text(response).strip()
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
    parsed = json.loads(text)
    if not isinstance(parsed, dict):
        raise ValueError("synthesis response must be a JSON object")
    return parsed


def _extract_text(response: Any) -> str:
    if isinstance(response, str):
        return response
    if isinstance(response, dict):
        content = response.get("content")
        if isinstance(content, str):
            return content
        if isinstance(content, list):
            parts = [
                str(item.get("text", ""))
                for item in content
                if isinstance(item, dict) and item.get("type") == "text"
            ]
            if parts:
                return "\n".join(parts)
        choices = response.get("choices") or []
        if choices:
            message = choices[0].get("message", {})
            content = message.get("content", "")
            if isinstance(content, str):
                return content
            if isinstance(content, list):
                return "\n".join(
                    str(item.get("text", ""))
                    for item in content
                    if isinstance(item, dict)
                )
    raise ValueError("could not extract text from provider response")


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
    return isinstance(intent, str) and (
        intent in {"workflow", "advisory", "decline", "chitchat"}
    ) and ("jobs" in result or "message" in result)


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
    env_project = os.environ.get("CHEMSMART_AGENT_PROJECT", "").strip()
    if env_project:
        return env_project
    try:
        from chemsmart.agent.provider_config import load_active_provider_config

        config = load_active_provider_config()
    except Exception:
        config = None
    if config is not None and config.project:
        return config.project.strip()
    return _single_available_project()


def _single_available_project() -> str:
    candidates: set[str] = set()
    for folder_name in ("gaussian", "orca"):
        folder = Path.home() / ".chemsmart" / folder_name
        if not folder.is_dir():
            continue
        for path in folder.iterdir():
            if path.is_file() and path.suffix.lower() in {".yaml", ".yml"}:
                candidates.add(path.stem)
    return next(iter(candidates)) if len(candidates) == 1 else ""


def _is_frontier_provider(provider: Any) -> bool:
    return str(getattr(provider, "name", "")).lower() in _FRONTIER_PROVIDER_NAMES


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
    return summaries.get(action, "The request was routed by the API intent gate.")


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
    if _program_has_project(tokens, program_index):
        return command
    updated = list(tokens)
    updated[program_index + 1 : program_index + 1] = ["-p", project]
    return shlex.join(updated)


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
    return any(token in {"-p", "--project"} for token in tokens[program_index + 1 :])


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
