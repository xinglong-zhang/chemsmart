"""Command preparation, routing, repair, and confirmation workflow."""

from __future__ import annotations

import json
import shlex
from typing import Any

import click
from rich.console import Console
from rich.panel import Panel

from chemsmart.agent.command_answerer import (
    ComposedAnswer,
    compose_command_answer,
    reason_missing_info,
)
from chemsmart.agent.harness.workflow_state import (
    current_workflow_state,
    select_project_from_request,
    select_workspace_project,
)
from chemsmart.agent.model_command_parser import parse_model_command
from chemsmart.agent.services.synthesis_support import (
    _build_decision_trace,
    _command_project_is_unresolved,
    _ensure_program_project,
    _extract_command_from_text,
    _is_project_name_missing,
    _request_needs_workspace_project,
    _selected_project_for_command,
    _parse_json_response,
)
from chemsmart.settings.workspace_project import resolve_workspace_project

JsonDict = dict[str, Any]
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


class SynthesisWorkflowMixin:
    """Stateful workflow methods mixed into the public synthesis facade."""

    def prepare_command(self, request: str) -> JsonDict:
        self._reset_preparation_state()
        clarification = self._resolve_workspace_context(request)
        if clarification is not None:
            return clarification
        routed = self._route_frontier_request(request)
        if routed is not None:
            return routed
        result = self.synthesize(request, enforce_intent=False)
        if result["status"] != "ready":
            return self._attach_harness_evidence(result)
        repaired = self._repair_ready_result(request, result)
        if repaired is None:
            return self._failed_preparation(request, result)
        if repaired.get("status") != "ready":
            return self._attach_harness_evidence(repaired)
        command = str(repaired.get("command") or "")
        self._remember_turn(
            request, command, assistant_message=self._last_raw_response
        )
        return self._attach_harness_evidence(repaired)

    def _reset_preparation_state(self) -> None:
        self._last_semantic_result = None
        self._last_intent_assertion = None
        self._last_decision_trace = None
        self._last_reasoning = ""

    def _resolve_workspace_context(self, request: str) -> JsonDict | None:
        requested = select_project_from_request(request)
        if requested and requested.get("selected"):
            self.default_project = str(requested["project"])
        elif not self.default_project:
            selected = current_workflow_state().project
            if selected is not None:
                self.default_project = selected.name
        status = resolve_workspace_project()
        if not self.default_project and status.loaded:
            self.default_project = status.project
            select_workspace_project(status.project, status.program)
        if self.default_project or status.loaded:
            return None
        if not _request_needs_workspace_project(request):
            return None
        return self._attach_decision_trace(
            {
                "status": "needs_clarification",
                "command": "",
                "explanation": (
                    "No workspace project YAML is loaded for this chemsmart "
                    "workspace. Build one with /init, then write it with "
                    "/write-project, or specify a workspace project YAML "
                    "with -p/--project."
                ),
                "confidence": "high",
                "missing_info": [
                    "Create or select a workspace project YAML before command synthesis."
                ],
                "alternatives": [],
                "yaml_check": {"loaded": False, "message": status.message},
            }
        )

    def _failed_preparation(self, request: str, result: JsonDict) -> JsonDict:
        missing = self._semantic_clarification()
        if missing:
            command = str(result.get("command") or "")
            if command:
                self._remember_turn(
                    request,
                    command,
                    assistant_message=self._last_raw_response,
                )
            return self._attach_harness_evidence(
                {
                    "status": "needs_clarification",
                    "command": command,
                    "explanation": (
                        "The command is almost ready but the chemsmart runtime "
                        "needs a bit more information to run it."
                    ),
                    "confidence": "medium",
                    "missing_info": missing,
                    "alternatives": [],
                }
            )
        return self._attach_harness_evidence(
            {
                "status": "infeasible",
                "command": "",
                "explanation": (
                    "Synthesized command failed validation and could not be repaired."
                ),
                "confidence": "low",
                "missing_info": [],
                "alternatives": [],
            }
        )

    def _route_frontier_request(self, request: str) -> JsonDict | None:
        if not self.enable_intent_router:
            return None
        decision = self._frontier_intent_decision(request)
        if decision is None:
            return None
        action = str(decision.get("action") or "synthesize_command")
        target = (
            _extract_command_from_text(request)
            or str(decision.get("target_command") or "").strip()
            or self._last_ready_command()
        )
        self._last_decision_trace = _build_decision_trace(
            decision=decision,
            action=action,
            request=request,
            target_command=target,
            default_project=self.default_project,
            last_command=self._last_ready_command(),
        )
        if action == "synthesize_command":
            return None
        if action in {"explain_command", "critique_command", "repair_command"}:
            return self._route_command_action(
                request, decision, action, target
            )
        if action == "clarify":
            return self._route_clarification(request, decision)
        return None

    def _route_command_action(
        self,
        request: str,
        decision: JsonDict,
        action: str,
        target: str,
    ) -> JsonDict:
        if not target:
            return self._attach_decision_trace(
                {
                    "status": "needs_clarification",
                    "command": "",
                    "explanation": "No previous chemsmart command is available.",
                    "confidence": "medium",
                    "missing_info": [
                        "Paste the chemsmart command you want me to explain or critique."
                    ],
                    "alternatives": [],
                    "action": action,
                }
            )
        target = _ensure_program_project(target, self.default_project)
        if self._last_decision_trace is not None:
            self._last_decision_trace["target_command"] = target
        semantic = None
        if action in {"critique_command", "repair_command"}:
            self._last_semantic_result = self._evaluate_semantics(target)
            semantic = self._last_semantic_result.to_dict()
        composed = compose_command_answer(
            self.provider,
            request=request,
            action=action,
            parsed=parse_model_command(target),
            semantic_summary=semantic,
            timeout_s=self.semantic_timeout_s,
        )
        self._record_answer_reasoning(composed)
        return self._attach_decision_trace(
            {
                "status": "informational",
                "command": target,
                "explanation": composed.answer,
                "confidence": str(decision.get("confidence") or "high"),
                "missing_info": [],
                "alternatives": [],
                "action": action,
                "project": self.default_project,
                "reasoning": list(composed.reasoning),
                "caveats": list(composed.caveats),
                "grounded": composed.grounded,
            }
        )

    def _route_clarification(
        self, request: str, decision: JsonDict
    ) -> JsonDict:
        missing = decision.get("missing_info")
        if not isinstance(missing, list):
            missing = []
        reasoning: list[str] = []
        composed = reason_missing_info(
            self.provider,
            request=request,
            hints=tuple(str(item) for item in missing),
            timeout_s=self.semantic_timeout_s,
        )
        if composed is not None:
            if composed.questions:
                missing = list(composed.questions)
            reasoning = list(composed.reasoning)
            if reasoning and self._last_decision_trace is not None:
                self._last_decision_trace["reasoning"] = reasoning
        if not missing:
            missing = ["What command or calculation should I work with?"]
        return self._attach_decision_trace(
            {
                "status": "needs_clarification",
                "command": "",
                "explanation": str(decision.get("decision_summary") or ""),
                "confidence": str(decision.get("confidence") or "medium"),
                "missing_info": missing,
                "alternatives": [],
                "action": "clarify",
                "reasoning": reasoning,
            }
        )

    def _record_answer_reasoning(self, composed: ComposedAnswer) -> None:
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

    def _attach_harness_evidence(self, result: JsonDict) -> JsonDict:
        normalized = dict(result)
        if self._last_semantic_result is not None:
            normalized.setdefault(
                "semantic", self._last_semantic_result.to_dict()
            )
        if self._last_intent_assertion is not None:
            normalized.setdefault(
                "intent_assertion", self._last_intent_assertion
            )
        return self._attach_decision_trace(normalized)

    def _frontier_intent_decision(self, request: str) -> JsonDict | None:
        messages = [
            {"role": "system", "content": _INTENT_ROUTER_SYSTEM_PROMPT},
            {
                "role": "user",
                "content": json.dumps(
                    {
                        "request": request,
                        "last_command": self._last_ready_command(),
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
        return (
            parsed
            if str(parsed.get("action") or "") in _ROUTER_ACTIONS
            else None
        )

    def _last_ready_command(self) -> str:
        for turn in reversed(self.memory.turns):
            command = str(turn.plan_rationale or "").strip()
            if command.startswith("chemsmart "):
                return command
        return ""

    def _repair_ready_result(
        self, current_request: str, result: JsonDict
    ) -> JsonDict | None:
        semantic_retried = False
        intent_retried = False
        validation_error = ""
        for attempt in range(self.max_retries + 1):
            if result["status"] != "ready":
                click.echo(
                    result.get("explanation")
                    or "Could not synthesize command."
                )
                return None
            result, retry, semantic_retried, intent_retried, error = (
                self._check_ready_result(
                    current_request,
                    result,
                    semantic_retried=semantic_retried,
                    intent_retried=intent_retried,
                )
            )
            if error:
                validation_error = error
            if retry is None:
                return result
            if attempt >= self.max_retries:
                return self._repair_exhausted(result, validation_error)
            result = self.synthesize(
                retry,
                intent_request=current_request,
                enforce_intent=False,
            )
        return None

    def _check_ready_result(
        self,
        request: str,
        result: JsonDict,
        *,
        semantic_retried: bool,
        intent_retried: bool,
    ) -> tuple[JsonDict, str | None, bool, bool, str]:
        command = str(result.get("command") or "")
        valid, error = self.validate_command(command)
        if not valid:
            retry = (
                f"Original request: {request}\nYour command failed validation: "
                f"{error}\nReturn a corrected legal command."
            )
            return result, retry, semantic_retried, intent_retried, error
        if self.semantic_gate:
            self._last_semantic_result = self._evaluate_semantics(command)
            if self._last_semantic_result.verdict == "reject":
                click.echo(
                    "Notice: synthesized command failed runtime semantic "
                    "validation; asking the model to repair it."
                )
                return (
                    result,
                    self._last_semantic_result.correction_prompt(request),
                    True,
                    intent_retried,
                    "",
                )
            if self._last_semantic_result.verdict == "warn":
                click.echo(
                    "Notice: synthesized command passed with runtime semantic warnings."
                )
        checked = self._enforce_request_intent(request, result)
        if checked.get("status") == "ready":
            if semantic_retried or intent_retried:
                click.echo(
                    "Notice: synthesized command was repaired after "
                    "deterministic harness validation."
                )
            return checked, None, semantic_retried, intent_retried, ""
        return (
            checked,
            self._intent_correction_prompt(request, checked),
            semantic_retried,
            True,
            "",
        )

    def _repair_exhausted(
        self, result: JsonDict, validation_error: str
    ) -> JsonDict | None:
        if result.get("intent_reject"):
            return result
        if self.semantic_gate and self._last_semantic_result:
            click.echo(
                "Synthesized command failed runtime semantic validation: "
                + "; ".join(
                    issue.message
                    for issue in self._last_semantic_result.issues
                )
            )
        else:
            click.echo(f"Synthesized command is invalid: {validation_error}")
        return None

    def confirm_and_execute(
        self, cmd: str, explanation: str, confidence: str
    ) -> None:
        text = (
            f"[bold]Command[/bold]\n{cmd}\n\n"
            f"[bold]Explanation[/bold]\n{explanation or '(none)'}\n\n"
            f"[bold]Confidence[/bold] {confidence}"
        )
        Console().print(Panel(text, title="ChemSmart synthesis"))
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
            self._run_subprocess(self._quiet_argv(shlex.split(cmd)))
        elif choice in {"N", "NO"}:
            click.echo("cancelled")
        elif choice in {"E", "EDIT"}:
            self._edit_and_confirm(cmd, explanation, confidence)
        elif choice in {"T", "TEST"}:
            self._test_submission(cmd)
        else:
            click.echo("cancelled")

    def _edit_and_confirm(
        self, cmd: str, explanation: str, confidence: str
    ) -> None:
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

    def _test_submission(self, cmd: str) -> None:
        tokens = shlex.split(cmd)
        if len(tokens) <= 1 or tokens[1] != "sub":
            click.echo("test only for sub")
            return
        tokens.insert(2, "--test")
        self._run_subprocess(self._quiet_argv(tokens))

    def _apply_default_project(self, result: JsonDict) -> JsonDict:
        if not self.default_project:
            return result
        normalized = dict(result)
        command = str(normalized.get("command") or "")
        if normalized.get("status") == "needs_clarification":
            return self._resolve_project_clarification(normalized, command)
        if normalized.get("status") != "ready":
            return result
        normalized["command"] = _ensure_program_project(
            command, self.default_project
        )
        alternatives = normalized.get("alternatives") or []
        if isinstance(alternatives, list):
            normalized["alternatives"] = [
                _ensure_program_project(str(item), self.default_project)
                if isinstance(item, str)
                else item
                for item in alternatives
            ]
        normalized["project"] = self.default_project
        return normalized

    def _resolve_project_clarification(
        self, normalized: JsonDict, command: str
    ) -> JsonDict:
        missing = [str(item) for item in normalized.get("missing_info") or []]
        project_missing = [
            item for item in missing if _is_project_name_missing(item)
        ]
        remaining = [
            item for item in missing if not _is_project_name_missing(item)
        ]
        selected = _selected_project_for_command(command, self.default_project)
        if (
            not project_missing
            or selected is None
            or not _command_project_is_unresolved(command)
        ):
            return normalized
        normalized["command"] = _ensure_program_project(command, selected.name)
        normalized["project"] = selected.name
        normalized["missing_info"] = remaining
        if not remaining:
            normalized["status"] = "ready"
            explanation = str(normalized.get("explanation") or "").strip()
            normalized["explanation"] = (
                f"{explanation} Using the selected workspace project "
                f"{selected.name}."
            ).strip()
        return normalized


__all__ = ["SynthesisWorkflowMixin"]
