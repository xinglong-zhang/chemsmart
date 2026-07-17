"""Interactive natural-language to ChemSmart CLI command synthesis."""

from __future__ import annotations

import json
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
from chemsmart.agent.harness.workflow_state import (
    current_workflow_state,
    select_project_from_request,
    select_workspace_project,
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
from chemsmart.agent.prompts.synthesis import build_synthesis_system_prompt
from chemsmart.agent.provider_adapter import extract_response_text
from chemsmart.agent.schema_prune import prune_schema_for_request
from chemsmart.agent.services.conversation_memory import (
    ConversationMemory,
    ConversationTurn,
)
from chemsmart.agent.services.synthesis_support import (
    _build_decision_trace,
    _clarification_from_semantic,
    _command_project_is_unresolved,
    _ensure_program_project,
    _extract_command_from_text,
    _extract_reasoning,
    _is_frontier_provider,
    _is_project_name_missing,
    _is_v8_spec,
    _normalize_result,
    _normalize_v8_spec,
    _parse_json_response,
    _request_needs_workspace_project,
    _selected_project_for_command,
    _validate_tokens_against_schema,
    resolve_default_project,
)
from chemsmart.settings.workspace_project import resolve_workspace_project

# Structural spec-invariant rule ids that are reliable to enforce at synthesis
# time. The query-heuristic rules (decline_contract / required_present / label)
# are intentionally excluded here to avoid false rejects on well-formed specs.
_RELIABLE_SPEC_REJECTS = frozenset(
    {
        "spec.kind.canonical",
        "spec.runtime_owned.leak",
        "spec.settings.allowed",
        "spec.ts.runtime_freq",
        "spec.ts.runtime_route",
        "spec.ts.bad_extra",
        "spec.jobs.empty",
    }
)


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


def quiet_chemsmart_argv(
    tokens: list[str],
    *,
    debug: bool = False,
) -> list[str]:
    """Inject ``--no-verbose`` after ``chemsmart`` unless debug is on."""

    if debug:
        return list(tokens)
    if not tokens or tokens[0] != "chemsmart":
        return list(tokens)
    if "--no-verbose" in tokens or "--verbose" in tokens:
        return list(tokens)
    argv = list(tokens)
    argv.insert(1, "--no-verbose")
    return argv


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
        self._last_reasoning = ""
        self._last_semantic_result: CommandSemanticResult | None = None
        self._last_intent_assertion: JsonDict | None = None
        self._last_decision_trace: JsonDict | None = None

    def synthesize(
        self,
        request: str,
        *,
        intent_request: str | None = None,
        enforce_intent: bool = True,
    ) -> JsonDict:
        """Return a structured synthesis result for ``request``.

        Invalid JSON is retried with a corrective message up to
        ``max_retries`` times.
        """

        messages = self._messages_for_request(request)
        last_error: Exception | None = None
        for attempt in range(self.max_retries + 1):
            response = self.provider.chat(messages)
            self._last_raw_response = extract_response_text(response)
            self._last_reasoning = _extract_reasoning(response)
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
                    result = self._apply_default_project(result)
                else:
                    result = self._apply_default_project(
                        _normalize_result(parsed)
                    )
                if enforce_intent:
                    return self._enforce_request_intent(
                        intent_request or request,
                        result,
                    )
                return result
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

    def _enforce_request_intent(
        self, request: str, result: JsonDict
    ) -> JsonDict:
        """Downgrade a ``ready`` command that does not preserve the request's
        intent (e.g. a coordinate scan silently collapsed to a plain ``opt``)
        to ``needs_clarification`` so the clarify loop repairs it. This mirrors
        the deterministic intent gate already applied on the tool-calling path
        (:mod:`chemsmart.agent.tools_command`) but was missing on the direct
        command-string synthesis path, letting silent-wrong commands surface as
        ready. The gate only asserts fields the request states explicitly, so
        under-specified requests are never over-rejected."""
        if result.get("status") != "ready":
            return result
        command = str(result.get("command") or "").strip()
        if not command:
            return result
        try:
            from chemsmart.agent.harness.intent import (
                IntentSpec,
                evaluate_intent,
            )

            intent = evaluate_intent(command, IntentSpec.from_request(request))
        except (
            Exception
        ):  # pragma: no cover - never fail synthesis on the gate
            return result
        self._last_intent_assertion = intent.to_dict()
        if intent.verdict != "reject":
            return result
        failed = list(intent.failed_rule_ids)
        result["intent_assertion"] = self._last_intent_assertion
        result["intent_reject"] = failed
        result["status"] = "needs_clarification"
        result["missing_info"] = failed or ["intent"]
        result["explanation"] = (
            (result.get("explanation") or "")
            + " The drafted command runs but does not preserve the requested "
            "intent ("
            + ", ".join(failed)
            + "); confirm the intended job before running."
        ).strip()
        return result

    def _intent_correction_prompt(
        self,
        original_request: str,
        result: JsonDict,
    ) -> str:
        assertion = self._last_intent_assertion or {}
        failed_rows = [
            row
            for row in assertion.get("assertions") or []
            if isinstance(row, dict) and row.get("status") == "fail"
        ]
        evidence = "\n".join(
            "- {id}: expected {expected!r}, observed {observed!r}".format(
                id=row.get("id"),
                expected=row.get("expected"),
                observed=row.get("observed"),
            )
            for row in failed_rows
        )
        return (
            f"Original request: {original_request}\n"
            "The command passed safe runtime validation but changed an "
            "explicitly stated user value or job intent.\n"
            f"Command: {result.get('command') or ''}\n"
            f"Intent mismatches:\n{evidence or '- unknown'}\n"
            "Return ONLY corrected JSON matching the synthesis schema. "
            "Preserve the explicit expected values above; do not invent "
            "any missing scientific setting."
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
        self._last_intent_assertion = None
        self._last_decision_trace = None
        self._last_reasoning = ""

        requested_project = select_project_from_request(request)
        if requested_project and requested_project.get("selected"):
            self.default_project = str(requested_project["project"])
        elif not self.default_project:
            selected = current_workflow_state().project
            if selected is not None:
                self.default_project = selected.name

        yaml_status = resolve_workspace_project()
        if not self.default_project and yaml_status.loaded:
            # A single workspace YAML is an unambiguous runtime default.  Keep
            # the command synthesizer and the observable workflow state in
            # sync so providers do not ask again for an already loaded project.
            self.default_project = yaml_status.project
            select_workspace_project(yaml_status.project, yaml_status.program)
        if (
            not self.default_project
            and not yaml_status.loaded
            and _request_needs_workspace_project(request)
        ):
            return self._attach_decision_trace(
                {
                    "status": "needs_clarification",
                    "command": "",
                    "explanation": (
                        "No workspace project YAML is loaded for this "
                        "chemsmart workspace. Build one with /init, then "
                        "write it with /write-project, or specify a workspace "
                        "project YAML with -p/--project."
                    ),
                    "confidence": "high",
                    "missing_info": [
                        "Create or select a workspace project YAML before command synthesis."
                    ],
                    "alternatives": [],
                    "yaml_check": {
                        "loaded": False,
                        "message": yaml_status.message,
                    },
                }
            )

        routed = self._route_frontier_request(request)
        if routed is not None:
            return routed

        # Runtime semantics and generated-input invariants must run before the
        # request-intent gate.  Otherwise a lossy command parser can reject a
        # valid command before the stronger runtime evidence exists, and a
        # repairable CLI failure never reaches the repair loop.
        result = self.synthesize(request, enforce_intent=False)
        if result["status"] != "ready":
            return self._attach_harness_evidence(result)

        repaired = self._repair_ready_result(request, result)
        if repaired is None:
            # If the runtime rejected the command because a required input
            # (e.g. charge/multiplicity) is missing, ask the user for it instead
            # of dead-ending on "could not be repaired".
            missing = _clarification_from_semantic(self._last_semantic_result)
            if missing:
                failed_command = str(result.get("command") or "")
                # Remember the near-ready turn so the user's short answer
                # (e.g. "charge 0 multiplicity 1") is merged with context.
                if failed_command:
                    self._remember_turn(
                        request,
                        failed_command,
                        assistant_message=self._last_raw_response,
                    )
                return self._attach_harness_evidence(
                    {
                        "status": "needs_clarification",
                        "command": failed_command,
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
                        "Synthesized command failed validation and could not be "
                        "repaired."
                    ),
                    "confidence": "low",
                    "missing_info": [],
                    "alternatives": [],
                }
            )
        if repaired.get("status") != "ready":
            return self._attach_harness_evidence(repaired)
        command = str(repaired.get("command") or "")
        self._remember_turn(
            request,
            command,
            assistant_message=self._last_raw_response,
        )
        return self._attach_harness_evidence(repaired)

    def _route_frontier_request(self, request: str) -> JsonDict | None:
        if not self.enable_intent_router:
            return None
        decision = self._frontier_intent_decision(request)
        if decision is None:
            return None
        action = str(decision.get("action") or "synthesize_command")
        target_command = (
            _extract_command_from_text(request)
            or str(decision.get("target_command") or "").strip()
            or self._last_ready_command()
        )
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
            return self._attach_decision_trace(
                {
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
                }
            )
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
            return self._attach_decision_trace(
                {
                    "status": "needs_clarification",
                    "command": "",
                    "explanation": summary,
                    "confidence": str(decision.get("confidence") or "medium"),
                    "missing_info": missing,
                    "alternatives": [],
                    "action": action,
                    "reasoning": reasoning,
                }
            )
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

    def _attach_harness_evidence(self, result: JsonDict) -> JsonDict:
        """Expose deterministic evidence without overwriting model fields."""

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
        intent_retried = False
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
                if semantic_result.verdict != "reject":
                    checked = self._enforce_request_intent(
                        current_request,
                        result,
                    )
                    if checked.get("status") != "ready":
                        result = checked
                        intent_retried = True
                        retry_request = self._intent_correction_prompt(
                            current_request,
                            result,
                        )
                    elif semantic_retried or intent_retried:
                        click.echo(
                            "Notice: synthesized command was repaired after "
                            "deterministic harness validation."
                        )
                        return checked
                    else:
                        return checked
            else:
                checked = self._enforce_request_intent(
                    current_request,
                    result,
                )
                if checked.get("status") == "ready":
                    return checked
                result = checked
                intent_retried = True
                retry_request = self._intent_correction_prompt(
                    current_request,
                    result,
                )

            if retry_request is None:
                return result
            if attempt >= self.max_retries:
                if result.get("intent_reject"):
                    return result
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
            # The repair prompt contains diagnostic text such as server names
            # and filenames. Preserve the original user request as the only
            # source of intent assertions so diagnostic context cannot become
            # a false user requirement.
            result = self.synthesize(
                retry_request,
                intent_request=current_request,
                enforce_intent=False,
            )
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
            subprocess.run(
                quiet_chemsmart_argv(shlex.split(cmd), debug=self.debug),
                check=False,
            )
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
                subprocess.run(
                    quiet_chemsmart_argv(tokens, debug=self.debug),
                    check=False,
                )
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
        return quiet_chemsmart_argv(tokens, debug=self.debug)

    def _messages_for_request(self, request: str) -> list[dict[str, str]]:
        selected_project = current_workflow_state().project
        workspace_program = (
            selected_project.program if selected_project is not None else None
        )
        prompt_schema = prune_schema_for_request(
            self.schema,
            request,
            workspace_program=workspace_program,
        )
        system_prompt = build_synthesis_system_prompt(prompt_schema)
        if (
            selected_project is not None
            and selected_project.name == self.default_project
            and Path(selected_project.path).is_file()
        ):
            system_prompt += (
                "\nRuntime-owned workspace context: the selected "
                f"{selected_project.program} project is "
                f"`{selected_project.name}`. For a matching program command, "
                "use this value for -p/--project and do not ask the user for "
                "the project name again."
            )
        messages = [
            {
                "role": "system",
                "content": system_prompt,
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

        if not self.default_project:
            return result
        normalized = dict(result)
        command = str(normalized.get("command") or "")

        if normalized.get("status") == "needs_clarification":
            missing = [
                str(item) for item in normalized.get("missing_info") or []
            ]
            project_missing = [
                item for item in missing if _is_project_name_missing(item)
            ]
            remaining = [
                item for item in missing if not _is_project_name_missing(item)
            ]
            selected = _selected_project_for_command(
                command,
                self.default_project,
            )
            if (
                not project_missing
                or selected is None
                or not _command_project_is_unresolved(command)
            ):
                return result
            normalized["command"] = _ensure_program_project(
                command,
                selected.name,
            )
            normalized["project"] = selected.name
            normalized["missing_info"] = remaining
            if not remaining:
                normalized["status"] = "ready"
                explanation = str(normalized.get("explanation") or "").strip()
                note = (
                    f"Using the selected workspace project {selected.name}."
                )
                normalized["explanation"] = (
                    f"{explanation} {note}".strip()
                )
            return normalized

        if normalized.get("status") != "ready":
            return result
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
