"""Interactive natural-language to ChemSmart CLI command synthesis."""

from __future__ import annotations

import json
import shlex
import subprocess
from pathlib import Path
from typing import Any

import click

from chemsmart.agent.cli_schema import build_chemsmart_cli_schema
from chemsmart.agent.harness.command_semantics import (
    CommandSemanticResult,
    evaluate_command_semantics,
)
from chemsmart.agent.harness.spec_invariants import check_spec
from chemsmart.agent.harness.workflow_state import (
    current_workflow_state,
)
from chemsmart.agent.kind_disambiguator import disambiguate
from chemsmart.agent.prompts.synthesis import build_synthesis_system_prompt
from chemsmart.agent.provider_adapter import extract_response_text
from chemsmart.agent.schema_prune import prune_schema_for_request
from chemsmart.agent.services.conversation_memory import (
    ConversationMemory,
    ConversationTurn,
)
from chemsmart.agent.services.synthesis_support import (
    _clarification_from_semantic,
    _extract_reasoning,
    _is_frontier_provider,
    _is_v8_spec,
    _normalize_result,
    _normalize_v8_spec,
    _parse_json_response,
    _validate_tokens_against_schema,
    resolve_default_project,
)
from chemsmart.agent.services.synthesis_workflow import SynthesisWorkflowMixin

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


class SynthesisSession(SynthesisWorkflowMixin):
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

    def _evaluate_semantics(self, command: str) -> CommandSemanticResult:
        """Evaluate through the public facade's patch-compatible boundary."""

        return evaluate_command_semantics(
            command,
            timeout_s=self.semantic_timeout_s,
        )

    def _semantic_clarification(self) -> list[str]:
        """Translate the latest semantic failure into missing information."""

        return _clarification_from_semantic(self._last_semantic_result)

    def _run_subprocess(self, tokens: list[str]) -> None:
        """Run a confirmed command through the public patch boundary."""

        subprocess.run(tokens, check=False)

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
