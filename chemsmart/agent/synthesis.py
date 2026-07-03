"""Interactive natural-language to ChemSmart CLI command synthesis."""

from __future__ import annotations

import json
import re
import shlex
import subprocess
from typing import Any

import click
from rich.console import Console
from rich.panel import Panel

from chemsmart.agent.cli_schema import build_chemsmart_cli_schema
from chemsmart.agent.harness.command_semantics import (
    CommandSemanticResult,
    evaluate_command_semantics,
)
from chemsmart.agent.harness.spec_invariants import check_spec
from chemsmart.agent.kind_disambiguator import disambiguate

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
        """

        if provider is None:
            from chemsmart.agent import providers

            provider = providers.get_provider()
        self.provider = provider
        self.schema = schema or build_chemsmart_cli_schema()
        self.max_retries = max_retries
        self.debug = debug
        self.semantic_gate = semantic_gate
        self.semantic_timeout_s = semantic_timeout_s
        self.memory = ConversationMemory()
        self._last_raw_response = ""
        self._last_semantic_result: CommandSemanticResult | None = None

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
                    result = _normalize_v8_spec(parsed)
                    _apply_spec_invariants(result, parsed, request)
                    return result
                return _normalize_result(parsed)
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

        result = self.synthesize(request)
        if result["status"] != "ready":
            return result

        repaired = self._repair_ready_result(request, result)
        if repaired is None:
            return {
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
        command = str(repaired.get("command") or "")
        self._remember_turn(
            request,
            command,
            assistant_message=self._last_raw_response,
        )
        return repaired

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


def _normalize_v8_spec(result: JsonDict) -> JsonDict:
    from chemsmart.agent.v8_adapter import adapt

    adapted = adapt(result, validate=True)
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
