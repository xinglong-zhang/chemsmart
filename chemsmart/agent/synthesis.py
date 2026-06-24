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
        """

        if provider is None:
            from chemsmart.agent import providers

            provider = providers.get_provider()
        self.provider = provider
        self.schema = schema or build_chemsmart_cli_schema()
        self.max_retries = max_retries
        self.debug = debug
        self.memory = ConversationMemory()
        # The model's verbatim output from the most recent synthesize() call, recorded as the
        # assistant turn for native multi-turn replay (see _messages_for_request). Kept off the
        # result dict so the public synthesize() return shape is unchanged.
        self._last_raw_response = ""

    def synthesize(self, request: str) -> JsonDict:
        """Return a structured synthesis result for ``request``.

        Invalid JSON is retried with a corrective message up to
        ``max_retries`` times.
        """

        messages = self._messages_for_request(request)
        last_error: Exception | None = None
        for attempt in range(self.max_retries + 1):
            response = self.provider.chat(messages)
            # Record the model's verbatim output for native multi-turn replay; kept on the
            # session (not the result dict) so synthesize()'s public shape is unchanged.
            self._last_raw_response = _extract_text(response)
            try:
                parsed = _parse_json_response(response)
                if _is_v8_spec(parsed):
                    return _normalize_v8_spec(parsed)
                return _normalize_result(parsed)
            except (TypeError, ValueError, json.JSONDecodeError) as exc:
                last_error = exc
                if attempt >= self.max_retries:
                    break
                messages.append(
                    {
                        "role": "assistant",
                        "content": _extract_text(response),
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

        command = str(result.get("command") or "")
        valid, error = self.validate_command(command)
        if not valid:
            retry_request = (
                f"Original request: {current_request}\n"
                f"Your command failed validation: {error}\n"
                "Return a corrected legal command."
            )
            result = self.synthesize(retry_request)
            if result["status"] != "ready":
                click.echo(
                    result.get("explanation")
                    or "Could not synthesize command."
                )
                return
            command = str(result.get("command") or "")
            valid, error = self.validate_command(command)
            if not valid:
                click.echo(f"Synthesized command is invalid: {error}")
                return

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
        """Build a NATIVE multi-turn message list: the system prompt, then each prior
        turn replayed as real ``user``/``assistant`` messages, then the current request.

        The assistant turns carry the model's verbatim SPEC (``assistant_message``), so a
        follow-up like "change the charge to -1" lets the model edit its own prior spec in
        context — the format it was effectively prompted with. Turns without an assistant
        message (e.g. clarification sub-turns, whose content is folded into a later request)
        are skipped to avoid duplicate/garbled user messages.
        """
        messages: list[dict[str, str]] = [
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
        self, request: str, result: str, assistant_message: str = ""
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


_V8_INTENTS = {"workflow", "advisory", "decline", "chitchat"}


def _is_v8_spec(parsed: JsonDict) -> bool:
    """A v8 spec-emission model returns {"intent": ..., "jobs"/"message": ...} instead of the
    {status, command, ...} synthesis shape."""
    return (
        isinstance(parsed, dict)
        and parsed.get("intent") in _V8_INTENTS
        and ("jobs" in parsed or "message" in parsed)
        and "status" not in parsed
    )


def _normalize_v8_spec(parsed: JsonDict) -> JsonDict:
    """Bridge a v8 job spec into the synthesis result shape via the deterministic v8 adapter
    (parse -> postprocess -> render chemsmart command -> validate)."""
    from chemsmart.agent import v8_adapter

    out = v8_adapter.adapt(parsed)
    if out["intent"] == "workflow":
        commands = out.get("commands") or []
        command = commands[0] if commands else ""
        # a chain renders as multiple commands; the lead command is executed, the rest are surfaced
        return _normalize_result(
            {
                "status": "ready" if command else "infeasible",
                "confidence": "high" if out.get("valid") else "medium",
                "command": command,
                "explanation": ("multi-step: " + " ; ".join(commands)) if len(commands) > 1 else "",
                "missing_info": [] if command else ["adapter produced no command"],
                "alternatives": commands[1:],
            }
        )
    status = "infeasible" if out["intent"] == "decline" else "needs_clarification"
    return _normalize_result(
        {
            "status": status,
            "confidence": "low",
            "command": "",
            "explanation": out.get("message") or "",
            "missing_info": [],
            "alternatives": [],
        }
    )


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
