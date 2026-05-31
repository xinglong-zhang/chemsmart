"""Tests for natural-language ChemSmart CLI command synthesis."""

from __future__ import annotations

import json
from typing import Any

import pytest
from click.testing import CliRunner

from chemsmart.agent.cli import agent
from chemsmart.agent.synthesis import SynthesisSession

READY = {
    "status": "ready",
    "command": "chemsmart sub gaussian ts -p oxetane -b def2-svp",
    "explanation": "Prepare a Gaussian transition-state job.",
    "confidence": "high",
    "missing_info": [],
    "alternatives": [],
}


class DummyProvider:
    """Deterministic provider double returning queued responses."""

    def __init__(self, responses: list[Any]) -> None:
        self.responses = list(responses)
        self.messages: list[list[dict[str, str]]] = []

    def chat(self, messages: list[dict[str, str]], **_kwargs: Any) -> Any:
        self.messages.append(messages)
        if not self.responses:
            raise AssertionError("provider response queue exhausted")
        return self.responses.pop(0)


def _json_response(payload: dict[str, Any]) -> str:
    return json.dumps(payload)


def test_synthesize_parses_happy_path_json() -> None:
    provider = DummyProvider([_json_response(READY)])
    session = SynthesisSession(provider=provider, schema={"subcommands": {}})

    result = session.synthesize("make a ts job")

    assert result == READY
    assert provider.messages[0][-1] == {
        "role": "user",
        "content": "make a ts job",
    }


def test_synthesize_retries_once_after_invalid_json() -> None:
    provider = DummyProvider(["not json", _json_response(READY)])
    session = SynthesisSession(
        provider=provider,
        schema={"subcommands": {}},
        max_retries=1,
    )

    result = session.synthesize("make a ts job")

    assert result["status"] == "ready"
    assert len(provider.messages) == 2
    assert (
        "previous response was invalid" in provider.messages[1][-1]["content"]
    )


def test_synthesize_raises_after_two_invalid_json_responses() -> None:
    provider = DummyProvider(["not json", "still not json"])
    session = SynthesisSession(
        provider=provider,
        schema={"subcommands": {}},
        max_retries=1,
    )

    with pytest.raises(ValueError, match="failed to return valid"):
        session.synthesize("make a ts job")

    assert len(provider.messages) == 2


def test_validate_command_accepts_legal_sub_gaussian_ts_command() -> None:
    session = SynthesisSession(provider=DummyProvider([]))

    ok, error = session.validate_command(
        "chemsmart sub gaussian ts -p oxetane -b def2-svp"
    )

    assert ok is True
    assert error == ""


def test_validate_command_rejects_invented_option() -> None:
    session = SynthesisSession(provider=DummyProvider([]))

    ok, error = session.validate_command("chemsmart sub gaussian ts -z bogus")

    assert ok is False
    assert "-z" in error


def test_validate_command_rejects_non_chemsmart_command() -> None:
    session = SynthesisSession(provider=DummyProvider([]))

    ok, error = session.validate_command("python -m chemsmart")

    assert ok is False
    assert "chemsmart" in error


def test_confirm_and_execute_inserts_test_only_for_sub(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = SynthesisSession(provider=DummyProvider([]))
    calls: list[list[str]] = []
    choices = iter(["T", "T"])

    monkeypatch.setattr(
        "click.prompt", lambda *_args, **_kwargs: next(choices)
    )
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.subprocess.run",
        lambda args, check=False: calls.append(args),
    )

    session.confirm_and_execute(
        "chemsmart sub gaussian ts -p oxetane",
        "explain",
        "medium",
    )
    session.confirm_and_execute("chemsmart config show", "explain", "medium")

    assert calls == [
        [
            "chemsmart",
            "--no-verbose",
            "sub",
            "--test",
            "gaussian",
            "ts",
            "-p",
            "oxetane",
        ]
    ]


def test_confirm_and_execute_injects_no_verbose_on_yes_by_default(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = SynthesisSession(provider=DummyProvider([]))
    calls: list[list[str]] = []

    monkeypatch.setattr("click.prompt", lambda *_args, **_kwargs: "Y")
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.subprocess.run",
        lambda args, check=False: calls.append(args),
    )

    session.confirm_and_execute(
        "chemsmart sub gaussian opt -p water -b 6-31g*",
        "explain",
        "high",
    )

    assert calls == [
        [
            "chemsmart",
            "--no-verbose",
            "sub",
            "gaussian",
            "opt",
            "-p",
            "water",
            "-b",
            "6-31g*",
        ]
    ]


def test_confirm_and_execute_skips_no_verbose_when_debug(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = SynthesisSession(provider=DummyProvider([]), debug=True)
    calls: list[list[str]] = []

    monkeypatch.setattr("click.prompt", lambda *_args, **_kwargs: "Y")
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.subprocess.run",
        lambda args, check=False: calls.append(args),
    )

    session.confirm_and_execute(
        "chemsmart sub gaussian opt -p water -b 6-31g*",
        "explain",
        "high",
    )

    assert calls == [
        [
            "chemsmart",
            "sub",
            "gaussian",
            "opt",
            "-p",
            "water",
            "-b",
            "6-31g*",
        ]
    ]


def test_quiet_argv_respects_existing_verbose_flags() -> None:
    session = SynthesisSession(provider=DummyProvider([]))
    argv_with_verbose = session._quiet_argv(
        ["chemsmart", "--verbose", "sub", "gaussian", "opt"]
    )
    argv_with_no_verbose = session._quiet_argv(
        ["chemsmart", "--no-verbose", "sub", "gaussian", "opt"]
    )

    assert argv_with_verbose == [
        "chemsmart",
        "--verbose",
        "sub",
        "gaussian",
        "opt",
    ]
    assert argv_with_no_verbose == [
        "chemsmart",
        "--no-verbose",
        "sub",
        "gaussian",
        "opt",
    ]


def test_multiturn_clarification_is_carried_in_memory(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    clarification = {
        "status": "needs_clarification",
        "command": "",
        "explanation": "Need method.",
        "confidence": "medium",
        "missing_info": ["What basis set?"],
        "alternatives": [],
    }
    provider = DummyProvider(
        [_json_response(clarification), _json_response(READY)]
    )
    session = SynthesisSession(provider=provider)

    monkeypatch.setattr("click.prompt", lambda *_args, **_kwargs: "def2-svp")
    monkeypatch.setattr(session, "confirm_and_execute", lambda *_args: None)

    session.run_interactive("optimize oxetane")

    assert len(provider.messages) == 2
    memory_message = provider.messages[1][1]["content"]
    assert "Conversation memory" in memory_message
    assert "What basis set?: def2-svp" in memory_message


def test_agent_ask_e2e_uses_synthesis_session(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    provider = DummyProvider([_json_response(READY)])
    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider",
        lambda: provider,
    )
    calls: list[list[str]] = []
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.subprocess.run",
        lambda args, check=False: calls.append(args),
    )

    result = CliRunner().invoke(agent, ["ask", "make a ts job"], input="N\n")

    assert result.exit_code == 0
    assert "chemsmart sub gaussian ts -p oxetane -b def2-svp" in result.output
    assert "Run? [Y]es/[E]dit/[N]o/[T]est" in result.output
    assert calls == []


def test_agent_ask_debug_flag_plumbs_to_synthesis_session(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    provider = DummyProvider([_json_response(READY)])
    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider",
        lambda: provider,
    )
    captured: dict[str, bool] = {}
    real_init = SynthesisSession.__init__

    def spy_init(self, *args: Any, **kwargs: Any) -> None:
        captured["debug"] = bool(kwargs.get("debug", False))
        real_init(self, *args, **kwargs)

    monkeypatch.setattr(SynthesisSession, "__init__", spy_init)
    calls: list[list[str]] = []
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.subprocess.run",
        lambda args, check=False: calls.append(args),
    )

    result = CliRunner().invoke(
        agent, ["ask", "--debug", "make a ts job"], input="Y\n"
    )

    assert result.exit_code == 0, result.output
    assert captured.get("debug") is True
    assert calls == [
        [
            "chemsmart",
            "sub",
            "gaussian",
            "ts",
            "-p",
            "oxetane",
            "-b",
            "def2-svp",
        ]
    ]
