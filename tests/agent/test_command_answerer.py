"""Unit tests for the tool-grounded command answerer."""

from __future__ import annotations

import json
from typing import Any

from chemsmart.agent.command_answerer import (
    build_command_facts,
    compose_command_answer,
    reason_missing_info,
)
from chemsmart.agent.model_command_parser import parse_model_command

_COMMAND = "chemsmart run gaussian -f examples/h2o.xyz -c 0 -m 1 -x B3LYP opt"


class _FakeProvider:
    """Provider double returning a queued response and recording calls."""

    def __init__(self, response: Any, *, name: str = "openai") -> None:
        self.name = name
        self._response = response
        self.calls: list[list[dict[str, str]]] = []

    def chat(self, messages: list[dict[str, str]], **_kwargs: Any) -> Any:
        self.calls.append(messages)
        if isinstance(self._response, Exception):
            raise self._response
        return self._response


class _RaisingTimeoutProvider:
    """Rejects the timeout_s kwarg to exercise the TypeError fallback path."""

    name = "openai"

    def __init__(self, response: Any) -> None:
        self._response = response
        self.calls = 0

    def chat(self, messages: list[dict[str, str]]) -> Any:  # no timeout_s
        self.calls += 1
        return self._response


def _parsed():
    return parse_model_command(_COMMAND)


def test_build_command_facts_includes_hard_fields() -> None:
    facts = build_command_facts(_parsed())
    assert facts["program"] == "gaussian"
    assert facts["job"] == "opt"
    assert facts["action"] == "run"
    assert facts["filename"] == "examples/h2o.xyz"
    assert facts["charge"] == "0"
    assert facts["multiplicity"] == "1"
    # None values (no project) are omitted, keeping the prompt compact.
    assert "project" not in facts


def test_compose_uses_grounded_model_answer() -> None:
    response = json.dumps(
        {
            "answer": "该命令用 gaussian 对 examples/h2o.xyz 做几何优化。",
            "facts_used": {
                "program": "gaussian",
                "job": "opt",
                "action": "run",
                "project": None,
                "filename": "examples/h2o.xyz",
                "charge": "0",
                "multiplicity": "1",
            },
            "reasoning": ["用户要求中文解释。"],
            "caveats": [],
        }
    )
    provider = _FakeProvider(response)

    composed = compose_command_answer(
        provider,
        request="이 명령 설명해줘. 중국어로.",
        action="explain_command",
        parsed=_parsed(),
    )

    assert composed.grounded is True
    assert composed.fallback_used is False
    assert composed.answer.startswith("该命令")
    assert composed.reasoning == ("用户要求中文解释。",)


def test_compose_falls_back_on_program_contradiction() -> None:
    response = json.dumps(
        {
            "answer": "This runs an orca calculation.",
            "facts_used": {"program": "orca", "job": "opt"},
            "reasoning": [],
            "caveats": [],
        }
    )
    provider = _FakeProvider(response)

    composed = compose_command_answer(
        provider,
        request="explain this",
        action="explain_command",
        parsed=_parsed(),
    )

    assert composed.grounded is False
    assert composed.fallback_used is True
    assert "grounding violation" in composed.fallback_reason
    # Deterministic fallback text is authoritative and English.
    assert "deterministic command parser:" in composed.answer


def test_compose_falls_back_on_prose_program_leak() -> None:
    # facts_used is clean but the prose names the wrong program only.
    response = json.dumps(
        {
            "answer": "This command runs an ORCA single point.",
            "facts_used": {"program": "gaussian"},
            "reasoning": [],
            "caveats": [],
        }
    )
    provider = _FakeProvider(response)

    composed = compose_command_answer(
        provider,
        request="explain this",
        action="explain_command",
        parsed=_parsed(),
    )

    assert composed.grounded is False
    assert "different program" in composed.fallback_reason


def test_compose_falls_back_on_provider_error() -> None:
    provider = _FakeProvider(RuntimeError("boom"))

    composed = compose_command_answer(
        provider,
        request="explain this",
        action="explain_command",
        parsed=_parsed(),
    )

    assert composed.fallback_used is True
    assert "provider error" in composed.fallback_reason


def test_compose_falls_back_on_malformed_json() -> None:
    provider = _FakeProvider("not json at all")

    composed = compose_command_answer(
        provider,
        request="explain this",
        action="explain_command",
        parsed=_parsed(),
    )

    assert composed.fallback_used is True
    assert "malformed response" in composed.fallback_reason


def test_compose_falls_back_on_empty_answer() -> None:
    provider = _FakeProvider(json.dumps({"answer": "", "facts_used": {}}))

    composed = compose_command_answer(
        provider,
        request="explain this",
        action="explain_command",
        parsed=_parsed(),
    )

    assert composed.fallback_used is True
    assert composed.fallback_reason == "empty answer"


def test_compose_handles_provider_without_timeout_kwarg() -> None:
    response = json.dumps(
        {
            "answer": "Runs a gaussian optimization on examples/h2o.xyz.",
            "facts_used": {"program": "gaussian", "job": "opt"},
        }
    )
    provider = _RaisingTimeoutProvider(response)

    composed = compose_command_answer(
        provider,
        request="explain this",
        action="explain_command",
        parsed=_parsed(),
    )

    assert composed.grounded is True
    assert provider.calls == 1


def test_reason_missing_info_returns_questions() -> None:
    response = json.dumps(
        {
            "present": ["file examples/h2o.xyz"],
            "missing": ["job type"],
            "questions": ["Which calculation do you want to run?"],
            "reasoning": ["No job keyword present in the request."],
        }
    )
    provider = _FakeProvider(response)

    cot = reason_missing_info(
        provider, request="do something with examples/h2o.xyz"
    )

    assert cot is not None
    assert cot.questions == ("Which calculation do you want to run?",)
    assert cot.missing == ("job type",)
    assert cot.reasoning == ("No job keyword present in the request.",)


def test_reason_missing_info_returns_none_on_malformed() -> None:
    provider = _FakeProvider("garbage")
    assert reason_missing_info(provider, request="hi") is None


def test_reason_missing_info_returns_none_when_empty() -> None:
    provider = _FakeProvider(json.dumps({"questions": [], "missing": []}))
    assert reason_missing_info(provider, request="hi") is None
