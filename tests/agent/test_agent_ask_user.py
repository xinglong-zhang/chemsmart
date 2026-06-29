from __future__ import annotations

import json
from typing import Any

from chemsmart.agent.core import AgentSession
from chemsmart.agent.loop import (
    ASK_USER_TOOL_DEF,
    ASK_USER_TOOL_NAME,
    ToolLoop,
)
from chemsmart.agent.permissions import RuntimePermissionMode
from chemsmart.agent.prompts.identity import build_system_prompt
from chemsmart.agent.registry import ToolRegistry

from ._agent_session_helpers import FakeProvider
from ._loop_helpers import (
    ScriptedRegistry,
    openai_final_response,
    openai_tool_call_response,
    tool_call,
)


class RecordingDecisionLog:
    def __init__(self) -> None:
        self.entries: list[dict[str, Any]] = []

    def write(
        self,
        kind: str,
        payload: dict[str, Any],
        rationale: str = "",
    ) -> None:
        self.entries.append(
            {
                "kind": kind,
                "payload": payload,
                "rationale": rationale,
            }
        )


class ExplodingPolicy:
    def __init__(self) -> None:
        self.calls = 0

    def resolve(self, request) -> Any:  # pragma: no cover - should not run
        self.calls += 1
        raise AssertionError(
            f"resolve should not be called for {request.name}"
        )


def test_ask_user_tool_def_is_available_for_provider_and_read_only_mode(
    tmp_path,
):
    loop = ToolLoop(
        provider=FakeProvider([]),
        registry=ToolRegistry.default(),
        handle_store=None,
        decision_log=RecordingDecisionLog(),
    )

    provider_defs = loop._tool_defs_for_provider("openai")
    read_only_defs = loop._tool_defs_for_mode(
        "openai",
        RuntimePermissionMode.READ_ONLY,
    )

    assert ASK_USER_TOOL_DEF in provider_defs
    assert ASK_USER_TOOL_DEF in read_only_defs
    assert [
        tool_def["function"]["name"]
        for tool_def in provider_defs
        if tool_def.get("type") == "function"
    ].count(ASK_USER_TOOL_NAME) == 1


def test_run_turn_returns_ask_user_and_skips_permission_resolution(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call(
                        "call_ask",
                        ASK_USER_TOOL_NAME,
                        {
                            "question": "which server?",
                            "options": ["chemnode1", "chemnode2", 7],
                        },
                    ),
                    tool_call(
                        "call_other",
                        "recommend_method",
                        {"task": "opt"},
                    ),
                )
            }
        ]
    )
    registry = ScriptedRegistry({"recommend_method": {"method": "b3lyp"}})
    decision_log = RecordingDecisionLog()
    policy = ExplodingPolicy()
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=None,
        decision_log=decision_log,
        policy=policy,
    )

    result = loop.run_turn(messages=[{"role": "user", "content": "Check."}])

    assert result["ask_user"] == {
        "question": "which server?",
        "options": ["chemnode1", "chemnode2"],
    }
    assert len(result["tool_requests"]) == 1
    assert len(result["tool_outcomes"]) == 1
    assert result["tool_outcomes"][0].status == "ask_user"
    assert policy.calls == 0
    assert registry.calls == []
    assert any(entry["kind"] == "ask_user" for entry in decision_log.entries)


def test_run_loop_surfaces_ask_user_question_and_keeps_turn_in_progress(
    tmp_path,
):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call(
                        "call_ask",
                        ASK_USER_TOOL_NAME,
                        {
                            "question": "which server?",
                            "options": ["chemnode1", "chemnode2"],
                        },
                    )
                )
            }
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ScriptedRegistry({"recommend_method": {"method": "b3lyp"}}),
        session_root=tmp_path,
    )

    result = session.run_loop("Check the queue.")

    assert result["ask_user_question"] == {
        "question": "which server?",
        "options": ["chemnode1", "chemnode2"],
    }
    assert result["assistant_output"] == ""
    assert result["completed_steps"] == 0
    assert session.conversation_history.turns[-1].status == "in_progress"
    assert session.state is not None
    assert session.state.pending_ask_user == result["ask_user_question"]
    assert session.state.pending_messages is not None


def test_run_loop_answer_continues_pending_ask_user_turn(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call(
                        "call_ask",
                        ASK_USER_TOOL_NAME,
                        {
                            "question": "which server?",
                            "options": ["chemnode1", "chemnode2"],
                        },
                    )
                )
            },
            {
                "__raw_response__": openai_final_response(
                    "Got it — I will check chemnode1."
                )
            },
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ScriptedRegistry({"recommend_method": {"method": "b3lyp"}}),
        session_root=tmp_path,
    )

    first = session.run_loop("Check the queue.")
    second = session.run_loop("chemnode1")

    assert first["ask_user_question"] is not None
    assert second["ask_user_question"] is None
    assert second["assistant_output"] == "Got it — I will check chemnode1."
    assert session.conversation_history.turns[-1].status == "completed"
    assert session.state is not None
    assert session.state.pending_ask_user is None
    assert session.state.pending_messages is None

    continued_messages = json.dumps(
        provider.calls[1]["messages"], sort_keys=True
    )
    assert "which server?" in continued_messages
    assert "chemnode1" in continued_messages

    prompt_context = session.conversation_history.prompt_context(
        current_turn_index=2
    )
    context_blob = json.dumps(prompt_context, sort_keys=True)
    assert "which server?" in context_blob
    assert "chemnode1" in context_blob


def test_identity_prompt_mentions_ask_user_for_truly_ambiguous_cases():
    prompt = build_system_prompt(
        registry=ToolRegistry.default(),
        stage_instructions="Follow the tool loop.",
    )

    assert "ask_user" in prompt
    assert "truly ambiguous" in prompt
