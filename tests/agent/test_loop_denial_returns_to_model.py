from __future__ import annotations

import json

from chemsmart.agent.core import DecisionLog
from chemsmart.agent.handles import HandleStore
from chemsmart.agent.loop import ToolLoop
from chemsmart.agent.permissions import PermissionMode, PermissionPolicy

from ._agent_session_helpers import FakeProvider
from ._loop_helpers import (
    ScriptedRegistry,
    anthropic_final_response,
    anthropic_tool_use,
    anthropic_tool_use_response,
    openai_final_response,
    openai_tool_call_response,
    tool_call,
)


def test_openai_denial_returns_tool_result_message(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call_1", "run_local", {"job": "job_1"})
                )
            },
            {"__raw_response__": openai_final_response("Denied.")},
        ]
    )
    loop = ToolLoop(
        provider=provider,
        registry=ScriptedRegistry({"run_local": {"ok": True}}),
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
        policy=PermissionPolicy(mode=PermissionMode.DRIVING),
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Run."}],
        tool_defs=[],
    )

    tool_message = next(
        message
        for message in provider.calls[1]["messages"]
        if message["role"] == "tool"
    )
    assert tool_message["role"] == "tool"
    assert tool_message["tool_call_id"] == "call_1"
    assert json.loads(tool_message["content"]) == {
        "ok": False,
        "error": {
            "type": "PermissionDenied",
            "message": "missing_yolo",
            "tool": "run_local",
        },
    }
    assert result["tool_outcomes"][0].status == "denied"


def test_anthropic_denial_returns_tool_result_block(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": anthropic_tool_use_response(
                    anthropic_tool_use(
                        "toolu_01A", "run_local", {"job": "job_1"}
                    )
                )
            },
            {"__raw_response__": anthropic_final_response("Denied.")},
        ]
    )
    provider.name = "anthropic"
    loop = ToolLoop(
        provider=provider,
        registry=ScriptedRegistry({"run_local": {"ok": True}}),
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
        policy=PermissionPolicy(mode=PermissionMode.DRIVING),
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Run."}],
        tool_defs=[],
    )

    tool_result_message = next(
        message
        for message in provider.calls[1]["messages"]
        if message["role"] == "user" and isinstance(message["content"], list)
    )
    assert tool_result_message["role"] == "user"
    assert len(tool_result_message["content"]) == 1
    block = tool_result_message["content"][0]
    assert block["type"] == "tool_result"
    assert block["tool_use_id"] == "toolu_01A"
    assert block["is_error"] is True
    assert json.loads(block["content"]) == {
        "ok": False,
        "error": {
            "type": "PermissionDenied",
            "message": "missing_yolo",
            "tool": "run_local",
        },
    }
    assert result["tool_outcomes"][0].status == "denied"
