from __future__ import annotations

from pydantic import create_model

from chemsmart.agent.core import DecisionLog
from chemsmart.agent.handles import HandleStore
from chemsmart.agent.loop import ToolLoop
from chemsmart.agent.permissions import RuntimePermissionMode
from chemsmart.agent.registry import ToolInputModel, ToolRegistry, ToolSpec
from chemsmart.agent.tool_protocol import RuntimeToolMetadata

from ._agent_session_helpers import FakeProvider
from ._loop_helpers import (
    ScriptedRegistry,
    openai_final_response,
    openai_tool_call_response,
    tool_call,
)


def test_tool_loop_runs_single_tool_then_final_assistant(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call_1", "recommend_method", {"task": "opt"}),
                    content="I'll recommend a method first.",
                )
            },
            {"__raw_response__": openai_final_response("Use B3LYP/6-31G*.")},
        ]
    )
    registry = ScriptedRegistry(
        {
            "recommend_method": {
                "method": "B3LYP/6-31G*",
                "why": "Balanced default",
            }
        }
    )
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Recommend a method."}],
        tool_defs=registry.openai_tool_defs(),
    )

    assert result["assistant_text"] == "Use B3LYP/6-31G*."
    assert result["stop_reason"] == "stop"
    assert result["model_steps"] == 2
    assert result["limit_reason"] is None
    assert [request.name for request in result["tool_requests"]] == [
        "recommend_method"
    ]
    assert [outcome.status for outcome in result["tool_outcomes"]] == ["ok"]
    assert registry.calls == [
        ("recommend_method", {"task": "opt"}),
    ]
    assert provider.calls[0]["tools"] == registry.openai_tool_defs()
    assert any(
        message["role"] == "tool" for message in provider.calls[1]["messages"]
    )


def test_tool_loop_mode_none_keeps_static_tool_schema(tmp_path):
    provider = FakeProvider(
        [{"__raw_response__": openai_final_response("No tools needed.")}]
    )
    registry = ToolRegistry(
        [
            _make_tool_spec(
                "read_only_tool",
                metadata=RuntimeToolMetadata(read_only=True),
            ),
            _make_tool_spec(
                "edit_safe_tool",
                metadata=RuntimeToolMetadata(edit_safe=True),
            ),
        ]
    )
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
    )

    tool_defs = registry.openai_tool_defs()
    result = loop.run_turn(
        messages=[{"role": "user", "content": "Just answer."}],
        tool_defs=tool_defs,
        mode=None,
    )

    assert result["assistant_text"] == "No tools needed."
    assert provider.calls[0]["tools"] == tool_defs
    assert provider.calls[0]["tools"] == registry.openai_tool_defs()


def test_tool_loop_read_only_mode_filters_to_read_safe_tools(tmp_path):
    provider = FakeProvider(
        [{"__raw_response__": openai_final_response("Read-only answer.")}]
    )
    registry = ToolRegistry(
        [
            _make_tool_spec(
                "read_only_tool",
                metadata=RuntimeToolMetadata(read_only=True),
            ),
            _make_tool_spec(
                "edit_safe_tool",
                metadata=RuntimeToolMetadata(edit_safe=True),
            ),
        ]
    )
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
    )

    loop.run_turn(
        messages=[{"role": "user", "content": "Stay read only."}],
        mode=RuntimePermissionMode.READ_ONLY,
    )

    tool_names = [
        tool_def["function"]["name"] for tool_def in provider.calls[0]["tools"]
    ]
    assert tool_names == ["read_only_tool", "ask_user"]


def _make_tool_spec(
    name: str,
    *,
    metadata: RuntimeToolMetadata,
) -> ToolSpec:
    def tool() -> dict[str, str]:
        """Test tool."""

        return {"tool": name}

    schema = create_model(
        f"{name.title().replace('_', '')}Input",
        __base__=ToolInputModel,
    )
    return ToolSpec(
        name=name,
        func=tool,
        input_schema=schema,
        metadata=metadata,
    )
