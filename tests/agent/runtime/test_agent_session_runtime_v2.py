from __future__ import annotations

import json

from pydantic import Field, create_model

from chemsmart.agent.core import AgentSession
from chemsmart.agent.registry import ToolInputModel, ToolRegistry, ToolSpec
from chemsmart.agent.tool_protocol import RuntimeToolMetadata

from .._agent_session_helpers import FakeProvider
from .._loop_helpers import (
    anthropic_final_response,
    anthropic_tool_use,
    anthropic_tool_use_response,
    openai_final_response,
    openai_tool_call_response,
    tool_call,
)


def test_active_runtime_limits_tools_and_writes_replayable_events(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call(
                        "call-1",
                        "synthesize_command",
                        {"request": "optimize water"},
                    )
                )
            },
            {"__raw_response__": openai_final_response("Command prepared.")},
        ]
    )
    command = "chemsmart run gaussian -p demo -f water.xyz -c 0 -m 1 opt"
    registry = _registry(
        _tool(
            "synthesize_command",
            lambda request: {
                "ok": True,
                "status": "ready",
                "command": command,
                "semantic": {"verdict": "ok"},
            },
            fields={"request": (str, Field(...))},
            read_only=True,
        ),
        _tool("build_job", lambda: {"ok": True}, read_only=True),
    )
    session = AgentSession(
        provider=provider,
        registry=registry,
        session_root=tmp_path / "sessions",
        runtime_v2="active",
    )

    result = session.run_loop("Prepare a Gaussian optimization for water.")

    exposed = [item["function"]["name"] for item in provider.calls[0]["tools"]]
    assert exposed == ["synthesize_command", "ask_user"]
    assert result["runtime_v2"]["mode"] == "active"
    assert result["runtime_v2"]["phase"] == "complete"
    assert session.session_dir is not None
    assert (session.session_dir / "runtime_events.jsonl").is_file()
    assert (session.session_dir / "runtime_state.json").is_file()
    assert session._runtime_controller is not None
    assert session._runtime_controller.state.previous_command == command
    metadata = json.loads(
        (session.session_dir / "session_metadata.json").read_text()
    )
    assert metadata["runtime_v2_mode"] == "active"
    assert metadata["runtime_v2_phase"] == "complete"
    assert metadata["runtime_v2_event_count"] > 0


def test_active_runtime_uses_same_contract_with_anthropic_tools(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": anthropic_tool_use_response(
                    anthropic_tool_use(
                        "call-1",
                        "synthesize_command",
                        {"request": "optimize water"},
                    )
                )
            },
            {"__raw_response__": anthropic_final_response("Prepared.")},
        ]
    )
    provider.name = "anthropic"
    registry = _registry(
        _tool(
            "synthesize_command",
            lambda request: {
                "ok": True,
                "status": "ready",
                "command": "chemsmart run gaussian -p demo -f water.xyz -c 0 -m 1 opt",
            },
            fields={"request": (str, Field(...))},
            read_only=True,
        )
    )
    session = AgentSession(
        provider=provider,
        registry=registry,
        session_root=tmp_path / "sessions",
        runtime_v2="active",
    )

    result = session.run_loop("Prepare a Gaussian optimization for water.")

    assert provider.calls[0]["tools"][0]["name"] == "synthesize_command"
    assert result["assistant_output"] == "Prepared."
    assert result["runtime_v2"]["phase"] == "complete"


def test_shadow_runtime_observes_legacy_tool_without_changing_behavior(
    tmp_path,
):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call-1", "build_job", {})
                )
            },
            {"__raw_response__": openai_final_response("Legacy path used.")},
        ]
    )
    registry = _registry(
        _tool(
            "synthesize_command",
            lambda request: {},
            fields={"request": (str, Field(...))},
            read_only=True,
        ),
        _tool("build_job", lambda: {"ok": True}, read_only=True),
    )
    session = AgentSession(
        provider=provider,
        registry=registry,
        session_root=tmp_path / "sessions",
        runtime_v2="shadow",
    )

    result = session.run_loop("Prepare a Gaussian optimization for water.")

    exposed = {item["function"]["name"] for item in provider.calls[0]["tools"]}
    assert "build_job" in exposed
    assert result["tool_outcomes"][0].status == "ok"
    assert result["runtime_v2"]["shadow_violations"] == [
        "runtime.tool.not_exposed"
    ]


def test_active_runtime_blocks_provider_call_to_hidden_tool(tmp_path):
    called = []
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call-1", "build_job", {})
                )
            },
            {"__raw_response__": openai_final_response("Call rejected.")},
        ]
    )
    registry = _registry(
        _tool(
            "synthesize_command",
            lambda request: {},
            fields={"request": (str, Field(...))},
            read_only=True,
        ),
        _tool(
            "build_job",
            lambda: called.append(True) or {"ok": True},
            read_only=True,
        ),
    )
    session = AgentSession(
        provider=provider,
        registry=registry,
        session_root=tmp_path / "sessions",
        runtime_v2="active",
    )

    result = session.run_loop("Prepare a Gaussian optimization for water.")

    assert called == []
    assert result["tool_outcomes"][0].status == "error"
    assert result["tool_outcomes"][0].error_type == "ToolExposureViolation"


def test_project_write_remains_approval_gated_in_active_runtime(tmp_path):
    called = []
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call(
                        "call-1",
                        "write_project_yaml",
                        {"project_name": "co2", "yaml_text": "gas: {}"},
                    )
                )
            },
            {"__raw_response__": openai_final_response("Write not approved.")},
        ]
    )

    def write_project_yaml(project_name, yaml_text):
        called.append((project_name, yaml_text))
        return {"ok": True}

    registry = _registry(
        _tool(
            "write_project_yaml",
            write_project_yaml,
            fields={
                "project_name": (str, Field(...)),
                "yaml_text": (str, Field(...)),
            },
            read_only=False,
        )
    )
    session = AgentSession(
        provider=provider,
        registry=registry,
        session_root=tmp_path / "sessions",
        runtime_v2="active",
    )

    result = session.run_loop("Write the validated project YAML now.")

    assert called == []
    assert result["tool_outcomes"][0].status == "denied"
    assert result["denials_count"] == 1


def test_runtime_off_keeps_legacy_session_artifacts(tmp_path):
    provider = FakeProvider(
        [{"__raw_response__": openai_final_response("No tools needed.")}]
    )
    session = AgentSession(
        provider=provider,
        registry=_registry(
            _tool(
                "synthesize_command",
                lambda request: {},
                fields={"request": (str, Field(...))},
                read_only=True,
            )
        ),
        session_root=tmp_path / "sessions",
        runtime_v2="off",
    )

    result = session.run_loop("Explain the current setup.")

    assert result["runtime_v2"] == {"mode": "off"}
    assert session.session_dir is not None
    assert not (session.session_dir / "runtime_events.jsonl").exists()


def _registry(*tools):
    return ToolRegistry(list(tools))


def _tool(name, func, *, fields=None, read_only):
    model = create_model(
        f"{name.title().replace('_', '')}Input",
        __base__=ToolInputModel,
        **(fields or {}),
    )
    return ToolSpec(
        name=name,
        func=func,
        input_schema=model,
        metadata=RuntimeToolMetadata(read_only=read_only),
    )
