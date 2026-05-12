from __future__ import annotations

import pytest
from pydantic import create_model

from chemsmart.agent.permissions import RuntimePermissionMode
from chemsmart.agent.registry import ToolInputModel, ToolRegistry, ToolSpec
from chemsmart.agent.tool_protocol import RuntimeToolMetadata
from chemsmart.jobs.gaussian.settings import GaussianJobSettings


def test_registry_round_trips_build_gaussian_settings_call():
    registry = ToolRegistry.default()

    result = registry.call(
        "build_gaussian_settings",
        {"functional": "B3LYP", "basis": "6-31G*", "charge": 0},
    )

    assert isinstance(result, GaussianJobSettings)
    assert result.functional == "B3LYP"
    assert result.basis == "6-31G*"

    tool_defs = registry.openai_tool_defs()
    names = [tool_def["function"]["name"] for tool_def in tool_defs]
    assert "build_gaussian_settings" in names
    assert "wizard_probe" in names
    assert "wizard_refresh" in names
    assert "wizard_verify" in names
    assert "wizard_write" in names
    build_job_def = next(
        tool_def
        for tool_def in tool_defs
        if tool_def["function"]["name"] == "build_job"
    )
    gaussian_settings_def = next(
        tool_def
        for tool_def in tool_defs
        if tool_def["function"]["name"] == "build_gaussian_settings"
    )
    assert build_job_def["function"]["parameters"]["properties"]["kind"] == {
        "type": "string",
        "enum": [
            "gaussian.opt",
            "gaussian.ts",
            "gaussian.freq",
            "gaussian.sp",
            "gaussian.irc",
            "gaussian.scan",
            "orca.opt",
            "orca.ts",
            "orca.freq",
            "orca.sp",
            "orca.irc",
            "orca.scan",
        ],
    }
    gaussian_settings_props = gaussian_settings_def["function"]["parameters"][
        "properties"
    ]
    assert "title" in gaussian_settings_props
    assert "freq" in gaussian_settings_props
    assert "numfreq" in gaussian_settings_props
    assert "additional_opt_options_in_route" in gaussian_settings_props
    assert "additional_route_parameters" in gaussian_settings_props


def test_registry_unknown_tool_name_raises_clearly():
    registry = ToolRegistry.default()

    with pytest.raises(ValueError) as excinfo:
        registry.call("does_not_exist", {})

    message = str(excinfo.value)
    assert "does_not_exist" in message
    assert "build_molecule" in message


def test_registry_forbids_unexpected_tool_arguments():
    registry = ToolRegistry.default()

    result = registry.call(
        "build_gaussian_settings",
        {
            "functional": "B3LYP",
            "basis": "6-31G*",
            "unexpected": "boom",
        },
    )

    assert result["ok"] is False
    assert result["error"]["type"] == "ValidationError"
    assert "unexpected" in result["error"]["message"]


def test_registry_exposes_machine_readable_tool_descriptions():
    registry = ToolRegistry.default()

    build_molecule_def = next(
        tool_def
        for tool_def in registry.openai_tool_defs()
        if tool_def["function"]["name"] == "build_molecule"
    )

    assert "Load one molecule" in build_molecule_def["function"]["description"]


def test_registry_default_registration_uses_empty_runtime_metadata():
    registry = ToolRegistry.default()

    tools = registry.list_tools()

    assert tools
    assert all(tool.metadata == RuntimeToolMetadata() for tool in tools)


def test_assemble_tool_pool_plan_mode_returns_empty_list():
    registry = ToolRegistry([_make_tool_spec("legacy_tool")])

    pool = registry.assemble_tool_pool(RuntimePermissionMode.PLAN)

    assert pool == []


@pytest.mark.parametrize(
    ("mode", "tool_name", "expected"),
    [
        (RuntimePermissionMode.READ_ONLY, "legacy_tool", True),
        (RuntimePermissionMode.READ_ONLY, "read_only_tool", True),
        (RuntimePermissionMode.READ_ONLY, "edit_safe_tool", False),
        (RuntimePermissionMode.READ_ONLY, "bypass_only_tool", False),
        (RuntimePermissionMode.ACCEPT_EDITS, "legacy_tool", True),
        (RuntimePermissionMode.ACCEPT_EDITS, "read_only_tool", True),
        (RuntimePermissionMode.ACCEPT_EDITS, "edit_safe_tool", True),
        (RuntimePermissionMode.ACCEPT_EDITS, "bypass_only_tool", False),
        (RuntimePermissionMode.BYPASS, "legacy_tool", True),
        (RuntimePermissionMode.BYPASS, "read_only_tool", True),
        (RuntimePermissionMode.BYPASS, "edit_safe_tool", True),
        (RuntimePermissionMode.BYPASS, "bypass_only_tool", True),
        (RuntimePermissionMode.PLAN, "legacy_tool", False),
        (RuntimePermissionMode.PLAN, "read_only_tool", False),
        (RuntimePermissionMode.PLAN, "edit_safe_tool", False),
        (RuntimePermissionMode.PLAN, "bypass_only_tool", False),
    ],
)
def test_assemble_tool_pool_runtime_mode_matrix(mode, tool_name, expected):
    registry = ToolRegistry(
        [
            _make_tool_spec("legacy_tool"),
            _make_tool_spec(
                "read_only_tool",
                metadata=RuntimeToolMetadata(read_only=True),
            ),
            _make_tool_spec(
                "edit_safe_tool",
                metadata=RuntimeToolMetadata(edit_safe=True),
            ),
            _make_tool_spec(
                "bypass_only_tool",
                metadata=RuntimeToolMetadata(
                    requires_mode=frozenset({RuntimePermissionMode.BYPASS})
                ),
            ),
        ]
    )

    names = {tool.name for tool in registry.assemble_tool_pool(mode)}

    assert (tool_name in names) is expected


def _make_tool_spec(
    name: str,
    *,
    metadata: RuntimeToolMetadata | None = None,
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
        metadata=metadata or RuntimeToolMetadata(),
    )
