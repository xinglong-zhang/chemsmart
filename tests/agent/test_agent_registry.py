from __future__ import annotations

import pytest

from chemsmart.agent.registry import ToolRegistry
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
