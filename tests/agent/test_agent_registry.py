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


def test_registry_unknown_tool_name_raises_clearly():
    registry = ToolRegistry.default()

    with pytest.raises(ValueError) as excinfo:
        registry.call("does_not_exist", {})

    message = str(excinfo.value)
    assert "does_not_exist" in message
    assert "build_molecule" in message
