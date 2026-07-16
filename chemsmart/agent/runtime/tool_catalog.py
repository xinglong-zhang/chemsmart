"""Phase-scoped tool exposure independent of model provider syntax."""

from __future__ import annotations

from enum import Enum
from typing import Any, Iterable

from pydantic import BaseModel, ConfigDict

from chemsmart.agent.runtime.contracts import ProviderRole, TaskPhase


class ToolExposure(str, Enum):
    DIRECT = "direct"
    DEFERRED = "deferred"
    HIDDEN = "hidden"


class ToolSelection(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    phase: TaskPhase
    provider_role: ProviderRole
    direct: tuple[str, ...]
    deferred: tuple[str, ...]
    hidden: tuple[str, ...]


_PHASE_TOOLS: dict[TaskPhase, tuple[str, ...]] = {
    TaskPhase.PROJECT: (
        "extract_project_protocol",
        "render_project_yaml",
        "validate_project_yaml",
        "critic_project_yaml",
        "search_basis_sets",
    ),
    TaskPhase.PROJECT_READ: (
        "read_project_yaml",
        "validate_project_yaml",
        "critic_project_yaml",
        "search_basis_sets",
    ),
    TaskPhase.PROJECT_WRITE: (
        "read_project_yaml",
        "validate_project_yaml",
        "critic_project_yaml",
        "write_project_yaml",
        "update_project_yaml",
    ),
    TaskPhase.SYNTHESIS: (
        "read_project_yaml",
        "synthesize_command",
        "repair_command",
    ),
    TaskPhase.VALIDATION: (
        "read_project_yaml",
        "validate_project_yaml",
        "critic_project_yaml",
        "synthesize_command",
        "repair_command",
    ),
    TaskPhase.REPAIR: (
        "repair_command",
        "read_project_yaml",
        "validate_project_yaml",
        "search_basis_sets",
        "synthesize_command",
    ),
    TaskPhase.EXECUTION: (
        "synthesize_command",
        "repair_command",
        "execute_chemsmart_command",
        "read_project_yaml",
    ),
    TaskPhase.DIAGNOSTICS: (
        "inspect_calculation",
        "read_project_yaml",
        "log_tail",
        "scheduler_query",
    ),
}
_SPECIALIST_TOOLS = ("synthesize_command", "repair_command")


class ToolCatalog:
    def __init__(self, registry: Any) -> None:
        self.registry = registry

    def select(
        self,
        *,
        phase: TaskPhase,
        provider_role: ProviderRole,
    ) -> ToolSelection:
        available = tuple(tool.name for tool in self.registry.list_tools())
        requested = (
            _SPECIALIST_TOOLS
            if provider_role is ProviderRole.SYNTHESIS_SPECIALIST
            else _PHASE_TOOLS.get(phase, ())
        )
        direct = tuple(name for name in requested if name in available)
        if len(direct) > 5:
            raise ValueError(
                "runtime phases may expose at most five real tools"
            )
        phase_capabilities = {
            name for names in _PHASE_TOOLS.values() for name in names
        }
        deferred = tuple(
            name
            for name in available
            if name not in direct and name in phase_capabilities
        )
        hidden = tuple(
            name
            for name in available
            if name not in direct and name not in deferred
        )
        return ToolSelection(
            phase=phase,
            provider_role=provider_role,
            direct=direct,
            deferred=deferred,
            hidden=hidden,
        )

    def provider_tool_defs(
        self,
        provider_name: str,
        selection: ToolSelection,
    ) -> list[dict[str, Any]]:
        tools = [
            tool
            for name in selection.direct
            if (tool := self.registry.get_tool(name)) is not None
        ]
        return self.registry.tool_defs_for_provider(provider_name, tools)

    @staticmethod
    def exposure_for(
        name: str,
        selection: ToolSelection,
    ) -> ToolExposure:
        if name in selection.direct:
            return ToolExposure.DIRECT
        if name in selection.deferred:
            return ToolExposure.DEFERRED
        return ToolExposure.HIDDEN


def filter_tool_names(
    names: Iterable[str],
    selection: ToolSelection,
) -> tuple[str, ...]:
    allowed = set(selection.direct)
    return tuple(name for name in names if name in allowed)


__all__ = [
    "ToolCatalog",
    "ToolExposure",
    "ToolSelection",
    "filter_tool_names",
]
