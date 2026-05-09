"""Public agent exports for chemsmart."""

from __future__ import annotations

from importlib import import_module

__all__ = [
    "AgentSession",
    "run_agent",
    "build_molecule",
    "build_gaussian_settings",
    "build_orca_settings",
    "build_job",
    "dry_run_input",
    "recommend_method",
    "validate_runtime",
    "run_local",
    "submit_hpc",
    "SubmitTransport",
    "LocalDryRunTransport",
    "SshQsubTransport",
    "MockTransport",
]

_EXPORTS = {
    "AgentSession": ("chemsmart.agent.core", "AgentSession"),
    "run_agent": ("chemsmart.agent.core", "run_agent"),
    "build_molecule": ("chemsmart.agent.tools", "build_molecule"),
    "build_gaussian_settings": (
        "chemsmart.agent.tools",
        "build_gaussian_settings",
    ),
    "build_orca_settings": (
        "chemsmart.agent.tools",
        "build_orca_settings",
    ),
    "build_job": ("chemsmart.agent.tools", "build_job"),
    "dry_run_input": ("chemsmart.agent.tools", "dry_run_input"),
    "recommend_method": ("chemsmart.agent.tools", "recommend_method"),
    "validate_runtime": ("chemsmart.agent.tools", "validate_runtime"),
    "run_local": ("chemsmart.agent.tools", "run_local"),
    "submit_hpc": ("chemsmart.agent.tools", "submit_hpc"),
    "SubmitTransport": ("chemsmart.agent.transport", "SubmitTransport"),
    "LocalDryRunTransport": (
        "chemsmart.agent.transport",
        "LocalDryRunTransport",
    ),
    "SshQsubTransport": ("chemsmart.agent.transport", "SshQsubTransport"),
    "MockTransport": ("chemsmart.agent.transport", "MockTransport"),
}


def __getattr__(name: str):
    try:
        module_name, attr_name = _EXPORTS[name]
    except KeyError as exc:  # pragma: no cover - import protocol guard
        raise AttributeError(name) from exc
    value = getattr(import_module(module_name), attr_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
