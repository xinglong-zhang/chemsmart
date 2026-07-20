"""Shared wizard probe dispatch and scalar normalization."""

from __future__ import annotations

import os

from chemsmart.agent.wizard.probe import (
    ALL_PROBE_SPECS,
    ProbeSpec,
    run_local_probe,
    run_ssh_probe,
)
from chemsmart.agent.wizard.topology import Topology


def run_probe(runner, topology: Topology, spec: ProbeSpec, **slots: str):
    if topology.mode == "A":
        return run_local_probe(runner, spec, **slots)
    if topology.mode == "B" and topology.host:
        return run_ssh_probe(runner, topology.host, spec, **slots)
    raise ValueError(f"Unsupported topology: {topology}")


def first_nonempty_line(text: str) -> str | None:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return None


def normalize_shell_value(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    if not stripped or stripped.startswith("$"):
        return None
    return stripped


def probe_env_values(
    runner,
    topology: Topology,
    names: list[str],
) -> list[str | None]:
    values: list[str | None] = []
    for name in names:
        result = run_probe(
            runner,
            topology,
            ALL_PROBE_SPECS["common.printenv_var"],
            env_name=name,
        )
        value = normalize_shell_value(first_nonempty_line(result.stdout))
        if value is None and topology.mode == "A":
            value = normalize_shell_value(os.environ.get(name))
        values.append(value)
    return values


__all__ = [
    "first_nonempty_line",
    "normalize_shell_value",
    "probe_env_values",
    "run_probe",
]
