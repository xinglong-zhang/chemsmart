"""Scratch-directory discovery for the wizard."""

from __future__ import annotations

import os
from dataclasses import dataclass

from chemsmart.agent.wizard.probe import (
    ALL_PROBE_SPECS,
    ProbeSpec,
    run_local_probe,
    run_ssh_probe,
)
from chemsmart.agent.wizard.topology import Topology


@dataclass(frozen=True)
class ScratchFinding:
    path: str | None
    source: str
    writable: bool
    candidates: list[tuple[str, str]]


_ENV_SOURCES = [
    ("SCRATCH", "env:SCRATCH"),
    ("WORK", "env:WORK"),
    ("TMPDIR", "env:TMPDIR"),
]
_HOME_SCRATCH = "~/scratch"


def discover_scratch(runner, topology: Topology) -> ScratchFinding:
    """Discover a likely scratch location and whether it is writable."""

    candidates = _discover_env_candidates(runner, topology)
    home_candidate = _discover_home_scratch_candidate(runner, topology)
    if home_candidate is not None:
        candidates.append(home_candidate)

    for source, path in candidates:
        actual_path = _actual_probe_path(runner, topology, source, path)
        if _is_writable(runner, topology, actual_path):
            return ScratchFinding(
                path=path,
                source=source,
                writable=True,
                candidates=candidates,
            )

    if candidates:
        source, path = candidates[0]
        return ScratchFinding(
            path=path,
            source=source,
            writable=False,
            candidates=candidates,
        )

    return ScratchFinding(
        path=None,
        source="none",
        writable=False,
        candidates=[],
    )


def _discover_env_candidates(
    runner,
    topology: Topology,
) -> list[tuple[str, str]]:
    candidates: list[tuple[str, str]] = []
    for name, source in _ENV_SOURCES:
        result = _run_probe(
            runner,
            topology,
            ALL_PROBE_SPECS["common.printenv_var"],
            env_name=name,
        )
        value = _normalize_shell_value(_first_nonempty_line(result.stdout))
        if value is None and topology.mode == "A":
            value = _normalize_shell_value(os.environ.get(name))
        if value:
            candidates.append((source, value))
    return candidates


def _discover_home_scratch_candidate(
    runner,
    topology: Topology,
) -> tuple[str, str] | None:
    result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["common.printenv_var"],
        env_name="HOME",
    )
    home = _normalize_shell_value(_first_nonempty_line(result.stdout))
    if home is None and topology.mode == "A":
        home = _normalize_shell_value(os.environ.get("HOME"))
    if not home:
        return None

    home_scratch = f"{home.rstrip('/')}/scratch"
    home_result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["scratch.test_dir_writable"],
        path=home_scratch,
    )
    if home_result.returncode == 0:
        return ("home:~/scratch", _HOME_SCRATCH)
    return None


def _actual_probe_path(
    runner, topology: Topology, source: str, path: str
) -> str:
    if source != "home:~/scratch":
        return path
    result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["common.printenv_var"],
        env_name="HOME",
    )
    home = _normalize_shell_value(_first_nonempty_line(result.stdout))
    if home is None and topology.mode == "A":
        home = _normalize_shell_value(os.environ.get("HOME"))
    if not home:
        return path
    return f"{home.rstrip('/')}/scratch"


def _is_writable(runner, topology: Topology, path: str) -> bool:
    result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["scratch.test_writable"],
        path=path,
    )
    return result.returncode == 0


def _run_probe(runner, topology: Topology, spec: ProbeSpec, **slots: str):
    if topology.mode == "A":
        return run_local_probe(runner, spec, **slots)
    if topology.mode == "B" and topology.host:
        return run_ssh_probe(runner, topology.host, spec, **slots)
    raise ValueError(f"Unsupported topology: {topology}")


def _first_nonempty_line(text: str) -> str | None:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return None


def _normalize_shell_value(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    if not stripped or stripped.startswith("$"):
        return None
    return stripped
