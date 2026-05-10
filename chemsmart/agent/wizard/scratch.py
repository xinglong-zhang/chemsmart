"""Scratch-directory discovery for the wizard."""

from __future__ import annotations

import os
import shlex
from dataclasses import dataclass

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
    home_result = _run_probe(
        runner,
        topology,
        ["test", "-d", _HOME_SCRATCH, "-a", "-w", _HOME_SCRATCH],
        remote_command="test -d ~/scratch -a -w ~/scratch",
    )
    if home_result.returncode == 0:
        candidates.append(("home:~/scratch", _HOME_SCRATCH))

    for source, path in candidates:
        if _is_writable(runner, topology, path):
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
    names = [name for name, _ in _ENV_SOURCES]
    remote_parts = [f'"${name}"' for name in names]
    result = _run_probe(
        runner,
        topology,
        ["printf", "%s\\n", *[f"${name}" for name in names]],
        remote_command=f"printf '%s\\n' {' '.join(remote_parts)}",
    )
    values = [
        _normalize_shell_value(line) for line in result.stdout.splitlines()
    ]
    if topology.mode == "A" and result.stdout.strip() and not any(values):
        values = [
            _normalize_shell_value(os.environ.get(name)) for name in names
        ]

    candidates: list[tuple[str, str]] = []
    for (_, source), value in zip(_ENV_SOURCES, values):
        if value:
            candidates.append((source, value))
    return candidates


def _is_writable(runner, topology: Topology, path: str) -> bool:
    result = _run_probe(
        runner,
        topology,
        ["test", "-w", path],
        remote_command=f"test -w {shlex.quote(path)}",
    )
    return result.returncode == 0


def _run_probe(
    runner,
    topology: Topology,
    command: list[str],
    remote_command: str | None = None,
):
    if topology.mode == "A":
        return runner.run_local(command)
    if topology.mode == "B" and topology.host:
        return runner.run_ssh(
            topology.host, remote_command or shlex.join(command)
        )
    raise ValueError(f"Unsupported topology: {topology}")


def _normalize_shell_value(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    if not stripped or stripped.startswith("$"):
        return None
    return stripped
