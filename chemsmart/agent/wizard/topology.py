"""Topology detection for scheduler probing."""

from __future__ import annotations

import re
from dataclasses import dataclass

from chemsmart.agent.wizard.probe import ALL_PROBE_SPECS, run_local_probe

_ENV_PATTERN = re.compile(r"^(SLURM_|PBS_|LSB_|SGE_)")


@dataclass(frozen=True)
class Topology:
    mode: str
    host: str | None
    evidence: list[str]


class NoTargetError(Exception):
    """Raised when no local or remote scheduler target can be inferred."""


_LOCAL_SCHEDULER_PROBES = [
    ALL_PROBE_SPECS["topology.sinfo"],
    ALL_PROBE_SPECS["topology.qstat"],
    ALL_PROBE_SPECS["topology.bqueues"],
    ALL_PROBE_SPECS["topology.qconf"],
]


def detect_topology(runner, ssh_host_hint: str | None = None) -> Topology:
    """Detect whether scheduler probes should run locally or over SSH."""

    evidence: list[str] = []
    env_result = run_local_probe(
        runner, ALL_PROBE_SPECS["common.printenv_all"]
    )
    if env_result.returncode == 0:
        matches = [
            line.split("=", 1)[0]
            for line in env_result.stdout.splitlines()
            if _ENV_PATTERN.match(line)
        ]
        if matches:
            evidence.extend(f"env:{name}" for name in matches)
            return Topology(mode="A", host="localhost", evidence=evidence)

    for spec in _LOCAL_SCHEDULER_PROBES:
        result = run_local_probe(runner, spec)
        if result.returncode == 0:
            evidence.append(f"local:{' '.join(spec.argv_template)}")
            return Topology(mode="A", host="localhost", evidence=evidence)

    if ssh_host_hint:
        evidence.append(f"ssh_host_hint:{ssh_host_hint}")
        return Topology(mode="B", host=ssh_host_hint, evidence=evidence)

    raise NoTargetError("No local scheduler evidence or SSH host hint found.")
