"""Project/account discovery for the wizard."""

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
class ProjectFinding:
    project: str | None
    source: str
    candidates: list[str]


def discover_project(
    runner,
    topology: Topology,
    scheduler: str | None,
) -> ProjectFinding:
    """Discover a likely scheduler account / project string."""

    scheduler_name = (scheduler or "").upper()
    if scheduler_name == "SLURM":
        env_finding = _discover_env_project(
            runner,
            topology,
            [
                ("SBATCH_ACCOUNT", "env:SBATCH_ACCOUNT"),
                ("SLURM_ACCOUNT", "env:SLURM_ACCOUNT"),
            ],
        )
        if env_finding.project is not None:
            return env_finding

        sacct_finding = _discover_sacctmgr_project(runner, topology)
        if sacct_finding.project is not None:
            return sacct_finding

        return _discover_groups_project(runner, topology)

    if scheduler_name == "PBS":
        env_finding = _discover_env_project(
            runner,
            topology,
            [("PBS_ACCOUNT", "env:PBS_ACCOUNT")],
        )
        if env_finding.project is not None:
            return env_finding
        return _discover_groups_project(runner, topology)

    if scheduler_name in {"LSF", "SGE"}:
        return _discover_groups_project(runner, topology)

    return _discover_groups_project(runner, topology)


def _discover_env_project(
    runner,
    topology: Topology,
    env_sources: list[tuple[str, str]],
) -> ProjectFinding:
    names = [name for name, _ in env_sources]
    values = _probe_env_values(runner, topology, names)
    candidates = [value for value in values if value]
    for value, (_, source) in zip(values, env_sources):
        if value:
            return ProjectFinding(
                project=value,
                source=source,
                candidates=candidates,
            )
    return ProjectFinding(project=None, source="none", candidates=candidates)


def _discover_sacctmgr_project(runner, topology: Topology) -> ProjectFinding:
    user_values = _probe_env_values(runner, topology, ["USER"])
    user = user_values[0] or _normalize_shell_value(os.environ.get("USER"))
    if not user:
        return ProjectFinding(project=None, source="none", candidates=[])

    seen_projects: list[str] = []
    result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["project.sacctmgr_show_user"],
        user=user,
    )
    if result.returncode != 0:
        return ProjectFinding(project=None, source="none", candidates=[])

    for line in result.stdout.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        default_account = _normalize_shell_value(stripped.split("|", 1)[0])
        if default_account:
            seen_projects.append(default_account)
    if seen_projects:
        return ProjectFinding(
            project=seen_projects[0],
            source="sacctmgr",
            candidates=seen_projects,
        )
    return ProjectFinding(project=None, source="none", candidates=[])


def _discover_groups_project(runner, topology: Topology) -> ProjectFinding:
    result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["project.groups"],
    )
    candidates = [
        _normalize_shell_value(token) for token in result.stdout.split()
    ]
    normalized = [token for token in candidates if token]
    if normalized:
        return ProjectFinding(
            project=normalized[0],
            source="groups",
            candidates=normalized,
        )
    return ProjectFinding(project=None, source="none", candidates=[])


def _probe_env_values(
    runner,
    topology: Topology,
    names: list[str],
) -> list[str | None]:
    values: list[str | None] = []
    for name in names:
        result = _run_probe(
            runner,
            topology,
            ALL_PROBE_SPECS["common.printenv_var"],
            env_name=name,
        )
        value = _normalize_shell_value(_first_nonempty_line(result.stdout))
        if value is None and topology.mode == "A":
            value = _normalize_shell_value(os.environ.get(name))
        values.append(value)
    return values


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
