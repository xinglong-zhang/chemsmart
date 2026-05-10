"""Project/account discovery for the wizard."""

from __future__ import annotations

import os
import shlex
from dataclasses import dataclass

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
    user = os.environ.get("USER", "$USER")
    commands = [
        (
            [
                "sacctmgr",
                "-n",
                "-p",
                "show",
                "user",
                "$USER",
                "format=DefaultAccount,Account",
            ],
            "sacctmgr -n -p show user $USER format=DefaultAccount,Account",
        ),
        (
            [
                "sacctmgr",
                "-n",
                "-p",
                "show",
                "user",
                user,
                "format=DefaultAccount,Account",
            ],
            (
                "sacctmgr -n -p show user "
                f"{shlex.quote(user)} format=DefaultAccount,Account"
            ),
        ),
    ]
    seen_projects: list[str] = []
    for command, remote_command in commands:
        result = _run_probe(
            runner,
            topology,
            command,
            remote_command=remote_command,
        )
        if result.returncode != 0:
            continue
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
        ["groups"],
        remote_command="groups",
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
    remote_parts = [f'"${name}"' for name in names]
    result = _run_probe(
        runner,
        topology,
        ["printf", "%s\\n", *[f"${name}" for name in names]],
        remote_command=f"printf '%s\\n' {' '.join(remote_parts)}",
    )
    values = [
        _normalize_shell_value(line)
        for line in result.stdout.splitlines()[: len(names)]
    ]
    while len(values) < len(names):
        values.append(None)

    if (
        topology.mode == "A"
        and result.stdout.strip()
        and all(value is None for value in values)
    ):
        return [_normalize_shell_value(os.environ.get(name)) for name in names]
    return values


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
