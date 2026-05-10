"""Software-environment discovery for the wizard."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

from chemsmart.agent.wizard.probe import (
    ALL_PROBE_SPECS,
    ProbeSpec,
    run_local_probe,
    run_ssh_probe,
)
from chemsmart.agent.wizard.topology import Topology


@dataclass(frozen=True)
class ModuleSystem:
    kind: str
    version: str | None


@dataclass(frozen=True)
class ProgramFinding:
    program: str
    exefolder: str | None
    source: str
    module_candidates: list[str]
    on_path: bool


@dataclass(frozen=True)
class SoftwareSurvey:
    module_system: ModuleSystem
    programs: dict[str, ProgramFinding]
    conda_base: str | None
    conda_env: str | None


def detect_module_system(runner, topology: Topology) -> ModuleSystem:
    """Detect the site module implementation, if any."""

    module_type = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["software.type_module"],
    )
    if module_type.returncode != 0:
        which_module = _run_probe(
            runner,
            topology,
            ALL_PROBE_SPECS["software.which_module"],
        )
        if which_module.returncode != 0:
            return ModuleSystem(kind="none", version=None)

    version_result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["software.module_version"],
    )
    merged_output = _merge_output(version_result)
    version_text = _first_nonempty_line(merged_output)
    kind = "lmod" if "Lmod" in merged_output else "envmodules"
    return ModuleSystem(kind=kind, version=version_text)


def find_program(
    runner,
    topology: Topology,
    key: str,
    exe_names: list[str],
    module_patterns: list[str],
) -> ProgramFinding:
    """Find a program on PATH first, then via module availability."""

    for exe_name in exe_names:
        resolved = _resolve_executable(runner, topology, exe_name)
        if resolved is not None:
            return ProgramFinding(
                program=key,
                exefolder=str(Path(resolved).parent),
                source="path",
                module_candidates=[],
                on_path=True,
            )

    module_candidates = _find_module_candidates(
        runner,
        topology,
        module_patterns,
    )
    if module_candidates:
        return ProgramFinding(
            program=key,
            exefolder=None,
            source="module",
            module_candidates=module_candidates,
            on_path=False,
        )

    return ProgramFinding(
        program=key,
        exefolder=None,
        source="none",
        module_candidates=[],
        on_path=False,
    )


def discover_conda(
    runner, topology: Topology
) -> tuple[str | None, str | None]:
    """Discover the active conda env and base path, if available."""

    env_value = _probe_env_values(runner, topology, ["CONDA_PREFIX"])[0]
    base_result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["software.conda_base"],
    )
    base_value = _first_nonempty_line(base_result.stdout)
    return _normalize_shell_value(base_value), env_value


def run_software_survey(runner, topology: Topology) -> SoftwareSurvey:
    """Collect module, program, and conda findings."""

    module_system = detect_module_system(runner, topology)
    conda_base, conda_env = discover_conda(runner, topology)
    programs = {
        "gaussian": find_program(
            runner,
            topology,
            "gaussian",
            ["g16", "g09"],
            ["gaussian", "g16", "g09"],
        ),
        "orca": find_program(
            runner,
            topology,
            "orca",
            ["orca"],
            ["orca"],
        ),
        "nciplot": find_program(
            runner,
            topology,
            "nciplot",
            ["nciplot"],
            ["nciplot", "nci"],
        ),
    }
    return SoftwareSurvey(
        module_system=module_system,
        programs=programs,
        conda_base=conda_base,
        conda_env=conda_env,
    )


def _resolve_executable(
    runner,
    topology: Topology,
    exe_name: str,
) -> str | None:
    command_result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["software.command_v"],
        exe_name=exe_name,
    )
    path_value = _first_nonempty_line(command_result.stdout)
    if not path_value and topology.mode == "A":
        which_result = _run_probe(
            runner,
            topology,
            ALL_PROBE_SPECS["software.which_exe"],
            exe_name=exe_name,
        )
        path_value = _first_nonempty_line(which_result.stdout)
    path_value = _normalize_shell_value(path_value)
    if not path_value:
        return None
    if not path_value.startswith(("/", "~")):
        return path_value

    readlink_result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["software.readlink"],
        path=path_value,
    )
    resolved = _normalize_shell_value(
        _first_nonempty_line(readlink_result.stdout)
    )
    return resolved or path_value


def _find_module_candidates(
    runner,
    topology: Topology,
    module_patterns: list[str],
) -> list[str]:
    result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["software.module_avail"],
    )
    if result.returncode != 0:
        return []

    pattern = re.compile("|".join(module_patterns), re.IGNORECASE)
    candidates = {
        _normalize_module_line(line)
        for line in _merge_output(result).splitlines()
        if pattern.search(line)
    }
    normalized = [candidate for candidate in candidates if candidate]
    return sorted(normalized, key=lambda name: (len(name), name))


def _normalize_module_line(line: str) -> str | None:
    candidate = line.strip()
    if not candidate:
        return None
    if candidate.endswith(":") or candidate.startswith("--"):
        return None
    if candidate.lower().startswith(("where:", "modulepath", 'use "module')):
        return None
    return candidate.split()[0]


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


def _merge_output(result) -> str:
    stdout = result.stdout.strip()
    stderr = result.stderr.strip()
    if stdout and stderr:
        return f"{stdout}\n{stderr}"
    return stdout or stderr


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
