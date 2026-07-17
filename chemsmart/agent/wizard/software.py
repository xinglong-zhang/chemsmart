"""Software-environment discovery for the wizard."""

from __future__ import annotations

import re
import shlex
from dataclasses import dataclass
from pathlib import Path

from chemsmart.agent.wizard.probe import ALL_PROBE_SPECS
from chemsmart.agent.wizard.probe_values import (
    first_nonempty_line as _first_nonempty_line,
)
from chemsmart.agent.wizard.probe_values import (
    normalize_shell_value as _normalize_shell_value,
)
from chemsmart.agent.wizard.probe_values import (
    probe_env_values as _probe_env_values,
)
from chemsmart.agent.wizard.probe_values import run_probe as _run_probe
from chemsmart.agent.wizard.topology import Topology

_WELL_KNOWN_CONDA_PATHS = (
    "~/miniforge3/bin/conda",
    "/opt/miniforge3/bin/conda",
    "~/anaconda3/bin/conda",
    "/opt/anaconda3/bin/conda",
    "~/miniconda3/bin/conda",
    "/opt/conda/bin/conda",
)
_MODULE_ENVVAR_PREFERENCES = {
    "gaussian": ("g16root", "GAUSS_EXEDIR"),
    "orca": ("ORCA_ROOT",),
    "nciplot": ("NCIPLOT_HOME",),
}
_LMOD_SETENV_PATTERN = re.compile(
    r'setenv\{\s*"(?P<name>[^"]+)"\s*,\s*"(?P<value>[^"]*)"\s*\}'
)
_LMOD_PREPEND_PATH_PATTERN = re.compile(
    r'prepend_path\{\s*"PATH"\s*,\s*"(?P<value>[^"]+)"'
)


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
class CondaEnvSurvey:
    base: str | None
    env_path: str | None
    env_name: str | None


@dataclass(frozen=True)
class SoftwareSurvey:
    module_system: ModuleSystem
    programs: dict[str, ProgramFinding]
    conda: CondaEnvSurvey


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
    verify_qchem: bool = False,
) -> ProgramFinding:
    """Find a program on PATH first, then via module availability.

    verify_qchem=True runs ``<exe> --version`` and requires "orca" +
    "version"/"release" in the output — this filters out the Linux
    accessibility tool ``/usr/bin/orca`` which shares the same name.
    """

    for exe_name in exe_names:
        resolved = _resolve_executable(runner, topology, exe_name)
        if resolved is not None:
            if verify_qchem and not _verify_orca_qchem(
                runner, topology, resolved
            ):
                continue
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
        exefolder = _probe_module_exefolder(
            runner,
            topology,
            key,
            module_candidates,
        )
        return ProgramFinding(
            program=key,
            exefolder=exefolder,
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
    runner,
    topology: Topology,
) -> CondaEnvSurvey:
    """Discover the active conda env and base path, if available."""

    env_value = _probe_env_values(runner, topology, ["CONDA_PREFIX"])[0]
    conda_path = _resolve_conda_executable(runner, topology)
    base = None
    env_path = env_value
    if conda_path is not None:
        base = _probe_conda_base(runner, topology, conda_path)
        if env_path is None:
            env_path = _select_conda_env_path(
                base=base,
                envs=_probe_conda_env_list(runner, topology, conda_path),
            )
    return CondaEnvSurvey(
        base=base,
        env_path=env_path,
        env_name=_derive_conda_env_name(base=base, env_path=env_path),
    )


def run_software_survey(runner, topology: Topology) -> SoftwareSurvey:
    """Collect module, program, and conda findings."""

    module_system = detect_module_system(runner, topology)
    conda = discover_conda(runner, topology)
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
            verify_qchem=True,
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
        conda=conda,
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


def _resolve_conda_executable(runner, topology: Topology) -> str | None:
    resolved = _resolve_executable(runner, topology, "conda")
    if resolved is not None:
        return resolved

    home = _probe_env_values(runner, topology, ["HOME"])[0]
    for candidate in _WELL_KNOWN_CONDA_PATHS:
        expanded = _expand_home_path(candidate, home=home, topology=topology)
        result = _run_probe(
            runner,
            topology,
            ALL_PROBE_SPECS["software.test_executable"],
            path=expanded,
        )
        if result.returncode == 0:
            return expanded
    return None


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
    candidate = candidate.split()[0]
    if candidate.endswith("/"):
        return None
    return candidate


def _probe_module_exefolder(
    runner,
    topology: Topology,
    program: str,
    module_candidates: list[str],
) -> str | None:
    for module_name in module_candidates:
        result = _run_probe(
            runner,
            topology,
            ALL_PROBE_SPECS["software.module_show"],
            module_name=module_name,
        )
        if result.returncode != 0:
            continue
        exefolder = _extract_exefolder_from_module_show(
            program=program,
            payload=_merge_output(result),
        )
        if exefolder is not None:
            return exefolder
    return None


def _extract_exefolder_from_module_show(
    program: str,
    payload: str,
) -> str | None:
    env_values, path_entries = _parse_module_show(payload)
    for env_name in _MODULE_ENVVAR_PREFERENCES.get(program, ()):
        value = _normalize_shell_value(env_values.get(env_name))
        if value is not None:
            return value
    for entry in path_entries:
        value = _normalize_shell_value(entry)
        if value is not None:
            return value
    return None


def _parse_module_show(payload: str) -> tuple[dict[str, str], list[str]]:
    env_values: dict[str, str] = {}
    path_entries: list[str] = []
    for raw_line in payload.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        lmod_env_match = _LMOD_SETENV_PATTERN.search(line)
        if lmod_env_match is not None:
            env_values[lmod_env_match.group("name")] = lmod_env_match.group(
                "value"
            )
            continue

        lmod_path_match = _LMOD_PREPEND_PATH_PATTERN.search(line)
        if lmod_path_match is not None:
            path_entries.append(lmod_path_match.group("value"))
            continue

        parsed = _parse_tcl_module_show_line(line)
        if parsed is None:
            continue

        directive, args = parsed
        if directive == "setenv" and len(args) >= 2:
            env_values[args[0]] = args[1]
            continue

        if (
            directive in {"prepend-path", "prepend_path"}
            and len(args) >= 2
            and args[0] == "PATH"
        ):
            path_entries.append(args[1])

    return env_values, path_entries


def _parse_tcl_module_show_line(
    line: str,
) -> tuple[str, list[str]] | None:
    if not line.startswith(("setenv ", "prepend-path ", "prepend_path ")):
        return None
    try:
        parts = shlex.split(line, posix=True)
    except ValueError:
        return None
    if not parts:
        return None
    return parts[0], parts[1:]


def _probe_conda_base(
    runner,
    topology: Topology,
    conda_path: str,
) -> str | None:
    spec_name = (
        "software.conda_base"
        if conda_path == "conda"
        else "software.conda_base_at_path"
    )
    slots = {} if conda_path == "conda" else {"conda_path": conda_path}
    result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS[spec_name],
        **slots,
    )
    return _normalize_shell_value(_first_nonempty_line(result.stdout))


def _probe_conda_env_list(
    runner,
    topology: Topology,
    conda_path: str,
) -> list[tuple[str | None, str, bool]]:
    if conda_path == "conda":
        return []

    result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["software.conda_env_list_at_path"],
        conda_path=conda_path,
    )
    merged_output = _merge_output(result)
    if result.returncode != 0 or not merged_output.strip():
        return []
    return _parse_conda_env_list(merged_output)


def _parse_conda_env_list(payload: str) -> list[tuple[str | None, str, bool]]:
    envs: list[tuple[str | None, str, bool]] = []
    for raw_line in payload.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = raw_line.split()
        if not parts:
            continue
        path_index = next(
            (
                index
                for index, value in enumerate(parts)
                if value.startswith(("/", "~"))
            ),
            None,
        )
        if path_index is None:
            continue
        path = _normalize_shell_value(parts[path_index])
        if path is None:
            continue
        prefix = parts[:path_index]
        active = "*" in prefix
        names = [value for value in prefix if value != "*"]
        name = names[0] if names else None
        envs.append((name, path, active))
    return envs


def _select_conda_env_path(
    base: str | None,
    envs: list[tuple[str | None, str, bool]],
) -> str | None:
    non_base = [entry for entry in envs if entry[1] != base]
    active_non_base = [path for _, path, active in non_base if active]
    if active_non_base:
        return active_non_base[0]

    named_chemsmart = [
        path for name, path, _ in non_base if name == "chemsmart"
    ]
    if named_chemsmart:
        return named_chemsmart[0]

    if len(non_base) == 1:
        return non_base[0][1]

    active_base = [path for _, path, active in envs if active]
    if active_base:
        return active_base[0]

    return None


def _expand_home_path(
    path: str,
    home: str | None,
    topology: Topology,
) -> str:
    if not path.startswith("~/"):
        return path
    if home:
        return str(Path(home) / path[2:])
    if topology.mode == "A":
        return str(Path.home() / path[2:])
    return path


def _merge_output(result) -> str:
    stdout = result.stdout.strip()
    stderr = result.stderr.strip()
    if stdout and stderr:
        return f"{stdout}\n{stderr}"
    return stdout or stderr


def _derive_conda_env_name(
    base: str | None,
    env_path: str | None,
) -> str | None:
    if env_path is None or env_path == base:
        return None

    env = Path(env_path)
    if env.parent.name != "envs":
        return None
    return env.name or None


def _verify_orca_qchem(runner, topology: Topology, exe_path: str) -> bool:
    """Return True only if exe_path is the quantum chemistry ORCA.

    The Linux screen reader (/usr/bin/orca) emits PyGI/GTK import warnings
    on stderr.  ORCA QChem prints "Program Version X.X.X - RELEASE -" or
    similar.  We reject any binary whose combined output contains the GTK
    screen-reader fingerprint.
    """
    try:
        import subprocess as _sp

        r = _sp.run(
            [exe_path, "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        combined = r.stdout + r.stderr
        # Fingerprint of the Linux accessibility screen reader
        if "PyGIWarning" in combined or "require_version" in combined:
            return False
        output_lower = combined.lower()
        # ORCA QChem version output contains "program version" or "release"
        return "orca" in output_lower and (
            "program version" in output_lower
            or "release" in output_lower
            or bool(__import__("re").search(r"\d+\.\d+\.\d+", combined))
        )
    except Exception:
        return False
