"""Live, non-computing conformance and environment bootstrapping.

These probes exist so an agent session can observe the checked-out source tree
and the exact local program environment before a model proposes work.  They do
not run an SCF cycle and they do not turn source presence into scientific
readiness.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Iterable

from click.testing import CliRunner

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.capabilities import (
    ProgramComponentConformanceReceiptV1,
    TrustedComputeEnvironmentReceiptV1,
    build_program_component_conformance_receipt,
    build_trusted_compute_environment_receipt,
)
from chemsmart.agent.cli_schema import LiveClickSchemaV1
from chemsmart.agent.commands import native_coordinate_options

_PYTHON_PROBE = r"""
import importlib
import importlib.metadata
import json
import pathlib
import sys

names = ("pyscf", "numpy", "h5py", "geometric", "gpu4pyscf", "cupy")
versions = {}
for name in names:
    try:
        module = importlib.import_module(name)
        version = getattr(module, "__version__", "")
        if not version:
            version = importlib.metadata.version(name)
        versions[name] = str(version or "available")
    except Exception:
        pass
print(json.dumps({
    "interpreter": str(pathlib.Path(sys.executable).resolve()),
    "versions": versions,
}, sort_keys=True))
"""


def probe_python_compute_environment(
    interpreter: str | Path,
    *,
    engine: str,
    timeout_seconds: float = 60.0,
) -> TrustedComputeEnvironmentReceiptV1:
    """Observe the exact Python used for PySCF CPU or GPU execution."""

    executable = Path(interpreter).expanduser().resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ContractError("PySCF compute interpreter is not executable")
    completed = subprocess.run(
        [str(executable), "-c", _PYTHON_PROBE],
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
        check=False,
    )
    if completed.returncode != 0:
        raise ContractError("PySCF compute environment probe failed")
    try:
        observed = json.loads(completed.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as exc:
        raise ContractError(
            "PySCF compute environment probe was malformed"
        ) from exc
    if Path(observed.get("interpreter", "")).resolve() != executable:
        raise ContractError(
            "PySCF compute interpreter identity changed during probe"
        )
    versions = {
        str(name).lower(): str(value)
        for name, value in dict(observed.get("versions") or {}).items()
    }
    required = {"pyscf", "numpy", "h5py"}
    if not required.issubset(versions):
        raise ContractError(
            "PySCF compute environment lacks required packages"
        )
    gpu_evidence = {
        "device_available": False,
        "gpu4pyscf_distribution": "gpu4pyscf" in versions,
        "gpu4pyscf_version": versions.get("gpu4pyscf", ""),
        "cupy_version": versions.get("cupy", ""),
    }
    return build_trusted_compute_environment_receipt(
        program="pyscf",
        engine=engine,
        compute_interpreter_sha256=file_sha256(executable),
        dependency_versions=versions,
        solver_evidence={"geometric": "geometric" in versions},
        gpu_evidence=gpu_evidence,
        source_probe_id=f"live-python-probe:{canonical_sha256(observed)}",
    )


#: The coordinate a conformance probe drives or holds.  The first two atoms of
#: the fixture geometry are enough to reach the same Click callback a real user
#: command reaches; nothing here computes, so the values only have to be legal.
_CONFORMANCE_COORDINATES = {
    "modred": {"constrained": [{"kind": "bond", "atoms": [1, 2]}]},
    "scan": {
        "scan": {
            "kind": "bond",
            "atoms": [1, 2],
            "start": 1.0,
            "stop": 2.0,
            "points": 3,
        }
    },
}


def bootstrap_program_conformance(
    *,
    program: str,
    engine: str,
    jobtypes: Iterable[str],
    input_path: str | Path,
    project_path: str | Path | None,
    server_path: str | Path | None = None,
    charge: int | None = None,
    multiplicity: int | None = None,
    registry_sha256: str,
    live_schema: LiveClickSchemaV1,
) -> ProgramComponentConformanceReceiptV1:
    """Exercise real Click fake previews before enabling an agent program.

    Both ``run`` and ``sub`` paths are parsed for every declared stage.  The
    resulting receipt enables planning and safe preview only; it deliberately
    carries no execution evidence.
    """

    input_path = Path(input_path).resolve()
    project = Path(project_path).resolve() if project_path else None
    server = Path(server_path).resolve() if server_path else None
    if (
        not input_path.is_file()
        or (project is not None and not project.is_file())
        or (server is not None and not server.is_file())
    ):
        raise ContractError(
            "conformance fixtures must be existing regular files"
        )
    if (charge is None) != (multiplicity is None):
        raise ContractError(
            "charge and multiplicity must be supplied together"
        )
    observations = []
    covered = []
    runner = CliRunner()
    for jobtype in tuple(sorted(set(str(item) for item in jobtypes))):
        schema_green = all(
            live_schema.has_path((target, program, jobtype))
            for target in ("run", "sub")
        )
        stage_green = schema_green
        stage_rows = []
        if schema_green:
            for target in ("run", "sub"):
                argv = [target]
                if server is not None:
                    argv.extend(("--server", str(server)))
                argv.extend(("--fake", "--no-scratch"))
                if target == "sub":
                    argv.append("--test")
                if program == "pyscf" and engine == "gpu":
                    argv.extend(("--num-gpus", "1"))
                argv.append(program)
                if project is not None:
                    argv.extend(("--project", str(project)))
                argv.extend(("--filename", str(input_path)))
                if charge is not None:
                    argv.extend(("--charge", str(charge)))
                    argv.extend(("--multiplicity", str(multiplicity)))
                if program == "pyscf":
                    argv.append("--gpu" if engine == "gpu" else "--no-gpu")
                argv.append(jobtype)
                # Multi-file stages need a second, existing fixture to reach
                # the same canonical Click callback as a real user command.
                # The live schema remains the source of the option spelling;
                # using the primary geometry as both ends is sufficient for a
                # non-computing syntax/materialization preview.
                if program == "orca" and jobtype == "neb":
                    job_command = live_schema.command(
                        (target, program, jobtype)
                    )
                    ending_option = (
                        job_command.option("ending_xyzfile")
                        if job_command is not None
                        else None
                    )
                    if ending_option is None:
                        stage_green = False
                    else:
                        argv.extend(
                            (ending_option.primary_flag, str(input_path))
                        )
                # A scan or a constrained optimisation is defined by the
                # coordinate it drives or holds, and its Click callback
                # refuses to build a job without one.  Render the fixture
                # coordinate through the same translation the production path
                # uses, so conformance exercises that translation rather than
                # a hand-written argv that could drift from it.
                if jobtype in {"modred", "scan"}:
                    job_command = live_schema.command(
                        (target, program, jobtype)
                    )
                    values = native_coordinate_options(
                        program, _CONFORMANCE_COORDINATES[jobtype]
                    )
                    for parameter_name, value in sorted(values.items()):
                        option = (
                            job_command.option(parameter_name)
                            if job_command is not None
                            else None
                        )
                        if option is None:
                            stage_green = False
                            break
                        argv.extend((option.primary_flag, str(value)))
                with runner.isolated_filesystem():
                    from chemsmart.cli.main import entry_point

                    result = runner.invoke(entry_point, argv)
                stage_rows.append(
                    {
                        "target": target,
                        "argv_shape": tuple(
                            (
                                "<input>"
                                if value == str(input_path)
                                else (
                                    "<project>"
                                    if project is not None
                                    and value == str(project)
                                    else (
                                        "<server>"
                                        if server is not None
                                        and value == str(server)
                                        else value
                                    )
                                )
                            )
                            for value in argv
                        ),
                        "exit_code": int(result.exit_code),
                        "exception": (
                            type(result.exception).__name__
                            if result.exception
                            else ""
                        ),
                    }
                )
                stage_green = stage_green and result.exit_code == 0
        observations.append(
            {
                "jobtype": jobtype,
                "schema_green": schema_green,
                "paths": tuple(stage_rows),
            }
        )
        if stage_green:
            covered.append(jobtype)
    fixture = {
        "program": program,
        "engine": engine,
        "input_sha256": file_sha256(input_path),
        "project_sha256": file_sha256(project) if project is not None else "",
        "server_sha256": file_sha256(server) if server is not None else "",
        "charge": charge,
        "multiplicity": multiplicity,
        "observations": tuple(observations),
    }
    fixture_sha = canonical_sha256(fixture)
    passed = bool(covered) and len(covered) == len(
        tuple(sorted(set(jobtypes)))
    )
    status = "passed" if passed else "failed"
    compiler_sha = canonical_sha256(
        {"schema": live_schema.schema_sha256, "paths": tuple(observations)}
    )
    preview_sha = canonical_sha256(fixture)
    # The fake CLI callback performs settings/materialisation preflight.  It
    # does not invoke the separate generated-artifact verifier: that verifier
    # runs later for the concrete node and must never be reported as observed
    # merely because its module is importable.
    preflight_sha = canonical_sha256(
        {"component": "program_node_preflight", "program": program}
    )
    return build_program_component_conformance_receipt(
        program=program,
        registry_sha256=registry_sha256,
        live_cli_schema_sha256=live_schema.schema_sha256,
        fixture_bundle_sha256=fixture_sha,
        covered_jobtypes=tuple(covered),
        covered_engines=(engine,) if passed else (),
        compiler_receipt_sha256=compiler_sha,
        preview_receipt_sha256=preview_sha,
        preflight_receipt_sha256=preflight_sha,
        verifier_receipt_sha256="",
        compiler_status=status,
        preview_status=status,
        preflight_status=status,
        verifier_status="not_observed",
    )


__all__ = [
    "bootstrap_program_conformance",
    "probe_python_compute_environment",
]
