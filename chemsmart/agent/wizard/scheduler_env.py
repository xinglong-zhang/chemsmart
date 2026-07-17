"""Auto-detect HPC scheduler installations and build an augmented env dict.

Called once per process (LRU-cached) so that ProbeRunner subprocesses always
receive a PATH and the scheduler-specific env vars they need, even when the
wizard is launched from a non-login shell or a shell whose profile does not
source the scheduler settings file.

Supported schedulers (detected in order of probe priority):
  * SGE / UGE (Sun/Univa Grid Engine)  — looks for util/arch + settings.sh
  * SLURM                               — looks for sinfo binary in common paths
  * PBS / Torque                        — looks for qstat -Q binary
"""

from __future__ import annotations

import logging
import os
import subprocess
from functools import lru_cache
from pathlib import Path

logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------- #
# Candidate paths                                                              #
# --------------------------------------------------------------------------- #

_SGE_ROOT_CANDIDATES: list[str] = [
    "/opt/sge",
    "/opt/uge",
    "/opt/univa",
    "/usr/sge",
    "/usr/local/sge",
    "/cm/shared/apps/sge",
    "/cm/shared/apps/uge",
    "/N/soft/sge",
    "/N/u/sge",
    "/scratch/sge",
    "/cm/local/apps/uge/current",
    "/cm/local/apps/sge/current",
    "/cm/shared/apps/uge/current",
    "/cm/shared/apps/sge/current",
    "/opt/gridengine",
]

_SLURM_BIN_CANDIDATES: list[str] = [
    "/usr/bin",
    "/usr/local/bin",
    "/usr/local/slurm/bin",
    "/opt/slurm/bin",
    "/cm/shared/apps/slurm/current/bin",
    "/apps/slurm/bin",
    "/cm/local/apps/slurm/current/bin",
    "/opt/slurm/current/bin",
    "/opt/ohpc/pub/apps/slurm/current/bin",
]

_PBS_BIN_CANDIDATES: list[str] = [
    "/usr/bin",
    "/usr/pbs/bin",
    "/opt/pbs/bin",
    "/opt/torque/bin",
    "/usr/local/pbs/bin",
    "/cm/shared/apps/torque/current/bin",
    "/opt/pbs/default/bin",
    "/var/spool/pbs/bin",
    "/cm/local/apps/pbspro/current/bin",
]

_GAUSSIAN_BIN_CANDIDATES: list[str] = [
    "/opt/gaussian/g16",
    "/opt/gaussian/g09",
    "/opt/g16",
    "/opt/g09",
    "/usr/local/gaussian/g16",
    "/usr/local/gaussian/g09",
    "/usr/local/g16",
    "/usr/local/g09",
    "/cm/shared/apps/gaussian/g16",
    "/cm/shared/apps/gaussian/g09",
    "/apps/gaussian/g16",
    "/apps/gaussian/g09",
    "/software/gaussian/g16",
    "/software/gaussian/g09",
]

_ORCA_BIN_CANDIDATES: list[str] = [
    "/opt/orca",
    "/opt/orca/bin",
    "/usr/local/orca",
    "/usr/local/orca/bin",
    "/cm/shared/apps/orca/current",
    "/apps/orca",
    "/software/orca",
]


# --------------------------------------------------------------------------- #
# SGE helpers                                                                  #
# --------------------------------------------------------------------------- #


def _sge_arch(sge_root: str) -> str | None:
    arch_bin = Path(sge_root) / "util" / "arch"
    if not arch_bin.exists():
        return None
    try:
        r = subprocess.run(
            [str(arch_bin)],
            capture_output=True,
            text=True,
            timeout=5,
        )
        arch = r.stdout.strip()
        return arch if arch else None
    except Exception:
        return None


def _find_sge() -> tuple[str, str, str] | None:
    """Return (sge_root, sge_cell, arch) or None."""
    # 1. Already set in current environment
    env_root = os.environ.get("SGE_ROOT", "").strip()
    env_cell = os.environ.get("SGE_CELL", "default").strip() or "default"
    if env_root and Path(env_root).is_dir():
        arch = _sge_arch(env_root)
        if arch:
            return env_root, env_cell, arch

    # 2. Probe candidate roots
    for root in _SGE_ROOT_CANDIDATES:
        root_path = Path(root)
        if not root_path.is_dir():
            continue
        arch = _sge_arch(root)
        if not arch:
            continue
        # Find cell (default first, then any subdirectory with common/settings.sh)
        for cell in ("default",):
            if (root_path / cell / "common" / "settings.sh").exists():
                return root, cell, arch
        for child in root_path.iterdir():
            if (child / "common" / "settings.sh").exists():
                return root, child.name, arch

    return None


def _parse_sge_settings(sge_root: str, sge_cell: str) -> dict[str, str]:
    """Parse KEY=value pairs from settings.sh without executing it."""
    settings = Path(sge_root) / sge_cell / "common" / "settings.sh"
    result: dict[str, str] = {}
    if not settings.exists():
        return result
    try:
        for raw in settings.read_text(errors="replace").splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            # Each meaningful line is one of:
            #   KEY=value; export KEY
            #   export KEY=value
            #   KEY=value
            # Split on ';' first to discard the trailing 'export KEY' clause
            line = line.split(";")[0].strip()
            line = line.removeprefix("export ").strip()
            if "=" not in line:
                continue
            key, _, val = line.partition("=")
            key = key.strip()
            val = val.strip().strip('"').strip("'")
            # Only accept clean identifier keys; skip shell variable expansions
            if key and all(c.isalnum() or c == "_" for c in key):
                if "$" not in val:
                    result[key] = val
    except OSError:
        pass
    return result


def _build_sge_env(base: dict[str, str]) -> dict[str, str]:
    env = dict(base)
    found = _find_sge()
    if not found:
        return env

    sge_root, sge_cell, arch = found
    settings = _parse_sge_settings(sge_root, sge_cell)

    # Inject SGE env vars (don't overwrite caller-set values)
    for key, val in settings.items():
        if key.startswith("SGE_") or key in ("DRMAA_LIBRARY_PATH",):
            env.setdefault(key, val)
    env.setdefault("SGE_ROOT", sge_root)
    env.setdefault("SGE_CELL", sge_cell)

    # Prepend SGE bin directories to PATH
    extra = [
        str(Path(sge_root) / "bin" / arch),
        str(Path(sge_root) / "bin"),
    ]
    current_path = env.get("PATH", "")
    path_parts = current_path.split(":")
    additions = [p for p in extra if p and p not in path_parts]
    if additions:
        env["PATH"] = ":".join(additions) + ":" + current_path

    logger.debug(
        "scheduler_env: SGE detected at %s (cell=%s arch=%s); "
        "prepended %s to PATH",
        sge_root,
        sge_cell,
        arch,
        additions,
    )
    return env


# --------------------------------------------------------------------------- #
# SLURM helpers                                                                #
# --------------------------------------------------------------------------- #


def _find_slurm_bin() -> str | None:
    for directory in _SLURM_BIN_CANDIDATES:
        if (Path(directory) / "sinfo").exists():
            return directory
    return None


def _build_slurm_env(base: dict[str, str]) -> dict[str, str]:
    env = dict(base)
    slurm_bin = _find_slurm_bin()
    if not slurm_bin:
        return env
    current_path = env.get("PATH", "")
    if slurm_bin not in current_path.split(":"):
        env["PATH"] = slurm_bin + ":" + current_path
        logger.debug(
            "scheduler_env: SLURM bin at %s prepended to PATH", slurm_bin
        )
    return env


# --------------------------------------------------------------------------- #
# PBS / Torque helpers                                                         #
# --------------------------------------------------------------------------- #


def _find_pbs_bin() -> str | None:
    for directory in _PBS_BIN_CANDIDATES:
        candidate = Path(directory) / "qstat"
        if candidate.exists():
            # Distinguish PBS qstat from SGE qstat by checking for -Q flag support
            try:
                r = subprocess.run(
                    [str(candidate), "-Q"],
                    capture_output=True,
                    text=True,
                    timeout=5,
                )
                # PBS qstat -Q lists queues; SGE qstat -Q shows job state
                if "Queue" in r.stdout or r.returncode == 0:
                    return directory
            except Exception:
                pass
    return None


def _build_pbs_env(base: dict[str, str]) -> dict[str, str]:
    env = dict(base)
    pbs_bin = _find_pbs_bin()
    if not pbs_bin:
        return env
    current_path = env.get("PATH", "")
    if pbs_bin not in current_path.split(":"):
        env["PATH"] = pbs_bin + ":" + current_path
        logger.debug("scheduler_env: PBS bin at %s prepended to PATH", pbs_bin)
    return env


# --------------------------------------------------------------------------- #
# Gaussian helpers                                                             #
# --------------------------------------------------------------------------- #


def _find_gaussian_bin() -> str | None:
    for directory in _GAUSSIAN_BIN_CANDIDATES:
        for exe in ("g16", "g09"):
            if (Path(directory) / exe).exists():
                return directory
    return None


def _build_gaussian_env(base: dict[str, str]) -> dict[str, str]:
    env = dict(base)
    gauss_bin = _find_gaussian_bin()
    if not gauss_bin:
        return env
    current_path = env.get("PATH", "")
    if gauss_bin not in current_path.split(":"):
        env["PATH"] = gauss_bin + ":" + current_path
        gauss_root = str(Path(gauss_bin).parent)
        env.setdefault("g16root", gauss_root)
        env.setdefault("g09root", gauss_root)
        env.setdefault("GAUSS_SCRDIR", "/tmp")
        logger.debug(
            "scheduler_env: Gaussian bin at %s prepended to PATH", gauss_bin
        )
    return env


# --------------------------------------------------------------------------- #
# ORCA (quantum chemistry) helpers                                             #
# --------------------------------------------------------------------------- #


def _find_orca_qchem_bin() -> str | None:
    """Find ORCA quantum chemistry binary (distinct from Linux orca screen reader)."""
    for directory in _ORCA_BIN_CANDIDATES:
        candidate = Path(directory) / "orca"
        if not candidate.exists():
            continue
        try:
            r = subprocess.run(
                [str(candidate), "--version"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            output = (r.stdout + r.stderr).lower()
            if "orca" in output and (
                "version" in output or "release" in output
            ):
                return directory
        except Exception:
            pass
    return None


def _build_orca_qchem_env(base: dict[str, str]) -> dict[str, str]:
    env = dict(base)
    orca_bin = _find_orca_qchem_bin()
    if not orca_bin:
        return env
    current_path = env.get("PATH", "")
    if orca_bin not in current_path.split(":"):
        env["PATH"] = orca_bin + ":" + current_path
        logger.debug(
            "scheduler_env: ORCA (qchem) bin at %s prepended to PATH", orca_bin
        )
    return env


# --------------------------------------------------------------------------- #
# Public API                                                                   #
# --------------------------------------------------------------------------- #


@lru_cache(maxsize=1)
def build_scheduler_env() -> dict[str, str]:
    """Return os.environ augmented with any auto-detected scheduler/program paths.

    Cached after the first call — safe to call repeatedly from probes.
    Detection order: SGE → SLURM → PBS → Gaussian → ORCA.
    """
    env = dict(os.environ)
    env = _build_sge_env(env)
    env = _build_slurm_env(env)
    env = _build_pbs_env(env)
    env = _build_gaussian_env(env)
    env = _build_orca_qchem_env(env)
    return env
